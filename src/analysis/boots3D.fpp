!=================================================================
      PROGRAM BOOTS3D
!=================================================================
! BOOTS3D code (part of the GHOST suite)
!
! Reads dataset, and computes dataset that is prolongated in 
! wavenumber space, and outputs to real space. The files
! to be read are specified on input, as is the linear size of the
! old grid; the output files are the originial files appended with
! the prolongation size 'PXXX-XXX-XXX'.
! This is the 'boostrap regridding' procedure (BOOTS).
!
! Note that the number of mpi tasks must be of the form 2^i * 
! nx_new / nx_old, with i .ne. 1. For example, for nx_new = 384 and
! nx_old = 128 possible num of tasks are 1, 6, 12, 24, 48, etc.
!
! NOTATION: index 'i' is 'x' 
!           index 'j' is 'y'
!           index 'k' is 'z'
!
! 2010 D. Rosenberg
!      NCAR
!
! 17 Nov 2010: Initial version
! 15 Mar 2019: Added support for anisotropic boxes and
!              error checking for number of tasks (M. Fontana)
!=================================================================

!
! Definitions for conditional compilation

! Modules

      USE fprecision
      USE commtypes
      USE mpivars
      USE filefmt
      USE iovar
      USE iompi
      USE grid
      USE fft
      USE ali
      USE var
      USE kes
      USE fftplans
      USE threads
      USE gutils
      IMPLICIT NONE

!
! Arrays for the fields

      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: C1t
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: B1

      REAL(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: vvt
      REAL(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: br
      REAL(KIND=GP)                                 :: fact, rmp, rmq, rms
      
!
! Auxiliary variables

      INTEGER :: nfiles,nmt,np,npt,ntprocs,npkeep
      INTEGER :: nxt,nyt,nzt
      INTEGER :: i,ib,ie,ind,itsta,itend,j,k,ktsta,ktend
      INTEGER :: istak,iendk,kstak,kendk
      INTEGER :: commtrunc, fh, groupworld, flags, grouptrunc
      INTEGER :: iExclude(3,1), iInclude(3,1)
 
      TYPE(IOPLAN)  :: planio, planiot
      TYPE(FFTPLAN) :: planrct

      CHARACTER(len=19)  :: suff
      CHARACTER(len=100) :: odir,idir
      CHARACTER(len=256) :: fname, fout, msg
      CHARACTER(len=1024):: fnlist

      NAMELIST / regrid / idir, odir, fnlist, iswap, oswap, nxt, nyt, nzt

! Initializes the MPI and I/O libraries considering all the tasks
! and the biggest possible array size (i.e. the new dimensions)
      CALL MPI_INIT(ierr)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
      CALL range(1,nx/2+1,nprocs,myrank,ista,iend)
      CALL range(1,nz,nprocs,myrank,ksta,kend)
      CALL io_init(myrank,(/nx,ny,nz/),ksta,kend,planio)
      idir   = '.'
      odir   = '.'
      fnlist = ''
      iswap  = 0
      oswap  = 0
      nxt    = 0
      nyt    = 0
      nzt    = 0
!
! Reads from the external file 'boots.txt' the 
! parameters that will be used to compute the transfer
!     idir   : directory for unformatted input (field components)
!     odir   : directory for unformatted output (prolongated data)
!     fnlist : file list to prolongate, separated by ';'
!     iswap  : do endian swap on input?
!     oswap  : do endian swap on output?
!     nxt    : original linear size of the old grid in x direction
!     nyt    : original linear size of the old grid in y direction
!     nzt    : original linear size of the old grid in z direction

      IF (myrank.eq.0) THEN
         OPEN(1,file='boots.inp',status='unknown',form="formatted")
         READ(1,NML=regrid)
         CLOSE(1)
      ENDIF
      CALL MPI_BCAST(idir  ,100 ,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(odir  ,100 ,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(fnlist,1024,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(iswap ,1   ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(oswap ,1   ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(nxt   ,1   ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(nyt   ,1   ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(nzt   ,1   ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)

      !
      !  n here refers to the prolongated grid, and nt
      !  is the grid size of the 'old grid'

      ! Check input
      IF ( nxt .GT. nx .OR. nxt .LT. 1 ) THEN
        IF ( myrank .eq. 0) PRINT*, 'MAIN: prolongation specification incorrect; input nxt must be less than Nx'
        CALL MPI_Finalize(ierr)
        STOP  
      ENDIF
      IF ( nyt .GT. ny .OR. nyt .LT. 1 ) THEN
        IF ( myrank .eq. 0) PRINT*, 'MAIN: prolongation specification incorrect; input nyt must be less than Ny'
        CALL MPI_Finalize(ierr)
        STOP
      ENDIF
      IF ( nzt .GT. nz .OR. nzt .LT. 1 ) THEN
        IF ( myrank .eq. 0) PRINT*, 'MAIN: prolongation specification incorrect; input nzt must be less than Nz'
        CALL MPI_Finalize(ierr)
        STOP
      ENDIF

      ! First check that nprocs is not 1, and that nx is not nxt.
      ! After that catch non-integer number of slices,
      ! then check quotient is power of two, finally exclude the 2^1 case
      IF ( nprocs .ne. 1 .AND. nxt .ne. nx .AND. &
            ( mod(nprocs*nxt, nx) .ne. 0 .OR. &
              IAND(nprocs*nxt/nx, nprocs*nxt/nx - 1) .ne. 0 .OR. &
              nprocs .eq. nx/nxt ) ) THEN
        IF ( myrank .eq. 0) PRINT*, 'MAIN: number of tasks must be of the form 2^i * (nx_new/nx_old) with i =/= 1'
        CALL MPI_Finalize(ierr)
        STOP
      ENDIF

      ! Suffix for the output
      WRITE(suff,'(a2,i5.5,a1,i5.5,a1,i5.5)') '_P', nx, "-", ny, "-", nz
     
! Determine size and create a new communicator for the truncated data:
      np      = nx / nprocs
      ntprocs = nxt / np
      nmt     = mod(nxt,np) 
      IF ( nmt .GT. 0 ) ntprocs = ntprocs + 1
      ntprocs = min(ntprocs, nprocs)

      CALL MPI_COMM_GROUP(MPI_COMM_WORLD, groupworld, ierr)
      commtrunc  = MPI_COMM_NULL
      grouptrunc = MPI_GROUP_NULL
      IF ( ntprocs .LT. nprocs ) THEN
        iExclude(1,1) = ntprocs
        iExclude(2,1) = nprocs-1
        iExclude(3,1) = 1
        CALL MPI_COMM_GROUP(MPI_COMM_WORLD, groupworld, ierr)
        commtrunc  = MPI_COMM_NULL
        CALL MPI_GROUP_RANGE_EXCL(groupworld, 1, iExclude, grouptrunc, ierr)   
        CALL MPI_COMM_GROUP(MPI_COMM_WORLD, groupworld, ierr)
        commtrunc  = MPI_COMM_NULL
        CALL MPI_COMM_CREATE(MPI_COMM_WORLD, grouptrunc, commtrunc, ierr)
      ELSE IF ( ntprocs .EQ. nprocs ) THEN
        CALL MPI_COMM_DUP(MPI_COMM_WORLD,commtrunc,ierr)
        CALL MPI_COMM_GROUP(MPI_COMM_WORLD,grouptrunc,ierr)
      ENDIF
      CALL range(1,nzt    ,ntprocs,myrank,ktsta,ktend)
      CALL range(1,nxt/2+1,ntprocs,myrank,itsta,itend)

! Allocate arrays and initialize I/O and FFTs

      ALLOCATE( vvt(nxt,nyt,ktsta:ktend) )
      ALLOCATE( C1t(nzt,nyt,itsta:itend) )
      ALLOCATE( B1(nz,ny,ista:iend) )
      ALLOCATE( kx(nx), ky(ny), kz(nz) )
      ALLOCATE( kn2(nz,ny,ista:iend) )
      ALLOCATE( br(nx,ny,ksta:kend) )

      IF ( myrank .LT. ntprocs ) THEN
        npkeep = nprocs; nprocs = ntprocs
        CALL io_init(myrank,(/nxt,nyt,nzt/),ktsta,ktend,planiot)
        nprocs = npkeep
      ENDIF

      flags  = FFTW_MEASURE

      CALL fftp3d_create_plan(plancr,(/nx,ny,nz/),FFTW_COMPLEX_TO_REAL, flags)

      IF ( myrank .LT. ntprocs ) THEN
        ! Temporarily redefine ksta,ista,kend,iend to accomodate
        ! for ntprocs tasks. Then return to the standard definition
        npkeep = nprocs; nprocs = ntprocs
        CALL range(1,nzt    ,ntprocs,myrank,ksta,kend)
        CALL range(1,nxt/2+1,ntprocs,myrank,ista,iend)
        CALL fftp3d_create_plan(planrct,(/nxt,nyt,nzt/),FFTW_REAL_TO_COMPLEX,flags)
        nprocs = npkeep
        CALL range(1,nx/2+1,nprocs,myrank,ista,iend)
        CALL range(1,nz,nprocs,myrank,ksta,kend)
      ENDIF

! Some constants for the FFT
!     kmax: maximum truncation for dealiasing
!     tiny: minimum truncation for dealiasing
      kmax =     1.0_GP/9.0_GP

! Populate wavenumber-associated arrays
      DO i = 1,nx/2
         kx(i) = real(i-1,kind=GP)
         kx(i+nx/2) = real(i-nx/2-1,kind=GP)
      END DO
      DO j = 1,ny/2
         ky(j) = real(j-1,kind=GP)
         ky(j+ny/2) = real(j-ny/2-1,kind=GP)
      END DO
      DO k = 1,nz/2
         kz(k) = real(k-1,kind=GP)
         kz(k+nz/2) = real(k-nz/2-1,kind=GP)
      END DO
      rmp = 1.0_GP/real(nx,kind=GP)**2
      rmq = 1.0_GP/real(ny,kind=GP)**2
      rms = 1.0_GP/real(nz,kind=GP)**2

!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
            DO k = 1,nz
               kn2(k,j,i) = rmp*kx(i)**2+rmq*ky(j)**2+rms*kz(k)**2
            END DO
         END DO
      END DO

      ib = 1
      ie = len(fnlist)
      DO WHILE ( len(trim(fnlist(ib:ie))) .GT. 0 )  

         ind = index(fnlist(ib:ie),";")
         IF ( ind .eq. 0 ) THEN
           fname = trim(adjustl(fnlist(ib:ie)))
           ib = ie + 1
         ELSE
           fname = trim(adjustl(fnlist(ib:(ib+ind-2))))
           ib = ib + ind
         ENDIF
!
! Read the external binary files in real space:
         IF ( myrank .LT. ntprocs ) THEN
           CALL MPI_FILE_OPEN(commtrunc,trim(idir) // '/' // trim(fname), &
                 MPI_MODE_RDONLY,MPI_INFO_NULL,fh,ioerr)
           IF ( ioerr .NE. MPI_SUCCESS ) THEN
             WRITE(*,*) 'main: Error opening file: ', trim(idir) // '/' // trim(fname)
             STOP
           ENDIF
           CALL MPI_FILE_SET_VIEW(fh,disp,GC_REAL,planiot%iotype,'native', &
               MPI_INFO_NULL,ioerr)
           IF ( ioerr .NE. MPI_SUCCESS ) THEN
             WRITE(*,*) 'main: Error with MPI_FILE_SET_VIEW'
             STOP
           ENDIF
           CALL MPI_FILE_READ_ALL(fh,vvt, &
               planiot%nx*planiot%ny*(planiot%kend-planiot%ksta+1),GC_REAL, &
               MPI_STATUS_IGNORE,ioerr)
           IF ( ioerr .NE. MPI_SUCCESS ) THEN
             WRITE(*,*) 'main: Error with MPI_FILE_READ_ALL'
             STOP
           ENDIF
           CALL MPI_FILE_CLOSE(fh,ioerr)

           IF ( iswap .NE. 0 ) THEN
             CALL rarray_byte_swap(vvt, nxt*nyt*(ktend-ktsta+1))
           ENDIF

         ENDIF

!
! Compute FT of variable on smaller grid:
      IF ( myrank .LT. ntprocs ) THEN
        ! Temporarily redefine ksta,ista,kend,iend to accomodate
        ! for ntprocs tasks. Then return to the standard definition
         npkeep = nprocs; nprocs = ntprocs
         CALL range(1,nzt    ,ntprocs,myrank,ksta,kend)
         CALL range(1,nxt/2+1,ntprocs,myrank,ista,iend)
         CALL fftp3d_real_to_complex(planrct,vvt,C1t,commtrunc)
         nprocs = npkeep
         CALL range(1,nx/2+1,nprocs,myrank,ista,iend)
         CALL range(1,nz,nprocs,myrank,ksta,kend)
      ENDIF

!
! Prolongate in Fourier space:
         fact = 1.0_GP/ \
            (real(nxt,kind=GP)*real(nyt,kind=GP)*real(nzt,kind=GP))
         B1 = 0.0
         DO i = itsta,itend
            DO j = 1,nyt/2
               DO k = 1,nzt/2
                  B1(k,j,i) = C1t(k,j,i) * fact
               END DO
               DO k = nz-nzt/2+1,nz
                  B1(k,j,i) = C1t(k-nz+nzt,j,i) * fact
               END DO
            END DO
            DO j = ny-nyt/2+1,ny
               DO k = 1,nzt/2
                  B1(k,j,i) = C1t(k,j-ny+nyt,i) * fact
               END DO
               DO k = nz-nzt/2+1,nz
                  B1(k,j,i) = C1t(k-nz+nzt,j-ny+nyt,i) * fact
               END DO
            END DO
         END DO

#if 1
! Spherically truncate prolongated spectrum in Fourier space:
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
            DO j = 1,ny
               DO k = 1,nz
                  IF (  kn2(k,j,i).GT.kmax ) B1(k,j,i) = 0.0
               END DO
            END DO
         END DO
#endif

! Compute inverse FFT of prolongated variable:
         CALL fftp3d_complex_to_real(plancr,B1,br,MPI_COMM_WORLD)
!
! Byte-swap on output:
         IF ( oswap .NE. 0 ) THEN
           CALL rarray_byte_swap(br, nx*ny*(kend-ksta+1))
         ENDIF
!
! Put to disk:
         fout = trim(odir) // '/' // trim(fname) // trim(suff)
         bmangle = 0
         CALL io_write(1,idir,fout,1,planio,br)
         bmangle = 1
         IF ( myrank .EQ. 0 ) THEN
           WRITE(*,*) 'main: ',trim(fout),' written.'
         ENDIF
      ENDDO ! end of file loop
!
      CALL fftp3d_destroy_plan(plancr)
      IF ( myrank .LT. ntprocs ) THEN
        CALL fftp3d_destroy_plan(planrct)
        CALL MPI_COMM_FREE(commtrunc, ierr)
        CALL MPI_GROUP_FREE(grouptrunc, ierr)
      ENDIF

      DEALLOCATE ( vvt )
      DEALLOCATE ( C1t )
      DEALLOCATE ( B1 )
      DEALLOCATE ( kx, ky, kz )
      DEALLOCATE ( kn2 )
      DEALLOCATE ( br )

      CALL MPI_FINALIZE(ierr)

      END PROGRAM BOOTS3D
