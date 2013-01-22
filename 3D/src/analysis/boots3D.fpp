!=================================================================
      PROGRAM BOOTS3D
!=================================================================
! BOOTS3D code (part of the GHOST suite)
!
! Reads dataset, and computes dataset that is prolongated in 
! wavenumber space, and outputs to real space. The files
! to be read are specified on input, as is the max linear size of the
! old grid; the output files are the originial files appended with with the 
! prolongation size 'PXXX'. This is the 'boostrap regridding'
! procedure (BOOTS).
!
! NOTATION: index 'i' is 'x' 
!           index 'j' is 'y'
!           index 'k' is 'z'
!
! 2010 D. Rosenberg
!      NCAR
!
! 17 Nov 2010: Initial version
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
! Arrays for the fields and structure functions

      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: C1t
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: B1


      REAL(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: vvt
      REAL(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: br
      REAL(KIND=GP)                                 :: fact
!
! Auxiliary variables

      REAL(KIND=GP)    :: kprol2, kprol
      INTEGER :: nfiles,nmt,np,npt,nt,ntprocs,npkeep
      INTEGER :: i,ib,ie,ind,iswap,oswap,itsta,itend,j,k,ktsta,ktend
      INTEGER :: istak,iendk,kstak,kendk
      INTEGER :: commtrunc, fh, groupworld, flags, grouptrunc, iExclude(3,1), iInclude(3,1)

      TYPE(IOPLAN)  :: planio, planiot
      TYPE(FFTPLAN) :: planrct

      CHARACTER(len=8)   :: suff
      CHARACTER(len=100) :: odir,idir
      CHARACTER(len=256) :: fname, fout, msg
      CHARACTER(len=1024):: fnlist
!
      NAMELIST / regrid / idir, odir, fnlist, iswap, oswap, nt

!
! Initializes the MPI and I/O libraries
      CALL MPI_INIT(ierr)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
      CALL range(1,n/2+1,nprocs,myrank,ista,iend)
      CALL range(1,n,nprocs,myrank,ksta,kend)
      CALL io_init(myrank,n,ksta,kend,planio)
      idir   = '.'
      odir   = '.'
      fnlist = ''
      iswap  = 0
      oswap  = 0
      nt     = 0
!
! Reads from the external file 'boots.txt' the 
! parameters that will be used to compute the transfer
!     idir   : directory for unformatted input (field components)
!     odir   : directory for unformatted output (prolongated data)
!     fnlist : file list to prolongate, separated by ';'
!     iswap  : do endian swap on input?
!     oswap  : do endian swap on output?
!     nt     : original linear size of the old grid

      IF (myrank.eq.0) THEN
         OPEN(1,file='boots.txt',status='unknown',form="formatted")
         READ(1,NML=regrid)
         CLOSE(1)
      ENDIF
      CALL MPI_BCAST(idir  ,100 ,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(odir  ,100 ,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(fnlist,1024,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(iswap ,1   ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(oswap ,1   ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(nt    ,1   ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)

      IF ( nt .GT. n .OR. nt .LT. 1 ) THEN
        STOP 'MAIN: prolongation specification incorrect; input nt must be less than N' 
      ENDIF

      WRITE(suff,'(a2,i5.5)') '_P', n
!
!  n, here refers to the prolongated grid, and nt
!  is the grid size of the 'old grid'
     
! Create new communicator for the truncated data:
      np      = n / nprocs
      ntprocs = nt / np
      nmt     = mod(nt,np) 
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
      CALL range(1,nt    ,ntprocs,myrank,ktsta,ktend)
      CALL range(1,nt/2+1,ntprocs,myrank,itsta,itend)


!     kprol  = real(n,kind=GP)/2.0 - 1.0
      kprol  = real(n,kind=GP)/3.0
      kprol2 = (real(kprol,kind=GP) )**2
      kmax   = (real(n,kind=GP)/3.)**2

      ALLOCATE( vvt(nt,nt,ktsta:ktend) )
      ALLOCATE( C1t(nt,nt,itsta:itend) )
      ALLOCATE( B1(n,n,ista:iend) )
      ALLOCATE( ka(n),ka2(n,n,ista:iend) )
      ALLOCATE( br(n,n,ksta:kend) )

      IF ( myrank .LT. ntprocs ) THEN
        npkeep = nprocs; nprocs = ntprocs
        CALL io_init(myrank,nt,ktsta,ktend,planiot)
        nprocs = npkeep
      ENDIF
!
      flags  = FFTW_MEASURE

      CALL fftp3d_create_plan(plancr,n,FFTW_COMPLEX_TO_REAL, flags)

      IF ( myrank .LT. ntprocs ) THEN
        npkeep = nprocs; nprocs = ntprocs
        CALL range(1,nt    ,ntprocs,myrank,ksta,kend)
        CALL range(1,nt/2+1,ntprocs,myrank,ista,iend)
        CALL fftp3d_create_plan(planrct,nt,FFTW_REAL_TO_COMPLEX,flags)
        nprocs = npkeep
        CALL range(1,n/2+1,nprocs,myrank,ista,iend)
        CALL range(1,n,nprocs,myrank,ksta,kend)
      ENDIF
!
      DO i = 1,n/2
         ka(i) = real(i-1,kind=GP)
         ka(i+n/2) = real(i-n/2-1,kind=GP)
      END DO

      DO i = ista,iend
         DO j = 1,n
            DO k = 1,n
               ka2(k,j,i) = ka(i)**2+ka(j)**2+ka(k)**2
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
               planiot%n*planiot%n*(planiot%kend-planiot%ksta+1),GC_REAL, &
               MPI_STATUS_IGNORE,ioerr)
           IF ( ioerr .NE. MPI_SUCCESS ) THEN
             WRITE(*,*) 'main: Error with MPI_FILE_READ_ALL'
             STOP
           ENDIF
           CALL MPI_FILE_CLOSE(fh,ioerr)

           IF ( iswap .NE. 0 ) THEN
             CALL rarray_byte_swap(vvt, nt*nt*(ktend-ktsta+1))
           ENDIF

         ENDIF

!
! Compute FT of variable on smaller grid:
      IF ( myrank .LT. ntprocs ) THEN
         npkeep = nprocs; nprocs = ntprocs
         CALL range(1,nt    ,ntprocs,myrank,ksta,kend)
         CALL range(1,nt/2+1,ntprocs,myrank,ista,iend)
         CALL fftp3d_real_to_complex(planrct,vvt,C1t,commtrunc)
         nprocs = npkeep
         CALL range(1,n/2+1,nprocs,myrank,ista,iend)
         CALL range(1,n,nprocs,myrank,ksta,kend)
      ENDIF

!
! Prolongate in Fourier space:
         fact = 1.0_GP/real(nt,kind=GP)**3
         B1 = 0.0
         DO i = itsta,itend
            DO j = 1,nt/2
               DO k = 1,nt/2
                  B1(k,j,i) = C1t(k,j,i) * fact
               END DO
               DO k = n-nt/2+1,n
                  B1(k,j,i) = C1t(k-n+nt,j,i) * fact
               END DO
            END DO
            DO j = n-nt/2+1,n
               DO k = 1,nt/2
                  B1(k,j,i) = C1t(k,j-n+nt,i) * fact
               END DO
               DO k = n-nt/2+1,n
                  B1(k,j,i) = C1t(k-n+nt,j-n+nt,i) * fact
               END DO
            END DO
         END DO
!
!
#if 1
! Spherically truncate prolongated spectrum in Fourier space:
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
            DO j = 1,n
               DO k = 1,n
                  IF (  ka2(k,j,i).GT.kprol2 ) B1(k,j,i) = 0.0
               END DO
            END DO
         END DO
#endif
!
! Compute inverse FT of prolongated variable:
         CALL fftp3d_complex_to_real(plancr,B1,br,MPI_COMM_WORLD)
!
! Byte-swap on output:
         IF ( oswap .NE. 0 ) THEN
           CALL rarray_byte_swap(br, n*n*(kend-ksta+1))
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
      DEALLOCATE ( ka)
      DEALLOCATE ( ka2)
      DEALLOCATE ( br )

      CALL MPI_FINALIZE(ierr)

      END PROGRAM BOOTS3D

