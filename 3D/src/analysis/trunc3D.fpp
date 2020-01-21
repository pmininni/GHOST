#include "fftw_wrappers.h"

!=================================================================
      PROGRAM TRUNC3D
!=================================================================
! TRUNC3D code (part of the GHOST suite)
!
! Reads dataset, and computes dataset that is truncated in 
! wavenumber space, and outputs to real space. The files
! to be read are specified on input, as is the truncation size;
! the output files are the originial files appended with with the 
! truncation size.
!
! When building either this or the 'boots3d' utillity, the Makefile.in
! should be modified with the largest dataset size. So, it would be
! the original data size for trunc3D.
!
! This tool is not guaranteed to work at all resolutions, or to
! truncate to some specific sizes. It also may not work depending
! on the number of processors used. See 'boots3d' for more details.
! Note that unlike 'boots3d', this tool does not check for a
! correct specification of the number of processors.
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
! Arrays for the fields

      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: C1
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: T1


      REAL(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: vv
      REAL(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: tr
      REAL(KIND=GP)                                 :: fact, rmp, rmq, rms

!
! Auxiliary variables

      INTEGER :: fftdir,flags,nfiles,nmt,np,npt,nxt,nyt,nzt,ntprocs,npkeep
      INTEGER :: i,ib,ie,ind,itsta,itend,j,k,ktsta,ktend
      INTEGER :: istak,iendk,kstak,kendk
      INTEGER :: mykrank,mytrank

      INTEGER :: commtrunc, fh, groupworld, grouptrunc, iExclude(3,1), iInclude(3,1)

      TYPE(IOPLAN)  :: planio, planiot
      TYPE(FFTPLAN) :: plancrt

      CHARACTER(len=8)   :: suff
      CHARACTER(len=100) :: odir,idir
      CHARACTER(len=256) :: fname, fout, msg
      CHARACTER(len=1024):: fnlist

      NAMELIST / regrid / idir, odir, fnlist, iswap, oswap, nxt, nyt, nzt

!
! Initializes the MPI and I/O libraries

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
! Reads from the external file 'trunc.txt' the 
! parameters that will be used to compute the transfer
!     idir   : directory for unformatted input (field components)
!     odir   : directory for unformatted output (truncated data)
!     fnlist : file list to truncate, separated by ';'
!     iswap  : do endian swap on input?
!     oswap  : do endian swap on output?
!     nxt    : truncated linear size in x direction
!     nyt    : truncated linear size in y direction
!     nzt    : truncated linear size in z direction

      IF (myrank.eq.0) THEN
         OPEN(1,file='trunc.txt',status='unknown')
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

! Check input
      IF ( nxt .GT. nx .OR. nxt .LT. 1 ) THEN
        IF ( myrank .eq. 0) &
	  PRINT*, 'MAIN: truncation specification incorrect; input nxt must be less than Nx'
        CALL MPI_Finalize(ierr)
        STOP
      ENDIF
      IF ( nyt .GT. ny .OR. nyt .LT. 1 ) THEN
        IF ( myrank .eq. 0) &
	  PRINT*, 'MAIN: truncation specification incorrect; input nyt must be less than Ny'
        CALL MPI_Finalize(ierr)
        STOP
      ENDIF
      IF ( nzt .GT. nz .OR. nzt .LT. 1 ) THEN
        IF ( myrank .eq. 0) &
	  PRINT*, 'MAIN: truncation specification incorrect; input nzt must be less than Nz'
        CALL MPI_Finalize(ierr)
        STOP
      ENDIF

! Suffix for the output
      WRITE(suff,'(a2,i5.5,a1,i5.5,a1,i5.5)') '_T', nxt, "-", nyt, "-", nzt

! Create new communicator for the truncated data:
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
        CALL MPI_GROUP_RANGE_EXCL(groupworld, 1, iExclude, grouptrunc, ierr)   
        CALL MPI_COMM_CREATE(MPI_COMM_WORLD, grouptrunc, commtrunc, ierr)
      ELSE IF ( ntprocs .EQ. nprocs ) THEN
        CALL MPI_COMM_DUP(MPI_COMM_WORLD,commtrunc,ierr)
        CALL MPI_COMM_GROUP(MPI_COMM_WORLD,grouptrunc,ierr)
      ENDIF
      CALL trrange(1,nz    ,nzt    ,nprocs,myrank,ktsta,ktend)
      CALL trrange(1,nx/2+1,nxt/2+1,nprocs,myrank,itsta,itend)

      ALLOCATE( vv(nx, ny, ksta:kend) )
      ALLOCATE( C1(nz, ny, ista:iend) )
      ALLOCATE( T1(nzt,nyt,itsta:itend) )
      ALLOCATE( kx(nx), ky(ny), kz(nz) )
      ALLOCATE( kn2(nz,ny,ista:iend) )
      ALLOCATE( kk2(nz,ny,ista:iend) )
      ALLOCATE( tr(nxt,nyt,ktsta:ktend) )

      IF ( myrank .LT. ntprocs ) THEN
        npkeep = nprocs; nprocs = ntprocs
        CALL io_init(myrank,(/nxt,nyt,nzt/),ktsta,ktend,planiot)
        nprocs = npkeep
      ENDIF
!
      flags  = FFTW_MEASURE
!
! Populate wavenumber-associated arrays
      kmax =     1.0_GP/9.0_GP
      
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
            DO k = 1,nz ! squared wavenumber normalized in the old grid
               kk2(k,j,i) = rmp*kx(i)**2+rmq*ky(j)**2+rms*kz(k)**2
            END DO
         END DO
      END DO
      
      rmp = 1.0_GP/real(nxt,kind=GP)**2
      rmq = 1.0_GP/real(nyt,kind=GP)**2
      rms = 1.0_GP/real(nzt,kind=GP)**2
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
            DO k = 1,nz ! squared wavenumber normalized in the new grid
               kn2(k,j,i) = rmp*kx(i)**2+rmq*ky(j)**2+rms*kz(k)**2
            END DO
         END DO
      END DO

!
! Create regular plan, and one for the truncated inverse transform:
      CALL fftp3d_create_plan(planrc,(/nx,ny,nz/),FFTW_REAL_TO_COMPLEX, flags)
      CALL trrange(1,nz    ,nzt    ,nprocs,myrank,ksta,kend)
      CALL trrange(1,nz/2+1,nzt/2+1,nprocs,myrank,ista,iend)
      IF ( myrank.eq.0 ) write(*,*)'main: creating trplan...'
      CALL fftp3d_create_trplan(plancrt,(/nx,ny,nz/),(/nxt,nyt,nzt/),FFTW_COMPLEX_TO_REAL,flags)
      CALL range(1,nx/2+1,nprocs,myrank,ista,iend)
      CALL range(1,nz    ,nprocs,myrank,ksta,kend)

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
         IF ( myrank .EQ. 0 ) THEN
           WRITE(*,*) 'main: Truncating file: '//trim(fname)//'...'
         ENDIF
!
! Read the external binary files in real space:
! 
         bmangle = 0
         CALL io_read(1,idir,fname,1,planio,vv)
         bmangle = 1

         IF ( iswap .NE. 0 ) THEN
           CALL rarray_byte_swap(vv, nx*ny*(kend-ksta+1))
         ENDIF

!
! Compute FT of variable:
         IF ( myrank.eq.0 ) write(*,*)'main: real_to_complex...'
         CALL fftp3d_real_to_complex(planrc,vv,C1,MPI_COMM_WORLD)
!
! Truncate in Fourier space:
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
            DO j = 1,ny
               DO k = 1,nz
                  IF (  kn2(k,j,i).GT.kmax ) C1(k,j,i) = 0.0
               END DO
            END DO
         END DO
!
         fact = 1.0_GP / (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
            DO j = 1,ny
               DO k = 1,nz
                  C1(k,j,i) = C1(k,j,i)*fact
               END DO
            END DO
         END DO
!
! Store in coeffs in reduced storage for inverse transform:
         T1 = 0.0
         DO i = itsta,itend
            DO j = 1,nyt/2
               DO k = 1,nzt/2
                  T1(k,j,i) = C1(k,j,i)
               END DO
               DO k = nzt/2+1,nzt
                  T1(k,j,i) = C1(k-nzt+nz,j,i)
               END DO
            END DO
            DO j = nyt/2+1,nyt
               DO k = 1,nzt/2
                  T1(k,j,i) = C1(k,j-nyt+ny,i)
               END DO
               DO k = nzt/2+1,nzt
                  T1(k,j,i) = C1(k-nzt+nz,j-nyt+ny,i)
               END DO
            END DO
         END DO
!
! Compute inverse FT of truncated variable:
!
         CALL trrange(1,nz    ,nzt    ,nprocs,myrank,ksta,kend)
         CALL trrange(1,nx/2+1,nxt/2+1,nprocs,myrank,ista,iend)
         IF ( myrank.eq.0 ) write(*,*)'main: complex_to_real...'
         CALL fftp3d_complex_to_real(plancrt,T1,tr,MPI_COMM_WORLD)
         IF ( myrank.eq.0 ) write(*,*)'main: complex_to_real done.'
         CALL range(1,nx/2+1,nprocs,myrank,ista,iend)
         CALL range(1,nz,    nprocs,myrank,ksta,kend)
!
! Put to disk:
         fout = trim(odir) // '/' // trim(fname) // trim(suff)
         IF ( myrank.eq.0 ) write(*,*)'main: writing to disk...'
         IF ( myrank .LT. ntprocs ) THEN
           IF ( iswap .NE. 0 ) THEN
             CALL rarray_byte_swap(tr, nxt*nyt*(ktend-ktsta+1))
           ENDIF
           CALL MPI_FILE_OPEN(commtrunc, fout, &
             MPI_MODE_CREATE+MPI_MODE_WRONLY, &
             MPI_INFO_NULL,fh,ioerr)
           CALL MPI_FILE_SET_VIEW(fh,disp,GC_REAL,planiot%iotype,'native', &
             MPI_INFO_NULL,ioerr)
           CALL MPI_FILE_WRITE_ALL(fh,tr, &
             planiot%nx*planiot%ny*(planiot%kend-planiot%ksta+1),GC_REAL,  &
             MPI_STATUS_IGNORE,ioerr)
           CALL MPI_FILE_CLOSE(fh,ioerr)
           IF ( myrank .EQ. 0 ) THEN
             WRITE(*,*) 'main: ',trim(fout),' written.'
           ENDIF
         ENDIF
!
      ENDDO
      IF ( myrank.eq.0 ) write(*,*)'main: cleaning up...'
      CALL fftp3d_destroy_plan(planrc)
      IF ( myrank .LT. ntprocs ) THEN
        CALL fftp3d_destroy_plan(plancrt)
      ENDIF
!
!
      IF ( commtrunc .NE. MPI_COMM_NULL ) THEN
        CALL MPI_COMM_FREE(commtrunc, ierr)
      endif
      IF ( grouptrunc .NE. MPI_GROUP_NULL ) THEN
        CALL MPI_GROUP_FREE(grouptrunc, ierr)
      ENDIF


      DEALLOCATE ( vv )
      DEALLOCATE ( kx,ky,kz )
      DEALLOCATE ( kn2,kk2  )
      DEALLOCATE ( C1 )
      DEALLOCATE ( T1 )
      DEALLOCATE ( tr )

      CALL MPI_FINALIZE(ierr)

      END PROGRAM TRUNC3D

