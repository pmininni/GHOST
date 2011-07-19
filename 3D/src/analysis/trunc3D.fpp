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
! When building either this of the 'boots3d' utillity, the Makefile.in
! should be modified with the largest dataset size. So, it would be
! the original data size for trunc3D.
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

      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: C1
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: T1


      REAL(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: vv
      REAL(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: tr
      REAL(KIND=GP)                                 :: fact


!
! Auxiliary variables

      REAL(KIND=GP)    :: ktrunc2, ktrunc
      INTEGER :: fftdir, flags, nfiles,nmt,np,npt,nt,ntprocs, npkeep
      INTEGER :: i,ib,ie,ind,iswap,oswap,itsta,itend,j,k,ktsta,ktend
      INTEGER :: istak,iendk,kstak,kendk
      INTEGER :: mykrank,mytrank

      INTEGER :: commtrunc, fh, groupworld, grouptrunc, iExclude(3,1), iInclude(3,1)

      TYPE(IOPLAN)  :: planio, planiot
      TYPE(FFTPLAN) :: plancrt

      CHARACTER(len=8)   :: suff
      CHARACTER(len=100) :: odir,idir
      CHARACTER(len=256) :: fname, fout, msg
      CHARACTER(len=1024):: fnlist

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
! Reads from the external file 'trunc.txt' the 
! parameters that will be used to compute the transfer
!     idir   : directory for unformatted input (field components)
!     odir   : directory for unformatted output (truncated data)
!     fnlist : file list to truncate, separated by ';'
!     iswap  : do endian swap on input?
!     oswap  : do endian swap on output?
!     nt     : truncation wavenumber

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
      CALL MPI_BCAST(nt    ,1   ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)

      IF ( nt .GT. n .OR. nt .LT. 1 ) THEN
        WRITE(*,*)'MAIN: truncation specification incorrect: Nmax=', n, ' nt=', nt
        STOP 
      ENDIF

      WRITE(suff,'(a2,i5.5)') '_T', nt

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
        CALL MPI_GROUP_RANGE_EXCL(groupworld, 1, iExclude, grouptrunc, ierr)   
        CALL MPI_COMM_CREATE(MPI_COMM_WORLD, grouptrunc, commtrunc, ierr)
      ELSE IF ( ntprocs .EQ. nprocs ) THEN
        CALL MPI_COMM_DUP(MPI_COMM_WORLD,commtrunc,ierr)
        CALL MPI_COMM_GROUP(MPI_COMM_WORLD,grouptrunc,ierr)
      ENDIF
      CALL trrange(1,n    ,nt    ,nprocs,myrank,ktsta,ktend)
      CALL trrange(1,n/2+1,nt/2+1,nprocs,myrank,itsta,itend)
      ktrunc  = real(nt,kind=GP)/2.0 - 1.0
      ktrunc2 = (real(ktrunc,kind=GP) )**2 

      ALLOCATE( vv(n,n,ksta:kend) )
      ALLOCATE( C1(n,n,ista:iend) )
      ALLOCATE( T1(nt,nt,itsta:itend) )
      ALLOCATE( ka(n),ka2(n,n,ista:iend) )
      ALLOCATE( tr(nt,nt,ktsta:ktend) )

      IF ( myrank .LT. ntprocs ) THEN
        npkeep = nprocs; nprocs = ntprocs
        CALL io_init(myrank,nt,ktsta,ktend,planiot)
        nprocs = npkeep
      ENDIF
!
      flags  = FFTW_MEASURE
!
      kmax = (real(n,kind=GP)/3.)**2
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
!
! Create regular plan, and one for the truncated inverse transform:
      CALL fftp3d_create_plan(planrc,n,FFTW_REAL_TO_COMPLEX, flags)
      CALL trrange(1,n    ,nt    ,nprocs,myrank,ksta,kend)
      CALL trrange(1,n/2+1,nt/2+1,nprocs,myrank,ista,iend)
      IF ( myrank.eq.0 ) write(*,*)'main: creating trplan...'
      CALL fftp3d_create_trplan(plancrt,n,nt,FFTW_COMPLEX_TO_REAL,flags)
      CALL range(1,n/2+1,nprocs,myrank,ista,iend)
      CALL range(1,n,nprocs,myrank,ksta,kend)

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
           CALL rarray_byte_swap(vv, n*n*(kend-ksta+1))
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
            DO j = 1,n
               DO k = 1,n
                  IF (  ka2(k,j,i).GT.ktrunc2 ) C1(k,j,i) = 0.0
               END DO
            END DO
         END DO
!
         fact = 1.0_GP / real(n,kind=GP)**3
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
            DO j = 1,n
               DO k = 1,n
                  C1(k,j,i) = C1(k,j,i)*fact
               END DO
            END DO
         END DO
!
! Store in coeffs in reduced storage for inverse transform:
         T1 = 0.0
         DO i = itsta,itend
            DO j = 1,nt/2
               DO k = 1,nt/2
                  T1(k,j,i) = C1(k,j,i)
               END DO
               DO k = nt/2+1,nt
                  T1(k,j,i) = C1(k-nt+n,j,i)
               END DO
            END DO
            DO j = nt/2+1,nt
               DO k = 1,nt/2
                  T1(k,j,i) = C1(k,j-nt+n,i)
               END DO
               DO k = nt/2+1,nt
                  T1(k,j,i) = C1(k-nt+n,j-nt+n,i)
               END DO
            END DO
         END DO
!
! Compute inverse FT of truncated variable:
!
         CALL trrange(1,n    ,nt    ,nprocs,myrank,ksta,kend)
         CALL trrange(1,n/2+1,nt/2+1,nprocs,myrank,ista,iend)
         IF ( myrank.eq.0 ) write(*,*)'main: complex_to_real...'
         CALL fftp3d_complex_to_real(plancrt,T1,tr,MPI_COMM_WORLD)
         IF ( myrank.eq.0 ) write(*,*)'main: complex_to_real done.'
         CALL range(1,n/2+1,nprocs,myrank,ista,iend)
         CALL range(1,n,nprocs,myrank,ksta,kend)
!
! Put to disk:
         fout = trim(odir) // '/' // trim(fname) // trim(suff)
         IF ( myrank.eq.0 ) write(*,*)'main: writing to disk...'
         IF ( myrank .LT. ntprocs ) THEN
           IF ( iswap .NE. 0 ) THEN
             CALL rarray_byte_swap(tr, nt*nt*(ktend-ktsta+1))
           ENDIF
           CALL MPI_FILE_OPEN(commtrunc, fout, &
             MPI_MODE_CREATE+MPI_MODE_WRONLY, &
             MPI_INFO_NULL,fh,ioerr)
           CALL MPI_FILE_SET_VIEW(fh,disp,GC_REAL,planiot%iotype,'native', &
             MPI_INFO_NULL,ioerr)
           CALL MPI_FILE_WRITE_ALL(fh,tr, &
             planiot%n*planiot%n*(planiot%kend-planiot%ksta+1),GC_REAL, &
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
      DEALLOCATE ( ka)
      DEALLOCATE ( ka2 )
      DEALLOCATE ( C1 )
      DEALLOCATE ( T1 )
      DEALLOCATE ( tr )

      CALL MPI_FINALIZE(ierr)

      END PROGRAM TRUNC3D

