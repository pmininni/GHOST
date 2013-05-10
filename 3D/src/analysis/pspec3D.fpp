!=================================================================
      PROGRAM PSPEC3D
!=================================================================
! PSPEC3D code (part of the GHOST suite)
!
! Reads real quantity, computes power spectrum, output to a file
!
!
! 2013 D. Rosenberg
!      NCAR
!
! 1o May 2013: Initial version
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

      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: vc


      REAL(KIND=GP)   , ALLOCATABLE, DIMENSION (:,:,:) :: rv
      REAL(KIND=GP)                                    :: tmp
!
! Auxiliary variables

      INTEGER :: i,ic,iir,ind,ir,iswap,it,maxfn,oswap,j,jc,jjc,k
      INTEGER :: nfiles

      TYPE(IOPLAN) :: planio
      CHARACTER(len=256)  :: odir,idir
      CHARACTER(len=256)  :: fnout
      CHARACTER(len=256)  :: fnin(256)
      CHARACTER(len=1024) :: fntmp
      CHARACTER(len=65537):: fnstr
!
      NAMELIST / pspec / idir, odir, fnstr, iswap, oswap

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
      iswap  = 0
      oswap  = 0
      maxfn  = 256;
!
! Reads from the external file 'vt`.txt' the 
! parameters that will be used to compute the transfer
!     idir   : directory for unformatted input (field components)
!     odir   : directory for unformatted output (prolongated data)
!     fnin   :  a ';'--separated list of file names to compute spectra of
!     iswap  : do endian swap on input?
!     oswap  : do endian swap on output?

      IF (myrank.eq.0) THEN
         OPEN(1,file='pspec.txt',status='unknown',form="formatted")
         READ(1,NML=pspec)
         CLOSE(1)
      ENDIF
      CALL MPI_BCAST(idir  ,256  ,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(odir  ,256  ,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(fnout ,256  ,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(fnstr ,65537,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(iswap ,1    ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(oswap ,1    ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)

!
      kmax   = (real(n,kind=GP)/3.)**2

      ALLOCATE( vc(n,n,ista:iend) )
      ALLOCATE( ka(n),ka2(n,n,ista:iend) )
      ALLOCATE( rv(n,n,ksta:kend) )
!

      CALL fftp3d_create_plan(planrc,n,FFTW_REAL_TO_COMPLEX, &
          FFTW_MEASURE)
      CALL fftp3d_create_plan(plancr,n,FFTW_COMPLEX_TO_REAL, &
          FFTW_MEASURE)
!
! Some constants for the FFT
!     kmax: maximum truncation for dealiasing
!     tiny: minimum truncation for dealiasing

      kmax = (REAL(n,KIND=GP)/3.)**2
      tiny = 1e-5

!
! Builds the wave number and the square wave 
! number matrixes

      DO i = 1,n/2
         ka(i) = REAL(i-1,KIND=GP)
         ka(i+n/2) = REAL(i-n/2-1,KIND=GP)
      END DO
      DO i = ista,iend
         DO j = 1,n
            DO k = 1,n
               ka2(k,j,i) = ka(i)**2+ka(j)**2+ka(k)**2
            END DO
         END DO
      END DO

      CALL parsestr(fnstr,';', fnin, maxfn, 256, nfiles) 

      tmp = 1./REAL(n,KIND=GP)**3
      DO it = 1,nfiles
            bmangle = 0
            CALL io_read(1,idir,fnin(it),1,planio,rv)
            bmangle = 1
! Byte-swap on input:
          IF ( iswap .NE. 0 ) THEN
            CALL rarray_byte_swap(rv, n*n*(kend-ksta+1))
          ENDIF
!
! Take FFT of component:
          CALL fftp3d_real_to_complex(planrc,rv,vc,MPI_COMM_WORLD)
!
! Compute power spectru and output it:
          fnout = 'sp' // trim(fnin(it)) // '.txt';
          fntmp = trim(odir) // '/' // trim(fnout)
          CALL pspectrum(vc,fntmp)
  
          IF ( myrank.EQ.0 ) THEN
            WRITE(*,*)'main: ', trim(fntmp), ' written.'
          ENDIF
      ENDDO     ! time index loop
!
      CALL fftp3d_destroy_plan(plancr)
      CALL fftp3d_destroy_plan(planrc)


      DEALLOCATE ( vc)
      DEALLOCATE ( rv)
      DEALLOCATE ( ka)
      DEALLOCATE ( ka2)

      CALL MPI_FINALIZE(ierr)

      END PROGRAM PSPEC3D


!*****************************************************************
      SUBROUTINE pspectrum(a,fnout)
!-----------------------------------------------------------------
!
! Computes the 'power' spectrum of specified 
! scalar quantity. The output is written to a 
! file by the first node.
!
! Parameters
!     a    : input matrix 
!     fnout: output file name 
!
      USE kes
      USE grid
      USE mpivars
      USE filefmt
      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(n/2+1)                     :: Ek
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,n,ista:iend) :: a
      CHARACTER(len=*), INTENT(IN)                           :: fnout
!
! Computes the energy and/or helicity spectra
      CALL pspectrumc(a,Ek)
!
! Exports the energy spectrum to a file
!
      IF (myrank.eq.0) THEN
         OPEN(1,file=trim(fnout))
         WRITE(1,20) Ek
   20    FORMAT( E23.15 ) 
         CLOSE(1)
      ENDIF
!
      RETURN
      END SUBROUTINE pspectrum


!*****************************************************************
      SUBROUTINE pspectrumc(a,Ektot)
!-----------------------------------------------------------------
!
! Computes the 'power' spectrum of a, returning it
!
! Parameters
!     a    : input matrix in the x-direction
!     Ektot: output power spectrum
!
      USE fprecision
      USE commtypes
      USE kes
      USE grid
      USE mpivars
!$    USE threads
      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(n/2+1) :: Ek
      DOUBLE PRECISION, INTENT(OUT), DIMENSION(n/2+1) :: Ektot
      DOUBLE PRECISION    :: tmq
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,n,ista:iend) :: a
      COMPLEX(KIND=GP), DIMENSION(n,n,ista:iend)             :: c1
      REAL(KIND=GP)       :: tmp
      INTEGER             :: i,j,k
      INTEGER             :: kmn

!
! Computes the kinetic energy spectrum
!
      tmp = 1.0_GP/real(n,kind=GP)**6
       DO i = 1,n/2+1
          Ek   (i) = 0.0D0
          Ektot(i) = 0.0D0
       END DO
       IF (ista.eq.1) THEN
!$omp parallel do private (k,kmn,tmq)
          DO j = 1,n
             DO k = 1,n
                kmn = int(sqrt(ka2(k,j,1))+.501)
                IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                   tmq = (abs(a(k,j,1))**2 )*tmp
!$omp atomic
                   Ek(kmn) = Ek(kmn)+tmq
                ENDIF
             END DO
          END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k,kmn,tmq)
          DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k,kmn,tmq)
             DO j = 1,n
                DO k = 1,n
                   kmn = int(sqrt(ka2(k,j,i))+.501)
                   IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                      tmq = 2*(abs(a(k,j,i))**2 )*tmp
!$omp atomic
                      Ek(kmn) = Ek(kmn)+tmq
                   ENDIF
                END DO
             END DO
          END DO
        ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k,kmn,tmq)
          DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,kmn,tmq)
             DO j = 1,n
                DO k = 1,n
                   kmn = int(sqrt(ka2(k,j,i))+.501)
                   IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                      tmq = 2*(abs(a(k,j,i))**2)*tmp
!$omp atomic
                      Ek(kmn) = Ek(kmn)+tmq
                   ENDIF
                END DO
             END DO
          END DO
       ENDIF
!
! Computes the reduction between nodes
!
       CALL MPI_ALLREDUCE(Ek,Ektot,n/2+1,MPI_DOUBLE_PRECISION,      &
                      MPI_SUM,MPI_COMM_WORLD,ierr)
!
       RETURN
       END SUBROUTINE pspectrumc
