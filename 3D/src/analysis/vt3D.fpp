!=================================================================
      PROGRAM VT3D
!=================================================================
! VT3D code (part of the GHOST suite)
!
! Reads velocity binaries and computes all or specified
! components of velocity gradient tensor in Fourier space, and
! then outputs these to disk in real space,
!
!
! 2011 D. Rosenberg
!      NCAR
!
! 17 Sep 2011: Initial version
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
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: dvc


      REAL(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: rv
      REAL(KIND=GP)                                 :: tmp
!
! Auxiliary variables

      INTEGER :: i,ic,iir,ind,ir,iswap,oswap,j,jc,jjc,k,irow(3),jcol(3),nrow,ncol,stat
      INTEGER :: fh

      TYPE(IOPLAN) :: planio
      CHARACTER(len=8)   :: pref
      CHARACTER(len=100) :: odir,idir
      CHARACTER(len=256) :: fout
      CHARACTER(len=100 ):: srow, scol
!
      NAMELIST / vt / idir, odir, srow, scol, iswap, oswap, stat

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
      srow   = '0'
      scol   = '0'
      iswap  = 0
      oswap  = 0
      irow   = 0
      jcol   = 0
!
! Reads from the external file 'vt`.txt' the 
! parameters that will be used to compute the transfer
!     idir   : directory for unformatted input (field components)
!     odir   : directory for unformatted output (prolongated data)
!     stat   : time index for which to compute VT
!     srow   : list of rows to do (';' separated)
!     scol   : list of columns to do (';' separated)
!     iswap  : do endian swap on input?
!     oswap  : do endian swap on output?

      IF (myrank.eq.0) THEN
         OPEN(1,file='vt.txt',status='unknown',form="formatted")
         READ(1,NML=vt)
         CLOSE(1)
      ENDIF
      CALL MPI_BCAST(idir  ,100 ,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(odir  ,100 ,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(srow  ,100 ,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(scol  ,100 ,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(iswap ,1   ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(oswap ,1   ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(stat  ,1   ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)

      pref = 'VT'
!
     
      kmax   = (real(n,kind=GP)/3.)**2

      ALLOCATE( vc(n,n,ista:iend) )
      ALLOCATE( dvc(n,n,ista:iend) )
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

      CALL parseind(srow, ';', irow, 3, nrow) 
      CALL parseind(scol, ';', jcol, 3, ncol) 

      WRITE(ext, fmtext) stat
      tmp = 1./REAL(n,KIND=GP)**3
      DO ir = 1, nrow
        iir = irow(ir)
! read in appropriate file:
        IF ( iir.EQ.1 ) THEN
          CALL io_read(1,idir,'vx',ext,planio,rv)
        ELSE IF ( iir.EQ.2 ) THEN
          CALL io_read(1,idir,'vy',ext,planio,rv)
        ELSE IF ( iir.EQ.3 ) THEN
          CALL io_read(1,idir,'vz',ext,planio,rv)
        ELSE
          WRITE(*,*) 'main: invalid row specified: ',i
          STOP 
        ENDIF
!
! Byte-swap on input:
        IF ( iswap .NE. 0 ) THEN
          CALL rarray_byte_swap(rv, n*n*(kend-ksta+1))
        ENDIF
!
! take FFT of component:
        CALL fftp3d_real_to_complex(planrc,rv,vc,MPI_COMM_WORLD)
        DO jc = 1, ncol
          jjc = jcol(jc)
          IF ( jjc.LT.1 .OR. jjc.GT.3 ) THEN
            WRITE(*,*) 'main: invalid column specified: ',j
            STOP 
          ENDIF
          CALL derivk3(vc, dvc, jjc)
          DO i = ista,iend
            DO j = 1,n
              DO k = 1,n
                dvc(k,j,i) = dvc(k,j,i)*tmp
              END DO
            END DO
          END DO
          CALL fftp3d_complex_to_real(plancr,dvc,rv,MPI_COMM_WORLD)
!
! Byte-swap on output:
          IF ( oswap .NE. 0 ) THEN
            CALL rarray_byte_swap(rv, n*n*(kend-ksta+1))
          ENDIF
          WRITE(fout,'(a2,i01,i01)') trim(pref),iir,jjc
          fout = trim(odir) // '/' // trim(fout) // '.' // ext // '.out' 
          bmangle = 0
          CALL io_write(1,odir,fout,0,planio,rv)
          bmangle = 1
          WRITE(*,*)'main: ', trim(fout), ' written.'
        ENDDO
      ENDDO
!
      CALL fftp3d_destroy_plan(plancr)
      CALL fftp3d_destroy_plan(planrc)


      DEALLOCATE ( vc,dvc)
      DEALLOCATE ( rv)
      DEALLOCATE ( ka)
      DEALLOCATE ( ka2)

      CALL MPI_FINALIZE(ierr)

      END PROGRAM VT3D

