!=================================================================
      PROGRAM TRIAD3D
!=================================================================
! TRIAD3D code (part of the GHOST suite)
!
! Numerically computes the triadic transfers for a fixed 
! giving shell q, a band of receiving shells k, and all 
! mediating shells p. A pseudo-spectral method is used to 
! compute spatial derivatives.
! To compile, you need the FFTW library installed on your 
! system. You should link with the FFTP subroutines and use 
! the FFTPLANS and MPIVARS modules (see the file 
! 'fftp_mod.f90').
!
! NOTATION: index 'i' is 'x' 
!           index 'j' is 'y'
!           index 'k' is 'z'
! 
! Conditional compilation options:
!           HD_SOL   builds the hydrodynamic transfer
!
! 2005 Alexandros Alexakis and Pablo D. Mininni.
!      National Center for Atmospheric Research.
!
! 15 Feb 2007: Main program for all transfers
! 21 Feb 2007: POSIX and MPI/IO support
! 10 Mar 2007: FFTW-2.x and FFTW-3.x support
!=================================================================

!
! Modules

      USE fprecision
      USE commtypes
      USE mpivars
      USE filefmt
      USE iovar
      USE grid
      USE fft
      USE ali
      USE var
      USE kes
      IMPLICIT NONE

!
! Arrays for the fields and the transfers

      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: vx,vy,vz,C1
      REAL(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:)    :: Rxx,Rxy,Rxz
      REAL(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:)    :: Ryx,Ryy,Ryz
      REAL(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:)    :: Rzx,Rzy,Rzz
      REAL(KIND=GP), ALLOCATABLE, DIMENSION (:,:)      :: uu

!
! Auxiliary variables

      REAL(KIND=GP)    :: tmp

      INTEGER :: stat
      INTEGER :: cold
      INTEGER :: kini
      INTEGER :: kmin
      INTEGER :: kstp
      INTEGER :: pmax
      INTEGER :: qval
      INTEGER :: i,j,k
      INTEGER :: kk,kp

      TYPE(IOPLAN) :: planio

      CHARACTER(len=100) :: odir,idir

!
! Initializes the MPI and I/O libraries

      CALL MPI_INIT(ierr)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
      CALL range(1,n/2+1,nprocs,myrank,ista,iend)
      CALL range(1,n,nprocs,myrank,ksta,kend)
      CALL io_init(myrank,n,ksta,kend,planio)

!
! Initializes the FFT library
! Use FFTW_ESTIMATE in short runs and FFTW_MEASURE 
! in long runs

      CALL fftp3d_create_plan(planrc,n,FFTW_REAL_TO_COMPLEX, &
          FFTW_MEASURE)
      CALL fftp3d_create_plan(plancr,n,FFTW_COMPLEX_TO_REAL, &
          FFTW_MEASURE)

!
! Allocates memory for distributed blocks

      ALLOCATE( vx(n,n,ista:iend) )
      ALLOCATE( vy(n,n,ista:iend) )
      ALLOCATE( vz(n,n,ista:iend) )
      ALLOCATE( C1(n,n,ista:iend) )
      ALLOCATE( ka(n), ka2(n,n,ista:iend) )
      ALLOCATE( Rxx(n,n,ksta:kend),Rxy(n,n,ksta:kend),Rxz(n,n,ksta:kend) )
      ALLOCATE( Ryx(n,n,ksta:kend),Ryy(n,n,ksta:kend),Ryz(n,n,ksta:kend) )
      ALLOCATE( Rzx(n,n,ksta:kend),Rzy(n,n,ksta:kend),Rzz(n,n,ksta:kend) )

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

!
! Reads from the external file 'triadtran.txt' the 
! parameters that will be used to compute the transfer
!     idir : directory for unformatted input (field components)
!     odir : directory for unformatted output (transfers)
!     stat : number of the file to analyze
!     cold : =0 restart a previous computation
!            =1 start a new computation
!     kmin : minimum value of k
!     kstp : maximum value of k
!     pmax : maximum value of p
!     qval : value of q

      IF (myrank.eq.0) THEN
         OPEN(1,file='triadtran.txt',status='unknown')
         READ(1,'(a100)') idir
         READ(1,'(a100)') odir
         READ(1,*) stat
         READ(1,*) cold
         READ(1,*) kmin
         READ(1,*) kstp
         READ(1,*) pmax
         READ(1,*) qval
         CLOSE(1)
      ENDIF
      CALL MPI_BCAST(idir,100,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(odir,100,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(stat,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(kmin,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(kstp,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(pmax,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(qval,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

!
! Reads the external binary files with the 
! fields components.

      WRITE(ext, fmtext) stat
      CALL io_read(1,idir,'vx',ext,planio,Rxx)
      CALL io_read(1,idir,'vy',ext,planio,Ryy)
      CALL io_read(1,idir,'vz',ext,planio,Rzz)
      CALL fftp3d_real_to_complex(planrc,Rxx,vx,MPI_COMM_WORLD)
      CALL fftp3d_real_to_complex(planrc,Ryy,vy,MPI_COMM_WORLD)
      CALL fftp3d_real_to_complex(planrc,Rzz,vz,MPI_COMM_WORLD)

!
! Allocates the arrays for the transfers

      ALLOCATE( uu(kmin:kstp,pmax) )

!
! If continuing a previous run (cold=0), reads 
! 'kini.txt' and the arrays with the transfers

 RI : IF (myrank.eq.0) THEN

      IF (cold.eq.0) THEN

         OPEN(1,file='kini.txt',status='unknown')
         READ(1,*) kini
         CLOSE(1)
         OPEN(1,file=trim(odir) // '/triadtran_uu.' &
              // ext // '.out' ,form='unformatted')
         READ(1) uu
         CLOSE(1)

      ELSE

         kini = kmin

      ENDIF

      ENDIF RI

      CALL MPI_BCAST(kini,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

!
! Computes the transfer terms
!
! Energy triadic interactions with fixed q:
!     uu(k,p) : -u_k.(u_p.grad u_q)

      CALL dshell(qval,vx,C1,1)
      CALL fftp3d_complex_to_real(plancr,C1,Rxx,MPI_COMM_WORLD)
      CALL dshell(qval,vx,C1,2)
      CALL fftp3d_complex_to_real(plancr,C1,Rxy,MPI_COMM_WORLD)
      CALL dshell(qval,vx,C1,3)
      CALL fftp3d_complex_to_real(plancr,C1,Rxz,MPI_COMM_WORLD)
      CALL dshell(qval,vy,C1,1)
      CALL fftp3d_complex_to_real(plancr,C1,Ryx,MPI_COMM_WORLD)
      CALL dshell(qval,vy,C1,2)
      CALL fftp3d_complex_to_real(plancr,C1,Ryy,MPI_COMM_WORLD)
      CALL dshell(qval,vy,C1,3)
      CALL fftp3d_complex_to_real(plancr,C1,Ryz,MPI_COMM_WORLD)
      CALL dshell(qval,vz,C1,1)
      CALL fftp3d_complex_to_real(plancr,C1,Rzx,MPI_COMM_WORLD)
      CALL dshell(qval,vz,C1,2)
      CALL fftp3d_complex_to_real(plancr,C1,Rzy,MPI_COMM_WORLD)
      CALL dshell(qval,vz,C1,3)
      CALL fftp3d_complex_to_real(plancr,C1,Rzz,MPI_COMM_WORLD)

      DO kk = kini,kstp
         DO kp = 1,pmax

            CALL triagtran(vx,vy,vz,Rxx,Rxy,Rxz,Ryx,Ryy,Ryz, &
                 Rzx,Rzy,Rzz,kk,kp-1,tmp)
            uu(kk,kp) = -tmp

         END DO

! Writes the results each 
! time a row is completed

 RO:     IF (myrank.eq.0) THEN

            OPEN(1,file=trim(odir) // '/triadtran_uu.' &
                 // ext // '.out' ,form='unformatted')
            WRITE(1) uu
            CLOSE(1)

            OPEN(1,file='kini.txt')
            WRITE(1,*) kk+1
            CLOSE(1)

         ENDIF RO

      END DO

!
! End of TRIAD3D

      CALL MPI_FINALIZE(ierr)
      CALL fftp3d_destroy_plan(plancr)
      CALL fftp3d_destroy_plan(planrc)
      DEALLOCATE( uu )
      DEALLOCATE( vx,vy,vz,C1 )
      DEALLOCATE( Rxx,Rxy,Rxz )
      DEALLOCATE( Ryx,Ryy,Ryz )
      DEALLOCATE( Rzx,Rzy,Rzz )

      END PROGRAM TRIAD3D
