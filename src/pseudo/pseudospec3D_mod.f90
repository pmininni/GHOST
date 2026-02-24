!=================================================================
! General MODULES for pseudospectral methods
!
! 2003 Pablo D. Mininni.
!      Department of Physics, 
!      Facultad de Ciencias Exactas y Naturales.
!      Universidad de Buenos Aires.
!      e-mail: mininni@df.uba.ar 
!=================================================================

!=================================================================

  MODULE filefmt
!
! Change the length of the string 'ext' to change the number 
! of characters used to number the binary files and files with 
! the spectra. The format fmtext should be consistent with the 
! length of the string, e.g. if len=5 then fmtext = '(i5.5)'.
      CHARACTER(len=4)      :: ext
      CHARACTER(len=6),SAVE :: fmtext = '(i4.4)'

  END MODULE filefmt
!=================================================================

  MODULE grid
!
! ni: number of points in the spatial grid in each direction
      INTEGER :: nx
      INTEGER :: ny
      INTEGER :: nz
    CONTAINS
!-----------------------------------------------------------------
      SUBROUTINE grid_init(infile_)
!
! Initializes the grid spatial resolution.
      USE commtypes
      USE mpivars
      IMPLICIT NONE

      CHARACTER(len=*), INTENT(IN) :: infile_
      NAMELIST/ grid / nx,ny,nz
      IF ( myrank .EQ. 0 ) THEN
         OPEN(1,file=infile_,status='unknown',form="formatted")
         READ(1,NML=grid)
         CLOSE(1)
      ENDIF
      CALL MPI_BCAST(nx,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(ny,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(nz,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      END SUBROUTINE grid_init

  END MODULE grid
!=================================================================

  MODULE kes
      USE fprecision
      REAL(KIND=GP), ALLOCATABLE, DIMENSION (:)        :: kx,ky,kz
      REAL(KIND=GP), TARGET, ALLOCATABLE, DIMENSION (:,:,:) :: kn2
      REAL(KIND=GP), POINTER, DIMENSION (:,:,:)             :: kk2
      INTEGER                                     :: nmax,nmaxperp
      SAVE

  END MODULE kes
!=================================================================

  MODULE ali
      USE fprecision
      REAL(KIND=GP) :: kmax
      REAL(KIND=GP) :: tiny
      REAL(KIND=GP) :: tinyf
      SAVE

  END MODULE ali
!=================================================================

  MODULE boxsize
!
! Lx, Ly, Lz, and Dk are the lengths of the sides of the box,
! and the width of Fourier shells.
!  Lx  : Length in x (in units of 2.pi, =1 gives a side of length 2.pi)
!  Ly  : Length in y
!  Lz  : Length in z
!  Dkk : Width of Fourier shells for 2D and 3D spectra
!        Default = min(1/Lx, 1/Ly, 1/Lz)
      USE fprecision
      REAL(KIND=GP)    :: Lx,Ly,Lz
      REAL(KIND=GP)    :: Dkx,Dky,Dkz,Dkk
      LOGICAL          :: anis   ! .TRUE. if in an anisotropic box
    CONTAINS
!-----------------------------------------------------------------
      SUBROUTINE box_init(infile_)
!
! Initializes the box size, and arrays with wavenumbers. If no
! parameters are present, the code attempts to initialize a cubic
! box with Lx=Ly=Lz=2.pi
      USE commtypes
      USE mpivars
      USE grid
      USE kes
      USE ali
!$    USE threads
      IMPLICIT NONE

      CHARACTER  (len=*), INTENT(IN) :: infile_
      CHARACTER(len=256)             :: iomsg
      REAL(KIND=GP)                  :: rmp,rmq,rms
      INTEGER                        :: i,j,k,ios
      NAMELIST / boxparams / Lx,Ly,Lz,Dkk

      Lx  = 1.0_GP; Ly  = 1.0_GP; Lz  = 1.0_GP
      Dkx = 1.0_GP; Dky = 1.0_GP; Dkz = 1.0_GP
      Dkk = 0.0_GP; anis= .FALSE.

      IF (myrank.eq.0) THEN
         OPEN(1,file=infile_,status='unknown',form="formatted")
         READ(1,NML=boxparams,iostat=ios,iomsg=iomsg)
	 IF (ios .eq. 0) THEN
	    anis = .TRUE.      ! Domain size defined in infile_
	 ELSE
	    PRINT *,'boxparams not found, attempting to use a 2.pi cubic box'
            IF ((nx.ne.ny).or.(ny.ne.nz)) anis = .TRUE.
	 ENDIF
         CLOSE(1)
      ENDIF
      CALL MPI_BCAST(Lx  ,1,GC_REAL    ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(Ly  ,1,GC_REAL    ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(Lz  ,1,GC_REAL    ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(Dkk ,1,GC_REAL    ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(anis,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      IF ( anis ) THEN
         Dkx = 1.0_GP/Lx
         Dky = 1.0_GP/Ly
         Dkz = 1.0_GP/Lz
      ENDIF
      IF (Dkk.lt.1e-5) Dkk = min(Dkx,Dky,Dkz)
      ! Allocates arrays with wavenumbers
      ALLOCATE( kx(nx), ky(ny), kz(nz) )
      ALLOCATE( kn2(nz,ny,ista:iend) )
      IF ( anis ) THEN
         ALLOCATE( kk2(nz,ny,ista:iend) )
      ELSE
         kk2 => kn2
      ENDIF
      ! Sets constants for the pseudospectral method
      !   kmax: maximum truncation for dealiasing
      !   tiny: minimum truncation for dealiasing
      kmax =     1.0_GP/9.0_GP
      nmax =     int(max(nx*Dkx,ny*Dky,nz*Dkz)/Dkk)
      nmaxperp = int(max(nx*Dkx,ny*Dky)/Dkk)
      IF ( .NOT. anis ) kmax = kmax*real(nx,kind=GP)**2
      tiny  = min(1e-5_GP ,.1_GP/(real(nmax,kind=GP)**2))
      tinyf = min(1e-15_GP,.1_GP/(real(nmax,kind=GP)**2))
      ! Initializez arrays with the wavenumbers and the
      ! square wavenumbers. At the end, kx, ky, and kz
      ! have wavenumbers with dimensions, kk2 has the
      ! squared wavenumbers with dimensions, and kn2 has
      ! the dimensionless and normalized squared
      ! wavenumbers used for dealiasing.
      DO i = 1,nx/2
         kx(i) = real(i-1,kind=GP)
         kx(i+nx/2) = real(i-nx/2-1,kind=GP)
      END DO
      IF (nx.eq.1) kx(1) = 0.0_GP ! Case with no extension in x
      DO j = 1,ny/2
         ky(j) = real(j-1,kind=GP)
         ky(j+ny/2) = real(j-ny/2-1,kind=GP)
      END DO
      IF (ny.eq.1) ky(1) = 0.0_GP ! Case with no extension in y
      DO k = 1,nz/2
         kz(k) = real(k-1,kind=GP)
         kz(k+nz/2) = real(k-nz/2-1,kind=GP)
      END DO
      IF (nz.eq.1) kz(1) = 0.0_GP ! Case with no extension in z
      IF ( anis ) THEN
         rmp = 1.0_GP/real(nx,kind=GP)**2
         rmq = 1.0_GP/real(ny,kind=GP)**2
         rms = 1.0_GP/real(nz,kind=GP)**2
      ELSE
         rmp = 1.0_GP
	 rmq = 1.0_GP
	 rms = 1.0_GP
      ENDIF
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
            DO k = 1,nz
               kn2(k,j,i) = rmp*kx(i)**2+rmq*ky(j)**2+rms*kz(k)**2
            END DO
         END DO
      END DO
      IF ( anis ) THEN
         kx = kx*Dkx
         ky = ky*Dky
         kz = kz*Dkz
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
            DO j = 1,ny
               DO k = 1,nz
                  kk2(k,j,i) = kx(i)**2+ky(j)**2+kz(k)**2
               END DO
            END DO
         END DO
      ENDIF
      END SUBROUTINE box_init

  END MODULE boxsize
!=================================================================

  MODULE fft
      USE fftplans
      TYPE(FFTPLAN) :: planrc, plancr
      SAVE

  END MODULE fft
!=================================================================

  MODULE var
      USE fprecision
      REAL(KIND=GP)    :: pi = 3.14159265358979323846_GP
      COMPLEX(KIND=GP) :: im = (0.0_GP,1.0_GP)
      SAVE

  END MODULE var
!=================================================================

  MODULE random
      USE var
      USE fprecision
    CONTAINS
!-----------------------------------------------------------------
      REAL(KIND=GP) FUNCTION randu(idum)
!
! Uniform distributed random numbers between -1 and 
! 1. The seed idum must be between 0 and the value 
! of mask
      INTEGER, PARAMETER :: iq=127773,ir=2836,mask=123459876
      INTEGER, PARAMETER :: ia=16807,im=2147483647
      INTEGER            :: k,idum
      
      REAL(KIND=GP), PARAMETER :: am=1./im
      idum = ieor(idum,mask)
      k = idum/iq
      idum = ia*(idum-k*iq)-ir*k
      IF (idum.lt.0) idum = idum+im
      randu = am*idum
      randu = (randu-.5)*2
      idum = ieor(idum,mask)
      END FUNCTION randu

!-----------------------------------------------------------------
      REAL(KIND=GP) FUNCTION randn(idum)
!
! Normally distributed random numbers with zero mean 
! and unit variance. The seed idum must be between 0 
! and the value of mask in randu.
      REAL(KIND=GP)      :: v1,v2,ran1
      REAL(KIND=GP)      :: fac,rsq
      REAL(KIND=GP),SAVE :: gset
      INTEGER       :: idum
      INTEGER, SAVE :: iset

      IF ((iset.ne.0).or.(iset.ne.1)) iset=0
      IF (idum.lt.0) iset=0
      IF (iset.eq.0) THEN
         rsq = 2.
         DO WHILE ((rsq.ge.1.).or.(rsq.eq.0.))
            v1 = randu(idum)
            v2 = randu(idum)
            rsq = v1**2+v2**2
         END DO
         fac = sqrt(-2.*log(rsq)/rsq)
         gset = v1*fac
         randn = v2*fac
         iset = 1
      ELSE
         randn = gset
         iset = 0
      ENDIF
      END FUNCTION randn

!-----------------------------------------------------------------
      SUBROUTINE randn_cmplx(x,y,idum)
!
! Normally distributed complex random numbers with 
! unit variance. The seed idum must be between 0 
! and the value of mask in randu.
      REAL(KIND=GP), INTENT (OUT) :: x,y
      REAL(KIND=GP)       :: fac,phi
      INTEGER, INTENT(IN) :: idum

      x = randu(idum)
      y = randu(idum)
      do while (y.eq.0)
         y = randu(idum)
      end do
      phi = 2.*pi*abs(x)
      fac = sqrt(-2.*log(abs(y)))
      x = fac*cos(phi)
      y = fac*sin(phi)
      END SUBROUTINE randn_cmplx

!-----------------------------------------------------------------
      SUBROUTINE prandom_seed(iseed)
!
! Computes a unique seed for each MPI task. Note: there is a limit
! of 2^(size(integer)) for the number of tasks that can have a
! unique seed. iseed must be an integer seed.
      USE mpivars
      INTEGER, INTENT(INOUT)   :: iseed
      INTEGER, ALLOCATABLE     :: iseed1(:)
      INTEGER                  :: j,k

      iseed     = mod(iseed+myrank,abs(huge(0)-iseed)-1)
      CALL random_seed(size=k)
      ALLOCATE (iseed1(k), source = iseed*[(j, j=0, k-1)])
      CALL random_seed(put=iseed1)
      DEALLOCATE(iseed1)
      END SUBROUTINE prandom_seed

!-----------------------------------------------------------------
      SUBROUTINE prandom_number(r)
!
! Generates a random number for this MPI task on the interval [0,1]. 
! Seed for this must be generated by using prandom_seed call.
! iseed must be an integer seed.
      REAL(KIND=GP), INTENT(OUT)   :: r
      CALL random_number(r)
      END SUBROUTINE prandom_number

  END MODULE random
!=================================================================
