!=================================================================
! GHOST computation of SGS training data
!
! 2024 D. Rosenberg
!      CIRA/ NOAA
!
!=================================================================
MODULE class_GSGS
 !    USE kes
 !    USE mpivars
      USE fprecision
      USE fftplans

      IMPLICIT NONE
      INCLUDE 'mpif.h' 


!     PRIVATE
      TYPE, PUBLIC :: GSGS
        PRIVATE
        ! Member data:
        INTEGER, DIMENSION(MPI_STATUS_SIZE)          :: istatus_
        INTEGER                                      :: myrank_,nprocs_
        INTEGER                                      :: nx,ny,nz
        INTEGER                                      :: ierr_
        INTEGER                                      :: comm_
        INTEGER                                      :: ista,iend,ksta,kend
        REAL(KIND=GP), ALLOCATABLE, DIMENSION(:,:,:) :: kk2
!       REAL(KIND=GP), TARGET     , ALLOCATABLE, DIMENSION(:,:,:) :: kn2
        REAL(KIND=GP), ALLOCATABLE, DIMENSION    (:) :: kx,ky,kz
        TYPE(FFTPLAN), POINTER                       :: plancr, planrc

        CHARACTER(len=MPI_MAX_ERROR_STRING)          :: serr_
      CONTAINS
        ! Public methods:
        PROCEDURE,PUBLIC :: GSGS_ctor
        FINAL            :: GSGS_dtor
        PROCEDURE,PUBLIC :: sgsv              => GSGS_Nuu
        PROCEDURE,PUBLIC :: sgsth             => GSGS_Ntheta
        PROCEDURE,PUBLIC :: derivk3           => GSGS_derivk3
        PROCEDURE,PUBLIC :: rotor3            => GSGS_rotor3
        PROCEDURE,PUBLIC :: gradre3           => GSGS_gradre3
        PROCEDURE,PUBLIC :: nonlhd3           => GSGS_nonlhd3
        PROCEDURE,PUBLIC :: prodre3           => GSGS_prodre3
        PROCEDURE,PUBLIC :: advect3           => GSGS_advect3
        PROCEDURE,PUBLIC :: project3          => GSGS_project3
      END TYPE GSGS

      PRIVATE :: GSGS_derivk3             , GSGS_rotor3
      PRIVATE :: GSGS_gradre3             , GSGS_nonlhd3
      PRIVATE :: GSGS_prodre3             , GSGS_advect3
      PRIVATE :: GSGS_Nuu                 , GSGS_Ntheta
      PRIVATE :: GSGS_project3


! Methods:
  CONTAINS

  SUBROUTINE GSGS_ctor(this, comm, ngrid, bnds, arbsz, Dk, plancr, planrc)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  Main explicit constructor
!  ARGUMENTS:
!    this    : 'this' class instance
!    comm    : MPI communicator
!    ngrid   : array of size 3 giving grid size
!    bnds    : array of [ista, iend, ksta, kend]
!    arbsz   : arbitrary size flag (0, 1)
!    Dk      : array of size 3 giving Fourier shell widths
!    plancr,
!     planrc : FFT plans
!-----------------------------------------------------------------
 !  USE var
 !  USE grid
 !  USE boxsize
 !  USE mpivars
 !  USE random
    USE commtypes
    USE fftplans

    IMPLICIT NONE
    CLASS(GSGS)      ,INTENT(INOUT)     :: this
    INTEGER          ,INTENT   (IN)     :: comm
    INTEGER          ,INTENT   (IN)     :: arbsz
    INTEGER          ,INTENT   (IN)     :: ngrid(3)
    INTEGER          ,INTENT   (IN)     :: bnds(4)
    INTEGER                             :: aniso,i,j,k,n(3)
    REAL(KIND=GP)    ,INTENT   (IN)     :: Dk(3)
    REAL(KIND=GP)                       :: rmp,rmq,rms
    TYPE(FFTPLAN)    ,INTENT   (IN), TARGET     :: plancr, planrc


    this%comm_ = comm
    CALL MPI_COMM_SIZE(this%comm_,this%nprocs_,this%ierr_)
    CALL MPI_COMM_RANK(this%comm_,this%myrank_,this%ierr_)

    this%nx = ngrid(1)
    this%ny = ngrid(2)
    this%nz = ngrid(3)
!   CALL range(1,this%nx/2+1,this%nprocs_,this%myrank_,this%ista,this%iend)
!   CALL range(1,this%nz,this%nprocs_,this%myrank_,this%ksta,this%kend)
    this%ista = bnds(1)
    this%iend = bnds(2)
    this%ksta = bnds(3)
    this%kend = bnds(4)


!   n = ngrid
!   CALL fftp3d_create_plan_comm(this%planrc,n,FFTW_REAL_TO_COMPLEX, &
!                                FFTW_ESTIMATE, this%comm_)
!   CALL fftp3d_create_plan_comm(this%plancr,n,FFTW_COMPLEX_TO_REAL, &
!                                FFTW_ESTIMATE, this%comm_)

    this%plancr => plancr
    this%planrc => planrc


    ALLOCATE( this%kx(this%nx), this%ky(this%ny), this%kz(this%nz) )
!   ALLOCATE( this%kn2(this%nz,this%ny,this%ista:this%iend) )
    ALLOCATE( this%kk2(this%nz,this%ny,this%ista:this%iend) )
    if ( arbsz .EQ. 1 ) THEN
      aniso = 1
    ELSE
      IF ((this%nx.ne.this%ny).or.(this%ny.ne.this%nz)) THEN
         aniso = 1
      ELSE
         aniso = 0
      ENDIF
    ENDIF

    DO i = 1,this%nx/2
       this%kx(i) = real(i-1,kind=GP)
       this%kx(i+this%nx/2) = real(i-this%nx/2-1,kind=GP)
     END DO
     DO j = 1,this%ny/2
        this%ky(j) = real(j-1,kind=GP)
        this%ky(j+this%ny/2) = real(j-this%ny/2-1,kind=GP)
     END DO
     IF (this%ny.eq.1) THEN
        this%ky(1) = 0.0_GP
     ENDIF
     DO k = 1,this%nz/2
        this%kz(k) = real(k-1,kind=GP)
        this%kz(k+this%nz/2) = real(k-this%nz/2-1,kind=GP)
     END DO
     IF (aniso.eq.1) THEN
        rmp = 1.0_GP/real(this%nx,kind=GP)**2
        rmq = 1.0_GP/real(this%ny,kind=GP)**2
        rms = 1.0_GP/real(this%nz,kind=GP)**2
     ELSE
        rmp = 1.0_GP
        rmq = 1.0_GP
        rms = 1.0_GP
     ENDIF

!$omp parallel do if (this%iend-this%ista.ge.nth) private (j,k)
     DO i = this%ista,this%iend
!$omp parallel do if (this%iend-this%ista.lt.nth) private (k)
        DO j = 1,this%ny
           DO k = 1,this%nz
              this%kk2(k,j,i) = rmp*this%kx(i)**2+rmq*this%ky(j)**2+rms*this%kz(k)**2
           END DO
        END DO
     END DO

     IF ( arbsz .eq. 1 ) THEN
       this%kx = this%kx*Dk(1)
       this%ky = this%ky*Dk(2)
       this%kz = this%kz*Dk(3)
     ENDIF

     IF (aniso.eq.1) THEN
!$omp parallel do if (this%iend-this%ista.ge.nth) private (j,k)
        DO i = this%ista,this%iend
!$omp parallel do if (this%iend-this%ista.lt.nth) private (k)
           DO j = 1,this%ny
              DO k = 1,this%nz
                 this%kk2(k,j,i) = this%kx(i)**2+this%ky(j)**2+this%kz(k)**2
              END DO
           END DO
        END DO
     ENDIF

  END SUBROUTINE GSGS_ctor
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GSGS_dtor(this)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  Main explicit destructor
!  ARGUMENTS:
!    this    : 'this' class instance
!-----------------------------------------------------------------

    IMPLICIT NONE
    TYPE(GSGS)   ,INTENT(INOUT)             :: this
    INTEGER                                 :: j

    IF ( ALLOCATED    (this%kk2) ) DEALLOCATE   (this%kk2)
!   IF ( ALLOCATED    (this%kn2) ) DEALLOCATE   (this%kn2)
    IF ( ALLOCATED    (this%kx)  ) DEALLOCATE   (this%kx)
    IF ( ALLOCATED    (this%ky)  ) DEALLOCATE   (this%ky)
    IF ( ALLOCATED    (this%kz)  ) DEALLOCATE   (this%kz)

    CALL fftp3d_destroy_plan(this%plancr)
    CALL fftp3d_destroy_plan(this%planrc)

  END SUBROUTINE GSGS_dtor
!-----------------------------------------------------------------
!-----------------------------------------------------------------


  SUBROUTINE GSGS_Nuu(this, vx, vy, vz, C1, C2, C3, idir, Nt)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  Computes nonlinear term for momentum
!  ARGUMENTS:
!    this    : 'this' class instance
!    vx,vy,vz: input vector
!    C1,C2,C3: tmp arrays
!    dir     : coord direction (1, 2, 3)
!    Nt      : output field
!-----------------------------------------------------------------

    IMPLICIT NONE
    class(GSGS)   ,INTENT(INOUT)             :: this
    INTEGER         , INTENT   (IN)         :: idir
    COMPLEX(KIND=GP), INTENT   (IN), DIMENSION(this%nz,this%ny,this%ista:this%iend) :: vx,vy,vz
    COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(this%nz,this%ny,this%ista:this%iend) :: C1,C2,C3
    COMPLEX(KIND=GP), INTENT  (OUT), DIMENSION(this%nz,this%ny,this%ista:this%iend) :: Nt

     CALL this%prodre3(vx,vy,vz,C1,C2,C3)
     CALL this%nonlhd3(C1,C2,C3,Nt,idir)


  END SUBROUTINE GSGS_Nuu
!-----------------------------------------------------------------
!-----------------------------------------------------------------

!*****************************************************************
  SUBROUTINE GSGS_Ntheta(this, vx, vy, vz, th, C1, Nt)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  Computes nonlinear term for energy
!  ARGUMENTS:
!    this    : 'this' class instance
!    vx,vy,vz: input vector
!    th      : input scalar
!    C1,C2,C3: tmp arrays
!    Nt      : output scalar
!-----------------------------------------------------------------

    IMPLICIT NONE
    class(GSGS)   ,INTENT(INOUT)             :: this
    COMPLEX(KIND=GP), INTENT   (IN), DIMENSION(this%nz,this%ny,this%ista:this%iend) :: vx,vy,vz,th
    COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(this%nz,this%ny,this%ista:this%iend) :: C1
    COMPLEX(KIND=GP), INTENT  (OUT), DIMENSION(this%nz,this%ny,this%ista:this%iend) :: Nt

    CALL this%advect3(vx,vy,vz,th,Nt)


  END SUBROUTINE GSGS_Ntheta
!-----------------------------------------------------------------
!-----------------------------------------------------------------

!*****************************************************************
      SUBROUTINE GSGS_derivk3(this,a,b,dir)
!-----------------------------------------------------------------
!
! Three-dimensional derivative of the matrix 'a'
!
! Parameters
!     a  : input matrix
!     b  : at the output contains the derivative da/dk_dir
!     dir: =1 derivative in the x-direction
!          =2 derivative in the y-direction
!          =3 derivative in the z-direction
!
   !  USE kes
   !  USE var
   !  USE grid
   !  USE mpivars
!$    USE threads
      IMPLICIT NONE

      CLASS(GSGS)   ,INTENT(INOUT)   :: this
      COMPLEX(KIND=GP), INTENT (IN), DIMENSION(this%nz,this%ny,this%ista:this%iend) :: a
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(this%nz,this%ny,this%ista:this%iend) :: b
      COMPLEX(KIND=GP) :: im = (0.0_GP,1.0_GP)
      INTEGER, INTENT(IN) :: dir
      INTEGER             :: i,j,k

!
! Derivative in the x-direction
!
      IF (dir.eq.1) THEN
!$omp parallel do if (this%iend-this%ista.ge.nth) private (j,k)
         DO i = this%ista,this%iend
!$omp parallel do if (this%iend-this%ista.lt.nth) private (k)
            DO j = 1,this%ny
               DO k = 1,this%nz
                  b(k,j,i) = im*this%kx(i)*a(k,j,i)
               END DO
            END DO
         END DO
!
! Derivative in the y-direction
!
      ELSE IF (dir.eq.2) THEN
!$omp parallel do if (this%iend-this%ista.ge.nth) private (j,k)
         DO i = this%ista,this%iend
!$omp parallel do if (this%iend-this%ista.lt.nth) private (k)
            DO j = 1,this%ny
               DO k = 1,this%nz
                  b(k,j,i) = im*this%ky(j)*a(k,j,i)
               END DO
            END DO
         END DO
!
! Derivative in the z-direction
!
      ELSE
!$omp parallel do if (this%iend-this%ista.ge.nth) private (j,k)
         DO i = this%ista,this%iend
!$omp parallel do if (this%iend-this%ista.lt.nth) private (k)
            DO j = 1,this%ny
               DO k = 1,this%nz
                  b(k,j,i) = im*this%kz(k)*a(k,j,i)
               END DO
            END DO
         END DO
      ENDIF

      RETURN
      END SUBROUTINE GSGS_derivk3
!-----------------------------------------------------------------
!-----------------------------------------------------------------


!*****************************************************************
      SUBROUTINE GSGS_rotor3(this,a,b,c,dir)
!-----------------------------------------------------------------
!
! Computes the curl of the vector field A in Fourier
! space. The needed components of the field A are 
! given by the matrixes a and b, and the order must 
! follow the right hand convention.
!
! Parameters
!     a  : input matrix
!     b  : input matrix
!     c  : at the output contains curl(A)_dir
!     dir: =1 computes the x-component
!          =2 computes the y-component
!          =3 computes the z-component
!
      USE fprecision
   !  USE grid
   !  USE mpivars
!$    USE threads
      IMPLICIT NONE

      CLASS(GSGS)   ,INTENT(INOUT)   :: this
      COMPLEX(KIND=GP), INTENT (IN), DIMENSION(this%nz,this%ny,this%ista:this%iend) :: a,b
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(this%nz,this%ny,this%ista:this%iend) :: c
      COMPLEX(KIND=GP), DIMENSION(this%nz,this%ny,this%ista:this%iend) :: c1,c2
      INTEGER, INTENT(IN) :: dir
      INTEGER             :: i,j,k

!
! Computes the x-component
!
      IF (dir.eq.1) THEN
         CALL this%derivk3(a,c1,3)
         CALL this%derivk3(b,c2,2)
!$omp parallel do if (this%iend-this%ista.ge.nth) private (j,k)
         DO i = this%ista,this%iend
!$omp parallel do if (this%iend-this%ista.lt.nth) private (k)
            DO j = 1,this%ny
               DO k = 1,this%nz
                  c(k,j,i) = c2(k,j,i)-c1(k,j,i)
               END DO
            END DO
         END DO
!
! Computes the y-component
!
      ELSE IF (dir.eq.2) THEN
         CALL this%derivk3(a,c1,3)
         CALL this%derivk3(b,c2,1)
!$omp parallel do if (this%iend-this%ista.ge.nth) private (j,k)
         DO i = this%ista,this%iend
!$omp parallel do if (this%iend-this%ista.lt.nth) private (k)
            DO j = 1,this%ny
               DO k = 1,this%nz
                  c(k,j,i) = c1(k,j,i)-c2(k,j,i)
               END DO
            END DO
         END DO
!
! Computes the z-component
!
      ELSE
         CALL this%derivk3(a,c1,2)
         CALL this%derivk3(b,c2,1)
!$omp parallel do if (this%iend-this%ista.ge.nth) private (j,k)
         DO i = this%ista,this%iend
!$omp parallel do if (this%iend-this%ista.lt.nth) private (k)
            DO j = 1,this%ny
               DO k = 1,this%nz
                  c(k,j,i) = c2(k,j,i)-c1(k,j,i)
               END DO
            END DO
         END DO
      ENDIF

      RETURN
      END SUBROUTINE GSGS_rotor3
!-----------------------------------------------------------------
!-----------------------------------------------------------------

!*****************************************************************
      SUBROUTINE GSGS_gradre3(this, a,b,c,d,e,f)
!-----------------------------------------------------------------
!
! Three-dimensional inner product A.grad(A) in 
! real space. The components of the field A are 
! given by the matrixes a, b and c
!
! Parameters
!     a: input matrix in the x-direction
!     b: input matrix in the y-direction
!     c: input matrix in the z-direction
!     d: product (A.grad)A_x in Fourier space [output]
!     e: product (A.grad)A_y in Fourier space [output]
!     f: product (A.grad)A_z in Fourier space [output]
!
      USE fprecision
      USE commtypes
  !   USE mpivars
  !   USE grid
      USE fftplans
!$    USE threads
      IMPLICIT NONE

      CLASS(GSGS)   ,INTENT(INOUT)   :: this
      COMPLEX(KIND=GP), INTENT (IN), DIMENSION(this%nz,this%ny,this%ista:this%iend) :: a,b,c
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(this%nz,this%ny,this%ista:this%iend) :: d,e,f
      COMPLEX(KIND=GP), DIMENSION(this%nz,this%ny,this%ista:this%iend) :: c1,c2
      COMPLEX(KIND=GP), DIMENSION(this%nz,this%ny,this%ista:this%iend) :: c3,c4
      REAL(KIND=GP), DIMENSION(this%nx,this%ny,this%ksta:this%kend)    :: r1,r2
      REAL(KIND=GP), DIMENSION(this%nx,this%ny,this%ksta:this%kend)    :: r3,r4
      REAL(KIND=GP), DIMENSION(this%nx,this%ny,this%ksta:this%kend)    :: rx,ry,rz
      REAL(KIND=GP)    :: tmp
      INTEGER :: i,j,k

!
! Computes (A_x.dx)A_dir
!
      c1 = a
      CALL this%derivk3(a,c2,1)
      CALL this%derivk3(b,c3,1)
      CALL this%derivk3(c,c4,1)
      CALL fftp3d_complex_to_real(this%plancr,c1,r1,this%comm_)
      CALL fftp3d_complex_to_real(this%plancr,c2,r2,this%comm_)
      CALL fftp3d_complex_to_real(this%plancr,c3,r3,this%comm_)
      CALL fftp3d_complex_to_real(this%plancr,c4,r4,this%comm_)

!$omp parallel do if (this%kend-this%ksta.ge.nth) private (j,i)
      DO k = this%ksta,this%kend
!$omp parallel do if (this%kend-this%ksta.lt.nth) private (i)
         DO j = 1,this%ny
            DO i = 1,this%nx
               rx(i,j,k) = r1(i,j,k)*r2(i,j,k)
               ry(i,j,k) = r1(i,j,k)*r3(i,j,k)
               rz(i,j,k) = r1(i,j,k)*r4(i,j,k)
            END DO
         END DO
      END DO
!
! Computes (A_y.dy)A_dir
!
      c1 = b
      CALL this%derivk3(a,c2,2)
      CALL this%derivk3(b,c3,2)
      CALL this%derivk3(c,c4,2)
      CALL fftp3d_complex_to_real(this%plancr,c1,r1,this%comm_)
      CALL fftp3d_complex_to_real(this%plancr,c2,r2,this%comm_)
      CALL fftp3d_complex_to_real(this%plancr,c3,r3,this%comm_)
      CALL fftp3d_complex_to_real(this%plancr,c4,r4,this%comm_)

!$omp parallel do if (this%kend-this%ksta.ge.nth) private (j,i)
      DO k = this%ksta,this%kend
!$omp parallel do if (this%kend-this%ksta.lt.nth) private (i)
         DO j = 1,this%ny
            DO i = 1,this%nx
               rx(i,j,k) = rx(i,j,k)+r1(i,j,k)*r2(i,j,k)
               ry(i,j,k) = ry(i,j,k)+r1(i,j,k)*r3(i,j,k)
               rz(i,j,k) = rz(i,j,k)+r1(i,j,k)*r4(i,j,k)
            END DO
         END DO
      END DO
!
! Computes (A_z.dz)A_dir
!
      c1 = c
      CALL this%derivk3(a,c2,3)
      CALL this%derivk3(b,c3,3)
      CALL this%derivk3(c,c4,3)
      CALL fftp3d_complex_to_real(this%plancr,c1,r1,this%comm_)
      CALL fftp3d_complex_to_real(this%plancr,c2,r2,this%comm_)
      CALL fftp3d_complex_to_real(this%plancr,c3,r3,this%comm_)
      CALL fftp3d_complex_to_real(this%plancr,c4,r4,this%comm_)

      tmp = 1.0_GP/ &
            (real(this%nx,kind=GP)*real(this%ny,kind=GP)*real(this%nz,kind=GP))**2
!$omp parallel do if (this%kend-this%ksta.ge.nth) private (j,i)
      DO k = this%ksta,this%kend
!$omp parallel do if (this%kend-this%ksta.lt.nth) private (i)
         DO j = 1,this%ny
            DO i = 1,this%nx
               rx(i,j,k) = (rx(i,j,k)+r1(i,j,k)*r2(i,j,k))*tmp
               ry(i,j,k) = (ry(i,j,k)+r1(i,j,k)*r3(i,j,k))*tmp
               rz(i,j,k) = (rz(i,j,k)+r1(i,j,k)*r4(i,j,k))*tmp
            END DO
         END DO
      END DO

      CALL fftp3d_real_to_complex(this%planrc,rx,d,this%comm_)
      CALL fftp3d_real_to_complex(this%planrc,ry,e,this%comm_)
      CALL fftp3d_real_to_complex(this%planrc,rz,f,this%comm_)

      RETURN
      END SUBROUTINE GSGS_gradre3
!-----------------------------------------------------------------
!-----------------------------------------------------------------

!*****************************************************************
      SUBROUTINE GSGS_prodre3(this,a,b,c,d,e,f)
!-----------------------------------------------------------------
!
! Three-dimensional cross product curl(A)xA in 
! real space. The components of the field A are 
! given by the matrixes a, b and c
!
! Parameters
!     a  : input matrix in the x-direction
!     b  : input matrix in the y-direction
!     c  : input matrix in the z-direction
!     d  : product [curl(A)xA]_x in Fourier space [output]
!     e  : product [curl(A)xA]_y in Fourier space [output]
!     f  : product [curl(A)xA]_z in Fourier space [output]
!
      USE fprecision
      USE commtypes
!     USE mpivars
!     USE grid
      USE fftplans
!$    USE threads
      IMPLICIT NONE

      CLASS(GSGS)   ,INTENT(INOUT)   :: this
      COMPLEX(KIND=GP), INTENT (IN), DIMENSION(this%nz,this%ny,this%ista:this%iend) :: a,b,c
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(this%nz,this%ny,this%ista:this%iend) :: d,e,f
      REAL(KIND=GP), DIMENSION(this%nx,this%ny,this%ksta:this%kend)    :: r1,r2
      REAL(KIND=GP), DIMENSION(this%nx,this%ny,this%ksta:this%kend)    :: r3,r4
      REAL(KIND=GP), DIMENSION(this%nx,this%ny,this%ksta:this%kend)    :: r5,r6
      REAL(KIND=GP), DIMENSION(this%nx,this%ny,this%ksta:this%kend)    :: r7
      REAL(KIND=GP)    :: tmp
      INTEGER :: i,j,k

!
! Computes curl(A)
!
      CALL this%rotor3(b,c,d,1)
      CALL this%rotor3(a,c,e,2)
      CALL this%rotor3(a,b,f,3)
      CALL fftp3d_complex_to_real(this%plancr,d,r1,this%comm_)
      CALL fftp3d_complex_to_real(this%plancr,e,r2,this%comm_)
      CALL fftp3d_complex_to_real(this%plancr,f,r3,this%comm_)

!
! Computes A
!
!$omp parallel do if (this%iend-this%ista.ge.nth) private (j,k)
      DO i = this%ista,this%iend
!$omp parallel do if (this%iend-this%ista.lt.nth) private (k)
         DO j = 1,this%ny
            DO k = 1,this%nz
               d(k,j,i) = a(k,j,i)
               e(k,j,i) = b(k,j,i)
               f(k,j,i) = c(k,j,i)
            END DO
         END DO
      END DO
      CALL fftp3d_complex_to_real(this%plancr,d,r4,this%comm_)
      CALL fftp3d_complex_to_real(this%plancr,e,r5,this%comm_)
      CALL fftp3d_complex_to_real(this%plancr,f,r6,this%comm_)
!
! Computes curl(A)xA
!
      tmp = 1.0_GP/ &
            (real(this%nx,kind=GP)*real(this%ny,kind=GP)*real(this%nz,kind=GP))**2
!$omp parallel do if (this%kend-this%ksta.ge.nth) private (j,i)
      DO k = this%ksta,this%kend
!$omp parallel do if (this%kend-this%ksta.lt.nth) private (i)
         DO j = 1,this%ny
            DO i = 1,this%nx
               r7(i,j,k) = (r2(i,j,k)*r6(i,j,k)-r5(i,j,k)*r3(i,j,k))*tmp
               r3(i,j,k) = (r3(i,j,k)*r4(i,j,k)-r6(i,j,k)*r1(i,j,k))*tmp
               r1(i,j,k) = (r1(i,j,k)*r5(i,j,k)-r4(i,j,k)*r2(i,j,k))*tmp
            END DO
         END DO
      END DO

      CALL fftp3d_real_to_complex(this%planrc,r7,d,this%comm_)
      CALL fftp3d_real_to_complex(this%planrc,r3,e,this%comm_)
      CALL fftp3d_real_to_complex(this%planrc,r1,f,this%comm_)

      RETURN
      END SUBROUTINE GSGS_prodre3
!-----------------------------------------------------------------
!-----------------------------------------------------------------

!*****************************************************************
      SUBROUTINE GSGS_nonlhd3(this,a,b,c,g,dir)
!-----------------------------------------------------------------
!
! Computes the nonlinear terms in 3D Navier-Stokes 
! equation. It takes the components of (v.grad)v or 
! curl(v)xv as input matrixes and computes 
! -(v.grad)v-grad(p) or -curl(v)xv-grad(p) in 
! Fourier space, with the pressure chosen to satisfy 
! the incompressibility condition.
!
! Parameters
!     a  : input matrix in the x-direction
!     b  : input matrix in the y-direction
!     c  : input matrix in the z-direction
!     g  : at the output contains the result in Fourier space
!     dir: =1 computes the x-component
!          =2 computes the y-component
!          =3 computes the z-component
!
      USE fprecision
   !  USE kes
   !  USE grid
   !  USE mpivars
      USE commtypes
!$    USE threads
      IMPLICIT NONE

      CLASS(GSGS)   ,INTENT(INOUT)   :: this
      COMPLEX(KIND=GP), INTENT (IN), DIMENSION(this%nz,this%ny,this%ista:this%iend) :: a,b,c
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(this%nz,this%ny,this%ista:this%iend) :: g
      INTEGER, INTENT(IN) :: dir
      INTEGER             :: i,j,k

!
! Computes the x-component
!
      IF (dir.eq.1) THEN
         IF (this%ista.eq.1) THEN
!$omp parallel do private (k)
            DO j = 1,this%ny
               DO k = 1,this%nz
                  g(k,j,1) = -a(k,j,1)
               END DO
            END DO
!$omp parallel do if (this%iend-2.ge.nth) private (j,k)
            DO i = 2,this%iend
!$omp parallel do if (this%iend-2.lt.nth) private (k)
               DO j = 1,this%ny
                  DO k = 1,this%nz
                     g(k,j,i) = -a(k,j,i)+this%kx(i)*(this%kx(i)*a(k,j,i) &
                       +this%ky(j)*b(k,j,i)+this%kz(k)*c(k,j,i))/this%kk2(k,j,i)
                  END DO
               END DO
            END DO
         ELSE
!$omp parallel do if (this%iend-this%ista.ge.nth) private (j,k)
            DO i = this%ista,this%iend
!$omp parallel do if (this%iend-this%ista.lt.nth) private (k)
               DO j = 1,this%ny
                  DO k = 1,this%nz
                     g(k,j,i) = -a(k,j,i)+this%kx(i)*(this%kx(i)*a(k,j,i) &
                       +this%ky(j)*b(k,j,i)+this%kz(k)*c(k,j,i))/this%kk2(k,j,i)
                  END DO
               END DO
            END DO
         ENDIF
!
! Computes the y-component
!
      ELSE IF (dir.eq.2) THEN
!$omp parallel do if (this%iend-this%ista.ge.nth) private (k)
         DO i = this%ista,this%iend
!$omp parallel do if (this%iend-this%ista.lt.nth)
            DO k = 1,this%nz
               g(k,1,i) = -b(k,1,i)
            END DO
         END DO
!$omp parallel do if (this%iend-this%ista.ge.nth) private (j,k)
         DO i = this%ista,this%iend
!$omp parallel do if (this%iend-this%ista.lt.nth) private (k)
            DO j = 2,this%ny
               DO k = 1,this%nz
                  g(k,j,i) = -b(k,j,i)+this%ky(j)*(this%kx(i)*a(k,j,i) &
                    +this%ky(j)*b(k,j,i)+this%kz(k)*c(k,j,i))/this%kk2(k,j,i)
               END DO
            END DO
         END DO
!
! Computes the z-component
!
      ELSE

!$omp parallel do if (this%iend-this%ista.ge.nth) private (j)
         DO i = this%ista,this%iend
!$omp parallel do if (this%iend-this%ista.lt.nth)
            DO j = 1,this%ny
               g(1,j,i) = -c(1,j,i)
            END DO
         END DO
!$omp parallel do if (this%iend-this%ista.ge.nth) private (j,k)
         DO i = this%ista,this%iend
!$omp parallel do if (this%iend-this%ista.lt.nth) private (k)
            DO j = 1,this%ny
               DO k = 2,this%nz
                  g(k,j,i) = -c(k,j,i)+this%kz(k)*(this%kx(i)*a(k,j,i) &
                    +this%ky(j)*b(k,j,i)+this%kz(k)*c(k,j,i))/this%kk2(k,j,i)
               END DO
            END DO
         END DO
      ENDIF
      RETURN
      END SUBROUTINE GSGS_nonlhd3
!-----------------------------------------------------------------
!-----------------------------------------------------------------

!*****************************************************************
      SUBROUTINE GSGS_advect3(this,a,b,c,d,e)
!-----------------------------------------------------------------
!
! Three-dimensional inner product -A.grad(B) in 
! real space. The components of the field A are 
! given by the arrays a, b and c, B is a scalar 
! quantity given by d.
!
! Parameters
!     a: input matrix in the x-direction
!     b: input matrix in the y-direction
!     c: input matrix in the z-direction
!     d: input matrix with the scalar
!     e: product (A.grad)B in Fourier space [output]
!
      USE fprecision
      USE commtypes
   !  USE mpivars
   !  USE grid
      USE fftplans
!$    USE threads
      IMPLICIT NONE

      CLASS(GSGS)   ,INTENT(INOUT)   :: this
      COMPLEX(KIND=GP), INTENT (IN), DIMENSION(this%nz,this%ny,this%ista:this%iend) :: a,b
      COMPLEX(KIND=GP), INTENT (IN), DIMENSION(this%nz,this%ny,this%ista:this%iend) :: c,d
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(this%nz,this%ny,this%ista:this%iend) :: e
      COMPLEX(KIND=GP), DIMENSION(this%nz,this%ny,this%ista:this%iend) :: c1,c2
      REAL(KIND=GP),    DIMENSION(this%nx,this%ny,this%ksta:this%kend) :: r1,r2
      REAL(KIND=GP),    DIMENSION(this%nx,this%ny,this%ksta:this%kend) :: r3
      REAL(KIND=GP)    :: tmp
      INTEGER :: i,j,k

!
! Computes (A_x.dx)B
!
      c1 = a
      CALL this%derivk3(d,c2,1)
      CALL fftp3d_complex_to_real(this%plancr,c1,r1,this%comm_)
      CALL fftp3d_complex_to_real(this%plancr,c2,r2,this%comm_)

!$omp parallel do if (this%kend-this%ksta.ge.nth) private (j,i)
      DO k = this%ksta,this%kend
!$omp parallel do if (this%kend-this%ksta.lt.nth) private (i)
         DO j = 1,this%ny
            DO i = 1,this%nx
               r3(i,j,k) = r1(i,j,k)*r2(i,j,k)
            END DO
         END DO
      END DO
!
! Computes (A_y.dy)B
!
      c1 = b
      CALL this%derivk3(d,c2,2)
      CALL fftp3d_complex_to_real(this%plancr,c1,r1,this%comm_)
      CALL fftp3d_complex_to_real(this%plancr,c2,r2,this%comm_)

!$omp parallel do if (this%kend-this%ksta.ge.nth) private (j,i)
      DO k = this%ksta,this%kend
!$omp parallel do if (this%kend-this%ksta.lt.nth) private (i)
         DO j = 1,this%ny
            DO i = 1,this%nx
               r3(i,j,k) = r3(i,j,k)+r1(i,j,k)*r2(i,j,k)
            END DO
         END DO
      END DO
!
! Computes (A_z.dz)B
!
      c1 = c
      CALL this%derivk3(d,c2,3)
      CALL fftp3d_complex_to_real(this%plancr,c1,r1,this%comm_)
      CALL fftp3d_complex_to_real(this%plancr,c2,r2,this%comm_)

! We need -A.grad(B)
      tmp = -1.0_GP/ &
            (real(this%nx,kind=GP)*real(this%ny,kind=GP)*real(this%nz,kind=GP))**2
!$omp parallel do if (this%kend-this%ksta.ge.nth) private (j,i)
      DO k = this%ksta,this%kend
!$omp parallel do if (this%kend-this%ksta.lt.nth) private (i)
         DO j = 1,this%ny
            DO i = 1,this%nx
               r3(i,j,k) = (r3(i,j,k)+r1(i,j,k)*r2(i,j,k))*tmp
            END DO
         END DO
      END DO

      CALL fftp3d_real_to_complex(this%planrc,r3,e,this%comm_)

      RETURN
      END SUBROUTINE GSGS_advect3


!*****************************************************************
      SUBROUTINE GSGS_project3(this,a,b,c,d,e,f)
!-----------------------------------------------------------------
!
! Project 3D vector field, (a,b,c) to incompressible space,
! returning projected components into input vector
!
! Parameters
!     a,b,c : input vector field
!     d,e,f : output field
!
      USE fprecision
      USE commtypes
!$    USE threads
      IMPLICIT NONE

      CLASS(GSGS)   ,INTENT(INOUT)   :: this
      COMPLEX(KIND=GP), INTENT  (IN), DIMENSION(this%nz,this%ny,this%ista:this%iend) :: a,b,c
      COMPLEX(KIND=GP), INTENT (OUT), DIMENSION(this%nz,this%ny,this%ista:this%iend) :: d,e,f
      REAL(KIND=GP)    :: tmp
      INTEGER :: i,j,k

        ! vx:
        IF (this%ista.eq.1) THEN
!$omp parallel do private (k)
            DO j = 1,this%ny
               DO k = 1,this%nz
                  d(k,j,1) = a(k,j,1)
               END DO
            END DO
!$omp parallel do if (this%iend-2.ge.nth) private (j,k)
            DO i = 2,this%iend
!$omp parallel do if (this%iend-2.lt.nth) private (k)
               DO j = 1,this%ny
                  DO k = 1,this%nz
                     d(k,j,i) = (1.0 - this%kx(i)*this%kx(i)/this%kk2(k,j,i) ) * a(k,j,i) &
                                     - this%kx(i)*this%ky(j)/this%kk2(k,j,i)   * b(k,j,i) &
                                     - this%kx(i)*this%kz(k)/this%kk2(k,j,i)   * c(k,j,i)  
                  END DO
               END DO
            END DO
         ELSE
!$omp parallel do if (this%iend-this%ista.ge.nth) private (j,k)
            DO i = this%ista,this%iend
!$omp parallel do if (this%iend-this%ista.lt.nth) private (k)
               DO j = 1,this%ny
                  DO k = 1,this%nz
                     d(k,j,i) = (1.0 - this%kx(i)*this%kx(i)/this%kk2(k,j,i) ) * a(k,j,i) &
                                     - this%kx(i)*this%ky(j)/this%kk2(k,j,i)   * b(k,j,i) &
                                     - this%kx(i)*this%kz(k)/this%kk2(k,j,i)   * c(k,j,i)  
                  END DO
               END DO
            END DO
         ENDIF

         ! vy:
!$omp parallel do if (this%iend-this%ista.ge.nth) private (k)
         DO i = this%ista,this%iend
!$omp parallel do if (this%iend-this%ista.lt.nth)
            DO k = 1,this%nz
               g(k,1,i) = b(k,1,i)
            END DO
         END DO
!$omp parallel do if (this%iend-this%ista.ge.nth) private (j,k)
         DO i = this%ista,this%iend
!$omp parallel do if (this%iend-this%ista.lt.nth) private (k)
            DO j = 2,this%ny
               DO k = 1,this%nz
                  e(k,j,i) =      - this%ky(j)*this%kx(i)/this%kk2(k,j,i)   * a(k,j,i) &
                           + (1.0 - this%ky(j)*this%ky(j)/this%kk2(k,j,i) ) * b(k,j,i) &
                                  - this%ky(j)*this%kz(k)/this%kk2(k,j,i)   * c(k,j,i)  
               END DO
            END DO
         END DO


         ! vz:
!$omp parallel do if (this%iend-this%ista.ge.nth) private (j)
         DO i = this%ista,this%iend
!$omp parallel do if (this%iend-this%ista.lt.nth)
            DO j = 1,this%ny
               f(1,j,i) = c(1,j,i)
            END DO
         END DO
!$omp parallel do if (this%iend-this%ista.ge.nth) private (j,k)
         DO i = this%ista,this%iend
!$omp parallel do if (this%iend-this%ista.lt.nth) private (k)
            DO j = 1,this%ny
               DO k = 2,this%nz
                  f(k,j,i) =      - this%kz(k)*this%kx(i)/this%kk2(k,j,i)   * a(k,j,i) &
                                  - this%kz(k)*this%ky(j)/this%kk2(k,j,i)   * b(k,j,i) &
                           + (1.0 - this%kz(k)*this%kz(k)/this%kk2(k,j,i) ) * c(k,j,i) 
               END DO
            END DO
         END DO


      END SUBROUTINE GSGS_project3
  
END MODULE class_GSGS
