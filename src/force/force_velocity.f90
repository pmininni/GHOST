! =====================================================================
! NAME       : force_velocity.f90
! DESCRIPTION: Forcing methods for all VelocityBase solver classes.
!              Each force needs an initial set up, and an update.
!              These forcing methods can be used to initialize and
!              update forcings on other solver classes, using
!              composition of operators, as long as no correlation
!              between the forces is needed.
!
! Forces avaliable:
!   null_fv    : Null forcing of the velocity field
!   tg_fv      : Taylor-Green forcing
!   abc_fv     : ABC forcing
!   random_fv  : Random forcing of the velocity field
!
! Update methods available:
!   constant_fv: Constant forcing (no method)
!   shift_fv   : Instantaneous phase shift
!   shuffle_fv : Slowly evolving random phases
!
! DATE       : 01/18/26 (PDM)
! =====================================================================

module force_velocity
  USE forcebase_mod

  IMPLICIT NONE

  ! ================= Forcing functions supported =====================
  type, extends(forceBase) :: forceNull_fv
    contains
      procedure :: init_GForce => init_nullfv
  end type forceNull_fv 
  type, extends(forceBase) :: forceTg_fv
    contains
      procedure :: init_GForce => init_tgfv
  end type forceTg_fv 
  type, extends(forceBase) :: forceAbc_fv
    contains
      procedure :: init_GForce => init_abcfv
  end type forceAbc_fv 
  type, extends(forceBase) :: forceRandom_fv
    contains
      procedure :: init_GForce => init_randomfv
  end type forceRandom_fv
! type, extends(forceBase) :: forceUserdef_fv
!   contains
!     procedure :: init_GForce => init_userdefv
! end type forceUserdef_fv

  ! ================= Update methods supported =======================
  type, extends(forceUpdt) :: shiftupdt_fv
    contains
      procedure :: update_GForce => update_shiftfv
  end type shiftupdt_fv 
  type, extends(forceUpdt) :: shuffleupdt_fv
    contains
      procedure :: update_GForce => update_shufflefv
  end type shuffleupdt_fv
! type, extends(forceUpdt) :: userupdt_fv
!   contains
!     procedure :: init_GForce => init_userupdtfv
! end type userupdt_fv

CONTAINS

  ! ===================================================================
  ! Forcing functions
  ! ===================================================================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Null forcing
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_nullfv(this,solver,state)
    use gstate_mod
    use hd_mod
    use grid
    use mpivars
!$  use threads
    implicit none
    
    class(forceNull_fv), intent   (in) :: this
    class(EquationBase), intent   (in) :: solver
    type   (GStateComp), intent(inout) :: state(:)
    integer                            :: i,j,k

    select type (solver)
    class is (VelocityBase)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
        DO j = 1,ny
          DO k = 1,nz
            state(solver%VELOCITY  )%ccomp(k,j,i) = 0.0_GP
            state(solver%VELOCITY+1)%ccomp(k,j,i) = 0.0_GP
            state(solver%VELOCITY+2)%ccomp(k,j,i) = 0.0_GP
          END DO
        END DO
      END DO
    class default
      error stop "This solver does not support velocity forcing"
    end select
  end subroutine init_nullfv


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Taylor-Green forcing
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   f0  : amplitude of the forcing
  !   kdn : minimum wave number (rounded to next integer)
  !   kup : maximum wave number (rounded to next integer)
  subroutine init_tgfv(this,solver,state)
    use gstate_mod
    use hd_mod
    use grid
    use boxsize
    use commtypes
    use fft
    use ali
    use var
    use pseudospec_norm
!$  use threads
    implicit none
    
    class  (forceTg_fv), intent   (in)       :: this
    class(EquationBase), intent   (in)       :: solver
    type   (GStateComp), intent(inout)       :: state(:)
    real(kind=GP), pointer, dimension(:,:,:) :: R1,R2
    real(kind=GP)                            :: f0,kdn,kup
    integer                                  :: i,j,k,ki
    logical                                  :: bret

    namelist/ tg_fv / f0,kdn,kup
    CALL solver%workspace_%get_real_tmp(R1,bret)
    CALL solver%workspace_%get_real_tmp(R2,bret)
    select type (solver)
    class is (VelocityBase)
    ! Read parameters from a namelist in the input file
    if ( myrank .eq. 0 ) then
      open(1,file=solver%infile_,status='unknown',form="formatted")
      read(1,NML=tg_fv)
      close(1)
    endif
    call mpi_bcast(f0 ,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(kup,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(kdn,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
    ! Generate a superposition of TG flows
    IF ( abs(Lx-Ly).gt.tiny ) THEN
      IF (myrank.eq.0) error stop "TG initial conditions require Lx=Ly"
    ENDIF
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
    DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
       DO j = 1,ny
          DO k = 1,nz
            state(solver%VELOCITY+2)%ccomp(k,j,i) = 0.0_GP !fz
          END DO
       END DO
    END DO
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
    DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
       DO j = 1,ny
          DO i = 1,nx
          R1(i,j,k) = 0.0_GP
          R2(i,j,k) = 0.0_GP
          DO ki = INT(kdn),INT(kup)
             R1(i,j,k) = R1(i,j,k)+SIN(2*pi*ki*(real(i,kind=GP)-1)/ &
                      real(nx,kind=GP))*COS(2*pi*ki*(real(j,kind=GP)-1)/ &
                      real(ny,kind=GP))*COS(2*pi*ki*(real(k,kind=GP)-1)/ &
                      real(nz,kind=GP))
             R2(i,j,k) = R2(i,j,k)-COS(2*pi*ki*(real(i,kind=GP)-1)/ &
                      real(nx,kind=GP))*SIN(2*pi*ki*(real(j,kind=GP)-1)/ &
                      real(ny,kind=GP))*COS(2*pi*ki*(real(k,kind=GP)-1)/ &
                      real(nz,kind=GP))
          END DO
          END DO
       END DO
    END DO
    CALL fftp3d_real_to_complex(planrc,R1, & !fx
                   state(solver%VELOCITY  )%ccomp,MPI_COMM_WORLD)
    CALL fftp3d_real_to_complex(planrc,R2, & !fy
                   state(solver%VELOCITY+1)%ccomp,MPI_COMM_WORLD)
    CALL normalize(state(solver%VELOCITY  )%ccomp, &
                   state(solver%VELOCITY+1)%ccomp, &
                   state(solver%VELOCITY+2)%ccomp, &
                   f0,1,MPI_COMM_WORLD)
    class default
      error stop "This solver does not support velocity forcing"
    end select
    CALL solver%workspace_%free_real_tmp(R1)
    CALL solver%workspace_%free_real_tmp(R2)
  end subroutine init_tgfv


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! ABC forcing
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   f0  : normalized amplitude of the forcing
  !   kdn : minimum wave number (rounded to next integer)
  !   kup : maximum wave number (rounded to next integer)
  !   A   : A amplitude (default =1)
  !   B   : B amplitude (default =1)
  !   C   : C amplitude (default =1)
  subroutine init_abcfv(this,solver,state)
    use gstate_mod
    use hd_mod
    use grid
    use boxsize
    use commtypes
    use fft
    use ali
    use var
    use pseudospec_norm
!$  use threads
    implicit none

    class (forceAbc_fv), intent   (in)       :: this
    class(EquationBase), intent   (in)       :: solver
    type   (GStateComp), intent(inout)       :: state(:)
    real(kind=GP), pointer, dimension(:,:,:) :: R1,R2,R3
    real(kind=GP)                            :: f0,kdn,kup
    real(kind=GP)                            :: A,B,C
    integer                                  :: i,j,k,ki
    logical                                  :: bret

    namelist/ abc_fv / f0,kdn,kup,A,B,C
    CALL solver%workspace_%get_real_tmp(R1,bret)
    CALL solver%workspace_%get_real_tmp(R2,bret)
    CALL solver%workspace_%get_real_tmp(R3,bret)
    select type (solver)
    class is (VelocityBase)
    ! Read parameters from a namelist in the input file
    A = 1.0_GP ; B = 1.0_GP ; C = 1.0_GP
    if ( myrank .eq. 0 ) then
      open(1,file=solver%infile_,status='unknown',form="formatted")
      read(1,NML=abc_fv)
      close(1)
    endif
    call mpi_bcast(f0 ,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(kup,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(kdn,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(A  ,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(B  ,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(C  ,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
    ! Generate a superposition of ABC flows
    if ( (abs(Lx-Ly).gt.tiny).or.(abs(Lx-Lz).gt.tiny) ) then
      if (myrank.eq.0) error stop 'ABC forcing requires Lx=Ly=Lz'
    endif
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
    DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
      DO j = 1,ny
        DO i = 1,nx
          R1(i,j,k) = 0.
          R2(i,j,k) = 0.
          R3(i,j,k) = 0.
          DO ki = INT(kdn),INT(kup)
          R1(i,j,k) = R1(i,j,k)+(B*COS(2*pi*ki*(real(j,kind=GP)-1)/ &
               real(ny,kind=GP))+C*SIN(2*pi*ki*(real(k,kind=GP)-1)/ &
               real(nz,kind=GP)))/ki**2
          R2(i,j,k) = R2(i,j,k)+(A*SIN(2*pi*ki*(real(i,kind=GP)-1)/ &
               real(nx,kind=GP))+C*COS(2*pi*ki*(real(k,kind=GP)-1)/ &
               real(nz,kind=GP)))/ki**2 
          R3(i,j,k) = R3(i,j,k)+(A*COS(2*pi*ki*(real(i,kind=GP)-1)/ &
               real(nx,kind=GP))+B*SIN(2*pi*ki*(real(j,kind=GP)-1)/ &
               real(ny,kind=GP)))/ki**2
          END DO
        END DO
      END DO
    END DO
    call fftp3d_real_to_complex(planrc,R1, & !fx
                   state(solver%VELOCITY  )%ccomp,MPI_COMM_WORLD)
    call fftp3d_real_to_complex(planrc,R2, & !fy
                   state(solver%VELOCITY+1)%ccomp,MPI_COMM_WORLD)
    call fftp3d_real_to_complex(planrc,R3, & !fz
                   state(solver%VELOCITY+2)%ccomp,MPI_COMM_WORLD)
    call normalize(state(solver%VELOCITY  )%ccomp, &
                   state(solver%VELOCITY+1)%ccomp, &
                   state(solver%VELOCITY+2)%ccomp, &
                   f0,1,MPI_COMM_WORLD)
    class default
      error stop "This solver does not support velocity field ICs"
    end select
    call solver%workspace_%free_real_tmp(R1)
    call solver%workspace_%free_real_tmp(R2)
    call solver%workspace_%free_real_tmp(R3)
  end subroutine init_abcfv

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Random forcing: Superposition of harmonic modes with 
  !! random phases correlated to give a tunable total 
  !! relative helicity [following Pouquet & Patterson, JFM 
  !! 1978] with a ~k^(-4) slope
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   f0   : velocity amplitude
  !   kdn  : minimum wave number (rounded to next integer)
  !   kup  : maximum wave number (rounded to next integer)
  !   alpha: sin(2*alpha) is the total relative helicity
  subroutine init_randomfv(this,solver,state)
    use gstate_mod
    use hd_mod
    use random
    use grid
    use boxsize
    use commtypes
    use status
    use fft
    use ali
    use kes
    use pseudospec_norm
!$  use threads
    implicit none

    class(forceRandom_fv), intent   (in)          :: this
    class  (EquationBase), intent   (in)          :: solver
    type     (GStateComp), intent(inout)          :: state(:)
    complex  (kind=GP), pointer, dimension(:,:,:) :: C1,C2,C3,C4
    complex  (kind=GP), pointer, dimension(:,:,:) :: C5,C6,C7,C8
    real     (kind=GP)                            :: f0,kdn,kup
    real     (kind=GP)                            :: alpha,a1,a2
    real     (kind=GP)                            :: dump,phase
    integer                                       :: i,j,k
    logical                                       :: bret

    namelist/ random_fv / f0,kdn,kup,alpha
    CALL solver%workspace_%get_complex_tmp(C1,bret)
    CALL solver%workspace_%get_complex_tmp(C2,bret)
    CALL solver%workspace_%get_complex_tmp(C3,bret)
    CALL solver%workspace_%get_complex_tmp(C4,bret)
    CALL solver%workspace_%get_complex_tmp(C5,bret)
    CALL solver%workspace_%get_complex_tmp(C6,bret)
    CALL solver%workspace_%get_complex_tmp(C7,bret)
    CALL solver%workspace_%get_complex_tmp(C8,bret)
    select type (solver)
    class is (VelocityBase)
    ! Read parameters from a namelist in the input file
    alpha = 0.0_GP
    if ( myrank .eq. 0 ) then
      open(1,file=solver%infile_,status='unknown',form="formatted")
      read(1,NML=random_fv)
      close(1)
    endif
    call mpi_bcast(f0   ,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(kup  ,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(kdn  ,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(alpha,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
    ! Generate the first random velocity field
    IF (ista.eq.1) THEN
      C1(1,1,1) = 0.
      C2(1,1,1) = 0.
      C3(1,1,1) = 0.
      DO j = 2,ny/2+1
        IF ((kk2(1,j,1).le.kup**2).and.(kk2(1,j,1).ge.kdn**2)) THEN
          dump = 1./sqrt(kk2(1,j,1))**4
          phase = 2*pi*randu(seed)
          C1(1,j,1) = (COS(phase)+im*SIN(phase))*dump
          C1(1,ny-j+2,1) = conjg(C1(1,j,1))
          phase = 2*pi*randu(seed)
          C2(1,j,1) = (COS(phase)+im*SIN(phase))*dump
          C2(1,ny-j+2,1) = conjg(C2(1,j,1))
          phase = 2*pi*randu(seed)
          C3(1,j,1) = (COS(phase)+im*SIN(phase))*dump
          C3(1,ny-j+2,1) = conjg(C3(1,j,1))
        ELSE
          C1(1,j,1) = 0.
          C1(1,ny-j+2,1) = 0.
          C2(1,j,1) = 0.
          C2(1,ny-j+2,1) = 0.
          C3(1,j,1) = 0.
          C3(1,ny-j+2,1) = 0.
        ENDIF
      END DO
      DO k = 2,nz/2+1
        IF ((kk2(k,1,1).le.kup**2).and.(kk2(k,1,1).ge.kdn**2)) THEN
          dump = 1./sqrt(kk2(k,1,1))**4
          phase = 2*pi*randu(seed)
          C1(k,1,1) = (COS(phase)+im*SIN(phase))*dump
          C1(nz-k+2,1,1) = conjg(C1(k,1,1))
          phase = 2*pi*randu(seed)
          C2(k,1,1) = (COS(phase)+im*SIN(phase))*dump
          C2(nz-k+2,1,1) = conjg(C2(k,1,1))
          phase = 2*pi*randu(seed)
          C3(k,1,1) = (COS(phase)+im*SIN(phase))*dump
          C3(nz-k+2,1,1) = conjg(C3(k,1,1))
        ELSE
          C1(k,1,1) = 0.
          C1(nz-k+2,1,1) = 0.
          C2(k,1,1) = 0.
          C2(nz-k+2,1,1) = 0.
          C3(k,1,1) = 0.
          C3(nz-k+2,1,1) = 0.
        ENDIF
      END DO
      DO j = 2,ny
        DO k = 2,nz/2+1
          IF ((kk2(k,j,1).le.kup**2).and.(kk2(k,j,1).ge.kdn**2)) THEN
            dump = 1./sqrt(kk2(k,j,1))**4
            phase = 2*pi*randu(seed)
            C1(k,j,1) = (COS(phase)+im*SIN(phase))*dump
            C1(nz-k+2,ny-j+2,1) = conjg(C1(k,j,1))
            phase = 2*pi*randu(seed)
            C2(k,j,1) = (COS(phase)+im*SIN(phase))*dump
            C2(nz-k+2,ny-j+2,1) = conjg(C2(k,j,1))
            phase = 2*pi*randu(seed)
            C3(k,j,1) = (COS(phase)+im*SIN(phase))*dump
            C3(nz-k+2,ny-j+2,1) = conjg(C3(k,j,1))
          ELSE
            C1(k,j,1) = 0.
            C1(nz-k+2,ny-j+2,1) = 0.
            C2(k,j,1) = 0.
            C2(nz-k+2,ny-j+2,1) = 0.
            C3(k,j,1) = 0.
            C3(nz-k+2,ny-j+2,1) = 0.
          ENDIF
        END DO
      END DO
      DO i = 2,iend
        DO j = 1,ny
          DO k = 1,nz
            IF ((kk2(k,j,i).le.kup**2).and.(kk2(k,j,i).ge.kdn**2)) THEN
              dump = 1./sqrt(kk2(k,j,i))**4
              phase = 2*pi*randu(seed)
              C1(k,j,i) = 2*(COS(phase)+im*SIN(phase))*dump
              phase = 2*pi*randu(seed)
              C2(k,j,i) = 2*(COS(phase)+im*SIN(phase))*dump
              phase = 2*pi*randu(seed)
              C3(k,j,i) = 2*(COS(phase)+im*SIN(phase))*dump
            ELSE
              C1(k,j,i) = 0.
              C2(k,j,i) = 0.
              C3(k,j,i) = 0.
            ENDIF
          END DO
        END DO
      END DO
    ELSE
      DO i = ista,iend
        DO j = 1,ny
          DO k = 1,nz
            IF ((kk2(k,j,i).le.kup**2).and.(kk2(k,j,i).ge.kdn**2)) THEN
              dump = 1./sqrt(kk2(k,j,i))**4
              phase = 2*pi*randu(seed)
              C1(k,j,i) = 2*(COS(phase)+im*SIN(phase))*dump
              phase = 2*pi*randu(seed)
              C2(k,j,i) = 2*(COS(phase)+im*SIN(phase))*dump
              phase = 2*pi*randu(seed)
              C3(k,j,i) = 2*(COS(phase)+im*SIN(phase))*dump
            ELSE
              C1(k,j,i) = 0.
              C2(k,j,i) = 0.
              C3(k,j,i) = 0.
            ENDIF
          END DO
        END DO
      END DO
    ENDIF
    CALL rotor3(C2,C3,C4,1)
    CALL rotor3(C1,C3,C5,2)
    CALL rotor3(C1,C2,C6,3)
    CALL normalize(C4,C5,C6,f0,1,MPI_COMM_WORLD)
    ! C4,C5,C6 are the components of the first random vector (say v1)
    ! Repeating this procedure we create the second random vector (v2)
    IF (ista.eq.1) THEN
      C1(1,1,1) = 0.
      C2(1,1,1) = 0.
      C3(1,1,1) = 0.
      DO j = 2,ny/2+1
        IF ((kk2(1,j,1).le.kup**2).and.(kk2(1,j,1).ge.kdn**2)) THEN
          dump = 1./sqrt(kk2(1,j,1))**4
          phase = 2*pi*randu(seed)
          C1(1,j,1) = (COS(phase)+im*SIN(phase))*dump
          C1(1,ny-j+2,1) = conjg(C1(1,j,1))
          phase = 2*pi*randu(seed)
          C2(1,j,1) = (COS(phase)+im*SIN(phase))*dump
          C2(1,ny-j+2,1) = conjg(C2(1,j,1))
          phase = 2*pi*randu(seed)
          C3(1,j,1) = (COS(phase)+im*SIN(phase))*dump
          C3(1,ny-j+2,1) = conjg(C3(1,j,1))
        ELSE
          C1(1,j,1) = 0.
          C1(1,ny-j+2,1) = 0.
          C2(1,j,1) = 0.
          C2(1,ny-j+2,1) = 0.
          C3(1,j,1) = 0.
          C3(1,ny-j+2,1) = 0.
        ENDIF
      END DO
      DO k = 2,nz/2+1
        IF ((kk2(k,1,1).le.kup**2).and.(kk2(k,1,1).ge.kdn**2)) THEN
          dump = 1./sqrt(kk2(k,1,1))**4
          phase = 2*pi*randu(seed)
          C1(k,1,1) = (COS(phase)+im*SIN(phase))*dump
          C1(nz-k+2,1,1) = conjg(C1(k,1,1))
          phase = 2*pi*randu(seed)
          C2(k,1,1) = (COS(phase)+im*SIN(phase))*dump
          C2(nz-k+2,1,1) = conjg(C2(k,1,1))
          phase = 2*pi*randu(seed)
          C3(k,1,1) = (COS(phase)+im*SIN(phase))*dump
          C3(nz-k+2,1,1) = conjg(C3(k,1,1))
        ELSE
          C1(k,1,1) = 0.
          C1(nz-k+2,1,1) = 0.
          C2(k,1,1) = 0.
          C2(nz-k+2,1,1) = 0.
          C3(k,1,1) = 0.
          C3(nz-k+2,1,1) = 0.
        ENDIF
      END DO
      DO j = 2,ny
        DO k = 2,nz/2+1
          IF ((kk2(k,j,1).le.kup**2).and.(kk2(k,j,1).ge.kdn**2)) THEN
            dump = 1./sqrt(kk2(k,j,1))**4
            phase = 2*pi*randu(seed)
            C1(k,j,1) = (COS(phase)+im*SIN(phase))*dump
            C1(nz-k+2,ny-j+2,1) = conjg(C1(k,j,1))
            phase = 2*pi*randu(seed)
            C2(k,j,1) = (COS(phase)+im*SIN(phase))*dump
            C2(nz-k+2,ny-j+2,1) = conjg(C2(k,j,1))
            phase = 2*pi*randu(seed)
            C3(k,j,1) = (COS(phase)+im*SIN(phase))*dump
            C3(nz-k+2,ny-j+2,1) = conjg(C3(k,j,1))
          ELSE
            C1(k,j,1) = 0.
            C1(nz-k+2,ny-j+2,1) = 0.
            C2(k,j,1) = 0.
            C2(nz-k+2,ny-j+2,1) = 0.
            C3(k,j,1) = 0.
            C3(nz-k+2,ny-j+2,1) = 0.
          ENDIF
        END DO
      END DO
      DO i = 2,iend
        DO j = 1,ny
          DO k = 1,nz
            IF ((kk2(k,j,i).le.kup**2).and.(kk2(k,j,i).ge.kdn**2)) THEN
              dump = 1./sqrt(kk2(k,j,i))**4
              phase = 2*pi*randu(seed)
              C1(k,j,i) = 2*(COS(phase)+im*SIN(phase))*dump
              phase = 2*pi*randu(seed)
              C2(k,j,i) = 2*(COS(phase)+im*SIN(phase))*dump
              phase = 2*pi*randu(seed)
              C3(k,j,i) = 2*(COS(phase)+im*SIN(phase))*dump
            ELSE
              C1(k,j,i) = 0.
              C2(k,j,i) = 0.
              C3(k,j,i) = 0.
            ENDIF
          END DO
        END DO
      END DO
    ELSE
      DO i = ista,iend
        DO j = 1,ny
          DO k = 1,nz
            IF ((kk2(k,j,i).le.kup**2).and.(kk2(k,j,i).ge.kdn**2)) THEN
              dump = 1./sqrt(kk2(k,j,i))**4
              phase = 2*pi*randu(seed)
              C1(k,j,i) = 2*(COS(phase)+im*SIN(phase))*dump
              phase = 2*pi*randu(seed)
              C2(k,j,i) = 2*(COS(phase)+im*SIN(phase))*dump
              phase = 2*pi*randu(seed)
              C3(k,j,i) = 2*(COS(phase)+im*SIN(phase))*dump
            ELSE
              C1(k,j,i) = 0.
              C2(k,j,i) = 0.
              C3(k,j,i) = 0.
            ENDIF
          END DO
        END DO
      END DO
    ENDIF
    CALL rotor3(C2,C3,C7,1)
    CALL rotor3(C1,C3,C8,2)
    CALL rotor3(C1,C2,C3,3)
    CALL normalize(C7,C8,C3,f0,1,MPI_COMM_WORLD)
    ! So far, v1 = (C4,C5,C6) and v2 = (C7,C8,C3) are two random 
    ! normalized vectors in Fourier space. Correlating vectors 
    ! v1 and v2 we create the force components fx, fy, and fz.
    a1=SIN(alpha)
    a2=COS(alpha)
    ! fx: components y and z of u2 are calculated in order to get 
    ! the x component of curl(u2)
    DO i = ista,iend
      DO j = 1,ny
        DO k = 1,nz
          IF ((kk2(k,j,i).le.kup**2).and.(kk2(k,j,i).ge.kdn**2)) THEN
            C1(k,j,i) = a1*C5(k,j,i)+a2*C8(k,j,i)
            C2(k,j,i) = a1*C6(k,j,i)+a2*C3(k,j,i)
          ENDIF
        END DO
      END DO
    END DO
    CALL rotor3(C1,C2,C1,1)
    DO i = ista,iend
      DO j = 1,ny
        DO k = 1,nz
          IF ((kk2(k,j,i).le.kup**2).and.(kk2(k,j,i).ge.kdn**2)) THEN
            dump = 1./sqrt(kk2(k,j,i))
            state(solver%VELOCITY  )%ccomp(k,j,i) = a2*C4(k,j,i) +  &
                        a1*C7(k,j,i) + C1(k,j,i)*dump
          ENDIF
        END DO
      END DO
    END DO
    ! fy: components x and z of u2 are calculated in order to get 
    ! the y component of curl(u2)
    DO i = ista,iend
      DO j = 1,ny
        DO k = 1,nz
          IF ((kk2(k,j,i).le.kup**2).and.(kk2(k,j,i).ge.kdn**2)) THEN
            C1(k,j,i) = a1*C4(k,j,i)+a2*C7(k,j,i)
            C2(k,j,i) = a1*C6(k,j,i)+a2*C3(k,j,i)
          ENDIF
        END DO
      END DO
    END DO
    CALL rotor3(C1,C2,C1,2)
    DO i = ista,iend
      DO j = 1,ny
        DO k = 1,nz
          IF ((kk2(k,j,i).le.kup**2).and.(kk2(k,j,i).ge.kdn**2)) THEN
            dump = 1./sqrt(kk2(k,j,i))
            state(solver%VELOCITY+1)%ccomp(k,j,i) = a2*C5(k,j,i) +  &
                        a1*C8(k,j,i) + C1(k,j,i)*dump
          ENDIF
        END DO
      END DO
    END DO
    ! fz: components x and y of u2 are calculated in order to get
    ! the z component of curl(u2)
    DO i = ista,iend
      DO j = 1,ny
        DO k = 1,nz
          IF ((kk2(k,j,i).le.kup**2).and.(kk2(k,j,i).ge.kdn**2)) THEN
            C1(k,j,i) = a1*C4(k,j,i)+a2*C7(k,j,i)
            C2(k,j,i) = a1*C5(k,j,i)+a2*C8(k,j,i)
          ENDIF
        END DO
      END DO
    END DO
    CALL rotor3(C1,C2,C1,3)
    DO i = ista,iend
      DO j = 1,ny
        DO k = 1,nz
          IF ((kk2(k,j,i).le.kup**2).and.(kk2(k,j,i).ge.kdn**2)) THEN
            dump = 1./sqrt(kk2(k,j,i))
            state(solver%VELOCITY+2)%ccomp(k,j,i) = a2*C6(k,j,i) +  &
                        a1*C3(k,j,i) + C1(k,j,i)*dump
          ENDIF
        END DO
      END DO
    END DO
    call normalize(state(solver%VELOCITY  )%ccomp, &
                   state(solver%VELOCITY+1)%ccomp, &
                   state(solver%VELOCITY+2)%ccomp, &
                   f0,1,MPI_COMM_WORLD)
    class default
      error stop "This solver does not support velocity forcing"
    end select
    call solver%workspace_%free_complex_tmp(C1)
    call solver%workspace_%free_complex_tmp(C2)
    call solver%workspace_%free_complex_tmp(C3)
    call solver%workspace_%free_complex_tmp(C4)
    call solver%workspace_%free_complex_tmp(C5)
    call solver%workspace_%free_complex_tmp(C6)
    call solver%workspace_%free_complex_tmp(C7)
    call solver%workspace_%free_complex_tmp(C8)
  end subroutine init_randomfv

  
  ! ===================================================================
  ! Update methods
  ! ===================================================================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Instantaneous phase shift
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine update_shiftfv(this, force, solver, state)
    use pseudospec_fluid
    use equationbase_mod
    use commtypes
    use random
    use status
    implicit none
    
    class(shiftupdt_fv),   intent(inout) :: this
    class   (forceBase),   intent   (in) :: force
    class(EquationBase),   intent   (in) :: solver
    type   (GStateComp),   intent(inout) :: state(:)
    complex(kind=GP)                     :: cdump
    real(kind=GP)                        :: phase

    select type (solver)
    class is (VelocityBase)
      if (timef.eq.fstep) then
        if (myrank.eq.0) phase = 2*pi*randu(seed)
        call MPI_BCAST(phase,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
        cdump = COS(phase)+im*SIN(phase)
        call phaseshift(state(solver%VELOCITY  )%ccomp,cdump)
        call phaseshift(state(solver%VELOCITY+1)%ccomp,cdump)
        call phaseshift(state(solver%VELOCITY+2)%ccomp,cdump)
      endif
    class default
      error stop "This solver does not support velocity forcing"
    end select
  end subroutine update_shiftfv


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Slowly evolving random phases. This update method
  !! requires randomnes in the generation of the forcing,
  !! i.e., it doesn't work with deterministic forcing as ABC.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine update_shufflefv(this, force, solver, state)
    use equationbase_mod
    use status
    use filefmt
    use fft
    use commtypes
    implicit none
    
    class(shuffleupdt_fv), intent(inout)             :: this
    class     (forceBase), intent   (in)             :: force
    class  (EquationBase), intent   (in)             :: solver
    type     (GStateComp), intent(inout)             :: state(:)
    complex     (kind=GP), pointer, dimension(:,:,:) :: C1,C2,C3
    real        (kind=GP), pointer, dimension(:,:,:) :: R1,R2,R3
    complex(kind=GP)                                 :: cdump
    real(kind=GP)                                    :: phase,rmp
    integer                                          :: i,j,k
    logical                                          :: bret

    ! If not initialized, initialize the method
    if ( .not. this%binit_ ) then
      call solver%workspace_%get_complex_tmp(this%fxold_,bret)
      call solver%workspace_%get_complex_tmp(this%fyold_,bret)
      call solver%workspace_%get_complex_tmp(this%fzold_,bret)
      call solver%workspace_%get_complex_tmp(this%fxnew_,bret)
      call solver%workspace_%get_complex_tmp(this%fynew_,bret)
      call solver%workspace_%get_complex_tmp(this%fznew_,bret)
      if ( stat .eq. 0 ) then ! Generate new forcing
         timef = fstep
      else                    ! Read previous forcing state
        write(ext, fmtext) tind
        call solver%workspace_%get_real_tmp(R1,bret)
        call solver%workspace_%get_real_tmp(R2,bret)
        call solver%workspace_%get_real_tmp(R3,bret)
        call io_read(1,solver%idir_,'fvxold',ext,solver%planio_,R1)
        call io_read(1,solver%idir_,'fvyold',ext,solver%planio_,R2)
        call io_read(1,solver%idir_,'fvzold',ext,solver%planio_,R3)
        call fftp3d_real_to_complex(planrc,R1,this%fxold_,MPI_COMM_WORLD)
        call fftp3d_real_to_complex(planrc,R2,this%fyold_,MPI_COMM_WORLD)
        call fftp3d_real_to_complex(planrc,R3,this%fzold_,MPI_COMM_WORLD)
        call io_read(1,solver%idir_,'fvxnew',ext,solver%planio_,R1)
        call io_read(1,solver%idir_,'fvynew',ext,solver%planio_,R2)
        call io_read(1,solver%idir_,'fvznew',ext,solver%planio_,R3)
        call fftp3d_real_to_complex(planrc,R1,this%fxnew_,MPI_COMM_WORLD)
        call fftp3d_real_to_complex(planrc,R2,this%fynew_,MPI_COMM_WORLD)
        call fftp3d_real_to_complex(planrc,R3,this%fznew_,MPI_COMM_WORLD)
        call solver%workspace_%free_real_tmp(R1)
        call solver%workspace_%free_real_tmp(R2)
        call solver%workspace_%free_real_tmp(R3)
      endif    
      this%binit_ = .TRUE.
    endif
    select type (solver)
    class is (VelocityBase)
      ! Generate new forcing states when the correlation time is reached
      if (timef.eq.fstep) then
        do i = ista,iend        ! Keeps a copy of the last forcing state
          do j = 1,ny
            do k = 1,nz
              this%fxold_(k,j,i) = state(solver%VELOCITY  )%ccomp(k,j,i)
              this%fyold_(k,j,i) = state(solver%VELOCITY+1)%ccomp(k,j,i)
              this%fzold_(k,j,i) = state(solver%VELOCITY+2)%ccomp(k,j,i)
            end do
          end do
        end do
        call force%init_GForce(solver,state)     ! Creates a new forcing
        do i = ista,iend                ! Copies the new forcing to fnew
          do j = 1,ny
            do k = 1,nz
              this%fxnew_(k,j,i) = state(solver%VELOCITY  )%ccomp(k,j,i)
              this%fynew_(k,j,i) = state(solver%VELOCITY+1)%ccomp(k,j,i)
              this%fznew_(k,j,i) = state(solver%VELOCITY+2)%ccomp(k,j,i)
            end do
          end do
        end do
      endif
      ! Slowly updates the forcing at every time step
      rmp = float(timef+1)/float(fstep)
      do i = ista,iend
        do j = 1,ny
          do k = 1,nz
            state(solver%VELOCITY  )%ccomp(k,j,i) = &
                 (1-rmp)*this%fxold_(k,j,i)+rmp*this%fxnew_(k,j,i)
            state(solver%VELOCITY+1)%ccomp(k,j,i) = &
                 (1-rmp)*this%fyold_(k,j,i)+rmp*this%fynew_(k,j,i)
            state(solver%VELOCITY+2)%ccomp(k,j,i) = &
                 (1-rmp)*this%fzold_(k,j,i)+rmp*this%fznew_(k,j,i)
          end do
        end do
      end do
    class default
      error stop "This solver does not support velocity forcing"
    end select
    ! If we just wrote binary files, we save the force state
    if ((timet.eq.0).and.(bench.eq.0)) then
      write(ext, fmtext) tind
      rmp = 1.0_GP/(real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))
      call solver%workspace_%get_complex_tmp(C1,bret)
      call solver%workspace_%get_complex_tmp(C2,bret)
      call solver%workspace_%get_complex_tmp(C3,bret)
      call solver%workspace_%get_real_tmp(R1,bret)
      call solver%workspace_%get_real_tmp(R2,bret)
      call solver%workspace_%get_real_tmp(R3,bret)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      do i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
        do j = 1,ny
          do k = 1,nz
            C1(k,j,i) = this%fxold_(k,j,i)*rmp
            C2(k,j,i) = this%fyold_(k,j,i)*rmp
            C3(k,j,i) = this%fzold_(k,j,i)*rmp
          end do
        end do
      end do
      call fftp3d_complex_to_real(plancr,C1,R1,MPI_COMM_WORLD)
      call fftp3d_complex_to_real(plancr,C2,R2,MPI_COMM_WORLD)
      call fftp3d_complex_to_real(plancr,C3,R3,MPI_COMM_WORLD)
      call io_write(1,solver%odir_,'fvxold',ext,solver%planio_,R1)
      call io_write(1,solver%odir_,'fvyold',ext,solver%planio_,R2)
      call io_write(1,solver%odir_,'fvzold',ext,solver%planio_,R3)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      do i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
        do j = 1,ny
          do k = 1,nz
            C1(k,j,i) = this%fxnew_(k,j,i)*rmp
            C2(k,j,i) = this%fynew_(k,j,i)*rmp
            C3(k,j,i) = this%fznew_(k,j,i)*rmp
          end do
        end do
      end do
      call fftp3d_complex_to_real(plancr,C1,R1,MPI_COMM_WORLD)
      call fftp3d_complex_to_real(plancr,C2,R2,MPI_COMM_WORLD)
      call fftp3d_complex_to_real(plancr,C3,R3,MPI_COMM_WORLD)
      call io_write(1,solver%odir_,'fvxnew',ext,solver%planio_,R1)
      call io_write(1,solver%odir_,'fvynew',ext,solver%planio_,R2)
      call io_write(1,solver%odir_,'fvznew',ext,solver%planio_,R3)
      call solver%workspace_%free_complex_tmp(C1)
      call solver%workspace_%free_complex_tmp(C2)
      call solver%workspace_%free_complex_tmp(C3)
      call solver%workspace_%free_real_tmp(R1)
      call solver%workspace_%free_real_tmp(R2)
      call solver%workspace_%free_real_tmp(R3)
    endif
  end subroutine update_shufflefv
  
end module force_velocity
