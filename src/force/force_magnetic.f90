! =====================================================================
! NAME       : force_magnetic.f90
! DESCRIPTION: Forcing methods for all MagneticBase solver classes.
!              Each force needs an initial set up, and an update.
!              These forcing methods can be used to initialize and
!              update forcings on other solver classes, using
!              composition of operators, as long as no correlation
!              between the forces is needed. For fv-fb correlated
!              cases, this file can provide methods to initialize
!              both forces simultaneously.
!
! Forces avaliable:
!   null_fb     : Null vector potential
!   random_fb   : Random (Pouquet-Patterson) electromotive forcing
!
! Update methods available:
!   constant_fb : Constant forcing (no method)
!   shift_fb    : Instantaneous phase shift (uncorrelated from fv)
!   shuffle_fb  : Slowly evolving random phases
!
! DATE       : 01/24/26 (PDM)
! =====================================================================

module force_magnetic
  USE forcebase_mod

  IMPLICIT NONE

  ! ================= Forcing functions supported =====================
  type, extends(forceBase) :: forceNull_fb
    contains
      procedure :: init_GForce => init_nullfb
  end type forceNull_fb 
  type, extends(forceBase) :: forceRandom_fb
    contains
      procedure :: init_GForce => init_randomfb
  end type forceRandom_fb 
! type, extends(forceBase) :: forceUserdef_fb
!   contains
!     procedure :: init_GForce => init_userdefb
! end type forceUserdef_fb

  ! ================= Update methods supported =======================
  type, extends(forceUpdt) :: shiftupdt_fb
    contains
      procedure ::   update_GForce => update_shiftfb
  end type shiftupdt_fb 
  type, extends(forceUpdt) :: shuffleupdt_fb
    contains
      procedure ::   update_GForce => update_shufflefb
  end type shuffleupdt_fb
! type, extends(forceUpdt) :: userupdt_fb
!   contains
!     procedure :: init_GForce => init_userupdtfb
! end type userupdt_fb

CONTAINS

  ! ===================================================================
  ! Forcing functions
  ! ===================================================================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Null electromotive forcing
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_nullfb(this,solver,state)
    use gstate_mod
    use mhd_mod
    use grid
    use mpivars
!$  use threads
    implicit none
    
    class(forceNull_fb), intent   (in) :: this
    class(EquationBase), intent   (in) :: solver
    type   (GStateComp), intent(inout) :: state(:)
    integer                            :: i,j,k

    select type (solver)
    class is (MagneticBase)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
        DO j = 1,ny
          DO k = 1,nz
            state(solver%MAGNETIC  )%ccomp(k,j,i) = 0.0_GP
            state(solver%MAGNETIC+1)%ccomp(k,j,i) = 0.0_GP
            state(solver%MAGNETIC+2)%ccomp(k,j,i) = 0.0_GP
          END DO
        END DO
      END DO
    class default
      error stop "This solver does not support electromotive forcing"
    end select
  end subroutine init_nullfb

  
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
  !   corr : correlation with the velocity forcing
  subroutine init_randomfb(this,solver,state)
    use gstate_mod
    use mhd_mod
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

    class(forceRandom_fb), intent   (in)          :: this
    class  (EquationBase), intent   (in)          :: solver
    type     (GStateComp), intent(inout)          :: state(:)
    complex  (kind=GP), pointer, dimension(:,:,:) :: C1,C2,C3,C4
    complex  (kind=GP), pointer, dimension(:,:,:) :: C5,C6,C7,C8
    real     (kind=GP)                            :: f0,kdn,kup
    real     (kind=GP)                            :: alpha,a1,a2
    real     (kind=GP)                            :: dump,phase
    real     (kind=GP)                            :: corr,rmp
    integer                                       :: i,j,k
    logical                                       :: bret

    namelist/ random_fb / f0,kdn,kup,alpha,corr
    CALL solver%workspace_%get_complex_tmp(C1,bret)
    CALL solver%workspace_%get_complex_tmp(C2,bret)
    CALL solver%workspace_%get_complex_tmp(C3,bret)
    CALL solver%workspace_%get_complex_tmp(C4,bret)
    CALL solver%workspace_%get_complex_tmp(C5,bret)
    CALL solver%workspace_%get_complex_tmp(C6,bret)
    CALL solver%workspace_%get_complex_tmp(C7,bret)
    CALL solver%workspace_%get_complex_tmp(C8,bret)
    select type (solver)
    class is (MagneticBase)
    ! Read parameters from a namelist in the input file
    alpha = 0.0_GP
    corr  = 0.0_GP
    if ( myrank .eq. 0 ) then
      open(1,file=solver%infile_,status='unknown',form="formatted")
      read(1,NML=random_fb)
      close(1)
    endif
    call mpi_bcast(f0   ,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(kup  ,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(kdn  ,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(alpha,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(corr ,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
    ! Generate the first random velocity field
    IF (ista.eq.1) THEN
      C1(1,1,1) = 0.
      C2(1,1,1) = 0.
      C3(1,1,1) = 0.
      DO j = 2,ny/2+1
        IF ((kk2(1,j,1).le.kup**2).and.(kk2(1,j,1).ge.kdn**2)) THEN
          dump = 1./sqrt(kk2(1,j,1))**5
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
          dump = 1./sqrt(kk2(k,1,1))**5
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
            dump = 1./sqrt(kk2(k,j,1))**5
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
              dump = 1./sqrt(kk2(k,j,i))**5
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
              dump = 1./sqrt(kk2(k,j,i))**5
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
    CALL normalize(C4,C5,C6,f0,0,MPI_COMM_WORLD)
    ! C4,C5,C6 are the components of the first random vector (say v1)
    ! Repeating this procedure we create the second random vector (v2)
    IF (ista.eq.1) THEN
      C1(1,1,1) = 0.
      C2(1,1,1) = 0.
      C3(1,1,1) = 0.
      DO j = 2,ny/2+1
        IF ((kk2(1,j,1).le.kup**2).and.(kk2(1,j,1).ge.kdn**2)) THEN
          dump = 1./sqrt(kk2(1,j,1))**5
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
          dump = 1./sqrt(kk2(k,1,1))**5
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
            dump = 1./sqrt(kk2(k,j,1))**5
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
              dump = 1./sqrt(kk2(k,j,i))**5
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
              dump = 1./sqrt(kk2(k,j,i))**5
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
    CALL normalize(C7,C8,C3,f0,0,MPI_COMM_WORLD)
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
            state(solver%MAGNETIC  )%ccomp(k,j,i) = a2*C4(k,j,i) +  &
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
            state(solver%MAGNETIC+1)%ccomp(k,j,i) = a2*C5(k,j,i) +  &
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
            state(solver%MAGNETIC+2)%ccomp(k,j,i) = a2*C6(k,j,i) +  &
                        a1*C3(k,j,i) + C1(k,j,i)*dump
          ENDIF
        END DO
      END DO
    END DO
    call normalize(state(solver%MAGNETIC  )%ccomp, &
                   state(solver%MAGNETIC+1)%ccomp, &
                   state(solver%MAGNETIC+2)%ccomp, &
                   f0,0,MPI_COMM_WORLD)
    if ( corr .gt. tiny ) then   ! Correlate with the velocity forcing
      DO i = ista,iend
        DO j = 1,ny
          DO k = 1,nz
            C1(k,j,i) = state(solver%VELOCITY  )%ccomp(k,j,i)
            C2(k,j,i) = state(solver%VELOCITY+1)%ccomp(k,j,i)
            C3(k,j,i) = state(solver%VELOCITY+2)%ccomp(k,j,i)
          END DO
        END DO
      END DO
      call normalize(C1,C2,C3,f0,1,MPI_COMM_WORLD)
      CALL rotor3(C2,C3,C4,1)
      CALL rotor3(C1,C3,C5,2)
      CALL rotor3(C1,C2,C6,3)
      rmp = sqrt(1-corr**2)
      DO i = ista,iend
        DO j = 1,ny
          DO k = 1,nz
            IF ((kk2(k,j,i).le.kup**2).and.(kk2(k,j,i).ge.kdn**2)) THEN
              dump = 1.0_GP/kk2(k,j,i)
              state(solver%MAGNETIC  )%ccomp(k,j,i) =          &
                   rmp*state(solver%MAGNETIC  )%ccomp(k,j,i) + &
                   corr*dump*C4(k,j,i)
              state(solver%MAGNETIC+1)%ccomp(k,j,i) =          &
                   rmp*state(solver%MAGNETIC+1)%ccomp(k,j,i) + &
                   corr*dump*C5(k,j,i)
              state(solver%MAGNETIC+2)%ccomp(k,j,i) =          &
                   rmp*state(solver%MAGNETIC+2)%ccomp(k,j,i) + &
                   corr*dump*C6(k,j,i)
            ENDIF
          END DO
        END DO
      END DO
      call normalize(state(solver%MAGNETIC  )%ccomp, &
                     state(solver%MAGNETIC+1)%ccomp, &
                     state(solver%MAGNETIC+2)%ccomp, &
                     f0,0,MPI_COMM_WORLD)
    endif
    class default
      error stop "This solver does not support electromotive forcing"
    end select
    call solver%workspace_%free_complex_tmp(C1)
    call solver%workspace_%free_complex_tmp(C2)
    call solver%workspace_%free_complex_tmp(C3)
    call solver%workspace_%free_complex_tmp(C4)
    call solver%workspace_%free_complex_tmp(C5)
    call solver%workspace_%free_complex_tmp(C6)
    call solver%workspace_%free_complex_tmp(C7)
    call solver%workspace_%free_complex_tmp(C8)
  end subroutine init_randomfb

  
  ! ===================================================================
  ! Update methods
  ! ===================================================================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Instantaneous phase shift
  !! This method destroys any cross correlation between forces
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine update_shiftfb(this, force, solver, state)
    use pseudospec_fluid
    use equationbase_mod
    use commtypes
    use random
    use status
    implicit none
    
    class(shiftupdt_fb),   intent(inout) :: this
    class   (forceBase),   intent   (in) :: force
    class(EquationBase),   intent   (in) :: solver
    type   (GStateComp),   intent(inout) :: state(:)
    complex(kind=GP)                     :: cdump
    real(kind=GP)                        :: phase

    select type (solver)
    class is (MagneticBase)
      if (timef.eq.fstep) then
        if (myrank.eq.0) phase = 2*pi*randu(seed)
        call MPI_BCAST(phase,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
        cdump = COS(phase)+im*SIN(phase)
        call phaseshift(state(solver%MAGNETIC  )%ccomp,cdump)
        call phaseshift(state(solver%MAGNETIC+1)%ccomp,cdump)
        call phaseshift(state(solver%MAGNETIC+2)%ccomp,cdump)
      endif
    class default
      error stop "This solver does not support electromotive forcing"
    end select
  end subroutine update_shiftfb


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Slowly evolving random phases. This update method
  !! requires randomnes in the generation of the forcing,
  !! but preserves cross correlations between forces.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine update_shufflefb(this, force, solver, state)
    use equationbase_mod
    use status
    use filefmt
    use fft
    use commtypes
    implicit none
    
    class(shuffleupdt_fb), intent(inout)             :: this
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
        call io_read(1,solver%idir_,'fbxold',ext,solver%planio_,R1)
        call io_read(1,solver%idir_,'fbyold',ext,solver%planio_,R2)
        call io_read(1,solver%idir_,'fbzold',ext,solver%planio_,R3)
        call fftp3d_real_to_complex(planrc,R1,this%fxold_,MPI_COMM_WORLD)
        call fftp3d_real_to_complex(planrc,R2,this%fyold_,MPI_COMM_WORLD)
        call fftp3d_real_to_complex(planrc,R3,this%fzold_,MPI_COMM_WORLD)
        call io_read(1,solver%idir_,'fbxnew',ext,solver%planio_,R1)
        call io_read(1,solver%idir_,'fbynew',ext,solver%planio_,R2)
        call io_read(1,solver%idir_,'fbznew',ext,solver%planio_,R3)
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
    class is (MagneticBase)
      ! Generate new forcing states when the correlation time is reached
      if (timef.eq.fstep) then
        do i = ista,iend        ! Keeps a copy of the last forcing state
          do j = 1,ny
            do k = 1,nz
              this%fxold_(k,j,i) = state(solver%MAGNETIC  )%ccomp(k,j,i)
              this%fyold_(k,j,i) = state(solver%MAGNETIC+1)%ccomp(k,j,i)
              this%fzold_(k,j,i) = state(solver%MAGNETIC+2)%ccomp(k,j,i)
            end do
          end do
        end do
        call force%init_GForce(solver,state)     ! Creates a new forcing
        do i = ista,iend                ! Copies the new forcing to fnew
          do j = 1,ny
            do k = 1,nz
              this%fxnew_(k,j,i) = state(solver%MAGNETIC  )%ccomp(k,j,i)
              this%fynew_(k,j,i) = state(solver%MAGNETIC+1)%ccomp(k,j,i)
              this%fznew_(k,j,i) = state(solver%MAGNETIC+2)%ccomp(k,j,i)
            end do
          end do
        end do
      endif
      ! Slowly updates the forcing at every time step
      rmp = float(timef+1)/float(fstep)
      do i = ista,iend
        do j = 1,ny
          do k = 1,nz
            state(solver%MAGNETIC  )%ccomp(k,j,i) = &
                 (1-rmp)*this%fxold_(k,j,i)+rmp*this%fxnew_(k,j,i)
            state(solver%MAGNETIC+1)%ccomp(k,j,i) = &
                 (1-rmp)*this%fyold_(k,j,i)+rmp*this%fynew_(k,j,i)
            state(solver%MAGNETIC+2)%ccomp(k,j,i) = &
                 (1-rmp)*this%fzold_(k,j,i)+rmp*this%fznew_(k,j,i)
          end do
        end do
      end do
    class default
      error stop "This solver does not support electromotive forcing"
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
      call io_write(1,solver%odir_,'fbxold',ext,solver%planio_,R1)
      call io_write(1,solver%odir_,'fbyold',ext,solver%planio_,R2)
      call io_write(1,solver%odir_,'fbzold',ext,solver%planio_,R3)
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
      call io_write(1,solver%odir_,'fbxnew',ext,solver%planio_,R1)
      call io_write(1,solver%odir_,'fbynew',ext,solver%planio_,R2)
      call io_write(1,solver%odir_,'fbznew',ext,solver%planio_,R3)
      call solver%workspace_%free_complex_tmp(C1)
      call solver%workspace_%free_complex_tmp(C2)
      call solver%workspace_%free_complex_tmp(C3)
      call solver%workspace_%free_real_tmp(R1)
      call solver%workspace_%free_real_tmp(R2)
      call solver%workspace_%free_real_tmp(R3)
    endif
  end subroutine update_shufflefb
  
end module force_magnetic
