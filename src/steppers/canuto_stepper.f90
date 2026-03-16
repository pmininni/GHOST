! =====================================================================
! NAME       : canuto_stepper.f90
!
!              Performs time stepping using the modified RK method 
!              of Canuto et al., Spectral Methods in Fluid
!              Dynamics. This is a low storage time stepper
!              for a specifiable order, though it is strictly
!              speaking of full truncation order only for
!              norder = 2, or for linear PDEs. Explicit number
!              of stages are not required. While norder > 2 may
!              not yield a truncation of that order in nonlinear
!              PDEs, it can still provide a benefit.
!
! INPUT FILE : Stepper looks for a "&stepper" namelist with:
!                sname   : 'TRADITIONAL' to use this time stepper
!                itype   : Ignored by this solver. The stepper type:
!                    =1 => Butcher type (norder = nstage, this solver)
!                    =2 => Mixed type (norder != nstage; not available)
!                    =3 => SSP type (strong stability preserving; not available)
!                norder  : Stepper order (>=1)
!                nstage  : No. stepper stages (ignored by this solver)
!
! DATE       : 2/8/26 (DLR)
! =====================================================================

module canuto_stepper_mod
  use gstepperbase_mod
  implicit none

  ! ================= Global parameters ===============================

  ! ================= Stepper traits ===================================

  ! ================= Stepper ==========================================
  ! Define class:
  type, extends(GStepperBase) :: CanutoStepper
    ! Member data (extends GStepperBase member data)
    type(GStepperTraits)      :: traits_                     ! GStepper traits
  contains
    procedure, public :: init                                ! initialize
    procedure, public :: step          => step_impl          ! take one timestep
    procedure, public :: pstep         => pstep_impl         ! take one part+field timestep
    procedure, public :: GStepper_ctor => CanutoStepper_ctor
    final             :: CanutoStepper_dtor
  end type CanutoStepper

contains

  ! ===================================================================
  ! Stepper-specific methods
  ! ===================================================================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Constructor
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine CanutoStepper_ctor(this, traits, workspace, solver, psolver)
    class(CanutoStepper), intent(inout)                   :: this
    type(GStepperTraits), intent(inout)                   :: traits
    type    (GWorkspace), intent(inout), target           :: workspace
    class (EquationBase), intent   (in), target           :: solver
    class (ParticleBase), intent   (in), target, optional :: psolver
    this%workspace_ => workspace
    this%solver_    => solver
    if (present(psolver)) then
      this%psolver_ => psolver
    else
      this%psolver_ => null()
    endif
    if (.not. associated(this%workspace_)) then
      stop 'CanutoStepper::CanutoStepper_ctor: Worskpace not associated'
    endif
    call this%init(traits)
  end subroutine CanutoStepper_ctor


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Destructor
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine CanutoStepper_dtor(this) 
    type(CanutoStepper), intent(inout) :: this
    if (associated(this%workspace_)) nullify(this%workspace_)
    if (associated(this%solver_))    nullify(this%solver_)
    if (associated(this%psolver_))   nullify(this%psolver_)
  end subroutine CanutoStepper_dtor


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Subroutine to initialize the stepper
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init(this, traits)
    class(CanutoStepper), intent (inout) :: this
    type(GStepperTraits), intent    (in) :: traits
    this%traits_ = traits;
  end subroutine init


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Implementation function to take one step
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine step_impl(this, time, uin, uf, dt, uout)
!$  use threads
    implicit none

    class(CanutoStepper), intent   (in) :: this
    type    (GStateComp), intent(inout) :: uin(:), uf(:), uout(:)
    real       (kind=GP), intent   (in) :: time, dt
    real       (kind=GP)                :: eff_dt
    integer                             :: i,j,k,o,state_size,ic
    logical                             :: bret
       
    if ( size(uin)  .ne. this%traits_%nstate &
     .or.size(uout) .ne. this%traits_%nstate  ) then
      stop 'CanutoStepper::step: Inconsistent input state'
    endif

    do o = this%traits_%norder,1,-1
      eff_dt = dt/real(o,kind=GP)
      call this%solver_%dudt(time, uout, uf, eff_dt, uout)
      do ic = 1,this%traits_%nstate
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
        do i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
          do j = 1,ny
            do k = 1,nz
              uout(ic)%ccomp(k,j,i) = uin(ic)%ccomp(k,j,i) + &
                              eff_dt*uout(ic)%ccomp(k,j,i)
            end do
          end do
        end do
      end do
    end do
  end subroutine step_impl


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Implementation function to take one step
  !! for particles + fields
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine pstep_impl(this, time, uin, upin, uf, dt, uout, upout)
!$  use threads
    use gpstate_mod
    implicit none

    class(CanutoStepper), intent   (in) :: this
    type    (GStateComp), intent(inout) :: uin (:), uf(:), uout(:)
    type   (GPStateComp), intent(inout) :: upin(:), upout(:)
    real       (kind=GP), intent   (in) :: time, dt
    real       (kind=GP)                :: eff_dt
    integer                             :: i,j,k,o,state_size,ic,ip,nparts
    logical                             :: bret
       
    if ( size(uin) .ne. this%traits_%nstate &
     .or.size(uout) .ne. this%traits_%nstate  ) then
      stop 'CanutoStepper::step: Inconsistent input state'
    endif
    
!!$    if ( .not. associated(this%pcallback_) ) then
!!$      stop 'CanutoStepper::step: RHS pcallback function not set'
!!$    endif
!!$
!!$    nparts = size(upin(1)%rcomp)
!!$
!!$    do o = this%traits_%norder,1,-1
!!$
!!$      eff_dt = dt/real(o,kind=GP)
!!$      call this%pcallback_(time, uout, upin, uf, eff_dt, uout, upout)
!!$
!!$      ! Update fields:
!!$      do ic = 1,this%traits_%nstate  ! for each state comp
!!$!$omp parallel do if (iend-ista.ge.nth) private (j,k)
!!$        do i = ista,iend
!!$!$omp parallel do if (iend-ista.lt.nth) private (k)
!!$          do j = 1,ny
!!$            do k = 1,nz
!!$              uout(ic)%ccomp(k,j,i) = uin(ic)%ccomp(k,j,i) + &
!!$                              eff_dt*uout(ic)%ccomp(k,j,i)
!!$            end do
!!$          end do
!!$        end do
!!$      end do
!!$
!!$      ! Update particles:
!!$      do ip = 1,this%traits_%npstate ! for each state comp
!!$        do k = 1, nparts
!!$          upout(ip)%rcomp(k) = upin(ip)%rcomp(k) + &
!!$                               eff_dt*upout(ip)%rcomp(k)
!!$        enddo
!!$      end do
!!$
!!$    end do ! end, o-loop
  end subroutine pstep_impl

end module canuto_stepper_mod
