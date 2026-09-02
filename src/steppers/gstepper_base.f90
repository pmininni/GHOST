! ===================================================================
! NAME       : gstepper_base.fpp
! DESCRIPTION: Forms base class for all time steppers
! DATE       : 2/8/26 (DLR)
! ===================================================================

module gstepperbase_mod
  use equationbase_mod, only: EquationBase, GWorkspace
  use particlebase_mod, only: ParticleBase
  use gpstate_mod
  use gstate_mod
  implicit none

  ! ================= Stepper traits ===================================
  type, public  :: GStepperTraits
    integer            :: itype       ! =1 Butcher type 
    integer            :: norder      ! Order
    integer            :: nstage      ! Number of stages
    integer            :: nstate  = 0 ! Number of field state components
    integer            :: npstate = 0 ! Number of particles state components
    logical            :: dopart  = .false. ! Use particle interfaces?   
    character(len=128) :: sname       ! Stepper name
  end type GStepperTraits

  ! ================= Base class for all steppers =======================
  ! Define an abstract base class
  type, abstract :: GStepperBase
      type   (GWorkspace), pointer   :: workspace_ => null()
      class(EquationBase), pointer   :: solver_    => null()
      class(ParticleBase), pointer   :: psolver_   => null()
      character (len=128)            :: infile_    ! config file name
    contains
      procedure(GStepper_ctor_interface), deferred :: GStepper_ctor ! Constructor
      procedure         (init_interface), deferred :: init          ! init method
      procedure         (step_interface), deferred :: step          ! fluid step method
      procedure        (pstep_interface), deferred :: pstep         ! part  step method
      procedure        (cstep_interface), deferred :: cstep         ! part+field step method
      generic              , public :: gstep => step, pstep, cstep  ! general method
  end type GStepperBase 
  
  abstract interface
    subroutine GStepper_ctor_interface(this, traits, workspace, solver, psolver)
      import :: GStepperBase, GStepperTraits, GWorkspace, EquationBase, ParticleBase
      class (GStepperBase), intent(inout)                   :: this
      type(GStepperTraits), intent(inout)                   :: traits
      type    (GWorkspace), intent(inout), target           :: workspace
      class (EquationBase), intent   (in), target           :: solver
      class (ParticleBase), intent   (in), target, optional :: psolver
    end subroutine GStepper_ctor_interface

    subroutine init_interface(this, traits) 
      import :: GStepperBase, GStepperTraits
      class (GStepperBase), intent(inout) :: this
      type(GStepperTraits), intent   (in) :: traits
    end subroutine init_interface

    subroutine step_interface(this, time, uin, uf, dt, uout)
      import :: GStepperBase, GStateComp, GP
      class(GStepperBase), intent(inout) :: this
      real      (kind=GP), intent   (in) :: time, dt
      type   (GStateComp), intent(inout), target :: uin(:),uf(:)
      type   (GStateComp), intent(inout), target :: uout(:) 
    end subroutine step_interface

    subroutine pstep_interface(this, time, uin, upin, dt, upout)
      import :: GStepperBase, GStateComp, GPStateComp, GP
      class(GStepperBase), intent(inout)                      :: this
      real      (kind=GP), intent   (in)                      :: time, dt
      type   (GStateComp), intent(inout)                      :: uin(:)
      type  (GPStateComp), intent(inout), target, allocatable :: upin(:), upout(:)
    end subroutine pstep_interface

    subroutine cstep_interface(this, time, uin, upin, uf, dt, uout, upout)
      import :: GStepperBase, GStateComp, GPStateComp, GP
      class(GStepperBase), intent(inout)                      :: this
      real      (kind=GP), intent   (in)                      :: time, dt
      type   (GStateComp), intent(inout)                      :: uin(:),uf(:), uout(:)
      type  (GPStateComp), intent(inout), target, allocatable :: upin(:), upout(:)
    end subroutine cstep_interface
  end interface

contains

  ! ===================================================================
  ! Concrete methods inherited by all steppers
  ! ===================================================================
  
end module gstepperbase_mod
