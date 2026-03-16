! ===================================================================
! NAME       : gstepper_base.fpp
! DESCRIPTION: Forms base class for all time steppers
! DATE       : 2/8/26 (DLR)
! ===================================================================

module gstepperbase_mod
  use equationbase_mod
  use particlebase_mod
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
  end type

  ! ================= Base class for all steppers =======================
  ! Define an abstract base class
  type, abstract :: GStepperBase
      type   (GWorkspace), pointer   :: workspace_
      class(EquationBase), pointer   :: solver_
      class(ParticleBase), pointer   :: psolver_
      character (len=128)            :: infile_
    contains
      procedure(GStepper_ctor_interface), deferred :: GStepper_ctor! Constructor
      procedure         (step_interface), deferred :: step         ! step method
      procedure        (pstep_interface), deferred :: pstep        ! part+field step method
      generic,                              public :: gstep => step, pstep ! general method
  end type GStepperBase
  
  abstract interface
    subroutine GStepper_ctor_interface(this, traits, workspace, solver, psolver)
      use equationbase_mod
      use particlebase_mod
      import :: GStepperBase, GStepperTraits
      class (GStepperBase), intent(inout)                   :: this
      type(GStepperTraits), intent(inout)                   :: traits
      type    (GWorkspace), intent(inout), target           :: workspace
      class (Equationbase), intent   (in), target           :: solver
      class (ParticleBase), intent   (in), target, optional :: psolver
    end subroutine GStepper_ctor_interface

!   subroutine init_interface(this) 
!     import :: GStepperBase
!     class (GStepperBase), intent (inout) :: this
!   end subroutine init_interface

    subroutine step_interface(this, time, uin, uf, dt, uout) 
      use gstate_mod
      use equationbase_mod
      import :: GStepperBase
      class(GStepperBase), intent   (in) :: this
      real      (kind=GP), intent   (in) :: time, dt
      type   (GStateComp), intent(inout) :: uin(:),uf(:)
      type   (GStateComp), intent(inout) :: uout(:) 
    end subroutine step_interface

    ! Define step function interface for 
    ! particles plus fields:
    subroutine pstep_interface(this, time, uin, upin, uf, dt, uout, upout) 
      use gstate_mod
      use gpstate_mod
      import :: GStepperBase
      class (GStepperBase), intent  (in) :: this
      real      (kind=GP), intent   (in) :: time, dt
      type   (GStateComp), intent(inout) :: uin(:),uf(:)
      type  (GPStateComp), intent(inout) :: upin(:)
      type   (GStateComp), intent(inout) :: uout(:) 
      type  (GPStateComp), intent(inout) :: upout(:) 
    end subroutine pstep_interface
  end interface

contains

  ! ===================================================================
  ! Concrete methods inherited by all steppers
  ! ===================================================================
  
end module gstepperbase_mod


