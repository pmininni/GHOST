! =====================================================================
! NAME       : canuto_stepper.f90
!
!              Performs time stepping using method of 
!              Canuto et al. Spectral Methods in Fluid
!              Dynamics. This is a low storage time stepper
!              for a specifiable order, though it is strictly
!              speaking of full truncation order only for
!              norder = 2. Explicit number of stages are
!              not required. While norder > 2 may not yield
!              a truncation of that order, it can still provide
!              benefit.
!
! INPUT FILE : Stepper looks for a "&stepper" namelist with:
!                norder  : Stepper order (arbitrary)
!
! DATE       : 2/8/26 (DLR)
! =====================================================================

module canuto_stepper_mod
  use iso_c_binding
  use class_GWorkspace3D
  use gstate_mod
  use stepperbase_mod


  implicit none

  ! Define callback function interface:
! abstract interface
!    subroutine dudt_interface(this, time, uin, uf, dt, dudt)   
!      use gstate_mod
!      import :: CanutoStepper
!      class(CanutoStepper), intent   (in)         :: this
!      real      (kind=GP), intent   (in)         :: time, dt
!      type   (GStateComp), intent(inout), target :: uin(:),uf(:)
!      type   (GStateComp), intent(inout)         :: dudt(:) 
!    end subroutine dudt_interface
! end interface

  ! ================= Global parameters ===============================

  ! ================= Stepper traits ===================================

  ! ================= Stepper ==========================================
  ! Define class:
  type, extends(GStepperBase) :: CanutoStepper
    ! Member data:
    type(GWorkspace), pointer     :: workspace_
    logical                       :: binit_=.false. ! is initialized?
    logical                       :: busing_butcher_=.true. ! using Butcher tableau?
    integer                       :: myrank_   ! MPI rank
    integer                       :: nprocs_   ! MPI rank 
    type(GStepperTraits)          :: traits_   ! GStepper traits


    procedure(dudt_interface), pointer, nopass :: callback_ => null()

  CONTAINS
    procedure, public :: init                        ! initialize
    procedure, public :: set_callback => set_callback_impl ! set RHS callback method
    procedure, public :: step         =>  step_impl  ! take one timestep
    procedure, public :: GStepper_ctor =>  CanutoStepper_ctor 
    final             :: CanutoStepper_dtor

  end type CanutoStepper

CONTAINS

  ! ===================================================================
  ! Stepper-specific methods
  ! ===================================================================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Constructor
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine CanutoStepper_ctor(this, traits, workspace)
    type(GStepperTraits), intent(inout), target :: traits
    type  (GWorkspace)  , intent(inout), target :: workspace

    this%workspace_ => workspace
    if (.not. associated(this%workspace_)) then
      stop 'CanutoStepper::CanutoStepper_ctor: Worskpace not associated'
    endif

    call this%init(traits)
  end subroutine CanutoStepper_ctor


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Destructor
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine CanutoStepper_dtor(this) 
    type  (CanutoStepper), intent(inout) :: this

    if (associated(this%workspace_))   nullify(this%workspace_)
    
  end subroutine CanutoStepper_dtor


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Subroutine to initialize the stepper
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init(this, traits)
!   use commtypes
    class(CanutoStepper), intent (inout) :: this
    type(GStepperTraits), intent    (in) :: traits

    this%traits_ = traits;

  end subroutine init

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Subroutine to set RHS callback function
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine set_callback_impl(this, fcn_callback)
    class (CanutoStepper), intent  (inout) :: this
    procedure(callback_interface), pointer :: fcn_callback

    this%callback_ => fcn_callback

  end subroutine set_callback_impl


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
       
    if ( size(uin) .ne. this%traits_%nstate &
     .or.size(uout) .ne. this%traits_%nstate  ) then
      stop 'CanutoStepper::step: Inconsistent input state'
    endif
    
    if ( .not. associated(this%callback_) ) then
      stop 'CanutoStepper::step: RHS callback function not set'
    endif

    do o = this%traits_%norder,1,-1
      eff_dt = dt/real(o,kind=GP)
      call this%callback(time, uout, uf, eff_dt, uout)
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


end module canuto_stepper_mod
