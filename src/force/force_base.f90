! ===================================================================
! NAME       : force_base.fpp
! DESCRIPTION: Forms base class for all forcing methods
! DATE       : 01/16/26 (PDM)
! ===================================================================

module forcebase_mod
  USE fprecision
  
  IMPLICIT NONE

  ! ================= Base class for all forcings  ==================
  ! Define an abstract base class
  type, abstract :: forceBase
    contains
      procedure(init_GForce_interface), deferred   :: init_GForce
  end type forceBase

  type, abstract :: forceUpdt
      complex(kind=GP), pointer, dimension(:,:,:)  :: fxold_ => null()
      complex(kind=GP), pointer, dimension(:,:,:)  :: fyold_ => null()
      complex(kind=GP), pointer, dimension(:,:,:)  :: fzold_ => null()
      complex(kind=GP), pointer, dimension(:,:,:)  :: fxnew_ => null()
      complex(kind=GP), pointer, dimension(:,:,:)  :: fynew_ => null()
      complex(kind=GP), pointer, dimension(:,:,:)  :: fznew_ => null()
      logical                                      :: binit_
      ! Set by the update method when it modifies the host copy of the
      ! forcing state in the current step (the main program then copies
      ! the state to the device). Reset by update_forcing before each
      ! call. Methods that change the forcing every step (shuffle) do
      ! the intermediate updates on the device copies and leave the host
      ! copies untouched on those steps.
      logical                                      :: changed_ = .false.
    contains
      procedure(update_GForce_interface), deferred :: update_GForce
      ! Second pass of the update (after all the methods of the chain
      ! ran update_GForce): blend of the old and new states of the
      ! slowly evolving (shuffle) methods. No-op by default.
      procedure                                    :: blend_GForce => blend_GForce_noop
  end type forceUpdt

  ! Concrete data type to do chain operation of forcing methods
  type           :: forceChain
      class(forceBase), allocatable        :: force
      class(forceUpdt), allocatable        :: update
  end type forceChain
 
  abstract interface
    subroutine init_GForce_interface(this, solver, state)
      USE equationbase_mod
      USE gstate_mod
      import :: forceBase
      class   (forceBase),   intent   (in) :: this
      class(EquationBase),   intent(inout) :: solver
      type   (GStateComp),   intent(inout) :: state(:)
    end subroutine
  end interface

  abstract interface
    subroutine update_GForce_interface(this, force, solver, state)
      USE equationbase_mod
      USE gstate_mod
      import :: forceBase
      import :: forceUpdt
      class   (forceUpdt),   intent(inout) :: this
      class   (forceBase),   intent   (in) :: force
      class(EquationBase),   intent(inout) :: solver
      type   (GStateComp),   intent(inout) :: state(:)
    end subroutine
  end interface

CONTAINS
  
  ! ================= Chain operators for ICs =======================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Initializes all forcing states from a list
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_forcing(chain, solver, state)
    USE equationbase_mod
    USE gstate_mod
    implicit none
    type   (forceChain), intent   (in) :: chain(:)
    class(EquationBase), intent(inout) :: solver
    type   (GStateComp), intent(inout) :: state(:)
    integer                            :: i
    do i = 1,size(chain)
      call chain(i)%force%init_GForce(solver,state)
    end do
  end subroutine init_forcing

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Updates all forcing states from a list
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine update_forcing(chain, solver, state)
    USE equationbase_mod
    USE gstate_mod
    USE status
    implicit none
    type   (forceChain), intent(inout) :: chain(:)
    class(EquationBase), intent(inout) :: solver
    type   (GStateComp), intent(inout) :: state(:)
    integer                            :: i
    ! First pass: new random states, phase shifts, etc. (host code, at
    ! the correlation time). A method generating a state correlated with
    ! another forcing (random_fb with corr > 0) sees the newly generated
    ! state of the latter, as in the original code
    do i = 1,size(chain)
      if ( allocated(chain(i)%update) ) then
        chain(i)%update%changed_ = .false.
        call chain(i)%update%update_GForce(chain(i)%force,solver,state)
      endif
    end do
    ! Second pass: blends of the slowly evolving methods (every step)
    do i = 1,size(chain)
      if ( allocated(chain(i)%update) ) then
        call chain(i)%update%blend_GForce(solver,state)
      endif
    end do
    if (timef.eq.fstep) timef = 0 ! Updates state counters
    timef = timef + 1
  end subroutine update_forcing

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Default second pass of the update methods: nothing to do
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine blend_GForce_noop(this, solver, state)
    USE equationbase_mod
    USE gstate_mod
    implicit none
    class   (forceUpdt),   intent(inout) :: this
    class(EquationBase),   intent(inout) :: solver
    type   (GStateComp),   intent(inout) :: state(:)
  end subroutine blend_GForce_noop

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Logical function, true if any update method modified the
  !! host copy of the forcing state in the last call to
  !! update_forcing: the state must then be copied to the
  !! device (GState_update_to) before the time step
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  logical function forcing_changed(chain) result(bchanged)
    implicit none
    type   (forceChain), intent(in) :: chain(:)
    integer                         :: i
    bchanged = .false.
    do i = 1,size(chain)
      if ( allocated(chain(i)%update) ) then
        if ( chain(i)%update%changed_ ) bchanged = .true.
      endif
    end do
  end function forcing_changed

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Blend of the old and new forcing states of the shuffle
  !! update methods, f = (1-rmp).fold + rmp.fnew. With ondevice
  !! the kernel runs on the device copies of the three arrays
  !! (offload builds; the host copy of f is left stale and is
  !! refreshed by the main program before the diagnostics),
  !! otherwise on the host copies. In host builds it always
  !! runs on the host.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine blend_forcing(f, fold, fnew, rmp, ondevice)
    USE pseudospec_fluid, only: saxpby_c
    USE gdevice, only: gdev_active
    USE grid
    USE mpivars
    implicit none
    complex(kind=GP), intent(inout), dimension(nz,ny,ista:iend) :: f
    complex(kind=GP), intent   (in), dimension(nz,ny,ista:iend) :: fold,fnew
    real   (kind=GP), intent   (in)                             :: rmp
    logical         , intent   (in)                             :: ondevice
    logical                                                     :: bsave
    bsave = gdev_active
    if ( ondevice ) gdev_active = .true.
    call saxpby_c(f,fold,1.0_GP-rmp,fnew,rmp)
    gdev_active = bsave
  end subroutine blend_forcing

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Logical function to keep track of static forcings:
  !! True if no forcing method has an update scheme, i.e. the
  !! forcing is constant in time (the arrays never change after
  !! init_forcing: the host copies are always current and the
  !! state never needs to be copied to or from a device again)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  logical function forcing_is_static(chain) result(bstatic)
    implicit none
    type   (forceChain), intent(in) :: chain(:)
    integer                         :: i
    bstatic = .true.
    do i = 1,size(chain)
      if ( allocated(chain(i)%update) ) bstatic = .false.
    end do
  end function forcing_is_static
  
end module forcebase_mod
