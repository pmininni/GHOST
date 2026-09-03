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
    contains
      procedure(update_GForce_interface), deferred :: update_GForce
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
      class(EquationBase),   intent   (in) :: solver
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
      class(EquationBase),   intent   (in) :: solver
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
    class(EquationBase), intent   (in) :: solver
    type   (GStateComp), intent(inout) :: state(:)
    integer                            :: i
    do i = 1,size(chain)
      call chain(i)%force%init_GForce(solver,state)
    end do
  end subroutine init_forcing


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Updates all forcing states from a list
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! True if no forcing method has an update scheme, i.e. the
  !! forcing is constant in time (the arrays never change after
  !! init_forcing and need not be copied to the device again)
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

  subroutine update_forcing(chain, solver, state)
    USE equationbase_mod
    USE gstate_mod
    USE status
    implicit none
    type   (forceChain), intent(inout) :: chain(:)
    class(EquationBase), intent   (in) :: solver
    type   (GStateComp), intent(inout) :: state(:)
    integer                            :: i
    do i = 1,size(chain)
      if ( allocated(chain(i)%update) ) then
        call chain(i)%update%update_GForce(chain(i)%force,solver,state)
      endif
    end do
    if (timef.eq.fstep) timef = 0 ! Updates state counters
    timef = timef + 1
  end subroutine update_forcing

end module forcebase_mod
