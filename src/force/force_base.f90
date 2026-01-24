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
      complex(kind=GP), pointer, dimension(:,:,:)  :: fxold_,fyold_,fzold_
      complex(kind=GP), pointer, dimension(:,:,:)  :: fxnew_,fynew_,fznew_
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
      type       (Gstate),   intent(inout) :: state(:)
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
      type       (Gstate),   intent(inout) :: state(:)
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
    type       (Gstate), intent(inout) :: state(:)
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
    class(EquationBase), intent   (in) :: solver
    type       (Gstate), intent(inout) :: state(:)
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
