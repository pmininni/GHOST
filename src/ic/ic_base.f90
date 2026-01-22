! ===================================================================
! NAME       : ic_base.fpp
! DESCRIPTION: Forms base class for all initial conditions (IC)
! DATE       : 01/16/26 (PDM)
! ===================================================================

module icbase_mod

  IMPLICIT NONE

  ! ================= Base class for all ICs  =======================
  ! Define an abstract base class
  type, abstract :: icBase
    contains
      procedure(init_GState_interface), deferred :: init_GState
  end type icBase

  ! Concrete data type to do chain operation of ICs
  type           :: icChain
    class(icBase), allocatable :: ic
  end type icchain
  
  abstract interface
    subroutine init_GState_interface(this, solver, state)
      USE equationbase_mod
      USE gstate_mod
      import :: icBase
      class      (icBase), intent   (in) :: this
      class(EquationBase), intent   (in) :: solver
      type       (Gstate), intent(inout) :: state(:)
    end subroutine
  end interface

contains

  ! ================= Chain operator for ICs ========================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Initializes all states from a list of ICs
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_allstates(chain, solver, state)
    USE equationbase_mod
    USE gstate_mod
    implicit none
    type      (icChain), intent   (in) :: chain(:)
    class(EquationBase), intent   (in) :: solver
    type       (Gstate), intent(inout) :: state(:)
    integer                            :: i
    do i = 1,size(chain)
      call chain(i)%ic%init_GState(solver,state)
    end do
  end subroutine init_allstates
  
end module icbase_mod


