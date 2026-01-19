! ===================================================================
! NAME       : ic_base.fpp
! DESCRIPTION: Forms base class for all initial conditions (IC)
! DATE       : 01/16/26 (PDM)
! ===================================================================

module icbase_mod
  USE gstate_mod

  IMPLICIT NONE

  ! ================= Base class for all ICs  =======================
  ! Define an abstract base class
  type, abstract :: icBase
    contains
      procedure(init_GState_interface), deferred :: init_GState
  end type icBase

! Right now we do just 1 IC, in the future we may loop through a
! list to reuse, e.g., velocity ICs from Navier-Stokes in MHD.
!  type           :: icChain
!    class(icBase), allocatable    :: list(:)
!  contains
!    procedure    :: init_chain
!  end type icchain
  
  abstract interface
    subroutine init_GState_interface(this, solver, state)
      USE equationbase_mod
      USE gstate_mod
      import :: icBase
      class      (icBase),   intent   (in) :: this
      class(EquationBase),   intent   (in) :: solver
      type       (Gstate),   intent(inout) :: state(:)
    end subroutine
  end interface

end module icbase_mod


