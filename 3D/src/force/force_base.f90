! ===================================================================
! NAME       : ic_base.fpp
! DESCRIPTION: Forms base class for all initial conditions (IC)
! DATE       : 01/16/26 (PDM)
! ===================================================================

module forcebase_mod
  USE gstate_mod

  IMPLICIT NONE

  ! ================= Base class for all ICs  =======================
  ! Define an abstract base class
  type, abstract :: forceBase
    contains
      procedure(init_GForce_interface), deferred :: init_GForce
      procedure, public         :: update_GForce => update_none
  end type forceBase

! Right now we do just 1 force, in the future we may loop through a
! list to reuse, e.g., mechanical forcing from Navier-Stokes in MHD.
!  type           :: forceChain
!    class(forceBase), allocatable    :: list(:)
!  contains
!    procedure    :: init_fchain
!  end type forcechain
  
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

CONTAINS
  
  ! ===================================================================
  ! Concrete methods inherited by all forcing schemes
  ! ===================================================================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Concrete method to keep the forcing as is (no update)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine update_none(this,solver,state)
    use equationbase_mod
    use gstate_mod
    class   (forceBase), intent   (in)       :: this
    class(EquationBase), intent   (in)       :: solver
    type       (Gstate), intent(inout)       :: state(:)
    return
  end subroutine update_none

  
end module forcebase_mod


