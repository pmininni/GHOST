! ===================================================================
! NAME       : icp_base.fpp
! DESCRIPTION: Base class for initial conditions of particles (ICP)
! DATE       : 03/21/26 (PDM)
! ===================================================================

module icpbase_mod

  IMPLICIT NONE

  ! ================= Base class for all ICPs  ======================
  ! Define an abstract base class
  type, abstract :: icpBase
    contains
      procedure(init_GPState_interface), deferred :: init_GPState
  end type icpBase

  ! Concrete data type to do chain operation of ICPs
  type           :: icpChain
    class(icpBase), allocatable :: icp
  end type icpChain
  
  abstract interface
    subroutine init_GPState_interface(this,pde,fluidstate,psolver,pstate)
      USE equationbase_mod
      USE particlebase_mod
      USE gpstate_mod
      USE gstate_mod
      import :: icpBase
      class     (icpBase),         intent   (in) :: this
      class(EquationBase),         intent   (in) :: pde
      type   (GStateComp), target, intent   (in) :: fluidstate(:) 
      class(ParticleBase),         intent(inout) :: psolver
      type  (GPStateComp),         intent(inout) :: pstate(:)
    end subroutine
  end interface

contains

  ! ================= Chain operator for ICPs =======================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Initializes all states from a list of ICPs
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_allpstates(chain,pde,fluidstate,psolver,pstate)
    USE equationbase_mod
    USE particlebase_mod
    USE gpstate_mod
    USE gstate_mod
    implicit none
    type     (icpChain), intent   (in) :: chain(:)
    class(EquationBase), intent   (in) :: pde
    type   (GStateComp), intent   (in) :: fluidstate(:) 
    class(ParticleBase), intent(inout) :: psolver
    type  (GPstateComp), intent(inout) :: pstate(:)
    integer                            :: i
    do i = 1,size(chain)
      call chain(i)%icp%init_GPState(pde,fluidstate,psolver,pstate)
    end do
  end subroutine init_allpstates
  
end module icpbase_mod
