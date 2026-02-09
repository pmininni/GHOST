! ===================================================================
! NAME       : gstepper_base.fpp
! DESCRIPTION: Forms base class for all time steppers
! DATE       : 2/8/26 (DLR)
! ===================================================================

module gstepperbase_mod
  use class_GWorkspace3D
  use gstate_mod

  implicit none

  ! ================= Stepper traits ===================================
  type, public  :: GStepperTraits
    integer            :: itype  = 1
    integer            :: norder = 2
    integer            :: nstage = 2
    integer            :: nstate = 0
    character(len=128) :: sname  = 'CANUTO' ! Currently: 'CANUTO' or 'GEXRK'
  end type


  ! ================= Base class for all steppers =======================
  ! Define an abstract base class
  type, abstract :: GStepperBase
      type(GWorkspace), pointer     :: workspace_
      integer                       :: myrank_   ! MPI rank
      integer                       :: nprocs_   ! MPI rank 
      character(len=128)            :: infile_
    contains
      procedure(GStepper_ctor_interface), deferred :: GStepper_ctor ! Constructor
!     procedure(init_interface),         deferred :: init        ! init method
      procedure(step_interface),         deferred :: step        ! step method
      procedure(dudt_interface),         deferred :: set_callback! step method
  end type GStepperBase

  abstract interface

     ! Define step constructor:
     subroutine GStepper_ctor_interface(this, traits, workspace)
       use class_GWorkspace3D
       import :: GStepperBase
       class (GStepperBase), intent(inout)        :: this
       type(GStepperTraits), intent(inout)        :: traits
       type   (GWorkspace), intent(inout), target:: workspace
     end subroutine GStepper_ctor_interface

!    subroutine init_interface(this) 
!      import :: GStepperBase
!      class (GStepperBase), intent (inout) :: this
!    end subroutine init_interface

     ! Define step function interface:
     subroutine step_interface(this, time, uin, uf, dt, uout) 
       use gstate_mod
       import :: GStepperBase
       class (GStepperBase), intent   (in)         :: this
       real      (kind=GP), intent   (in)         :: time, dt
       type   (GStateComp), intent(inout), target :: uin(:),uf(:)
       type   (GStateComp), intent(inout)         :: uout(:) 
     end subroutine step_interface

     ! Define callback function interface:
     subroutine dudt_interface(this, time, uin, uf, dt, dudt)
       use gstate_mod
       import :: GStepperBase
       class (GStepperBase), intent   (in)         :: this
       real      (kind=GP), intent   (in)         :: time, dt
       type   (GStateComp), intent(inout), target :: uin(:),uf(:)
       type   (GStateComp), intent(inout)         :: dudt(:)
     end subroutine dudt_interface

  end interface

CONTAINS

  ! ===================================================================
  ! Concrete methods inherited by all steppers
  ! ===================================================================

  
end module gstepperbase_mod


