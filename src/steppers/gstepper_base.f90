! ===================================================================
! NAME       : gstepper_base.fpp
! DESCRIPTION: Forms base class for all time steppers
! DATE       : 2/8/26 (DLR)
! ===================================================================

module gstepperbase_mod
!  use class_GWorkspace3D
  use equationbase_mod
!  use gstate_mod
  implicit none

  ! ================= Stepper traits ===================================
  type, public  :: GStepperTraits
    integer            :: itype   = 1 ! TRADITIONAL
    integer            :: norder  = 2 ! order 2
    integer            :: nstage  = 2 ! arbitrary
    integer            :: nstate  = 0 ! no. continuum state components
    integer            :: npstate = 0 ! no. particle state components
    character(len=128) :: sname  = 'TRADITIONAL' ! Currently: 'TRADITIONAL' or 'GEXRK'
  end type

  ! ================= Base class for all steppers =======================
  ! Define an abstract base class
  type, abstract :: GStepperBase
      type(GWorkspace), pointer     :: workspace_
      integer                       :: myrank_   ! MPI rank
      integer                       :: nprocs_   ! MPI rank 
      logical                       :: doparts_ = .false. ! use particle interfaces?   
      character     (len=128)       :: infile_
      procedure(dudt_interface), pointer, nopass :: callback_ => null()
    contains
      procedure(GStepper_ctor_interface), deferred :: GStepper_ctor! Constructor
      procedure(step_interface),          deferred :: step         ! step method
!     procedure(pstep_interface),         deferred :: pstep        ! part+field step method
      procedure(set_callback_interface),  deferred :: set_callback ! callback to RHS method
!     procedure(pdudt_interface),         deferred :: set_pcallback! part+field RHS method
  end type GStepperBase

  abstract interface

     ! Define step constructor:
     subroutine GStepper_ctor_interface(this, traits, workspace, nparts)
       use class_GWorkspace3D
       import :: GStepperBase, GStepperTraits
       class (GStepperBase), intent(inout)           :: this
       type(GStepperTraits), intent(inout)           :: traits
       type    (GWorkspace), intent(inout), target   :: workspace
       integer             , intent   (in), optional :: nparts
!!$
!!$       if ( nparts .gt. 0 ) then
!!$         this%doparts_ = .true.
!!$         if ( workspace%get_nparts() .ne. nparts ) then
!!$           stop 'GStepperBase_ctor:: Incompatible workspace' 
!!$         endif
!!$       endif
     end subroutine GStepper_ctor_interface

!    subroutine init_interface(this) 
!      import :: GStepperBase
!      class (GStepperBase), intent (inout) :: this
!    end subroutine init_interface

     ! Define step function interface:
     subroutine step_interface(this, solver, time, uin, uf, dt, uout) 
       use gstate_mod
       use equationbase_mod
       import :: GStepperBase
       class(GStepperBase), intent   (in) :: this
       class(EquationBase), intent   (in) :: solver
       real      (kind=GP), intent   (in) :: time, dt
       type   (GStateComp), intent(inout) :: uin(:),uf(:)
       type   (GStateComp), intent(inout) :: uout(:) 
     end subroutine step_interface

!!$     ! Define step function interface for 
!!$     ! particles plus fields:
!!$     subroutine pstep_interface(this, time, uin, upin, uf, dt, uout, upout) 
!!$       use gstate_mod
!!$       use gpstate_mod
!!$       import :: GStepperBase
!!$       class (GStepperBase), intent  (in) :: this
!!$       real      (kind=GP), intent   (in) :: time, dt
!!$       type   (GStateComp), intent(inout) :: uin(:),uf(:)
!!$       type  (GPStateComp), intent(inout) :: upin(:)
!!$       type   (GStateComp), intent(inout) :: uout(:) 
!!$       type  (GPStateComp), intent(inout) :: upout(:) 
!!$     end subroutine pstep_interface

     ! Define callback function interface:
     subroutine set_callback_interface(this, fcn_callback)
       use gstate_mod
       import :: GStepperBase
       class (GStepperBase), intent(inout) :: this
       procedure(dudt_interface)   :: fcn_callback
     end subroutine set_callback_interface

!!$     ! Define callback function particle + field interface:
!!$     subroutine pdudt_interface(this, time, uin, upin, uf, dt, dudt, pdudt)
!!$       use gstate_mod
!!$       use gpstate_mod
!!$       import :: GStepperBase
!!$       class (GStepperBase), intent(inout)         :: this
!!$       real      (kind=GP),  intent   (in)         :: time, dt
!!$       type   (GStateComp),  intent(inout), target :: uin(:),uf(:)
!!$       type  (GPStateComp),  intent(inout), target :: upin(:)
!!$       type   (GStateComp),  intent(inout)         :: dudt(:)
!!$       type  (GPStateComp),  intent(inout)         :: pdudt(:)
!!$     end subroutine pdudt_interface

  end interface

contains

  ! ===================================================================
  ! Concrete methods inherited by all steppers
  ! ===================================================================
  
end module gstepperbase_mod


