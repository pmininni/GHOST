! =====================================================================
! NAME       : canuto_stepper.f90
!
!              Performs time stepping using the modified RK method 
!              of Canuto et al., Spectral Methods in Fluid
!              Dynamics. This is a low storage time stepper
!              for a specifiable order, though it is strictly
!              speaking of full truncation order only for
!              norder = 2, or for linear ODEs/PDEs. Explicit number
!              of stages are not required. While norder > 2 may
!              not yield a truncation of that order in nonlinear
!              PDEs, it can still provide a benefit. Particles in
!              this stepper are sync'd at every substepper stage.
!
! INPUT FILE : Stepper looks for a "&stepper" namelist with:
!                sname   : 'TRADITIONAL' to use this time stepper
!                itype   : Ignored by this solver. The stepper type:
!                    =1 => Butcher type (norder = nstage, this solver)
!                    =2 => Mixed type (norder != nstage; not available)
!                    =3 => SSP type (strong stability preserving; not available)
!                norder  : Stepper order (>=1)
!                nstage  : No. stepper stages (ignored by this solver)
!
! DATE       : 2/8/26 (DLR)
! =====================================================================

module canuto_stepper_mod
  use gstepperbase_mod
  use gdevice, only: gdev_active
  implicit none

  ! ================= Global parameters ===============================

  ! ================= Stepper traits ===================================

  ! ================= Stepper ==========================================
  ! Define class:
  type, extends(GStepperBase) :: CanutoStepper
    ! Member data (extends GStepperBase member data)
    type(GStepperTraits)      :: traits_                     ! GStepper traits
  contains
    procedure, public :: init          => init_impl          ! initialize
    procedure, public :: step          => step_impl          ! take one field timestep
    procedure, public :: pstep         => pstep_impl         ! take one part  timestep
    procedure, public :: cstep         => cstep_impl         ! take one part+field timestep
    procedure, public :: GStepper_ctor => CanutoStepper_ctor
    final             :: CanutoStepper_dtor
  end type CanutoStepper

contains

  ! ===================================================================
  ! Stepper-specific methods
  ! ===================================================================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Constructor
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine CanutoStepper_ctor(this, traits, workspace, solver, psolver)
    class(CanutoStepper), intent(inout)                   :: this
    type(GStepperTraits), intent(inout)                   :: traits
    type    (GWorkspace), intent(inout), target           :: workspace
    class (EquationBase), intent   (in), target           :: solver
    class (ParticleBase), intent   (in), target, optional :: psolver
    this%workspace_ => workspace
    this%solver_    => solver
    if (present(psolver)) then
      this%psolver_ => psolver
    else
      this%psolver_ => null()
    endif
    if (.not. associated(this%workspace_)) then
      stop 'CanutoStepper::CanutoStepper_ctor: Worskpace not associated'
    endif
    call this%init(traits)
  end subroutine CanutoStepper_ctor


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Destructor
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine CanutoStepper_dtor(this) 
    type(CanutoStepper), intent(inout) :: this
    if (associated(this%workspace_)) nullify(this%workspace_)
    if (associated(this%solver_))    nullify(this%solver_)
    if (associated(this%psolver_))   nullify(this%psolver_)
  end subroutine CanutoStepper_dtor


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Subroutine to initialize the stepper
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_impl(this, traits)
    class(CanutoStepper), intent(inout) :: this
    type(GStepperTraits), intent   (in) :: traits
    this%traits_ = traits;
  end subroutine init_impl


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Implementation function to take one step of fields
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine step_impl(this, time, uin, uf, dt, uout)
    implicit none

    class(CanutoStepper), intent(inout) :: this
    type    (GStateComp), intent(inout), target :: uin(:), uf(:), uout(:)
    real       (kind=GP), intent   (in) :: time, dt
    real       (kind=GP)                :: eff_dt
    integer                             :: i,j,k,o,state_size,ic
    logical                             :: bret
    complex    (kind=GP), pointer       :: po(:,:,:), pi(:,:,:)
       
    if ( size (uin) .ne. this%traits_%nstate &
     .or.size(uout) .ne. this%traits_%nstate  ) then
      stop 'CanutoStepper::step: Inconsistent input state'
    endif

    do o = this%traits_%norder,1,-1
      eff_dt = dt/real(o,kind=GP)
      ! We assume that at input, uout has a copy of uin
      call this%solver_%dudt(time, uout, uf, eff_dt, uout)
      ! The components are addressed through pointers: indexing
      ! uout(ic)%ccomp inside a target region does not work with flang
      do ic = 1,this%traits_%nstate
        po => uout(ic)%ccomp
        pi => uin (ic)%ccomp
#if defined(GHOST_GPU)
!$omp target teams distribute parallel do collapse(3) if(target: gdev_active)
        do i = ista,iend
          do j = 1,ny
            do k = 1,nz
#else
!$omp parallel do collapse(2) private (k)
        do i = ista,iend
          do j = 1,ny
            do concurrent (k=1:nz)
#endif
              po(k,j,i) = pi(k,j,i) + eff_dt*po(k,j,i)
            end do
          end do
        end do
      end do
    end do
  end subroutine step_impl


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Implementation function to take one step of the particles.
  !! This time stepper does not perform computation of feedback
  !! forces on the fluid (one-way coupled), and evolves the
  !! particles in a fixed velocity field during the substepping 
  !! stages. This method is mainly inteded for manual
  !! integration of multiple sets of (same type) one-way
  !! coupled particles.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine pstep_impl(this, time, uin, upin, dt, upout)
    use gpstate_mod
    implicit none

    class(CanutoStepper), intent(inout)         :: this
    type    (GStateComp), intent(inout)         :: uin (:)
    type   (GPStateComp), intent(inout), target, allocatable :: upin(:), upout(:)
    real       (kind=GP), intent   (in)         :: time, dt
    real       (kind=GP)                        :: eff_dt
    integer                                     :: i,j,k,o,state_size
    integer                                     :: ic,ip,nparts
    logical                                     :: bret
       
    if ( size (uin) .ne. this%traits_%nstate ) then
      stop 'CanutoStepper::pstep: Inconsistent input state'
    endif
    if ( size(upin) .ne. this%traits_%npstate &
    .or.size(upout) .ne. this%traits_%npstate  ) then
      stop 'CanutoStepper::pstep: Inconsistent particle state'
    endif
    if ( .not. this%traits_%dopart ) then
      stop 'CanutoStepper::pstep: The particle stepper was not initialized'
    endif
    if ( this%psolver_%hasfeedback_ ) then
      stop 'CanutoStepper::pstep: Particles are not one-way coupled'
    endif
    nparts = size(upin(1)%rcomp)

    do o = this%traits_%norder,1,-1
      eff_dt = dt/real(o,kind=GP)
      ! We compute dpdt; we assume that at input, upout has a copy of upin
      call this%psolver_%dpdt(time, this%solver_, uin, upout, eff_dt, upout)
      ! Update particles:
!$omp parallel do collapse(2)
      do ip = 1,this%traits_%npstate
        do k = 1,this%psolver_%nparts_
          upout(ip)%rcomp(k) = upin(ip)%rcomp(k) + eff_dt*upout(ip)%rcomp(k)
        end do
      end do
      ! Sync particles at each stage
      call this%psolver_%end_stage(upin,upout)
    end do ! end, o-loop
  end subroutine pstep_impl


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Implementation function to take one *coupled* step for
  !! particles + fields. Feedback forces are computed if 
  !! needed. Both particles and fields are evolved
  !! simultaneously in the substepping stages.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine cstep_impl(this, time, uin, upin, uf, dt, uout, upout)
    use gpstate_mod
    implicit none

    class(CanutoStepper), intent(inout)         :: this
    type    (GStateComp), intent(inout)         :: uin (:), uf(:), uout(:)
    type    (GStateComp), allocatable           :: fdbk(:)
    type   (GPStateComp), intent(inout), target, allocatable :: upin(:), upout(:)
    real       (kind=GP), intent   (in)         :: time, dt
    real       (kind=GP)                        :: eff_dt
    integer                                     :: i,j,k,o,state_size
    integer                                     :: ic,ip,nparts
    logical                                     :: bret
       
    if ( size (uin) .ne. this%traits_%nstate &
     .or.size(uout) .ne. this%traits_%nstate  ) then
      stop 'CanutoStepper::cstep: Inconsistent input state'
    endif
    if ( size(upin) .ne. this%traits_%npstate &
    .or.size(upout) .ne. this%traits_%npstate  ) then
      stop 'CanutoStepper::cstep: Inconsistent particle state'
    endif
    if ( .not. this%traits_%dopart ) then
      stop 'CanutoStepper::cstep: The particle stepper was not initialized'
    endif
    if ( this%psolver_%hasfeedback_ ) then  ! We alloc arrays if we have feedback
      call GState_alloc(fdbk, this%traits_%nstate)
    endif
    nparts = size(upin(1)%rcomp)

    do o = this%traits_%norder,1,-1
      eff_dt = dt/real(o,kind=GP)
      ! If needed we compute the feedback of the particles in the fluid 
      if ( this%psolver_%hasfeedback_ ) call this%psolver_%feedback(upout,fdbk)
      ! We compute dpdt first; we assume that at input, upout has a copy of upin
      call this%psolver_%dpdt(time, this%solver_, uout, upout, eff_dt, upout)
      ! We assume that at input, uout has a copy of uin
      call this%solver_%dudt (time, uout, uf, eff_dt, uout)
      ! If needed we update the forcing to include the feedback
      if ( this%psolver_%hasfeedback_ ) then
        call this%psolver_%feedback(upout,fdbk)
        do ic = 1,this%traits_%nstate
          uout(ic)%ccomp = uout(ic)%ccomp + fdbk(ic)%ccomp
        end do
      endif
      ! Update fields and particles inside one parallel region, so the
      ! team is forked once per stage for both updates:
!$omp parallel private (k)
!$omp do collapse(3)
      do ic = 1,this%traits_%nstate
        do i = ista,iend
          do j = 1,ny
            do concurrent (k=1:nz)
              uout(ic)%ccomp(k,j,i) = uin(ic)%ccomp(k,j,i) + eff_dt*uout(ic)%ccomp(k,j,i)
            end do
          end do
        end do
      end do
!$omp end do
!$omp do collapse(2)
      do ip = 1,this%traits_%npstate
        do k = 1,this%psolver_%nparts_
          upout(ip)%rcomp(k) = upin(ip)%rcomp(k) + eff_dt*upout(ip)%rcomp(k)
        end do
      end do
!$omp end do
!$omp end parallel
      ! Sync particles at each stage
      call this%psolver_%end_stage(upin,upout)
    end do ! end, o-loop
    if ( this%psolver_%hasfeedback_ ) call GState_dealloc(fdbk)
  end subroutine cstep_impl

end module canuto_stepper_mod
