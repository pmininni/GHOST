! =====================================================================
! NAME       : gexrk_stepper.f90
!  
!              Performs time stepping using explicit RK with user-specified
!              order and number of stages. Currently stepper supports
!              only explicit RK schemes where norder=nstage, but the 
!              class is designed to eventually accommodate both 
!              GEXRK_MIXED (itype=2) and GEXRK_SSP (itype=3) types.
!              Particles in this stepper are sync'd only at the end of
!              a full time step.
!              
! INPUT FILE : Stepper looks for a "&stepper" namelist with:
!                sname   : 'GEXRK' to use this time stepper
!                itype   : Stepper type:
!                    =1 => Butcher type (norder = nstage)
!                    =2 => Mixed type (norder != nstage; not available)
!                    =3 => SSP type (strong stability preserving; not available)
!                norder  : Stepper order (1-4 currently)
!                nstage  : No. stepper stages (1-4 currently)
!
! DATE       : 2/2/26 (DLR)
! =====================================================================

module gexrk_stepper_mod
  use gstepperbase_mod
  use iso_c_binding
  implicit none

  ! ================= Global parameters ===============================
  enum, bind(c)
    enumerator :: GEXRK_BUTCHER = 1, GEXRK_MIXED = 2, GEXRK_SSP = 3
  end enum

  ! ================= Stepper traits ===================================
  
  ! ================= Stepper ==========================================
  ! Define class:
  type, extends(GStepperBase) :: GExRKStepper
    ! Member data (extends GStepperBase member data)
    type(GStepperTraits)   :: traits_                 ! GStepper traits
    logical                :: busing_butcher_= .true. ! using Butcher tableau?
    real    (kind=GP), allocatable, dimension  (:) :: alpha_, c_
    real    (kind=GP), allocatable, dimension(:,:) :: beta_
    type (GCStateArr), allocatable, dimension  (:) :: K_
    type (GStateComp), allocatable, dimension  (:) :: utmp_
    type (GPStateArr), allocatable, dimension  (:) :: pK_
    type(GPStateComp), allocatable, dimension  (:) :: putmp_
  contains
    procedure, public  :: init          => init_impl   ! initialize
    procedure, public  :: step          => step_impl   ! take one field timestep
    procedure, public  :: pstep         => pstep_impl  ! take one part  timestep
    procedure, public  :: cstep         => cstep_impl  ! take one part+field timestep
    procedure, public  :: GStepper_ctor => GEXRKStepper_ctor 
    final              :: GExRKStepper_dtor
    procedure, private :: init_butcher   ! initialize Butcher table
    procedure, private :: init_mixed     ! initialize mixed type
    procedure, private :: init_ssp       ! initialize for SSP
    procedure, private :: step_butcher   ! take a field Butcher step
    procedure, private :: step_mixed     ! take a field mixed step
    procedure, private :: step_ssp       ! take a field SSP step
    procedure, private :: pstep_butcher  ! take a part  Butcher step
    procedure, private :: pstep_mixed    ! take a part  mixed step
    procedure, private :: pstep_ssp      ! take a part  SSP step
    procedure, private :: cstep_butcher  ! take a part+field Butcher step
    procedure, private :: cstep_mixed    ! take a part+field mixed step
    procedure, private :: cstep_ssp      ! take a part+field SSP step
    procedure, private :: alloc_tmp      ! allocate tmp arrays
    procedure, private :: dealloc_tmp    ! deallocate tmp arrays
  end type GExRKStepper

contains

  ! ===================================================================
  ! Stepper-specific methods
  ! ===================================================================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Constructor
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine GExRKStepper_ctor(this, traits, workspace, solver, psolver)
    class (GExRKStepper), intent(inout)                   :: this
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
      stop 'GExRKStepper::GExRKStepper_ctor: Worskpace not associated'
    endif
    call this%init(traits)
    call this%alloc_tmp()
  end subroutine GExRKStepper_ctor


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Destructor
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine GExRKStepper_dtor(this) 
    type  (GExRKStepper), intent(inout) :: this
    if (associated(this%workspace_)) nullify(this%workspace_)
    if (associated(this%solver_))    nullify(this%solver_)
    if (associated(this%psolver_))   nullify(this%psolver_)
    if ( allocated(this%alpha_) )    deallocate(this%alpha_)
    if ( allocated (this%beta_) )    deallocate (this%beta_)
    if ( allocated    (this%c_) )    deallocate    (this%c_)
    call this%dealloc_tmp()        ! deallocate tmp arrays    
  end subroutine GExRKStepper_dtor


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Subroutine to initialize the stepper
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_impl(this, traits)
    class (GExRKStepper), intent(inout) :: this
    type(GStepperTraits), intent   (in) :: traits
    this%traits_ = traits
    ! Validate norder, nstage:
    if ( this%traits_%norder .le. 0 .or. this%traits_%nstage .le. 0 ) then
      stop 'GExRKStepper::init: Invalid norder/nstage'
    endif
  end subroutine init_impl


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Subroutine to initialize Butcher table
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_butcher(this)
    class(GExRKStepper), intent (inout) :: this
    this%busing_butcher_ = .true.
    ! When norder = nstage, we can use the following
    ! Butcher tableau formulation:
    if (allocated(this%alpha_)) deallocate(this%alpha_)
    if (allocated (this%beta_)) deallocate (this%beta_)
    if (allocated    (this%c_)) deallocate    (this%c_)
    allocate(this%alpha_(this%traits_%norder))
    allocate(this%c_    (this%traits_%norder))
    allocate(this%beta_ (this%traits_%norder, this%traits_%norder) )
    this%alpha_ = 0.0_GP
    this%beta_  = 0.0_GP
    this%c_     = 0.0_GP
    select case ( this%traits_%norder )
      case (1)
        this%alpha_(1) = 0.0_GP  ;  this%c_(1) = 1.0_GP
      case (2)
        this%alpha_  (1) = 0.0_GP;  this%c_(1) = 0.5_GP
        this%alpha_  (2) = 1.0_GP;  this%c_(2) = 0.5_GP
        this%beta_ (2,1) = 1.0_GP
      case (3)
        this%alpha_  (1) = 0.0_GP       ;  this%c_(1) = 0.25_GP
        this%alpha_  (2) = 1.0_GP/3.0_GP;  this%c_(2) =  0.0_GP
        this%alpha_  (3) = 2.0_GP/3.0_GP;  this%c_(3) = 0.75_GP
        this%beta_ (2,1) = 1.0_GP/3.0_GP
        this%beta_ (3,2) = 2.0_GP/3.0_GP
      case (4)
        this%alpha_  (1) = 0.0_GP       ;  this%c_(1) = 1.0_GP/6.0_GP
        this%alpha_  (2) = 0.5_GP       ;  this%c_(2) = 1.0_GP/3.0_GP
        this%alpha_  (3) = 0.5_GP       ;  this%c_(3) = 1.0_GP/3.0_GP
        this%alpha_  (4) = 1.0_GP       ;  this%c_(4) = 1.0_GP/6.0_GP
        this%beta_ (2,1) = 0.5_GP
        this%beta_ (3,2) = 0.5_GP
        this%beta_ (4,3) = 1.0_GP 
      case default
        stop 'GExRKStepper::init_butcher: Invalid norder/nstage'
    end select 
  end subroutine init_butcher

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Subroutine to initialize for norder != nstage methods
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_mixed(this)
    class(GExRKStepper), intent(inout) :: this
    stop 'GExRKStepper::init_mixed: GEXRK_MIXED not yet supported!'
  end subroutine init_mixed

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Subroutine to initialize for strong stability
  !! preserving methods
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_ssp(this)
    class(GExRKStepper), intent(inout) :: this
    stop 'GExRKStepper::init_ssp: GEXRK_SSP not yet supported!'
  end subroutine init_ssp


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Subroutine to allocate set of tmp arrays
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine alloc_tmp(this)
    use gstate_mod
    class(GExRKStepper), intent(inout) :: this

    call this%dealloc_tmp()    ! deallocate tmp arrays
    call GState_alloc(this%utmp_, this%traits_%nstate)    
    select case ( this%traits_%itype )
      case (GEXRK_BUTCHER)
        call this%init_butcher()
        call GCStateArr_alloc(this%K_,this%traits_%nstage,this%traits_%nstate)
      case (GEXRK_MIXED)
        call this%init_mixed()
      case (GEXRK_SSP)
        call this%init_ssp()
      case default
        stop 'GExRKStepper::init: Invalid stepper type'
    end select
  end subroutine alloc_tmp


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Subroutine to deallocate set of tmp pointer arrays
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine dealloc_tmp(this)
    class(GExRKStepper), intent(inout) :: this
    integer                            :: i

    if ( allocated( this%utmp_) ) call GState_dealloc (this%utmp_ )
    if ( allocated(this%putmp_) ) call GPState_dealloc(this%putmp_)
    if ( allocated(    this%K_) ) call GCStateArr_dealloc(this%K_ )
    if ( allocated(   this%pK_) ) call GPStateArr_dealloc(this%PK_)
  end subroutine dealloc_tmp


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Implementation function to take one step of fields
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine step_impl(this, time, uin, uf, dt, uout)
    implicit none

    class(GExRKStepper), intent(inout) :: this
    type   (GStateComp), intent(inout) :: uin(:), uf(:), uout(:)
    real      (kind=GP), intent   (in) :: time, dt
       
    if ( size(uin) .ne. this%traits_%nstate &
     .or.size(uout) .ne. this%traits_%nstate  ) then
      stop 'GExRKStepper::step: Inconsistent input state'
    endif

    select case ( this%traits_%itype )
      case (GEXRK_BUTCHER)
        call this%step_butcher(time, uin, uf, dt, uout)
      case (GEXRK_MIXED)
        call this%step_mixed  (time, uin, uf, dt, uout)
      case (GEXRK_SSP)
        call this%step_ssp    (time, uin, uf, dt, uout)
      case default
        stop 'GExRKStepper::step: Invalid stepper type'
    end select
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

    class(GExRKStepper), intent(inout)         :: this
    type   (GStateComp), intent(inout)         :: uin (:)
    type  (GPStateComp), intent(inout), target, allocatable :: upin(:), upout(:)
    real      (kind=GP), intent   (in)         :: time, dt

    if ( size (uin) .ne. this%traits_%nstate ) then
      stop 'GExRKStepper::pstep: Inconsistent input state'
    endif
    if ( size(upin) .ne. this%traits_%npstate &
    .or.size(upout) .ne. this%traits_%npstate  ) then
      stop 'GExRKStepper::pstep: Inconsistent particle state'
    endif
    if ( .not. this%traits_%dopart ) then
      stop 'GExRKStepper::pstep: The particle stepper was not initialized'
    endif
    if ( this%psolver_%hasfeedback_ ) then
      stop 'GExRKStepper::pstep: Particles are not one-way coupled'
    endif

    select case ( this%traits_%itype )
      case (GEXRK_BUTCHER)
        call this%pstep_butcher(time, uin, upin, dt, upout)
      case (GEXRK_MIXED)
        call this%pstep_mixed  (time, uin, upin, dt, upout)
      case (GEXRK_SSP)
        call this%pstep_ssp    (time, uin, upin, dt, upout)
      case default
        stop 'GExRKStepper::pstep: Invalid stepper type'
    end select
  end subroutine pstep_impl


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Implementation function to take one *coupled* step for
  !! particles + fields. Feedback forces are computed if 
  !! needed. Both particles and fields are evolved
  !! simultaneously in the substepping stages.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine cstep_impl(this, time, uin, upin, uf, dt, uout, upout)
    implicit none

    class(GExRKStepper), intent(inout)         :: this
    type   (GStateComp), intent(inout)         :: uin (:), uf(:), uout(:)
    type   (GStateComp), allocatable           :: fdbk(:)
    type  (GPStateComp), intent(inout), target, allocatable :: upin(:), upout(:)
    real      (kind=GP), intent   (in)         :: time, dt
       
    if ( size(uin) .ne. this%traits_%nstate &
     .or.size(uout) .ne. this%traits_%nstate  ) then
      stop 'GExRKStepper::cstep: Inconsistent input state'
    endif
    if ( size(upin) .ne. this%traits_%npstate &
     .or.size(upout) .ne. this%traits_%npstate  ) then
      stop 'GExRKStepper::cstep: Inconsistent particle state'
    endif
    if ( .not. this%traits_%dopart ) then
      stop 'GExRKStepper::cstep: The particle stepper was not initialized'
    endif

    select case ( this%traits_%itype )
      case (GEXRK_BUTCHER)
        call this%cstep_butcher(time, uin, upin, uf, dt, uout, upout)
      case (GEXRK_MIXED)
        call this%cstep_mixed  (time, uin, upin, uf, dt, uout, upout)
      case (GEXRK_SSP)
        call this%cstep_ssp    (time, uin, upin, uf, dt, uout, upout)
      case default
        stop 'GExRKStepper::cstep: Invalid stepper type'
    end select
  end subroutine cstep_impl


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Function to take one GEXRK_BUTCHER step of fields
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine step_butcher(this, time, uin, uf, dt, uout)
    use pseudospec_fluid
    use gstate_mod
    implicit none

    class(GExRKStepper), intent(inout) :: this
    type   (GStateComp), intent(inout) :: uin(:), uf(:), uout(:)
    real      (kind=GP), intent   (in) :: time, dt
    real      (kind=GP)                :: tt, eff_dt
    integer                            :: j,m,n
    ! alpha_ : time fractions for each stage
    ! beta_  : stage coefficient matrix 
    ! c_     : weights for final combination
    ! K_     : stage data
   
    ! Compute stage data:
    do m = 1, this%traits_%nstage
      this%utmp_ = uin                 ! set temp state
      do j = 1, m-1                    ! utmp = utmp + h beta K_j
        do n = 1, this%traits_%nstate  ! set utmp components
          call saxpby_c(this%utmp_(n)%ccomp, this%utmp_(n)%ccomp,    &
               1.0_GP, this%K_(j)%cstate(n)%ccomp, this%beta_(m,j)*dt)
        enddo
      enddo ! j-loop
      tt = time + this%alpha_(m) * dt  ! dudt called AFTER j-loop
      eff_dt    = this%alpha_(m) * dt
      call this%solver_%dudt(tt,this%utmp_,uf,eff_dt,this%K_(m)%cstate) 
    enddo ! stage m loop

    ! Combine stages to get step update:
    uout = uin
    do m = 1, this%traits_%nstage  
      do n = 1, this%traits_%nstate    ! uout = uout + h * c_m * K_m
        call saxpby_c(uout(n)%ccomp, uout(n)%ccomp, 1.0_GP,          &
             this%K_(m)%cstate(n)%ccomp, this%c_(m)*dt)
      enddo
    enddo ! m-loop
  end subroutine step_butcher


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Function to take one particle GEXRK_BUTCHER pstep
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine pstep_butcher(this, time, uin, upin, dt, upout)
    use class_GPartComm
    use pseudospec_fluid
    implicit none

    class(GExRKStepper), intent(inout) :: this
    type   (GStateComp), intent(inout) :: uin (:)
    type  (GPStateComp), intent(inout), target, allocatable :: upin(:), upout(:)
    real      (kind=GP), intent   (in) :: time, dt
    real      (kind=GP)                :: tt, eff_dt
    integer                            :: j,m,n,nparts
    ! alpha_ : time fractions for each stage
    ! beta_  : stage coefficient matrix 
    ! c_     : weights for final combination
    ! pK_    : particle stage data

    ! We allocate tmp arrays for particles here as we need the particle state size
    nparts = this%psolver_%partbuff_
    call GPState_alloc(this%putmp_,this%traits_%npstate,nparts)
    call GPStateArr_alloc(this%pK_,this%traits_%nstage,this%traits_%npstate,nparts)
    
    ! Compute stage data:
    do m = 1, this%traits_%nstage
      this%putmp_ = upin               ! set temp state
      do j = 1, m-1                    ! putmp = putmp + h beta pK_j
        eff_dt = this%beta_(m,j) * dt
        do n = 1, this%traits_%npstate ! set putmp components
          this%putmp_(n)%rcomp = this%putmp_(n)%rcomp+this%pK_(j)%rpstate(n)%rcomp*eff_dt
        enddo
      enddo ! j-loop
      tt = time + this%alpha_(m) * dt  ! dpdt called AFTER j-loop
      eff_dt    = this%alpha_(m) * dt
      call this%psolver_%dpdt(tt,this%solver_,uin,this%putmp_,eff_dt,this%pK_(m)%rpstate)
    enddo ! stage m loop

    ! Combine stages to get step update
    upout = upin
    do m = 1, this%traits_%nstage  
      eff_dt = this%c_(m) * dt
      do n = 1, this%traits_%npstate   ! upout = upout + h * c_m * pK_m
        upout(n)%rcomp = upout(n)%rcomp + this%pK_(m)%rpstate(n)%rcomp * eff_dt
      enddo
    enddo ! m-loop

    ! We now can deallocate the tmp arrays, and sync/resize upin, upout.
    ! Note that this is the only sync done, if the particles move too fast
    ! during the substepping stages (dt too long) this method may fail.
    call GPState_dealloc(this%putmp_)
    call GPStateArr_dealloc(this%pK_)
    call this%psolver_%end_stage(upin,upout)
  end subroutine pstep_butcher


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Function to take one part+field GEXRK_BUTCHER step
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine cstep_butcher(this, time, uin, upin, uf, dt, uout, upout)
    use class_GPartComm
    use pseudospec_fluid
    implicit none

    class(GExRKStepper), intent(inout) :: this
    type   (GStateComp), intent(inout) :: uin(:), uf(:), uout(:)
    type   (GStateComp), allocatable   :: fdbk(:)
    type  (GPStateComp), intent(inout), target, allocatable :: upin(:), upout(:)
    real      (kind=GP), intent   (in) :: time, dt
    real      (kind=GP)                :: tt, eff_dt
    integer                            :: j,m,n,ic,nparts
    ! alpha_ : time fractions for each stage
    ! beta_  : stage coefficient matrix 
    ! c_     : weights for final combination
    ! K_     : stage data
    ! pK_    : particle stage data
    ! fdbk   : feedback forces

    ! We allocate fdbk if we have feedback on the fluid
    if ( this%psolver_%hasfeedback_ ) then
      call GState_alloc(fdbk, this%traits_%nstate)
    endif

    ! We allocate tmp arrays for particles here as we need the particle state size
    nparts = this%psolver_%partbuff_
    call GPState_alloc(this%putmp_,this%traits_%npstate,nparts)
    call GPStateArr_alloc(this%pK_,this%traits_%nstage,this%traits_%npstate,nparts)
    
    ! Compute stage data:
    do m = 1, this%traits_%nstage
      this%utmp_  = uin                ! set temp states
      this%putmp_ = upin
      do j = 1, m-1
        do n = 1, this%traits_%nstate  ! set utmp  = utmp  + h beta K_j
          call saxpby_c(this%utmp_(n)%ccomp,  this%utmp_(n)%ccomp,   &
               1.0_GP, this%K_(j)%cstate(n)%ccomp,   this%beta_(m,j)*dt)
        enddo
        eff_dt  = this%beta_(m,j) * dt
        do n = 1, this%traits_%npstate ! set putmp = putmp + h beta pK_j
          this%putmp_(n)%rcomp = this%putmp_(n)%rcomp+this%pK_(j)%rpstate(n)%rcomp*eff_dt
        enddo
      enddo ! j-loop
      tt = time + this%alpha_(m) * dt  ! feedback and d/dt called AFTER j-loop
      eff_dt    = this%alpha_(m) * dt
      ! If needed we include the feedback of the particles in the fluid in uf
      if ( this%psolver_%hasfeedback_ ) then
        call this%psolver_%feedback(this%putmp_,fdbk)
        do ic = 1,this%traits_%nstate
          fdbk(ic)%ccomp = fdbk(ic)%ccomp + uf(ic)%ccomp
        end do
        call this%solver_%dudt(tt,this%utmp_,fdbk,eff_dt,this%K_(m)%cstate)
      else
        call this%solver_%dudt(tt,this%utmp_,  uf,eff_dt,this%K_(m)%cstate)
      end if
      call this%psolver_%dpdt(tt,this%solver_,uin,this%putmp_,eff_dt,this%pK_(m)%rpstate)
    enddo ! stage m loop
    
    ! Combine stages to get step update
    uout  = uin
    upout = upin
    do m = 1, this%traits_%nstage  
      do n = 1, this%traits_%nstate    ! uout  = uout  + h * c_m * K_m
        call saxpby_c(uout(n)%ccomp, uout(n)%ccomp, 1.0_GP,          &
             this%K_(m)%cstate(n)%ccomp, this%c_(m)*dt)
     enddo
     eff_dt = this%c_(m) * dt
     do n = 1, this%traits_%npstate   ! upout = upout + h * c_m * pK_m
        upout(n)%rcomp = upout(n)%rcomp + this%pK_(m)%rpstate(n)%rcomp * eff_dt
      enddo
    enddo ! m-loop
    
    ! We now can deallocate the tmp arrays, and sync/resize upin, upout.
    ! Note that this is the only sync done, if the particles move too fast
    ! during the substepping stages (dt too long) this method may fail.
    call GPState_dealloc(this%putmp_)
    call GPStateArr_dealloc(this%pK_)
    if ( this%psolver_%hasfeedback_ ) call GState_dealloc(fdbk)
    call this%psolver_%end_stage(upin,upout)    
  end subroutine cstep_butcher


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Function to take one GEXRK_MIXED step
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine step_mixed(this, time, uin, uf, dt, uout)
    implicit none

    class(GExRKStepper), intent   (in) :: this
    type   (GStateComp), intent(inout) :: uin(:), uf(:), uout(:)
    real      (kind=GP), intent   (in) :: time, dt
    stop 'GExRKStepper::step_mixed: GEXRK_MIXED not yet supported!'
  end subroutine step_mixed


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Function to take one particle GEXRK_MIXED pstep
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine pstep_mixed(this, time, uin, upin, dt, upout)
    implicit none

    class(GExRKStepper), intent   (in) :: this
    type   (GStateComp), intent(inout) :: uin (:)
    type  (GPStateComp), intent(inout), target, allocatable :: upin(:), upout(:)
    real      (kind=GP), intent   (in) :: time, dt
    stop 'GExRKStepper::pstep_mixed: GEXRK_MIXED not yet supported!'
  end subroutine pstep_mixed


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Function to take one part+field GEXRK_MIXED step
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine cstep_mixed(this, time, uin, upin, uf, dt, uout, upout)
    implicit none

    class(GExRKStepper), intent   (in) :: this
    type   (GStateComp), intent(inout) :: uin(:), uf(:), uout(:)
    type  (GPStateComp), intent(inout), target, allocatable :: upin(:), upout(:)
    real      (kind=GP), intent   (in) :: time, dt
    stop 'GExRKStepper::cstep_mixed: GEXRK_MIXED not yet supported!'
  end subroutine cstep_mixed


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Function to take one GEXRK_SSP step
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine step_ssp(this, time, uin, uf, dt, uout)
    implicit none

    class(GExRKStepper), intent   (in) :: this
    type   (GStateComp), intent(inout) :: uin(:), uf(:), uout(:)
    real      (kind=GP), intent   (in) :: time, dt
    stop 'GExRKStepper::step_ssp: GEXRK_SSP not yet supported!'
  end subroutine step_ssp

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Function to take one particle GEXRK_SSP pstep
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine pstep_ssp(this, time, uin, upin, dt, upout)
    implicit none

    class(GExRKStepper), intent   (in) :: this
    type   (GStateComp), intent(inout) :: uin (:)
    type  (GPStateComp), intent(inout), target, allocatable :: upin(:), upout(:)
    real      (kind=GP), intent   (in) :: time, dt
    stop 'GExRKStepper::pstep_ssp: GEXRK_SSP not yet supported!'
  end subroutine pstep_ssp


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Function to take one part+field GEXRK_SSP step
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine cstep_ssp(this, time, uin, upin, uf, dt, uout, upout)
    implicit none

    class(GExRKStepper), intent   (in) :: this
    type   (GStateComp), intent(inout) :: uin(:), uf(:), uout(:)
    type  (GPStateComp), intent(inout), target, allocatable :: upin(:), upout(:)
    real      (kind=GP), intent   (in) :: time, dt
    stop 'GExRKStepper::cstep_ssp: GEXRK_SSP not yet supported!'
  end subroutine cstep_ssp

end module gexrk_stepper_mod
