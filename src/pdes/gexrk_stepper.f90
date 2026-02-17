! =====================================================================
! NAME       : gexrk.f90
!  
!              Performs time stepping using explicit RK 
!              with user-specified order and number of
!              stages
!              
! INPUT FILE : Stepper looks for a "&stepper" namelist with:
!                itype   : Stepper type:
!                    =1 => Butcher type (norder = nstage)
!                    =2 => Mixed type (norder != nstage; not available)
!                    =3 => SSP type (strong stability preserving; not available)
!                norder  : Stepper order (1-4 currently)
!                nstage  : No. stepper stages (1-4 currently)
!
!              Currently stepper supports only explicit RK
!              schemes where norder=nstage, but the class 
!              is designed to eventually accommodate both 
!              GEXRK_MIXED (itype=2), and GEXRK_SSP (itype=3) types.
!
! DATE       : 2/2/26 (DLR)
! =====================================================================

module gexrk_mod
  use iso_c_binding
  use class_GWorkspace3D
  use gstate_mod
  use stepperbase_mod


  implicit none

  ! Define callback function interface:
! abstract interface
!    subroutine dudt_interface(this, time, uin, uf, dt, dudt)   
!      use gstate_mod
!      import :: GExRKStepper
!      class(GExRKStepper), intent   (in)         :: this
!      real      (kind=GP), intent   (in)         :: time, dt
!      type   (GStateComp), intent(inout), target :: uin(:),uf(:)
!      type   (GStateComp), intent(inout)         :: dudt(:) 
!    end subroutine dudt_interface
! end interface

  ! ================= Global parameters ===============================
  enum, bind(c)
    enumerator :: GEXRK_BUTCHER = 1, GEXRK_MIXED = 2, GEXRK_SSP = 3
  end enum

  ! ================= Stepper ==========================================
  ! Define class:
  type, extends(StepperBase) :: GExRKStepper
    ! Member data:
    type(GWorkspace), pointer     :: workspace_
    logical                       :: binit_=.false. ! is initialized?
    logical                       :: busing_butcher_=.true. ! using Butcher tableau?
!   integer                       :: myrank_   ! MPI rank
!   integer                       :: nprocs_   ! MPI rank 
    type(GStepperTraits)          :: traits_
    real   (kind=GP),        , allocatable, &
                               dimension    (:) :: alpha_, c_
    real   (kind=GP),        , allocatable, &
                               dimension  (:,:) :: beta_
    type   (GCState)         , allocatable, &
                               dimension    (:) :: K_
    type   (GPState),        , allocatable, &
                               dimension    (:) :: pK_
    type(GStateComp)         , allocatable, &
                               dimension    (:) :: utmp_
    type(GPStateComp)         , allocatable, &
                               dimension    (:) :: putmp_

    procedure (dudt_interface), pointer, nopass :: callback_  => null()
    procedure(pdudt_interface), pointer, nopass :: pcallback_ => null()

  CONTAINS
    procedure, public :: init                        ! initialize
    procedure, public :: set_callback  => set_callback_impl
                                                     ! set RHS callback method
    procedure, public :: set_pcallback => set_callback_impl
                                                     ! set RHS part callback method
    interface step
      procedure, public :: step         =>  step_impl, pstep_impl  ! take one timestep
    end interface

    procedure, public :: Stepper_ctor =>  GEXRKStepper_ctor 
    final             :: GExRKStepper_dtor

    procedure, private:: init_butcher   ! initialize Butcher table
    procedure, private:: init_mixed     ! initialize mixed type
    procedure, private:: init_ssp       ! initialize for SSP
    procedure, private:: step_butcher   ! take a Butcher step
    procedure, private:: step_mixed     ! take a mixed step
    procedure, private:: step_ssp       ! take a SSP step
    procedure, private:: pstep_butcher  ! take a part+field Butcher step
    procedure, private:: pstep_mixed    ! take a part+field mixed step
    procedure, private:: pstep_ssp      ! take a part+field SSP step
!   procedure, private:: set_tmp        ! set data from tmp pool
!   procedure, private:: free_tmp       ! free data from tmp pool
    procedure, private:: alloc_tmp      ! allocate tmp arrays
    procedure, private:: dealloc_tmp    ! deallocate tmp arrays
  end type GExRKStepper

CONTAINS

  ! ===================================================================
  ! Stepper-specific methods
  ! ===================================================================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Constructor
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine GExRKStepper_ctor(this, traits, workspace, nparts)
    use iovar
    class (GExRKStepper), intent(inout)          :: this
    type(GStepperTraits), intent(inout), target  :: traits
    type  (GWorkspace)  , intent(inout), target  :: workspace
    integer             , intent(inout), optional:: nparts

    this%workspace_ => workspace
    if (.not. associated(this%workspace_)) then
      stop 'GExRKStepper::GExRKStepper_ctor: Worskpace not associated'
    endif
    
    ! Call StepperBase constructor:
    this%GStepperBase%GStepper_ctor_interface(traits, workspace, nparts)

    call this%init(traits)
  end subroutine GExRKStepper_ctor


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Destructor
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine GExRKStepper_dtor(this) 
    type  (GExRKStepper), intent(inout) :: this

    if ( associated(this%workspace_) )   nullify(this%workspace_)
    if ( allocated(this%alpha_) ) deallocate(this%alpha_)
    if ( allocated (this%beta_) ) deallocate (this%beta_)
    if ( allocated    (this%c_) ) deallocate    (this%c_)


    call this%dealloc_tmp() ! deallocate tmp arrays
    
  end subroutine GExRKStepper_dtor


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Subroutine to initialize the stepper
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init(this, traits)
!   use commtypes
    class (GExRKStepper), intent (inout) :: this
    type(GStepperTraits), intent    (in) :: traits
    interger                             :: i

    this%traits_ = traits

    ! Validate norder, nstage:
    if ( this%traits_%norder .le. 0 .or. this%traits_%nstage .le. 0 ) then
      stop 'GExRKStepper::init: Invalid norder/nstage'
    endif

    this%traits_%itype  = itype
    this%traits_%norder = norder
    this%traits_%nstage = nstage

    call MPI_COMM_SIZE(MPI_COMM_WORLD,this%nprocs_,ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,this%myrank_,ierr)


    ! Not sure if it's best to set tmp from pool
    ! for lifetime of this object, or if we
    ! should make following call each time
    ! step method is called. For now, we keep
    ! it while object is in scope by setting it here:
!   call this%set_tmp()

  end subroutine init

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Subroutine to initialize Butcher table
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_butcher(this)
    class  (GExRKStepper), intent (inout) :: this

    this%busing_butcher_ = .true.

    ! When norder = nstage, we can use the following
    ! Butcher tableau formulation:
    if (allocated(this%alpha_)) deallocate(this%alpha_)
    if (allocated (this%beta_)) deallocate (this%beta_)
    if (allocated    (this%c_)) deallocate    (this%c_)

    allocate(this%alpha_(this%traits_.norder))
    allocate(this%c_    (this%traits_.norder))
    allocate(this%beta_ (this%traits_.norder, &
                         this%traits_.norder) )

    this%alpha_ = 0.0_GP
    this%beta_  = 0.0_GP
    this%c_     = 0.0_GP

    select case ( this%traits_%norder )
      case 1:
        this%alpha_(1) = 0.0_GP;  this%c_(1) = 1.0_GP; 
      case 2:
        this%alpha_  (1) = 0.0_GP;  this%c_(1) = 0.5_GP; 
        this%alpha_  (2) = 1.0_GP;  this%c_(2) = 0.5_GP; 
        this%beta_ (2,1) = 1.0_GP; 
      case 3:
        this%alpha_  (1) = 0.0_GP       ;  this%c_(1) = 0.25_GP; 
        this%alpha_  (2) = 1.0_GP/3.0_GP;  this%c_(2) = 0.0_GP; 
        this%alpha_  (3) = 2.0_GP/3.0_GP;  this%c_(3) = 0.75_GP; 
        this%beta_ (2,1) = 1.0_GP/3.0_GP; 
        this%beta_ (3,2) = 2.0_GP/3.0_GP; 
      case 4:
        this%alpha_  (1) = 0.0_GP       ;  this%c_(1) = 1.0_GP/6.0_GP; 
        this%alpha_  (2) = 0.5_GP       ;  this%c_(2) = 1.0_GP/3.0_GP;
        this%alpha_  (3) = 0.5_GP       ;  this%c_(3) = 1.0_GP/3.0_GP;
        this%alpha_  (4) = 1.0_GP       ;  this%c_(4) = 1.0_GP/6.0_GP;
        this%beta_ (2,1) = 0.5_GP; 
        this%beta_ (3,2) = 0.5_GP
        this%beta_ (4,3) = 1.0_GP; 
      case default:
        stop 'GExRKStepper::init_butcher: Invalid norder/nstage'
    end select 

  end subroutine init_butcher

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Subroutine to initialize for norder != nstage methods
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_mixed(this)
    class  (GExRKStepper), intent (inout) :: this

    stop 'GExRKStepper::init_mixed: GEXRK_MIXED not yet supported!'
  end subroutine init_mixed

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Subroutine to initialize for strong stability
  !! preserving methods
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_ssp(this)
    class  (GExRKStepper), intent (inout) :: this

    stop 'GExRKStepper::init_ssp: GEXRK_SSP not yet supported!'
  end subroutine init_ssp

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !! Subroutine to allow 'volatile' set of tmp
! !! data from workspace
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine set_tmp(this)
!   class  (GExRKStepper), intent (inout) :: this
!   logical                               :: bret
!   integer                               :: istage,icomp

!   if ( this%traits_%nstate .le. 0 ) then
!     stop 'GExRKStepper::set_tmp: Invalid nstate'
!   endif 

!   if ( this%traits_%norder .le. 0 .or. this%traits_%nstage .le. 0 ) then
!     stop 'GExRKStepper::set_tmp: Invalid norder/nstage'
!   endif 


!   if ( this%busing_butcher_ ) then
!     if ( .not. allocated(this%K_) ) then
!       stop 'GExRKStepper::set_tmp: volatile not allocated'
!     endif
!     do istage = 1, this%traits_%nstage   
!       do icomp = 1, this%traits_%nstate
!         CALL this%workspace_%get_complex_tmp(this%K_(istage)%cstate(icomp),bret)
!         if ( .not. bret  ) then
!           stop 'GExRKStepper::set_tmp: Workspace field failure'
!         endif
!       enddo
!     enddo
!   endif

!   ! Set particle tmp space, if using particles:
!   if ( this%doparts_ ) then
!     if ( this%traits_%npstate .le. 0 ) then
!       stop 'GExRKStepper::set_tmp: Invalid npstate'
!     endif 
!     if ( .not. allocated(this%pK_) ) then
!       stop 'GExRKStepper::set_tmp: volatile not allocated'
!     endif
!     do istage = 1, this%traits_%nstage   
!       do icomp = 1, this%traits_%npstate
!         CALL this%workspace_%get_pcomp_tmp(this%pK_(istage)%rpstate(icomp),bret)
!         if ( .not. bret  ) then
!           stop 'GExRKStepper::set_tmp: Workspace particle failure'
!         endif
!       enddo
!     enddo
!   endif
!   

! end subroutine set_tmp


! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !! Subroutine to free 'volatile' set of tmp
! !! data in workspace
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine free_tmp(this)
!   class  (GExRKStepper), intent (inout) :: this
!   logical                               :: bret
!   integer                               :: istage,istate

!   do istage = 1, this%traits_%nstage   
!     do istate = 1, this%traits_%nstate   
!       CALL this%workspace_%free_complex_tmp(this%K_(istage)%cstate(istate))
!     enddo
!   enddo

!   do istage = 1, this%traits_%nstage   
!     do istate = 1, this%traits_%npstate   
!       CALL this%workspace_%free_complex_tmp(this%pK_(istage)%rpstate(istate))
!     enddo
!   enddo

! end subroutine free_tmp


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Subroutine to allocate set of tmp arrays
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine alloc_tmp(this)
    class  (GExRKStepper), intent (inout) :: this
    logical                               :: bret
    integer                               :: i


!   call this%free_tmp()    ! free workspace pointers
    call this%dealloc_tmp() ! deallocate tmp pointer arrays

    if ( allocated(this%utmp_) ) then
      GState_dealloc(this%utmp_)
    endif
    call GState_alloc(this%utmp_, this%traits_%nstate)

    if ( this%doparts_ ) then
      if ( allocated(this%putmp_) ) deallocate(this%putmp_)
      GPState_dealloc(this%putmp_)
      call GPState_alloc(this%putmp_, this%traits_%npstate)
    endif

    select case ( this%traits_%itype )
      case GEXRK_BUTCHER:
        call this%init_butcher()
        allocate(this%K_(this%traits_%nstage))
        do i = 1, this%traits_%nstage
          call  GState_alloc(this%K_(i)%cstate, this%traits_%nstate)
        enddo
        if ( this%doparts_ ) then
          allocate(this%pK_(this%traits_%nstage))
          do i = 1, this%traits_%nstage
            call  GPState_alloc(this%pK_(i)%rpstate, this%traits_%npstate)
          enddo
        endif
      case GEXRK_MIXED:
        call this%init_mixed()
      case GEXRK_SSP:
        call this%init_ssp()
      case default:
        stop 'GExRKStepper::init: Invalid stepper type'
    end select


  end subroutine alloc_tmp

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Subroutine to deallocate set of tmp pointer arrays
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine dealloc_tmp(this)
    class  (GExRKStepper), intent (inout) :: this
    logical                               :: bret
    integer                               :: i

!   this%free_tmp(); ! free workspace
    call  GState_dealloc (this%utmp_)
    call  GPState_dealloc(this%putmp_)

    if ( allocated(this%K_) ) then
      do i = 1, this%traits_%nstage   
        if ( allocated(this%K_(i)) ) deallocate(this%K_(i))
      enddo
      deallocate(this%K_)
    enddo

    if ( allocated(this%pK_) ) then
      do i = 1, this%traits_%nstage   
        if ( allocated(this%pK_(i)) ) deallocate(this%pK_(i))
      enddo
      deallocate(this%pK_))
    enddo

  end subroutine dealloc_tmp


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Subroutine to set RHS callback function
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine set_callback_impl(this, fcn_callback)
    class  (GExRKStepper), intent  (inout) :: this
    procedure(callback_interface), pointer :: fcn_callback

    this%callback_ => fcn_callback

  end subroutine set_callback_impl


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Subroutine to set RHS part+field callback function
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine set_pcallback_impl(this, fcn_callback)
    class  (GExRKStepper), intent  (inout) :: this
    procedure(pcallback_interface), pointer :: fcn_callback

    this%pcallback_ => fcn_callback

  end subroutine set_callback_impl


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Implementation function to take one step
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine step_impl(this, time, uin, uf, dt, uout)
!$  use threads
    implicit none

    class (GExRKStepper), intent   (in) :: this
    type    (GStateComp), intent(inout) :: uin(:), uf(:), uout(:)
    real       (kind=GP), intent   (in) :: time, dt
    logical                             :: bret
       
    if ( size(uin) .ne. this%traits_%nstate &
     .or.size(uout) .ne. this%traits_%nstate  ) then
      stop 'GExRKStepperi::step: Inconsistent input state'
    endif
    
    if ( .not. associated(this%callback_) ) then
      stop 'GExRKStepperi::step: RHS callback function not set'
    endif

    select case ( this%traits_%itype )
      case GEXRK_BUTCHER:
        call this%step_butcher(time, uin, uf, dt, uout)
      case GEXRK_MIXED:
        call this%step_mixed(time, uin, uf, dt, uout)
      case GEXRK_SSP:
        call this%step_ssp(time, uin, uf, dt, uout)
      case default:
        stop 'GExRKStepper::step: Invalid stepper type'
    end select

  end subroutine step_impl

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Implementation function to take one part+field step
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine pstep_impl(this, time, uin, puin, uf, dt, uout, puout)
!$  use threads
    implicit none

    class (GExRKStepper), intent   (in) :: this
    type    (GStateComp), intent(inout) :: uin(:), uf(:), uout(:)
    type   (GPStateComp), intent(inout) :: puin(:), puout(:)
    real       (kind=GP), intent   (in) :: time, dt
    logical                             :: bret
       
    if ( size(uin) .ne. this%traits_%nstate &
     .or.size(uout) .ne. this%traits_%nstate  ) then
      stop 'GExRKStepperi::step: Inconsistent input state'
    endif
    if ( size(puin) .ne. this%traits_%npstate &
     .or.size(puout) .ne. this%traits_%npstate  ) then
      stop 'GExRKStepperi::step: Inconsistent particle state'
    endif
    
    if ( .not. associated(this%pcallback_) ) then
      stop 'GExRKStepperi::pstep: RHS pcallback function not set'
    endif

    select case ( this%traits_%itype )
      case GEXRK_BUTCHER:
        call this%step_butcher(time, uin, puin, uf, dt, uout, puout)
      case GEXRK_MIXED:
        call this%step_mixed(time, uin, puin, uf, dt, uout, puout)
      case GEXRK_SSP:
        call this%step_ssp(time, uin, puin, uf, dt, uout, puout)
      case default:
        stop 'GExRKStepper::step: Invalid stepper type'
    end select

  end subroutine pstep_impl


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Function to take one GEXRK_BUTCHER step
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine step_butcher(this, time, uin, uf, dt, uout)
!$  use threads
    implicit none

    class (GExRKStepper), intent   (in) :: this
    type    (GStateComp), intent(inout) :: uin(:), uf(:), uout(:)
!   complex(kind=GP), pointer, &
!                      dimension(:,:,:) :: sum
    real       (kind=GP), intent   (in) :: time, dt
    real       (kind=GP)                :: tt
    logical                             :: bret
    integer                             :: j,m,n

!   CALL this%workspace_%get_complex_tmp(sum,bret)

    ! alpha_ : time fractions for each stage
    ! beta_  : stage coefficient matrix 
    ! c_     : weights for final combination
    ! K_     : stage data: 
    ! pK_    : particle stage data: 
   
    ! Compute stage data:
    do m = 1, this%traits_%nstage

      do n = 1, this%traits_%nstate  ! set temp state
        this%utmp_(n) = uin(n)
      enddo
      do j = 1, m-1  ! utmp = utmp + h beta K_j
        do n = 1, this%traits_%nstate  ! set utmp
          call saxpby_c(this%utmp_(n), this%utmp_(n), 1.0_GP, this%K_(j)%cstate(n), this%beta_(m,j)*dt)
        enddo
        tt = time + this%alpha_(m) * dt
        this%callback_( tt, this%utmp_, uf, dt, K_(m)%cstate ); 
      enddo ! j-loop

    enddo ! stage m loop

    ! Combine stages to get step update:
    do n = 1, this%traits_%nstate  ! uout = uin
      uout(n) = uin(n)
    enddo

    do m = 1, this%traits_%nstage  

      do n = 1, this%traits_%nstate ! uout = uout + h * c_ * K:
        call saxpby_c(uout(n), uout(n), 1.0_GP, this%K_(m)%cstate(n), this%c_(m,j)*dt)
      enddo

    enddo ! m-loop

!   CALL this%workspace_%free_complex_tmp(sum)

  end subroutine step_butcher

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Function to take one part+field GEXRK_BUTCHER step
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine pstep_butcher(this, time, uin, puin, uf, dt, uout, puout)
!$  use threads
    implicit none

    class (GExRKStepper), intent   (in) :: this
    type    (GStateComp), intent(inout) :: uin(:), uf(:), uout(:)
    type   (GPStateComp), intent(inout) :: puin(:), puout(:)
    complex(kind=GP), pointer, &
                       dimension(:,:,:) :: sum
    real       (kind=GP), intent   (in) :: time, dt
    real       (kind=GP)                :: tt
    logical                             :: bret
    integer                             :: j,m,n


    ! alpha_ : time fractions for each stage
    ! beta_  : stage coefficient matrix 
    ! c_     : weights for final combination
    ! K_     : stage data: 
    ! pK_    : particle stage data: 
   
    ! Compute stage data:
    do m = 1, this%traits_%nstage

      do n = 1, this%traits_%nstate  ! set temp state
        this%utmp_(n) = uin(n)
      enddo
      do n = 1, this%traits_%npstate  ! set ptemp state
        this%putmp_(n) = puin(n)
      enddo
      do j = 1, m-1  ! utmp = utmp + h beta K_j
        do n = 1, this%traits_%nstate  ! set utmp
          call saxpby_c(this%utmp_(n), this%utmp_(n), 1.0_GP, this%K_(j)%cstate(n), this%beta_(m,j)*dt)
        enddo
        do n = 1, this%traits_%npstate  ! set utmp
          call saxpby_c(this%putmp_(n), this%putmp_(n), 1.0_GP, this%pK_(j)%rpstate(n), this%beta_(m,j)*dt)
        enddo
        tt = time + this%alpha_(m) * dt
        this%pcallback_( tt, this%utmp_, this%putmp_, uf, dt, K_(m)%cstate, pK_(m)%rpstate); 
      enddo ! j-loop

    enddo ! stage m loop

    ! Combine stages to get step update:
    do n = 1, this%traits_%nstate  ! uout = uin
      uout(n) = uin(n)
    enddo
    do n = 1, this%traits_%npstate  ! uout = uin
      puout(n) = puin(n)
    enddo

    do m = 1, this%traits_%nstage  

      do n = 1, this%traits_%nstate ! uout = uout + h * c_ * K:
        call saxpby_c(uout(n), uout(n), 1.0_GP, this%K_(m)%cstate(n), this%c_(m,j)*dt)
      enddo
      do n = 1, this%traits_%npstate ! uout = uout + h * c_ * K:
        call saxpby_c(puout(n), puout(n), 1.0_GP, this%pK_(m)%rpstate(n), this%c_(m,j)*dt)
      enddo

    enddo ! m-loop


  end subroutine pstep_butcher

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Function to take one GEXRK_MIXED step
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine step_mixed(this, time, uin, uf, dt, uout)
!$  use threads
    implicit none

    class (GExRKStepper), intent   (in) :: this
    type    (GStateComp), intent(inout) :: uin(:), uf(:), uout(:)
    real       (kind=GP), intent   (in) :: time, dt
    logical                             :: bret

    stop 'GExRKStepper::step_mixed: GEXRK_MIXED not yet supported!'
  end subroutine step_mixed

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Function to take one part+field GEXRK_MIXED step
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine pstep_mixed(this, time, uin, puin, uf, dt, uout, puout)
!$  use threads
    implicit none

    class (GExRKStepper), intent   (in) :: this
    type    (GStateComp), intent(inout) :: uin(:), uf(:), uout(:)
    type   (GPStateComp), intent(inout) :: puin(:), puout(:)
    real       (kind=GP), intent   (in) :: time, dt
    logical                             :: bret

    stop 'GExRKStepper::pstep_mixed: GEXRK_MIXED not yet supported!'
  end subroutine pstep_mixed


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Function to take one GEXRK_SSP step
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine step_ssp(this, time, uin, uf, dt, uout)
!$  use threads
    implicit none

    class (GExRKStepper), intent   (in) :: this
    type    (GStateComp), intent(inout) :: uin(:), uf(:), uout(:)
    real       (kind=GP), intent   (in) :: time, dt
    logical                             :: bret

    stop 'GExRKStepper::step_ssp: GEXRK_SSP not yet supported!'
  end subroutine step_ssp

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Function to take one part+field GEXRK_SSP step
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine pstep_ssp(this, time, uin, puin, uf, dt, uout, puout)
!$  use threads
    implicit none

    class (GExRKStepper), intent   (in) :: this
    type    (GStateComp), intent(inout) :: uin(:), uf(:), uout(:)
    type   (GPStateComp), intent(inout) :: puin(:), puout(:)
    real       (kind=GP), intent   (in) :: time, dt
    logical                             :: bret

    stop 'GExRKStepper::pstep_ssp: GEXRK_SSP not yet supported!'
  end subroutine pstep_ssp


end module gexrk_mod
