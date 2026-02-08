! =====================================================================
! NAME       : gexrk.f90
!              using explicit RK of specified order and number of
!              stages
!
! INPUT FILE : Stepper looks for a "&stepper" namelist with:
!                norder  : Stepper order (1-4 currently)
!                nstage  : No. stepper stages (1-4 currently)
!
!              Currently stepper supports only explicit RK
!              schemes where norder=nstage, but this may change
!              in the futurel this object is designed up to 
!              accommodate this case eventually.
!
! DATE       : 2/2/26 (DLR)
! =====================================================================

module gexrk_mod
  use iso_c_binding
  use class_GWorkspace3D
  use gstate_mod
  use equationbase_mod


  implicit none

  ! Define callback function interface:
  abstract interface
     subroutine dudt_interface(this, time, uin, uf, dt, dudt)   
       use gstate_mod
       import :: GExRKStepper
       class(GExRKStepper), intent   (in)         :: this
       real      (kind=GP), intent   (in)         :: time, dt
       type   (GStateComp), intent(inout), target :: uin(:),uf(:)
       type   (GStateComp), intent(inout)         :: dudt(:) 
     end subroutine dudt_interface
  end interface

  ! ================= Global parameters ===============================
  enum, bind(c)
    enumerator :: GEXRK_BUTCHER = 1, GEXRK_MIXED = 2, GEXRK_SSP = 3
  end enum

  ! ================= Stepper traits ===================================
  type, public  :: GEXRKTraits
!   logical       :: bSSP       = .FALSE. ! strong stability flag?
!   logical       :: blowstore  = .FALSE. ! low storage version, if possible?
    integer       :: itype      = GEXRK_BUTCHER
    integer       :: norder     = 2  
    integer       :: nstage     = 2 
  end type


  ! ================= Stepper ==========================================
  ! Define class:
  type :: GExRKStepper
    ! Member data:
    type(GWorkspace), pointer     :: workspace_
    logical                       :: binit_=.false. ! is initialized?
    logical                       :: busing_butcher_=.true. ! using Butcher tableau?
    integer                       :: nstate_=0 ! no. state members
    integer                       :: myrank_   ! MPI rank
    integer                       :: nprocs_   ! MPI rank 
    type  (GEXRKTraits)           :: traits_
    real   (kind=GP),        , allocatable, &
                               dimension    (:) :: alpha_, c_
    real   (kind=GP),        , allocatable, &
                               dimension  (:,:) :: beta_
    type   (GCState), pointer, allocatable, &
                               dimension    (:) :: K_

    procedure(dudt_interface), pointer, nopass :: callback_ => null()
  CONTAINS
    procedure, public :: set_order     ! set order, nstages
    procedure, public :: init          ! initialize
    procedure, public :: set_callback  ! set RHS callback method
    procedure, public :: step          ! take one timestep
    procedure, public :: GEXRKStepper_ctor 
    final             :: GExRKStepper_dtor

    procedure, private:: init_butcher   ! initialize Butcher table
    procedure, private:: init_mixed     ! initialize mixed type
    procedure, private:: init_ssp       ! initialize for SSP
    procedure, private:: step_butcher   ! take a Butcher step
    procedure, private:: step_mixed     ! take a mixed step
    procedure, private:: step_ssp       ! take a SSP step
    procedure, private:: set_tmp        ! set data from tmp pool
    procedure, private:: free_tmp       ! free data from tmp pool
  end type GExRKStepper

CONTAINS

  ! ===================================================================
  ! Stepper-specific methods
  ! ===================================================================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Constructor
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine GExRKStepper_ctor(this, infile, workspace, nstate)
    use iovar
    class  (GExRKStepper), intent(inout)     :: this
    type(GWorkspace) , intent(inout), target :: workspace
    character(len=*) , intent   (in)         :: infile

    this%workspace_ => workspace
    if (.not. associated(this%workspace_)) then
      stop 'GExRKStepper::GExRKStepper_ctor: Worskpace not associated'
    endif

    call this%init(infile, nstate)
  end subroutine GExRKStepper_ctor


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Destructor
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine GExRKStepper_dtor(this) 
    type  (GExRKStepper), intent(inout) :: this

    if (associated(this%workspace_))   nullify(this%workspace_)
    if (allocated(this%alpha_)) deallocate(this%alpha_)
    if (allocated (this%beta_)) deallocate (this%beta_)
    if (allocated    (this%c_)) deallocate    (this%c_)

    call this%free_tmp() ! must call this before deleting K_:
    if (allocated    (this%K_)) deallocate    (this%K_)
    
  end subroutine GExRKStepper_dtor


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Subroutine to initialize the stepper
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init(this, infile, nstate)
!   use commtypes
    class  (GExRKStepper), intent (inout) :: this
    integer          , intent   (in)      :: nstate
    character(len=*) , intent   (in)      :: infile

    ! Temporary data to read from namelists:
    integer                    :: itype,norder,nstage
    integer                    :: ierr
    
    ! Required namelists:
    namelist/ stepper    / norder, nstage

    if ( nstate .le. 0 ) then
      stop 'GExRKStepper::init: Invalid nstate'
    endif
    this%nstate_    = nstate    ! no. state components
    this%infile_    = infile    ! input file

    ! Get stepper traits from input file:
    itype    = GEXRK_BUTCHER
    norder   = 2
    nstage   = 2

    if ( this%myrank_ .eq. 0 ) then
      open(1,file=this%infile_,status='unknown',form="formatted")
      read(1,NML=HD)
      close(1)
    endif
    call mpi_bcast(itype    ,1 ,MPI_INTEGER,    0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(norder   ,1 ,MPI_INTEGER,    0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(nstage   ,1 ,MPI_INTEGER,    0,MPI_COMM_WORLD,ierr)
    ! Validate norder, nstage:
    if ( norder .le. 0 .or. nstage .le. 0 ) then
      stop 'GExRKStepper::init: Invalid norder/nstage'
    endif

    this%traits_%itype  = itype
    this%traits_%norder = norder
    this%traits_%nstage = nstage

    call MPI_COMM_SIZE(MPI_COMM_WORLD,this%nprocs_,ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,this%myrank_,ierr)

    select case ( this%traits_%itype )
      case GEXRK_BUTCHER:
        this%busing_butcher_ = .true.
        call this%init_butcher()
      case GEXRK_MIXED:
        call this%init_mixed()
      case GEXRK_SSP:
        call this%init_ssp()
      case default:
        stop 'GExRKStepper::init: Invalid stepper type'
    end select

    if (allocated(this%K_)) deallocate(this%K_)
    allocate(this%K_(this%nstage_))

    ! Not sure if it's best to set tmp from pool
    ! for lifetime of this object, or if we
    ! should make following call each time
    ! step method is called. For now, we keep
    ! it while object is in scope by setting it here:
    call this%set_tmp()

  end subroutine init

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Subroutine to initialize Butcher table
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_butcher(this)
    class  (GExRKStepper), intent (inout) :: this

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

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Subroutine to allow 'volatile' set of tmp
  !! data from workspace
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine set_tmp(this)
    class  (GExRKStepper), intent (inout) :: this
    logical                               :: bret
    integer                               :: istage,istate

    if ( this%nstate_ .le. 0 ) then
      stop 'GExRKStepper::set_tmp: Invalid nstate'
    endif 

    if ( this%norder_ .le. 0 .or. this%nstage_ .le. 0 ) then
      stop 'GExRKStepper::set_tmp: Invalid norder/nstage'
    endif 


    if ( this%busing_butcher_ ) then
      if ( .not. allocated(this%K_) ) then
        stop 'GExRKStepper::set_tmp: volatile not allocated'
      endif
      do istage = 1, this%nstage_   
        do istate = 1, this%nstate_   
          CALL this%workspace_%get_complex_tmp(this%K_(istage)(istate),bret)
          if ( .not. bret  ) then
            stop 'GExRKStepper::set_tmp: Workspace failure'
          endif
        enddo
      enddo
    endif

  end subroutine set_tmp


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Subroutine to free 'volatile' set of tmp
  !! data in workspace
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine free_tmp(this)
    class  (GExRKStepper), intent (inout) :: this
    logical                               :: bret
    integer                               :: istage,istate

    do istage = 1, this%nstage_   
      do istate = 1, this%nstate_   
        CALL this%workspace_%free_complex_tmp(this%K_(istage)(istate))
      enddo
    enddo

  end subroutine free_tmp

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Subroutine to set RHS callback function
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine set_callback(this, fcn_callback)
    class  (GExRKStepper), intent  (inout) :: this
    procedure(callback_interface), pointer :: fcn_callback

    this%callback_ => fcn_callback

  end subroutine set_callback


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Function to take one step
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine step(this, time, uin, uf, dt, uout)
!$  use threads
    implicit none

    class (GExRKStepper), intent   (in) :: this
    type    (GStateComp), intent(inout) :: uin(:), uf(:), uout(:)
    real       (kind=GP), intent   (in) :: time, dt
    logical                             :: bret
       
    if ( size(uin) .ne. this%nstate_ &
     .or.size(uout) .ne. this%nstate_  ) then
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

  end subroutine step


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Function to take one GEXRK_BUTCHER step
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine step_butcher(this, time, uin, uf, dt, uout)
!$  use threads
    implicit none

    class (GExRKStepper), intent   (in) :: this
    type    (GStateComp), intent(inout) :: uin(:), uf(:), uout(:)
    complex(kind=GP), pointer, &
                       dimension(:,:,:) :: sum
    real       (kind=GP), intent   (in) :: time, dt
    real       (kind=GP)                :: tt
    logical                             :: bret
    integer                             :: j,m,n

    CALL this%workspace_%get_complex_tmp(sum,bret)

    ! alpha_ : time fractions for each stage
    ! beta_  : stage coefficient matrix 
    ! c_     : weights for final combination
    ! K_     : stage data: 
   
    ! Compute stage data:
    do m = 1, this%nstage_

      do n = 1, this%nstate_  ! set temp state
        utmp(n) = uin(n)
      enddo
      do j = 1, m-1  ! utmp = utmp + h beta K_j
        do n = 1, this%nstate_  ! set utmp
          call saxpby_c(utmp(n), utmp(n), 1.0_GP, this%K_(j)(n), this%beta_(m,j)*dt)
        enddo
        tt = time + this%alpha_(m) * dt
        this%callback_( tt, utmp, uf, dt, K_(m) ); 
      enddo ! j-loop

    enddo ! stage m loop

    ! Combine stages to get step update:
    do n = 1, this%nstate_  ! uout = uin
      uout(n) = uin(n)
    enddo

    do m = 1, this%nstage_  

      do n = 1, this%nstate_ ! uout = uout + h * c_ * K:
        call saxpby_c(uout(n), uout(n), 1.0_GP, this%K_(m)(n), this%c_(m,j)*dt)
      enddo

    enddo ! m-loop

    CALL this%workspace_%free_complex_tmp(sum)

  end subroutine step_butcher


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


end module gexrk_mod
