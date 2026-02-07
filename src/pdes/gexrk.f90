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
  USE class_GWorkspace3D
  USE gstate_mod


  IMPLICIT NONE

  ! ================= Stepper traits ===================================
  type, public  :: GEXRKTraits
!   logical       :: bSSP       = .FALSE. ! strong stability flag?
!   logical       :: blowstore  = .FALSE. ! low storage version, if possible?
    integer       :: norder     = 2  
    integer       :: nstage     = 2 
  end type


  ! ================= Stepper ==========================================
  ! Define class:
  type :: GExRKStepper
    ! Member data:
    type(GWorkspace), pointer     :: workspace_
    integer                       :: nstate_   ! no. state members
    integer                       :: myrank_   ! MPI rank
    integer                       :: nprocs_   ! MPI rank 
    logical                       :: binit_=.false. ! is initialized?
    type  (GEXRKTraits)           :: traits_
    real   (kind=GP),        , allocatable, &
                               dimension    (:) :: alpha_, c_
    real   (kind=GP),        , allocatable, &
                               dimension  (:,:) :: beta_
    complex(kind=GP), pointer, dimension(:,:,:) :: K_

  CONTAINS
    procedure, public :: set_order     ! set order, nstages
    procedure, public :: init          ! initialize
    procedure, public :: step          ! take one timestep
    procedure, public :: GEXRKStepper_ctor 
    final             :: GExRKStepper_dtor

    procedure, private:: init_butcher  ! initialize Butcher table
    procedure, private:: set_tmp       ! set data from tmp pool
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

    this%infile_    =  infile    ! input file
    this%workspace_ => workspace

    this%nstate_    = nstate ! must be set before init is called!

    call this%init()
  end subroutine GExRKStepper_ctor


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Destructor
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine GExRKStepper_dtor(this) 
    type  (GExRKStepper), intent(inout) :: this

!   if (associated(this%workspace_))   nullify(this%workspace_)
    if (allocated(this%alpha_)) deallocate(this%alpha_)
    if (allocated (this%beta_)) deallocate (this%beta_)
    if (allocated    (this%c_)) deallocate    (this%c_)

    
  end subroutine GExRKStepper_dtor


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Subroutine to initialize the stepper
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init(this)
!   USE commtypes
    class  (GExRKStepper), intent (inout) :: this

    ! Temporary data to read from namelists:
    integer                    :: norder,nstage
    integer                    :: ierr
    
    ! Required namelists:
    namelist/ stepper    / norder, nstage

    ! Get stepper traits from input file:
    norder   = 2
    nstage   = 2
!   bSSP     = .false.
    if ( this%myrank_ .eq. 0 ) then
      open(1,file=this%infile_,status='unknown',form="formatted")
      read(1,NML=HD)
      close(1)
    endif
    call mpi_bcast(norder   ,1 ,MPI_INTEGER,    0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(nstage   ,1 ,MPI_INTEGER,    0,MPI_COMM_WORLD,ierr)
!   call mpi_bcast(bSSP     ,1 ,MPI_LOGICAL,    0,MPI_COMM_WORLD,ierr)

    if ( norder .ne. nstage ) then
      stop 'GExRKStepper::init: For now, norder must equal nstage'
    endif

    this%traits_%norder = norder
    this%traits_%nstage = nstage
!   this%traits_%bSSP   = bSSP

    call MPI_COMM_SIZE(MPI_COMM_WORLD,this%nprocs_,ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,this%myrank_,ierr)

    call this%init_butcher()

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
    if      ( this%traits_%norder .eq. 1 ) then
    else if ( this%traits_%norder .eq. 1 ) then
    endif
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
    end select 

  end subroutine init_butcher

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Subroutine to allow volatile set of tmp
  !! data from workspace
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine set_tmp(this)
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
    if      ( this%traits_%norder .eq. 1 ) then
    else if ( this%traits_%norder .eq. 1 ) then
    endif
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
    end select 

  end subroutine init_butcher


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Function to take one step
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine step(this, time, uin, uf, dt, dudt) 
    use ali
    use kes
    use var
    use grid
    use mpivars
!$  use threads
    implicit none

    class (GExRKStepper), intent   (in)             :: this
    real   (kind=GP), intent   (in)             :: time, dt
    type    (GState), intent(inout), target     :: uin(:),uf(:)
    type    (GState), intent(inout)             :: dudt(:) 
    complex(kind=GP), pointer, dimension(:,:,:) :: fx,fy,fz,vx,vy,vz
    complex(kind=GP), pointer, dimension(:,:,:) :: fth,th
    complex(kind=GP), pointer, dimension(:,:,:) :: C1,C2,C3,C4,C5,C6
    complex(kind=GP), pointer, dimension(:,:,:) :: C7,C8
    integer                                     :: i,j,k
    logical                                     :: bret
       
    if ( size(uin) .ne. this.nstate_ ) then
      stop 'GExRKStepperi::step: Inconsistent input state'
    endif

    CALL this%workspace_%get_complex_tmp(C1,bret)
    CALL this%workspace_%get_complex_tmp(C2,bret)
    CALL this%workspace_%get_complex_tmp(C3,bret)
    CALL this%workspace_%get_complex_tmp(C4,bret)
    CALL this%workspace_%get_complex_tmp(C5,bret)
    CALL this%workspace_%get_complex_tmp(C6,bret)
    CALL this%workspace_%get_complex_tmp(C7,bret)
    CALL this%workspace_%get_complex_tmp(C8,bret)

    vx  => uin(this%VELOCITY  )%ccomp
    vy  => uin(this%VELOCITY+1)%ccomp
    vz  => uin(this%VELOCITY+2)%ccomp
    th  => uin(this%TEMP)      %ccomp
    fx  => uf (this%VELOCITY  )%ccomp
    fy  => uf (this%VELOCITY+1)%ccomp
    fz  => uf (this%VELOCITY+2)%ccomp
    fth => uf (this%TEMP)      %ccomp
      
    CALL this%workspace_%free_complex_tmp(C1)
    CALL this%workspace_%free_complex_tmp(C2)
    CALL this%workspace_%free_complex_tmp(C3)
    CALL this%workspace_%free_complex_tmp(C4)
    CALL this%workspace_%free_complex_tmp(C5)
    CALL this%workspace_%free_complex_tmp(C6)
    CALL this%workspace_%free_complex_tmp(C7)
    CALL this%workspace_%free_complex_tmp(C8)

  end subroutine step


end module gexrk_mod
