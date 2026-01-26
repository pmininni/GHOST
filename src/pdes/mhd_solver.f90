! =====================================================================
! NAME       : mhd_solver.f90
! DESCRIPTION: Forms class for incompressible MHD solver, computing:
!
!              dv/dt +  w x v  = j x (b+B_0) - Grad p + nu Del^2 v
!              da/dt = (v - eps j) x (b+B_0) - Grad phi + eta Del^2 a 
!              ds_i/dt + v.Grad s_i = kappa_i Del^2 s_i 
!                                             i = 1, ..., numpassive
!              State ordering is:
!                v1, v2, v3, s1, s2, ..., s_numpassive
!
!              State sector ids are:
!                VELOCITY (VELOCITY+1, VELOCITY+2)
!                MAGNETIC (MAGNETIC+1, MAGNETIC+2)
!                PASSIVE  ( PASSIVE+1,  PASSIVE+2, ...)
!
! INPUT FILE : For solver='MHD', looks for a "&MHD" namelist with:
!              nu      : fluid kinematic viscosity
!              eta     : magnetic diffusivity
!              doB0    : do mean magnetic field, = .TRUE. or .FALSE.
!              dohall  : do Hall physics         = .TRUE. or .FALSE.
!              epsilon : amplitude of the Hall term
!              B0x     : amplitude of the guide field along x
!              B0y     : amplitude of the guide field along y
!              B0z     : amplitude of the guide field along z
!              npassive: number of passive scalars (default=0)
!
!              For npassive > 0, looks for a "&passive" namelist with:
!              kappa   : vector with npassive diffusivities
!
! DATE       : 01/17/26 (PDM)
! =====================================================================

module mhd_mod
  USE equationbase_mod
  USE gstate_mod

  IMPLICIT NONE

  ! ================= Solver traits ===================================
  type, public  :: NHTraits
    logical       :: doB0         = .FALSE. ! guide field flag
    logical       :: dohall       = .FALSE. ! compute hall term
    real(kind=GP) :: nu           = 0.0_GP  ! dissipation
    real(kind=GP) :: eta          = 0.0_GP  ! magnetic diffusivity
    real(kind=GP) :: epsilon      = 0.0_GP  ! Ion inertial length scale
    real(kind=GP), allocatable :: kappa(:)  ! scalar diffusivities
    real(kind=GP)              :: B0(3)     ! Guide field
  end type

  ! ================= Global parameters ===============================
  integer, parameter, public   :: MAXPASSIVE = 20 ! max # passive scalars

  ! ================= Solver ==========================================
  ! Define class:
  type, extends(MagneticBase) :: MHDSolver 
    ! Member data:
    logical           :: binit_=.false. ! is initialized?
    type  (NHTraits)  :: traits_

  CONTAINS
    procedure, public :: init          =>          init_impl ! init method
    procedure, public :: dudt          =>          dudt_impl ! RHS method
    procedure, public :: global        =>        global_impl ! Writes global qtys
    procedure, public :: spectra       =>       spectra_impl ! Writes spectra
    procedure, public :: state_size    =>    state_size_impl ! state size
    procedure, public :: sstate2istate => sstate2istate_impl ! state names
    procedure, public :: get_sstate    =>    get_sstate_impl ! get state name list
    procedure, public :: Solver_ctor   =>     MHDSolver_ctor ! constructor
    final             :: MHDSolver_dtor
  end type MHDSolver

CONTAINS

  ! ===================================================================
  ! Solver initialization, this is where parameter files are read
  ! ===================================================================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Subroutine to initialize the solver
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_impl(this)
    USE commtypes
    class  (MHDSolver), intent (inout) :: this

    ! Temporary data to read from namelists:
    logical                    :: doB0
    logical                    :: dohall
    integer                    :: npassive
    integer                    :: ierr
    real(kind=GP)              :: nu, eta, epsilon
    real(kind=GP)              :: B0x, B0y, B0z
    real(kind=GP), allocatable :: kappa(:)

    ! Required namelists:
    namelist/ MHD      / nu, eta, doB0, B0x, B0y, B0z, npassive
    namelist/ MHD      / dohall, epsilon, npassive
    namelist/ passive  / kappa

    call MPI_COMM_SIZE(MPI_COMM_WORLD,this%nprocs_,ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,this%myrank_,ierr)

    ! Get trait variables from input file:
    doB0     = .FALSE.
    dohall   = .FALSE.
    nu       = 0.0
    eta      = 0.0
    npassive = 0
    B0x = 0.0_GP; B0y = 0.0_GP; B0z = 0.0_GP
    if ( this%myrank_ .eq. 0 ) then
      open(1,file=this%infile_,status='unknown',form="formatted")
      read(1,NML=MHD)
      close(1)
    endif
    call mpi_bcast(nu       ,1 ,GC_REAL,    0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(eta      ,1 ,GC_REAL,    0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(epsilon  ,1 ,GC_REAL,    0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(doB0     ,1 ,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(dohall   ,1 ,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(B0x      ,1 ,GC_REAL,    0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(B0y      ,1 ,GC_REAL,    0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(B0z      ,1 ,GC_REAL,    0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(npassive ,1 ,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    this%numpassive_ = npassive
    if ( npassive .gt. 0 ) then
      allocate(kappa(npassive))
      if ( this%myrank_ .eq. 0 ) then
        if ( npassive .gt. MAXPASSIVE ) stop 'Max # of passive scalars exceeded'
        open(1,file=this%infile_,status='unknown',form="formatted")
        read(1,NML=passive)
        close(1)
      endif
      call mpi_bcast(kappa,npassive,GC_REAL,0,MPI_COMM_WORLD,ierr)
    endif

    ! Set traits from inputfile data:
    this%traits_%      doB0 = doB0
    this%traits_%    dohall = dohall
    this%traits_%        nu = nu
    this%traits_%       eta = eta
    this%traits_%   epsilon = epsilon
    this%traits_%        B0 = (/B0x,B0y,B0z/)
    if ( npassive .gt. 0 ) then
      if ( allocated(this%traits_%kappa) ) then
        deallocate(this%traits_%kappa);
      endif
      allocate(this%traits_%kappa(npassive))
      this%traits_%kappa = kappa
      deallocate(kappa)
    endif
   
    if ( dohall ) call this%workspace_%add_complex_entries(3)
    this%order_   = 2                        ! Time stepping order
    this%nd_      = 3                        ! 3d
    this%nc_      = this%nd_                 ! # field components
    this%VELOCITY = 1                        ! start of vel sector
    this%MAGNETIC = this%VELOCITY + this%nc_ ! start of mag sector
    this%PASSIVE  = this%MAGNETIC + this%nc_ ! start of scalar sector

    allocate(this%sstate_(this%state_size()))
    call this%get_sstate(this%sstate_)

    this%binit_ = .true.  
  end subroutine init_impl

  
  ! ===================================================================
  ! Computation of RHS, the solver equations are defined here
  ! ===================================================================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Function to compute RHS with guide field and hall terms
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine dudt_impl(this, time, uin, uf, dt, dudt) 
    use pseudospec_magnetic
    use ali
    use kes
    use var
    use grid
    use mpivars
!$  use threads
    implicit none

    class(MHDSolver), intent   (in)             :: this
    real   (kind=GP), intent   (in)             :: time, dt
    type    (GState), intent(inout), target     :: uin(:),uf(:)
    type    (GState), intent(inout)             :: dudt(:) 
    complex(kind=GP), pointer, dimension(:,:,:) :: fx,fy,fz,vx,vy,vz
    complex(kind=GP), pointer, dimension(:,:,:) :: mx,my,mz,ax,ay,az
    complex(kind=GP), pointer, dimension(:,:,:) :: C1,C2,C3,C4,C5,C6
    complex(kind=GP), pointer, dimension(:,:,:) :: C7,C8,C9,C10,C11
    complex(kind=GP), pointer, dimension(:,:,:) :: C12,C13,C14,C15
    real   (kind=GP)                            :: nu,eta,ep
    real   (kind=GP)                            :: b0x,b0y,b0z
    integer                                     :: i,j,k
    logical                                     :: bret

    nu  = this%traits_%nu
    eta = this%traits_%eta
    ep  = this%traits_%epsilon

    CALL this%workspace_%get_complex_tmp(C1,bret)
    CALL this%workspace_%get_complex_tmp(C2,bret)
    CALL this%workspace_%get_complex_tmp(C3,bret)
    CALL this%workspace_%get_complex_tmp(C4,bret)
    CALL this%workspace_%get_complex_tmp(C5,bret)
    CALL this%workspace_%get_complex_tmp(C6,bret)
    CALL this%workspace_%get_complex_tmp(C7,bret)
    CALL this%workspace_%get_complex_tmp(C8,bret)
    CALL this%workspace_%get_complex_tmp(C9,bret)
    CALL this%workspace_%get_complex_tmp(C10,bret)
    CALL this%workspace_%get_complex_tmp(C11,bret)
    CALL this%workspace_%get_complex_tmp(C12,bret)
    if ( this%traits_%dohall ) then
       CALL this%workspace_%get_complex_tmp(C13,bret)
       CALL this%workspace_%get_complex_tmp(C14,bret)
       CALL this%workspace_%get_complex_tmp(C15,bret)
    endif

    vx => uin(this%VELOCITY  )%ccomp
    vy => uin(this%VELOCITY+1)%ccomp
    vz => uin(this%VELOCITY+2)%ccomp
    fx => uf (this%VELOCITY  )%ccomp
    fy => uf (this%VELOCITY+1)%ccomp
    fz => uf (this%VELOCITY+2)%ccomp
    ax => uin(this%MAGNETIC  )%ccomp
    ay => uin(this%MAGNETIC+1)%ccomp
    az => uin(this%MAGNETIC+2)%ccomp
    mx => uf (this%MAGNETIC  )%ccomp
    my => uf (this%MAGNETIC+1)%ccomp
    mz => uf (this%MAGNETIC+2)%ccomp

    call rotor3(ay,az,C1,1)         ! b = curl(a)
    call rotor3(ax,az,C2,2)
    call rotor3(ax,ay,C3,3)
    if ( this%traits_%doB0 ) then   ! b = b + B_0
      if (myrank.eq.0) then
        b0x = this%traits_%B0(1)
        b0y = this%traits_%B0(2)
        b0z = this%traits_%B0(3)
        C1(1,1,1) = b0x*real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP)
        C2(1,1,1) = b0y*real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP)
        C3(1,1,1) = b0z*real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP)
      endif
    endif
    call prodre3(vx,vy,vz,C4,C5,C6) ! w x v
    call prodre3(C1,C2,C3,C7,C8,C9) ! j x b
    call nonlin3(C4,C5,C6,C7,C8,C9,C10,1) ! [- w x v + j x b - Grad p]_x
    call nonlin3(C4,C5,C6,C7,C8,C9,C11,2) ! [- w x v + j x b - Grad p]_y
    call nonlin3(C4,C5,C6,C7,C8,C9,C12,3) ! [- w x v + j x b - Grad p]_z
    call laplak3(ax,C4)             ! Del^2 ax
    call laplak3(ay,C5)             ! Del^2 ay
    call laplak3(az,C6)             ! Del^2 az
    if ( this%traits_%dohall ) then ! electron velocity: v_e = v - epsilon j
      do i = ista,iend               
        do j = 1,ny
          do k = 1,nz
            C7(k,j,i) = vx(k,j,i)+ep*C4(k,j,i)
            C8(k,j,i) = vy(k,j,i)+ep*C5(k,j,i)
            C9(k,j,i) = vz(k,j,i)+ep*C6(k,j,i)
          end do
        end do
      end do
      call vector3(C7,C8,C9,C1,C2,C3,C13,C14,C15) ! v_e x b (hall)
      call gauge3(C13,C14,C15,C1,1) ! [v_e x b - Grad phi]_x
      call gauge3(C13,C14,C15,C2,2) ! [v_e x b - Grad phi]_y
      call gauge3(C13,C14,C15,C3,3) ! [v_e x b - Grad phi]_z
    else
      call vector3(vx,vy,vz,C1,C2,C3,C7,C8,C9)    ! v x b (no hall)
      call gauge3(C7, C8, C9, C1,1) ! [v_e x b - Grad phi]_x
      call gauge3(C7, C8, C9, C2,2) ! [v_e x b - Grad phi]_y
      call gauge3(C7, C8, C9, C3,3) ! [v_e x b - Grad phi]_z
    endif
    call laplak3(vx,C7)             ! Del^2 vx
    call laplak3(vy,C8)             ! Del^2 vy
    call laplak3(vz,C9)             ! Del^2 vz

!$omp parallel do if (iend-ista.ge.nth) private (j,k)
    do i = ista,iend
!$omp parallel do if (iend-ista.ge.nth) private (k)
    do j = 1,ny
    do k = 1,nz
      if ((kn2(k,j,i).le.kmax).and.(kn2(k,j,i).ge.tiny)) then
        dudt(this%VELOCITY  )%ccomp(k,j,i) = nu*C7(k,j,i) + C10(k,j,i) + fx(k,j,i)
        dudt(this%VELOCITY+1)%ccomp(k,j,i) = nu*C8(k,j,i) + C11(k,j,i) + fy(k,j,i)
        dudt(this%VELOCITY+2)%ccomp(k,j,i) = nu*C9(k,j,i) + C12(k,j,i) + fz(k,j,i)
        dudt(this%MAGNETIC  )%ccomp(k,j,i) = eta*C4(k,j,i) + C1(k,j,i) + mx(k,j,i)
        dudt(this%MAGNETIC+1)%ccomp(k,j,i) = eta*C5(k,j,i) + C2(k,j,i) + my(k,j,i)
        dudt(this%MAGNETIC+2)%ccomp(k,j,i) = eta*C6(k,j,i) + C3(k,j,i) + mz(k,j,i)
      else
        dudt(this%VELOCITY  )%ccomp(k,j,i) = 0.0_GP
        dudt(this%VELOCITY+1)%ccomp(k,j,i) = 0.0_GP
        dudt(this%VELOCITY+2)%ccomp(k,j,i) = 0.0_GP
        dudt(this%MAGNETIC  )%ccomp(k,j,i) = 0.0_GP
        dudt(this%MAGNETIC+1)%ccomp(k,j,i) = 0.0_GP
        dudt(this%MAGNETIC+2)%ccomp(k,j,i) = 0.0_GP
      endif
    enddo
    enddo
    enddo

    CALL this%workspace_%free_complex_tmp(C1)
    CALL this%workspace_%free_complex_tmp(C2)
    CALL this%workspace_%free_complex_tmp(C3)
    CALL this%workspace_%free_complex_tmp(C4)
    CALL this%workspace_%free_complex_tmp(C5)
    CALL this%workspace_%free_complex_tmp(C6)
    CALL this%workspace_%free_complex_tmp(C7)
    CALL this%workspace_%free_complex_tmp(C8)
    CALL this%workspace_%free_complex_tmp(C9)
    CALL this%workspace_%free_complex_tmp(C10)
    CALL this%workspace_%free_complex_tmp(C11)
    CALL this%workspace_%free_complex_tmp(C12)
    if ( this%traits_%dohall ) then
       CALL this%workspace_%free_complex_tmp(C13)
       CALL this%workspace_%free_complex_tmp(C14)
       CALL this%workspace_%free_complex_tmp(C15)
    endif

    ! Compute passive scalars:
    call this%rhs_passive(uin, uf, this%traits_%kappa, dudt)

  end subroutine dudt_impl


  ! ===================================================================
  ! Computation of global quantities and spectra
  ! ===================================================================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Function to compute and write global quantities
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine global_impl(this, uin, uf, t) 
    use pseudospec_mhd
    use pseudospec_phd
    use status
    implicit none

    class(MHDSolver), intent(in)                :: this
    type    (GState), intent(in), target        :: uin(:), uf(:)
    integer         , intent(in)                :: t
    complex(kind=GP), pointer, dimension(:,:,:) :: fx,fy,fz,vx,vy,vz
    complex(kind=GP), pointer, dimension(:,:,:) :: mx,my,mz,ax,ay,az
    double precision                            :: eps,epm
    real   (kind=GP)                            :: rmp,rmq
    integer                                     :: i

    vx => uin(this%VELOCITY  )%ccomp
    vy => uin(this%VELOCITY+1)%ccomp
    vz => uin(this%VELOCITY+2)%ccomp
    fx => uf (this%VELOCITY  )%ccomp
    fy => uf (this%VELOCITY+1)%ccomp
    fz => uf (this%VELOCITY+2)%ccomp
    ax => uin(this%MAGNETIC  )%ccomp
    ay => uin(this%MAGNETIC+1)%ccomp
    az => uin(this%MAGNETIC+2)%ccomp
    mx => uf (this%MAGNETIC  )%ccomp
    my => uf (this%MAGNETIC+1)%ccomp
    mz => uf (this%MAGNETIC+2)%ccomp
    if ( .not. this%traits_%dohall ) then
      CALL mhdcheck(vx,vy,vz,ax,ay,az,t,dt,this%traits_%epsilon,1,1,0)
    else
      CALL mhdcheck(vx,vy,vz,ax,ay,az,t,dt,this%traits_%epsilon,1,2,0)
    endif
    CALL cross(vx,vy,vz,fx,fy,fz,eps,1)
    CALL cross(ax,ay,az,mx,my,mz,epm,0)
    CALL maxabs(vx,vy,vz,rmp,0)
    CALL maxabs(ax,ay,az,rmq,1)
    IF (myrank.eq.0) THEN
      OPEN(1,file='injection.txt',position='append')
      WRITE(1,FMT='(E13.6,E22.14,E22.14)') (t-1)*dt,eps,epm
      CLOSE(1)
      OPEN(1,file='maximum.txt',position='append')
      WRITE(1,FMT='(E13.6,E13.6,E13.6)') (t-1)*dt,rmp,rmq
      CLOSE(1)
    ENDIF
    do i = this%PASSIVE, this%PASSIVE+this%numpassive_-1
      call pscheck(uin(i)%ccomp,uf(i)%ccomp,t,dt,trim(this%sstate_(i)))
    end do
  end subroutine global_impl

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Function to compute and write spectra
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine spectra_impl(this, uin)
    use pseudospec_aniso
    use pseudospec_scalar
    use pseudospec_anisca
    use filefmt
    use status
    implicit none

    class(MHDSolver), intent(in)                :: this
    type    (GState), intent(in), target        :: uin(:)
    complex(kind=GP), pointer, dimension(:,:,:) :: vx,vy,vz
    complex(kind=GP), pointer, dimension(:,:,:) :: ax,ay,az
    complex(kind=GP), pointer, dimension(:,:,:) :: C1,C2,C3
    integer                                     :: i,j,k
    logical                                     :: bret

    WRITE(ext, fmtext) sind
    vx => uin(this%VELOCITY  )%ccomp
    vy => uin(this%VELOCITY+1)%ccomp
    vz => uin(this%VELOCITY+2)%ccomp
    ax => uin(this%MAGNETIC  )%ccomp
    ay => uin(this%MAGNETIC+1)%ccomp
    az => uin(this%MAGNETIC+2)%ccomp
    CALL spectrum(vx,vy,vz,ext,1,1)
    CALL spectrum(ax,ay,az,ext,0,1)
    if ( this%traits_%doB0 ) then
      CALL specpara(vx,vy,vz,ext,1,1)
      CALL specpara(ax,ay,az,ext,0,1)
      CALL specperp(vx,vy,vz,ext,1,1)
      CALL specperp(ax,ay,az,ext,0,1)
    endif
    if (  this%traits_%dohall ) then
      call this%workspace_%get_complex_tmp(C1,bret)
      call this%workspace_%get_complex_tmp(C2,bret)
      call this%workspace_%get_complex_tmp(C3,bret)
      DO i = ista,iend
        DO j = 1,ny
          DO k = 1,nz
            C1(k,j,i) = ax(k,j,i)+this%traits_%epsilon*vx(k,j,i)
            C2(k,j,i) = ay(k,j,i)+this%traits_%epsilon*vy(k,j,i)
            C3(k,j,i) = az(k,j,i)+this%traits_%epsilon*vz(k,j,i)
          END DO
        END DO
      END DO
      CALL spectrum(C1,C2,C3,ext,2,1)
      if ( this%traits_%doB0 ) then
        CALL specpara(C1,C2,C3,ext,2,1)
        CALL specperp(C1,C2,C3,ext,2,1)
      endif
      call this%workspace_%free_complex_tmp(C1)
      call this%workspace_%free_complex_tmp(C2)
      call this%workspace_%free_complex_tmp(C3)
    endif
    if ( this%numpassive_ .gt. 0) then
      do i = this%PASSIVE, this%PASSIVE+this%numpassive_-1
        call spectrsc(uin(i)%ccomp,ext,0,trim(this%sstate_(i)))
        if ( this%traits_%doB0 ) then
          call specscpa(uin(i)%ccomp,ext,0,trim(this%sstate_(i)))
          call specscpe(uin(i)%ccomp,ext,0,trim(this%sstate_(i)))
        endif
      end do
    endif

  end subroutine spectra_impl
  

  ! ===================================================================
  ! Solver specific methods
  ! ===================================================================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Constructor
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine MHDSolver_ctor(this, infile, workspace, plan)
    use iovar
    class (MHDSolver), intent(inout)         :: this
    type(GWorkspace) , intent(inout), target :: workspace
    type    (ioplan) , intent(inout), target :: plan
    character(len=*) , intent   (in)         :: infile
    this%infile_    =  infile    ! input file
    this%workspace_ => workspace
    this%planio_    => plan
    call this%init();
  end subroutine MHDSolver_ctor


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Destructor
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine MHDSolver_dtor(this) 
    type  (MHDSolver), intent(inout) :: this
    if (associated(this%workspace_))   nullify(this%workspace_)
    if (associated(this%planio_))      nullify(this%planio_)
    if (allocated(this%sstate_))       deallocate(this%sstate_)
    if (allocated(this%traits_%kappa)) deallocate(this%traits_%kappa)
  end subroutine MHDSolver_dtor


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Convert input state name to index in state vector
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine sstate2istate_impl(this, sstate, istate) 
    class(MHDSolver)              , intent   (in) :: this
    character(len=8)              , intent   (in) :: sstate(:)
    integer         , allocatable , intent(inout) :: istate(:)
    integer                                       :: i,j
    if ( size(sstate) .ne. size(istate) ) then
      stop 'MHDSolver::sstate2istate_impl: Incompatible sstate and istate'
    endif  
    do i = 1, size(sstate)
      istate(i) = -1 ! return unusable index
      do j = 1, size(this%sstate_)
        if ( sstate(i) .eq. this%sstate_(j) ) then
          istate(i) = j
        endif
      enddo
    enddo
  end subroutine sstate2istate_impl
  

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Get state variable names
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_sstate_impl(this, sstate) 
    class (MHDSolver)             , intent   (in) :: this
    character (len=8), allocatable, intent(inout) :: sstate(:)
    character(len=100)                            :: snum
    character(len=1)                              :: comp(3)
    integer                                       :: j
    comp = ['x', 'y', 'z']
    do j = this%VELOCITY,this%VELOCITY+this%nc_-1
       sstate(j) = 'v' // comp(j-this%VELOCITY+1)
    enddo
    do j = this%MAGNETIC,this%MAGNETIC+this%nc_-1
       sstate(j) = 'a' // comp(j-this%MAGNETIC+1)
    enddo
    do j = this%PASSIVE,this%PASSIVE+this%numpassive_-1
       write(snum,'(I0)') j-this%PASSIVE+1
       sstate(j) = 's' // trim(snum)
    enddo
  end subroutine get_sstate_impl
  

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Function to compute number of state members (equations)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PURE function state_size_impl(this) result(num)
    class(MHDSolver), intent(in) :: this
    integer                      :: num
    num = this%nc_               ! # vel. components
    num = num + this%nc_         ! # vec. potential components
    num = num + this%numpassive_ ! # scalars
  end function state_size_impl

end module mhd_mod
