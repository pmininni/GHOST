! =====================================================================
! NAME       : cmhd_solver.f90
! DESCRIPTION: Forms class for compressible MHD solver, computing:
!
!              dv/dt + w x v = - Grad(v^2/2 + h) + (j x B)/(amach^2 rho)
!                              + [nu Del^2 v + nu2 Grad(Div v)]/rho
!              da/dt = v_e x B - Grad phi + eta Del^2 a
!              drho/dt = - Div(rho v)
!              ds_i/dt + v.Grad s_i = kappa_i Del^2 s_i
!                                             i = 1, ..., numpassive
!
!              where B = curl(a) + B_0 is the total magnetic field,
!              j = curl(B) the current density, v_e = v - eps j/rho
!              the electron velocity when the Hall term is used
!              (v_e = v otherwise), and h the enthalpy of a polytropic
!              gas with p ~ rho^gamma, h = cp1 rho^(gamma-1)/2 with
!              cp1 = 2/[(gamma-1) smach^2], so that the sound speed at
!              rho = 1 is 1/smach. The Alfven speed at rho = 1 and
!              B = 1 is 1/amach. The mean modes (k = 0) of all fields
!              evolve: the mean density is part of the state. The
!              mean momentum is conserved, but not the mean velocity.
!
!              State ordering is:
!                v1, v2, v3, a1, a2, a3, rho, s1, s2, ..., s_numpassive
!
!              State sector ids are:
!                VELOCITY (VELOCITY+1, VELOCITY+2)     : velocity sector
!                MAGNETIC (MAGNETIC+1, MAGNETIC+2)     : magnetic sector
!                DENSITY                               : mass density
!                PASSIVE  ( PASSIVE+1,  PASSIVE+2, ...): passive scalars
!
!              The forcing state has a component for the density,
!              which is ignored (the continuity equation has no
!              source).
!
! INPUT FILE : For solver='CMHD', looks for a "&CMHD" namelist with:
!              fidir   : changes class binary input  dir (default: idir)
!              fodir   : changes class binary output dir (default: odir)
!              todir   : changes the class TXT output dir (default: '')
!              nu      : fluid kinematic viscosity
!              nu2     : second (bulk) viscosity for the Grad(Div v) term
!                        (see Stokeshyp, default=0)
!              Stokeshyp: = 0 uses nu2 + nu/3 as the bulk viscosity
!                              (default)
!                         > 0 Stokes hypothesis, nu2 = -2 nu/3
!                         < 0 uses nu2 as given
!              eta     : magnetic diffusivity
!              smach   : sound Mach number (sound speed 1/smach at rho=1)
!              amach   : Alfvenic Mach number (Alfven speed 1/amach for
!                        unit field and density)
!              gamma   : polytropic index of the equation of state
!                        (default=5/3)
!              doB0    : do mean magnetic field, = .TRUE. or .FALSE.
!              dohall  : do Hall physics         = .TRUE. or .FALSE.
!              epsilon : amplitude of the Hall term
!              B0x     : amplitude of the guide field along x
!              B0y     : amplitude of the guide field along y
!              B0z     : amplitude of the guide field along z
!              spectlod: spectral output level of detail (in [1,3]):
!                          1: All 1d spectra, KE fluxes
!                          2: 2D spectra, directional spectra,
!                             fluxes, helicity flux, if drot=true
!                          3: KE Fourier modes
!              npassive: number of passive scalars (default=0)
!
!              For npassive > 0, looks for a "&passive" namelist with:
!              kappa   : vector with npassive diffusivities
!
! DATE       : 09/11/26 (PDM)
! =====================================================================

module cmhd_mod
  USE equationbase_mod
  USE gstate_mod

  IMPLICIT NONE

  ! ================= Solver traits ===================================
  type, public  :: CMHDTraits
    logical       :: doB0         = .FALSE. ! guide field flag
    logical       :: dohall       = .FALSE. ! compute hall term
    integer       :: spectlod     = 1       ! standard level of spectra detail
    integer       :: stokeshyp    = 0       ! bulk viscosity option
    real(kind=GP) :: nu           = 0.0_GP  ! dissipation
    real(kind=GP) :: nu2          = 0.0_GP  ! bulk viscosity (as used)
    real(kind=GP) :: eta          = 0.0_GP  ! magnetic diffusivity
    real(kind=GP) :: epsilon      = 0.0_GP  ! Ion inertial length scale
    real(kind=GP) :: smach        = 1.0_GP  ! sound Mach number
    real(kind=GP) :: amach        = 1.0_GP  ! Alfvenic Mach number
    real(kind=GP) :: gam          = 5.0_GP/3.0_GP ! polytropic index
    real(kind=GP) :: gam1         = 2.0_GP/3.0_GP ! gamma - 1
    real(kind=GP) :: cp1          = 0.0_GP  ! 2/(gam1 smach^2)
    real(kind=GP) :: cp2          = 1.0_GP  ! 1/amach^2
    real(kind=GP), allocatable :: kappa(:)  ! scalar diffusivities
    real(kind=GP)              :: B0(3)     ! Guide field
  end type

  ! ================= Global parameters ===============================
  integer, parameter, public   :: MAXPASSIVE = 20 ! max # passive scalars

  ! ================= Solver ==========================================
  ! Define class:
  type, extends(CompMagneticBase) :: CMHDSolver
    ! Member data:
    logical           :: binit_ = .false. ! is initialized?
    type (CMHDTraits) :: traits_

  CONTAINS
    procedure, public :: init          =>          init_impl ! init method
    procedure, public :: dudt          =>          dudt_impl ! RHS method
    procedure, public :: global        =>        global_impl ! Writes global qtys
    procedure, public :: spectra       =>       spectra_impl ! Writes spectra
    procedure, public :: state_size    =>    state_size_impl ! state size
    procedure, public :: sstate2istate => sstate2istate_impl ! state names
    procedure, public :: get_sstate    =>    get_sstate_impl ! get state name list
    procedure, public :: Solver_ctor   =>    CMHDSolver_ctor ! constructor
    final             :: CMHDSolver_dtor
  end type CMHDSolver

CONTAINS

  ! ===================================================================
  ! Solver initialization, this is where parameter files are read
  ! ===================================================================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Subroutine to initialize the solver
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_impl(this)
    USE commtypes
    use status
    class(CMHDSolver), intent (inout) :: this

    ! Temporary data to read from namelists:
    logical                    :: doB0
    logical                    :: dohall
    integer                    :: npassive
    integer                    :: spectlod
    integer                    :: Stokeshyp
    integer                    :: ierr
    real(kind=GP)              :: nu, nu2, eta, epsilon
    real(kind=GP)              :: smach, amach, gamma
    real(kind=GP)              :: B0x, B0y, B0z
    real(kind=GP), allocatable :: kappa(:)
    character(len=128)         :: fidir, fodir, todir

    ! Required namelists:
    namelist/ CMHD    / fidir, fodir, todir
    namelist/ CMHD    / nu, nu2, Stokeshyp, eta, smach, amach, gamma
    namelist/ CMHD    / doB0, B0x, B0y, B0z
    namelist/ CMHD    / dohall, epsilon, npassive, spectlod
    namelist/ passive / kappa

    call MPI_COMM_SIZE(MPI_COMM_WORLD,this%nprocs_,ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,this%myrank_,ierr)

    ! Get I/O and trait variables from input file:
    fidir     = idir ! Set the default to status idir
    fodir     = odir ! Set the default to status odir
    todir     = '.'  ! Set the default to the current dir
    doB0      = .FALSE.
    dohall    = .FALSE.
    spectlod  = 1 ! standard lod
    Stokeshyp = 0
    nu        = 0.0_GP
    nu2       = 0.0_GP
    eta       = 0.0_GP
    epsilon   = 0.0_GP
    smach     = 1.0_GP
    amach     = 1.0_GP
    gamma     = 5.0_GP/3.0_GP
    npassive  = 0
    B0x = 0.0_GP; B0y = 0.0_GP; B0z = 0.0_GP
    if ( this%myrank_ .eq. 0 ) then
      open(1,file=this%infile_,status='unknown',form="formatted")
      read(1,NML=CMHD)
      close(1)
    endif
    call MPI_BCAST(fidir    ,128,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(fodir    ,128,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(todir    ,128,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(nu       ,1  ,GC_REAL      ,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(nu2      ,1  ,GC_REAL      ,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(Stokeshyp,1  ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(eta      ,1  ,GC_REAL      ,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(smach    ,1  ,GC_REAL      ,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(amach    ,1  ,GC_REAL      ,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(gamma    ,1  ,GC_REAL      ,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(epsilon  ,1  ,GC_REAL      ,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(doB0     ,1  ,MPI_LOGICAL  ,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(dohall   ,1  ,MPI_LOGICAL  ,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(B0x      ,1  ,GC_REAL      ,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(B0y      ,1  ,GC_REAL      ,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(B0z      ,1  ,GC_REAL      ,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(npassive ,1  ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(spectlod ,1  ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
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

    ! Set I/O and traits from inputfile data:
    this%idir_  = fidir ! If present in &CMHD, replaces the class default idir
    this%odir_  = fodir ! If present in &CMHD, replaces the class default odir
    this%todir_ = todir ! If present in &CMHD, replaces the class default todir
    this%traits_%     doB0 = doB0
    this%traits_%   dohall = dohall
    this%traits_% spectlod = spectlod
    this%traits_%stokeshyp = Stokeshyp
    this%traits_%       nu = nu
    this%traits_%      eta = eta
    this%traits_%  epsilon = epsilon
    this%traits_%    smach = smach
    this%traits_%    amach = amach
    this%traits_%      gam = gamma
    this%traits_%     gam1 = gamma - 1.0_GP
    this%traits_%      cp1 = 2.0_GP/((gamma - 1.0_GP)*smach*smach)
    this%traits_%      cp2 = 1.0_GP/(amach*amach)
    this%traits_%       B0 = (/B0x,B0y,B0z/)
    ! Bulk viscosity: the Grad(Div v) coefficient is nu2 + nu/3 with
    ! Stokeshyp = 0, -2 nu/3 (Stokes hypothesis) with Stokeshyp > 0,
    ! and nu2 as given with Stokeshyp < 0
    if ( Stokeshyp .gt. 0 ) then
      this%traits_%nu2 = -2.0_GP*nu/3.0_GP
    else if ( Stokeshyp .eq. 0 ) then
      this%traits_%nu2 = nu2 + nu/3.0_GP
    else
      this%traits_%nu2 = nu2
    endif
    if ( npassive .gt. 0 ) then
      if ( allocated(this%traits_%kappa) ) then
        deallocate(this%traits_%kappa);
      endif
      allocate(this%traits_%kappa(npassive))
      this%traits_%kappa = kappa
      deallocate(kappa)
    endif

    this%nd_      = 3                        ! 3d
    this%nc_      = this%nd_                 ! # field components
    this%VELOCITY = 1                        ! start of vel sector
    this%MAGNETIC = this%VELOCITY + this%nc_ ! start of mag sector
    this%DENSITY  = this%MAGNETIC + this%nc_ ! start of density sector
    this%PASSIVE  = this%DENSITY  + 1        ! start of scalar sector

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
    use pseudospec_compressible
    use ali
    use kes
    use var
    use grid
    use mpivars
    use gdevice, only: gdev_active
    implicit none

    class(CMHDSolver), intent(inout) :: this
    real   (kind=GP), intent   (in)             :: time, dt
    type(GStateComp), intent(inout), target     :: uin(:),uf(:)
    type(GStateComp), intent(inout), target     :: dudt(:)
    complex(kind=GP), pointer, dimension(:,:,:) :: fx,fy,fz,vx,vy,vz
    complex(kind=GP), pointer, dimension(:,:,:) :: mx,my,mz,ax,ay,az,rho
    complex(kind=GP), pointer, dimension(:,:,:) :: C1,C2,C3,C4,C5,C6
    complex(kind=GP), pointer, dimension(:,:,:) :: C7,C8,C9,C10,C11
    complex(kind=GP), pointer, dimension(:,:,:) :: C12,C13,C14,C15,C16
    complex(kind=GP), pointer, dimension(:,:,:) :: dvx,dvy,dvz,dax,day,daz
    complex(kind=GP), pointer, dimension(:,:,:) :: drho
    real   (kind=GP)                            :: nu,nu2,eta,ep,cp1,cp2,gam1
    real   (kind=GP)                            :: b0x,b0y,b0z
    integer                                     :: i,j,k
    logical                                     :: bret

    if ( .not. this%binit_ ) then
      stop 'CMHDSolver::dudt: Solver not initialized'
    endif

    nu   = this%traits_%nu
    nu2  = this%traits_%nu2
    eta  = this%traits_%eta
    ep   = this%traits_%epsilon
    cp1  = this%traits_%cp1
    cp2  = this%traits_%cp2
    gam1 = this%traits_%gam1

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
    CALL this%workspace_%get_complex_tmp(C13,bret)
    CALL this%workspace_%get_complex_tmp(C14,bret)
    CALL this%workspace_%get_complex_tmp(C15,bret)
    CALL this%workspace_%get_complex_tmp(C16,bret)

    vx  => uin(this%VELOCITY  )%ccomp
    vy  => uin(this%VELOCITY+1)%ccomp
    vz  => uin(this%VELOCITY+2)%ccomp
    fx  => uf (this%VELOCITY  )%ccomp
    fy  => uf (this%VELOCITY+1)%ccomp
    fz  => uf (this%VELOCITY+2)%ccomp
    ax  => uin(this%MAGNETIC  )%ccomp
    ay  => uin(this%MAGNETIC+1)%ccomp
    az  => uin(this%MAGNETIC+2)%ccomp
    mx  => uf (this%MAGNETIC  )%ccomp
    my  => uf (this%MAGNETIC+1)%ccomp
    mz  => uf (this%MAGNETIC+2)%ccomp
    rho => uin(this%DENSITY   )%ccomp
    dvx  => dudt(this%VELOCITY  )%ccomp
    dvy  => dudt(this%VELOCITY+1)%ccomp
    dvz  => dudt(this%VELOCITY+2)%ccomp
    dax  => dudt(this%MAGNETIC  )%ccomp
    day  => dudt(this%MAGNETIC+1)%ccomp
    daz  => dudt(this%MAGNETIC+2)%ccomp
    drho => dudt(this%DENSITY   )%ccomp

    call rotor3(ay,az,C1,1)         ! B = curl(a)
    call rotor3(ax,az,C2,2)
    call rotor3(ax,ay,C3,3)
    if ( this%traits_%doB0 ) then   ! B = B + B_0
      if (myrank.eq.0) then
        b0x = this%traits_%B0(1)*real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP)
        b0y = this%traits_%B0(2)*real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP)
        b0z = this%traits_%B0(3)*real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP)
        call setmode3(C1,1,1,1,cmplx(b0x,0.0_GP,kind=GP))
        call setmode3(C2,1,1,1,cmplx(b0y,0.0_GP,kind=GP))
        call setmode3(C3,1,1,1,cmplx(b0z,0.0_GP,kind=GP))
      endif
    endif
    call prodre3(C1,C2,C3,C4,C5,C6) ! j x B
    call divide(rho,C4,C5,C6)       ! (j x B)/rho
    call vector3(vx,vy,vz,C1,C2,C3,C7,C8,C9)    ! v x B
    if ( this%traits_%dohall ) then ! v_e x B = v x B - eps (j x B)/rho
      call saxpby_c(C7,C7,1.0_GP,C4,-ep)
      call saxpby_c(C8,C8,1.0_GP,C5,-ep)
      call saxpby_c(C9,C9,1.0_GP,C6,-ep)
    endif
    call gauge3(C7,C8,C9,C10,1)     ! [v_e x B - Grad phi]_x
    call gauge3(C7,C8,C9,C11,2)     ! [v_e x B - Grad phi]_y
    call gauge3(C7,C8,C9,C12,3)     ! [v_e x B - Grad phi]_z
    call prodre3(vx,vy,vz,C1,C2,C3) ! w x v (B is no longer needed)
    call saxpby_c(C1,C4,cp2,C1,-1.0_GP) ! (j x B)/(amach^2 rho) - w x v
    call saxpby_c(C2,C5,cp2,C2,-1.0_GP)
    call saxpby_c(C3,C6,cp2,C3,-1.0_GP)
    call laplak3(ax,C4)             ! Del^2 ax
    call laplak3(ay,C5)             ! Del^2 ay
    call laplak3(az,C6)             ! Del^2 az
    call gradpress(cp1,gam1,rho,vx,vy,vz,C7,C8,C9) ! Grad(v^2/2 + h)
    call divrhov(rho,vx,vy,vz,0,C16)   ! Div(rho v)
    call vdiss(nu,nu2,vx,vy,vz,C13,C14,C15) ! nu Del^2 v + nu2 Grad(Div v)
    call divide(rho,C13,C14,C15)    ! viscous term / rho

    ! The components of dudt are addressed through pointers (see the
    ! HD solver). All modes with kn2 <= kmax evolve, including k = 0:
    ! the mean density is part of the state, and the mean velocity
    ! is not conserved in a compressible flow (the mean momentum is)
#if defined(GHOST_GPU)
!$omp target teams distribute parallel do collapse(3) if(target: gdev_active)
#else
!$omp parallel do collapse(2) private (k)
#endif
    do i = ista,iend
      do j = 1,ny
        do k = 1,nz
          if (kn2(k,j,i).le.kmax) then
            dvx (k,j,i) = C1(k,j,i) - C7(k,j,i) + C13(k,j,i) + fx(k,j,i)
            dvy (k,j,i) = C2(k,j,i) - C8(k,j,i) + C14(k,j,i) + fy(k,j,i)
            dvz (k,j,i) = C3(k,j,i) - C9(k,j,i) + C15(k,j,i) + fz(k,j,i)
            dax (k,j,i) = eta*C4(k,j,i) + C10(k,j,i) + mx(k,j,i)
            day (k,j,i) = eta*C5(k,j,i) + C11(k,j,i) + my(k,j,i)
            daz (k,j,i) = eta*C6(k,j,i) + C12(k,j,i) + mz(k,j,i)
            drho(k,j,i) = -C16(k,j,i)
          else
            dvx (k,j,i) = 0.0_GP
            dvy (k,j,i) = 0.0_GP
            dvz (k,j,i) = 0.0_GP
            dax (k,j,i) = 0.0_GP
            day (k,j,i) = 0.0_GP
            daz (k,j,i) = 0.0_GP
            drho(k,j,i) = 0.0_GP
          endif
        end do
      end do
    end do

    CALL this%workspace_%free_complex_tmp(C16)
    CALL this%workspace_%free_complex_tmp(C15)
    CALL this%workspace_%free_complex_tmp(C14)
    CALL this%workspace_%free_complex_tmp(C13)
    CALL this%workspace_%free_complex_tmp(C12)
    CALL this%workspace_%free_complex_tmp(C11)
    CALL this%workspace_%free_complex_tmp(C10)
    CALL this%workspace_%free_complex_tmp(C9)
    CALL this%workspace_%free_complex_tmp(C8)
    CALL this%workspace_%free_complex_tmp(C7)
    CALL this%workspace_%free_complex_tmp(C6)
    CALL this%workspace_%free_complex_tmp(C5)
    CALL this%workspace_%free_complex_tmp(C4)
    CALL this%workspace_%free_complex_tmp(C3)
    CALL this%workspace_%free_complex_tmp(C2)
    CALL this%workspace_%free_complex_tmp(C1)

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
    use pseudospec_compr
    use pseudospec_phd
    use status
    implicit none

    class(CMHDSolver), intent(in)                :: this
    type (GStateComp), intent(in), target        :: uin(:), uf(:)
    integer          , intent(in)                :: t
    complex (kind=GP), pointer, dimension(:,:,:) :: fx,fy,fz,vx,vy,vz
    complex (kind=GP), pointer, dimension(:,:,:) :: mx,my,mz,ax,ay,az,rho
    double precision                             :: eps,epm
    real    (kind=GP)                            :: rmp,rmq
    integer                                      :: i

    vx  => uin(this%VELOCITY  )%ccomp
    vy  => uin(this%VELOCITY+1)%ccomp
    vz  => uin(this%VELOCITY+2)%ccomp
    fx  => uf (this%VELOCITY  )%ccomp
    fy  => uf (this%VELOCITY+1)%ccomp
    fz  => uf (this%VELOCITY+2)%ccomp
    ax  => uin(this%MAGNETIC  )%ccomp
    ay  => uin(this%MAGNETIC+1)%ccomp
    az  => uin(this%MAGNETIC+2)%ccomp
    mx  => uf (this%MAGNETIC  )%ccomp
    my  => uf (this%MAGNETIC+1)%ccomp
    mz  => uf (this%MAGNETIC+2)%ccomp
    rho => uin(this%DENSITY   )%ccomp
    ! Kinetic and internal energy with the mass density
    CALL energycompr(this%traits_%gam1,this%traits_%cp1,rho,vx,vy,vz,t,dt, &
                     this%todir_)
    ! Energies, helicities and divergences of v and b (<v^2> without rho)
    if ( .not. this%traits_%dohall ) then
      CALL mhdcheck(vx,vy,vz,ax,ay,az,t,dt,this%traits_%epsilon,1,1,1,this%todir_)
    else
      CALL mhdcheck(vx,vy,vz,ax,ay,az,t,dt,this%traits_%epsilon,1,2,1,this%todir_)
    endif
    CALL cross(vx,vy,vz,fx,fy,fz,eps,1)
    CALL cross(ax,ay,az,mx,my,mz,epm,0)
    CALL maxabs(vx,vy,vz,rmp,0)
    CALL maxabs(ax,ay,az,rmq,1)
    IF (myrank.eq.0) THEN
      OPEN(1,file=trim(this%todir_) // '/injection.txt',position='append')
      WRITE(1,FMT='(E13.6,E22.14,E22.14)') (t-1)*dt,eps,epm
      CLOSE(1)
      OPEN(1,file=trim(this%todir_) // '/maximum.txt',position='append')
      WRITE(1,FMT='(E13.6,E13.6,E13.6)'  ) (t-1)*dt,rmp,rmq
      CLOSE(1)
    ENDIF
    do i = this%PASSIVE, this%PASSIVE+this%numpassive_-1
      call pscheck(uin(i)%ccomp,uf(i)%ccomp,t,dt,this%todir_,trim(this%sstate_(i)))
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

    class(CMHDSolver), intent(in)                :: this
    type (GStateComp), intent(in), target        :: uin(:)
    complex (kind=GP), pointer, dimension(:,:,:) :: vx,vy,vz
    complex (kind=GP), pointer, dimension(:,:,:) :: ax,ay,az,rho
    complex (kind=GP), pointer, dimension(:,:,:) :: C1,C2,C3
    integer                                      :: i,j,k
    logical                                      :: bret

    WRITE(ext, fmtext) sind
    vx  => uin(this%VELOCITY  )%ccomp
    vy  => uin(this%VELOCITY+1)%ccomp
    vz  => uin(this%VELOCITY+2)%ccomp
    ax  => uin(this%MAGNETIC  )%ccomp
    ay  => uin(this%MAGNETIC+1)%ccomp
    az  => uin(this%MAGNETIC+2)%ccomp
    rho => uin(this%DENSITY   )%ccomp
    CALL spectrum(vx,vy,vz,this%todir_,ext,1,1)
    CALL spectrum(ax,ay,az,this%todir_,ext,0,1)
    CALL spectrsc(rho,this%todir_,ext,-1)          ! rhospect.XXX.txt
    if ( this%traits_%doB0 ) then
      CALL specpara(vx,vy,vz,this%todir_,ext,1,1)
      CALL specpara(ax,ay,az,this%todir_,ext,0,1)
      CALL specperp(vx,vy,vz,this%todir_,ext,1,1)
      CALL specperp(ax,ay,az,this%todir_,ext,0,1)
      CALL specscpa(rho,this%todir_,ext,-1)
      CALL specscpe(rho,this%todir_,ext,-1)
    endif
    if (  this%traits_%dohall ) then ! generalized helicity spectrum
      call this%workspace_%get_complex_tmp(C1,bret)
      call this%workspace_%get_complex_tmp(C2,bret)
      call this%workspace_%get_complex_tmp(C3,bret)
      !$omp parallel do collapse(2) private (k)
      DO i = ista,iend
        DO j = 1,ny
          DO k = 1,nz
            C1(k,j,i) = ax(k,j,i)+this%traits_%epsilon*vx(k,j,i)
            C2(k,j,i) = ay(k,j,i)+this%traits_%epsilon*vy(k,j,i)
            C3(k,j,i) = az(k,j,i)+this%traits_%epsilon*vz(k,j,i)
          END DO
        END DO
      END DO
      CALL spectrum(C1,C2,C3,this%todir_,ext,2,1)
      if ( this%traits_%doB0 ) then
        CALL specpara(C1,C2,C3,this%todir_,ext,2,1)
        CALL specperp(C1,C2,C3,this%todir_,ext,2,1)
      endif
      call this%workspace_%free_complex_tmp(C1)
      call this%workspace_%free_complex_tmp(C2)
      call this%workspace_%free_complex_tmp(C3)
    endif
    if ( this%numpassive_ .gt. 0) then
      do i = this%PASSIVE, this%PASSIVE+this%numpassive_-1
        call spectrsc(uin(i)%ccomp,this%todir_,ext,0,trim(this%sstate_(i)))
        if ( this%traits_%doB0 ) then
          call specscpa(uin(i)%ccomp,this%todir_,ext,0,trim(this%sstate_(i)))
          call specscpe(uin(i)%ccomp,this%todir_,ext,0,trim(this%sstate_(i)))
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
  subroutine CMHDSolver_ctor(this, infile, workspace, plan)
    use iovar
    class(CMHDSolver), intent(inout)         :: this
    type(GWorkspace) , intent(inout), target :: workspace
    type    (ioplan) , intent(inout), target :: plan
    character(len=*) , intent   (in)         :: infile
    this%infile_    =  infile    ! input file
    this%workspace_ => workspace
    this%planio_    => plan
    call this%init();
  end subroutine CMHDSolver_ctor


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Destructor
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine CMHDSolver_dtor(this)
    type  (CMHDSolver), intent(inout) :: this
    if (associated(this%workspace_))   nullify(this%workspace_)
    if (associated(this%planio_))      nullify(this%planio_)
    if (allocated(this%sstate_))       deallocate(this%sstate_)
    if (allocated(this%traits_%kappa)) deallocate(this%traits_%kappa)
  end subroutine CMHDSolver_dtor


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Convert input state name to index in state vector
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine sstate2istate_impl(this, sstate, istate)
    class(CMHDSolver), intent   (in) :: this
    character (len=8), intent   (in) :: sstate(:)
    integer          , intent(inout) :: istate(:)
    integer                          :: i,j
    if ( size(sstate) .ne. size(istate) ) then
      stop 'CMHDSolver::sstate2istate_impl: Incompatible sstate and istate'
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
    class (CMHDSolver), intent   (in) :: this
    character  (len=8), intent(inout) :: sstate(:)
    character(len=100)                :: snum
    character  (len=1)                :: comp(3)
    integer                           :: j
    comp = ['x', 'y', 'z']
    do j = this%VELOCITY,this%VELOCITY+this%nc_-1
       sstate(j) = 'v' // comp(j-this%VELOCITY+1)
    enddo
    do j = this%MAGNETIC,this%MAGNETIC+this%nc_-1
       sstate(j) = 'a' // comp(j-this%MAGNETIC+1)
    enddo
    sstate(this%DENSITY) = 'rho'
    do j = this%PASSIVE,this%PASSIVE+this%numpassive_-1
       write(snum,'(I0)') j-this%PASSIVE+1
       sstate(j) = 's' // trim(snum)
    enddo
  end subroutine get_sstate_impl


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Function to compute number of state members (equations)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PURE function state_size_impl(this) result(num)
    class(CMHDSolver), intent(in) :: this
    integer                       :: num
    num = this%nc_               ! # vel. components
    num = num + this%nc_         ! # vec. potential components
    num = num + 1                ! # mass density
    num = num + this%numpassive_ ! # scalars
  end function state_size_impl

end module cmhd_mod
