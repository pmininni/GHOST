! =====================================================================
! NAME       : moist_solver.f90
! DESCRIPTION: Forms class for moist Boussinesq solvers, computing:
!
!              dv/dt + (w + 2 Omega) x v = - Grad p + 
!                          + xmom (Nu bu Hu + Ns bs Hs) \hat{z} + nu Del^2 v
!              dbu/dt + v.Grad bu  = - xtemp Nu v_z + bkappa Del^2 bu
!              dbs/dt + v.Grad bs  = - xtemp Ns v_z + bkappa Del^2 bs
!
!              ds_i/dt + v.Grad s_i = kappa_i Del^2 s_i 
!                                             i = 1, ..., numpassive
!              Div v = 0,
!               
!              where Hu is the heaviside fuction, Hs = 1-Hu, Ns=bvsat*bvfreq
!              and Nu = bvuns*bvfreq are buoyancy frequencies, nu is the
!              viscosity, bkappa is a diffusivity, and kappa is a vector with
!              the diffusivities of the passive scalar. bu and bs have units
!              of velocity, and the buoyancy with units of acceleration is
!              obtained by multiplying them by their buoyancy frequencies.
!
!              State ordering is:
!                v1, v2, v3, th1, th2, s1, s2, ..., s_numpassive
!
!              State sector ids are:
!                VELOCITY (VELOCITY+1, VELOCITY+2)   : momentum sector
!                ACTIVESC (ACTIVESC+1)               : active scalar sector
!                PASSIVE  (PASSIVE+1, PASSIVE+2, ...): passive scalar sector
!
! INPUT FILE : For solver='MOIST', looks for a "&MOIST" namelist with:
!              fidir   : changes class binary input  dir (default: idir)
!              fodir   : changes class binary output dir (default: odir)
!              todir   : changes the class TXT output dir (default: '')
!              nu      : fluid kinematic viscosity
!              bkappa  : active scalar diffusivity
!              bvuns   : Unsaturated multiplicative coefficient
!              bvsat   : Saturated multiplicative coefficient
!              bvfreq  : Base Brunt-Vaisala frequency
!              xmom    : factor multiplying buoyancy force in momentum eq.
!              xtemp   : factor multiplying vertical veloc in temp. eq.
!              dorot   : do rotation = .TRUE. or .FALSE.
!              omegax  : amplitude of the uniform rotation along x
!              omegay  : amplitude of the uniform rotation along y
!              omegaz  : amplitude of the uniform rotation along z
!              npassive: number of passive scalars (default=0)
!              spectlod: spectral output level of detail (in [1,4]):
!                          1: All 1d spectra, KE and PE fluxes
!                          2: 2D spectra, helicity flluxes
!                          3: KE Fourier modes
!                          4: PV spectra, horizontally-avged quantities
!
!              For npassive > 0, looks for a "&passive" namelist with:
!              kappa   : vector with npassive diffusivities
!
! DATE       : 3/30/26 (JBG)
! =====================================================================

module moist_mod
  USE equationbase_mod
  USE gstate_mod

  implicit none

  ! ================= Solver traits ===================================
  type, public    :: MOISTTraits
    logical       :: dorot         = .FALSE. ! rotation flag
    integer       :: spectlod      = 1       ! standard level of spectra detail 
    real(kind=GP) :: nu            = 0.0_GP  ! dissipation
    real(kind=GP) :: bkappa        = 0.0_GP  ! diffusivity
    real(kind=GP) :: bvuns         = 0.0_GP  ! unsaturated coef
    real(kind=GP) :: bvsat         = 0.0_GP  ! saturated coef
    real(kind=GP) :: xmom          = 1.0_GP  ! buoy mom coef    
    real(kind=GP) :: xtemp         = 1.0_GP  ! buoy temp coef
    real(kind=GP) :: bvfreq        = 0.0_GP  ! Brunt-Vaisala freq
    real(kind=GP), allocatable :: kappa(:)   ! diffusivities
    real(kind=GP)              :: omega(3)   ! rotation vector
  end type

  ! ================= Global parameters ===============================
  integer, parameter, public   :: MAXPASSIVE = 20 ! max # passive scalars

  ! ================= Solver ==========================================
  ! Define class:
  type, extends(ActiveScalarBase) :: MOISTSolver 
    ! Member data:
    logical              :: binit_  = .false. ! is initialized?
    type  (MOISTTraits)  :: traits_
  contains
    procedure, public :: init          =>          init_impl ! init method
    procedure, public :: dudt          =>          dudt_impl ! RHS method
    procedure, public :: global        =>        global_impl ! Writes global qtys
    procedure, public :: spectra       =>       spectra_impl ! Writes spectra
    procedure, public :: state_size    =>    state_size_impl ! state size
    procedure, public :: sstate2istate => sstate2istate_impl ! state names
    procedure, public :: get_sstate    =>    get_sstate_impl ! get state name list
    procedure, public :: Solver_ctor   =>   MOISTSolver_ctor ! constructor
    final             :: MOISTSolver_dtor
  end type MOISTSolver

contains

  ! ===================================================================
  ! Solver initialization, this is where parameter files are read
  ! ===================================================================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Subroutine to initialize the solver
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_impl(this)
    USE commtypes
    use status
    class  (MOISTSolver), intent (inout) :: this

    ! Temporary data to read from namelists:
    logical                    :: dorot
    integer                    :: npassive
    integer                    :: spectlod
    integer                    :: ierr
    real(kind=GP)              :: nu, bkappa, bvuns, bvsat, bvfreq
    real(kind=GP)              :: xmom, xtemp
    real(kind=GP)              :: omegax, omegay, omegaz
    real(kind=GP), allocatable :: kappa(:)
    character(len=128)         :: fidir, fodir, todir

    ! Required namelists:
    namelist/ MOIST   / fidir, fodir, todir
    namelist/ MOIST   / nu, bkappa, bvuns, bvsat, bvfreq, xmom, xtemp
    namelist/ MOIST   / dorot, omegax, omegay, omegaz, npassive, spectlod
    namelist/ passive / kappa

    call MPI_COMM_SIZE(MPI_COMM_WORLD,this%nprocs_,ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,this%myrank_,ierr)

    ! Get I/O and trait variables from input file:
    fidir    = idir ! Set the default to status idir
    fodir    = odir ! Set the default to status odir
    todir    = '.'  ! Set the default to the current dir
    dorot    = .FALSE.
    spectlod = 1 ! standard lod
    nu       = 0.0_GP
    bkappa   = 0.0_GP
    bvuns    = 0.0_GP
    bvsat    = 0.0_GP
    bvfreq   = 0.0_GP
    xmom     = 1.0_GP
    xtemp    = 1.0_GP
    npassive = 0
    omegax   = 0.0_GP; omegay = 0.0_GP; omegaz = 0.0_GP
    if ( this%myrank_ .eq. 0 ) then
      open(1,file=this%infile_,status='unknown',form="formatted")
      read(1,NML=MOIST)
      close(1)
    endif
    call MPI_BCAST(fidir    ,128,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(fodir    ,128,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(todir    ,128,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(nu       ,1  ,GC_REAL      ,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(bkappa   ,1  ,GC_REAL      ,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(bvuns    ,1  ,GC_REAL      ,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(bvsat    ,1  ,GC_REAL      ,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(bvfreq   ,1  ,GC_REAL      ,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(xmom     ,1  ,GC_REAL      ,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(xtemp    ,1  ,GC_REAL      ,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(dorot    ,1  ,MPI_LOGICAL  ,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(omegax   ,1  ,GC_REAL      ,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(omegay   ,1  ,GC_REAL      ,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(omegaz   ,1  ,GC_REAL      ,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(npassive ,1  ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(spectlod ,1  ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)

    this%numpassive_ = npassive
    this%numactivesc_  = 2
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
    this%idir_  = fidir ! If present in &HD, replaces the class default idir
    this%odir_  = fodir ! If present in &HD, replaces the class default odir
    this%todir_ = todir ! If present in &HD, replaces the class default todir
    this%traits_%   dorot = dorot
    this%traits_%spectlod = spectlod
    this%traits_%      nu = nu
    this%traits_%  bkappa = bkappa
    this%traits_%   bvuns = bvuns
    this%traits_%   bvsat = bvsat
    this%traits_%  bvfreq = bvfreq
    this%traits_%    xmom = xmom
    this%traits_%   xtemp = xtemp   
    this%traits_%   omega = (/omegax,omegay,omegaz/)
    if ( npassive .gt. 0 ) then
      if ( allocated(this%traits_%kappa) ) then
        deallocate(this%traits_%kappa);
      endif
      allocate(this%traits_%kappa(npassive))
      this%traits_%kappa = kappa
      deallocate(kappa)
    endif

    this%nd_        = 3                                 ! 3d
    this%nc_        = this%nd_                          ! #vel field components
    this%VELOCITY   = 1                                 ! start of vel sector
    this%ACTIVESC   = this%VELOCITY + this%nc_          ! Start of active scalar sector
    this%PASSIVE    = this%ACTIVESC + this%numactivesc_ ! Start of passive scalar sector

    allocate(this%sstate_(this%state_size()))
    call this%get_sstate(this%sstate_)
    this%binit_ = .true.  
  end subroutine init_impl

  
  ! ===================================================================
  ! Computation of RHS, the solver equations are defined here
  ! ===================================================================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Function to compute RHS with rotation
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine dudt_impl(this, time, uin, uf, dt, dudt) 
    use pseudospec_scalar
    use ali
    use kes
    use var
    use grid
    use mpivars
    use commtypes
    use fft
!$  use threads
    implicit none

    class (MOISTSolver), intent   (in)          :: this
    real   (kind=GP), intent   (in)             :: time, dt
    type(GStateComp), intent(inout), target     :: uin(:),uf(:)
    type(GStateComp), intent(inout)             :: dudt(:)
    complex(kind=GP), pointer, dimension(:,:,:) :: fx,fy,fz,vx,vy,vz
    complex(kind=GP), pointer, dimension(:,:,:) :: fth1, fth2,th1,th2
    complex(kind=GP), pointer, dimension(:,:,:) :: C1,C2,C3,C4,C5,C6
    complex(kind=GP), pointer, dimension(:,:,:) :: C7,C8
    real   (kind=GP), pointer, dimension(:,:,:) :: R1,R2
    real   (kind=GP)                            :: bkappa,bvuns, bvsat,bvfreq,nu
    real   (kind=GP)                            :: xmom,xtemp
    real   (kind=GP)                            :: tmp
    real   (kind=GP)                            :: omegax,omegay,omegaz
    integer                                     :: i,j,k
    logical                                     :: bret

    if ( .not. this%binit_ ) then
      stop 'MOISTSolver::dudt: Solver not initialized'
    endif
       
    nu     = this%traits_%nu
    bkappa = this%traits_%bkappa
    bvfreq = this%traits_%bvfreq
    bvuns  = this%traits_%bvuns
    bvsat  = this%traits_%bvsat
    xmom   = this%traits_%xmom  * this%traits_%bvfreq
    xtemp  = this%traits_%xtemp * this%traits_%bvfreq

    call this%workspace_%get_complex_tmp(C1,bret)
    call this%workspace_%get_complex_tmp(C2,bret)
    call this%workspace_%get_complex_tmp(C3,bret)
    call this%workspace_%get_complex_tmp(C4,bret)
    call this%workspace_%get_complex_tmp(C5,bret)
    call this%workspace_%get_complex_tmp(C6,bret)
    call this%workspace_%get_complex_tmp(C7,bret)
    call this%workspace_%get_complex_tmp(C8,bret)
    call this%workspace_%get_real_tmp(R1,bret)
    call this%workspace_%get_real_tmp(R2,bret)

    vx   => uin(this%VELOCITY  )%ccomp
    vy   => uin(this%VELOCITY+1)%ccomp
    vz   => uin(this%VELOCITY+2)%ccomp
    th1  => uin(this%ACTIVESC  )%ccomp
    th2  => uin(this%ACTIVESC+1)%ccomp
    fx   => uf (this%VELOCITY  )%ccomp
    fy   => uf (this%VELOCITY+1)%ccomp
    fz   => uf (this%VELOCITY+2)%ccomp
    fth1 => uf (this%ACTIVESC  )%ccomp
    fth2 => uf (this%ACTIVESC+1)%ccomp

    call prodre3(vx,vy,vz,C4,C5,C6)                    ! w x v
    if ( this%traits_%dorot ) then
      omegax = this%traits_%omega(1)
      omegay = this%traits_%omega(2)
      omegaz = this%traits_%omega(3)
      call saxpby_c(C1, vz, 2*omegay, vy, -2.0*omegaz) ! 2 Omega x v
      call saxpby_c(C2, vx, 2*omegaz, vz, -2.0*omegax)
      call saxpby_c(C3, vy, 2*omegax, vx, -2.0*omegay)
!$omp parallel do collapse(2) private (k)
      do i = ista,iend
      do j = 1,ny
      do concurrent (k=1:nz)
         C4(k,j,i) = C4(k,j,i) + C1(k,j,i) ! (w x v + 2 Omega x v)_x
         C5(k,j,i) = C5(k,j,i) + C2(k,j,i) ! (w x v + 2 Omega x v)_y
         C6(k,j,i) = C6(k,j,i) + C3(k,j,i) ! (w x v + 2 Omega x v)_z
      end do
      end do
      end do
    endif

    C7 = th1
    call fftp3d_complex_to_real(plancr, C7, R1, MPI_COMM_WORLD)
    C8 = th2
    call fftp3d_complex_to_real(plancr, C8, R2, MPI_COMM_WORLD)
    tmp = 1.0_GP/(real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))
!$omp parallel do collapse(2) private (i)
    do k = ksta,kend
    do j = 1,ny
    do concurrent (i=1:nx)   ! Buoyancy force w/Heaviside
      if ( (bvuns*R1(i,j,k)).gt.(bvsat*R2(i,j,k)) ) then
        R1(i,j,k) = tmp*bvuns*xmom*R1(i,j,k)
      else
        R1(i,j,k) = tmp*bvsat*xmom*R2(i,j,k)
      endif
    end do
    end do
    end do
  
    call fftp3d_real_to_complex(planrc, R1, C7, MPI_COMM_WORLD)
!$omp parallel do collapse(2) private (k)
    do i = ista,iend
    do j = 1,ny                            ! NL term in z + Buoyancy
    do concurrent (k=1:nz)                 ! It becomes negative as it changes
      C6(k,j,i) = C6(k,j,i) + C7(k,j,i)    ! sign after the call to nonlhd
    end do
    end do
    end do

    call nonlhd3(C4,C5,C6,C1,1)   ! -[(w + 2 Omega) x v + Grad p]_x
    call nonlhd3(C4,C5,C6,C2,2)   ! -[(w + 2 Omega) x v + Grad p]_y
    call nonlhd3(C4,C5,C6,C3,3)   ! -[(w + 2 Omega) x v + Grad p]_z
    call advect3(vx,vy,vz,th1,C7) ! -(v.Grad) th1
    call advect3(vx,vy,vz,th2,C8) ! -(v.Grad) th2

!$omp parallel do collapse(2) private (k)
    do i = ista,iend
    do j = 1,ny
    do concurrent (k=1:nz)   ! heat 'currrents'
      C7(k,j,i) = C7(k,j,i) + bvuns*xtemp*vz(k,j,i)
      C8(k,j,i) = C8(k,j,i) + bvsat*xtemp*vz(k,j,i)
    end do
    end do
    end do

    call laplak3(vx,C4)           ! Del^2 vx
    call laplak3(vy,C5)           ! Del^2 vy
    call laplak3(vz,C6)           ! Del^2 vz
    call laplak3(th1,th1)         ! Del^2 th1
    call laplak3(th2,th2)         ! Del^2 th2

!$omp parallel do collapse(2) private (k)
    do i = ista,iend
    do j = 1,ny
    do concurrent (k=1:nz)
      if ((kn2(k,j,i).le.kmax).and.(kn2(k,j,i).ge.tiny)) then
        dudt(this%VELOCITY  )%ccomp(k,j,i) = nu*C4(k,j,i) + C1(k,j,i) + fx(k,j,i)
        dudt(this%VELOCITY+1)%ccomp(k,j,i) = nu*C5(k,j,i) + C2(k,j,i) + fy(k,j,i)
        dudt(this%VELOCITY+2)%ccomp(k,j,i) = nu*C6(k,j,i) + C3(k,j,i) + fz(k,j,i)
        dudt(this%ACTIVESC  )%ccomp(k,j,i) = bkappa*th1(k,j,i) + C7(k,j,i) + fth1(k,j,i)
        dudt(this%ACTIVESC+1)%ccomp(k,j,i) = bkappa*th2(k,j,i) + C8(k,j,i) + fth2(k,j,i)
      else
        dudt(this%VELOCITY  )%ccomp(k,j,i) = 0.0_GP
        dudt(this%VELOCITY+1)%ccomp(k,j,i) = 0.0_GP
        dudt(this%VELOCITY+2)%ccomp(k,j,i) = 0.0_GP
        dudt(this%ACTIVESC  )%ccomp(k,j,i) = 0.0_GP
        dudt(this%ACTIVESC+1)%ccomp(k,j,i) = 0.0_GP
      endif
    end do
    end do
    end do

    call this%workspace_%free_complex_tmp(C1)
    call this%workspace_%free_complex_tmp(C2)
    call this%workspace_%free_complex_tmp(C3)
    call this%workspace_%free_complex_tmp(C4)
    call this%workspace_%free_complex_tmp(C5)
    call this%workspace_%free_complex_tmp(C6)
    call this%workspace_%free_complex_tmp(C7)
    call this%workspace_%free_complex_tmp(C8)
    call this%workspace_%free_real_tmp(R1)
    call this%workspace_%free_real_tmp(R2)

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
    use pseudospec_hd
    use pseudospec_phd
    use status
    implicit none

    class(MOISTSolver), intent(in)              :: this
    type  (GStateComp), intent(in), target      :: uin(:), uf(:)
    integer, intent(in)                         :: t
    complex(kind=GP), pointer, dimension(:,:,:) :: fx, fy, fz
    complex(kind=GP), pointer, dimension(:,:,:) :: fth1, fth2
    complex(kind=GP), pointer, dimension(:,:,:) :: vx, vy, vz
    complex(kind=GP), pointer, dimension(:,:,:) :: th1, th2
    complex(kind=GP), pointer, dimension(:,:,:) :: c1, c2, c3
    real(kind=GP)                               :: rmp1, rmq1, rmq2
    integer                                     :: i
    logical                                     :: bret

    call this%workspace_%get_complex_tmp(c1, bret)
    call this%workspace_%get_complex_tmp(c2, bret)
    call this%workspace_%get_complex_tmp(c3, bret)

    vx   => uin(this%VELOCITY  )%ccomp
    vy   => uin(this%VELOCITY+1)%ccomp
    vz   => uin(this%VELOCITY+2)%ccomp
    th1  => uin(this%ACTIVESC  )%ccomp
    th2  => uin(this%ACTIVESC+1)%ccomp
    fx   => uf (this%VELOCITY  )%ccomp
    fy   => uf (this%VELOCITY+1)%ccomp
    fz   => uf (this%VELOCITY+2)%ccomp
    fth1 => uf (this%ACTIVESC  )%ccomp
    fth2 => uf (this%ACTIVESC+1)%ccomp

    call hdcheck(vx,vy,vz,fx,fy,fz,t,dt,1,1,this%todir_)
    call pscheck(th1,fth1,t,dt,this%todir_,trim(this%sstate_(this%ACTIVESC  )))
    call pscheck(th2,fth2,t,dt,this%todir_,trim(this%sstate_(this%ACTIVESC+1)))

    call maxabs(vx,vy,vz,rmp1,0) ! max |vorticity|
    call derivk3(th1, c1, 1)
    call derivk3(th1, c2, 2)
    call derivk3(th1, c3, 3)
    call maxabs(c1,c2,c3,rmq1,2) ! max |Grad bu|
    call derivk3(th2,c1,1)
    call derivk3(th2,c2,2)
    call derivk3(th2,c3,3)
    call maxabs(c1,c2,c3,rmq2,2) ! max |Grad bs|
    if (myrank.eq.0) then
      open(1,file=trim(this%todir_) // '/maximum.txt',position='append')
      write(1,FMT='(E13.6,E13.6,E13.6,E13.6)') (t-1)*dt,rmp1,rmq1,rmq2
      close(1)
    endif

    do i = this%PASSIVE, this%PASSIVE+this%numpassive_-1
      call pscheck(uin(i)%ccomp,uf(i)%ccomp,t,dt,this%todir_,trim(this%sstate_(i)))
    end do   

    call this%workspace_%free_complex_tmp(c1)
    call this%workspace_%free_complex_tmp(c2)
    call this%workspace_%free_complex_tmp(c3)
  end subroutine global_impl


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Function to compute and write spectra
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine spectra_impl(this, uin) 
    use pseudospec_aniso
    use pseudospec_scalar
    use pseudospec_anisca
    use pseudospec_bouss
    use filefmt
    use status
    use iovar
    implicit none

    class (MOISTSolver), intent(in)             :: this
    type(GStateComp), intent(in), target        :: uin(:)
    complex(kind=GP), pointer, dimension(:,:,:) :: vx,vy,vz,th1,th2
    complex(kind=GP), pointer, dimension(:,:,:) :: c1,c2,c3
    real   (kind=GP)                            :: bvfreq,omegaz
    integer                                     :: i
    logical                                     :: bret

    bvfreq = this%traits_%bvfreq 
    omegaz = this%traits_%omega(3) 

    write(ext, fmtext) sind
    vx  => uin(this%VELOCITY  )%ccomp
    vy  => uin(this%VELOCITY+1)%ccomp
    vz  => uin(this%VELOCITY+2)%ccomp
    th1 => uin(this%ACTIVESC  )%ccomp
    th2 => uin(this%ACTIVESC+1)%ccomp
    call spectrum(vx,vy,vz,this%todir_,ext,1,1)
    call this%workspace_%get_complex_tmp(c1,bret)
    call this%workspace_%get_complex_tmp(c2,bret)
    call this%workspace_%get_complex_tmp(c3,bret)
    call gradre3(vx,vy,vz,c1,c2,c3)                      ! Computes v.Grad(v)
    call entrans(vx,vy,vz,-c1,-c2,-c3,this%todir_,ext,1) ! Writes the energy flux
    
    call specpara(vx,vy,vz,this%todir_,ext,1,1)
    call specperp(vx,vy,vz,this%todir_,ext,1,1)
    call spectrsc(th1,this%todir_,ext,0,trim(this%sstate_(this%ACTIVESC  )))
    call specscpa(th1,this%todir_,ext,0,trim(this%sstate_(this%ACTIVESC  )))
    call specscpe(th1,this%todir_,ext,0,trim(this%sstate_(this%ACTIVESC  )))
    call spectrsc(th2,this%todir_,ext,1,trim(this%sstate_(this%ACTIVESC+1)))
    call specscpa(th2,this%todir_,ext,1,trim(this%sstate_(this%ACTIVESC+1)))
    call specscpe(th2,this%todir_,ext,1,trim(this%sstate_(this%ACTIVESC+1)))
    call entpara(vx,vy,vz,-c1,-c2,-c3,this%todir_,ext,1) ! Writes energy flux
    call entperp(vx,vy,vz,-c1,-c2,-c3,this%todir_,ext,1) ! Writes energy flux

    ! Write helicity fluxes, 2D spectra:
    if ( this%traits_%spectlod .ge. 2 ) then
      call heltrans(vx,vy,vz,-c1,-c2,-c3,this%todir_,ext,1)
      ! Write 2D spectra:
      call spec2D  (vx,vy,vz,ext,this%odir_,1,1)
      call specsc2D(th1,ext,this%odir_,0,trim(this%sstate_(this%ACTIVESC  )))
      call specsc2D(th2,ext,this%odir_,0,trim(this%sstate_(this%ACTIVESC+1)))
    endif

    ! Write Fourier modes:
    if ( this%traits_%spectlod .ge. 3 ) then
      call write_fourier(vx, trim(this%sstate_(this%VELOCITY  )),ext,this%odir_)
      call write_fourier(vy, trim(this%sstate_(this%VELOCITY+1)),ext,this%odir_)
      call write_fourier(vz, trim(this%sstate_(this%VELOCITY+2)),ext,this%odir_)
      call write_fourier(th1,trim(this%sstate_(this%ACTIVESC  )),ext,this%odir_)
      call write_fourier(th2,trim(this%sstate_(this%ACTIVESC+1)),ext,this%odir_)
    endif

    call this%workspace_%free_complex_tmp(c1)
    call this%workspace_%free_complex_tmp(c2)
    call this%workspace_%free_complex_tmp(c3)
    if ( this%numpassive_ .gt. 0) then
      do i = this%PASSIVE, this%PASSIVE+this%numpassive_-1
        call spectrsc(uin(i)%ccomp,this%todir_,ext,0,trim(this%sstate_(i)))
        if ( this%traits_%dorot ) then
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
  subroutine MOISTSolver_ctor(this, infile, workspace, plan)
    use iovar
    class(MOISTSolver), intent(inout)        :: this
    type(GWorkspace) , intent(inout), target :: workspace
    type(ioplan)     , intent(inout), target :: plan
    character(len=*), intent(in)             :: infile
    this%infile_    =  infile
    this%workspace_ => workspace
    this%planio_    => plan
    call this%init()
  end subroutine MOISTSolver_ctor


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Destructor
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine MOISTSolver_dtor(this)
    type(MOISTSolver), intent(inout) :: this
    if (associated(this%workspace_))   nullify(this%workspace_)
    if (associated(this%planio_))      nullify(this%planio_)
    if (allocated(this%sstate_))       deallocate(this%sstate_)
    if (allocated(this%traits_%kappa)) deallocate(this%traits_%kappa)
  end subroutine MOISTSolver_dtor


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Convert input state name to index in state vector
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine sstate2istate_impl(this, sstate, istate)
    class(MOISTSolver), intent(in)            :: this
    character(len=8), intent(in)              :: sstate(:)
    integer, allocatable, intent(inout)       :: istate(:)
    integer                                   :: i, j
    if (.not. allocated(istate)) then
      allocate(istate(size(sstate)))
    else if (size(istate) /= size(sstate)) then
      deallocate(istate)
      allocate(istate(size(sstate)))
    end if
    do i = 1, size(sstate)
      istate(i) = -1
      do j = 1, size(this%sstate_)
        if (trim(sstate(i)) == trim(this%sstate_(j))) then
          istate(i) = j
          exit
        end if
      end do
    end do
  end subroutine sstate2istate_impl


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Get state variable names
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_sstate_impl(this, sstate)
    class(MOISTSolver), intent   (in) :: this
    character  (len=8), intent(inout) :: sstate(:)
    character  (len=8)                :: name
    character  (len=1)                :: comp(3)
    integer                           :: j
    comp = ['x', 'y', 'z']
    do j = 1, this%nc_
      sstate(this%VELOCITY + j - 1) = 'v' // comp(j)
    end do
    sstate(this%ACTIVESC  ) = 'bu'
    sstate(this%ACTIVESC+1) = 'bs'
    ! Passive scalars
    do j = 1, this%numpassive_
      write(name,'("s",I0)') j
      sstate(this%PASSIVE + j - 1) = trim(name)
    end do
  end subroutine get_sstate_impl


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Function to compute number of state members (equations)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  pure function state_size_impl(this) result(num)
    class(MOISTSolver), intent(in) :: this
    integer                        :: num
    num = this%nc_ + this%numactivesc_ + this%numpassive_
  end function state_size_impl
end module moist_mod
