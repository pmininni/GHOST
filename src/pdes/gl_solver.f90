! =====================================================================
! NAME       : gl_solver.f90
! DESCRIPTION: Forms class for the Ginzburg-Landau solver (the
!              advective real Ginzburg-Landau equation, ARGL, with or
!              without rotation and a trapping potential), computing:
!
!              dz/dt = omegag.z - beta.|z|^2 z + alpha.Lap(z) - i v.grad(z)
!                      - |v|^2.z/(4.alpha) + f
!                      - i omegaz.(x dz/dy - y dz/dx) - V(x,y).z
!
!              for the complex order parameter z = zre + i zim, with
!              alpha = c.xi/sqrt(2), omegag = c/(xi.sqrt(2)) and beta =
!              omegag/rho0 (c: speed of sound, xi: coherence length, rho0:
!              density at infinity). v is an advective velocity (constant
!              in time, from the initial conditions), f a thermal forcing
!              (the forcing state), V(x,y) = V0.(x^2 + y^2) a trapping
!              potential with V0 = m.w0^2/(2.hbar) (from the forcing
!              factory), and omegaz the rotation rate of the frame. The
!              real and imaginary parts are solved separately:
!                dzre/dt - alpha.Lap(zre) = omegag.zre - beta.|z|^2 zre
!                  + v.grad(zim) - |v|^2.zre/(4.alpha) + fre
!                  + omegaz.(x.dzim/dy - y.dzim/dx) - V.zre
!                dzim/dt - alpha.Lap(zim) = omegag.zim - beta.|z|^2 zim
!                  - v.grad(zre) - |v|^2.zim/(4.alpha) + fim
!                  - omegaz.(x.dzre/dy - y.dzre/dx) - V.zim
!
!              TIME STEPPING: the equations are relaxation equations
!              with a stiff diffusion term, and should be integrated with
!              first order implicit Euler for the Laplacian. To this end
!              dudt returns the effective time derivative of the
!              semi-implicit Euler step,
!                dudt = [(z + dt.N(z))/(1 + alpha.k^2.dt) - z]/dt,
!              with N the explicit terms, so that the traditional stepper
!              of order 1 (z + dt.dudt) performs the implicit step. This
!              solver requires sname='TRADITIONAL' and norder=1 (checked
!              at initialization). In finite temperature runs (kttherm >
!              0) the mass is renormalized to rho0 after each step (or,
!              with a trap, the chemical potential omegag is adjusted to
!              match the mean density rhom).
!
!              State ordering is:
!                zre, zim
!
!              State sector ids are:
!                ZFUNC (ZFUNC+1)
!
! INPUT FILE : For solver='GL', looks for a "&GL" namelist with:
!              fidir   : changes class binary input  dir (default: idir)
!              fodir   : changes class binary output dir (default: odir)
!              todir   : changes the class TXT output dir (default: '')
!              cspeed  : speed of sound
!              lambda  : coherence length
!              rho0    : density at infinity (or at the zero of the
!                        potential, default=1)
!              kttherm : k.T of the thermal bath (>0: finite temperature
!                        run with renormalization of the mass, default=0)
!              rhom    : mean density for finite temperature runs with a
!                        trap (default=rho0)
!              V0      : amplitude of the trapping potential (default=0;
!                        the potential itself is set by the forcing method
!                        'cyltrap_vq')
!              docflow : .true.=counterflow run, the |v|^2 term is dropped
!                        (default=.false.)
!              dorot   : do rotation, = .TRUE. or .FALSE.
!              omegaz  : rotation rate along z (rotating frame)
!              spectlod: spectral output level of detail (in [1,3])
!
!              The advective velocity v is set by the initial conditions
!              of the order parameter (e.g., 'abc_z', 'tg_z', 'ring_z'),
!              the thermal forcing by 'thermal_fq', and the trapping
!              potential by 'cyltrap_vq' (see ic_quantum.f90 and
!              force_quantum.f90).
!
! DATE       : 09/11/26 (PDM)
! =====================================================================

module gl_mod
  use equationbase_mod
  use gstate_mod
  use gdevice, only: gdev_active

  implicit none

  ! ================= Solver traits ===================================
  type, public  :: GLTraits
    logical       :: dorot    = .FALSE. ! rotation flag
    logical       :: docflow  = .FALSE. ! counterflow (no |v|^2 term)
    integer       :: spectlod = 1       ! standard level of spectra detail
    real(kind=GP) :: cspeed   = 1.0_GP  ! speed of sound
    real(kind=GP) :: lambda   = 1.0_GP  ! coherence length
    real(kind=GP) :: rho0     = 1.0_GP  ! density at infinity
    real(kind=GP) :: kttherm  = 0.0_GP  ! k.T of the thermal bath
    real(kind=GP) :: rhom     = 1.0_GP  ! mean density (trap, finite T)
    real(kind=GP) :: V0       = 0.0_GP  ! trapping potential amplitude
    real(kind=GP) :: omegaz   = 0.0_GP  ! rotation rate
  end type

  ! ================= Global parameters ===============================
  integer, parameter, public   :: MAXPASSIVE = 0 ! No passive scalars

  ! ================= Solver ==========================================
  ! Define class:
  type, extends(QuantumBase) :: GLSolver
    ! Member data:
    logical           :: binit_  = .false. ! is initialized?
    type  (GLTraits)  :: traits_
  CONTAINS
    procedure, public :: init          =>          init_impl ! init method
    procedure, public :: dudt          =>          dudt_impl ! RHS method
    procedure, public :: global        =>        global_impl ! Writes global qtys
    procedure, public :: spectra       =>       spectra_impl ! Writes spectra
    procedure, public :: state_size    =>    state_size_impl ! state size
    procedure, public :: sstate2istate => sstate2istate_impl ! state names
    procedure, public :: get_sstate    =>    get_sstate_impl ! get state name list
    procedure, public :: sync_device   =>   sync_device_impl ! aux. arrays to device
    procedure, public :: Solver_ctor   =>      GLSolver_ctor ! constructor
    final             :: GLSolver_dtor
  end type GLSolver

CONTAINS

  ! ===================================================================
  ! Solver initialization, this is where parameter files are read
  ! ===================================================================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Subroutine to initialize the solver
  !! Reads the &GL namelist, sets the constants of the equation
  !! and the sector indices, and checks the time stepper
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_impl(this)
    use commtypes
    use status
    class  (GLSolver), intent (inout) :: this

    ! Temporary data to read from namelists:
    logical                    :: dorot, docflow
    integer                    :: spectlod
    integer                    :: ierr, norder, nstage, itype
    real(kind=GP)              :: cspeed, lambda, rho0, kttherm
    real(kind=GP)              :: rhom, V0, omegaz
    character(len=128)         :: fidir, fodir, todir, sname

    ! Required namelists:
    namelist/ GL      / fidir, fodir, todir, cspeed, lambda, rho0
    namelist/ GL      / kttherm, rhom, V0, docflow, dorot, omegaz, spectlod
    namelist/ stepper / sname, itype, norder, nstage

    call MPI_COMM_SIZE(MPI_COMM_WORLD,this%nprocs_,ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,this%myrank_,ierr)

    ! Get I/O and trait variables from input file:
    fidir    = idir ! Set the default to status idir
    fodir    = odir ! Set the default to status odir
    todir    = '.'  ! Set the default to the current dir
    dorot    = .FALSE.
    docflow  = .FALSE.
    spectlod = 1 ! standard lod
    cspeed   = 1.0_GP
    lambda   = 1.0_GP
    rho0     = 1.0_GP
    kttherm  = 0.0_GP
    rhom     = -1.0_GP ! defaults to rho0 below
    V0       = 0.0_GP
    omegaz   = 0.0_GP
    sname    = 'TRADITIONAL'; itype = 1; norder = 2; nstage = 2
    if ( this%myrank_ .eq. 0 ) then
      open(1,file=this%infile_,status='unknown',form="formatted")
      read(1,NML=GL)
      close(1)
      open(1,file=this%infile_,status='unknown',form="formatted")
      read(1,NML=stepper)
      close(1)
    endif
    call MPI_BCAST(fidir    ,128,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(fodir    ,128,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(todir    ,128,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(sname    ,128,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(cspeed   ,1  ,GC_REAL      ,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(lambda   ,1  ,GC_REAL      ,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(rho0     ,1  ,GC_REAL      ,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(kttherm  ,1  ,GC_REAL      ,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(rhom     ,1  ,GC_REAL      ,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(V0       ,1  ,GC_REAL      ,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(omegaz   ,1  ,GC_REAL      ,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(dorot    ,1  ,MPI_LOGICAL  ,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(docflow  ,1  ,MPI_LOGICAL  ,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(spectlod ,1  ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(norder   ,1  ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
    if ( rhom .lt. 0.0_GP ) rhom = rho0

    ! The semi-implicit Euler step is built inside dudt (see the header)
    if ( (trim(adjustl(sname)).ne.'TRADITIONAL') .or. (norder.ne.1) ) then
      if ( this%myrank_ .eq. 0 ) then
        write(*,*) 'GLSolver::init: the GL solver requires sname = TRADITIONAL', &
                   ' and norder = 1 in the stepper namelist'
      endif
      stop
    endif

    ! Set I/O and traits from inputfile data:
    this%idir_  = fidir ! If present in &GL, replaces the class default idir
    this%odir_  = fodir ! If present in &GL, replaces the class default odir
    this%todir_ = todir ! If present in &GL, replaces the class default todir
    this%traits_%   dorot = dorot
    this%traits_% docflow = docflow
    this%traits_%spectlod = spectlod
    this%traits_%  cspeed = cspeed
    this%traits_%  lambda = lambda
    this%traits_%    rho0 = rho0
    this%traits_% kttherm = kttherm
    this%traits_%    rhom = rhom
    this%traits_%      V0 = V0
    this%traits_%  omegaz = omegaz
    ! Constants of the equation (shared with the ICs and the forcing)
    this%alpha_  = cspeed*lambda/sqrt(2.0_GP)
    this%omegag_ = cspeed/(lambda*sqrt(2.0_GP))
    this%beta_   = this%omegag_/rho0
    this%rho0_   = rho0
    this%V0_     = V0
    this%omegaz_ = omegaz
    this%dorot_  = dorot

    this%ZFUNC = 1                            ! start of the order parameter
    allocate(this%sstate_(this%state_size()))
    call this%get_sstate(this%sstate_)
    this%binit_ = .true.
  end subroutine init_impl


  ! ===================================================================
  ! Computation of RHS, the solver equations are defined here
  ! ===================================================================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Function to compute the effective RHS: dudt = (z_new - z)/dt
  !! with z_new the semi-implicit Euler update (see the header).
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine dudt_impl(this, time, uin, uf, dt, dudt)
    use pseudospec_fluid
    use pseudospec_scalar, only: advect3
    use pseudospec_quantum
    use ali
    use kes
    use var
    use grid
    use mpivars
    implicit none

    class (GLSolver), intent(inout)             :: this
    real   (kind=GP), intent   (in)             :: time, dt
    type(GStateComp), intent(inout), target     :: uin(:),uf(:)
    type(GStateComp), intent(inout), target     :: dudt(:)
    complex(kind=GP), pointer, dimension(:,:,:) :: zre,zim,fre,fim,dre,dim
    complex(kind=GP), pointer, dimension(:,:,:) :: C1,C2,C3,C4,C5,C6,C7,C8
    real   (kind=GP), pointer, dimension(:,:,:) :: R1
    real   (kind=GP)                            :: alpha,beta,omegag,rmq
    double precision                            :: mass
    integer                                     :: iadv,irot
    logical                                     :: bret

    if ( .not. this%binit_ ) then
      stop 'GLSolver::dudt: Solver not initialized'
    endif

    alpha  = this%alpha_
    beta   = this%beta_
    omegag = this%omegag_
    iadv = 0; if ( this%hasadv_ ) iadv = 1
    irot = 0; if ( this%dorot_  ) irot = 1

    call this%workspace_%get_complex_tmp(C1,bret)
    call this%workspace_%get_complex_tmp(C2,bret)
    call this%workspace_%get_complex_tmp(C3,bret)
    call this%workspace_%get_complex_tmp(C4,bret)
    call this%workspace_%get_complex_tmp(C5,bret)
    call this%workspace_%get_complex_tmp(C6,bret)
    if ( irot.eq.1 ) then
      call this%workspace_%get_complex_tmp(C7,bret)
      call this%workspace_%get_complex_tmp(C8,bret)
    endif
    call this%workspace_%get_real_tmp(R1,bret)

    zre => uin (this%ZFUNC  )%ccomp
    zim => uin (this%ZFUNC+1)%ccomp
    fre => uf  (this%ZFUNC  )%ccomp
    fim => uf  (this%ZFUNC+1)%ccomp
    dre => dudt(this%ZFUNC  )%ccomp
    dim => dudt(this%ZFUNC+1)%ccomp

    ! Nonlinear and potential terms: (omegag/beta - |z|^2 - |v|^2/(4
    ! alpha beta) - V/beta).z, multiplied by beta below. The square of
    ! the advective velocity and the potential are class arrays, used
    ! through explicit-shape dummies in the kernels.
    rmq = real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP)*omegag/beta
    call squareabs(zre,zim,R1,1)
    call gl_potential(R1,this%vsq_,this%vpot_,rmq, &
                      (this%hasadv_.and..not.this%traits_%docflow),this%haspot_)
    call nonlgpe(R1,zre,C3)
    call nonlgpe(R1,zim,C4)
    ! Advection: -v.grad(z)
    if ( iadv.eq.1 ) then
      call advect3(this%vx_,this%vy_,this%vz_,zre,C5)
      call advect3(this%vx_,this%vy_,this%vz_,zim,C6)
    endif
    ! Rotation: omegaz.(x dz/dy - y dz/dx) for zre (C7) and zim (C8)
    if ( irot.eq.1 ) then
      call derivk3(zre,C7,1)                ! dzre/dx
      call nonlgpe(this%vliny_,C7,C1)       ! omegaz y dzre/dx
      call derivk3(zre,C7,2)                ! dzre/dy
      call nonlgpe(this%vlinx_,C7,C2)       ! omegaz x dzre/dy
      call saxpby_c(C7,C2,1.0_GP,C1,-1.0_GP)
      call derivk3(zim,C8,1)
      call nonlgpe(this%vliny_,C8,C1)
      call derivk3(zim,C8,2)
      call nonlgpe(this%vlinx_,C8,C2)
      call saxpby_c(C8,C2,1.0_GP,C1,-1.0_GP)
    endif
    ! Semi-implicit Euler step into C1 (zre) and C2 (zim)
    if ( irot.eq.1 ) then
      call gl_step(zre,zim,C3,C4,C5,C6,C7,C8,fre,fim,C1,C2,alpha,beta,dt,iadv,irot)
    else
      call gl_step(zre,zim,C3,C4,C5,C6,C5,C6,fre,fim,C1,C2,alpha,beta,dt,iadv,irot)
    endif
    ! Renormalization in finite temperature runs: the mass is set to
    ! rho0, or with a trap the chemical potential is adjusted to match
    ! the mean density rhom
    if ( this%traits_%kttherm .gt. 0.0_GP ) then
      call gpe_mass(C1,C2,mass)
      if ( this%V0_ .le. tiny ) then
        rmq = sqrt(omegag/beta)/sqrt(real(mass,kind=GP))
        call scal3(C1,rmq)
        call scal3(C2,rmq)
      else
        this%omegag_ = this%omegag_ - 150.0_GP*(real(mass,kind=GP)-this%traits_%rhom)
      endif
    endif
    ! Effective time derivative for the stepper
    call saxpby_c(dre,C1,1.0_GP/dt,zre,-1.0_GP/dt)
    call saxpby_c(dim,C2,1.0_GP/dt,zim,-1.0_GP/dt)

    call this%workspace_%free_real_tmp(R1)
    if ( irot.eq.1 ) then
      call this%workspace_%free_complex_tmp(C8)
      call this%workspace_%free_complex_tmp(C7)
    endif
    call this%workspace_%free_complex_tmp(C6)
    call this%workspace_%free_complex_tmp(C5)
    call this%workspace_%free_complex_tmp(C4)
    call this%workspace_%free_complex_tmp(C3)
    call this%workspace_%free_complex_tmp(C2)
    call this%workspace_%free_complex_tmp(C1)
  end subroutine dudt_impl


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Kernel: r = rmq - r [- vsq] [- vpot] in real space (the
  !! argument of the nonlinear term). Module procedure with
  !! explicit-shape dummies so that the class arrays can be used
  !! in the device kernel.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine gl_potential(r,vsq,vpot,rmq,useadv,usepot)
    use grid
    use mpivars
    implicit none
    real(kind=GP), intent(inout), dimension(nx,ny,ksta:kend) :: r
    real(kind=GP), intent   (in), dimension(nx,ny,ksta:kend) :: vsq,vpot
    real(kind=GP), intent   (in)                             :: rmq
    logical      , intent   (in)                             :: useadv,usepot
    real(kind=GP)                                            :: ca,cp
    integer                                                  :: i,j,k
    ca = 0.0_GP; if ( useadv ) ca = 1.0_GP
    cp = 0.0_GP; if ( usepot ) cp = 1.0_GP
#if defined(GHOST_GPU)
!$omp target teams distribute parallel do collapse(3) if(target: gdev_active)
#else
!$omp parallel do collapse(2) private (i)
#endif
    do k = ksta,kend
      do j = 1,ny
        do i = 1,nx
          r(i,j,k) = rmq - r(i,j,k) - ca*vsq(i,j,k) - cp*vpot(i,j,k)
        end do
      end do
    end do
  end subroutine gl_potential


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Kernel: semi-implicit Euler step of the two components,
  !!   nre = (zre + dt (beta.c3 - c6 + fre + c8))/(1 + alpha k^2 dt)
  !!   nim = (zim + dt (beta.c4 + c5 + fim - c7))/(1 + alpha k^2 dt)
  !! for the modes with kn2 <= kmax (zero otherwise). The advection
  !! (c5,c6) and rotation (c7,c8) terms are included when iadv,
  !! irot = 1.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine gl_step(zre,zim,c3,c4,c5,c6,c7,c8,fre,fim,nre,nim,alpha,beta,dt,iadv,irot)
    use grid
    use mpivars
    use kes
    use ali
    implicit none
    complex(kind=GP), intent (in), dimension(nz,ny,ista:iend) :: zre,zim,c3,c4
    complex(kind=GP), intent (in), dimension(nz,ny,ista:iend) :: c5,c6,c7,c8
    complex(kind=GP), intent (in), dimension(nz,ny,ista:iend) :: fre,fim
    complex(kind=GP), intent(out), dimension(nz,ny,ista:iend) :: nre,nim
    real   (kind=GP), intent (in)                             :: alpha,beta,dt
    integer         , intent (in)                             :: iadv,irot
    real   (kind=GP)                                          :: ca,cr,rmp
    integer                                                   :: i,j,k
    ca = real(iadv,kind=GP)
    cr = real(irot,kind=GP)
#if defined(GHOST_GPU)
!$omp target teams distribute parallel do collapse(3) if(target: gdev_active) private(rmp)
#else
!$omp parallel do collapse(2) private (k,rmp)
#endif
    do i = ista,iend
      do j = 1,ny
        do k = 1,nz
          if (kn2(k,j,i).le.kmax) then
            rmp = 1.0_GP/(1.0_GP+alpha*kk2(k,j,i)*dt)
            nre(k,j,i) = (zre(k,j,i)+dt*(beta*c3(k,j,i)-ca*c6(k,j,i) &
                         +fre(k,j,i)+cr*c8(k,j,i)))*rmp
            nim(k,j,i) = (zim(k,j,i)+dt*(beta*c4(k,j,i)+ca*c5(k,j,i) &
                         +fim(k,j,i)-cr*c7(k,j,i)))*rmp
          else
            nre(k,j,i) = 0.0_GP
            nim(k,j,i) = 0.0_GP
          endif
        end do
      end do
    end do
  end subroutine gl_step


  ! ===================================================================
  ! Computation of global quantities and spectra
  ! ===================================================================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Function to compute and write global quantities
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine global_impl(this, uin, uf, t)
    use pseudospec_gpe
    use pseudospec_rgpe
    use status
    implicit none

    class (GLSolver), intent(in)                :: this
    type(GStateComp), intent(in), target        :: uin(:), uf(:)
    integer         , intent(in)                :: t
    complex(kind=GP), pointer, dimension(:,:,:) :: zre,zim

    zre => uin(this%ZFUNC  )%ccomp
    zim => uin(this%ZFUNC+1)%ccomp
    call gpecheck(zre,zim,this%alpha_,this%beta_,this%omegag_,t,dt,this%todir_)
    call momentum(zre,zim,this%alpha_,t,dt,this%todir_)
    call gpehelicity(zre,zim,this%alpha_,this%beta_,this%omegag_,t,dt,this%todir_)
    if ( this%haspot_ ) then
      call trapenergy(zre,zim,this%vpot_,this%alpha_,this%beta_,t,dt,this%todir_)
    endif
    if ( this%dorot_ ) then
      call rotenergy(zre,zim,this%vlinx_,this%vliny_,this%alpha_,t,dt,this%todir_)
    endif
  end subroutine global_impl


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Function to compute and write spectra
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine spectra_impl(this, uin)
    use pseudospec_gpe
    use filefmt
    use status
    implicit none

    class (GLSolver), intent(in)                :: this
    type(GStateComp), intent(in), target        :: uin(:)
    complex(kind=GP), pointer, dimension(:,:,:) :: zre,zim

    write(ext, fmtext) sind
    zre => uin(this%ZFUNC  )%ccomp
    zim => uin(this%ZFUNC+1)%ccomp
    call gpemassspec(zre,zim,this%todir_,ext)
    call gperealspec(zre,zim,this%alpha_,this%beta_,this%omegag_,this%todir_,ext)
    call gpehelspec (zre,zim,this%alpha_,this%beta_,this%omegag_,this%todir_,ext)
    call gpemomtspec(zre,zim,this%alpha_,this%todir_,ext)
  end subroutine spectra_impl


  ! ===================================================================
  ! Solver specific methods
  ! ===================================================================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Auxiliary arrays to the device: computes |v|^2/(4 alpha beta)
  !! (dealiased) from the advective velocity set by the initial
  !! conditions, and copies all the auxiliary arrays to the
  !! device. Called by the main program after the initial
  !! conditions and the forcing are set (on the host).
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine sync_device_impl(this)
    use pseudospec_fluid
    use grid
    use mpivars
    use kes
    use ali
    use fft
    use commtypes
    implicit none
    class (GLSolver), intent(inout)             :: this
    complex(kind=GP), pointer, dimension(:,:,:) :: C1
    real   (kind=GP), pointer, dimension(:,:,:) :: R1,R2,R3
    real   (kind=GP)                            :: rmp
    integer                                     :: i,j,k
    logical                                     :: bret

    if ( this%hasadv_ ) then
      call this%workspace_%get_complex_tmp(C1,bret)
      call this%workspace_%get_real_tmp(R1,bret)
      call this%workspace_%get_real_tmp(R2,bret)
      call this%workspace_%get_real_tmp(R3,bret)
      call copy3(this%vx_,C1)
      call fftp3d_complex_to_real(plancr,C1,R1,MPI_COMM_WORLD)
      call copy3(this%vy_,C1)
      call fftp3d_complex_to_real(plancr,C1,R2,MPI_COMM_WORLD)
      call copy3(this%vz_,C1)
      call fftp3d_complex_to_real(plancr,C1,R3,MPI_COMM_WORLD)
      rmp = 1.0_GP/(4*this%alpha_*this%beta_*(real(nx,kind=GP)* &
            real(ny,kind=GP)*real(nz,kind=GP))**2)
!$omp parallel do collapse(2) private (i)
      do k = ksta,kend
        do j = 1,ny
          do i = 1,nx
            this%vsq_(i,j,k) = (R1(i,j,k)**2+R2(i,j,k)**2+R3(i,j,k)**2)*rmp
          end do
        end do
      end do
      ! Dealiases |v|^2 (the result is not normalized)
      call fftp3d_real_to_complex(planrc,this%vsq_,C1,MPI_COMM_WORLD)
!$omp parallel do collapse(2) private (k)
      do i = ista,iend
        do j = 1,ny
          do k = 1,nz
            if (kn2(k,j,i).gt.kmax) C1(k,j,i) = 0.0_GP
          end do
        end do
      end do
      call fftp3d_complex_to_real(plancr,C1,this%vsq_,MPI_COMM_WORLD)
      call this%workspace_%free_real_tmp(R3)
      call this%workspace_%free_real_tmp(R2)
      call this%workspace_%free_real_tmp(R1)
      call this%workspace_%free_complex_tmp(C1)
    endif
    call this%aux_to_device()
  end subroutine sync_device_impl


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Constructor: reads the parameters and allocates the
  !! auxiliary arrays (zero by default)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine GLSolver_ctor(this, infile, workspace, plan)
    use iovar
    class  (GLSolver), intent(inout)         :: this
    type(GWorkspace) , intent(inout), target :: workspace
    type(ioplan)     , intent(inout), target :: plan
    character(len=*) , intent   (in)         :: infile
    this%infile_    =  infile    ! input file
    this%workspace_ => workspace
    this%planio_    => plan
    call this%init()
    call this%alloc_aux()
  end subroutine GLSolver_ctor


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Destructor
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine GLSolver_dtor(this)
    type  (GLSolver), intent(inout) :: this
    if (associated(this%workspace_))   nullify(this%workspace_)
    if (associated(this%planio_))      nullify(this%planio_)
    if (allocated(this%sstate_))       deallocate(this%sstate_)
    call this%free_aux()
  end subroutine GLSolver_dtor


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Convert input state name to index in state vector
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine sstate2istate_impl(this, sstate, istate)
    class (GLSolver), intent   (in) :: this
    character(len=8), intent   (in) :: sstate(:)
    integer         , intent(inout) :: istate(:)
    integer                         :: i,j
    if ( size(sstate) .ne. size(istate) ) then
      stop 'GLSolver::sstate2istate_impl: Incompatible sstate and istate'
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
    class  (GLSolver), intent   (in) :: this
    character (len=8), intent(inout) :: sstate(:)
    sstate(this%ZFUNC  ) = 'zre'
    sstate(this%ZFUNC+1) = 'zim'
  end subroutine get_sstate_impl


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Function to compute number of state members (equations)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PURE function state_size_impl(this) result(num)
    class(GLSolver), intent(in) :: this
    integer                     :: num
    num = 2                      ! real and imaginary parts of z
  end function state_size_impl

end module gl_mod
