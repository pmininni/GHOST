! =====================================================================
! NAME       : gpe_solver.f90
! DESCRIPTION: Forms class for the Gross-Pitaevskii solver (GPE, with
!              or without rotation and a trapping potential), computing:
!
!              dz/dt = i [omegag.z - beta.|z|^2 z + alpha.Lap(z)
!                         - i omegaz.(x dz/dy - y dz/dx) - V(x,y).z]
!
!              for the complex order parameter z = zre + i zim, with
!              alpha = c.xi/sqrt(2), omegag = c/(xi.sqrt(2)) and beta =
!              omegag/rho0 (c: speed of sound, xi: coherence length, rho0:
!              density at infinity). V(x,y) = V0.(x^2 + y^2) is a trapping
!              potential with V0 = m.w0^2/(2.hbar) (w0 the trapping
!              frequency, set by the forcing factory), and omegaz the
!              rotation rate of the frame. The real and imaginary parts
!              are solved separately:
!                dzre/dt = - omegag.zim + beta.|z|^2 zim - alpha.Lap(zim)
!                          + omegaz.(x.dzre/dy - y.dzre/dx) + V.zim
!                dzim/dt =   omegag.zre - beta.|z|^2 zre + alpha.Lap(zre)
!                          + omegaz.(x.dzim/dy - y.dzim/dx) - V.zre
!
!              The equations are Hamiltonian and are integrated with the
!              explicit Runge-Kutta steppers (TRADITIONAL or GEXRK). The
!              forcing state is not used by this solver (the potential
!              is the only external field, set with 'cyltrap_vq'); use
!              'null_fq,constant_fq' in the forcing namelist.
!
!              State ordering is:
!                zre, zim
!
!              State sector ids are:
!                ZFUNC (ZFUNC+1)
!
! INPUT FILE : For solver='GPE', looks for a "&GPE" namelist with:
!              fidir   : changes class binary input  dir (default: idir)
!              fodir   : changes class binary output dir (default: odir)
!              todir   : changes the class TXT output dir (default: '')
!              cspeed  : speed of sound
!              lambda  : coherence length
!              rho0    : density at infinity (or at the zero of the
!                        potential, default=1)
!              V0      : amplitude of the trapping potential (default=0;
!                        the potential itself is set by the forcing method
!                        'cyltrap_vq')
!              dorot   : do rotation, = .TRUE. or .FALSE.
!              omegaz  : rotation rate along z (rotating frame)
!              dotrans : .true. computes the energy transfer functions
!                        at the times of the spectra (default=.false.):
!                        the spectra of the state advanced one time step
!                        with a second order Runge-Kutta step are
!                        compared with the spectra of the state
!              spectlod: spectral output level of detail (in [1,3])
!
!              The initial conditions of the order parameter are the
!              same as for the GL solver ('read_z', 'abc_z', 'ring_z',
!              etc., see ic_quantum.f90); the advective velocity they
!              define is ignored by this solver. The trapping potential
!              is set by 'cyltrap_vq' (see force_quantum.f90).
!
! DATE       : 09/11/26 (PDM)
! =====================================================================

module gpe_mod
  use equationbase_mod
  use gstate_mod
  use gdevice, only: gdev_active

  implicit none

  ! ================= Solver traits ===================================
  type, public  :: GPETraits
    logical       :: dorot    = .FALSE. ! rotation flag
    logical       :: dotrans  = .FALSE. ! energy transfer functions
    integer       :: spectlod = 1       ! standard level of spectra detail
    real(kind=GP) :: cspeed   = 1.0_GP  ! speed of sound
    real(kind=GP) :: lambda   = 1.0_GP  ! coherence length
    real(kind=GP) :: rho0     = 1.0_GP  ! density at infinity
    real(kind=GP) :: V0       = 0.0_GP  ! trapping potential amplitude
    real(kind=GP) :: omegaz   = 0.0_GP  ! rotation rate
  end type

  ! ================= Global parameters ===============================
  integer, parameter, public   :: MAXPASSIVE = 0 ! No passive scalars

  ! ================= Solver ==========================================
  ! Define class:
  type, extends(QuantumBase) :: GPESolver
    ! Member data:
    logical           :: binit_  = .false. ! is initialized?
    type  (GPETraits) :: traits_
  CONTAINS
    procedure, public :: init          =>          init_impl ! init method
    procedure, public :: dudt          =>          dudt_impl ! RHS method
    procedure, public :: global        =>        global_impl ! Writes global qtys
    procedure, public :: spectra       =>       spectra_impl ! Writes spectra
    procedure, public :: state_size    =>    state_size_impl ! state size
    procedure, public :: sstate2istate => sstate2istate_impl ! state names
    procedure, public :: get_sstate    =>    get_sstate_impl ! get state name list
    procedure, public :: sync_device   =>   sync_device_impl ! aux. arrays to device
    procedure, public :: Solver_ctor   =>     GPESolver_ctor ! constructor
    final             :: GPESolver_dtor
  end type GPESolver

CONTAINS

  ! ===================================================================
  ! Solver initialization, this is where parameter files are read
  ! ===================================================================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Subroutine to initialize the solver
  !! Reads the &GPE namelist, sets the constants of the equation
  !! and the sector indices
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_impl(this)
    use commtypes
    use status
    class  (GPESolver), intent (inout) :: this

    ! Temporary data to read from namelists:
    logical                    :: dorot, dotrans
    integer                    :: spectlod
    integer                    :: ierr
    real(kind=GP)              :: cspeed, lambda, rho0, V0, omegaz
    character(len=128)         :: fidir, fodir, todir

    ! Required namelists:
    namelist/ GPE     / fidir, fodir, todir, cspeed, lambda, rho0
    namelist/ GPE     / V0, dorot, omegaz, dotrans, spectlod

    call MPI_COMM_SIZE(MPI_COMM_WORLD,this%nprocs_,ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,this%myrank_,ierr)

    ! Get I/O and trait variables from input file:
    fidir    = idir ! Set the default to status idir
    fodir    = odir ! Set the default to status odir
    todir    = '.'  ! Set the default to the current dir
    dorot    = .FALSE.
    dotrans  = .FALSE.
    spectlod = 1 ! standard lod
    cspeed   = 1.0_GP
    lambda   = 1.0_GP
    rho0     = 1.0_GP
    V0       = 0.0_GP
    omegaz   = 0.0_GP
    if ( this%myrank_ .eq. 0 ) then
      open(1,file=this%infile_,status='unknown',form="formatted")
      read(1,NML=GPE)
      close(1)
    endif
    call MPI_BCAST(fidir    ,128,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(fodir    ,128,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(todir    ,128,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(cspeed   ,1  ,GC_REAL      ,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(lambda   ,1  ,GC_REAL      ,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(rho0     ,1  ,GC_REAL      ,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(V0       ,1  ,GC_REAL      ,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(omegaz   ,1  ,GC_REAL      ,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(dorot    ,1  ,MPI_LOGICAL  ,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(dotrans  ,1  ,MPI_LOGICAL  ,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(spectlod ,1  ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)

    ! Set I/O and traits from inputfile data:
    this%idir_  = fidir ! If present in &GPE, replaces the class default idir
    this%odir_  = fodir ! If present in &GPE, replaces the class default odir
    this%todir_ = todir ! If present in &GPE, replaces the class default todir
    this%traits_%   dorot = dorot
    this%traits_% dotrans = dotrans
    this%traits_%spectlod = spectlod
    this%traits_%  cspeed = cspeed
    this%traits_%  lambda = lambda
    this%traits_%    rho0 = rho0
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
  !! Function to compute the RHS of the GPE
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine dudt_impl(this, time, uin, uf, dt, dudt)
    implicit none

    class (GPESolver), intent(inout)            :: this
    real   (kind=GP), intent   (in)             :: time, dt
    type(GStateComp), intent(inout), target     :: uin(:),uf(:)
    type(GStateComp), intent(inout), target     :: dudt(:)
    complex(kind=GP), pointer, dimension(:,:,:) :: zre,zim,dre,dim

    if ( .not. this%binit_ ) then
      stop 'GPESolver::dudt: Solver not initialized'
    endif

    zre => uin (this%ZFUNC  )%ccomp
    zim => uin (this%ZFUNC+1)%ccomp
    dre => dudt(this%ZFUNC  )%ccomp
    dim => dudt(this%ZFUNC+1)%ccomp
    call gpe_rhs(this,zre,zim,dre,dim)
  end subroutine dudt_impl


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Computes the RHS of the two components of the GPE (see the
  !! header) for the order parameter (zre,zim), into (dre,dim).
  !! Shared by dudt and by the computation of the transfer
  !! functions.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine gpe_rhs(this, zre, zim, dre, dim)
    use pseudospec_fluid
    use pseudospec_quantum
    use grid
    use mpivars
    implicit none

    class (GPESolver), intent   (in)                          :: this
    complex(kind=GP) , intent   (in), dimension(nz,ny,ista:iend) :: zre,zim
    complex(kind=GP) , intent(inout), dimension(nz,ny,ista:iend) :: dre,dim ! may alias zre,zim
    complex(kind=GP) , pointer, dimension(:,:,:) :: C1,C2,C3,C4,C5,C6,C7,C8
    real   (kind=GP) , pointer, dimension(:,:,:) :: R1
    integer                                      :: irot
    logical                                      :: bret

    irot = 0; if ( this%dorot_ ) irot = 1

    call this%workspace_%get_complex_tmp(C3,bret)
    call this%workspace_%get_complex_tmp(C4,bret)
    call this%workspace_%get_complex_tmp(C5,bret)
    call this%workspace_%get_complex_tmp(C6,bret)
    if ( irot.eq.1 ) then
      call this%workspace_%get_complex_tmp(C1,bret)
      call this%workspace_%get_complex_tmp(C2,bret)
      call this%workspace_%get_complex_tmp(C7,bret)
      call this%workspace_%get_complex_tmp(C8,bret)
    endif
    call this%workspace_%get_real_tmp(R1,bret)

    ! Nonlinear and potential terms: (|z|^2 + V/beta).z, multiplied
    ! by beta below. The potential is a class array, used through
    ! explicit-shape dummies in the kernel.
    call squareabs(zre,zim,R1,1)
    if ( this%haspot_ ) call gpe_addpot(R1,this%vpot_)
    call nonlgpe(R1,zre,C3)
    call nonlgpe(R1,zim,C4)
    call laplak3(zre,C5)                    ! Lap(zre)
    call laplak3(zim,C6)                    ! Lap(zim)
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
      call gpe_rhs_kernel(zre,zim,C3,C4,C5,C6,C7,C8,dre,dim, &
                          this%alpha_,this%beta_,this%omegag_,irot)
    else
      call gpe_rhs_kernel(zre,zim,C3,C4,C5,C6,C5,C6,dre,dim, &
                          this%alpha_,this%beta_,this%omegag_,irot)
    endif

    call this%workspace_%free_real_tmp(R1)
    if ( irot.eq.1 ) then
      call this%workspace_%free_complex_tmp(C8)
      call this%workspace_%free_complex_tmp(C7)
      call this%workspace_%free_complex_tmp(C2)
      call this%workspace_%free_complex_tmp(C1)
    endif
    call this%workspace_%free_complex_tmp(C6)
    call this%workspace_%free_complex_tmp(C5)
    call this%workspace_%free_complex_tmp(C4)
    call this%workspace_%free_complex_tmp(C3)
  end subroutine gpe_rhs


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Kernel: r = r + vpot in real space. Module procedure with
  !! explicit-shape dummies so that the class array can be used
  !! in the device kernel.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine gpe_addpot(r,vpot)
    use grid
    use mpivars
    implicit none
    real(kind=GP), intent(inout), dimension(nx,ny,ksta:kend) :: r
    real(kind=GP), intent   (in), dimension(nx,ny,ksta:kend) :: vpot
    integer                                                  :: i,j,k
#if defined(GHOST_GPU)
!$omp target teams distribute parallel do collapse(3) if(target: gdev_active)
#else
!$omp parallel do collapse(2) private (i)
#endif
    do k = ksta,kend
      do j = 1,ny
        do i = 1,nx
          r(i,j,k) = r(i,j,k) + vpot(i,j,k)
        end do
      end do
    end do
  end subroutine gpe_addpot


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Kernel: RHS of the two components,
  !!   dre = -omegag.zim + beta.c4 - alpha.c6 [+ c7]
  !!   dim =  omegag.zre - beta.c3 + alpha.c5 [+ c8]
  !! for the modes with kn2 <= kmax (zero otherwise). The rotation
  !! terms (c7,c8) are included when irot = 1. The steppers can
  !! call dudt in place (dre and dim overwriting zre and zim), so
  !! the terms with zre and zim are first accumulated in the
  !! temporaries c3 and c4 (overwritten), and the output is then
  !! built from the temporaries only.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine gpe_rhs_kernel(zre,zim,c3,c4,c5,c6,c7,c8,dre,dim,alpha,beta,omegag,irot)
    use grid
    use mpivars
    use kes
    use ali
    implicit none
    complex(kind=GP), intent   (in), dimension(nz,ny,ista:iend) :: zre,zim
    complex(kind=GP), intent(inout), dimension(nz,ny,ista:iend) :: c3,c4
    complex(kind=GP), intent   (in), dimension(nz,ny,ista:iend) :: c5,c6,c7,c8
    complex(kind=GP), intent  (out), dimension(nz,ny,ista:iend) :: dre,dim
    real   (kind=GP), intent   (in)                             :: alpha,beta,omegag
    integer         , intent   (in)                             :: irot
    real   (kind=GP)                                            :: cr
    integer                                                     :: i,j,k
    cr = real(irot,kind=GP)
#if defined(GHOST_GPU)
!$omp target teams distribute parallel do collapse(3) if(target: gdev_active)
#else
!$omp parallel do collapse(2) private (k)
#endif
    do i = ista,iend
      do j = 1,ny
        do k = 1,nz
          c4(k,j,i) = -omegag*zim(k,j,i)+beta*c4(k,j,i)-alpha*c6(k,j,i)
          c3(k,j,i) =  omegag*zre(k,j,i)-beta*c3(k,j,i)+alpha*c5(k,j,i)
        end do
      end do
    end do
#if defined(GHOST_GPU)
!$omp target teams distribute parallel do collapse(3) if(target: gdev_active)
#else
!$omp parallel do collapse(2) private (k)
#endif
    do i = ista,iend
      do j = 1,ny
        do k = 1,nz
          if (kn2(k,j,i).le.kmax) then
            dre(k,j,i) = c4(k,j,i)+cr*c7(k,j,i)
            dim(k,j,i) = c3(k,j,i)+cr*c8(k,j,i)
          else
            dre(k,j,i) = 0.0_GP
            dim(k,j,i) = 0.0_GP
          endif
        end do
      end do
    end do
  end subroutine gpe_rhs_kernel


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

    class (GPESolver), intent(in)               :: this
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
  !! Function to compute and write spectra (and the transfer
  !! functions if dotrans)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine spectra_impl(this, uin)
    use pseudospec_gpe
    use pseudospec_rgpe
    use filefmt
    use status
    implicit none

    class (GPESolver), intent(in)               :: this
    type(GStateComp), intent(in), target        :: uin(:)
    complex(kind=GP), pointer, dimension(:,:,:) :: zre,zim

    write(ext, fmtext) sind
    zre => uin(this%ZFUNC  )%ccomp
    zim => uin(this%ZFUNC+1)%ccomp
    call gpemassspec(zre,zim,this%todir_,ext)
    call gperealspec(zre,zim,this%alpha_,this%beta_,this%omegag_,this%todir_,ext)
    call gpehelspec (zre,zim,this%alpha_,this%beta_,this%omegag_,this%todir_,ext)
    call gpemomtspec(zre,zim,this%alpha_,this%todir_,ext)
    if ( this%dorot_ ) then
      call gperealspecperp(zre,zim,this%alpha_,this%beta_,this%omegag_,this%todir_,ext)
    endif
    if ( this%traits_%dotrans ) call gpe_transfer(this,zre,zim,ext)
  end subroutine spectra_impl


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Energy transfer functions: the spectra of the state
  !! advanced one time step (with a second order Runge-Kutta
  !! step) are compared with the spectra of the current state,
  !! T(k) = dE(k)/dt (as in the old code, where the spectra
  !! before and after a time step were used).
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine gpe_transfer(this, zre, zim, nmb)
    use pseudospec_fluid
    use pseudospec_gpe
    use grid
    use kes
    use mpivars
    use status
    implicit none

    class (GPESolver), intent(in)                              :: this
    complex(kind=GP) , intent(in), dimension(nz,ny,ista:iend)  :: zre,zim
    character(len=*) , intent(in)                              :: nmb
    complex(kind=GP) , pointer, dimension(:,:,:)               :: Z1,Z2,K1,K2
    double precision , dimension(nmax/2+1)                     :: io,qo,ko,co
    double precision , dimension(nmax/2+1)                     :: in,qn,kn,cn
    logical                                                    :: bret

    call this%workspace_%get_complex_tmp(Z1,bret)
    call this%workspace_%get_complex_tmp(Z2,bret)
    call this%workspace_%get_complex_tmp(K1,bret)
    call this%workspace_%get_complex_tmp(K2,bret)
    call gperealspecc(zre,zim,this%alpha_,this%beta_,this%omegag_,io,qo,ko,co)
    call gpe_rhs(this,zre,zim,K1,K2)                ! K1 = f(z)
    call saxpby_c(Z1,zre,1.0_GP,K1,0.5_GP*dt)       ! z + dt/2 f(z)
    call saxpby_c(Z2,zim,1.0_GP,K2,0.5_GP*dt)
    call gpe_rhs(this,Z1,Z2,K1,K2)                  ! K2 = f(z + dt/2 f(z))
    call saxpby_c(Z1,zre,1.0_GP,K1,dt)              ! z(t+dt)
    call saxpby_c(Z2,zim,1.0_GP,K2,dt)
    call gperealspecc(Z1,Z2,this%alpha_,this%beta_,this%omegag_,in,qn,kn,cn)
    call gperealtrans(dt,io,qo,ko,co,in,qn,kn,cn,this%todir_,nmb)
    call this%workspace_%free_complex_tmp(K2)
    call this%workspace_%free_complex_tmp(K1)
    call this%workspace_%free_complex_tmp(Z2)
    call this%workspace_%free_complex_tmp(Z1)
  end subroutine gpe_transfer


  ! ===================================================================
  ! Solver specific methods
  ! ===================================================================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Auxiliary arrays to the device (the potential and the
  !! linear ramps for the rotation). Called by the main program
  !! after the initial conditions and the forcing are set (on
  !! the host).
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine sync_device_impl(this)
    implicit none
    class (GPESolver), intent(inout) :: this
    call this%aux_to_device()
  end subroutine sync_device_impl


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Constructor: reads the parameters and allocates the
  !! auxiliary arrays (zero by default)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine GPESolver_ctor(this, infile, workspace, plan)
    use iovar
    class  (GPESolver), intent(inout)        :: this
    type(GWorkspace) , intent(inout), target :: workspace
    type(ioplan)     , intent(inout), target :: plan
    character(len=*) , intent   (in)         :: infile
    this%infile_    =  infile    ! input file
    this%workspace_ => workspace
    this%planio_    => plan
    call this%init()
    call this%alloc_aux()
  end subroutine GPESolver_ctor


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Destructor
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine GPESolver_dtor(this)
    type  (GPESolver), intent(inout) :: this
    if (associated(this%workspace_))   nullify(this%workspace_)
    if (associated(this%planio_))      nullify(this%planio_)
    if (allocated(this%sstate_))       deallocate(this%sstate_)
    call this%free_aux()
  end subroutine GPESolver_dtor


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Convert input state name to index in state vector
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine sstate2istate_impl(this, sstate, istate)
    class (GPESolver), intent   (in) :: this
    character(len=8), intent   (in) :: sstate(:)
    integer         , intent(inout) :: istate(:)
    integer                         :: i,j
    if ( size(sstate) .ne. size(istate) ) then
      stop 'GPESolver::sstate2istate_impl: Incompatible sstate and istate'
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
    class  (GPESolver), intent   (in) :: this
    character (len=8), intent(inout) :: sstate(:)
    sstate(this%ZFUNC  ) = 'zre'
    sstate(this%ZFUNC+1) = 'zim'
  end subroutine get_sstate_impl


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Function to compute number of state members (equations)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PURE function state_size_impl(this) result(num)
    class(GPESolver), intent(in) :: this
    integer                      :: num
    num = 2                      ! real and imaginary parts of z
  end function state_size_impl

end module gpe_mod
