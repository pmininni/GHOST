! =====================================================================
! NAME       : bouss_solver.f90
! DESCRIPTION: Forms class for incompressible Boussinesq solvers, computing:
!
!              dv/dt + (w + 2 Omega) x v = - Grad p - N theta \hat{z} + nu Del^2 v
!              dtheta/dt + v.Grad theta  = - N v_z + kappa  Del^2 theta
!
!              ds_i/dt + v.Grad s_i = kappa_i Del^2 s_i 
!                                             i = 1, ..., numpassive
!              Div v = 0,
!
!              where N is Brunt-Vaisala frequency, nu, kappa are
!              viscostiy and diffusivity.
!
!
!              State ordering is:
!                v1, v2, v3, theta, s1, s2, ..., s_numpassive
!
!              State sector ids are:
!                VELOCITY (VELOCITY+1, VELOCITY+2)   : momentum sector
!                TEMP                                : temperature sector
!                PASSIVE  (PASSIVE+1, PASSIVE+2, ...): passive scalar sector
!
! INPUT FILE : For solver='BOUSS', looks for a "&BOUSS" namelist with:
!              nu      : fluid kinematic viscosity
!              bkappa  : active scalar diffusivity
!              bvfreq  : Brunt-Vaisala frequency
!              xmom    : factor multiplying buoyancy tendency 
!                        in momentum eqn
!              xtemp   : factor multiplying buoyancy tendency i
!                        in temperature eqn
!              dorot   : do rotation, = .TRUE. or .FALSE.
!              doparts : use particles = .TRUE. or .FALSE.
!              omegax  : amplitude of the uniform rotation along x
!              omegay  : amplitude of the uniform rotation along y
!              omegaz  : amplitude of the uniform rotation along z
!              npassive: number of passive scalars (default=0)
!              spectlod: spectral output level of detail (in [1,4]):
!                          1: All 1d spectra, KE and PE fluxes
!                          2: 2D spectra, helicity flluxes
!                          3: KE Fourier modes
!                          4: PV spectra, horizontally-avged quantities

!              For npassive > 0, looks for a "&passive" namelist with:
!              kappa   : vector with npassive diffusivities
!
! DATE       : 2/2/26 (JBG)
! =====================================================================

module moist_mod
  USE equationbase_mod
  USE gstate_mod

  implicit none

  ! ================= Solver traits ===================================
  type, public  :: BSTraits
    logical       :: dorot        = .FALSE. ! rotation flag
    integer       :: spectlod     = 1       ! standard level of spectra detail 
    real(kind=GP) :: nu           = 0.0_GP  ! dissipation
    real(kind=GP) :: bvfreq       = 0.0_GP  ! 
    real(kind=GP) :: cbu         = 0.0_GP  ! 
    real(kind=GP) :: cbs         = 0.0_GP  ! 
    real(kind=GP) :: kappau         = 0.0_GP  ! 
    real(kind=GP) :: kappas         = 0.0_GP  ! 
    real(kind=GP), allocatable :: kappa(:)  ! diffusivities
    real(kind=GP)              :: omega(3)  ! rotation vector
  end type

  ! ================= Global parameters ===============================
  integer, parameter, public   :: MAXPASSIVE = 20 ! max # passive scalars

  ! ================= Solver ==========================================
  ! Define class:
  type, extends(VelocityBase) :: MOISTSolver 
    ! Member data:
    logical           :: binit_=.false. ! is initialized?
    type  (BSTraits)  :: traits_
  contains
    procedure, public :: init          =>          init_impl ! init method
    procedure, public :: dudt          =>          dudt_impl ! RHS method
    procedure, public :: pdudt         =>         pdudt_impl ! part+field RHS method
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
    class  (MOISTSolver), intent (inout) :: this

    ! Temporary data to read from namelists:
    logical                    :: dorot
    ! logical                    :: doparts
    integer                    :: npassive
    integer                    :: spectlod
    integer                    :: ierr
    real(kind=GP)              :: nu, bvfreq, cbu, cbs    
    real(kind=GP)              :: kappau, kappas
    real(kind=GP)              :: omegax, omegay, omegaz
    real(kind=GP), allocatable :: kappa(:)

    ! Required namelists:
    namelist /MOIST/ nu, bvfreq, cbu, cbs, kappau, kappas, &
                    dorot, omegax, omegay, omegaz, npassive, spectlod
    namelist /passive/ kappa

    call MPI_COMM_SIZE(MPI_COMM_WORLD,this%nprocs_,ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,this%myrank_,ierr)

    ! Get trait variables from input file:
    dorot    = .FALSE.
    ! doparts  = .FALSE.
    spectlod = 1 ! standard lod
    nu       = 0.0_GP
    bvfreq   = 0.0_GP
    cbu      = 0.0_GP
    cbs      = 0.0_GP
    kappau   = 0.0_GP
    kappas   = 0.0_GP
    npassive = 0
    omegax   = 0.0_GP; omegay = 0.0_GP; omegaz = 0.0_GP

    if ( this%myrank_ .eq. 0 ) then
      open(1,file=this%infile_,status='unknown',form="formatted")
      read(1,NML=MOIST)
      close(1)
    endif

    call mpi_bcast(nu      ,1,GC_REAL    ,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(bvfreq  ,1,GC_REAL    ,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(cbu     ,1,GC_REAL    ,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(cbs     ,1,GC_REAL    ,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(kappau  ,1,GC_REAL    ,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(kappas  ,1,GC_REAL    ,0,MPI_COMM_WORLD,ierr)

    call mpi_bcast(dorot   ,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)

    call mpi_bcast(omegax  ,1,GC_REAL    ,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(omegay  ,1,GC_REAL    ,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(omegaz  ,1,GC_REAL    ,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(npassive,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(spectlod,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

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

    ! Set traits structure from inputfile data:
    this%traits_%  dorot = dorot
    this%traits_%     nu = nu
    this%traits_% bvfreq = bvfreq 
    this%traits_%    cbu = cbu 
    this%traits_%    cbs = cbs 
    this%traits_% kappau = kappau 
    this%traits_% kappas = kappas
    this%traits_%  omega = (/omegax,omegay,omegaz/)

    if ( npassive .gt. 0 ) then
      if ( allocated(this%traits_%kappa) ) then
        deallocate(this%traits_%kappa);
      endif
      allocate(this%traits_%kappa(npassive))
      this%traits_%kappa = kappa
      deallocate(kappa)
    endif

    this%traits_%spectlod =  spectlod

    this%order_   = 2                        ! Time stepping order
    this%nd_      = 3                        ! 3d
    this%nc_      = this%nd_                 ! #vel field components
    this%VELOCITY = 1                        ! start of vel sector
    this%BUOYU    = this%VELOCITY+this%nc_   ! start of buoyancies sector
    this%BUOYS    = this%BUOYU + 1           ! start of buoyancies sector
    this%PASSIVE  = this%BUOYS + 1         ! start of scalar sector

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
    use pseudospec_fluid
    use ali
    use kes
    use var
    use grid
    use mpivars
  !$ use threads
    implicit none
  
    class (MOISTSolver), intent(in)            :: this
    real   (kind=GP), intent(in)               :: time, dt
    type(GStateComp), intent(inout), target    :: uin(:), uf(:)
    type(GStateComp), intent(inout)            :: dudt(:)
  
    complex(kind=GP), pointer, dimension(:,:,:) :: fx, fy, fz, vx, vy, vz
    complex(kind=GP), pointer, dimension(:,:,:) :: fbu, fbs, bu, bs
    complex(kind=GP), pointer, dimension(:,:,:) :: C1, C2, C3, C4, C5, C6
    complex(kind=GP), pointer, dimension(:,:,:) :: C7, C8
  
    real(kind=GP) :: nu
    real(kind=GP) :: bvfreq, cbu, cbs
    real(kind=GP) :: kappau, kappas
    real(kind=GP) :: nufreq, nsfreq
    real(kind=GP) :: omegax, omegay, omegaz
    integer       :: i, j, k
    logical       :: bret
  
    if (.not. this%binit_) then
      stop 'MOISTSolver::dudt: Solver not initialized'
    endif
  
    nu     = this%traits_%nu
    bvfreq = this%traits_%bvfreq
    cbu    = this%traits_%cbu
    cbs    = this%traits_%cbs
    kappau = this%traits_%kappau
    kappas = this%traits_%kappas
  
    nufreq = bvfreq * cbu
    nsfreq = bvfreq * cbs
  
    call this%workspace_%get_complex_tmp(C1,bret)
    call this%workspace_%get_complex_tmp(C2,bret)
    call this%workspace_%get_complex_tmp(C3,bret)
    call this%workspace_%get_complex_tmp(C4,bret)
    call this%workspace_%get_complex_tmp(C5,bret)
    call this%workspace_%get_complex_tmp(C6,bret)
    call this%workspace_%get_complex_tmp(C7,bret)
    call this%workspace_%get_complex_tmp(C8,bret)
  
    vx  => uin(this%VELOCITY  )%ccomp
    vy  => uin(this%VELOCITY+1)%ccomp
    vz  => uin(this%VELOCITY+2)%ccomp
    bu  => uin(this%BUOYU     )%ccomp
    bs  => uin(this%BUOYS     )%ccomp
  
    fx  => uf (this%VELOCITY  )%ccomp
    fy  => uf (this%VELOCITY+1)%ccomp
    fz  => uf (this%VELOCITY+2)%ccomp
    fbu => uf (this%BUOYU     )%ccomp
    fbs => uf (this%BUOYS     )%ccomp
  
    call prodre3(vx,vy,vz,C4,C5,C6)     ! w x v
  
    if (this%traits_%dorot) then
      omegax = this%traits_%omega(1)
      omegay = this%traits_%omega(2)
      omegaz = this%traits_%omega(3)
  
      call saxpby_c(C1, vz, 2*omegay, vy, -2.0_GP*omegaz)
      call saxpby_c(C2, vx, 2*omegaz, vz, -2.0_GP*omegax)
      call saxpby_c(C3, vy, 2*omegax, vx, -2.0_GP*omegay)
  
  !$omp parallel do if (iend-ista.ge.nth) private(j,k)
      do i = ista, iend
  !$omp parallel do if (iend-ista.lt.nth) private(k)
        do j = 1, ny
          do k = 1, nz
            C4(k,j,i) = C4(k,j,i) + C1(k,j,i)
            C5(k,j,i) = C5(k,j,i) + C2(k,j,i)
            C6(k,j,i) = C6(k,j,i) + C3(k,j,i)
          enddo
        enddo
      enddo
    endif
  
    C7 = bu
    call fftp3d_complex_to_real(plancr, C7, R1, MPI_COMM_WORLD)
  
    C8 = bs
    call fftp3d_complex_to_real(plancr, C8, R2, MPI_COMM_WORLD)
  
    tmp = 1.0_GP/(real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))
  
  !$omp parallel do if (kend-ksta.ge.nth) private(j,i)
    do k = ksta, kend
  !$omp parallel do if (kend-ksta.lt.nth) private(i)
      do j = 1, ny
        do i = 1, nx
          if (R1(i,j,k) .gt. R2(i,j,k)) then
            R1(i,j,k) = tmp*xmom*R1(i,j,k)
          else
            R1(i,j,k) = tmp*xmom*R2(i,j,k)
          endif
        enddo
      enddo
    enddo
  
    call fftp3d_real_to_complex(planrc, R1, C7, MPI_COMM_WORLD)
  
    DO i = ista,iend               ! Buoyancy force
  !$omp parallel do if (iend-ista.lt.nth) private(k)
       DO j = 1,ny
          DO k = 1,nz
             C6(k,j,i) = C6(k,j,i)+C7(k,j,i)
          END DO
       END DO
    END DO
    CALL nonlhd3(C4,C5,C6,C1,1)     ! -[(w + 2 Omega) x v + Grad p]_x
    CALL nonlhd3(C4,C5,C6,C2,2)     ! -[(w + 2 Omega) x v + Grad p]_y
    CALL nonlhd3(C4,C5,C6,C3,3)     ! -[(w + 2 Omega) x v + Grad p]_z
    CALL advect3(vx,vy,vz,bu,C7)
    CALL advect3(vx,vy,vz,bs,C8)
  
  !$omp parallel do if (iend-ista.ge.nth) private(j,k)
    do i = ista, iend
  !$omp parallel do if (iend-ista.lt.nth) private(k)
      do j = 1, ny
        do k = 1, nz
          C7(k,j,i) = C7(k,j,i) + xtemp*vz(k,j,i)*nufreq
          C8(k,j,i) = C8(k,j,i) + xtemp*vz(k,j,i)*nsfreq
        enddo
      enddo
    enddo
  
    call laplak3(vx, vx)
    call laplak3(vy, vy)
    call laplak3(vz, vz)
    call laplak3(bu, bu)
    call laplak3(bs, bs)
  
  !$omp parallel do if (iend-ista.ge.nth) private(j,k)
    do i = ista, iend
  !$omp parallel do if (iend-ista.lt.nth) private(k)
      do j = 1, ny
        do k = 1, nz
          if ((kn2(k,j,i) .le. kmax) .and. (kn2(k,j,i) .ge. tiny)) then
            dudt(this%VELOCITY  )%ccomp(k,j,i) = nu    *vx(k,j,i) + C1(k,j,i) + fx(k,j,i)
            dudt(this%VELOCITY+1)%ccomp(k,j,i) = nu    *vy(k,j,i) + C2(k,j,i) + fy(k,j,i)
            dudt(this%VELOCITY+2)%ccomp(k,j,i) = nu    *vz(k,j,i) + C3(k,j,i) + fz(k,j,i)
            dudt(this%BUOYU     )%ccomp(k,j,i) = kappau*bu(k,j,i) + C7(k,j,i) + fbu(k,j,i)
            dudt(this%BUOYS     )%ccomp(k,j,i) = kappas*bs(k,j,i) + C8(k,j,i) + fbs(k,j,i)
          else
            dudt(this%VELOCITY  )%ccomp(k,j,i) = 0.0_GP
            dudt(this%VELOCITY+1)%ccomp(k,j,i) = 0.0_GP
            dudt(this%VELOCITY+2)%ccomp(k,j,i) = 0.0_GP
            dudt(this%BUOYU     )%ccomp(k,j,i) = 0.0_GP
            dudt(this%BUOYS     )%ccomp(k,j,i) = 0.0_GP
          endif
        enddo
      enddo
    enddo
  
    call this%workspace_%free_complex_tmp(C1)
    call this%workspace_%free_complex_tmp(C2)
    call this%workspace_%free_complex_tmp(C3)
    call this%workspace_%free_complex_tmp(C4)
    call this%workspace_%free_complex_tmp(C5)
    call this%workspace_%free_complex_tmp(C6)
    call this%workspace_%free_complex_tmp(C7)
    call this%workspace_%free_complex_tmp(C8)
  
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

    class (MOISTSolver), intent(in)             :: this
    type(GStateComp), intent(in), target        :: uin(:), uf(:)
    integer         , intent(in)                :: t
    complex(kind=GP), pointer, dimension(:,:,:) :: fbu,fbs,fx,fy,fz,vx,vy,vz,bu,bs
    complex(kind=GP), pointer, dimension(:,:,:) :: c1,c2,c3
    real   (kind=GP)                            :: rmp,rmq_bu,rmq_bs
    integer                                     :: i
    logical                                     :: bret

    call this%workspace_%get_complex_tmp(C1,bret)
    call this%workspace_%get_complex_tmp(C2,bret)
    call this%workspace_%get_complex_tmp(C3,bret)

    vx  => uin(this%VELOCITY  )%ccomp
    vy  => uin(this%VELOCITY+1)%ccomp
    vz  => uin(this%VELOCITY+2)%ccomp
    bu  => uin(this%BUOYU     )%ccomp
    bs  => uin(this%BUOYS     )%ccomp
    fx  => uf (this%VELOCITY  )%ccomp
    fy  => uf (this%VELOCITY+1)%ccomp
    fz  => uf (this%VELOCITY+2)%ccomp
    fbu => uf (this%BUOYU     )%ccomp
    fbs => uf (this%BUOYS     )%ccomp

    CALL hdcheck(vx,vy,vz,fx,fy,fz,t,dt,1,1)
    ! CALL tbouss(vx,vy,vz,bu,t,dt,0.0_GP,bvfreq)  TODO subrutine tbouss not defined
    CALL maxabs(vx,vy,vz,rmp,0)
    ! CALL mpscheck2(th1,fs1,th2,fs2,t,dt)     TODO subrutine mpsckeck2 not defined

    call pscheck(bu,fbu,t,dt)
    call derivk3(bu,c1,1)
    call derivk3(bu,c2,2)
    call derivk3(bu,c3,3)
    call maxabs(c1,c2,c3,rmq_bu,2)

    call pscheck(bs,fbs,t,dt)
    call derivk3(bs,c1,1)
    call derivk3(bs,c2,2)
    call derivk3(bs,c3,3)
    call maxabs(c1,c2,c3,rmq_bs,2)

    IF (myrank.eq.0) THEN
      open(1,file='maximum.txt',position='append')
      write(1,FMT='(E13.6,E13.6,E13.6,E13.6)') (t-1)*dt,rmp,rmq_bu,rmq_bs
      close(1)
    endif

    do i = this%PASSIVE, this%PASSIVE+this%numpassive_-1
      call pscheck(uin(i)%ccomp,uf(i)%ccomp,t,dt,trim(this%sstate_(i)))
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
    use filefmt
    use status
    use iovar
    implicit none

    class (MOISTSolver), intent(in)             :: this
    type(GStateComp), intent(in), target        :: uin(:)
    complex(kind=GP), pointer, dimension(:,:,:) :: vx,vy,vz,th
    complex(kind=GP), pointer, dimension(:,:,:) :: c1,c2,c3
    real   (kind=GP)                            :: bvfreq,omegaz
    integer                                     :: i
    logical                                     :: bret

    bvfreq = this%traits_%bvfreq 
    omegaz = this%traits_% omega(3) 

!   write(ext, fmtext) sind
    vx => uin(this%VELOCITY  )%ccomp
    vy => uin(this%VELOCITY+1)%ccomp
    vz => uin(this%VELOCITY+2)%ccomp
    th => uin(this%TEMP)      %ccomp
    call spectrum(vx,vy,vz,ext,1,1)
    call this%workspace_%get_complex_tmp(c1,bret)
    call this%workspace_%get_complex_tmp(c2,bret)
    call this%workspace_%get_complex_tmp(c3,bret)
    call gradre3(vx,vy,vz,c1,c2,c3)            ! Computes v.Grad(v)
    call entrans(vx,vy,vz,-c1,-c2,-c3,ext,1)   ! Writes the energy flux

    call specpara(vx,vy,vz,ext,1,1)
    call specperp(vx,vy,vz,ext,1,1)
    call spectrsc(th,ext,0)
    call specscpa(th,ext,0)
    call specscpe(th,ext,0)

    call entpara(vx,vy,vz,-c1,-c2,-c3,ext,1) ! Writes the energy flux
    call entperp(vx,vy,vz,-c1,-c2,-c3,ext,1) ! Writes the energy fluxq

    ! Write helicity fluxes, 2D spectra:
    if ( this%traits_%spectlod .ge. 2 ) then
      call heltrans(vx,vy,vz,-c1,-c2,-c3,ext,1)
      call helpara(vx,vy,vz,-c1,-c2,-c3,ext,1)
      call helperp(vx,vy,vz,-c1,-c2,-c3,ext,1)
      ! Write 2D spectra:
      call spec2D(vx,vy,vz,ext,odir,1,1)
      call specsc2D(th,ext,odir,0)
    endif

    ! Write Fourier modes:
    if ( this%traits_%spectlod .ge. 3 ) then
      call write_fourier(vx,'vx',ext,odir)
      call write_fourier(vy,'vy',ext,odir)
      call write_fourier(vz,'vz',ext,odir)
      call write_fourier(th,'th',ext,odir)
    endif

    ! Write PV spectra, horizontally-averaged data:
    if ( this%traits_%spectlod .ge. 4 ) then
      call spectpv(vx,vy,vz,th,ext)
      call havgwrite(0,'shear'  ,ext,vx,vy,vz,th,omegaz,bvfreq) ! shear
      call havgwrite(1,'tgradz' ,ext,vx,vy,vz,th,omegaz,bvfreq) ! dtheta/dz
      call havgwrite(2,'hawdtdz',ext,vx,vy,vz,th,omegaz,bvfreq) ! u_z*dtheta/dz
      call havgwrite(3,'hahke'  ,ext,vx,vy,vz,th,omegaz,bvfreq) ! hor. k.e.
      call havgwrite(4,'havke'  ,ext,vx,vy,vz,th,omegaz,bvfreq) ! vert. k.e.
      call havgwrite(5,'haphel' ,ext,vx,vy,vz,th,omegaz,bvfreq) ! perp. helicity
      call havgwrite(6,'haomzt' ,ext,vx,vy,vz,th,omegaz,bvfreq) ! ometa_z*theta
      call havgwrite(7,'hapv2'  ,ext,vx,vy,vz,th,omegaz,bvfreq) ! pot'l vorticity^2
      call havgwrite(8,'hasuph' ,ext,vx,vy,vz,th,omegaz,bvfreq) ! super-helicity
      call havgwrite(9,'hari'   ,ext,vx,vy,vz,th,omegaz,bvfreq) ! Richardson no.

    endif


    call this%workspace_%free_complex_tmp(c1)
    call this%workspace_%free_complex_tmp(c2)
    call this%workspace_%free_complex_tmp(c3)
    if ( this%numpassive_ .gt. 0) then
      do i = this%PASSIVE, this%PASSIVE+this%numpassive_-1
        call spectrsc(uin(i)%ccomp,ext,0,trim(this%sstate_(i)))
        if ( this%traits_%dorot ) then
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
  subroutine MOISTSolver_ctor(this, infile, workspace, plan)
    use iovar
    class  (MOISTSolver), intent(inout)      :: this
    type(GWorkspace) , intent(inout), target :: workspace
    type(ioplan)     , intent(inout), target :: plan
    character(len=*) , intent   (in)         :: infile
    this%infile_    =  infile    ! input file
    this%workspace_ => workspace
    this%planio_    => plan
    call this%init()
  end subroutine BOUSSSolver_ctor


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Destructor
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine MOISTSolver_dtor(this) 
    type  (MOISTSolver), intent(inout) :: this
    if (associated(this%workspace_))   nullify(this%workspace_)
    if (associated(this%planio_))      nullify(this%planio_)
    if (allocated(this%sstate_))       deallocate(this%sstate_)
    if (allocated(this%traits_%kappa)) deallocate(this%traits_%kappa)
  end subroutine BOUSSSolver_dtor


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Convert input state name to index in state vector
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine sstate2istate_impl(this, sstate, istate) 
    class (BOUSSSolver)              , intent   (in) :: this
    character(len=8)              , intent   (in) :: sstate(:)
    integer         , allocatable , intent(inout) :: istate(:)
    integer                                       :: i,j
    if ( size(sstate) .ne. size(istate) ) then
      stop 'BOUSSSolver::sstate2istate_impl: Incompatible sstate and istate'
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
    class  (BOUSSSolver)             , intent   (in) :: this
    character (len=8), allocatable, intent(inout) :: sstate(:)
    character(len=100)                            :: snum
    character(len=1)                              :: comp(3)
    integer                                       :: j

    if ( size(sstate) .lt. this%state_size() ) then
      deallocate(sstate)
      allocate(sstate,this%state_size())
    endif

    comp = ['x', 'y', 'z']
    do j = this%VELOCITY,this%VELOCITY+this%nc_-1
       sstate(j) = 'v' // comp(j-this%VELOCITY+1)
    enddo
    sstate(this%TEMP) = 'th' 
    do j = this%PASSIVE,this%PASSIVE+this%numpassive_-1
       write(snum,'(I0)') j-this%PASSIVE+1
       sstate(j) = 's' // trim(snum)
    enddo
  end subroutine get_sstate_impl
  

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Function to compute number of state members (equations)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PURE function state_size_impl(this) result(num)
    class(BOUSSSolver), intent(in) :: this
    integer                     :: num
    num = this%nc_ + 1           ! # vel. components + theta
    num = num + this%numpassive_ ! # passive scalars
  end function state_size_impl

end module bouss_mod
