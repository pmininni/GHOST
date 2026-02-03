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
!              omegax  : amplitude of the uniform rotation along x
!              omegay  : amplitude of the uniform rotation along y
!              omegaz  : amplitude of the uniform rotation along z
!              npassive: number of passive scalars (default=0)

!              For npassive > 0, looks for a "&passive" namelist with:
!              kappa   : vector with npassive diffusivities
!
! DATE       : 2/2/26 (DLR)
! =====================================================================

module bouss_mod
  USE equationbase_mod
  USE gstate_mod

  IMPLICIT NONE

  ! ================= Solver traits ===================================
  type, public  :: NHTraits
    logical       :: dorot        = .FALSE. ! rotation flag
    real(kind=GP) :: nu           = 0.0_GP  ! dissipation
    real(kind=GP) :: bkappa       = 0.0_GP  ! diffusivity
    real(kind=GP) :: bvfreq       = 0.0_GP  ! Brunt-Vaisala freq
    real(kind=GP), allocatable :: kappa(:)  ! diffusivities
    real(kind=GP)              :: omega(3)  ! rotation vector
  end type

  ! ================= Global parameters ===============================
  integer, parameter, public   :: MAXPASSIVE = 20 ! max # passive scalars

  ! ================= Solver ==========================================
  ! Define class:
  type, extends(VelocityBase) :: BOUSSSolver 
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
    procedure, public :: Solver_ctor   =>      BOUSSSolver_ctor ! constructor
    final             :: BOUSSSolver_dtor
  end type BOUSSSolver

CONTAINS

  ! ===================================================================
  ! Solver initialization, this is where parameter files are read
  ! ===================================================================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Subroutine to initialize the solver
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_impl(this)
    USE commtypes
    class  (BOUSSSolver), intent (inout) :: this

    ! Temporary data to read from namelists:
    logical                    :: dorot
    integer                    :: npassive
    integer                    :: ierr
    real(kind=GP)              :: nu, bkappa, bvfreq, omegax, omegay, omegaz
    real(kind=GP), allocatable :: kappa(:)

    ! Required namelists:
    namelist/ BOUSS      / nu, bkappa, bvfreq, xmom, xtamp, &
                           dorot, omegax, omegay, omegaz, npassive
    namelist/ passive    / kappa

    call MPI_COMM_SIZE(MPI_COMM_WORLD,this%nprocs_,ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,this%myrank_,ierr)

    ! Get trait variables from input file:
    dorot    = .FALSE.
    nu       = 0.0
    bkappa   = 0.0
    bvfreq   = 0.0
    xmom     = 1.0
    xtemp    = 1.0
    npassive = 0
    omegax   = 0.0_GP; omegay = 0.0_GP; omegaz = 0.0_GP
    if ( this%myrank_ .eq. 0 ) then
      open(1,file=this%infile_,status='unknown',form="formatted")
      read(1,NML=HD)
      close(1)
    endif
    call mpi_bcast(nu       ,1 ,GC_REAL,    0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(bkappa   ,1 ,GC_REAL,    0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(bvfreq   ,1 ,GC_REAL,    0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(xmom     ,1 ,GC_REAL,    0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(xtemp    ,1 ,GC_REAL,    0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(dorot    ,1 ,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(omegax   ,1 ,GC_REAL,    0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(omegay   ,1 ,GC_REAL,    0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(omegaz   ,1 ,GC_REAL,    0,MPI_COMM_WORLD,ierr)
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

    ! Set traits structure from inputfile data:
    this%traits_% dorot = dorot
    this%traits_%    nu = nu
    this%traits_%bkappa = bkappa
    this%traits_%bvfreq = bvfreq
    this%traits_% omega = (/omegax,omegay,omegaz/)
    if ( npassive .gt. 0 ) then
      if ( allocated(this%traits_%kappa) ) then
        deallocate(this%traits_%kappa);
      endif
      allocate(this%traits_%kappa(npassive))
      this%traits_%kappa = kappa
      deallocate(kappa)
    endif

    this%order_   = 2                        ! Time stepping order
    this%nd_      = 3                        ! 3d
    this%nc_      = this%nd_                 ! #vel field components
    this%VELOCITY = 1                        ! start of vel sector
    this%TEMP     = this%VELOCITY+this%nc_   ! start of temperature sector
    this%PASSIVE  = this%TEMP + 1            ! start of scalar sector

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
!$  use threads
    implicit none

    class (BOUSSSolver), intent   (in)             :: this
    real   (kind=GP), intent   (in)             :: time, dt
    type    (GState), intent(inout), target     :: uin(:),uf(:)
    type    (GState), intent(inout)             :: dudt(:) 
    complex(kind=GP), pointer, dimension(:,:,:) :: fx,fy,fz,vx,vy,vz
    complex(kind=GP), pointer, dimension(:,:,:) :: fth,th
    complex(kind=GP), pointer, dimension(:,:,:) :: C1,C2,C3,C4,C5,C6
    complex(kind=GP), pointer, dimension(:,:,:) :: C7,C8
    real   (kind=GP)                            :: bkappa,bvfreq,nu
    real   (kind=GP)                            :: xmom,xtemp
    real   (kind=GP)                            :: omegax,omegay,omegaz
    integer                                     :: i,j,k
    logical                                     :: bret
       
    nu     = this%traits_%nu
    bkappa = this%traits_%bkappa
    bvfreq = this%traits_%bvfreq
    xmom   = this%traits_%xmom  * this%bvfreq
    xtemp  = this%traits_%xtemp * this%traits_%bvfreq

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
      
    call prodre3(vx,vy,vz,C4,C5,C6)                    ! w x v
    if ( this%traits_%dorot ) then
      omegax = this%traits_%omega(1)
      omegay = this%traits_%omega(2)
      omegaz = this%traits_%omega(3)
      call saxpby_c(C1, vz, 2*omegay, vy, -2.0*omegaz) ! 2 Omega x v
      call saxpby_c(C2, vx, 2*omegaz, vz, -2.0*omegax)
      call saxpby_c(C3, vy, 2*omegax, vx, -2.0*omegay)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      do i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
      do j = 1,ny
      do k = 1,nz
         C4(k,j,i) = C4(k,j,i) + C1(k,j,i) ! (w x v + 2 Omega x v)_x
         C5(k,j,i) = C5(k,j,i) + C2(k,j,i) ! (w x v + 2 Omega x v)_y
         C6(k,j,i) = C6(k,j,i) + C3(k,j,i) ! (w x v + 2 Omega x v)_z
      end do
      end do
      end do
    endif

!$omp parallel do if (iend-ista.ge.nth) private (j,k)
    do i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
    do j = 1,ny
    do k = 1,nz
       C6(k,j,i) = C6(k,j,i) - xmom*th(k,j,i) ! buoyancy term
    end do
    end do
    end do

    call nonlhd3(C4,C5,C6,C1,1)  ! -[(w + 2 Omega) x v + Grad p]_x
    call nonlhd3(C4,C5,C6,C2,2)  ! -[(w + 2 Omega) x v + Grad p]_y
    call nonlhd3(C4,C5,C6,C3,3)  ! -[(w + 2 Omega) x v + Grad p]_z
    call laplak3(vx,C4)          ! Del^2 vx
    call laplak3(vy,C5)          ! Del^2 vy
    call laplak3(vz,C6)          ! Del^2 vz

    CALL advect3(vx,vy,vz,th,C7) ! -(v.Grad) th
    call laplak3(th,C8)          ! Del^2 th

!$omp parallel do if (iend-ista.ge.nth) private (j,k)
    DO i = ista,iend               ! heat 'currrent':
!$omp parallel do if (iend-ista.lt.nth) private (k)
      DO j = 1,ny
        DO k = 1,nz
          C7(k,j,i) = C7(k,j,i) + xtemp*vz(k,j,i) ! add N vz term
        END DO
      END DO
    END DO


!$omp parallel do if (iend-ista.ge.nth) private (j,k)
    do i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
    do j = 1,ny
    do k = 1,nz
      if ((kn2(k,j,i).le.kmax).and.(kn2(k,j,i).ge.tiny)) then
        dudt(this%VELOCITY  )%ccomp(k,j,i) = nu*C4(k,j,i) + C1(k,j,i) + fx(k,j,i)
        dudt(this%VELOCITY+1)%ccomp(k,j,i) = nu*C5(k,j,i) + C2(k,j,i) + fy(k,j,i)
        dudt(this%VELOCITY+2)%ccomp(k,j,i) = nu*C6(k,j,i) + C3(k,j,i) + fz(k,j,i)
        dudt(this%TEMP)%ccomp(k,j,i) = bkappa*C8(k,j,i) + C7(k,j,i) + fth(k,j,i)
      else
        dudt(this%VELOCITY  )%ccomp(k,j,i) = 0.0_GP
        dudt(this%VELOCITY+1)%ccomp(k,j,i) = 0.0_GP
        dudt(this%VELOCITY+2)%ccomp(k,j,i) = 0.0_GP
        dudt(this%TEMP)      %ccomp(k,j,i) = 0.0_GP
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

    class (BOUSSSolver), intent(in)                :: this
    type    (GState), intent(in), target        :: uin(:), uf(:)
    integer         , intent(in)                :: t
    complex(kind=GP), pointer, dimension(:,:,:) :: fx,fy,fz,vx,vy,vz,th
    real   (kind=GP)                            :: rmp
    integer                                     :: i

    vx => uin(this%VELOCITY  )%ccomp
    vy => uin(this%VELOCITY+1)%ccomp
    vz => uin(this%VELOCITY+2)%ccomp
    th => uin(this%TEMP)      %ccomp
    fx => uf (this%VELOCITY  )%ccomp
    fy => uf (this%VELOCITY+1)%ccomp
    fz => uf (this%VELOCITY+2)%ccomp
    CALL hdcheck(vx,vy,vz,fx,fy,fz,t,dt,1,0)
    CALL maxabs(vx,vy,vz,rmp,0)
    IF (myrank.eq.0) THEN
      OPEN(1,file='maximum.txt',position='append')
      WRITE(1,FMT='(E13.6,E13.6)') (t-1)*dt,rmp
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
    use iovar
    implicit none

    class (BOUSSSolver), intent(in)                :: this
    type    (GState), intent(in), target        :: uin(:)
    complex(kind=GP), pointer, dimension(:,:,:) :: vx,vy,vz,th
    complex(kind=GP), pointer, dimension(:,:,:) :: c1,c2,c3
    integer                                     :: i
    logical                                     :: bret

    WRITE(ext, fmtext) sind
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
    ! Uncomment the following line to compute the helicity flux
    ! call heltrans(vx,vy,vz,-c1,-c2,-c3,ext,1)
    if ( this%traits_%dorot ) then
      call specpara(vx,vy,vz,ext,1,1)
      call specperp(vx,vy,vz,ext,1,1)
      call entpara(vx,vy,vz,-c1,-c2,-c3,ext,1) ! Writes the energy flux
      call entperp(vx,vy,vz,-c1,-c2,-c3,ext,1) ! Writes the energy fluxq
      ! The following two lines compute anisotropic helicity fluxes
      ! call helpara(vx,vy,vz,-c1,-c2,-c3,ext,1)
      ! call helperp(vx,vy,vz,-c1,-c2,-c3,ext,1)
      ! Uncomment the following line to compute 2D spectra
      ! CALL spec2D(vx,vy,vz,ext,odir,1,1)
      ! Uncomment the following lines to compute spatio-temporal spectra
      ! CALL write_fourier(vx,'vx',ext,odir)
      ! CALL write_fourier(vy,'vy',ext,odir)
      ! CALL write_fourier(vz,'vz',ext,odir)
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
  subroutine BOUSSSolver_ctor(this, infile, workspace, plan)
    use iovar
    class  (BOUSSSolver), intent(inout)         :: this
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
  subroutine BOUSSSolver_dtor(this) 
    type  (BOUSSSolver), intent(inout) :: this
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
