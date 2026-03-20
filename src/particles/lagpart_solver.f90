! =====================================================================
! NAME       : lagpart_solver.f90
! DESCRIPTION: Forms solver class for Lagrangian particles, computing:
!
!              dx/dt = v(x,t)
!
!              State ordering is:
!                x1 (x2, x3)
!
!              State sector ids are:
!                POSITION (POSITION+1, POSITION+2)
!
! INPUT FILE : For particles='lagpart', looks for "&lagpart" namelist:
!
!
! DATE       : 19/02/26 (PDM)
! =====================================================================

module lagpart_mod
  use particlebase_mod
  use gpstate_mod

  implicit none

  ! ================= Solver traits ===================================
  type, public  :: NHTraits
  end type

  ! ================= Global parameters ===============================
  
  ! ================= Solver ==========================================
  ! Define class:
  type, extends(VelocParticleBase) :: GPart 
    ! Member data:
    logical           :: binit_ = .false. ! is initialized?
    type  (NHTraits)  :: traits_
  CONTAINS
    procedure, public :: init          =>          init_impl ! init method
    procedure, public :: dpdt          =>          dpdt_impl ! part RHS method
!   procedure, public :: end_step      =>      end_step_impl ! finalizes time evolution
    procedure, public :: feedback      =>      null_feedback ! feedback in the fluid
    procedure, public :: write_pstate  =>  write_pstate_impl ! write states
    procedure, public :: state_size    =>    state_size_impl ! state size
!   procedure, public :: sstate2istate => sstate2istate_impl ! state names
!   procedure, public :: get_sstate    =>    get_sstate_impl ! get state name list
    procedure, public :: part_ctor     =>         GPart_ctor ! constructor
    final             :: GPart_dtor
  end type GPart

CONTAINS

  ! ===================================================================
  ! Solver initialization, this is where parameter files are read
  ! ===================================================================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Subroutine to initialize the solver
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_impl(this)
    use commtypes
    class  (GPart), intent (inout) :: this
    this%POSITION = 1                        ! start of position sector
  end subroutine init_impl
  
  ! ===================================================================
  ! Computation of RHS, the solver equations are defined here
  ! ===================================================================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Function to compute the rhs of the equations of motion
  !! of Lagrangian particles.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE dpdt_impl(this, time, pde, fluidstate, pstate, dt, dpdtout)
    use equationbase_mod
    use fft
    IMPLICIT NONE
    class       (GPart),         intent(inout) :: this
    class(EquationBase),         intent   (in) :: pde
    real      (kind=GP),         intent   (in) :: time, dt
    type   (GStateComp), target, intent   (in) :: fluidstate(:)
    type  (GPStateComp),         intent   (in) :: pstate(:)
    type  (GPStateComp),         intent(inout) :: dpdtout(:)
    complex   (KIND=GP), pointer,DIMENSION(:,:,:) :: velc
    real      (KIND=GP), pointer,DIMENSION(:,:,:) :: velr,tmp1,tmp2
    real      (kind=GP)                           :: rmp
    integer                                       :: i,j,k,m
    logical                                       :: bret
    
    CALL GTStart(this%htimers_(GPTIME_STEP))
    call this%workspace_%get_complex_tmp(velc,bret)
    call this%workspace_%get_real_tmp   (velr,bret)
    call this%workspace_%get_real_tmp   (tmp1,bret)
    call this%workspace_%get_real_tmp   (tmp2,bret)

    select type (pde)
    class is (VelocityBase)
      ! Find F(u*):
      do m = pde%VELOCITY, pde%VELOCITY+pde%nc_-1
      !$omp parallel do if (iend-ista.ge.nth) private (j,k)
        do i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
          do j = 1,ny
            do k = 1,nz
              velc(k,j,i) = fluidstate(m)%ccomp(k,j,i)*rmp
            end do
          end do
        end do
        call fftp3d_complex_to_real(plancr,velc,velr,MPI_COMM_WORLD)
        call EulerToLag(this,dpdtout(m)%rcomp,this%nparts_,velr,.true.,tmp1,tmp2)
      end do
    class default
      stop "lagpart: This solver does not support pdes without a velocity field"
    end select
 
    call this%workspace_%free_complex_tmp(velc)
    call this%workspace_%free_real_tmp   (velr)
    call this%workspace_%free_real_tmp   (tmp1)
    call this%workspace_%free_real_tmp   (tmp2)
    CALL GTAcc(this%htimers_(GPTIME_STEP))    
  END SUBROUTINE dpdt_impl

    
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Functions to compute fluid and particle couplings
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine null_feedback(this, pstate, feedback)
    class     (GPart), intent   (in) :: this
    type(GPStateComp), intent   (in) :: pstate
    type (GStateComp), intent(inout) :: feedback
    return
  end subroutine null_feedback
  
  ! ===================================================================
  ! Output methods 
  ! ===================================================================
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Function to compute and write particle states
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine write_pstate_impl(this, time, pde, fluidstate, pstate)
    use equationbase_mod
    use particlebase_mod
    
    class       (GPart),         intent(inout) :: this
    class(EquationBase),         intent   (in) :: pde
    real      (kind=GP),         intent   (in) :: time
    type   (GStateComp), target, intent   (in) :: fluidstate(:)
    type  (GPStateComp),         intent   (in) :: pstate(:)
    real      (kind=GP) :: rmp
    integer             :: i,j,k

!!$    ! Set the Lagrangian velocities so output doesn't
!!$    ! give 0: 
!!$    rmp = 1.0_GP/(real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))
!!$!$omp parallel do if (iend-ista.ge.nth) private (j,k)
!!$    DO i = ista,iend
!!$!$omp parallel do if (iend-ista.lt.nth) private (k)
!!$      DO j = 1,ny
!!$        DO k = 1,nz
!!$          C7(k,j,i) = vx(k,j,i)*rmp
!!$        END DO
!!$      END DO
!!$    END DO
!!$    CALL fftp3d_complex_to_real(plancr,C7,R1,MPI_COMM_WORLD)
!!$!$omp parallel do if (iend-ista.ge.nth) private (j,k)
!!$    DO i = ista,iend
!!$!$omp parallel do if (iend-ista.lt.nth) private (k)
!!$      DO j = 1,ny
!!$        DO k = 1,nz
!!$          C7(k,j,i) = vy(k,j,i)*rmp
!!$        END DO
!!$      END DO
!!$    END DO
!!$    CALL fftp3d_complex_to_real(plancr,C7,R2,MPI_COMM_WORLD)
!!$!$omp parallel do if (iend-ista.ge.nth) private (j,k)
!!$      DO i = ista,iend
!!$!$omp parallel do if (iend-ista.lt.nth) private (k)
!!$        DO j = 1,ny
!!$          DO k = 1,nz
!!$            C7(k,j,i) = vz(k,j,i)*rmp
!!$          END DO
!!$        END DO
!!$      END DO
!!$    CALL fftp3d_complex_to_real(plancr,C7,R3,MPI_COMM_WORLD)
!!$    CALL lagpart%SetLagVec(R1,R2,R3,.true.,R4,R5)
!!$
!!$    timep = 0
!!$    pind = pind+1
!!$    WRITE(lgext,lgfmtext) pind
!!$
!!$!!!!!! Write Lagrangian vorticity components: !!!!!!
!!$!
!!$!$omp parallel do if (iend-ista.ge.nth) private (j,k)
!!$    DO i = ista,iend
!!$!$omp parallel do if (iend-ista.lt.nth) private (k)
!!$      DO j = 1,ny
!!$        DO k = 1,nz
!!$          C1(k,j,i) = vx(k,j,i)*rmp
!!$          C2(k,j,i) = vy(k,j,i)*rmp
!!$          C3(k,j,i) = vz(k,j,i)*rmp
!!$        END DO
!!$      END DO
!!$    END DO
!!$    CALL rotor3(C2,C3,C4,1)
!!$    CALL rotor3(C1,C3,C5,2)
!!$    CALL rotor3(C1,C2,C6,3)
!!$    CALL fftp3d_complex_to_real(plancr,C4,R1,MPI_COMM_WORLD)
!!$    CALL fftp3d_complex_to_real(plancr,C5,R2,MPI_COMM_WORLD)
!!$    CALL fftp3d_complex_to_real(plancr,C6,R3,MPI_COMM_WORLD)
!!$    CALL lagpart%SetLagVec(R1,R2,R3,.false.,R4,R5)
!!$    CALL lagpart%io_write_vec(1,odir,'wlg'  ,lgext,(t-1)*dt)
!!$
!!$    nwpart = nwpart + 1
  end subroutine write_pstate_impl
        
  ! ===================================================================
  ! Solver specific methods
  ! ===================================================================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Constructor
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE GPart_ctor(this,infile, workspace, pstate, pstate_cpy)
    USE var
    USE grid
    USE boxsize
    USE mpivars
    USE pstatus
    USE commtypes
    USE fftplans
    USE pstatus
    USE random
    IMPLICIT NONE
    CLASS     (GPart), intent(inout)                      :: this
    type (GWorkspace), intent(inout),              target :: workspace
    type(GPStateComp), intent(inout), allocatable, target :: pstate(:), pstate_cpy(:)
    character(len=*) , intent   (in)                      :: infile
    INTEGER                             :: disp(3),lens(3),types(3),szreal
    INTEGER                             :: tsta,tend,num_components
    INTEGER                             :: j,nc
    logical                             :: bret

    this%infile_      =  infile    ! input file
    this%workspace_   => workspace
    call pstatus_init(this%infile_)
    this%nc_          = 3          ! fixed for now
    this%nparts_      = 0 
    this%npartsm_     = 0 
    this%nvdb_        = 0
    this%comm_        = MPI_COMM_WORLD
    this%maxparts_    = maxparts
    this%nd_(1)       = nx
    this%nd_(2)       = ny
    this%nd_(3)       = nz
    this%delta_(1)    = 2*pi*Lx/real(nx,kind=GP)
    this%delta_(2)    = 2*pi*Ly/real(ny,kind=GP)
    this%delta_(3)    = 2*pi*Lz/real(nz,kind=GP)
    this%invdel_(1)   = real(nx,kind=GP)/(2*pi*Lx)
    this%invdel_(2)   = real(ny,kind=GP)/(2*pi*Ly)
    this%invdel_(3)   = real(nz,kind=GP)/(2*pi*Lz)
    this%seedfile_    = lgseedfile
    this%iinterp_     = 3          ! fixed for now
    this%itorder_     = 2
    this%intorder_    = max(intorder,1)
    this%iseed_       = 1000
    this%istep_       = 0   
    this%iexchtype_   = ilgexchtype
    this%iouttype_    = ilgouttype
    this%bcollective_ = ilgcoll
    this%itimetype_   = GT_WTIME
    this%wrtunit_     = ilgwrtunit
    CALL prandom_seed(this%iseed_)
    CALL MPI_COMM_SIZE(this%comm_,this%nprocs_,this%ierr_)
    CALL MPI_COMM_RANK(this%comm_,this%myrank_,this%ierr_)
 
    IF (this%iexchtype_.EQ.GPEXCHTYPE_VDB) THEN
      this%partbuff_ = maxparts  
    ELSE IF (this%iexchtype_.EQ.GPEXCHTYPE_NN) THEN
      this%partbuff_      = 1 + (maxparts - 1)/this%nprocs_
      this%partchunksize_ = (this%partbuff_ + 9)/10
      this%partbuff_      =  this%partbuff_ + this%partchunksize_
      IF ((this%bcollective_.EQ.0).AND.(this%myrank_.EQ.0)) THEN
        this%partbuff_   = maxparts
      END IF
      this%stepcounter_ = 0
    END IF
   
    ! Initialize timers (get handles):
    DO j = 1, GPMAXTIMERS
      CALL GTInitHandle(this%htimers_(j),this%itimetype_)
      IF ( this%htimers_(j).EQ.GTNULLHANDLE ) THEN
        WRITE(*,*) 'GPart_ctor: Not enough timers available'
        STOP
      ENDIF
    ENDDO

    ! Initialize communicators
    CALL this%gpcomm_%GPartComm_ctor(GPCOMM_INTRFC_SF,this%partbuff_, &
         this%nd_,this%intorder_-1,this%comm_,this%htimers_(GPTIME_COMM))
    CALL this%gpcomm_%SetCacheParam(csize,nstrip)
    CALL this%gpcomm_%Init()

    this%libnds_(1,1) = 1
    this%libnds_(1,2) = nx
    this%lxbnds_(1,1) = 0.0_GP
    this%lxbnds_(1,2) = real(nx,kind=GP)
    this%libnds_(2,1) = 1
    this%libnds_(2,2) = ny
    this%lxbnds_(2,1) = 0.0_GP
    this%lxbnds_(2,2) = real(ny,kind=GP)
    this%libnds_(3,1) = ksta 
    this%libnds_(3,2) = kend 
    this%lxbnds_(3,1) = real(ksta-1,kind=GP)          !- 0.50_GP
    this%lxbnds_(3,2) = real(kend-1,kind=GP) + 1.0_GP !0.50_GP
    CALL range(1,nx,nprocs,myrank,tsta,tend) !Bounds of transposed real array
    this%tibnds_(1,1) = 1
    this%tibnds_(1,2) = nz
    this%tibnds_(2,1) = 1
    this%tibnds_(2,2) = ny
    this%tibnds_(3,1) = tsta 
    this%tibnds_(3,2) = tend
    DO j = 1,3
      this%gext_ (j) = real(this%nd_(j),kind=GP)
    ENDDO

    ! Call init
    call this%init()

    ! Instantiate interp operation. Remember that a valid timer 
    ! handle must be passed:
    CALL this%intop_%GPSplineInt_ctor(3,this%nd_,this%libnds_,this%lxbnds_, &
         this%tibnds_,this%intorder_,this%partbuff_,this%gpcomm_,&
         this%htimers_(GPTIME_DATAEX),this%htimers_(GPTIME_TRANSP))

    ! Create part. d.b. structure type for time evolution and I/O
    CALL MPI_TYPE_SIZE(GC_REAL,szreal,this%ierr_)
    ALLOCATE(this%id_      (this%partbuff_))
    ALLOCATE(this%tmpint_  (this%partbuff_))
    num_components = this%state_size()
    CALL GPState_alloc(pstate    , num_components, this%partbuff_)
    CALL GPState_alloc(pstate_cpy, num_components, this%partbuff_)
    this%px_ => pstate(this%POSITION  )%rcomp
    this%py_ => pstate(this%POSITION+1)%rcomp
    this%pz_ => pstate(this%POSITION+2)%rcomp
    call this%workspace_%set_nparts(this%partbuff_)
    call this%workspace_%get_pcomp_tmp(this%lvx_,bret)    
    call this%workspace_%get_pcomp_tmp(this%lvy_,bret)    
    call this%workspace_%get_pcomp_tmp(this%lvz_,bret)    
    ALLOCATE(this%ptmp0_ (3,this%partbuff_))
    IF ( this%iexchtype_.EQ.GPEXCHTYPE_VDB ) THEN
      ALLOCATE(this%vdb_ (3,this%partbuff_))
    ENDIF
  END SUBROUTINE GPart_ctor

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Destructor
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE GPart_dtor(this)
    IMPLICIT NONE
    TYPE(GPart)   ,INTENT(INOUT) :: this
    integer                      :: j
    IF ( ALLOCATED    (this%id_) ) DEALLOCATE    (this%id_)
    IF ( ALLOCATED(this%tmpint_) ) DEALLOCATE(this%tmpint_)
    IF ( ALLOCATED   (this%idm_) ) DEALLOCATE   (this%idm_)
    IF ( ALLOCATED   (this%vdb_) ) DEALLOCATE   (this%vdb_)
    IF ( ALLOCATED (this%ptmp0_) ) DEALLOCATE (this%ptmp0_)
    IF ( ASSOCIATED   (this%px_) ) NULLIFY       (this%px_)
    IF ( ASSOCIATED   (this%py_) ) NULLIFY       (this%py_)
    IF ( ASSOCIATED   (this%pz_) ) NULLIFY       (this%pz_)
    ! Destroy timers:
    DO j = 1, GPMAXTIMERS
      CALL GTFree(this%htimers_(j))
    ENDDO
  END SUBROUTINE GPart_dtor


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Convert input state name to index in state vector
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Get state variable names
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Function to compute number of state members (equations)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PURE function state_size_impl(this) result(num)
    class(GPart), intent(in)    :: this
    integer                     :: num
    num = this%nc_              ! # spatial coordinates
  end function state_size_impl
  
end module lagpart_mod
