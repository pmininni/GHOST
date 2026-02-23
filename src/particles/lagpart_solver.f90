! =====================================================================
! NAME       : lagpart_solver.f90
! DESCRIPTION: Forms solver class for Lagrangian particles, computing:
!
!              dx/dt = v(x,t)
!
!              State ordering is:
!                x1, x2, x3
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
    logical           :: binit_=.false. ! is initialized?
    type  (NHTraits)  :: traits_
  CONTAINS
    procedure, public :: init          =>          init_impl ! init method
    procedure, public :: dvpdt         =>         dvpdt_impl ! part RHS method
    procedure, public :: state_size    =>    state_size_impl ! state size
    procedure, public :: sstate2istate => sstate2istate_impl ! state names
    procedure, public :: get_sstate    =>    get_sstate_impl ! get state name list
    procedure, public :: Solver_ctor   =>         GPart_ctor ! constructor
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





    
  ! ===================================================================
  ! Computation of RHS, the solver equations are defined here
  ! ===================================================================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Function to compute RHS
!  METHOD     : Step
!  DESCRIPTION: Carries out one stage of explicit RK-like time
!               inegration step.  Intended for explicit step within 
!               an outer stepper method of the form:
!
!               X = X_0 + dt * V(X(t),t) * xk,
!       
!               Note that the vx, vy, vz, will be overwritten here.
!  ARGUMENTS  :
!    this    : 'this' class instance
!    vz,vy,vz: compoments of velocity field, in real space, partially
!              updated, possibly. These will be overwritten!
!    dt      : integration timestep
!    xk      : multiplicative RK time stage factor
!    tmpX    : REAL temp arrays the same size as vx, vy, vz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE GPart_Step(this, vx, vy, vz, dt, xk, tmp1, tmp2, tmp3)
    USE fprecision
    USE commtypes
    USE mpivars
    USE grid

    IMPLICIT NONE
    CLASS(GPart) ,INTENT(INOUT)                            :: this
    INTEGER                                                :: i,j
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(nx,ny,ksta:kend) :: vx,vy,vz
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(nx,ny,ksta:kend) :: tmp1,tmp2,tmp3
    REAL(KIND=GP),INTENT   (IN)                            :: dt,xk
    REAL(KIND=GP)                                          :: dtfact
    REAL(KIND=GP),ALLOCATABLE  ,DIMENSION              (:) :: lid,gid

    CALL GTStart(this%htimers_(GPTIME_STEP))

    ! Find F(u*):
    ! ... x:
    CALL GPart_R3toR3(this,tmp3,vx) ! Want vx intact to use later
    CALL GPart_EulerToLag(this,this%lvx_,this%nparts_,tmp3,.true.,tmp1,tmp2)
    ! ux* <-- ux + dt * F(U*)*xk:
    dtfact = dt*xk*this%invdel_(1)
!$omp parallel do
    DO j = 1, this%nparts_
      this%px_(j) = this%ptmp0_(1,j) + dtfact*this%lvx_(j)
    ENDDO
    
    !
    ! ... y:
    ! Exchange bdy data for velocities, so that we
    ! can perform local interpolations:
    CALL GPart_R3toR3(this,tmp3,vy) ! Want vy intact to use later
    CALL GPart_EulerToLag(this,this%lvy_,this%nparts_,tmp3,.false.,tmp1,tmp2)
    ! uy* <-- uy + dt * F(U*)*xk:
    dtfact = dt*xk*this%invdel_(2)
!$omp parallel do
    DO j = 1, this%nparts_
      this%py_(j) = this%ptmp0_(2,j) + dtfact*this%lvy_(j)
    ENDDO

    ! ... z:
    ! Exchange bdy data for velocities, so that we
    ! can perform local interpolations:
    CALL GPart_R3toR3(this,tmp3,vz) ! Want vz intact to use later
    CALL GPart_EulerToLag(this,this%lvz_,this%nparts_,tmp3,.false.,tmp1,tmp2)
    ! uz* <-- uz + dt * F(U*)*xk:
    dtfact = dt*xk*this%invdel_(3)
!$omp parallel do
    DO j = 1, this%nparts_
      this%pz_(j) = this%ptmp0_(3,j) + dtfact*this%lvz_(j)
    ENDDO

    ! Enforce periodicity in x-y only:
    CALL GPart_MakePeriodicP(this,this%px_,this%py_,this%pz_,this%nparts_,3)

    CALL GTAcc(this%htimers_(GPTIME_STEP))

!   ALLOCATE  (lid(this%maxparts_))
!   ALLOCATE  (gid(this%maxparts_))
!   lid = 0
!   gid = 0
!   do j=1,this%maxparts_
!     lid(this%id_(j)+1) = 1
!   enddo
!   call mpi_allreduce(lid,gid,this%maxparts_,MPI_INTEGER,MPI_SUM,this%comm_,this%ierr_)
!   if ( this%myrank_.eq. 0 ) then
!     do j=1,this%maxparts_
!       if ( gid(j) .gt.1 ) then
!         write(*,*)this%myrank_,': StepRKK: multiplicity > 1: id=',j-1,' p=', &
!                   this%vdb_(1,j),this%vdb_(2,j),this%vdb_(3,j)
!       else if ( gid(j).eq.0 ) then
!         write(*,*)this%myrank_,': StepRKK: particle missing: id=',j-1
!       endif    
!     enddo
!   endif
!
!   DEALLOCATE(lid,gid)

  ! At this point, the vx, vy, vz should be intact:
  CALL GPart_EndStageRKK(this,vx,vy,vz,xk,tmp1,tmp2)

  END SUBROUTINE GPart_StepRKK
!-----------------------------------------------------------------
!-----------------------------------------------------------------
  

    
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Functions to compute fluid and particle couplings
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  

  ! ===================================================================
  ! Output methods 
  ! ===================================================================
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Function to compute and write particle states
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Set the Lagrangian velocities so output doesn't
! give 0: 
!
           rmp = 1.0_GP/(real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
           DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
             DO j = 1,ny
               DO k = 1,nz
                 C7(k,j,i) = vx(k,j,i)*rmp
               END DO
             END DO
           END DO
           CALL fftp3d_complex_to_real(plancr,C7,R1,MPI_COMM_WORLD)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
           DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
             DO j = 1,ny
               DO k = 1,nz
                 C7(k,j,i) = vy(k,j,i)*rmp
               END DO
             END DO
           END DO
           CALL fftp3d_complex_to_real(plancr,C7,R2,MPI_COMM_WORLD)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
           DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
             DO j = 1,ny
               DO k = 1,nz
                 C7(k,j,i) = vz(k,j,i)*rmp
               END DO
             END DO
           END DO
           CALL fftp3d_complex_to_real(plancr,C7,R3,MPI_COMM_WORLD)
           CALL lagpart%SetLagVec(R1,R2,R3,.true.,R4,R5)

           timep = 0
           pind = pind+1
           WRITE(lgext,lgfmtext) pind

!       NOTE: if dopacc > 0 (set in solver), write out position and velocity corresp. 
!             to acceleration time stamp; change name to reflect that these lag by one 
!             dt from the current time:
           CALL lagpart%io_write_pdb  (1,odir,'xlg'  ,lgext,(t-1)*dt)
           CALL lagpart%io_write_vec  (1,odir,'vlg'  ,lgext,(t-1)*dt)
           CALL lagpart%io_write_pdbm1(1,odir,'xlgm1',lgext,(t-2)*dt)
           CALL lagpart%io_write_vecm1(1,odir,'vlgm1',lgext,(t-2)*dt)
!!!!!! Write internal Lagrangian acceleration components: !!!!!!
!      NOTE: if dopacc > 0, then time centering of 'vlgm1', xlgm1' are
!            the same as 'algm1', while the other quantities are at the
!            most recent time.] 
           CALL lagpart%io_write_acc(tbeta,1,odir,'algm1',lgext,(t-2)*dt)
!
!!!!!! Write Lagrangian vorticity components: !!!!!!
!
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
           DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
              DO j = 1,ny
                 DO k = 1,nz
                    C1(k,j,i) = vx(k,j,i)*rmp
                    C2(k,j,i) = vy(k,j,i)*rmp
                    C3(k,j,i) = vz(k,j,i)*rmp
                 END DO
              END DO
           END DO
           CALL rotor3(C2,C3,C4,1)
           CALL rotor3(C1,C3,C5,2)
           CALL rotor3(C1,C2,C6,3)
           CALL fftp3d_complex_to_real(plancr,C4,R1,MPI_COMM_WORLD)
           CALL fftp3d_complex_to_real(plancr,C5,R2,MPI_COMM_WORLD)
           CALL fftp3d_complex_to_real(plancr,C6,R3,MPI_COMM_WORLD)
           CALL lagpart%SetLagVec(R1,R2,R3,.false.,R4,R5)
           CALL lagpart%io_write_vec(1,odir,'wlg'  ,lgext,(t-1)*dt)

           nwpart = nwpart + 1

!
!!!!!! Write nonlinear terms:!!!!!!
!
           if ( .false. ) then
           CALL prodre3(vx,vy,vz,C4,C5,C6)
           CALL fftp3d_complex_to_real(plancr,C4,R4,MPI_COMM_WORLD)
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
           DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
             DO j = 1,ny
                DO i = 1,nz
                   R4(i,j,k) = R4(i,j,k)*rmp
                END DO
             END DO
           END DO
           CALL lagpart%io_write_euler(1,odir,'v1nllg'  ,lgext,(t-1)*dt,R4,.false.,R2,R3)

           CALL fftp3d_complex_to_real(plancr,C5,R4,MPI_COMM_WORLD)
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
           DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
             DO j = 1,ny
                DO i = 1,nz
                   R4(i,j,k) = R4(i,j,k)*rmp
                END DO
             END DO
           END DO
           CALL lagpart%io_write_euler(1,odir,'v2nllg'  ,lgext,(t-1)*dt,R4,.false.,R2,R3)

           CALL fftp3d_complex_to_real(plancr,C6,R4,MPI_COMM_WORLD)
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
           DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
             DO j = 1,ny
                DO i = 1,nz
                   R4(i,j,k) = R4(i,j,k)*rmp
                END DO
             END DO
           END DO
           CALL lagpart%io_write_euler(1,odir,'v3nllg'  ,lgext,(t-1)*dt,R4,.false.,R2,R3)
           endif

!
!!!!!!! Write strain-rate tensor components: !!!!!!
!
           CALL derivk3(vx,C1,1)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
           DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
              DO j = 1,ny
                 DO k = 1,nz
                    C1(k,j,i) = C1(k,j,i)*rmp
                 END DO
              END DO
           END DO
           CALL fftp3d_complex_to_real(plancr,C1,R1,MPI_COMM_WORLD)
           CALL lagpart%io_write_euler(1,odir,'s11'  ,lgext,(t-1)*dt,R1,.false.,R2,R3)
           CALL derivk3(vx,C1,2)
           CALL derivk3(vy,C2,1)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
           DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
              DO j = 1,ny
                 DO k = 1,nz
                    C1(k,j,i) = 0.5*(C1(k,j,i)+C2(k,j,i))*rmp
                 END DO
              END DO
           END DO
           CALL fftp3d_complex_to_real(plancr,C1,R1,MPI_COMM_WORLD)
           CALL lagpart%io_write_euler(1,odir,'s12',lgext,(t-1)*dt,R1,.false.,R2,R3)
           CALL derivk3(vx,C1,3)
           CALL derivk3(vz,C2,1)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
           DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
              DO j = 1,ny
                 DO k = 1,nz
                    C1(k,j,i) = 0.5*(C1(k,j,i)+C2(k,j,i))*rmp
                 END DO
              END DO
           END DO
           CALL fftp3d_complex_to_real(plancr,C1,R1,MPI_COMM_WORLD)
           CALL lagpart%io_write_euler(1,odir,'s13',lgext,(t-1)*dt,R1,.false.,R2,R3)
           CALL derivk3(vy,C1,2)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
           DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
              DO j = 1,ny
                 DO k = 1,nz
                    C1(k,j,i) = C1(k,j,i)*rmp
                 END DO
              END DO
           END DO
           CALL fftp3d_complex_to_real(plancr,C1,R1,MPI_COMM_WORLD)
           CALL lagpart%io_write_euler(1,odir,'s22',lgext,(t-1)*dt,R1,.false.,R2,R3)
           CALL derivk3(vy,C1,3)
           CALL derivk3(vz,C2,2)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
           DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
              DO j = 1,ny
                 DO k = 1,nz
                    C1(k,j,i) = 0.5*(C1(k,j,i)+C2(k,j,i))*rmp
                 END DO
              END DO
           END DO
           CALL fftp3d_complex_to_real(plancr,C1,R1,MPI_COMM_WORLD)
           CALL lagpart%io_write_euler(1,odir,'s23',lgext,(t-1)*dt,R1,.false.,R2,R3)

  
  ! ===================================================================
  ! Solver specific methods
  ! ===================================================================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Constructor
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE GPart_ctor(this,comm,mparts,inittype,iinterp,intorder,iexchtyp, &
                        iouttyp,bcoll,csize,nstrip,intacc,wrtunit)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  Main explicit constructor
!  ARGUMENTS:
!    this    : 'this' class instance
!    comm    : MPI communicator
!    mparts  : max no. particles allowed on grid. If inittype = GPINIT_USERLOC
!              then this is computed from the user file, and can
!              be any value.
!    inittype: GPINIT-typed quantity to give the type of particle
!              initialization
!    iinterp : GPINTRP-type quantity providing the interpolation
!              scheme
!    intorder: for variable-order (e.g., Lagrange) interpolation, 
!              the order (1, 2, 3...). Sets the number of 'ghost' zones
!              of data transferred between MPI tasks.
!    iexchtyp: format for exchanging particle data: GPEXCHTYPE_NN=nearest neighbor, 
!              suggesting that data is passed between neighboring
!              tasks only; GPEXCHTYPE_VDB form suggests that all tasks retain a copy
!              of the 'master' particle d.b., so that passing between
!              neighbors is not necessary. But the VDB form does require
!              a global reduction, and may be more expensive.
!    iouttup : output type: 0==binary, 1=ASCII
!    bcoll   : if doing binary I/O, do collective (==1); or not (==0)
!    csize   : cache size param for local transposes
!    nstrip  : 'strip-mining' size for local transposes
!    intacc  : compute acceleration internally to class (==1); or not (==0). Storage
!              allocated only if intacc==1.
!    wrtunit : (optional) write particle positions in box units (==1) (i.e.,
!              x,y,z in [0,2.pi]), or in grid units (==0) (x,y,z in [0,N]).
!-----------------------------------------------------------------
    USE var
    USE grid
    USE boxsize
    USE mpivars
    USE commtypes
    USE random

    IMPLICIT NONE
    CLASS(GPart)     ,INTENT(INOUT)     :: this
    INTEGER          ,INTENT   (IN)     :: bcoll,comm,mparts
    INTEGER          ,INTENT   (IN)     :: csize,nstrip,intacc
    INTEGER, OPTIONAL,INTENT   (IN)     :: wrtunit
    INTEGER                             :: disp(3),lens(3),types(3),szreal
    INTEGER          ,INTENT   (IN)     :: iexchtyp,iinterp,inittype
    INTEGER          ,INTENT   (IN)     :: intorder,iouttyp
    INTEGER                             :: tsta,tend
    INTEGER                             :: j,nc

    this%nparts_   = 0 
    this%npartsm_  = 0 
    this%nvdb_     = 0
    this%comm_     = comm
    this%maxparts_ = mparts
    this%nd_(1)    = nx
    this%nd_(2)    = ny
    this%nd_(3)    = nz
    this%delta_(1) = 2*pi*Lx/real(nx,kind=GP)
    this%delta_(2) = 2*pi*Ly/real(ny,kind=GP)
    this%delta_(3) = 2*pi*Lz/real(nz,kind=GP)
    this%invdel_(1)= real(nx,kind=GP)/(2*pi*Lx)
    this%invdel_(2)= real(ny,kind=GP)/(2*pi*Ly)
    this%invdel_(3)= real(nz,kind=GP)/(2*pi*Lz)
    this%seedfile_ = 'gploc.dat'
    this%iinterp_  = 3          ! fixed for now
    this%inittype_ = inittype
    this%itorder_  = 2
    this%intorder_ = max(intorder,1)
    this%iseed_    = 1000
    this%istep_    = 0   
    this%iexchtype_   =  iexchtyp
    this%iouttype_    =  iouttyp
    this%bcollective_ =  bcoll
    this%itimetype_   =  GT_WTIME
    this%intacc_      =  intacc
    IF ( present(wrtunit) ) THEN
       this%wrtunit_  =  wrtunit
    ELSE
       this%wrtunit_  =  0
    ENDIF

    CALL prandom_seed(this%iseed_)
    IF ( this%intorder_ .NE. 3 ) THEN
      WRITE(*,*) 'GPart::ctor: Only 3rd order allowed for now' 
    ENDIF

    CALL MPI_COMM_SIZE(this%comm_,this%nprocs_,this%ierr_)
    CALL MPI_COMM_RANK(this%comm_,this%myrank_,this%ierr_)
 
    IF (iexchtyp.EQ.GPEXCHTYPE_VDB) THEN
      this%partbuff_ = mparts  
    ELSE IF (iexchtyp.EQ.GPEXCHTYPE_NN) THEN
      this%partbuff_      = 1 + (mparts - 1)/this%nprocs_
      this%partchunksize_ = (this%partbuff_ + 9)/10
      this%partbuff_      =  this%partbuff_ + this%partchunksize_
      IF ((bcoll.EQ.0).AND.(this%myrank_.EQ.0)) THEN
        this%partbuff_   = mparts
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

    CALL this%gpcomm_%GPartComm_ctor(GPCOMM_INTRFC_SF,this%partbuff_, &
         this%nd_,this%intorder_-1,this%comm_,this%htimers_(GPTIME_COMM))
    CALL this%gpcomm_%SetCacheParam(csize,nstrip)
    CALL this%gpcomm_%Init()

    this%libnds_(1,1) = 1  ;
    this%libnds_(1,2) = nx ;
    this%lxbnds_(1,1) = 0.0_GP
    this%lxbnds_(1,2) = real(nx,kind=GP)
    this%libnds_(2,1) = 1  ;
    this%libnds_(2,2) = ny ;
    this%lxbnds_(2,1) = 0.0_GP
    this%lxbnds_(2,2) = real(ny,kind=GP)
    this%libnds_(3,1) = ksta ; 
    this%libnds_(3,2) = kend ; 
    this%lxbnds_(3,1) = real(ksta-1,kind=GP)          !- 0.50_GP
    this%lxbnds_(3,2) = real(kend-1,kind=GP) + 1.0_GP !0.50_GP
    CALL range(1,nx,nprocs,myrank,tsta,tend) !Bounds of transposed real array
    this%tibnds_(1,1) = 1  ;
    this%tibnds_(1,2) = nz ;
    this%tibnds_(2,1) = 1  ;
    this%tibnds_(2,2) = ny ;
    this%tibnds_(3,1) = tsta ; 
    this%tibnds_(3,2) = tend ; 

    DO j = 1,3
      this%gext_ (j) = real(this%nd_(j),kind=GP)
    ENDDO

    ! Instantiate interp operation. Remember that a valid timer 
    ! handle must be passed:
    CALL this%intop_%GPSplineInt_ctor(3,this%nd_,this%libnds_,this%lxbnds_, &
         this%tibnds_,this%intorder_,this%partbuff_,this%gpcomm_,&
         this%htimers_(GPTIME_DATAEX),this%htimers_(GPTIME_TRANSP))

    ! Create part. d.b. structure type for I/O
    CALL MPI_TYPE_SIZE(GC_REAL,szreal,this%ierr_)

    ALLOCATE(this%id_      (this%partbuff_))
    ALLOCATE(this%tmpint_  (this%partbuff_))
    IF ( this%intacc_.EQ.1 ) THEN
    ALLOCATE(this%idm_     (this%partbuff_))
    ENDIF
    ALLOCATE(this%px_      (this%partbuff_))
    ALLOCATE(this%py_      (this%partbuff_))
    ALLOCATE(this%pz_      (this%partbuff_))
    ALLOCATE(this%lvx_     (this%partbuff_))
    ALLOCATE(this%lvy_     (this%partbuff_))
    ALLOCATE(this%lvz_     (this%partbuff_))
    ALLOCATE(this%ptmp0_ (3,this%partbuff_))
    ALLOCATE(this%ptmp1_ (3,this%partbuff_))
    IF ( this%iexchtype_.EQ.GPEXCHTYPE_VDB ) THEN
      ALLOCATE(this%gptmp0_ (3,this%partbuff_))
      ALLOCATE(this%vdb_ (3,this%partbuff_))
    ENDIF
    IF ( this%intacc_.EQ. 1 ) THEN
      ALLOCATE(this%vk0_  (3,this%partbuff_))
      ALLOCATE(this%vk1_  (3,this%partbuff_))
      ALLOCATE(this%vk2_  (3,this%partbuff_))
      ALLOCATE(this%xk1_  (3,this%partbuff_))
      ALLOCATE(this%ptmp2_(3,this%partbuff_))
    ENDIF
    ALLOCATE(this%ltmp0_ (this%partbuff_))
    ALLOCATE(this%ltmp1_ (this%partbuff_))

  END SUBROUTINE GPart_ctor
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Destructor
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE GPart_dtor(this)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  Main explicit destructor
!  ARGUMENTS:
!    this    : 'this' class instance
!-----------------------------------------------------------------

    IMPLICIT NONE
    TYPE(GPart)   ,INTENT(INOUT)             :: this
    INTEGER                                  :: j

    IF ( ALLOCATED    (this%id_) ) DEALLOCATE    (this%id_)
    IF ( ALLOCATED(this%tmpint_) ) DEALLOCATE(this%tmpint_)
    IF ( ALLOCATED   (this%idm_) ) DEALLOCATE   (this%idm_)
    IF ( ALLOCATED    (this%px_) ) DEALLOCATE    (this%px_)
    IF ( ALLOCATED    (this%py_) ) DEALLOCATE    (this%py_)
    IF ( ALLOCATED    (this%pz_) ) DEALLOCATE    (this%pz_)
    IF ( ALLOCATED   (this%lvx_) ) DEALLOCATE   (this%lvx_)
    IF ( ALLOCATED   (this%lvy_) ) DEALLOCATE   (this%lvy_)
    IF ( ALLOCATED   (this%lvz_) ) DEALLOCATE   (this%lvz_)
    IF ( ALLOCATED (this%ptmp0_) ) DEALLOCATE (this%ptmp0_)
    IF ( ALLOCATED (this%gptmp0_)) DEALLOCATE(this%gptmp0_)
    IF ( ALLOCATED (this%ptmp1_) ) DEALLOCATE (this%ptmp1_)
    IF ( ALLOCATED (this%ptmp2_) ) DEALLOCATE (this%ptmp2_)
    IF ( ALLOCATED   (this%vdb_) ) DEALLOCATE   (this%vdb_)
    IF ( ALLOCATED (this%ltmp0_) ) DEALLOCATE (this%ltmp0_)
    IF ( ALLOCATED (this%ltmp1_) ) DEALLOCATE (this%ltmp1_)
    IF ( ALLOCATED   (this%lvy_) ) DEALLOCATE   (this%lvy_)
    IF ( ALLOCATED   (this%lvz_) ) DEALLOCATE   (this%lvz_)
    IF ( ALLOCATED  (this%fpid_) ) DEALLOCATE  (this%fpid_)
    IF ( ALLOCATED   (this%vk0_) ) DEALLOCATE   (this%vk0_)
    IF ( ALLOCATED   (this%vk1_) ) DEALLOCATE   (this%vk1_)
    IF ( ALLOCATED   (this%vk2_) ) DEALLOCATE   (this%vk2_)
    IF ( ALLOCATED   (this%xk1_) ) DEALLOCATE   (this%xk1_)

    ! Destroy timers:
    DO j = 1, GPMAXTIMERS
      CALL GTFree(this%htimers_(j))
    ENDDO
  
  END SUBROUTINE GPart_dtor
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Convert input state name to index in state vector
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Get state variable names
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Function to compute number of state members (equations)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
