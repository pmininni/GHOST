! =====================================================================
! NAME       : inerpart_solver.f90
! DESCRIPTION: Forms solver class for inertial particles, computing:
!
!              dx/dt   = v_p
!              dv_p/dt = (1/tau) * (u(x,t) - v_p) - grav * z_hat
!                     or f(Re_p)/tau (u(x,t) - v_p)
!                     if nonlinear drag is used (see Clift & Gauvin 1971),
!                     where f(Re_p) = 1 + 0.15 Re_p^0.687
!                                   + (Re_p/57.14)/(1 + 4.25e4 Re_p^(-1.16)),
!                     Re_p = (18 tau gamma/nu)^(1/2) |u - v_p|, and
!                     |u - v_p| = [(ux-vpx)^2+(uy-vpy)^2+(uz-vpz)^2]^(1/2).
!
!              State ordering is:
!                x1 (x2, x3), v1 (v2, v3)
!
!              State sector ids are:
!                POSITION (POSITION+1, POSITION+2)
!                VELOCITY (VELOCITY+1, VELOCITY+2)
!
! INPUT FILE : For psolver='inerpart', looks for a "&inerpart" namelist with:
!              pidir   : changes class binary input  dir (default: status idir)
!              podir   : changes class binary output dir (default: status odir)
!              tau     : Stokes time (particle response time)
!              dograv  : .false.=no gravity [default], .true.=gravity
!              grav    : gravity acceleration (z-direction, positive
!                        downward, default=0)
!              donldrag: nonlinear drag (NLD) (.false.=linear [default],
!                        .true.=Clift & Gauvin 1971 correction).
!                        When true, nu is obtained from the PDE solver.
!              gamma   : mass ratio m_f/m_p (default=1, only for NLD)
!              partlod : particle output level of detail (default=1):
!                         1: position (xlg), part veloc (vip), fluid veloc (vlg)
!                         2: Lagrangian vorticity (wlg),
!                            strain-rate tensor (s11,s12,s13,s22,s23)
!
! DATE       : 04/02/26 (BLE)
! =====================================================================

module inerpart_mod
  use particlebase_mod
  use gpstate_mod

  implicit none

  ! ================= Solver traits ===================================
  type, public  :: InerTraits
    integer       :: partlod  = 1       ! particle output level of detail
    real(kind=GP) :: tau      = 1.0_GP  ! Stokes time
    real(kind=GP) :: grav     = 0.0_GP  ! gravity (z, positive downward)
    real(kind=GP) :: gamma    = 1.0_GP  ! mass ratio m_f/m_p
    real(kind=GP) :: nu       = 0.0_GP  ! fluid kinematic viscosity
    real(kind=GP) :: invtau   = 1.0_GP  ! precomputed 1/tau
    real(kind=GP) :: nld_rep  = 0.0_GP  ! precomputed NLD const (=18 tau gamma/nu)
    logical       :: dograv   = .false. ! .false.=no gravity, .true.=grav
    logical       :: donldrag = .false. ! .false.=linear drag, .true.=NLD
  end type

  ! ================= Solver ==========================================
  type, extends(VelocParticleBase) :: Ipart
    ! Member data:
    type (InerTraits) :: traits_
  CONTAINS
    procedure, public :: init          =>          init_impl
    procedure, public :: dpdt          =>          dpdt_impl
    procedure, public :: end_stage     =>     end_stage_impl
    procedure, public :: feedback      =>      null_feedback
    procedure, public :: write_pstate  =>  write_pstate_impl
    procedure, public :: state_size    =>    state_size_impl
    procedure, public :: part_ctor     =>         Ipart_ctor
    final             :: Ipart_dtor
  end type Ipart

CONTAINS

  ! ===================================================================
  ! Solver initialization, this is where parameter files are read
  ! ===================================================================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Subroutine to initialize the solver.
  !! Reads the &inerpart namelist and sets sector indices.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_impl(this)
    use commtypes
    class      (Ipart), intent (inout) :: this
    real     (kind=GP)                 :: tau, grav, gamma
    integer                            :: ierr, partlod
    logical                            :: dograv, donldrag
    character(len=128)                 :: pidir, podir
    namelist/ inerpart / pidir,podir,partlod,tau,grav,gamma,dograv,donldrag

    this%POSITION = 1
    this%VELOCITY = this%POSITION + this%nc_
    ! Defaults
    pidir    = this%idir_ ! Set to the pde class idir at ctor
    podir    = this%odir_ ! Set to the pde class odir at ctor
    tau      = 1.0_GP
    grav     = 0.0_GP
    gamma    = 1.0_GP
    dograv   = .false.
    donldrag = .false.
    partlod  = 1
    if ( this%myrank_ .eq. 0 ) then
      open(1,file=this%infile_,status='unknown',form="formatted")
      read(1,NML=inerpart)
      close(1)
    endif
    call MPI_BCAST(pidir   ,128,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(podir   ,128,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(tau     ,1  ,GC_REAL      ,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(grav    ,1  ,GC_REAL      ,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(gamma   ,1  ,GC_REAL      ,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(dograv  ,1  ,MPI_LOGICAL  ,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(donldrag,1  ,MPI_LOGICAL  ,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(partlod ,1  ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)

    this%idir_ = pidir ! If present in &inerpart, replaces the class default idir
    this%odir_ = podir ! If present in &inerpart, replaces the class default odir
    this%sstate_pos_ = 'xlg' ! state name of positions
    this%sstate_lag_ = 'vlg' ! state name of Lagrangian velocities
    this%sstate_vel_ = 'vip' ! state name of particles velocities
    this%traits_%     tau = tau
    this%traits_%    grav = grav
    this%traits_%   gamma = gamma
    this%traits_%  dograv = dograv
    this%traits_%donldrag = donldrag
    this%traits_%  invtau = 1.0_GP / tau
    this%traits_% partlod = partlod
  end subroutine init_impl

  ! ===================================================================
  ! Computation of RHS, the solver equations are defined here
  ! ===================================================================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Function to compute the rhs of the equations of motion
  !! of inertial particles:
  !!   dx/dt   = v_p
  !!   dv_p/dt = (1/tau)*(u - v_p) - grav * z_hat (or nonlinear-drag)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE dpdt_impl(this, time, pde, fluidstate, pstate, dt, dpdtout)
    use equationbase_mod
    use fft
!$  use threads
    IMPLICIT NONE
    class       (Ipart),             intent(inout) :: this
    class(EquationBase),             intent   (in) :: pde
    real      (kind=GP),             intent   (in) :: time, dt
    type   (GStateComp),             intent   (in) :: fluidstate(:)
    type  (GPStateComp), target ,    intent   (in) :: pstate(:)
    type  (GPStateComp),             intent(inout) :: dpdtout(:)
    complex   (KIND=GP), pointer, dimension(:,:,:) :: velc
    real      (KIND=GP), pointer, dimension(:,:,:) :: velr,tmp1,tmp2
    real      (kind=GP)                            :: rmp, invtau, grav
    real      (kind=GP)                            :: cdrag, rep2, dx, dy, dz
    real      (kind=GP)                            :: rep2_coef
    integer                                        :: i,j,k
    logical                                        :: bret

    CALL GTStart(this%htimers_(GPTIME_STEP))
    call this%workspace_%get_complex_tmp(velc,bret)
    call this%workspace_%get_real_tmp   (velr,bret)
    call this%workspace_%get_real_tmp   (tmp1,bret)
    call this%workspace_%get_real_tmp   (tmp2,bret)
    call this%AssignLagPos(pstate) ! We assign px_,py_,pz_ to the pstate

    invtau = this%traits_%invtau
    grav   = this%traits_%grav

    select type (pde)
    class is (VelocityBase)
      if (this%nc_ .ne. pde%nc_) then
        stop "Inerpart: # of components of the particles and pdes must be equal"
      endif

      ! IMPORTANT: pstate and dpdtout may alias the same array (some steppers
      ! can pass upout for both). We must be careful about ordering:
      !   1. Interpolate fluid velocity with positions still intact
      !   2. Write position RHS (overwrites pstate(POSITION), but we're done)
      !   3. Compute velocity RHS (reads particle velocity before overwrite)

      ! Step 1: Interpolate fluid velocity to lvx_, lvy_, lvz_
      rmp = 1.0_GP/(real(this%nd_(1),kind=GP)*real(this%nd_(2),kind=GP)* &
                    real(this%nd_(3),kind=GP))
      !$omp parallel do collapse(2) private (k)
      do i = ista,iend
        do j = 1,ny
          do concurrent (k=1:nz)
            velc(k,j,i) = fluidstate(pde%VELOCITY  )%ccomp(k,j,i)*rmp
          end do
        end do
      end do
      call fftp3d_complex_to_real(plancr,velc,velr,MPI_COMM_WORLD)
      call this%EulerToLag(this%lvx_,this%nparts_,velr,.true. ,tmp1,tmp2)
      !$omp parallel do collapse(2) private (k)
      do i = ista,iend
        do j = 1,ny
          do concurrent (k=1:nz)
            velc(k,j,i) = fluidstate(pde%VELOCITY+1)%ccomp(k,j,i)*rmp
          end do
        end do
      end do
      call fftp3d_complex_to_real(plancr,velc,velr,MPI_COMM_WORLD)
      call this%EulerToLag(this%lvy_,this%nparts_,velr,.false.,tmp1,tmp2)
      !$omp parallel do collapse(2) private (k)
      do i = ista,iend
        do j = 1,ny
          do concurrent (k=1:nz)
            velc(k,j,i) = fluidstate(pde%VELOCITY+2)%ccomp(k,j,i)*rmp
          end do
        end do
      end do
      call fftp3d_complex_to_real(plancr,velc,velr,MPI_COMM_WORLD)
      call this%EulerToLag(this%lvz_,this%nparts_,velr,.false.,tmp1,tmp2)

      rep2_coef = this%traits_%nld_rep
      ! One parallel region for the position, drag and gravity loops:
      ! the team is forked once, and every thread takes the same IF branches
!$omp parallel private(dx,dy,dz,rep2,cdrag)
      ! Step 2: Position RHS: dx/dt = v_p
      ! (may overwrite pstate(POSITION) depending on the stepper call)
!$omp do
      do j = 1,this%nparts_
        dpdtout(this%POSITION  )%rcomp(j)=pstate(this%VELOCITY  )%rcomp(j)*this%invdel_(1)
        dpdtout(this%POSITION+1)%rcomp(j)=pstate(this%VELOCITY+1)%rcomp(j)*this%invdel_(2)
        dpdtout(this%POSITION+2)%rcomp(j)=pstate(this%VELOCITY+2)%rcomp(j)*this%invdel_(3)
      end do

      ! Step 3: Velocity RHS: dv_p/dt = cdrag/tau * (u - v_p)
      ! (reads pstate(VELOCITY) before it may be overwritten by some steppers)
      if ( .not. this%traits_%donldrag ) then
        ! Linear Stokes drag
!$omp do
        do j = 1,this%nparts_
          dpdtout(this%VELOCITY  )%rcomp(j) = &
            (this%lvx_(j) - pstate(this%VELOCITY  )%rcomp(j)) * invtau
          dpdtout(this%VELOCITY+1)%rcomp(j) = &
            (this%lvy_(j) - pstate(this%VELOCITY+1)%rcomp(j)) * invtau
          dpdtout(this%VELOCITY+2)%rcomp(j) = &
            (this%lvz_(j) - pstate(this%VELOCITY+2)%rcomp(j)) * invtau
        end do
      else
        ! Nonlinear drag (Clift & Gauvin 1971):
        !   cdrag = 1 + 0.15 * Re_p^0.687
        !             + Re_p^2.16 / (57.14 * (Re_p^1.16 + 4.25e4))
        !   where Re_p^2 = rep2_coef * |u - v_p|^2,
        !         rep2_coef = 18 tau gamma / nu
!$omp do
        do j = 1,this%nparts_
          dx = this%lvx_(j) - pstate(this%VELOCITY  )%rcomp(j)
          dy = this%lvy_(j) - pstate(this%VELOCITY+1)%rcomp(j)
          dz = this%lvz_(j) - pstate(this%VELOCITY+2)%rcomp(j)
          rep2 = rep2_coef * (dx**2 + dy**2 + dz**2)
          cdrag = 1.0_GP + 0.15_GP * rep2 ** 0.3435_GP                  &
                + rep2 ** 1.08_GP / (57.14_GP * (rep2 ** 0.58_GP + 4.25e4_GP))
          dpdtout(this%VELOCITY  )%rcomp(j) = dx * cdrag * invtau
          dpdtout(this%VELOCITY+1)%rcomp(j) = dy * cdrag * invtau
          dpdtout(this%VELOCITY+2)%rcomp(j) = dz * cdrag * invtau
        end do
      endif

      ! Add gravity in z-direction
      if ( this%traits_%dograv ) then
!$omp do
        do j = 1,this%nparts_
          dpdtout(this%VELOCITY+2)%rcomp(j) = &
            dpdtout(this%VELOCITY+2)%rcomp(j) - grav
        end do
      endif
!$omp end parallel
    class default
      stop "Inerpart: This solver does not support pdes without a velocity field"
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
    class     (Ipart), intent   (in) :: this
    type(GPStateComp), intent   (in) :: pstate(:)
    type (GStateComp), intent(inout) :: feedback(:)
    return
  end subroutine null_feedback


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Functions to sync particles after doing a time step.
  !! Syncs only initial and final states, does not sync
  !! tmp arrays that may be used by the stepper.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine end_stage_impl(this, upin, upout)
    use gpstate_mod
    implicit none
    class     (IPart), intent(inout) :: this
    type(GPStateComp), intent(inout) :: upin (:) ! state at t0
    type(GPStateComp), intent(inout) :: upout(:) ! state after sub-stage
    integer                          :: j, ng

    ! Branch 1: Nearest-neighbour (NN) exchange -------------------------
    if (this%iexchtype_ .EQ. GPEXCHTYPE_NN) then
      ! We first enforce periodicity in x-y only
      CALL this%MakePeriodicP(upout(this%POSITION  )%rcomp,            &
                              upout(this%POSITION+1)%rcomp,            &
                              upout(this%POSITION+2)%rcomp, this%nparts_, 3)
      ! We identify particles that have left the local slab in z
      CALL GTStart(this%htimers_(GPTIME_COMM))
      CALL this%gpcomm_%IdentifyExchV(this%id_,                        &
               upout(this%POSITION+2)%rcomp,                           &
               this%nparts_, ng, this%lxbnds_(3,1),  this%lxbnds_(3,2))
      ! We resize internal buffers if the exchange set is too large
      if (ng .GT. this%partbuff_) then
        WRITE(*,'(A,I0,A,I0,A,I0,A,I0)') &
          'EndStage: Rank ', this%myrank_, ' resizing: nparts=', ng,   &
          ' | partbuff=', this%partbuff_, ' --> ', this%partbuff_ +    &
          (1 + (ng - this%partbuff_) / this%partchunksize_) * this%partchunksize_
        this%partbuff_ = this%partbuff_ + (1 + (ng - this%partbuff_) / &
               this%partchunksize_) * this%partchunksize_
        call this%ResizeArrays   (this%partbuff_,.true.)
        call GPState_resize(upin ,this%partbuff_)
        call GPState_resize(upout,this%partbuff_)
      end if
      ! Exchange current positions (upout) across MPI tasks
      CALL this%gpcomm_%PartExchangeV(this%id_,                        &
               upout(this%POSITION  )%rcomp,                           &
               upout(this%POSITION+1)%rcomp,                           &
               upout(this%POSITION+2)%rcomp,                           &
               this%nparts_,this%lxbnds_(3,1),this%lxbnds_(3,2), GPEXCH_INIT)
      ! Exchange the current velocity (upout)
      CALL this%gpcomm_%PartExchangeV(this%id_,                        &
               upout(this%VELOCITY  )%rcomp,                           &
               upout(this%VELOCITY+1)%rcomp,                           &
               upout(this%VELOCITY+2)%rcomp,                           &
               this%nparts_,this%lxbnds_(3,1),this%lxbnds_(3,2), GPEXCH_UPDT)
      ! Exchange the previous positions (upin)
      CALL this%gpcomm_%PartExchangeV(this%id_,                        &
               upin (this%POSITION  )%rcomp,                           &
               upin (this%POSITION+1)%rcomp,                           &
               upin (this%POSITION+2)%rcomp,                           &
               this%nparts_,this%lxbnds_(3,1),this%lxbnds_(3,2), GPEXCH_UPDT)
      ! Exchange the previous velocity (upin)
      CALL this%gpcomm_%PartExchangeV(this%id_,                        &
               upin (this%VELOCITY  )%rcomp,                           &
               upin (this%VELOCITY+1)%rcomp,                           &
               upin (this%VELOCITY+2)%rcomp,                           &
               this%nparts_,this%lxbnds_(3,1),this%lxbnds_(3,2), GPEXCH_END )
      CALL GTAcc(this%htimers_(GPTIME_COMM))
      ! x-y periodicity already enforced, we enforce z in upout and upin
      CALL this%MakePeriodicZ(upout(this%POSITION+2)%rcomp,            &
                              upin (this%POSITION+2)%rcomp, this%nparts_)
      ! Buffer shrink
      if (this%stepcounter_ .GE. GPSWIPERATE) then
        if ((this%bcollective_ .EQ. 1) .OR. (this%myrank_ .NE. 0)) then
          ng = this%partbuff_ - this%nparts_
          ng = this%partbuff_ - (ng/this%partchunksize_-1)*this%partchunksize_
          if (ng .LT. this%partbuff_) then
            WRITE(*,'(A,I0,A,I0,A,I0,A,I0)') 'EndStage: Rank ',        &
              this%myrank_, ' shrinking: nparts=', this%nparts_,       &
              ' | partbuff=', this%partbuff_, ' --> ', ng
            this%partbuff_ = ng
            call this%ResizeArrays   (this%partbuff_,.false.)
            call GPState_resize(upin ,this%partbuff_)
            call GPState_resize(upout,this%partbuff_)
          end if
        end if
        this%stepcounter_ = 1
      else
        this%stepcounter_ = this%stepcounter_ + 1
      end if
    end if  ! GPEXCHTYPE_NN
    ! Branch 2: Voxel Database (VDB) exchange ---------------------------
    if (this%iexchtype_ .EQ. GPEXCHTYPE_VDB) then
      ! Enforce x-y-z periodicity on updated positions
      CALL this%MakePeriodicP(upout(this%POSITION  )%rcomp,            &
                              upout(this%POSITION+1)%rcomp,            &
                              upout(this%POSITION+2)%rcomp, this%nparts_, 7)
      ! Consistency check
      if (.NOT. this%PartNumConsistent(this%nparts_)) then
        if (this%myrank_ .EQ. 0) then
          WRITE(*,*) 'Inerpart EndStage (VDB): inconsistent particle count'
          print *,this%nparts_,this%maxparts_
        end if
      end if
      ! Sync global VDB for the current positions (upout)
      CALL GTStart(this%htimers_(GPTIME_COMM))
      CALL this%gpcomm_%VDBSynch(this%vdb_,this%maxparts_,this%id_,    &
               upout(this%POSITION  )%rcomp,                           &
               upout(this%POSITION+1)%rcomp,                           &
               upout(this%POSITION+2)%rcomp,                           &
               this%nparts_,this%ptmp0_)
      ! Sync current velocities (unpout) into gptmp0_ and extract locally
      CALL this%gpcomm_%VDBSynch(this%gptmp0_,this%maxparts_,this%id_, &
               upout(this%VELOCITY  )%rcomp,                           &
               upout(this%VELOCITY+1)%rcomp,                           &
               upout(this%VELOCITY+2)%rcomp,                           &
               this%nparts_,this%ptmp0_)
      CALL this%CopyLocalWrk(upout(this%VELOCITY  )%rcomp,             &
                             upout(this%VELOCITY+1)%rcomp,             &
                             upout(this%VELOCITY+2)%rcomp,             &
                             this%vdb_,this%gptmp0_,this%maxparts_)
      ! Sync previous velocities (upin) into gptmp0_ and extract locally
      CALL this%gpcomm_%VDBSynch(this%gptmp0_,this%maxparts_,this%id_, &
               upin (this%VELOCITY  )%rcomp,                           &
               upin (this%VELOCITY+1)%rcomp,                           &
               upin (this%VELOCITY+2)%rcomp,                           &
               this%nparts_,this%ptmp0_)
      CALL this%CopyLocalWrk(upin (this%VELOCITY  )%rcomp,             &
                             upin (this%VELOCITY+1)%rcomp,             &
                             upin (this%VELOCITY+2)%rcomp,             &
                             this%vdb_,this%gptmp0_,this%maxparts_)
      ! Sync gptmp VDB for the previous positions (upin)
      CALL this%gpcomm_%VDBSynch(this%gptmp0_,this%maxparts_,this%id_, &
               upin (this%POSITION  )%rcomp,                           &
               upin (this%POSITION+1)%rcomp,                           &
               upin (this%POSITION+2)%rcomp,                           &
               this%nparts_,this%ptmp0_)
      CALL GTAcc(this%htimers_(GPTIME_COMM))
      ! Get local particles (upout and upin) based on positions
      CALL this%GetLocalWrk_aux(this%id_,                              &
               upout(this%POSITION  )%rcomp,                           &
               upout(this%POSITION+1)%rcomp,                           &
               upout(this%POSITION+2)%rcomp,                           &
               upin (this%POSITION  )%rcomp,                           &
               upin (this%POSITION+1)%rcomp,                           &
               upin (this%POSITION+2)%rcomp,                           &
               this%nparts_,this%vdb_,this%gptmp0_,this%maxparts_)
      ! Global particle-count sanity check
      CALL MPI_ALLREDUCE(this%nparts_, ng, 1, MPI_INTEGER,             &
                         MPI_SUM, this%comm_, this%ierr_)
      if (this%myrank_ .EQ. 0 .AND. ng .NE. this%maxparts_) then
        WRITE(*,*) 'Inerpart EndStage (VDB): inconsistent d.b.: expected: ', &
                   this%maxparts_, '; found: ', ng
        CALL this%ascii_write_lag(1, '.', trim(this%sstate_pos_) // 'err',   &
             '000', 0.0_GP, this%maxparts_,this%vdb_)
        STOP
      end if
    end if  ! GPEXCHTYPE_VDB
  end subroutine end_stage_impl


  ! ===================================================================
  ! Output methods
  ! ===================================================================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Function to compute and write particle states.
  !! Writes: positions (xlg), fluid velocity (vlg),
  !!         inertial particle velocity (vip)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine write_pstate_impl(this, time, pde, fluidstate, pstate)
    use equationbase_mod
    use particlebase_mod
    use pseudospec_fluid
    use status
    use pstatus
    use fft
!$  use threads
    class       (IPart),             intent(inout) :: this
    class(EquationBase),             intent   (in) :: pde
    type   (GStateComp),             intent   (in) :: fluidstate(:)
    type  (GPStateComp), target ,    intent   (in) :: pstate(:)
    real      (kind=GP),             intent   (in) :: time
    complex   (kind=GP), pointer, dimension(:,:,:) :: velc, velc2
    real      (kind=GP), pointer, dimension(:,:,:) :: velr,tmp1,tmp2
    real      (kind=GP)                            :: rmp
    integer                                        :: i,j,k
    logical                                        :: bret

    call this%workspace_%get_complex_tmp(velc,bret)
    call this%workspace_%get_real_tmp   (velr,bret)
    call this%workspace_%get_real_tmp   (tmp1,bret)
    call this%workspace_%get_real_tmp   (tmp2,bret)
    call this%AssignLagPos(pstate)

    select type (pde)
    class is (VelocityBase)
      rmp = 1.0_GP/(real(this%nd_(1),kind=GP)*real(this%nd_(2),kind=GP)* &
                    real(this%nd_(3),kind=GP))
      ! Interpolate fluid velocity to particle positions
!$omp parallel do collapse(2) private (k)
      DO i = ista,iend
        DO j = 1,ny
          DO CONCURRENT (k=1:nz)
            velc(k,j,i) = fluidstate(pde%VELOCITY  )%ccomp(k,j,i)*rmp
          END DO
        END DO
      END DO
      CALL fftp3d_complex_to_real(plancr,velc,velr,MPI_COMM_WORLD)
      CALL this%EulerToLag(this%lvx_,this%nparts_,velr,.true. ,tmp1,tmp2)
!$omp parallel do collapse(2) private (k)
      DO i = ista,iend
        DO j = 1,ny
          DO CONCURRENT (k=1:nz)
            velc(k,j,i) = fluidstate(pde%VELOCITY+1)%ccomp(k,j,i)*rmp
          END DO
        END DO
      END DO
      CALL fftp3d_complex_to_real(plancr,velc,velr,MPI_COMM_WORLD)
      CALL this%EulerToLag(this%lvy_,this%nparts_,velr,.false.,tmp1,tmp2)
!$omp parallel do collapse(2) private (k)
      DO i = ista,iend
        DO j = 1,ny
          DO CONCURRENT (k=1:nz)
            velc(k,j,i) = fluidstate(pde%VELOCITY+2)%ccomp(k,j,i)*rmp
          END DO
        END DO
      END DO
      CALL fftp3d_complex_to_real(plancr,velc,velr,MPI_COMM_WORLD)
      CALL this%EulerToLag(this%lvz_,this%nparts_,velr,.false.,tmp1,tmp2)
      ! Write positions and Lagrangian fluid velocity
      WRITE(lgext,lgfmtext) pind
      CALL this%io_write_pdb(1,this%odir_,trim(this%sstate_pos_),lgext,time)
      CALL this%io_write_vec(1,this%odir_,trim(this%sstate_lag_),lgext,time)
      ! Write inertial particle velocity
      do j = 1,this%nparts_
        this%lvx_(j) = pstate(this%VELOCITY  )%rcomp(j)
        this%lvy_(j) = pstate(this%VELOCITY+1)%rcomp(j)
        this%lvz_(j) = pstate(this%VELOCITY+2)%rcomp(j)
      end do
      CALL this%io_write_vec(1,this%odir_,trim(this%sstate_vel_),lgext,time)
! partlod >= 2: write Lagrangian vorticity and strain-rate tensor
      if ( this%traits_%partlod .ge. 2 ) then
! Write Lagrangian vorticity components
        CALL rotor3(fluidstate(pde%VELOCITY+1)%ccomp, &
                    fluidstate(pde%VELOCITY+2)%ccomp, velc, 1)
!$omp parallel do collapse(2) private (k)
        DO i = ista,iend
          DO j = 1,ny
            DO CONCURRENT (k=1:nz)
              velc(k,j,i) = velc(k,j,i)*rmp
            END DO
          END DO
        END DO
        CALL fftp3d_complex_to_real(plancr,velc,velr,MPI_COMM_WORLD)
        CALL this%EulerToLag(this%lvx_,this%nparts_,velr,.false.,tmp1,tmp2)
        CALL rotor3(fluidstate(pde%VELOCITY  )%ccomp, &
                    fluidstate(pde%VELOCITY+2)%ccomp, velc, 2)
!$omp parallel do collapse(2) private (k)
        DO i = ista,iend
          DO j = 1,ny
            DO CONCURRENT (k=1:nz)
              velc(k,j,i) = velc(k,j,i)*rmp
            END DO
          END DO
        END DO
        CALL fftp3d_complex_to_real(plancr,velc,velr,MPI_COMM_WORLD)
        CALL this%EulerToLag(this%lvy_,this%nparts_,velr,.false.,tmp1,tmp2)
        CALL rotor3(fluidstate(pde%VELOCITY  )%ccomp, &
                    fluidstate(pde%VELOCITY+1)%ccomp, velc, 3)
!$omp parallel do collapse(2) private (k)
        DO i = ista,iend
          DO j = 1,ny
            DO CONCURRENT (k=1:nz)
              velc(k,j,i) = velc(k,j,i)*rmp
            END DO
          END DO
        END DO
        CALL fftp3d_complex_to_real(plancr,velc,velr,MPI_COMM_WORLD)
        CALL this%EulerToLag(this%lvz_,this%nparts_,velr,.false.,tmp1,tmp2)
        CALL this%io_write_vec(1,this%odir_,'wlg',lgext,time)
! Write strain-rate tensor components
        call this%workspace_%get_complex_tmp(velc2,bret)
        ! S11 = dv_x/dx
        CALL derivk3(fluidstate(pde%VELOCITY  )%ccomp, velc, 1)
!$omp parallel do collapse(2) private (k)
        DO i = ista,iend
          DO j = 1,ny
            DO CONCURRENT (k=1:nz)
              velc(k,j,i) = velc(k,j,i)*rmp
            END DO
          END DO
        END DO
        CALL fftp3d_complex_to_real(plancr,velc,velr,MPI_COMM_WORLD)
        CALL this%io_write_euler(1,this%odir_,'s11',lgext,time,velr,.false.,tmp1,tmp2)
        ! S12 = 0.5*(dv_x/dy + dv_y/dx)
        CALL derivk3(fluidstate(pde%VELOCITY  )%ccomp, velc,  2)
        CALL derivk3(fluidstate(pde%VELOCITY+1)%ccomp, velc2, 1)
!$omp parallel do collapse(2) private (k)
        DO i = ista,iend
          DO j = 1,ny
            DO CONCURRENT (k=1:nz)
              velc(k,j,i) = 0.5_GP*(velc(k,j,i)+velc2(k,j,i))*rmp
            END DO
          END DO
        END DO
        CALL fftp3d_complex_to_real(plancr,velc,velr,MPI_COMM_WORLD)
        CALL this%io_write_euler(1,this%odir_,'s12',lgext,time,velr,.false.,tmp1,tmp2)
        ! S13 = 0.5*(dv_x/dz + dv_z/dx)
        CALL derivk3(fluidstate(pde%VELOCITY  )%ccomp, velc,  3)
        CALL derivk3(fluidstate(pde%VELOCITY+2)%ccomp, velc2, 1)
!$omp parallel do collapse(2) private (k)
        DO i = ista,iend
          DO j = 1,ny
            DO CONCURRENT (k=1:nz)
              velc(k,j,i) = 0.5_GP*(velc(k,j,i)+velc2(k,j,i))*rmp
            END DO
          END DO
        END DO
        CALL fftp3d_complex_to_real(plancr,velc,velr,MPI_COMM_WORLD)
        CALL this%io_write_euler(1,this%odir_,'s13',lgext,time,velr,.false.,tmp1,tmp2)
        ! S22 = dv_y/dy
        CALL derivk3(fluidstate(pde%VELOCITY+1)%ccomp, velc, 2)
!$omp parallel do collapse(2) private (k)
        DO i = ista,iend
          DO j = 1,ny
            DO CONCURRENT (k=1:nz)
              velc(k,j,i) = velc(k,j,i)*rmp
            END DO
          END DO
        END DO
        CALL fftp3d_complex_to_real(plancr,velc,velr,MPI_COMM_WORLD)
        CALL this%io_write_euler(1,this%odir_,'s22',lgext,time,velr,.false.,tmp1,tmp2)
        ! S23 = 0.5*(dv_y/dz + dv_z/dy)
        CALL derivk3(fluidstate(pde%VELOCITY+1)%ccomp, velc,  3)
        CALL derivk3(fluidstate(pde%VELOCITY+2)%ccomp, velc2, 2)
!$omp parallel do collapse(2) private (k)
        DO i = ista,iend
          DO j = 1,ny
            DO CONCURRENT (k=1:nz)
              velc(k,j,i) = 0.5_GP*(velc(k,j,i)+velc2(k,j,i))*rmp
            END DO
          END DO
        END DO
        CALL fftp3d_complex_to_real(plancr,velc,velr,MPI_COMM_WORLD)
        CALL this%io_write_euler(1,this%odir_,'s23',lgext,time,velr,.false.,tmp1,tmp2)
        call this%workspace_%free_complex_tmp(velc2)
      endif
    class default
      stop "inerpart: This solver does not support pdes without a velocity field"
    end select

    call this%workspace_%free_complex_tmp(velc)
    call this%workspace_%free_real_tmp   (velr)
    call this%workspace_%free_real_tmp   (tmp1)
    call this%workspace_%free_real_tmp   (tmp2)
  end subroutine write_pstate_impl


  ! ===================================================================
  ! Solver specific methods
  ! ===================================================================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Constructor
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE Ipart_ctor(this,infile, pde, workspace, pstate, pstate_cpy)
    USE equationbase_mod
    USE hd_mod,    ONLY: HDSolver
    USE mhd_mod,   ONLY: MHDSolver
    USE moist_mod, ONLY: MOISTSolver
    USE bouss_mod, ONLY: BOUSSSolver
    USE var
    USE grid
    USE boxsize
    USE mpivars
    USE pstatus
    USE commtypes
    USE fftplans
    USE pstatus
    USE status
    USE random
    IMPLICIT NONE
    CLASS       (IPart), intent(inout)              :: this
    class(EquationBase), intent   (in)              :: pde
    type   (GWorkspace), intent(inout), target      :: workspace
    type  (GPStateComp), intent(inout), allocatable :: pstate(:), pstate_cpy(:)
    character   (len=*), intent   (in)              :: infile
    integer                                         :: disp(3),lens(3),types(3)
    integer                                         :: tsta,tend,num_components
    integer                                         :: j,nc,szreal
    logical                                         :: bret

    this%infile_      =  infile
    this%idir_        =  pde%idir_ ! input  directory, same as in the pde class
    this%odir_        =  pde%odir_ ! output directory, same as in the pde class
    this%workspace_   => workspace
    this%hasfeedback_ = .false.
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
    CALL this%SetRandSeed(seed)
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
        WRITE(*,*) 'Inerpart_ctor: Not enough timers available'
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
    this%lxbnds_(3,1) = real(ksta-1,kind=GP)
    this%lxbnds_(3,2) = real(kend-1,kind=GP) + 1.0_GP
    CALL range(1,nx,nprocs,myrank,tsta,tend)
    this%tibnds_(1,1) = 1
    this%tibnds_(1,2) = nz
    this%tibnds_(2,1) = 1
    this%tibnds_(2,2) = ny
    this%tibnds_(3,1) = tsta
    this%tibnds_(3,2) = tend
    DO j = 1,3
      this%gext_ (j) = real(this%nd_(j),kind=GP)
    ENDDO

    ! Call init (reads &inerpart namelist, sets POSITION/VELOCITY indices)
    call this%init()

    ! If nonlinear drag is enabled, get nu from the PDE solver
    if ( this%traits_%donldrag ) then
      select type (pde)
        type is (HDSolver)
          this%traits_%nu = pde%traits_%nu
        type is (MHDSolver)
          this%traits_%nu = pde%traits_%nu
        type is (BOUSSSolver)
          this%traits_%nu = pde%traits_%nu
        type is (MOISTSolver)
          this%traits_%nu = pde%traits_%nu
        class default
          stop "Inerpart_ctor: nonlinear drag not supported by this PDE solver"
      end select
      this%traits_%nld_rep = 18.0_GP * this%traits_%tau * &
        this%traits_%gamma / this%traits_%nu
    endif

    ! Instantiate interp operation
    CALL this%intop_%GPSplineInt_ctor(3,this%nd_,this%libnds_,this%lxbnds_, &
         this%tibnds_,this%intorder_,this%partbuff_,this%gpcomm_,&
         this%htimers_(GPTIME_DATAEX),this%htimers_(GPTIME_TRANSP))

    ! Allocate particle arrays
    CALL MPI_TYPE_SIZE(GC_REAL,szreal,this%ierr_)
    ALLOCATE(this%id_      (this%partbuff_))
    ALLOCATE(this%tmpint_  (this%partbuff_))
    ALLOCATE(this%ptmp0_ (3,this%partbuff_))
    ALLOCATE(this%gptmp0_(3,this%partbuff_))
    num_components = this%state_size()  ! Returns 6
    CALL GPState_alloc(pstate    , num_components, this%partbuff_)
    CALL GPState_alloc(pstate_cpy, num_components, this%partbuff_)
    call this%workspace_%set_nparts(this%partbuff_)
    call this%workspace_%init_pcomp_arrays(this%partbuff_)
    call this%workspace_%get_pcomp_tmp(this%lvx_,bret)
    call this%workspace_%get_pcomp_tmp(this%lvy_,bret)
    call this%workspace_%get_pcomp_tmp(this%lvz_,bret)
    IF ( this%iexchtype_.EQ.GPEXCHTYPE_VDB ) THEN
      ALLOCATE(this%vdb_   (3,this%partbuff_))
    ENDIF
  END SUBROUTINE Ipart_ctor


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Destructor
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE Ipart_dtor(this)
    IMPLICIT NONE
    TYPE(Ipart),INTENT(INOUT) :: this
    integer                   :: j
    IF ( ALLOCATED    (this%id_) ) DEALLOCATE    (this%id_)
    IF ( ALLOCATED(this%tmpint_) ) DEALLOCATE(this%tmpint_)
    IF ( ALLOCATED   (this%idm_) ) DEALLOCATE   (this%idm_)
    IF ( ALLOCATED   (this%vdb_) ) DEALLOCATE   (this%vdb_)
    IF ( ALLOCATED (this%ptmp0_) ) DEALLOCATE (this%ptmp0_)
    IF ( ALLOCATED(this%gptmp0_) ) DEALLOCATE(this%gptmp0_)
    IF ( ASSOCIATED   (this%px_) ) NULLIFY       (this%px_)
    IF ( ASSOCIATED   (this%py_) ) NULLIFY       (this%py_)
    IF ( ASSOCIATED   (this%pz_) ) NULLIFY       (this%pz_)
    IF ( ASSOCIATED(this%lvx_) ) CALL this%workspace_%free_pcomp_tmp(this%lvx_)
    IF ( ASSOCIATED(this%lvy_) ) CALL this%workspace_%free_pcomp_tmp(this%lvy_)
    IF ( ASSOCIATED(this%lvz_) ) CALL this%workspace_%free_pcomp_tmp(this%lvz_)
    DO j = 1, GPMAXTIMERS
      CALL GTFree(this%htimers_(j))
    ENDDO
  END SUBROUTINE Ipart_dtor


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Function to compute number of state members (equations)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PURE function state_size_impl(this) result(num)
    class(Ipart), intent(in) :: this
    integer                  :: num
    num = 2 * this%nc_          ! 3 positions + 3 velocities
  end function state_size_impl

end module inerpart_mod
