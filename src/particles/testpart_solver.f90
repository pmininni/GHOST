! =====================================================================
! NAME       : testpart_solver.f90
! DESCRIPTION: Forms solver class for test particles (charged particles
!              that do not act back on the fields), computing:
!
!              dx/dt   = v_p
!              dv_p/dt = gyrof * [ E(x(t)) + v_p x B(x(t)) ]
!                      = gyrof * [ (v_p - u_e) x B + eta j ]
!
!                     where gyrof = q/m is the gyrofrequency for a unit
!                     magnetic field, B = curl(a) + B_0 is the total
!                     magnetic field, j = curl(B) the current density,
!                     eta the magnetic diffusivity, and E the electric
!                     field from Ohm's law, E = -u_e x B + eta j. The
!                     velocity u_e is the fluid velocity u, or the
!                     electron velocity u_e = u - dii j if dokinelv is
!                     .TRUE. (Hall-MHD). The magnetic diffusivity and the
!                     guide field are those of the MHD solver.
!
!              State ordering is:
!                x1 (x2, x3), v1 (v2, v3)
!
!              State sector ids are:
!                POSITION (POSITION+1, POSITION+2)
!                VELOCITY (VELOCITY+1, VELOCITY+2)
!
! INPUT FILE : For psolver='testpart', looks for a "&testpart" namelist with:
!              pidir   : changes class binary input  dir (default: status idir)
!              podir   : changes class binary output dir (default: status odir)
!              gyrof   : gyrofrequency for a unit magnetic field (q/m,
!                        default=1)
!              dii     : ion inertial length scale (Hall-MHD epsilon,
!                        default=0, only used if dokinelv is .TRUE.)
!              dokinelv: .false.=E computed with the fluid velocity [default],
!                        .true.=E computed with the electron velocity
!                        u_e = u - dii j (Hall-MHD correction)
!              partlod : particle output level of detail (default=1):
!                         1: position (xlg), fluid velocity (vlg), test
!                            particle velocity (vtp), magnetic field (blg)
!                            and current density (jlg) at the particles
!                         2: Lagrangian vorticity (wlg),
!                            strain-rate tensor (s11,s12,s13,s22,s23)
!              The initial velocity of the particles is set by the
!              particle initial conditions (e.g., 'thermal_v' with a
!              "&thermal_v" namelist and the thermal speed vtherm).
!
! NOTE       : For compressible MHD solvers (not available yet) the old
!              code also supported the electron pressure correction to
!              the electric field (ambipolar diffusion), selected with
!              the namelist flag dokinelp:
!                E = -u_e x B + eta j - (dii/2) grad(p)/rho
!              with u_e = u - dii j/rho. The places where dokinelp and
!              the density enter are marked with "COMPRESSIBLE" comments
!              in this file: the traits and the namelist (init_impl),
!              the electron velocity (dpdt_impl) and the current density
!              term (tpart_current), and the selection of the solver
!              traits in the constructor (Tpart_ctor).
!
! DATE       : 09/10/26 (PDM)
! =====================================================================

module testpart_mod
  use particlebase_mod
  use gpstate_mod
  use pseudospec_fluid, only: copy3, scal3, saxpby_c, rotor3, derivk3, &
                              laplak3, setmode3
  use gmem
  use gdevice, only: gdev_active

  implicit none

  ! ================= Solver traits ===================================
  type, public  :: TestTraits
    integer       :: partlod  = 1       ! particle output level of detail
    real(kind=GP) :: gyrof    = 1.0_GP  ! gyrofrequency (q/m)
    real(kind=GP) :: dii      = 0.0_GP  ! ion inertial length (Hall epsilon)
    real(kind=GP) :: eta      = 0.0_GP  ! magnetic diffusivity (from the pde)
    real(kind=GP) :: gyeta    = 0.0_GP  ! precomputed gyrof*eta
    real(kind=GP) :: B0(3)    = 0.0_GP  ! guide field (from the pde)
    logical       :: doB0     = .false. ! guide field flag (from the pde)
    logical       :: dokinelv = .false. ! .true.=electron velocity in E
    ! COMPRESSIBLE: add here the traits of compressible MHD solvers,
    ! e.g., dokinelp (.true.=electron pressure correction to E) and the
    ! equation of state parameters needed to compute grad(p)/rho (in the
    ! old code, cp1 = dii*betae/2 and gam1 = gamma-1 from the solver).
  end type

  ! ================= Solver ==========================================
  type, extends(VelocParticleBase) :: Tpart
    ! Member data:
    type (TestTraits) :: traits_
  CONTAINS
    procedure, public :: init          =>          init_impl
    procedure, public :: dpdt          =>          dpdt_impl
    procedure, public :: end_stage     =>     end_stage_impl
    procedure, public :: feedback      =>      null_feedback
    procedure, public :: write_pstate  =>  write_pstate_impl
    procedure, public :: state_size    =>    state_size_impl
    procedure, public :: part_ctor     =>         Tpart_ctor
    final             :: Tpart_dtor
  end type Tpart

CONTAINS

  ! ===================================================================
  ! Solver initialization, this is where parameter files are read
  ! ===================================================================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Subroutine to initialize the solver.
  !! Reads the &testpart namelist and sets sector indices.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_impl(this)
    use commtypes
    class      (Tpart), intent (inout) :: this
    real     (kind=GP)                 :: gyrof, dii
    integer                            :: ierr, partlod
    logical                            :: dokinelv
    character(len=128)                 :: pidir, podir
    ! COMPRESSIBLE: add dokinelp (logical, default .false.) to the
    ! namelist, its default, MPI_BCAST and trait below, next to dokinelv.
    ! It must be rejected (or ignored with a warning) when the pde is
    ! not compressible, since grad(p)/rho is not available then.
    namelist/ testpart / pidir,podir,partlod,gyrof,dii,dokinelv

    this%POSITION = 1
    this%VELOCITY = this%POSITION + this%nc_
    ! Defaults
    pidir    = this%idir_ ! Set to the pde class idir at ctor
    podir    = this%odir_ ! Set to the pde class odir at ctor
    partlod  = 1
    gyrof    = 1.0_GP
    dii      = 0.0_GP
    dokinelv = .false.
    if ( this%myrank_ .eq. 0 ) then
      open(1,file=this%infile_,status='unknown',form="formatted")
      read(1,NML=testpart)
      close(1)
    endif
    call MPI_BCAST(pidir   ,128,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(podir   ,128,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(partlod ,1  ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(gyrof   ,1  ,GC_REAL      ,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(dii     ,1  ,GC_REAL      ,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(dokinelv,1  ,MPI_LOGICAL  ,0,MPI_COMM_WORLD,ierr)

    this%idir_ = pidir ! If present in &testpart, replaces the class default idir
    this%odir_ = podir ! If present in &testpart, replaces the class default odir
    this%sstate_pos_ = 'xlg' ! state name of positions
    this%sstate_lag_ = 'vlg' ! state name of Lagrangian velocities
    this%sstate_vel_ = 'vtp' ! state name of particles velocities
    this%traits_% partlod = partlod
    this%traits_%   gyrof = gyrof
    this%traits_%     dii = dii
    this%traits_%dokinelv = dokinelv
  end subroutine init_impl

  ! ===================================================================
  ! Computation of RHS, the solver equations are defined here
  ! ===================================================================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Function to compute the rhs of the equations of motion
  !! of test particles:
  !!   dx/dt   = v_p
  !!   dv_p/dt = gyrof [ (v_p - u_e) x B + eta j ]
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE dpdt_impl(this, time, pde, fluidstate, pstate, dt, dpdtout)
    use equationbase_mod
    use fft
    IMPLICIT NONE
    class       (Tpart),             intent(inout) :: this
    class(EquationBase),             intent   (in) :: pde
    real      (kind=GP),             intent   (in) :: time, dt
    type   (GStateComp), target ,    intent   (in) :: fluidstate(:)
    type  (GPStateComp), target ,    intent   (in) :: pstate(:)
    type  (GPStateComp), target ,    intent(inout) :: dpdtout(:)
    complex   (KIND=GP), pointer, dimension(:,:,:) :: velc,velc2,vc,ac
    real      (KIND=GP), pointer, dimension(:,:,:) :: velr,tmp1,tmp2
    real      (kind=GP), pointer, dimension(:)     :: lbx,lby,lbz,lfx,lfy,lfz
    real      (kind=GP), pointer, dimension(:)     :: dpx,dpy,dpz,dvx,dvy,dvz
    real      (kind=GP), pointer, dimension(:)     :: pvx,pvy,pvz
    real      (kind=GP)                            :: rmp
    integer                                        :: m
    logical                                        :: bret

    CALL GTStart(this%htimers_(GPTIME_STEP))
    call this%workspace_%get_complex_tmp(velc,bret)
    if ( this%traits_%dokinelv ) call this%workspace_%get_complex_tmp(velc2,bret)
    call this%workspace_%get_real_tmp   (velr,bret)
    call this%workspace_%get_real_tmp   (tmp1,bret)
    call this%workspace_%get_real_tmp   (tmp2,bret)
    ! Particle-sized temporaries for the magnetic field and the current
    ! density at the particles (lvx_,lvy_,lvz_ hold the fluid velocity)
    call this%workspace_%get_pcomp_tmp  (lbx ,bret)
    call this%workspace_%get_pcomp_tmp  (lby ,bret)
    call this%workspace_%get_pcomp_tmp  (lbz ,bret)
    call this%workspace_%get_pcomp_tmp  (lfx ,bret)
    call this%workspace_%get_pcomp_tmp  (lfy ,bret)
    call this%workspace_%get_pcomp_tmp  (lfz ,bret)
    call this%AssignLagPos(pstate) ! We assign px_,py_,pz_ to the pstate

    select type (pde)
    class is (MagneticBase)
      if (this%nc_ .ne. pde%nc_) then
        stop "Testpart: # of components of the particles and pdes must be equal"
      endif

      ! IMPORTANT: pstate and dpdtout may alias the same array (some steppers
      ! can pass upout for both). We must be careful about ordering:
      !   1. Interpolate u_e, B and j with positions still intact
      !   2. Write position RHS (overwrites pstate(POSITION), but we're done)
      !   3. Compute velocity RHS (reads particle velocity before overwrite)

      ! Step 1: Interpolate the fluid (or electron) velocity to lvx_,
      ! lvy_, lvz_, the magnetic field to lbx, lby, lbz, and the current
      ! density to lfx, lfy, lfz. Only the first interpolation updates
      ! the interpolation points.
      rmp = 1.0_GP/(real(this%nd_(1),kind=GP)*real(this%nd_(2),kind=GP)* &
                    real(this%nd_(3),kind=GP))
      do m = 1,3
        vc => fluidstate(pde%VELOCITY+m-1)%ccomp
        if ( this%traits_%dokinelv ) then  ! u_e = u - dii j = u + dii Del^2 a
          ! COMPRESSIBLE: the electron velocity is u_e = u - dii j/rho.
          ! For compressible solvers divide -j (velc2) by the density
          ! (in the old code: divide(th,C14,C15,C16), a real-space
          ! product with the FFTs it requires) before adding it to u.
          ac => fluidstate(pde%MAGNETIC+m-1)%ccomp
          CALL laplak3(ac,velc2)
          CALL saxpby_c(velc,vc,rmp,velc2,this%traits_%dii*rmp)
        else
          CALL copy3(vc,velc)
          CALL scal3(velc,rmp)
        endif
        if (m.eq.1) then
          call tpart_c2lag(this,velc,this%lvx_,.true. ,velr,tmp1,tmp2)
        else if (m.eq.2) then
          call tpart_c2lag(this,velc,this%lvy_,.false.,velr,tmp1,tmp2)
        else
          call tpart_c2lag(this,velc,this%lvz_,.false.,velr,tmp1,tmp2)
        endif
      end do
      call tpart_magnetic(this,pde,fluidstate,lbx,lby,lbz,velc,velr,tmp1,tmp2)
      call tpart_current (this,pde,fluidstate,lfx,lfy,lfz,velc,velr,tmp1,tmp2)

      ! Steps 2 and 3 in one kernel over the particles: position RHS
      ! dx/dt = v_p and velocity RHS dv_p/dt = gyrof [(v_p-u_e) x B + eta j].
      ! pstate and dpdtout may alias (some steppers pass upout for
      ! both): each particle reads its velocity before writing its RHS.
      ! Pointers and a kernel (tpart_rhs) are used to help offloading.
      dpx => dpdtout(this%POSITION  )%rcomp
      dpy => dpdtout(this%POSITION+1)%rcomp
      dpz => dpdtout(this%POSITION+2)%rcomp
      dvx => dpdtout(this%VELOCITY  )%rcomp
      dvy => dpdtout(this%VELOCITY+1)%rcomp
      dvz => dpdtout(this%VELOCITY+2)%rcomp
      pvx => pstate (this%VELOCITY  )%rcomp
      pvy => pstate (this%VELOCITY+1)%rcomp
      pvz => pstate (this%VELOCITY+2)%rcomp
      call tpart_rhs(this%nparts_,this%lvx_,this%lvy_,this%lvz_,lbx,lby,lbz,  &
                     lfx,lfy,lfz,pvx,pvy,pvz,dpx,dpy,dpz,dvx,dvy,dvz,          &
                     this%invdel_,this%traits_%gyrof,this%traits_%gyeta)
    class default
      stop "Testpart: This solver does not support pdes without a magnetic field"
    end select

    call this%workspace_%free_pcomp_tmp  (lfz)
    call this%workspace_%free_pcomp_tmp  (lfy)
    call this%workspace_%free_pcomp_tmp  (lfx)
    call this%workspace_%free_pcomp_tmp  (lbz)
    call this%workspace_%free_pcomp_tmp  (lby)
    call this%workspace_%free_pcomp_tmp  (lbx)
    call this%workspace_%free_real_tmp   (tmp2)
    call this%workspace_%free_real_tmp   (tmp1)
    call this%workspace_%free_real_tmp   (velr)
    if ( this%traits_%dokinelv ) call this%workspace_%free_complex_tmp(velc2)
    call this%workspace_%free_complex_tmp(velc)
    CALL GTAcc(this%htimers_(GPTIME_STEP))
  END SUBROUTINE dpdt_impl


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Functions to compute fluid and particle couplings
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine null_feedback(this, pstate, feedback)
    class     (Tpart), intent   (in) :: this
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
    class     (TPart), intent(inout) :: this
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
          WRITE(*,*) 'Testpart EndStage (VDB): inconsistent particle count'
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
        WRITE(*,*) 'Testpart EndStage (VDB): inconsistent d.b.: expected: ', &
                   this%maxparts_, '; found: ', ng
        CALL this%ascii_write_lag(1, '.', trim(this%sstate_pos_) // 'err',   &
             '000', 0.0_GP, this%maxparts_,this%vdb_)
        STOP
      end if
    end if  ! GPEXCHTYPE_VDB
  end subroutine end_stage_impl


  ! ===================================================================
  ! Internal routines: fields at the particles and RHS kernel
  ! ===================================================================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Transforms the (already normalized) Fourier field velc to
  !! real space in velr and interpolates it to the particles
  !! in lag. Contents of velc and velr are lost.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine tpart_c2lag(this,velc,lag,doupdate,velr,tmp1,tmp2)
    use fft
    use grid
    use mpivars
    implicit none
    class     (Tpart), intent(inout)                          :: this
    complex(kind=GP), intent(inout), dimension(nz,ny,ista:iend) :: velc
    real   (kind=GP), intent(inout), dimension(nx,ny,ksta:kend) :: velr,tmp1,tmp2
    real   (kind=GP), intent(inout), dimension(*)               :: lag
    logical         , intent   (in)                             :: doupdate

    call fftp3d_complex_to_real(plancr,velc,velr,MPI_COMM_WORLD)
    call this%EulerToLag(lag,this%nparts_,velr,doupdate,tmp1,tmp2)
  end subroutine tpart_c2lag


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Computes the total magnetic field B = curl(a) + B_0 and
  !! interpolates it to the particles in lbx, lby, lbz. The
  !! interpolation points must be already updated.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine tpart_magnetic(this,pde,fluidstate,lbx,lby,lbz,velc,velr,tmp1,tmp2)
    use equationbase_mod
    use grid
    use mpivars
    implicit none
    class      (Tpart), intent(inout)                          :: this
    class(MagneticBase), intent   (in)                          :: pde
    type   (GStateComp), intent   (in), target                  :: fluidstate(:)
    complex(kind=GP), intent(inout), dimension(nz,ny,ista:iend) :: velc
    real   (kind=GP), intent(inout), dimension(nx,ny,ksta:kend) :: velr,tmp1,tmp2
    real   (kind=GP), intent(inout), dimension(*)               :: lbx,lby,lbz
    complex(kind=GP), pointer, dimension(:,:,:)                 :: ax,ay,az
    real   (kind=GP)                                            :: rmp,b0
    integer                                                     :: m

    rmp = 1.0_GP/(real(this%nd_(1),kind=GP)*real(this%nd_(2),kind=GP)* &
                  real(this%nd_(3),kind=GP))
    ax => fluidstate(pde%MAGNETIC  )%ccomp
    ay => fluidstate(pde%MAGNETIC+1)%ccomp
    az => fluidstate(pde%MAGNETIC+2)%ccomp
    do m = 1,3
      if (m.eq.1) then
        CALL rotor3(ay,az,velc,1)
      else if (m.eq.2) then
        CALL rotor3(ax,az,velc,2)
      else
        CALL rotor3(ax,ay,velc,3)
      endif
      if ( this%traits_%doB0 .and. (this%myrank_.eq.0) ) then ! b = b + B_0
        b0 = this%traits_%B0(m)/rmp
        CALL setmode3(velc,1,1,1,cmplx(b0,0.0_GP,kind=GP))
      endif
      CALL scal3(velc,rmp)
      if (m.eq.1) then
        call tpart_c2lag(this,velc,lbx,.false.,velr,tmp1,tmp2)
      else if (m.eq.2) then
        call tpart_c2lag(this,velc,lby,.false.,velr,tmp1,tmp2)
      else
        call tpart_c2lag(this,velc,lbz,.false.,velr,tmp1,tmp2)
      endif
    end do
  end subroutine tpart_magnetic


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Computes the current density j = curl(B) = -Del^2 a and
  !! interpolates it to the particles in lfx, lfy, lfz. The
  !! interpolation points must be already updated.
  !!
  !! COMPRESSIBLE: with the electron pressure correction (dokinelp)
  !! the dissipative part of E is eta j - (dii/2) grad(p)/rho. The
  !! simplest way to add it here is to interpolate, instead of j,
  !! the field f = j - (dii/(2 eta)) grad(p)/rho (or, if eta = 0,
  !! to pass gyrof and gyeta = gyrof*eta separately and add a third
  !! term in tpart_rhs). In the old code grad(p)/rho was computed
  !! by gradpstate(cp1,gam1,th,C11,C12,C13) from the density th and
  !! the equation of state, and combined in Fourier space as
  !! C14 = -gyrof*eta*C14 - 0.5*gyrof*dii*C11 (C14 = -j) before the
  !! inverse FFT. Both branches (dokinelp true/false) need a check
  !! that the pde is compressible; write_pstate_impl uses this
  !! routine to write jlg, so keep the plain current density
  !! available there (e.g., with an optional argument).
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine tpart_current(this,pde,fluidstate,lfx,lfy,lfz,velc,velr,tmp1,tmp2)
    use equationbase_mod
    use grid
    use mpivars
    implicit none
    class      (Tpart), intent(inout)                          :: this
    class(MagneticBase), intent   (in)                          :: pde
    type   (GStateComp), intent   (in), target                  :: fluidstate(:)
    complex(kind=GP), intent(inout), dimension(nz,ny,ista:iend) :: velc
    real   (kind=GP), intent(inout), dimension(nx,ny,ksta:kend) :: velr,tmp1,tmp2
    real   (kind=GP), intent(inout), dimension(*)               :: lfx,lfy,lfz
    complex(kind=GP), pointer, dimension(:,:,:)                 :: ac
    real   (kind=GP)                                            :: rmp
    integer                                                     :: m

    rmp = 1.0_GP/(real(this%nd_(1),kind=GP)*real(this%nd_(2),kind=GP)* &
                  real(this%nd_(3),kind=GP))
    do m = 1,3
      ac => fluidstate(pde%MAGNETIC+m-1)%ccomp
      CALL laplak3(ac,velc)    ! Del^2 a = -j
      CALL scal3(velc,-rmp)
      if (m.eq.1) then
        call tpart_c2lag(this,velc,lfx,.false.,velr,tmp1,tmp2)
      else if (m.eq.2) then
        call tpart_c2lag(this,velc,lfy,.false.,velr,tmp1,tmp2)
      else
        call tpart_c2lag(this,velc,lfz,.false.,velr,tmp1,tmp2)
      endif
    end do
  end subroutine tpart_current


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Internal kernel to compute right-hand side of n particles
  !!   dx/dt   = v_p/delta         (positions in grid units)
  !!   dv_p/dt = gyrof [ (v_p - u_e) x B + eta j ]
  !! with u_e in lv*, B in lb*, j in lf*, and gyeta = gyrof*eta.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine tpart_rhs(n,lvx,lvy,lvz,lbx,lby,lbz,lfx,lfy,lfz,pvx,pvy,pvz, &
                       dpx,dpy,dpz,dvx,dvy,dvz,invdel,gyrof,gyeta)
    implicit none
    integer      , intent(in)    :: n
    real(kind=GP), intent(in)    :: lvx(n),lvy(n),lvz(n),lbx(n),lby(n),lbz(n)
    real(kind=GP), intent(in)    :: lfx(n),lfy(n),lfz(n),pvx(n),pvy(n),pvz(n)
    real(kind=GP), intent(inout) :: dpx(n),dpy(n),dpz(n),dvx(n),dvy(n),dvz(n)
    real(kind=GP), intent(in)    :: invdel(3),gyrof,gyeta
    real(kind=GP)                :: vx,vy,vz,wx,wy,wz,bx,by,bz
    integer                      :: j
#if defined(GHOST_GPU)
!$omp target teams distribute parallel do if(target: gdev_active) &
!$omp   private(vx,vy,vz,wx,wy,wz,bx,by,bz)
#else
!$omp parallel do private(vx,vy,vz,wx,wy,wz,bx,by,bz)
#endif
    do j = 1,n
      vx = pvx(j); vy = pvy(j); vz = pvz(j)
      bx = lbx(j); by = lby(j); bz = lbz(j)
      wx = vx - lvx(j)                       ! v_p - u_e
      wy = vy - lvy(j)
      wz = vz - lvz(j)
      dpx(j) = vx*invdel(1)
      dpy(j) = vy*invdel(2)
      dpz(j) = vz*invdel(3)
      dvx(j) = gyeta*lfx(j) + gyrof*(wy*bz - wz*by)
      dvy(j) = gyeta*lfy(j) + gyrof*(wz*bx - wx*bz)
      dvz(j) = gyeta*lfz(j) + gyrof*(wx*by - wy*bx)
    end do
  end subroutine tpart_rhs


  ! ===================================================================
  ! Output methods
  ! ===================================================================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Function to compute and write particle states.
  !! Writes: positions (xlg), fluid velocity (vlg), test
  !!         particle velocity (vtp), magnetic field (blg) and
  !!         current density (jlg) at the particles
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine write_pstate_impl(this, time, pde, fluidstate, pstate)
    use equationbase_mod
    use particlebase_mod
    use pseudospec_fluid
    use status
    use pstatus
    use fft
    class       (TPart),             intent(inout) :: this
    class(EquationBase),             intent   (in) :: pde
    type   (GStateComp), target ,    intent   (in) :: fluidstate(:)
    type  (GPStateComp), target ,    intent   (in) :: pstate(:)
    real      (kind=GP),             intent   (in) :: time
    complex   (kind=GP), pointer, dimension(:,:,:) :: velc, velc2, vc
    real      (kind=GP), pointer, dimension(:,:,:) :: velr,tmp1,tmp2
    real      (kind=GP)                            :: rmp
    integer                                        :: m
    logical                                        :: bret,wasdev

    call this%workspace_%get_complex_tmp(velc,bret)
    call this%workspace_%get_real_tmp   (velr,bret)
    call this%workspace_%get_real_tmp   (tmp1,bret)
    call this%workspace_%get_real_tmp   (tmp2,bret)
    call this%AssignLagPos(pstate)
    ! The interpolations run on the device copies (the fluid state is
    ! current there); the I/O routines switch to the host copies
    wasdev = gdev_active
    gdev_active = .TRUE.

    select type (pde)
    class is (MagneticBase)
      rmp = 1.0_GP/(real(this%nd_(1),kind=GP)*real(this%nd_(2),kind=GP)* &
                    real(this%nd_(3),kind=GP))
      ! Interpolate fluid velocity to particle positions
      do m = 1,3
        vc => fluidstate(pde%VELOCITY+m-1)%ccomp
        CALL copy3(vc,velc)
        CALL scal3(velc,rmp)
        if (m.eq.1) then
          call tpart_c2lag(this,velc,this%lvx_,.true. ,velr,tmp1,tmp2)
        else if (m.eq.2) then
          call tpart_c2lag(this,velc,this%lvy_,.false.,velr,tmp1,tmp2)
        else
          call tpart_c2lag(this,velc,this%lvz_,.false.,velr,tmp1,tmp2)
        endif
      end do
      ! Write positions and Lagrangian fluid velocity
      WRITE(lgext,lgfmtext) pind
      CALL this%io_write_pdb(1,this%odir_,trim(this%sstate_pos_),lgext,time)
      CALL this%io_write_vec(1,this%odir_,trim(this%sstate_lag_),lgext,time)
      ! Write test particle velocity
      call gcopy(this%lvx_, pstate(this%VELOCITY  )%rcomp)
      call gcopy(this%lvy_, pstate(this%VELOCITY+1)%rcomp)
      call gcopy(this%lvz_, pstate(this%VELOCITY+2)%rcomp)
      CALL this%io_write_vec(1,this%odir_,trim(this%sstate_vel_),lgext,time)
      ! Write magnetic field (including the guide field) at the particles
      call tpart_magnetic(this,pde,fluidstate,this%lvx_,this%lvy_,this%lvz_, &
                          velc,velr,tmp1,tmp2)
      CALL this%io_write_vec(1,this%odir_,'blg',lgext,time)
      ! Write current density at the particles
      call tpart_current (this,pde,fluidstate,this%lvx_,this%lvy_,this%lvz_, &
                          velc,velr,tmp1,tmp2)
      CALL this%io_write_vec(1,this%odir_,'jlg',lgext,time)
! partlod >= 2: write Lagrangian vorticity and strain-rate tensor
      if ( this%traits_%partlod .ge. 2 ) then
! Write Lagrangian vorticity components
        CALL rotor3(fluidstate(pde%VELOCITY+1)%ccomp, &
                    fluidstate(pde%VELOCITY+2)%ccomp, velc, 1)
        CALL scal3(velc,rmp)
        call tpart_c2lag(this,velc,this%lvx_,.false.,velr,tmp1,tmp2)
        CALL rotor3(fluidstate(pde%VELOCITY  )%ccomp, &
                    fluidstate(pde%VELOCITY+2)%ccomp, velc, 2)
        CALL scal3(velc,rmp)
        call tpart_c2lag(this,velc,this%lvy_,.false.,velr,tmp1,tmp2)
        CALL rotor3(fluidstate(pde%VELOCITY  )%ccomp, &
                    fluidstate(pde%VELOCITY+1)%ccomp, velc, 3)
        CALL scal3(velc,rmp)
        call tpart_c2lag(this,velc,this%lvz_,.false.,velr,tmp1,tmp2)
        CALL this%io_write_vec(1,this%odir_,'wlg',lgext,time)
! Write strain-rate tensor components
        call this%workspace_%get_complex_tmp(velc2,bret)
        ! S11 = dv_x/dx
        CALL derivk3(fluidstate(pde%VELOCITY  )%ccomp, velc, 1)
        CALL scal3(velc,rmp)
        CALL fftp3d_complex_to_real(plancr,velc,velr,MPI_COMM_WORLD)
        CALL this%io_write_euler(1,this%odir_,'s11',lgext,time,velr,.false.,tmp1,tmp2)
        ! S12 = 0.5*(dv_x/dy + dv_y/dx)
        CALL derivk3(fluidstate(pde%VELOCITY  )%ccomp, velc,  2)
        CALL derivk3(fluidstate(pde%VELOCITY+1)%ccomp, velc2, 1)
        CALL saxpby_c(velc,velc,0.5_GP*rmp,velc2,0.5_GP*rmp)
        CALL fftp3d_complex_to_real(plancr,velc,velr,MPI_COMM_WORLD)
        CALL this%io_write_euler(1,this%odir_,'s12',lgext,time,velr,.false.,tmp1,tmp2)
        ! S13 = 0.5*(dv_x/dz + dv_z/dx)
        CALL derivk3(fluidstate(pde%VELOCITY  )%ccomp, velc,  3)
        CALL derivk3(fluidstate(pde%VELOCITY+2)%ccomp, velc2, 1)
        CALL saxpby_c(velc,velc,0.5_GP*rmp,velc2,0.5_GP*rmp)
        CALL fftp3d_complex_to_real(plancr,velc,velr,MPI_COMM_WORLD)
        CALL this%io_write_euler(1,this%odir_,'s13',lgext,time,velr,.false.,tmp1,tmp2)
        ! S22 = dv_y/dy
        CALL derivk3(fluidstate(pde%VELOCITY+1)%ccomp, velc, 2)
        CALL scal3(velc,rmp)
        CALL fftp3d_complex_to_real(plancr,velc,velr,MPI_COMM_WORLD)
        CALL this%io_write_euler(1,this%odir_,'s22',lgext,time,velr,.false.,tmp1,tmp2)
        ! S23 = 0.5*(dv_y/dz + dv_z/dy)
        CALL derivk3(fluidstate(pde%VELOCITY+1)%ccomp, velc,  3)
        CALL derivk3(fluidstate(pde%VELOCITY+2)%ccomp, velc2, 2)
        CALL saxpby_c(velc,velc,0.5_GP*rmp,velc2,0.5_GP*rmp)
        CALL fftp3d_complex_to_real(plancr,velc,velr,MPI_COMM_WORLD)
        CALL this%io_write_euler(1,this%odir_,'s23',lgext,time,velr,.false.,tmp1,tmp2)
        call this%workspace_%free_complex_tmp(velc2)
      endif
    class default
      stop "Testpart: This solver does not support pdes without a magnetic field"
    end select

    call this%workspace_%free_complex_tmp(velc)
    call this%workspace_%free_real_tmp   (velr)
    call this%workspace_%free_real_tmp   (tmp1)
    call this%workspace_%free_real_tmp   (tmp2)
    gdev_active = wasdev
  end subroutine write_pstate_impl


  ! ===================================================================
  ! Solver specific methods
  ! ===================================================================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Constructor
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE Tpart_ctor(this,infile, pde, workspace, pstate, pstate_cpy)
    USE equationbase_mod
    USE mhd_mod,   ONLY: MHDSolver
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
    CLASS       (TPart), intent(inout)              :: this
    class(EquationBase), intent   (in)              :: pde
    type   (GWorkspace), intent(inout), target      :: workspace
    type  (GPStateComp), intent(inout), allocatable :: pstate(:), pstate_cpy(:)
    character   (len=*), intent   (in)              :: infile
    integer                                         :: tsta,tend,num_components
    integer                                         :: j,szreal
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
        WRITE(*,*) 'Testpart_ctor: Not enough timers available'
        STOP
      ENDIF
    ENDDO

    ! Initialize communicators
    CALL this%gpcomm_%GPartComm_ctor(GPCOMM_INTRFC_SF,this%partbuff_, &
         this%nd_,this%intorder_-1,this%comm_,this%htimers_(GPTIME_COMM))
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

    ! Call init (reads &testpart namelist, sets POSITION/VELOCITY indices)
    call this%init()

    ! The magnetic diffusivity and the guide field are those of the pde.
    ! The traits are private to each solver class, so each magnetic
    ! solver needs its own "type is" clause here (and a USE of its
    ! module above). To add a solver: copy eta, doB0 and B0 from its
    ! traits, and any other trait the force needs.
    ! COMPRESSIBLE: for compressible MHD solvers also copy here the
    ! parameters of the equation of state needed for grad(p)/rho (see
    ! the COMPRESSIBLE notes in TestTraits and tpart_current), and set a
    ! trait flagging that the density is available (the electron
    ! velocity and the electron pressure correction need it). If the
    ! solver stores b instead of a, tpart_magnetic and tpart_current
    ! must also branch on that trait (rotor3/laplak3 assume a).
    select type (pde)
    type is (MHDSolver)
      this%traits_% eta = pde%traits_%eta
      this%traits_%doB0 = pde%traits_%doB0
      this%traits_%  B0 = pde%traits_%B0
    class default
      stop "Testpart_ctor: eta and B_0 are only known for the MHD solver"
    end select
    this%traits_%gyeta = this%traits_%gyrof*this%traits_%eta

    ! Instantiate interp operation
    CALL this%intop_%GPSplineInt_ctor(3,this%nd_,this%libnds_,this%lxbnds_, &
         this%tibnds_,this%intorder_,this%partbuff_,this%gpcomm_,&
         this%htimers_(GPTIME_DATAEX),this%htimers_(GPTIME_TRANSP))

    ! Allocate particle arrays
    CALL MPI_TYPE_SIZE(GC_REAL,szreal,this%ierr_)
    CALL galloc(this%id_    ,this%partbuff_)
    CALL galloc(this%tmpint_,this%partbuff_)
    CALL galloc(this%iwrk_  ,this%partbuff_)
    CALL galloc(this%ptmp0_ ,3,this%partbuff_)
    CALL galloc(this%gptmp0_,3,this%partbuff_)
    num_components = this%state_size()  ! Returns 6
    CALL GPState_alloc(pstate    , num_components, this%partbuff_)
    CALL GPState_alloc(pstate_cpy, num_components, this%partbuff_)
    call this%workspace_%set_nparts(this%partbuff_)
    call this%workspace_%init_pcomp_arrays(this%partbuff_)
    call this%workspace_%get_pcomp_tmp(this%lvx_,bret)
    call this%workspace_%get_pcomp_tmp(this%lvy_,bret)
    call this%workspace_%get_pcomp_tmp(this%lvz_,bret)
    IF ( this%iexchtype_.EQ.GPEXCHTYPE_VDB ) THEN
      CALL galloc(this%vdb_,3,this%partbuff_)
    ENDIF
  END SUBROUTINE Tpart_ctor


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Destructor
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE Tpart_dtor(this)
    IMPLICIT NONE
    TYPE(Tpart),INTENT(INOUT) :: this
    integer                   :: j
    CALL gfree(this%id_)
    CALL gfree(this%tmpint_)
    CALL gfree(this%iwrk_)
    CALL gfree(this%vdb_)
    CALL gfree(this%ptmp0_)
    CALL gfree(this%gptmp0_)
    IF ( ASSOCIATED   (this%px_) ) NULLIFY       (this%px_)
    IF ( ASSOCIATED   (this%py_) ) NULLIFY       (this%py_)
    IF ( ASSOCIATED   (this%pz_) ) NULLIFY       (this%pz_)
    IF ( ASSOCIATED(this%lvx_) ) CALL this%workspace_%free_pcomp_tmp(this%lvx_)
    IF ( ASSOCIATED(this%lvy_) ) CALL this%workspace_%free_pcomp_tmp(this%lvy_)
    IF ( ASSOCIATED(this%lvz_) ) CALL this%workspace_%free_pcomp_tmp(this%lvz_)
    DO j = 1, GPMAXTIMERS
      CALL GTFree(this%htimers_(j))
    ENDDO
  END SUBROUTINE Tpart_dtor


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Function to compute number of state members (equations)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PURE function state_size_impl(this) result(num)
    class(Tpart), intent(in) :: this
    integer                  :: num
    num = 2 * this%nc_          ! 3 positions + 3 velocities
  end function state_size_impl

end module testpart_mod
