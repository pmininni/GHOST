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
!                     electron velocity u_e = u - dii j (u - dii j/rho in
!                     compressible solvers) if dokinelv is .TRUE.
!                     (Hall-MHD). In compressible solvers the electron
!                     pressure can be added to E (dokinelp):
!                       E = -u_e x B + eta j - (dii/2) Grad(p)/rho
!                     with Grad(p)/rho = Grad(h), h the enthalpy of the
!                     polytropic gas of the solver. The magnetic
!                     diffusivity, the guide field and the equation of
!                     state are those of the fluid solver (MHD or CMHD).
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
!                        u_e = u - dii j/rho (Hall-MHD correction; rho=1
!                        in incompressible solvers)
!              dokinelp: .false.=Ohmic electric field eta j [default],
!                        .true.=adds the electron pressure to E,
!                        - (dii/2) Grad(p)/rho (compressible solvers only)
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
    logical       :: dokinelp = .false. ! .true.=electron pressure in E
    ! Compressible solvers: density available, equation of state
    logical       :: compressible = .false.
    real(kind=GP) :: cp1      = 0.0_GP  ! enthalpy h = cp1 rho^gam1/2
    real(kind=GP) :: gam1     = 0.0_GP  ! gamma - 1 (from the pde)
  end type

  ! ================= Solver ==========================================
  type, extends(VelocParticleBase) :: Tpart
    ! Member data:
    type (TestTraits) :: traits_
  CONTAINS
    procedure, public :: init          =>          init_impl
    procedure, public :: dpdt          =>          dpdt_impl
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
    logical                            :: dokinelv, dokinelp
    character(len=128)                 :: pidir, podir
    namelist/ testpart / pidir,podir,partlod,gyrof,dii,dokinelv,dokinelp

    this%POSITION = 1
    this%VELOCITY = this%POSITION + this%nc_
    ! Defaults
    pidir    = this%idir_ ! Set to the pde class idir at ctor
    podir    = this%odir_ ! Set to the pde class odir at ctor
    partlod  = 1
    gyrof    = 1.0_GP
    dii      = 0.0_GP
    dokinelv = .false.
    dokinelp = .false.
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
    call MPI_BCAST(dokinelp,1  ,MPI_LOGICAL  ,0,MPI_COMM_WORLD,ierr)

    this%idir_ = pidir ! If present in &testpart, replaces the class default idir
    this%odir_ = podir ! If present in &testpart, replaces the class default odir
    this%sstate_pos_ = 'xlg' ! state name of positions
    this%sstate_lag_ = 'vlg' ! state name of Lagrangian velocities
    this%sstate_vel_ = 'vtp' ! state name of particles velocities
    this%traits_% partlod = partlod
    this%traits_%   gyrof = gyrof
    this%traits_%     dii = dii
    this%traits_%dokinelv = dokinelv
    this%traits_%dokinelp = dokinelp
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
    use pseudospec_compressible, only: divide
    use fft
    IMPLICIT NONE
    class       (Tpart),             intent(inout) :: this
    class(EquationBase),             intent   (in) :: pde
    real      (kind=GP),             intent   (in) :: time, dt
    type   (GStateComp), target ,    intent   (in) :: fluidstate(:)
    type  (GPStateComp), target ,    intent   (in) :: pstate(:)
    type  (GPStateComp), target ,    intent(inout) :: dpdtout(:)
    complex   (KIND=GP), pointer, dimension(:,:,:) :: velc,velc2,vc,ac,rho
    complex   (KIND=GP), pointer, dimension(:,:,:) :: C1,C2,C3
    real      (KIND=GP), pointer, dimension(:,:,:) :: velr,tmp1,tmp2
    real      (kind=GP), pointer, dimension(:)     :: lbx,lby,lbz,lfx,lfy,lfz
    real      (kind=GP), pointer, dimension(:)     :: dpx,dpy,dpz,dvx,dvy,dvz
    real      (kind=GP), pointer, dimension(:)     :: pvx,pvy,pvz
    real      (kind=GP)                            :: rmp,cf
    integer                                        :: m
    logical                                        :: bret,dorho,dolap

    ! The density enters the electron velocity (dokinelv) and the
    ! electron pressure (dokinelp) in the compressible solvers; the
    ! incompressible electron velocity only needs Del^2 a
    dorho = this%traits_%compressible .and. &
            (this%traits_%dokinelv .or. this%traits_%dokinelp)
    dolap = this%traits_%dokinelv .and. .not.this%traits_%compressible

    CALL GTStart(this%htimers_(GPTIME_STEP))
    call this%workspace_%get_complex_tmp(velc,bret)
    if ( dolap ) call this%workspace_%get_complex_tmp(velc2,bret)
    if ( dorho ) then
      call this%workspace_%get_complex_tmp(C1,bret)
      call this%workspace_%get_complex_tmp(C2,bret)
      call this%workspace_%get_complex_tmp(C3,bret)
    endif
    call this%workspace_%get_real_tmp   (velr,bret)
    call this%workspace_%get_real_tmp   (tmp1,bret)
    call this%workspace_%get_real_tmp   (tmp2,bret)
    ! Particle-sized temporaries for the magnetic field and the current
    ! density (or the dissipative electric field) at the particles
    ! (lvx_,lvy_,lvz_ hold the fluid velocity)
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
      rho => null()
      select type (pde)
      class is (CompMagneticBase)
        rho => fluidstate(pde%DENSITY)%ccomp
      end select

      ! IMPORTANT: pstate and dpdtout may alias the same array (some steppers
      ! can pass upout for both). We must be careful about ordering:
      !   1. Interpolate u_e, B and j with positions still intact
      !   2. Write position RHS (overwrites pstate(POSITION), but we're done)
      !   3. Compute velocity RHS (reads particle velocity before overwrite)

      ! Step 1: Interpolate the fluid (or electron) velocity to lvx_,
      ! lvy_, lvz_, the magnetic field to lbx, lby, lbz, and the current
      ! density (or the dissipative electric field) to lfx, lfy, lfz.
      ! Only the first interpolation updates the interpolation points.
      rmp = 1.0_GP/(real(this%nd_(1),kind=GP)*real(this%nd_(2),kind=GP)* &
                    real(this%nd_(3),kind=GP))
      if ( this%traits_%dokinelv .and. this%traits_%compressible ) then
        CALL laplak3(fluidstate(pde%MAGNETIC  )%ccomp,C1) ! Del^2 a = -j
        CALL laplak3(fluidstate(pde%MAGNETIC+1)%ccomp,C2)
        CALL laplak3(fluidstate(pde%MAGNETIC+2)%ccomp,C3)
        CALL divide(rho,C1,C2,C3)                         ! -j/rho
      endif
      do m = 1,3
        vc => fluidstate(pde%VELOCITY+m-1)%ccomp
        if ( this%traits_%dokinelv ) then
          if ( this%traits_%compressible ) then ! u_e = u - dii j/rho
            if (m.eq.1) then
              ac => C1
            else if (m.eq.2) then
              ac => C2
            else
              ac => C3
            endif
            CALL saxpby_c(velc,vc,rmp,ac,this%traits_%dii*rmp)
          else                                  ! u_e = u - dii j = u + dii Del^2 a
            ac => fluidstate(pde%MAGNETIC+m-1)%ccomp
            CALL laplak3(ac,velc2)
            CALL saxpby_c(velc,vc,rmp,velc2,this%traits_%dii*rmp)
          endif
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
      if ( this%traits_%dokinelp ) then
        ! Dissipative electric field eta j - (dii/2) Grad(h) in lf*,
        ! multiplied by gyrof in the kernel
        call tpart_edissip(this,pde,fluidstate,rho,lfx,lfy,lfz,velc,C1,C2,C3, &
                           velr,tmp1,tmp2)
        cf = this%traits_%gyrof
      else
        ! Current density in lf*, multiplied by gyrof*eta in the kernel
        call tpart_current(this,pde,fluidstate,lfx,lfy,lfz,velc,velr,tmp1,tmp2)
        cf = this%traits_%gyeta
      endif

      ! Steps 2 and 3 in one kernel over the particles: position RHS
      ! dx/dt = v_p and velocity RHS dv_p/dt = gyrof [(v_p-u_e) x B + E_d].
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
                     this%invdel_,this%traits_%gyrof,cf)
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
    if ( dorho ) then
      call this%workspace_%free_complex_tmp(C3)
      call this%workspace_%free_complex_tmp(C2)
      call this%workspace_%free_complex_tmp(C1)
    endif
    if ( dolap ) call this%workspace_%free_complex_tmp(velc2)
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


  ! ===================================================================
  ! Internal routines: RHS kernel and fields at the particles
  ! ===================================================================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Internal kernel to compute right-hand side of n particles
  !!   dx/dt   = v_p/delta         (positions in grid units)
  !!   dv_p/dt = gyrof [ (v_p - u_e) x B ] + gyeta lf
  !! with u_e in lv*, B in lb*, and in lf* either the current density
  !! j (then gyeta = gyrof*eta) or the dissipative electric field
  !! E_d = eta j - (dii/2) Grad(h) (then gyeta = gyrof).
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
    class    (Tpart), intent(inout)                             :: this
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
    class       (Tpart), intent(inout)                          :: this
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
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine tpart_current(this,pde,fluidstate,lfx,lfy,lfz,velc,velr,tmp1,tmp2)
    use equationbase_mod
    use grid
    use mpivars
    implicit none
    class       (Tpart), intent(inout)                          :: this
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
  !! Computes the dissipative part of the electric field of the
  !! compressible solvers with the electron pressure,
  !!   E_d = eta j - (dii/2) Grad(p)/rho = -eta Del^2 a - (dii/2) Grad(h)
  !! with h = cp1 rho^gam1/2 the enthalpy of the polytropic gas,
  !! and interpolates it to the particles in lfx, lfy, lfz. The
  !! interpolation points must be already updated. gx, gy, gz are
  !! field-sized temporaries for Grad(h).
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine tpart_edissip(this,pde,fluidstate,rho,lfx,lfy,lfz,velc,gx,gy,gz, &
                           velr,tmp1,tmp2)
    use equationbase_mod
    use pseudospec_compressible, only: gradpstate
    use grid
    use mpivars
    implicit none
    class       (Tpart), intent(inout)                          :: this
    class(MagneticBase), intent   (in)                          :: pde
    type   (GStateComp), intent   (in), target                  :: fluidstate(:)
    complex(kind=GP), intent   (in), dimension(nz,ny,ista:iend) :: rho
    complex(kind=GP), intent(inout), dimension(nz,ny,ista:iend) :: velc,gx,gy,gz
    real   (kind=GP), intent(inout), dimension(nx,ny,ksta:kend) :: velr,tmp1,tmp2
    real   (kind=GP), intent(inout), dimension(*)               :: lfx,lfy,lfz
    complex(kind=GP), pointer, dimension(:,:,:)                 :: ac
    real   (kind=GP)                                            :: rmp,ceta,cdii
    integer                                                     :: m

    rmp  = 1.0_GP/(real(this%nd_(1),kind=GP)*real(this%nd_(2),kind=GP)* &
                   real(this%nd_(3),kind=GP))
    ceta = -this%traits_%eta*rmp          ! eta j = -eta Del^2 a
    cdii = -0.5_GP*this%traits_%dii*rmp   ! -(dii/2) Grad(h)
    CALL gradpstate(this%traits_%cp1,this%traits_%gam1,rho,gx,gy,gz) ! Grad(h)
    do m = 1,3
      ac => fluidstate(pde%MAGNETIC+m-1)%ccomp
      CALL laplak3(ac,velc)               ! Del^2 a = -j
      if (m.eq.1) then
        CALL saxpby_c(velc,velc,ceta,gx,cdii)
        call tpart_c2lag(this,velc,lfx,.false.,velr,tmp1,tmp2)
      else if (m.eq.2) then
        CALL saxpby_c(velc,velc,ceta,gy,cdii)
        call tpart_c2lag(this,velc,lfy,.false.,velr,tmp1,tmp2)
      else
        CALL saxpby_c(velc,velc,ceta,gz,cdii)
        call tpart_c2lag(this,velc,lfz,.false.,velr,tmp1,tmp2)
      endif
    end do
  end subroutine tpart_edissip


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
    USE cmhd_mod,  ONLY: CMHDSolver
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
    this%sclass_      =  'testpart'
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
         this%nd_,GPSI_NZGHOST,this%comm_,this%htimers_(GPTIME_COMM))
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

    ! The magnetic diffusivity, the guide field and (for compressible
    ! solvers) the equation of state are those of the pde. The traits
    ! are private to each solver class, so each magnetic solver needs
    ! its own "type is" clause here (and a USE of its module above);
    ! tpart_magnetic and tpart_current assume the solver stores the
    ! vector potential a.
    select type (pde)
    type is (MHDSolver)
      this%traits_%eta  = pde%traits_%eta
      this%traits_%doB0 = pde%traits_%doB0
      this%traits_%B0   = pde%traits_%B0
      this%traits_%compressible = .false.
    type is (CMHDSolver)
      this%traits_%eta  = pde%traits_%eta
      this%traits_%doB0 = pde%traits_%doB0
      this%traits_%B0   = pde%traits_%B0
      this%traits_%cp1  = pde%traits_%cp1
      this%traits_%gam1 = pde%traits_%gam1
      this%traits_%compressible = .true.
    class default
      stop "Testpart_ctor: the traits of this magnetic solver are not known"
    end select
    this%traits_%gyeta = this%traits_%gyrof*this%traits_%eta
    if ( this%traits_%dokinelp .and. .not.this%traits_%compressible ) then
      if ( this%myrank_ .eq. 0 ) then
        WRITE(*,*) 'Testpart_ctor: dokinelp (electron pressure) requires', &
                   ' a compressible solver (CMHD)'
      endif
      stop
    endif

    ! Instantiate interp operation
    CALL this%intop_%GPSplineInt_ctor(3,this%nd_,this%libnds_,this%lxbnds_, &
         this%tibnds_,this%partbuff_,this%gpcomm_,&
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
