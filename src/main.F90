!=================================================================
      PROGRAM MAIN
!=================================================================
! GHOST code: Geophysical High Order Suite for Turbulence
!
! Numerically integrates several fluid dynamics equations
! in 2 and 3 dimensions with periodic boundary conditions
! and external forcing. A pseudo-spectral method is used to
! compute spatial derivatives, while different time steppings
! can be used to evolve the system in the time domain.
!
! Notation: index 'i' is 'x'
!           index 'j' is 'y'
!           index 'k' is 'z'
!
! 2003 Pablo D. Mininni.
!      Department of Physics,
!      Facultad de Ciencias Exactas y Naturales.
!      Universidad de Buenos Aires.
!      e-mail: mininni@df.uba.ar
!
! References:
! Mininni PD, Rosenberg DL, Reddy R, Pouquet A.; P.Comp.37, 123 (2011)
! Rosenberg DL, Mininni PD, Reddy R, Pouquet A.: Atmosph.11, 178 (2020)
!=================================================================

! Modules
      USE commtypes
      USE filefmt
      USE iovar
      USE fft
      USE threads
      USE gdevice
      USE boxsize
      USE status
      USE pstatus
      USE gtimer
      USE gbench
      USE ic_factory
      USE force_factory
      USE equation_factory
      USE particle_factory
      USE icp_factory
      USE stepper_factory
      IMPLICIT NONE

! Arrays for the field and particle states, workspace, I/O, and solver classes
      TYPE   (GStateComp), ALLOCATABLE :: field(:),field_nxt(:),force (:)
      TYPE  (GPStateComp), ALLOCATABLE :: part (:),part_nxt (:)
      TYPE   (GWorkspace), TARGET      :: workspace
      TYPE       (ioplan)              :: planio
      CLASS(EquationBase), ALLOCATABLE :: fluid
      CLASS(ParticleBase), ALLOCATABLE :: particle
      CLASS(GStepperBase), ALLOCATABLE :: stepper
      CLASS     (icChain), ALLOCATABLE :: iclist(:)
      CLASS  (forceChain), ALLOCATABLE :: forcemethod(:)
      CLASS    (icpChain), ALLOCATABLE :: icplist(:)

! Auxiliary variables
      REAL(KIND=GP) :: time
      INTEGER       :: t, num_components
      LOGICAL       :: fstatic

! Initialization
! Initializes the MPI and I/O libraries
      CALL MPI_INIT_THREAD(MPI_THREAD_FUNNELED,provided,ierr)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)

! Initializes the grid. This must be done early to have nx, ny, nz.
      CALL grid_init('parameter.inp')

! Binds this MPI task to a GPU in offload builds (no-op otherwise).
! GHOST assumes the number of MPI tasks per node is a multiple of the
! number of GPUs per node; the user must ensure this condition.
      CALL device_init(myrank)

! Initialization of fluid and particles integration parameters
      CALL status_init ('parameter.inp')
      CALL pstatus_init('parameter.inp')

! I/O initialization
      CALL range(1,nx/2+1,nprocs,myrank,ista,iend)
      CALL range(1,nz,nprocs,myrank,ksta,kend)
      CALL io_init(myrank,(/nx,ny,nz/),ksta,kend,planio)

! Now we can initialize the PDE methods and read the particles status
      fluid        = init_pdes_from_file('parameter.inp')
      if (dopart) particle = init_particles_from_file('parameter.inp')
      CALL workspace%initialize_pool(NUMTMPREAL,NUMTMPCOMP,NUMTMPPART)
      CALL fluid%Solver_ctor('parameter.inp',workspace,planio)
      num_components = fluid%state_size()
      CALL GState_alloc(field    , num_components)
      CALL GState_alloc(field_nxt, num_components)
      CALL GState_alloc(force    , num_components)
      iclist       = init_ic_from_file(     'parameter.inp')
      forcemethod  = init_forcing_from_file('parameter.inp',workspace)

! Initialization of the numerical domain
      CALL box_init('parameter.inp')

! Initialization of the particles
      if (dopart) then
         CALL particle%part_ctor('parameter.inp',fluid,workspace,part,part_nxt)
         icplist = init_icp_from_file('parameter.inp')
      endif

! Initializes the FFT library. This must be done at this
! stage as it requires status and benchmark initialization.
! Use FFTW_ESTIMATE or FFTW_MEASURE in short runs
! Use FFTW_PATIENT or FFTW_EXHAUSTIVE in very long runs
! FFTW_MEASURE is used below: at 128^3 on two MPI ranks it makes the FFTs
! about 23% faster than FFTW_ESTIMATE, 13 to 18% of the whole time step.
      nth = 1
!$    nth = omp_get_max_threads()
!$    CALL fftp3d_init_threads(ierr)
      CALL GTBenchInit(bench,ihcpu1,ihomp1,ihwtm1,ihcpu2,ihomp2,ihwtm2)
      IF (bench.eq.2) THEN
         CALL GTStart(ihcpu2); CALL GTStart(ihomp2); CALL GTStart(ihwtm2)
      ENDIF
      CALL fftp3d_create_plan(planrc,(/nx,ny,nz/),FFTW_REAL_TO_COMPLEX, &
                             FFTW_MEASURE)
      CALL fftp3d_create_plan(plancr,(/nx,ny,nz/),FFTW_COMPLEX_TO_REAL, &
                             FFTW_MEASURE)
      IF (bench.eq.2) THEN
         CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
         CALL GTStop(ihcpu2);  CALL GTStop(ihomp2);  CALL GTStop(ihwtm2)
      ENDIF

! Initial states
 IC : IF (stat.eq.0) THEN                 ! If stat=0 we start a new run
        ini  = 1
        sind = 0                          ! index for the spectrum
        tind = 0                          ! index for the binaries
        pind = 0                          ! index for the particles
        timet = tstep
        timep = pstep
        timec = cstep
        times = sstep
        timef = 0
      ELSE                    ! If stat.ne.0 a previous run is continued
        ini = int((stat-1)*tstep) + 1
        tind = int(stat)
        sind = int(real(ini,kind=GP)/real(sstep,kind=GP)+1)
        pind = int((stat-1)*lgmult+1)
        timet = 0
        timep = 0
        timec = int(modulo(float(ini-1),float(cstep)))
        times = int(modulo(float(ini-1),float(sstep)))
        timef = int(modulo(float(ini-1),float(fstep)))
      ENDIF IC
      CALL init_allstates(iclist,fluid,field)
      CALL GState_update_to(field)  ! Device copies (no-op in host builds)
      CALL GState_copy(field_nxt,field) ! nxt is used by I/O and all steppers
      CALL init_forcing(forcemethod,fluid,force)
      CALL GState_update_to(force)
      fstatic = forcing_is_static(forcemethod)
      if (dopart) then
#if defined(GHOST_GPU)
         PRINT *,'Particles are not supported yet in offload builds'
         STOP
#endif
         CALL init_allpstates(icplist,fluid,field,particle,part)
         if (size(part(1)%rcomp) .ne. size(part_nxt(1)%rcomp)) then
            call GPState_resize(part_nxt,particle%partbuff_) ! We resize part_nxt
         endif
         part_nxt = part ! We also update part_nxt
      endif

! Sets up the time stepper
      if (dopart) then
        stepper = build_stepper_from_file('parameter.inp',workspace,fluid,particle)
      else
        stepper = build_stepper_from_file('parameter.inp',workspace,fluid)
      endif

! Time integration scheme starts here. In offload builds the fields
! are worked on the device inside the time step (gdev_active set) and
! on their host copies elsewhere: the host copies of the fields are
! refreshed before any output, and the forcing, which is computed on
! the host, is copied to the device after each update.
! If we are doing a benchmark, we measure cputime before
! starting. We also re-inititialize the fftp timers.
      IF (bench.eq.1) THEN
         ffttime  = 0.D00; tratime  = 0.0D0; comtime  = 0.D00; tottime  = 0.0D0
         CALL GTStart(ihcpu1); CALL GTStart(ihomp1); CALL GTStart(ihwtm1)
      ENDIF

 RK : DO t = ini,step
         time = (t-1)*dt
! Refreshes the host copies of the fields if any output is due
         IF (((timet.eq.tstep).or.(timec.eq.cstep).or.(times.eq.sstep) &
              .or.(dopart.and.(timep.eq.pstep))).and.(bench.eq.0)) THEN
            CALL GState_update_from(field_nxt)
         ENDIF
! Every 'tstep' steps, stores the fields in binary files
         IF ((timet.eq.tstep).and.(bench.eq.0)) THEN
            timet = 0
            tind = tind+1
            CALL fluid%write_states(field_nxt, planio)
         ENDIF

! Every 'pstep' steps, stores the particle states
         if (dopart) then
            IF ((timep.eq.pstep).and.(bench.eq.0)) THEN
               timep = 0
               pind = pind+1
               CALL particle%write_pstate(time,fluid,field_nxt,part_nxt)
            ENDIF
         endif

! Every 'cstep' steps writes global quantities
         IF ((timec.eq.cstep).and.(bench.eq.0)) THEN
            timec = 0
            CALL fluid%global(field_nxt, force, t)
         ENDIF

! Every 'sstep' steps writes spectra
         IF ((times.eq.sstep).and.(bench.eq.0)) THEN
            times = 0
            sind = sind+1
            CALL fluid%spectra(field_nxt)
         ENDIF

! Time evolution
         CALL update_forcing(forcemethod,fluid,force)
         IF (.not.fstatic) CALL GState_update_to(force)
         gdev_active = .TRUE.
         if (dopart) then
            CALL GState_copy(field,field_nxt)
            part  = part_nxt
            CALL stepper%gstep(time, field, part, force, dt, field_nxt, part_nxt)
         else
            CALL GState_copy(field,field_nxt)
            CALL stepper%gstep(time, field, force, dt, field_nxt)
         endif
         gdev_active = .FALSE.
         timet = timet+1; timep = timep+1; timec = timec+1; times = times+1
      END DO RK

! Finishes and writes the benchmark results
      IF (bench.gt.0) THEN
         CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
         CALL GTStop(ihcpu1); CALL GTStop(ihomp1); CALL GTStop(ihwtm1)
      ENDIF
      CALL GTBenchReport(bench,myrank,ini,step,ihcpu1,ihomp1,ihwtm1,             &
               ihcpu2,ihomp2,ihwtm2)
      IF (dopart) THEN
        rbal = rbal + GetLoadBal(particle) ! Get load balancing
        CALL GTBenchReportParticles(bench, myrank, ini, step, maxparts, rbal,    &
             [ GetTime(particle,GPTIME_STEP),   GetTime(particle,GPTIME_COMM),   &
               GetTime(particle,GPTIME_SPLINE), GetTime(particle,GPTIME_TRANSP), &
               GetTime(particle,GPTIME_DATAEX), GetTime(particle,GPTIME_INTERP), &
               GetTime(particle,GPTIME_PUPDATE) ], 'gpbenchmark.txt')
      ENDIF

! End of main
      CALL GTFree(ihcpu1)
      CALL GTFree(ihomp1)
      CALL GTFree(ihwtm1)
      CALL GTFree(ihcpu2)
      CALL GTFree(ihomp2)
      CALL GTFree(ihwtm2)
      CALL MPI_FINALIZE(ierr)
      CALL fftp3d_destroy_plan(plancr)
      CALL fftp3d_destroy_plan(planrc)
      
      END PROGRAM MAIN
