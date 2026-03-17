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
      USE offloading
      USE boxsize
      USE status
      USE pstatus
      USE gtimer
      USE ic_factory
      USE force_factory
      USE equation_factory
      USE particle_factory
      USE stepper_factory
      IMPLICIT NONE

! Arrays for the field and particle states, workspace, I/O, and solver classes
      TYPE(GStateComp),    ALLOCATABLE :: field(:),field_nxt(:),force (:)
      TYPE(GPStateComp),   ALLOCATABLE :: part (:),part_nxt (:)
      TYPE(GWorkspace)                 :: workspace
      TYPE(IOPLAN)                     :: planio
      CLASS(EquationBase), ALLOCATABLE :: fluid
      CLASS(ParticleBase), ALLOCATABLE :: particle
      CLASS(GStepperBase), ALLOCATABLE :: stepper
      CLASS(icChain),      ALLOCATABLE :: iclist(:)
      CLASS(forceChain),   ALLOCATABLE :: forcemethod(:)

! Auxiliary variables
      INTEGER :: t
      INTEGER :: num_components
      INTEGER :: idevice, iret, ncuda, ngcuda, ppn
      INTEGER :: ihcpu1,ihcpu2
      INTEGER :: ihomp1,ihomp2
      INTEGER :: ihwtm1,ihwtm2
      LOGICAL :: bbenchexist

! Initialization
! Initializes the MPI and I/O libraries
      CALL MPI_INIT_THREAD(MPI_THREAD_FUNNELED,provided,ierr)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)

! Initializes the grid. This must be done early to have nx, ny, nz.
      CALL grid_init('parameter.inp')

! Initialization of offloading to GPUs using OpenMP (this is independent
! of CUDA initialization in systems with NVIDIA GPUs). GHOST
! assumes the number of MPI jobs in each node is equal to the
! number of GPUs available in the node. The user must ensure this
! condition is fulfilled.
#if defined(DO_HYBRIDoffl)
      CALL init_offload(myrank,numdev,hostdev,targetdev)
#endif

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
      if (dopart) CALL particle%part_ctor('parameter.inp',workspace,part,part_nxt)

! Initializes the FFT library. This must be done at
! this stage as it requires the variable "bench" to
! be properly initialized.
! Use FFTW_ESTIMATE or FFTW_MEASURE in short runs
! Use FFTW_PATIENT or FFTW_EXHAUSTIVE in long runs
      nth = 1
!$    nth = omp_get_max_threads()
#if !defined(DEF_GHOST_CUDA_)
!$    CALL fftp3d_init_threads(ierr)
#endif
      IF (bench.eq.2) THEN
         CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
         CALL GTInitHandle(ihcpu2,GT_CPUTIME)
         CALL GTInitHandle(ihomp2,GT_OMPTIME)
         CALL GTInitHandle(ihwtm2,GT_WTIME)
         CALL GTStart(ihcpu2)
         CALL GTStart(ihomp2)
         CALL GTStart(ihwtm2)
      ENDIF
      CALL fftp3d_create_plan(planrc,(/nx,ny,nz/),FFTW_REAL_TO_COMPLEX, &
                             FFTW_ESTIMATE)
      CALL fftp3d_create_plan(plancr,(/nx,ny,nz/),FFTW_COMPLEX_TO_REAL, &
                             FFTW_ESTIMATE)
      IF (bench.eq.2) THEN
         CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
         CALL GTStop(ihcpu2)
         CALL GTStop(ihomp2)
         CALL GTStop(ihwtm2)
      ENDIF

! Initial states
 IC : IF (stat.eq.0) THEN                 ! If stat=0 we start a new run
        ini  = 1
        sind = 0                          ! index for the spectrum
        tind = 0                          ! index for the binaries
!       pind = 0                          ! index for the particles
        timet = tstep
        timep = pstep
        timec = cstep
        times = sstep
        timef = 0
      ELSE                    ! If stat.ne.0 a previous run is continued
        ini = int((stat-1)*tstep) + 1
        tind = int(stat)
        sind = int(real(ini,kind=GP)/real(sstep,kind=GP)+1)
!       pind = int((stat-1)*lgmult+1)
        timet = 0
        timep = 0
        timec = int(modulo(float(ini-1),float(cstep)))
        times = int(modulo(float(ini-1),float(sstep)))
        timef = int(modulo(float(ini-1),float(fstep)))
      ENDIF IC
      CALL init_allstates(iclist,fluid,field)
      CALL init_forcing(forcemethod,fluid,force)
!     if (dopart) CALL init_particles()

! Sets up the time stepper
      if (dopart) then
        stepper = build_stepper_from_file('parameter.inp',workspace,fluid,particle)
      else
        stepper = build_stepper_from_file('parameter.inp',workspace,fluid)
      endif

! Time integration scheme starts here.
! If we are doing a benchmark, we measure
! cputime before starting.
      IF (bench.eq.1) THEN
         CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
         CALL GTInitHandle(ihcpu1,GT_CPUTIME)
         CALL GTInitHandle(ihomp1,GT_OMPTIME)
         CALL GTInitHandle(ihwtm1,GT_WTIME)
         ffttime  = 0.D00 ! re-inititialize fftp timers
         tratime  = 0.0D0
         comtime  = 0.D00
         tottime  = 0.0D0
#if defined(DEF_GHOST_CUDA_)
         memtime  = 0.0D0
         asstime  = 0.D00
#endif
         CALL GTStart(ihcpu1)
         CALL GTStart(ihomp1)
         CALL GTStart(ihwtm1)
      ENDIF

 RK : DO t = ini,step

! Every 'tstep' steps, stores the fields in binary files
         IF ((timet.eq.tstep).and.(bench.eq.0)) THEN
            timet = 0
            tind = tind+1
	    CALL fluid%write_states(field, planio)
	 ENDIF

! Every 'pstep' steps, stores the particle states
         if (dopart) then
           IF ((timep.eq.pstep).and.(bench.eq.0)) THEN
             timep = 0
	     pind = pind+1
  	     CALL particle%write_pstate(time,fluid,field,part)
	   ENDIF
	 endif

! Every 'cstep' steps writes global quantities
         IF ((timec.eq.cstep).and.(bench.eq.0)) THEN
            timec = 0
	    CALL fluid%global(field, force, t)
         ENDIF

! Every 'sstep' steps writes spectra
         IF ((times.eq.sstep).and.(bench.eq.0)) THEN
            times = 0
            sind = sind+1
            CALL fluid%spectra(field)
         ENDIF

! Time evolution
         CALL update_forcing(forcemethod,fluid,force)
         if (dopart) then
           CALL stepper%gstep(time, field, part, force, dt, field_nxt, part_nxt)
           field = field_nxt
	   part  = part_nxt
	 else
           CALL stepper%gstep(time, field, force, dt, field_nxt)
	   field = field_nxt
	 endif
         timet = timet+1
         timep = timep+1
         timec = timec+1
         times = times+1

      END DO RK

!
! End of Runge-Kutta

! Computes the benchmark

      IF (bench.gt.0) THEN
         CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
         CALL GTStop(ihcpu1)
         CALL GTStop(ihomp1)
         CALL GTStop(ihwtm1)
         inquire( file='benchmark.txt', exist=bbenchexist )
         IF (myrank.eq.0) THEN
            OPEN(1,file='benchmark.txt',position='append')
#if defined(DEF_GHOST_CUDA_)
            IF ( .NOT. bbenchexist ) THEN
               WRITE(1,*) &
	       '# nx ny nz nsteps nprocs nth nstrm TCPU TOMP TWTIME TFFT TTRA TCOM TMEM TASS TTOT'
            ENDIF
            WRITE(1,*) nx,ny,nz,(step-ini+1),nprocs,nth, &
                       nstreams                        , &
                       GTGetTime(ihcpu1)/(step-ini+1)  , &
                       GTGetTime(ihomp1)/(step-ini+1)  , &
                       GTGetTime(ihwtm1)/(step-ini+1)  , &
                       ffttime/(step-ini+1), tratime/(step-ini+1), &
                       comtime/(step-ini+1), memtime/(step-ini+1), &
                       asstime/(step-ini+1), tottime/(step-ini+1)
            WRITE(*,*) 'wtime=', GTGetTime(ihwtm1)/(step-ini+1),   &
	               ' fft=', ffttime/(step-ini+1),    &
		       ' transp=',tratime/(step-ini+1),  &
		       ' comm=',comtime/(step-ini+1),    &
		       ' mem=', memtime/(step-ini+1),    &
		       ' ttot=',tottime/(step-ini+1)
#else
            IF ( .NOT. bbenchexist ) THEN
               WRITE(1,*) &
	       '# nx ny nz nsteps nprocs nth TCPU TOMP TWTIME TFFT TTRA TCOM TTOT'
            ENDIF
            WRITE(1,*) nx,ny,nz,(step-ini+1),nprocs,nth, &
                       GTGetTime(ihcpu1)/(step-ini+1),   &
                       GTGetTime(ihomp1)/(step-ini+1),   &
                       GTGetTime(ihwtm1)/(step-ini+1),   &
                       ffttime/(step-ini+1), tratime/(step-ini+1), &
                       comtime/(step-ini+1), tottime/(step-ini+1)
            WRITE(*,*) 'wtime=', GTGetTime(ihwtm1)/(step-ini+1),   &
	               ' fft=', ffttime/(step-ini+1),    &
		       ' transp=',tratime/(step-ini+1),  &
		       ' comm=',comtime/(step-ini+1),    &
		       ' mem=',0.0, ' ttot=',tottime/(step-ini+1)
#endif
            IF (bench.eq.2) THEN
               WRITE(1,*) 'FFTW: Create_plan = ',      &
                       GTGetTime(ihcpu2)/(step-ini+1), &
                       GTGetTime(ihomp2)/(step-ini+1), &
                       GTGetTime(ihwtm2)/(step-ini+1)
            ENDIF
            CLOSE(1)
         ENDIF
      ENDIF

!
! End of MAIN3D

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
