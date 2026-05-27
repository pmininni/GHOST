!=================================================================
! gpython_lib.fpp - Python Interface for GHOST
! 2026 Patricio Clark di Leoni
!=================================================================
MODULE gpython_lib
  USE, INTRINSIC :: iso_c_binding
  USE commtypes
  USE filefmt
  USE iovar
  USE fft
  USE threads
  USE offloading
  USE boxsize
  USE status
  USE pstatus
  USE ic_factory
  USE force_factory
  USE equation_factory
  USE particle_factory
  USE icp_factory
  USE stepper_factory
  IMPLICIT NONE

  ! --- Global State (Persists between Python calls) ---
  TYPE   (GStateComp), ALLOCATABLE, TARGET :: field(:), field_nxt(:), force(:)
  TYPE   (GWorkspace)                      :: workspace
  TYPE       (IOPLAN)                      :: planio
  CLASS(EquationBase), ALLOCATABLE         :: fluid
  CLASS(GStepperBase), ALLOCATABLE         :: stepper
  CLASS     (icChain), ALLOCATABLE         :: iclist(:)
  CLASS  (forceChain), ALLOCATABLE         :: forcemethod(:)
  REAL      (KIND=GP) :: time
  INTEGER             :: current_step = 0
  INTEGER             :: num_components

CONTAINS

  !=================================================================
  ! Initialize GHOST
  !=================================================================
  SUBROUTINE c_ghost_init() BIND(C, name='ghost_init')
    CHARACTER(len=64) :: file = 'parameter.inp'

    ! MPI Init
    CALL MPI_INIT_THREAD(MPI_THREAD_FUNNELED, provided, ierr)
    CALL MPI_COMM_SIZE  (MPI_COMM_WORLD, nprocs, ierr)
    CALL MPI_COMM_RANK  (MPI_COMM_WORLD, myrank, ierr)

    ! Grid & Status
    CALL grid_init  (trim(file))
    CALL status_init(trim(file))

    ! I/O Setup
    CALL range(1, nx/2+1, nprocs, myrank, ista, iend)
    CALL range(1, nz    , nprocs, myrank, ksta, kend)
    CALL io_init(myrank, (/nx,ny,nz/), ksta, kend, planio)

    ! PDE & Arrays
    fluid = init_pdes_from_file(trim(file))
    CALL workspace%initialize_pool(NUMTMPREAL, NUMTMPCOMP)
    CALL fluid%Solver_ctor(trim(file), workspace, planio)

    num_components = fluid%state_size()
    CALL GState_alloc(field    , num_components)
    CALL GState_alloc(field_nxt, num_components)
    CALL GState_alloc(force    , num_components)

    iclist      = init_ic_from_file     (trim(file))
    forcemethod = init_forcing_from_file(trim(file),workspace)

    ! Initialization of the numerical domain
    CALL box_init(trim(file))

    ! FFT Init (Crucial: Must match main.fpp logic)
    CALL fftp3d_create_plan(planrc, (/nx,ny,nz/), FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE)
    CALL fftp3d_create_plan(plancr, (/nx,ny,nz/), FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE)

    ! Initialize fluid states and the stepper
    CALL init_allstates(iclist,fluid,field)
    field_nxt = field
    CALL init_forcing(forcemethod,fluid,force)
    stepper = build_stepper_from_file(trim(file),workspace,fluid)
  END SUBROUTINE c_ghost_init

  !=================================================================
  ! Run steps
  !=================================================================
  SUBROUTINE c_ghost_run(num_steps) BIND(C, name='ghost_run')
    INTEGER(C_INT), VALUE :: num_steps
    INTEGER :: i
    REAL(KIND=GP) :: dt_val

    DO i = 1, num_steps
       time = (current_step-1)*dt
       CALL update_forcing(forcemethod,fluid,force)
       field = field_nxt
       CALL stepper%gstep(time, field, force, dt, field_nxt)
       current_step = current_step + 1
    END DO
  END SUBROUTINE c_ghost_run

  !=================================================================
  ! Finalize
  !=================================================================
  SUBROUTINE c_ghost_finalize() BIND(C, name='ghost_finalize')
    CALL fftp3d_destroy_plan(plancr)
    CALL fftp3d_destroy_plan(planrc)
    CALL MPI_FINALIZE(ierr)
  END SUBROUTINE c_ghost_finalize

END MODULE gpython_lib
