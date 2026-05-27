!=================================================================
! gpython_lib.fpp - Python Interface for GHOST
! For the moment, it provides access to the fluid solvers only.
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
  TYPE    (GStateComp), ALLOCATABLE, TARGET :: field(:),field_nxt(:),force(:)
  TYPE(GStateRealComp), ALLOCATABLE, TARGET :: realfield(:)
  TYPE    (GWorkspace)                      :: workspace
  TYPE        (ioplan)                      :: planio
  CLASS (EquationBase), ALLOCATABLE         :: fluid
  CLASS (GStepperBase), ALLOCATABLE         :: stepper
  CLASS      (icChain), ALLOCATABLE         :: iclist(:)
  CLASS   (forceChain), ALLOCATABLE         :: forcemethod(:)
  REAL       (KIND=GP) :: time
  INTEGER              :: current_step = 0
  INTEGER              :: num_components

CONTAINS

  !=================================================================
  ! Initialize GHOST
  !=================================================================
  SUBROUTINE ghost_init(c_file) BIND(C)
    TYPE           (C_PTR), VALUE   :: c_file
    CHARACTER(kind=C_CHAR), POINTER :: tmp(:)
    CHARACTER    (len=256)          :: file
    INTEGER                         :: i

    IF (.NOT. C_ASSOCIATED(c_file)) THEN
      file = 'parameter.inp'
    ELSE
      CALL c_f_pointer(c_file, tmp, [256])
      file = ''
      DO i = 1, SIZE(tmp)
        IF (tmp(i) == C_NULL_CHAR) EXIT
        file(i:i) = tmp(i)
      END DO
    END IF

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
    CALL GState_alloc    (field    , num_components)
    CALL GState_alloc    (field_nxt, num_components)
    CALL GState_alloc    (force    , num_components)
    CALL GStateReal_alloc(realfield, num_components)

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
  END SUBROUTINE ghost_init

  !=================================================================
  ! Run steps
  !=================================================================
  SUBROUTINE ghost_run(num_steps) BIND(C)
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
  END SUBROUTINE ghost_run

  !=================================================================
  ! Finalize
  !=================================================================
  SUBROUTINE ghost_finalize() BIND(C)
    CALL fftp3d_destroy_plan(plancr)
    CALL fftp3d_destroy_plan(planrc)
    CALL MPI_FINALIZE(ierr)
  END SUBROUTINE ghost_finalize

  !=================================================================
  ! Helpers
  !=================================================================
  ! Returns the pointer to the latest field component at position num_field
  FUNCTION ghost_get_complex_field(num_field) RESULT(ptr) BIND(C)
    INTEGER(C_INT), VALUE :: num_field
    TYPE   (C_PTR)        :: ptr
    ptr = C_LOC(field_nxt(num_field)%ccomp)
  END FUNCTION ghost_get_complex_field
  
  ! Returns the pointer to the forcing component at position num_field
  FUNCTION ghost_get_complex_forcing(num_field) RESULT(ptr) BIND(C)
    INTEGER(C_INT), VALUE :: num_field
    TYPE   (C_PTR)        :: ptr
    ptr = C_LOC(force(num_field)%ccomp)
  END FUNCTION ghost_get_complex_forcing

  ! Returns the pointer to the Fourier-transformed (real) latest field
  ! component at position num_field.
  FUNCTION ghost_get_real_field(num_field) RESULT(ptr) BIND(C)
    INTEGER(C_INT)  , VALUE                     :: num_field
    TYPE   (C_PTR)                              :: ptr
    REAL   (kind=GP)                            :: rmp
    COMPLEX(kind=GP), POINTER, DIMENSION(:,:,:) :: ctmp
    LOGICAL                                     :: bret
    CALL workspace%get_complex_tmp(ctmp,bret)
    rmp = 1.0_GP/(real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))
    ctmp = field_nxt(num_field)%ccomp * rmp
    CALL fftp3d_complex_to_real(plancr,ctmp,realfield(num_field)%rcomp,MPI_COMM_WORLD)
    ptr = C_LOC(realfield(num_field)%rcomp)
  END FUNCTION ghost_get_real_field
  
  ! Return complex array sizes 
  FUNCTION ghost_get_complex_size(num_comp) RESULT(sz) BIND(C)
    INTEGER(C_INT), VALUE :: num_comp
    INTEGER(C_INT)        :: sz
    sz = SIZE(field(1)%ccomp,num_comp)
  END FUNCTION ghost_get_complex_size

  ! Return real array sizes 
  FUNCTION ghost_get_real_size(num_comp) RESULT(sz) BIND(C)
    INTEGER(C_INT), VALUE :: num_comp
    INTEGER(C_INT)        :: sz
    sz = SIZE(realfield(1)%rcomp,num_comp)
  END FUNCTION ghost_get_real_size

END MODULE gpython_lib
