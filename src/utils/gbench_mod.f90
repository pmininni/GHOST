!=================================================================
! GHOST suite benchmark auxiliary module
! Provides subroutines to write timer outputs from the main code
!
! 2026 Pablo Mininni
! Departamento de Fisica, FCEN, UBA
!=================================================================
module gbench
  use gtimer
  implicit none

  ! Time handles for benchmarks in the main code
  REAL          :: rbal
  INTEGER       :: ihcpu1,ihcpu2
  INTEGER       :: ihomp1,ihomp2
  INTEGER       :: ihwtm1,ihwtm2

contains

  subroutine GTBenchInit(benchlv,hcpu1,homp1,hwtm1,hcpu2,homp2,hwtm2)
!-----------------------------------------------------------------
! DESCRIPTION: Initialises all benchmark timer handles and resets
! the global module-level timing accumulators (ffttime, tratime,
! comtime, and tottime). Called once in main before the
! time-integration loop. Caller is responsible for GTStart calls
! after this call.
!
! ARGUMENTS:
!   benchlv           : benchmark level (0=off, 1=loop only, 2=+FFTW plan)
!   hcpu1,homp1,hwtm1 : handles for the main integration loop
!   hcpu2,homp2,hwtm2 : handles for the FFTW plan (benchlv=2 only)
!-----------------------------------------------------------------
    USE commtypes
    USE mpivars
    USE fft,       ONLY: ffttime, tratime, comtime, tottime
    IMPLICIT NONE 
    INTEGER, INTENT(IN)    :: benchlv
    INTEGER, INTENT(INOUT) :: hcpu1, homp1, hwtm1
    INTEGER, INTENT(INOUT) :: hcpu2, homp2, hwtm2

    IF (benchlv .EQ. 2) THEN
       CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
       CALL GTInitHandle(hcpu2, GT_CPUTIME)
       CALL GTInitHandle(homp2, GT_OMPTIME)
       CALL GTInitHandle(hwtm2, GT_WTIME)
    ENDIF
    IF (benchlv .GE. 1) THEN
       CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
       CALL GTInitHandle(hcpu1, GT_CPUTIME)
       CALL GTInitHandle(homp1, GT_OMPTIME)
       CALL GTInitHandle(hwtm1, GT_WTIME)
       ffttime = 0.D0
       tratime = 0.D0
       comtime = 0.D0
       tottime = 0.D0
       rbal    = 0.0
    ENDIF
    end subroutine GTBenchInit
 

!-----------------------------------------------------------------
  subroutine GTBenchReport(benchlv, myrank, ini, step,       &
                           hcpu1, homp1, hwtm1,              &
                           hcpu2, homp2, hwtm2)
!-----------------------------------------------------------------
! DESCRIPTION: Writes the fluid benchmark results to 'benchmark.txt'
! and to stdout, then frees all six handles.  Called once in main
! after the time-integration loop has finished and GTStop has been
! called by the caller.
!
! ARGUMENTS:
!   benchlv           : benchmark level
!   myrank            : MPI rank (only rank 0 writes)
!   ini, step         : first and last time-step indices
!   hcpu1,homp1,hwtm1 : handles for the main loop
!   hcpu2,homp2,hwtm2 : handles for the FFTW plan (benchlv=2 only)
!-----------------------------------------------------------------
    USE fft,     ONLY: ffttime, tratime, comtime, tottime
    USE grid,    ONLY: nx, ny, nz
    USE mpivars, ONLY: nprocs
    USE threads, ONLY: nth
    IMPLICIT NONE
    INTEGER, INTENT(IN)    :: benchlv, myrank, ini, step
    INTEGER, INTENT(INOUT) :: hcpu1, homp1, hwtm1
    INTEGER, INTENT(INOUT) :: hcpu2, homp2, hwtm2
    LOGICAL :: bbenchexist
    INTEGER :: nsteps
 
    IF (benchlv .EQ. 0) RETURN
 
    nsteps = step - ini + 1
    inquire(file='benchmark.txt', exist=bbenchexist)
    IF (myrank .EQ. 0) THEN
       OPEN(1, file='benchmark.txt', position='append')
       IF (.NOT. bbenchexist) THEN
          WRITE(1,*) &
           '# nx ny nz nsteps nprocs nth' //                 &
           ' TCPU TOMP TWTIME TFFT TTRA TCOM TTOT'
       ENDIF
       WRITE(1,*) nx, ny, nz, nsteps, nprocs, nth,           &
                  GTGetTime(hcpu1)/nsteps,                   &
                  GTGetTime(homp1)/nsteps,                   &
                  GTGetTime(hwtm1)/nsteps,                   &
                  ffttime/nsteps, tratime/nsteps,            &
                  comtime/nsteps, tottime/nsteps
       WRITE(*,*) 'wtime=', GTGetTime(hwtm1)/nsteps,         &
                  ' fft=',    ffttime/nsteps,                &
                  ' transp=', tratime/nsteps,                &
                  ' comm=',   comtime/nsteps,                &
                  ' mem=',    0.0,                           &
                  ' ttot=',   tottime/nsteps
       IF (benchlv .EQ. 2) THEN
          WRITE(1,*) 'FFTW: Create_plan = ',                 &
                     GTGetTime(hcpu2)/nsteps,                &
                     GTGetTime(homp2)/nsteps,                &
                     GTGetTime(hwtm2)/nsteps
          ENDIF
          CLOSE(1)
      ENDIF
      CALL GTFree(hcpu1); CALL GTFree(homp1); CALL GTFree(hwtm1)
      CALL GTFree(hcpu2); CALL GTFree(homp2); CALL GTFree(hwtm2)
    end subroutine GTBenchReport
 
 
!-----------------------------------------------------------------
    subroutine GTBenchReportParticles(benchlv, myrank, ini,  &
                                      step, nparts, rbal,    &
                                      ptimes, filename)
!-----------------------------------------------------------------
! DESCRIPTION: Writes particle benchmark results to a named file
! (e.g. 'gpbenchmark.txt').  Intentionally decoupled from GTBenchReport
! so it can be called for different particles independently, and without
! introducing a dependency on the particle class hierarchy in gtimer.
! The caller extracts timing values from the particle object via
! particle%GetTime(GPTIME_*) and packs them into ptimes(:) before
! calling this routine.  The GPTIME_* index constants are defined in
! the particle module; the array is indexed [1:7] here in the same
! order they appear in the benchmark file header.
!
! ARGUMENTS:
!   myrank   : MPI rank (only rank 0 writes)
!   ini,step : first and last time-step indices
!   nparts   : number of particles (maxparts or partpcell*nx*ny*nz)
!   rbal     : cumulative load-balance sum over all steps (raw,
!              as accumulated with rbal = rbal + particle%GetLoadBal()
!              in the RK loop); divided by nsteps here
!   ptimes   : raw (un-normalised) particle timer values, packed as:
!                ptimes(1) = GetTime(GPTIME_STEP)
!                ptimes(2) = GetTime(GPTIME_COMM)
!                ptimes(3) = GetTime(GPTIME_SPLINE)
!                ptimes(4) = GetTime(GPTIME_TRANSP)
!                ptimes(5) = GetTime(GPTIME_DATAEX)
!                ptimes(6) = GetTime(GPTIME_INTERP)
!                ptimes(7) = GetTime(GPTIME_PUPDATE)
!   filename : output file name
!-----------------------------------------------------------------
    USE grid,    ONLY: nx, ny, nz
    USE mpivars, ONLY: nprocs
    USE threads, ONLY: nth
    IMPLICIT NONE
    INTEGER,          INTENT(IN) :: benchlv, myrank
    INTEGER,          INTENT(IN) :: ini, step, nparts
    REAL,             INTENT(IN) :: rbal
    DOUBLE PRECISION, INTENT(IN) :: ptimes(7)
    CHARACTER(LEN=*), INTENT(IN) :: filename
    LOGICAL :: bbenchexist
    INTEGER :: nsteps

    IF (benchlv .EQ. 0) RETURN

    nsteps = step - ini + 1
    IF (myrank .EQ. 0) THEN
       inquire(file=trim(filename), exist=bbenchexist)
       OPEN(1, file=trim(filename), position='append')
       IF (.NOT. bbenchexist) THEN
          WRITE(1,*) &
           '# nx ny nz nparts rbal nsteps nprocs nth' //     &
           ' TRK TCOMM TSPL TTRANS TDEX TINT TUPD'
       ENDIF
       WRITE(1,*) nx, ny, nz, nparts, rbal/nsteps,           &
                  nsteps, nprocs, nth,                       &
                  ptimes(1)/nsteps,                          &
                  ptimes(2)/nsteps,                          &
                  ptimes(3)/nsteps,                          &
                  ptimes(4)/nsteps,                          &
                  ptimes(5)/nsteps,                          &
                  ptimes(6)/nsteps,                          &
                  ptimes(7)/nsteps
       CLOSE(1)
    ENDIF
    end subroutine GTBenchReportParticles

end module gbench
    
