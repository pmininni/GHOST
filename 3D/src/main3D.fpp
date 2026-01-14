!=================================================================
      PROGRAM MAIN
!=================================================================
! GHOST code: Geophysical High Order Suite for Turbulence
!
! Numerically integrates several fluid dynamics equations
! in 2 and 3 dimensions with periodic boundary conditions
! and external forcing. A pseudo-spectral method is used to
! compute spatial derivatives, while different time steppings
! can be usedd to evolve the system in the time domain.
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
! Mininni PD, Gomez DO, Mahajan SM; Astrophys. J. 619, 1019 (2005)
! Gomez DO, Mininni PD, Dmitruk P; Phys. Scripta T116, 123 (2005)
! Gomez DO, Mininni PD, Dmitruk P; Adv. Sp. Res. 35, 899 (2005)
!=================================================================

!
! Definitions for conditional compilation
#include "ghost3D.h"

!
! Modules

      USE fprecision
      USE commtypes
      USE mpivars
      USE filefmt
      USE iovar
      USE grid
      USE fft
      USE ali
      USE var
      USE kes
      USE order
      USE random
      USE threads
      USE offloading
      USE boxsize
      USE gtimer
      USE fftplans
      USE dns
      USE equation_factory
!     USE class_GPart

      IMPLICIT NONE

!
! Arrays for the fields and external forcings

      TYPE(GState), ALLOCATABLE, TARGET :: field(:),field_nxt(:),force(:)
      TYPE(GWorkspace)                  :: workspace
      CLASS(EquationBase), ALLOCATABLE  :: pde

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Hack to initialize fields, should be removed later
#ifdef VELOC_
      COMPLEX(KIND=GP), POINTER, DIMENSION (:,:,:) :: vx,vy,vz
      COMPLEX(KIND=GP), POINTER, DIMENSION (:,:,:) :: fx,fy,fz
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! Temporal data storage arrays
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: C1,C2,C3
      REAL(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:)    :: R1,R2,R3

!
! Auxiliary variables

      COMPLEX(KIND=GP) :: cdump,jdump
      COMPLEX(KIND=GP) :: cdumq,jdumq
      COMPLEX(KIND=GP) :: cdumr,jdumr
      DOUBLE PRECISION :: tmp,tmq,tmr
      DOUBLE PRECISION :: eps,epm

      REAL(KIND=GP)    :: dt,nu,mu
      REAL(KIND=GP)    :: kup,kdn
      REAL(KIND=GP)    :: rmp,rmq,rms
      REAL(KIND=GP)    :: rmt,rm1,rm2
      REAL(KIND=GP)    :: dump
      REAL(KIND=GP)    :: stat
      REAL(KIND=GP)    :: f0,u0
      REAL(KIND=GP)    :: phase,ampl,cort,time
      REAL(KIND=GP)    :: fparam0,fparam1,fparam2,fparam3,fparam4
      REAL(KIND=GP)    :: fparam5,fparam6,fparam7,fparam8,fparam9
      REAL(KIND=GP)    :: vparam0,vparam1,vparam2,vparam3,vparam4
      REAL(KIND=GP)    :: vparam5,vparam6,vparam7,vparam8,vparam9

      INTEGER :: idevice, iret, ncuda, ngcuda, ppn
      INTEGER :: ini,step
      INTEGER :: num_components
      INTEGER :: tstep,cstep
      INTEGER :: sstep,fstep
      INTEGER :: bench,trans
      INTEGER :: outs,mean
      INTEGER :: seed,rand
      INTEGER :: mult
      INTEGER :: t,o
      INTEGER :: i,j,k
      INTEGER :: ki,kj,kk
      INTEGER :: pind,tind,sind
      INTEGER :: timet,timec
      INTEGER :: times,timef
      INTEGER :: timep,pstep,lgmult
      INTEGER :: ihcpu1,ihcpu2
      INTEGER :: ihomp1,ihomp2
      INTEGER :: ihwtm1,ihwtm2
      TYPE(IOPLAN)          :: planio
      CHARACTER(len=100)    :: odir,idir
      LOGICAL               :: bbenchexist

!
! Namelists for the input files

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Status, parameter should be read by some init that could be part of a module (grid?)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      NAMELIST / status / idir,odir,stat,mult,bench,outs,mean,trans,iswap
      NAMELIST / parameter / dt,step,tstep,sstep,cstep,rand,cort,seed
#if defined(VELOC_) || defined(ADVECT_)
      NAMELIST / velocity / f0,u0,kdn,kup,nu,fparam0,fparam1,fparam2
      NAMELIST / velocity / fparam3,fparam4,fparam5,fparam6,fparam7
      NAMELIST / velocity / fparam8,fparam9,vparam0,vparam1,vparam2
      NAMELIST / velocity / vparam3,vparam4,vparam5,vparam6,vparam7
      NAMELIST / velocity / vparam8,vparam9
#endif

!
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

! Initialization of I/O
      CALL range(1,nx/2+1,nprocs,myrank,ista,iend)
      CALL range(1,nz,nprocs,myrank,ksta,kend)
      CALL io_init(myrank,(/nx,ny,nz/),ksta,kend,planio)

! Now we can initialize the PDE method
     pde = init_pdes_from_file('parameter.inp')
     CALL workspace%initialize_pool(NUMTMPREAL,NUMTMPCOMP)
     CALL pde%Solver_ctor('parameter.inp',workspace,NUMFIELDS)
     num_components = pde%state_size()
     CALL GState_alloc(field    , num_components)
     CALL GState_alloc(field_nxt, num_components)
     CALL GState_alloc(force    , num_components)

! Initialization of the numerical domain
     CALL box_init('parameter.inp')

     ALLOCATE( C1(nz,ny,ista:iend), C2(nz,ny,ista:iend), C3(nz,ny,ista:iend) )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Temporary hack to initialize fields
#ifdef VELOC_
      vx => field(1)%ccomp
      vy => field(2)%ccomp
      vz => field(3)%ccomp
      fx => force(1)%ccomp
      fy => force(2)%ccomp
      fz => force(3)%ccomp
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ALLOCATE( R1(nx,ny,ksta:kend), R2(nx,ny,ksta:kend), R3(nx,ny,ksta:kend) )

!
! The following lines read the file 'parameter.inp'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Status, parameter and maybe the resolution (in boxparams) should
!!! be read by some init that could be part of a module (grid?)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Reads general configuration flags from the namelist
! 'status' on the external file 'parameter.inp'
!     idir : directory for unformatted input
!     odir : directory for unformatted output
!     stat : = 0 starts a new run
!            OR  gives the number of the file used to continue a run
!     mult : time step multiplier
!     bench: = 0 production run
!            = 1 benchmark run (no I/O)
!            = 2 higher level benchmark run (+time to create plans)
!     outs : = 0 writes velocity [and vect potential (MAGFIELD_)]
!            = 1 writes vorticity [and magnetic field (MAGFIELD_)]
!            = 2 writes current density (MAGFIELD_)
!     mean : = 0 skips mean field computation
!            = 1 performs mean field computation
!     trans: = 0 skips energy transfer computation
!            = 1 performs energy transfer computation
!     iswap  = 0 does nothing to restart binary data
!            = 1 does a byte swap of restart binary data

      IF (myrank.eq.0) THEN
         OPEN(1,file='parameter.inp',status='unknown',form="formatted")
         READ(1,NML=status)
         CLOSE(1)
      ENDIF
      CALL MPI_BCAST(idir,100,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(odir,100,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(stat,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(mult,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(bench,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(outs,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(mean,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(trans,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(iswap,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

! Reads parameters that will be used to control the
! time integration from the namelist 'parameter' on
! the external file 'parameter.inp'
!     dt   : time step size
!     step : total number of time steps to compute
!     tstep: number of steps between binary output
!     sstep: number of steps between power spectrum output
!     cstep: number of steps between output of global quantities
!     rand : = 0 constant force
!            = 1 random phases
!            = 2 slowly varying random phases (only for the velocity and
!                magnetic forcings)
!            = 3 user-defined forcing scheme
!     cort : time correlation of the external forcing
!     seed : seed for the random number generator

      IF (myrank.eq.0) THEN
         OPEN(1,file='parameter.inp',status='unknown',form="formatted")
         READ(1,NML=parameter)
         CLOSE(1)
         dt = dt/real(mult,kind=GP)
         step = step*mult
         tstep = tstep*mult
         sstep = sstep*mult
         cstep = cstep*mult
         fstep = int(cort/dt)
      ENDIF
      CALL MPI_BCAST(dt,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(step,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(tstep,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(sstep,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(cstep,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(fstep,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(rand,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(seed,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      Lx  = 1.0_GP
      Ly  = 1.0_GP
      Lz  = 1.0_GP
      Dkx = 1.0_GP
      Dky = 1.0_GP
      Dkz = 1.0_GP
      Dkk = 0.0_GP

#if defined(VELOC_) || defined(ADVECT_)
! Reads parameters for the velocity field from the
! namelist 'velocity' on the external file 'parameter.inp'
!     f0   : amplitude of the mechanical forcing
!     u0   : amplitude of the initial velocity field
!     kdn  : minimum wave number in v/mechanical forcing
!     kup  : maximum wave number in v/mechanical forcing
!     nu   : kinematic viscosity
!     fparam0-9 : ten real numbers to control properties of
!            the mechanical forcing
!     vparam0-9 : ten real numbers to control properties of
!            the initial conditions for the velocity field

      IF (myrank.eq.0) THEN
         OPEN(1,file='parameter.inp',status='unknown',form="formatted")
         READ(1,NML=velocity)
         CLOSE(1)
      ENDIF
      CALL MPI_BCAST(f0,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(u0,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(kdn,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(kup,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(nu,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(fparam0,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(fparam1,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(fparam2,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(fparam3,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(fparam4,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(fparam5,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(fparam6,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(fparam7,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(fparam8,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(fparam9,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(vparam0,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(vparam1,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(vparam2,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(vparam3,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(vparam4,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(vparam5,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(vparam6,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(vparam7,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(vparam8,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(vparam9,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
#endif

! Initializes arrays and constants for the pseudospectral method

! Some constants for the FFT
!     kmax: maximum truncation for dealiasing
!     tiny: minimum truncation for dealiasing

      kmax =     1.0_GP/9.0_GP
      nmax =     int(max(nx*Dkx,ny*Dky,nz*Dkz)/Dkk)
      nmaxperp = int(max(nx*Dkx,ny*Dky)/Dkk)
#ifndef DEF_ARBSIZE_
      IF ( .NOT. anis )  kmax = kmax*real(nx,kind=GP)**2
#endif
#ifdef EDQNM_
      kmax = (real(n,kind=GP)/2.0_GP-0.5_GP)**2
#endif
      tiny  = min(1e-5_GP ,.1_GP/(real(nmax,kind=GP)**2))
      tinyf = min(1e-15_GP,.1_GP/(real(nmax,kind=GP)**2))

! Builds arrays with the wavenumbers and the
! square wavenumbers. At the end, kx, ky, and kz
! have wavenumbers with dimensions, kk2 has the
! squared wavenumbers with dimensions, and kn2 has
! the dimensionless and normalized squared
! wavenumbers used for dealiasing.

      DO i = 1,nx/2
         kx(i) = real(i-1,kind=GP)
         kx(i+nx/2) = real(i-nx/2-1,kind=GP)
      END DO
      IF (nx.eq.1) THEN
         kx(1) = 0.0_GP
      END IF
      DO j = 1,ny/2
         ky(j) = real(j-1,kind=GP)
         ky(j+ny/2) = real(j-ny/2-1,kind=GP)
      END DO
      IF (ny.eq.1) THEN
         ky(1) = 0.0_GP
      ENDIF
      DO k = 1,nz/2
         kz(k) = real(k-1,kind=GP)
         kz(k+nz/2) = real(k-nz/2-1,kind=GP)
      END DO
      IF ( anis ) THEN
         rmp = 1.0_GP/real(nx,kind=GP)**2
         rmq = 1.0_GP/real(ny,kind=GP)**2
         rms = 1.0_GP/real(nz,kind=GP)**2
      ELSE
         rmp = 1.0_GP
	 rmq = 1.0_GP
	 rms = 1.0_GP
      ENDIF
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
            DO k = 1,nz
               kn2(k,j,i) = rmp*kx(i)**2+rmq*ky(j)**2+rms*kz(k)**2
            END DO
         END DO
      END DO
#ifdef DEF_ARBSIZE_
      kx = kx*Dkx
      ky = ky*Dky
      kz = kz*Dkz
#endif
      IF ( anis ) THEN
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
            DO j = 1,ny
               DO k = 1,nz
                  kk2(k,j,i) = kx(i)**2+ky(j)**2+kz(k)**2
               END DO
            END DO
         END DO
      ENDIF

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

!
! Sets the external forcing

#ifdef VELOC_
      INCLUDE 'initialfv.f90'           ! mechanical forcing
#endif

! If stat=0 we start a new run.
! Generates initial conditions for the fields and particles.
 IC : IF (stat.eq.0) THEN

      ini  = 1
      sind = 0                          ! index for the spectrum
      tind = 0                          ! index for the binaries
      pind = 0                          ! index for the particles
      timet = tstep
      timec = cstep
      times = sstep
      timep = pstep
#if defined(VELOC_) || defined (ADVECT_)
      INCLUDE 'initialv.f90'            ! initial velocity
#endif

      ELSE

! If stat.ne.0 a previous run is continued

      ini = int((stat-1)*tstep) + 1
      tind = int(stat)
      sind = int(real(ini,kind=GP)/real(sstep,kind=GP)+1)
      pind = int((stat-1)*lgmult+1)
      WRITE(ext, fmtext) tind
      timet = 0
      timep = 0
      times = int(modulo(float(ini-1),float(sstep)))
      timec = int(modulo(float(ini-1),float(cstep)))
      timef = int(modulo(float(ini-1),float(fstep)))

      CALL io_read(1,idir,'vx',ext,planio,R1)
      CALL io_read(1,idir,'vy',ext,planio,R2)
      CALL io_read(1,idir,'vz',ext,planio,R3)
      CALL fftp3d_real_to_complex(planrc,R1,vx,MPI_COMM_WORLD)
      CALL fftp3d_real_to_complex(planrc,R2,vy,MPI_COMM_WORLD)
      CALL fftp3d_real_to_complex(planrc,R3,vz,MPI_COMM_WORLD)

      ENDIF IC

!
! Time integration scheme starts here.
! Does ord iterations of Runge-Kutta. If
! we are doing a benchmark, we measure
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

! Every 'tstep' steps, stores the fields
! in binary files

         IF ((timet.eq.tstep).and.(bench.eq.0)) THEN
            timet = 0
            tind = tind+1
            WRITE(ext, fmtext) tind
#ifdef VELOC_
            rmp = 1.0_GP/ &
	          (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))
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
            CALL fftp3d_complex_to_real(plancr,C1,R1,MPI_COMM_WORLD)
            CALL fftp3d_complex_to_real(plancr,C2,R2,MPI_COMM_WORLD)
            CALL fftp3d_complex_to_real(plancr,C3,R3,MPI_COMM_WORLD)
            CALL io_write(1,odir,'vx',ext,planio,R1)
            CALL io_write(1,odir,'vy',ext,planio,R2)
            CALL io_write(1,odir,'vz',ext,planio,R3)
#endif
         ENDIF

! Every 'cstep' steps, generates external files 
! with global quantities. If mean=1 also updates 
! the mean fields.

         IF ((timec.eq.cstep).and.(bench.eq.0)) THEN
            timec = 0
            INCLUDE GLOBALOUTPUT_
         ENDIF

! Every 'sstep' steps, generates external files 
! with the power spectrum.

         IF ((times.eq.sstep).and.(bench.eq.0)) THEN
            times = 0
            sind = sind+1
            WRITE(ext, fmtext) sind
            INCLUDE SPECTROUTPUT_
         ENDIF

! Runge-Kutta
         CALL pde%timestep(time, field, force, dt, field_nxt)
         field = field_nxt

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
      DEALLOCATE( R1,R2,R3 )
      DEALLOCATE( C1,C2,C3 )
      DEALLOCATE( kx,ky,kz )
      IF ( anis ) THEN
         DEALLOCATE( kk2 )
      ELSE
         NULLIFY( kk2 )
      ENDIF
      DEALLOCATE( kn2 )
      
      END PROGRAM MAIN
