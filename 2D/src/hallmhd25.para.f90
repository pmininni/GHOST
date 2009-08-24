!=================================================================
      PROGRAM HALLMHD25D
!=================================================================
! HALLMHD25D code
!
! Numerically integrates the incompressible Hall-MHD 
! equations in 2.5 dimensions using the streamfunction 
! and vector potential formulation with an external force. 
! A pseudo-spectral method is used to compute spatial 
! derivatives, while variable order Runge-Kutta method 
! is used to evolve the system in time domain.
! To compile, you need the FFTW library installed on 
! your system. You should link with the FFTP subroutines
! and use the FFTPLANS and MPIVARS modules (see the file 
! 'fftp_mod.f90').
!
! NOTATION: index 'i' is 'x' 
!           index 'j' is 'y'
!
! 2003 Pablo D. Mininni.
!      Department of Physics, 
!      Facultad de Ciencias Exactas y Naturales.
!      Universidad de Buenos Aires.
!      e-mail: mininni@df.uba.ar 
!=================================================================

      USE mpivars
      USE fft
      USE ali
      USE var
      USE kes
      USE grid
      USE random
      IMPLICIT NONE

!
! Integration parameters
!     ord  : order of the Runge-Kutta method used

      INTEGER, PARAMETER :: ord = 2

!
! streamfunction, vector potential, z component 
! of the fields and external force matrixes

      COMPLEX, ALLOCATABLE, DIMENSION (:,:) :: ps,az
      COMPLEX, ALLOCATABLE, DIMENSION (:,:) :: vz,bz
      COMPLEX, ALLOCATABLE, DIMENSION (:,:) :: fk,fz
      COMPLEX, ALLOCATABLE, DIMENSION (:,:) :: mk,mz

!
! Temporal data storing matrixes

      COMPLEX, ALLOCATABLE, DIMENSION (:,:) :: C1,C2,C3
      COMPLEX, ALLOCATABLE, DIMENSION (:,:) :: C4,C5,C6
      COMPLEX, ALLOCATABLE, DIMENSION (:,:) :: C7,C8,C9
      COMPLEX, ALLOCATABLE, DIMENSION (:,:) :: C10,C11
      REAL, ALLOCATABLE, DIMENSION (:,:)    :: R1,R2

!
! Some auxiliary matrixes

      REAL*8 :: tmp,tmq
      REAL*8 :: tmr,tms
      REAL   :: dt,nu,mu
      REAL   :: kup,kdn
      REAL   :: mkup,mkdn
      REAL   :: stat,dump
      REAL   :: f0,z0,m0,n0
      REAL   :: p0,u0,a0,b0
      REAL   :: time,ttime
      REAL   :: cputime1,cputime2
      REAL   :: corr,cort
      REAL   :: phase
      REAL   :: ep

      INTEGER :: ini
      INTEGER :: step
      INTEGER :: tstep
      INTEGER :: cstep
      INTEGER :: sstep
      INTEGER :: fstep
      INTEGER :: bench
      INTEGER :: mult
      INTEGER :: adap
      INTEGER :: seed
      INTEGER :: t,o
      INTEGER :: i,j
      INTEGER :: ki,kj
      INTEGER :: ic,id,iu
      INTEGER :: jc,jd,ju
      INTEGER :: timet,timec
      INTEGER :: times,timef

      CHARACTER     :: c,d,u
      CHARACTER*3   :: node,ext
      CHARACTER*100 :: ldir

!
! Initializes the MPI library

      CALL MPI_INIT(ierr)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
      ic = 48+int(myrank/100)
      id = 48+int(myrank/10)-int(myrank/100)*10
      iu = 48+int(myrank)-int(myrank/10)*10
      c = char(ic)
      d = char(id)
      u = char(iu)
      node = c // d // u

!
! Allocates memory for distributed blocks

      CALL range(1,n/2+1,nprocs,myrank,ista,iend)
      CALL range(1,n,nprocs,myrank,jsta,jend)

      ALLOCATE( R1(n,jsta:jend), R2(n,jsta:jend)  )
      ALLOCATE( C1(n,ista:iend), C2(n,ista:iend)  )
      ALLOCATE( C3(n,ista:iend), C4(n,ista:iend)  )
      ALLOCATE( C5(n,ista:iend), C6(n,ista:iend)  )
      ALLOCATE( C7(n,ista:iend), C8(n,ista:iend)  )
      ALLOCATE( C9(n,ista:iend), C10(n,ista:iend) )
      ALLOCATE( C11(n,ista:iend))
      ALLOCATE( ps(n,ista:iend), vz(n,ista:iend)  )
      ALLOCATE( az(n,ista:iend), bz(n,ista:iend)  ) 
      ALLOCATE( fk(n,ista:iend), fz(n,ista:iend)  )
      ALLOCATE( mk(n,ista:iend), mz(n,ista:iend)  )
      ALLOCATE( ka(n), ka2(n,ista:iend) )

!
! Reads from the external file 'status.txt'
! the status of a previous run (if any)
!     stat : last output of a previous run
!     mult : time step multiplier
!     adap : adaptive time step factor
!     bench: performs a benchmark run

      IF (myrank.eq.0) THEN
         OPEN(1,file='status.txt',status='unknown')
         READ(1,*) stat
         READ(1,*) mult
         READ(1,*) adap
         READ(1,*) bench
         CLOSE(1)
      ENDIF
      CALL MPI_BCAST(stat,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(mult,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(adap,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(bench,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

!
! Reads from the external file 'parameter.txt' the 
! parameters that will be used during the integration
!     dt   : time step size
!     step : total number of time steps to compute
!     tstep: number of steps between I/O writing
!     sstep: number of steps between power spectrum I/O
!     cstep: number of steps between information output
!     f0   : amplitude of the external mechanic force
!     z0   : amplitude of the external mechanic force in z
!     m0   : amplitude of the external magnetic force
!     n0   : amplitude of the external magnetic force in z
!     p0   : amplitude of the initial streamfunction
!     u0   : amplitude of the velocity field in z
!     a0   : amplitude of the initial vector potential
!     b0   : amplitude of the initial magnetic field in z
!     kdn  : minimum wave number in the external force
!     kup  : maximum wave number in the external force
!     mkdn : minimum wave number in the seed
!     mkup : maximum wave number in the seed
!     nu   : kinematic viscosity
!     mu   : magnetic difusivity
!     corr : initial correlation u.b/|u||b|
!     cort : time correlation of the external force
!     seed : seed for random phase initial conditions
!     ldir : local directory for I/O

      IF (myrank.eq.0) THEN
         OPEN(1,file='parameter.txt',status='unknown')
         READ(1,*) dt
         READ(1,*) step
         READ(1,*) tstep
         READ(1,*) sstep
         READ(1,*) cstep
         READ(1,*) f0
         READ(1,*) z0
         READ(1,*) m0
         READ(1,*) n0
         READ(1,*) p0
         READ(1,*) u0
         READ(1,*) a0
         READ(1,*) b0
         READ(1,*) kdn
         READ(1,*) kup
         READ(1,*) mkdn
         READ(1,*) mkup
         READ(1,*) nu
         READ(1,*) mu
         READ(1,*) corr
         READ(1,*) cort
         READ(1,*) seed
         READ(1,'(a100)') ldir
         CLOSE(1)
         dt = dt/float(mult)
         step = step*mult
         tstep = tstep*mult
         sstep = sstep*mult
         cstep = cstep*mult
         fstep = int(cort/dt)
         ttime = step*dt
      ENDIF
      CALL MPI_BCAST(dt,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(step,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(tstep,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(sstep,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(cstep,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(fstep,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(ttime,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(f0,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(z0,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(m0,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(n0,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(p0,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(u0,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(a0,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(b0,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(kdn,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(kup,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(mkdn,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(mkup,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(nu,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(mu,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(corr,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(cort,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(seed,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(ldir,100,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)

!
! Reads from the external file 'hall.txt'
! the value of the Hall parameters
!     ep  : amplitude of the Hall effect

      IF (myrank.eq.0) THEN
         OPEN(1,file='hall.txt',status='unknown')
         READ(1,*) ep
         CLOSE(1)
      ENDIF
      CALL MPI_BCAST(ep,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)

!
! Some numerical constants

      ic = 48
      id = 48
      iu = 48
      jc = 48
      jd = 48
      ju = 48

!
! Some constants for the FFT
!     kmax: maximum truncation for dealiasing
!     tiny: minimum truncation for dealiasing

      kmax = (float(n)/3.)**2
      tiny = 1e-5

!
! Builds the wave number and the square wave 
! number matrixes

      DO i = 1,n/2
         ka(i) = float(i-1)
         ka(i+n/2) = float(i-n/2-1)
      END DO
      DO i = ista,iend
         DO j = 1,n
            ka2(j,i) = ka(i)**2+ka(j)**2
         END DO
      END DO

!
! Initializes the FFT library
! Use FFTW_ESTIMATE in short runs and FFTW_MEASURE 
! in long runs

      CALL fftp2d_create_plan(planrc,n,FFTW_REAL_TO_COMPLEX, &
                             FFTW_MEASURE)
      CALL fftp2d_create_plan(plancr,n,FFTW_COMPLEX_TO_REAL, &
                             FFTW_MEASURE)

!
! Sets the initial conditions. If the first line of the 
! file 'status.txt' has a number equal to zero, the initial 
! velocity and magnetic fields are given by a superposition 
! of modes in Fourier space. In other case, a previous run 
! is continued from the step indicated by the number found 
! in this file.

      IF (stat.eq.0) THEN

         ini = 1
         timet = tstep
         timec = cstep
         times = sstep
         timef = fstep

         IF (ista.eq.1) THEN
            ps(1,1) = 0.
            vz(1,1) = 0.
            DO j = 2,n/2+1
               IF ((ka2(j,1).le.kup**2).and.(ka2(j,1).ge.kdn**2)) THEN
                  phase = 2*pi*randu(seed)
                  ps(j,1) = (COS(phase)+im*SIN(phase))/sqrt(ka2(j,1))
                  ps(n-j+2,1) = conjg(ps(j,1))
                  phase = 2*pi*randu(seed)
                  vz(j,1) = (COS(phase)+im*SIN(phase))/sqrt(ka2(j,1))
                  vz(n-j+2,1) = conjg(vz(j,1))
               ELSE
                  ps(j,1) = 0.
                  ps(n-j+2,1) = 0.
                  vz(j,1) = 0.
                  vz(n-j+2,1) = 0.
               ENDIF
            END DO
            DO j = 1,n
               DO i = 2,iend
                  IF ((ka2(j,i).le.kup**2).and.(ka2(j,i).ge.kdn**2)) THEN
                     phase = 2*pi*randu(seed)
                     ps(j,i) = 2*(COS(phase)+im*SIN(phase))/sqrt(ka2(j,i))
                     phase = 2*pi*randu(seed)
                     vz(j,i) = 2*(COS(phase)+im*SIN(phase))/sqrt(ka2(j,i))
                  ELSE
                     ps(j,i) = 0.
                     vz(j,i) = 0.
                  ENDIF
               END DO
            END DO
         ELSE
             DO j = 1,n
               DO i = ista,iend
                  IF ((ka2(j,i).le.kup**2).and.(ka2(j,i).ge.kdn**2)) THEN
                     phase = 2*pi*randu(seed)
                     ps(j,i) = 2*(COS(phase)+im*SIN(phase))/sqrt(ka2(j,i))
                     phase = 2*pi*randu(seed)
                     vz(j,i) = 2*(COS(phase)+im*SIN(phase))/sqrt(ka2(j,i))
                  ELSE
                     ps(j,i) = 0.
                     vz(j,i) = 0.
                  ENDIF
               END DO
            END DO
         ENDIF
         CALL energy(ps,tmp,1)
         CALL energy(vz,tmq,2)
         CALL MPI_BCAST(tmp,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
         CALL MPI_BCAST(tmq,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
         ps = ps*sqrt(p0)/sqrt(tmp)
         vz = vz*sqrt(u0)/sqrt(tmq)

         IF (ista.eq.1) THEN
            az(1,1) = 0.
            bz(1,1) = 0.
            DO j = 2,n/2+1
               IF ((ka2(j,1).le.mkup**2).and.(ka2(j,1).ge.mkdn**2)) THEN
                  phase = 2*pi*randu(seed)
                  az(j,1) = (COS(phase)+im*SIN(phase))/sqrt(ka2(j,1))
                  az(n-j+2,1) = conjg(az(j,1))
                  phase = 2*pi*randu(seed)
                  bz(j,1) = (COS(phase)+im*SIN(phase))/sqrt(ka2(j,1))
                  bz(n-j+2,1) = conjg(bz(j,1))
               ELSE
                  az(j,1) = 0.
                  az(n-j+2,1) = 0.
                  bz(j,1) = 0.
                  bz(n-j+2,1) = 0.
               ENDIF
            END DO
            DO j = 1,n
               DO i = 2,iend
                  IF ((ka2(j,i).le.mkup**2).and.(ka2(j,i).ge.mkdn**2)) THEN
                     phase = 2*pi*randu(seed)
                     az(j,i) = 2*(COS(phase)+im*SIN(phase))/sqrt(ka2(j,i))
                     phase = 2*pi*randu(seed)
                     bz(j,i) = 2*(COS(phase)+im*SIN(phase))/sqrt(ka2(j,i))
                  ELSE
                     az(j,i) = 0.
                     bz(j,i) = 0.
                  ENDIF
               END DO
            END DO
         ELSE
             DO j = 1,n
               DO i = ista,iend
                  IF ((ka2(j,i).le.mkup**2).and.(ka2(j,i).ge.mkdn**2)) THEN
                     phase = 2*pi*randu(seed)
                     az(j,i) = 2*(COS(phase)+im*SIN(phase))/sqrt(ka2(j,i))
                     phase = 2*pi*randu(seed)
                     bz(j,i) = 2*(COS(phase)+im*SIN(phase))/sqrt(ka2(j,i))
                  ELSE
                     az(j,i) = 0.
                     bz(j,i) = 0.
                  ENDIF
               END DO
            END DO
         ENDIF
         CALL energy(az,tmp,1)
         CALL energy(bz,tmq,2)
         CALL MPI_BCAST(tmp,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
         CALL MPI_BCAST(tmq,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
         az = az*sqrt(a0)/sqrt(tmp)
         bz = bz*sqrt(b0)/sqrt(tmq)
         IF ((corr.gt.tiny).and.(p0.gt.tiny)) THEN
            az = corr*sqrt(a0)*ps/sqrt(p0)+sqrt(1.-corr**2)*az
         ENDIF
         IF ((corr.gt.tiny).and.(u0.gt.tiny)) THEN
            bz = corr*sqrt(a0)*vz/sqrt(u0)+sqrt(1.-corr**2)*bz
         ENDIF

      ELSE

         ini = int((stat-1)*tstep)
         dump = float(ini)/float(sstep)+1
         times = 0
         timet = 0
         timec = 0
         timef = fstep
         jc = 48+int(dump/100)
         jd = 48+int(dump/10)-int(dump/100)*10
         ju = 48+int(dump)-int(dump/10)*10
         ic = 48+int(stat/100)
         id = 48+int(stat/10)-int(stat/100)*10
         iu = 48+int(stat)-int(stat/10)*10
         c = char(ic)
         d = char(id)
         u = char(iu)

         OPEN(1,file='hall25Dps.' // node // '.' // c // d // u // &
              '.out',form='unformatted')
         READ(1) R1
         CLOSE(1)
         OPEN(1,file='hall25Daz.' // node // '.' // c // d // u // &
              '.out',form='unformatted')
         READ(1) R2
         CLOSE(1)

         CALL fftp2d_real_to_complex(planrc,R1,ps,MPI_COMM_WORLD)
         CALL fftp2d_real_to_complex(planrc,R2,az,MPI_COMM_WORLD)

         OPEN(1,file='hall25Dvz.' // node // '.' // c // d // u // &
              '.out',form='unformatted')
         READ(1) R1
         CLOSE(1)
         OPEN(1,file='hall25Dbz.' // node // '.' // c // d // u // &
              '.out',form='unformatted')
         READ(1) R2
         CLOSE(1)

         CALL fftp2d_real_to_complex(planrc,R1,vz,MPI_COMM_WORLD)
         CALL fftp2d_real_to_complex(planrc,R2,bz,MPI_COMM_WORLD)

      ENDIF

!
! Time integration scheme starts here
! Uses Runge-Kutta of order 'ord'

      time = (ini-1)*dt
      IF (bench.eq.1) THEN
         CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
         CALL CPU_Time(cputime1)
      ENDIF
 RK : DO t = ini,step

! Sets the external force. Every 'fsteps'
! the phases are changed randomly.

         IF (timef.eq.fstep) THEN
            timef = 0
            IF (ista.eq.1) THEN
               fk(1,1) = 0.
               fz(1,1) = 0.
               mk(1,1) = 0.
               mz(1,1) = 0.
               DO j = 2,n/2+1
                  IF ((ka2(j,1).le.kup**2).and.(ka2(j,1).ge.kdn**2)) THEN
                     phase = 2*pi*randu(seed)
                     fk(j,1) = (COS(phase)+im*SIN(phase))/sqrt(ka2(j,1))
                     fk(n-j+2,1) = conjg(fk(j,1))
                     phase = 2*pi*randu(seed)
                     fz(j,1) = (COS(phase)+im*SIN(phase))/sqrt(ka2(j,1))
                     fz(n-j+2,1) = conjg(fz(j,1))
                  ELSE
                     fk(j,1) = 0.
                     fk(n-j+2,1) = 0.
                     fz(j,1) = 0.
                     fz(n-j+2,1) = 0.
                  ENDIF
                  IF ((ka2(j,1).le.mkup**2).and.(ka2(j,1).ge.mkdn**2)) THEN
                     phase = 2*pi*randu(seed)
                     mk(j,1) = (COS(phase)+im*SIN(phase))/sqrt(ka2(j,1))
                     mk(n-j+2,1) = conjg(mk(j,1))
                     phase = 2*pi*randu(seed)
                     mz(j,1) = (COS(phase)+im*SIN(phase))/sqrt(ka2(j,1))
                     mz(n-j+2,1) = conjg(mz(j,1))
                  ELSE
                     mk(j,1) = 0.
                     mk(n-j+2,1) = 0.
                     mz(j,1) = 0.
                     mz(n-j+2,1) = 0.
                  ENDIF
               END DO
               DO i = 2,iend
                  DO j = 1,n
                     IF ((ka2(j,i).le.kup**2).and.(ka2(j,i).ge.kdn**2)) THEN
                        phase = 2*pi*randu(seed)
                        fk(j,i) = 2*(COS(phase)+im*SIN(phase))/sqrt(ka2(j,i))
                        phase = 2*pi*randu(seed)
                        fz(j,i) = 2*(COS(phase)+im*SIN(phase))/sqrt(ka2(j,i))
                     ELSE
                        fk(j,i) = 0.
                        fz(j,i) = 0.
                     ENDIF
                     IF ((ka2(j,i).le.mkup**2).and.(ka2(j,i).ge.mkdn**2)) THEN
                        phase = 2*pi*randu(seed)
                        mk(j,i) = 2*(COS(phase)+im*SIN(phase))/sqrt(ka2(j,i))
                        phase = 2*pi*randu(seed)
                        mz(j,i) = 2*(COS(phase)+im*SIN(phase))/sqrt(ka2(j,i))
                     ELSE
                        mk(j,i) = 0.
                        mz(j,i) = 0.
                     ENDIF
                  END DO
               END DO
            ELSE
               DO i = ista,iend
                  DO j = 1,n
                     IF ((ka2(j,i).le.kup**2).and.(ka2(j,i).ge.kdn**2)) THEN
                        phase = 2*pi*randu(seed)
                        fk(j,i) = 2*(COS(phase)+im*SIN(phase))/sqrt(ka2(j,i))
                        phase = 2*pi*randu(seed)
                        fz(j,i) = 2*(COS(phase)+im*SIN(phase))/sqrt(ka2(j,i))
                     ELSE
                        fk(j,i) = 0.
                        fz(j,i) = 0.
                     ENDIF
                     IF ((ka2(j,i).le.mkup**2).and.(ka2(j,i).ge.mkdn**2)) THEN
                        phase = 2*pi*randu(seed)
                        mk(j,i) = 2*(COS(phase)+im*SIN(phase))/sqrt(ka2(j,i))
                        phase = 2*pi*randu(seed)
                        mz(j,i) = 2*(COS(phase)+im*SIN(phase))/sqrt(ka2(j,i))
                     ELSE
                        mk(j,i) = 0.
                        mz(j,i) = 0.
                     ENDIF
                  END DO
               END DO
            ENDIF
            CALL energy(fk,tmp,1)
            CALL energy(fz,tmq,1)
            CALL energy(mk,tmr,1)
            CALL energy(mz,tms,1)
            CALL MPI_BCAST(tmp,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            CALL MPI_BCAST(tmq,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            CALL MPI_BCAST(tmr,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            CALL MPI_BCAST(tms,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            tmp = sqrt(f0)/sqrt(tmp)
            tmq = sqrt(z0)/sqrt(tmq)
            tmr = sqrt(m0)/sqrt(tmr)
            tms = sqrt(n0)/sqrt(tms)
            DO i = ista,iend
               DO j = 1,n
                  fk(j,i) = fk(j,i)*tmp
                  fz(j,i) = fz(j,i)*tmq
                  mk(j,i) = mk(j,i)*tmr
                  mz(j,i) = mz(j,i)*tms
               END DO
            END DO
         ENDIF

! Every 'cstep' steps, generates external files 
! to check consistency and convergency. See the 
! mhdcheck25 subroutine for details.

         IF ((timec.eq.cstep).and.(bench.eq.0)) THEN
            timec = 0
            CALL mhdcheck25(ps,az,vz,bz,fk,mk,fz,mz,time,1)
         ENDIF

! Every 'sstep' steps, generates external files 
! with the power spectrum, and checks if adaptive 
! time step adjustment is needed

         IF ((times.eq.sstep).and.(bench.eq.0)) THEN
            times = 0
            ju = ju+1
            IF (ju.eq.58) THEN
               ju = 48
               jd = jd+1
            ENDIF
            IF (jd.eq.58) THEN
               jd = 48
               jc = jc+1
            ENDIF
            c = char(jc)
            d = char(jd)
            u = char(ju)
            ext = c // d // u
            CALL spectrum25(ps,vz,ext,1)
            CALL spectrum25(az,bz,ext,0)
            IF ((adap.gt.0).and.(stat.gt.0)) THEN
               CALL energy25(ps,vz,tmp,1)
               CALL energy25(az,bz,tmq,1)
               IF (myrank.eq.0) THEN
                  dt = 1./(float(adap)*float(mult)*sqrt(max(tmp,tmq))*n)
                  timef = int(timef*cort/(dt*fstep))
                  fstep = int(cort/dt)
                  IF (fstep.eq.0) THEN
                     fstep = 1
                     timef = 0
                  ENDIF
               ENDIF
               CALL MPI_BCAST(dt,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
               CALL MPI_BCAST(timef,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
               CALL MPI_BCAST(fstep,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            ENDIF
            stat = 1
         ENDIF

! Every 'tstep' steps, stores the results of the integration

         IF ((timet.eq.tstep).and.(bench.eq.0)) THEN
            timet = 0
            iu = iu+1
            IF (iu.eq.58) THEN
               iu = 48
               id = id+1
            ENDIF
            IF (id.eq.58) THEN
               id = 48
               ic = ic+1
            ENDIF
            c = char(ic)
            d = char(id)
            u = char(iu)

            DO i = ista,iend
               DO j = 1,n
                  C1(j,i) = ps(j,i)/float(n)**2
                  C2(j,i) = az(j,i)/float(n)**2
                  C3(j,i) = vz(j,i)/float(n)**2
                  C4(j,i) = bz(j,i)/float(n)**2
               END DO
            END DO

            CALL fftp2d_complex_to_real(plancr,C1,R1,MPI_COMM_WORLD)
            CALL fftp2d_complex_to_real(plancr,C2,R2,MPI_COMM_WORLD)

            OPEN(1,file=trim(ldir) // '/hall25Dps.' // node // '.' &
                 // c // d // u // '.out',form='unformatted')
            WRITE(1) R1
            CLOSE(1)
            OPEN(1,file=trim(ldir) // '/hall25Daz.' // node // '.' &
                 // c // d // u // '.out',form='unformatted')
            WRITE(1) R2
            CLOSE(1)

            CALL fftp2d_complex_to_real(plancr,C3,R1,MPI_COMM_WORLD)
            CALL fftp2d_complex_to_real(plancr,C4,R2,MPI_COMM_WORLD)

            OPEN(1,file=trim(ldir) // '/hall25Dvz.' // node // '.' &
                 // c // d // u // '.out',form='unformatted')
            WRITE(1) R1
            CLOSE(1)
            OPEN(1,file=trim(ldir) // '/hall25Dbz.' // node // '.' &
                 // c // d // u // '.out',form='unformatted')
            WRITE(1) R2
            CLOSE(1)
         ENDIF

         timet = timet+1
         times = times+1
         timec = timec+1
         timef = timef+1
         time = time+dt

! Runge-Kutta step 1
! Copies the streamfunction, the vector potential, and the z
! components of the fields into the auxiliary matrixes C1-C4

         DO i = ista,iend
            DO j = 1,n
               C1(j,i) = ps(j,i)
               C2(j,i) = az(j,i)
               C3(j,i) = vz(j,i)
               C4(j,i) = bz(j,i)
            END DO
         END DO

! Runge-Kutta step 2

         DO o = ord,2,-1

         CALL laplak2(C1,C5)
         CALL laplak2(C2,C6)
         CALL poisson(C1,C2,C7)
         CALL poisson(C1,C3,C8)
         CALL poisson(C1,C4,C9)
         CALL poisson(C2,C3,C10)
         CALL poisson(C2,C4,C11)
         CALL poisson(C1,C5,C1)
         CALL poisson(C2,C6,C2)
         CALL laplak2(C3,C3)
         CALL laplak2(C4,C4)

         DO i = ista,iend
            DO j = 1,n

            IF ((ka2(j,i).le.kmax).and.(ka2(j,i).ge.tiny)) THEN
               C1(j,i) = ps(j,i)+dt*(nu*C5(j,i)+(C2(j,i)-C1(j,i))    &
              /ka2(j,i)+fk(j,i))/float(o)
               C2(j,i) = az(j,i)+dt*(mu*C6(j,i)+C7(j,i)+ep*C11(j,i)  &
              +mk(j,i))/float(o)
               C3(j,i) = vz(j,i)+dt*(nu*C3(j,i)+C8(j,i)-C11(j,i)     &
              +fz(j,i))/float(o)
               C4(j,i) = bz(j,i)+dt*(mu*C4(j,i)+C9(j,i)-C10(j,i)     &
              -ep*C2(j,i)+mz(j,i))/float(o)
            ELSE
               C1(j,i) = 0.
               C2(j,i) = 0.
               C3(j,i) = 0.
               C4(j,i) = 0.
            ENDIF

            END DO
         END DO
         
         END DO

! Runge-Kutta step 3
! Copies the result from the auxiliary matrixes into ps, az, vz, bz

         o = 1

         CALL laplak2(C1,C5)
         CALL laplak2(C2,C6)
         CALL poisson(C1,C2,C7)
         CALL poisson(C1,C3,C8)
         CALL poisson(C1,C4,C9)
         CALL poisson(C2,C3,C10)
         CALL poisson(C2,C4,C11)
         CALL poisson(C1,C5,C1)
         CALL poisson(C2,C6,C2)
         CALL laplak2(C3,C3)
         CALL laplak2(C4,C4)

         DO i = ista,iend
            DO j = 1,n

            IF ((ka2(j,i).le.kmax).and.(ka2(j,i).ge.tiny)) THEN
               ps(j,i) = ps(j,i)+dt*(nu*C5(j,i)+(C2(j,i)-C1(j,i))    &
              /ka2(j,i)+fk(j,i))/float(o)
               az(j,i) = az(j,i)+dt*(mu*C6(j,i)+C7(j,i)+ep*C11(j,i)  &
              +mk(j,i))/float(o)
               vz(j,i) = vz(j,i)+dt*(nu*C3(j,i)+C8(j,i)-C11(j,i)     &
              +fz(j,i))/float(o)
               bz(j,i) = bz(j,i)+dt*(mu*C4(j,i)+C9(j,i)-C10(j,i)     &
              -ep*C2(j,i)+mz(j,i))/float(o)
            ELSE
               ps(j,i) = 0.
               az(j,i) = 0.
               vz(j,i) = 0.
               bz(j,i) = 0.
            ENDIF

            END DO
         END DO

         IF (time.gt.ttime) EXIT RK

      END DO RK

! Computes the benchmark
      IF (bench.eq.1) THEN
         CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
         CALL CPU_Time(cputime2)
         IF (myrank.eq.0) THEN
            OPEN(1,file='benchmark.txt',position='append')
            WRITE(1,*) nprocs,(cputime2-cputime1)/(step-ini+1)
            CLOSE(1)
         ENDIF
      ENDIF

!
! End of HALLMHD25

      CALL MPI_FINALIZE(ierr)
      CALL fftp2d_destroy_plan(plancr)
      CALL fftp2d_destroy_plan(planrc)
      DEALLOCATE( R1,R2 )
      DEALLOCATE( ps,az,vz,bz )
      DEALLOCATE( fk,fz,mk,mz )
      DEALLOCATE( C1,C2,C3,C4,C5 )
      DEALLOCATE( C6,C7,C8,C9 )
      DEALLOCATE( C10,C11 )
      DEALLOCATE( ka,ka2 )

      END PROGRAM HALLMHD25D
