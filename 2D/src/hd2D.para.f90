!=================================================================
      PROGRAM HD2D
!=================================================================
! HD2D code
!
! Numerically integrates the incompressible HD equations 
! in 2 dimensions using the streamfunction formulation 
! with an external force. 
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
! 2004 Pablo D. Mininni.
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
      IMPLICIT NONE

!
! Integration parameters
!     ord  : order of the Runge-Kutta method used

      INTEGER, PARAMETER :: ord = 2
      INTEGER :: ini
      INTEGER :: step
      INTEGER :: tstep
      INTEGER :: cstep
      INTEGER :: sstep

!
! streamfunction, vector potential, z component 
! of the fields and external force matrixes

      COMPLEX, ALLOCATABLE, DIMENSION (:,:) :: ps
      COMPLEX, ALLOCATABLE, DIMENSION (:,:) :: fk

!
! Temporal data storing matrixes

      COMPLEX, ALLOCATABLE, DIMENSION (:,:) :: C1
      COMPLEX, ALLOCATABLE, DIMENSION (:,:) :: C5
      REAL, ALLOCATABLE, DIMENSION (:,:)    :: R1

!
! Some auxiliary matrixes

      REAL :: dt,nu
      REAL :: kup,kdn
      REAL :: amus,amuc
      REAL :: stat,dump
      REAL :: f0,u0
      REAL :: time

      INTEGER :: mult
      INTEGER :: t,o
      INTEGER :: i,j
      INTEGER :: ki,kj
      INTEGER :: ic,id,iu
      INTEGER :: jc,jd,ju
      INTEGER :: timet,timec,times

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

      ALLOCATE( R1(n,jsta:jend) )
      ALLOCATE( C1(n,ista:iend) )
      ALLOCATE( C5(n,ista:iend) )
      ALLOCATE( ps(n,ista:iend) )
      ALLOCATE( fk(n,ista:iend) )
      ALLOCATE( ka(n), ka2(n,ista:iend) )

!
! Reads from the external file 'status.txt'
! the status of a previous run (if any)
!     stat: last output of a previous run
!     mult: time step multiplier

      IF (myrank.eq.0) THEN
         OPEN(1,file='status.txt',status='unknown')
         READ(1,*) stat
         READ(1,*) mult
         CLOSE(1)
      ENDIF
      CALL MPI_BCAST(stat,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)

!
! Reads from the external file 'parameter.txt' the 
! parameters that will be used during the integration
!     dt   : time step size
!     step : total number of time steps to compute
!     tstep: number of steps between I/O writing
!     sstep: number of steps between power spectrum I/O
!     cstep: number of steps between information output
!     f0   : amplitude of the external kinetic force
!     u0   : amplitude of the initial streamfunction
!     kdn  : minimum wave number in the external force
!     kup  : maximum wave number in the external force
!     nu   : kinematic viscosity
!     amus : amplitude of the sin terms in the force
!     amuc : amplitude of the cos terms in the force
!     ldir : local directory for I/O

      IF (myrank.eq.0) THEN
         OPEN(1,file='parameter.txt',status='unknown')
         READ(1,*) dt
         READ(1,*) step
         READ(1,*) tstep
         READ(1,*) sstep
         READ(1,*) cstep
         READ(1,*) f0
         READ(1,*) u0
         READ(1,*) kdn
         READ(1,*) kup
         READ(1,*) nu
         READ(1,*) amus
         READ(1,*) amuc
         READ(1,'(a100)') ldir
         CLOSE(1)
         dt = dt/float(mult)
         step = step*mult
         tstep = tstep*mult
         sstep = sstep*mult
         cstep = cstep*mult
      ENDIF
      CALL MPI_BCAST(dt,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(step,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(tstep,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(sstep,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(cstep,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(f0,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(u0,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(kdn,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(kup,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(nu,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(amus,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(amuc,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(ldir,100,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)

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
! velocity field is given by the problem of the merger of two 
! positive vortices and a negative vortex (Schneider et al. 
! Teor. Comp. Fluid Dyn. 9, 191). In other case, a previous 
! run is continued from the step indicated by the number 
! found in this file.

      IF (stat.eq.0) THEN

         ini = 1
         timet = tstep
         timec = cstep
         times = sstep

         DO j = jsta,jend
            DO i = 1,n
               R1(i,j) = amus*(exp(-((2*pi*(float(i)-1)/float(n)-3*pi/4)**2 &
                       +(2*pi*(float(j)-1)/float(n)-pi)**2)/amuc**2)        &
                       +exp(-((2*pi*(float(i)-1)/float(n)-5*pi/4)**2        &
                       +(2*pi*(float(j)-1)/float(n)-pi)**2)/amuc**2)        &
                       -.5*exp(-((2*pi*(float(i)-1)/float(n)-5*pi/4)**2     &
                       +(2*pi*(float(j)-1)/float(n)-pi*(1+1./(2*sqrt(2.)))) &
                       **2)/amuc**2))
            END DO
         END DO
         CALL fftp2d_real_to_complex(planrc,R1,ps,MPI_COMM_WORLD)
         DO i = ista,iend
            DO j = 1,n
               IF ((ka2(j,i).le.kmax).and.(ka2(j,i).ge.tiny)) THEN
                  fk(j,i) = f0*nu*ps(j,i)
                  ps(j,i) = u0*ps(j,i)/ka2(j,i)
               ELSE
                  fk(j,i) = 0.
                  ps(j,i) = 0.
               ENDIF
            END DO
         END DO

      ELSE

         ini = int((stat-1)*tstep)
         dump = float(ini)/float(sstep)+1
         times = 0
         timet = 0
         timec = 0
         jc = 48+int(dump/100)
         jd = 48+int(dump/10)-int(dump/100)*10
         ju = 48+int(dump)-int(dump/10)*10
         ic = 48+int(stat/100)
         id = 48+int(stat/10)-int(stat/100)*10
         iu = 48+int(stat)-int(stat/10)*10
         c = char(ic)
         d = char(id)
         u = char(iu)

         OPEN(1,file='hd2Dps.' // node // '.' // c // d // u // &
              '.out',form='unformatted')
         READ(1) R1
         CLOSE(1)

         CALL fftp2d_real_to_complex(planrc,R1,ps,MPI_COMM_WORLD)

      ENDIF

!
! Time integration scheme starts here
! Uses Runge-Kutta of order 'ord'

      time = (ini-1)*dt
 RK : DO t = ini,step

! Every 'cstep' steps, generates external files 
! to check consistency and convergency. See the 
! hdcheck subroutine for details.

         IF (timec.eq.cstep) THEN
            timec = 0
            CALL hdcheck(ps,fk,time)
         ENDIF

! Every 'sstep' steps, generates external files 
! with the power spectrum

         IF (times.eq.sstep) THEN
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
            CALL spectrum(ps,ext,1)
         ENDIF

! Every 'tstep' steps, stores the results of the integration

         IF (timet.eq.tstep) THEN
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
               END DO
            END DO

            CALL fftp2d_complex_to_real(plancr,C1,R1,MPI_COMM_WORLD)

            OPEN(1,file=trim(ldir) // '/hd2Dps.' // node // '.' &
                 // c // d // u // '.out',form='unformatted')
            WRITE(1) R1
            CLOSE(1)
         ENDIF

         timet = timet+1
         times = times+1
         timec = timec+1
         time = time+dt

! Runge-Kutta step 1
! Copies the streamfunction into the auxiliary matrix C1

         DO i = ista,iend
            DO j = 1,n
               C1(j,i) = ps(j,i)
            END DO
         END DO

! Runge-Kutta step 2

         DO o = ord,2,-1

         CALL laplak2(C1,C5)
         CALL poisson(C1,C5,C1)

         DO i = ista,iend
            DO j = 1,n

            IF ((ka2(j,i).le.kmax).and.(ka2(j,i).ge.tiny)) THEN
               C1(j,i) = ps(j,i)+dt*(nu*C5(j,i)-C1(j,i)/ka2(j,i)   &
              +fk(j,i))/float(o)
            ELSE
               C1(j,i) = 0.
            ENDIF

            END DO
         END DO
         
         END DO

! Runge-Kutta step 3
! Copies the result from the auxiliary matrixes into ps, az

         o = 1

         CALL laplak2(C1,C5)
         CALL poisson(C1,C5,C1)

         DO i = ista,iend
            DO j = 1,n

            IF ((ka2(j,i).le.kmax).and.(ka2(j,i).ge.tiny)) THEN
               ps(j,i) = ps(j,i)+dt*(nu*C5(j,i)-C1(j,i)/ka2(j,i)   &
              +fk(j,i))/float(o)
            ELSE
               ps(j,i) = 0.
            ENDIF

            END DO
         END DO

      END DO RK

!
! End of Runge-Kutta

      CALL MPI_FINALIZE(ierr)
      CALL fftp2d_destroy_plan(plancr)
      CALL fftp2d_destroy_plan(planrc)
      DEALLOCATE( R1 )
      DEALLOCATE( ps,fk )
      DEALLOCATE( C1,C5 )
      DEALLOCATE( ka,ka2 )

      END PROGRAM HD2D
