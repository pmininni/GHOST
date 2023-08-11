
#ifdef PROTH_SOL
#define SCALAR_
#endif

#ifdef BOUSS_SOL
#define BOUSSINESQ_
#define SCALAR_
#endif

#ifdef ROTBOUSS_SOL
#define BOUSSINESQ_
#define SCALAR_
#define ROTATION_
#endif

#if !defined(PROTH_SOL) && !defined(BOUSS_SOL) && !defined(ROTBOUSS_SOL)
#error "Invalid solver specified"
#endif

#if defined(SCALAR_)
  #error "SCALAR defined"
#endif


      PROGRAM SHEAR3D
!=================================================================
! SHEAR3D code (part of the GHOST suite)
!
! Reads velocity binaries and computes all required
! components of strain rate tensor in Fourier space, computes
! eigenvalues, and isotropic spectra of the eigenvallue fields.
! Note: this utility is _strictly_ for incompressible flows!
!
! Additionally, if isolve > 0 , the utility will also compute the
! eigenvector at each spatial location corresponding to the max
! |eigenvalue|. This gives what might be termed the 'principle-principle'
! axis representing the direction of the largest amount of 'shear'
! This computation relied heavily on the following article:
!
!  'Eigensystems for 3x3 Symmetric Matrices (Revisited)'
!   David Eberly, Geometric Tools, LLC
!   http://www.geometrictools.com
!   (c) 1998-2012, May 2011 
!
! If jpdf=1, jpoint and 1d pdfs of energy dissipation and 
! entrophy density, energy diss and lambda, energy diss and (relative?)
! helicity will be computed and written, as will a variety of other 
! 1d and joint PDFs.
!
!
! 2013 D. Rosenberg
!      ORNL/NCCS
!
! 28 May 2013: Initial version
!=================================================================

!
! Definitions for conditional compilation

! Modules

      USE fprecision
      USE commtypes
      USE mpivars
      USE filefmt
      USE iovar
      USE iompi
      USE grid
      USE fft
      USE ali
      USE var
      USE kes
      USE fftplans
      USE random
      USE threads
      USE gutils
      USE gtimer
      IMPLICIT NONE

!
! Arrays for the fields and structure functions

      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: vx,vy,vz,th
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: ctmp,sij
      COMPLEX(KIND=GP)                                 :: cdump,jdump


      REAL(KIND=GP)    :: dt,nu,mu
      REAL(KIND=GP)    :: kup,kdn
      REAL(KIND=GP)    :: rmp,rmq,rms
      REAL(KIND=GP)    :: rmt,rm1,rm2
      REAL(KIND=GP)    :: dump
      REAL(KIND=GP)    :: stat
      REAL(KIND=GP)    :: f0,u0
      REAL(KIND=GP)    :: phase,ampl,cort
      REAL(KIND=GP)    :: fparam0,fparam1,fparam2,fparam3,fparam4
      REAL(KIND=GP)    :: fparam5,fparam6,fparam7,fparam8,fparam9
      REAL(KIND=GP)    :: vparam0,vparam1,vparam2,vparam3,vparam4
      REAL(KIND=GP)    :: vparam5,vparam6,vparam7,vparam8,vparam9

      
      REAL   (KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: lamb,R1,R2,R3,R4,R5,R6
      REAL   (KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: evx,evy,evz
      REAL   (KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: beta,pv,rig,rign
      REAL   (KIND=GP)                                 :: btrunc,sg,sl,tmp
      REAL   (KIND=GP)                                 :: ensmin,ensmax,dismin,dismax
      REAL   (KIND=GP)                                 :: krmin,krmax
      REAL   (KIND=GP)                                 :: krmin2,krmax2 
      REAL   (KIND=GP)                                 :: ktmin ,ktmax 
      REAL   (KIND=GP)                                 :: ktmin2,ktmax2
      REAL   (KIND=GP)                                 :: fmin(2),fmax(2)
#ifdef BOUSSINESQ_
      REAL   (KIND=GP)                                 :: bvfreq,xmom,xtemp
#endif
#ifdef ROTATION_
      REAL   (KIND=GP)                                 :: omegax,omegay,omegaz
#endif
      REAL   (KIND=GP)                                 :: sheb(3),sup,suz,sth,fup,fuz,fth
#ifdef SCALAR_
      REAL(KIND=GP)    :: kappa
      REAL(KIND=GP)    :: c0,s0
      REAL(KIND=GP)    :: cparam0,cparam1,cparam2,cparam3,cparam4
      REAL(KIND=GP)    :: cparam5,cparam6,cparam7,cparam8,cparam9
      REAL(KIND=GP)    :: sparam0,sparam1,sparam2,sparam3,sparam4
      REAL(KIND=GP)    :: sparam5,sparam6,sparam7,sparam8,sparam9
#endif
#if defined(SCALAR_) || defined(MULTISCALAR_)
      REAL(KIND=GP)    :: skup,skdn
#endif


!
! Auxiliary variables
      INTEGER :: ini,step
      INTEGER :: tstep,cstep
      INTEGER :: sstep,fstep
      INTEGER :: bench,trans
      INTEGER :: outs,mean
      INTEGER :: seed,rand
      INTEGER :: mult

      INTEGER :: i,ic,iir,ind,ir,isolve,istrain,it,j,jc,jjc,k
      INTEGER :: demean,dolog,ilamb,iobin,inorm,istat(4096),jpdf,nstat
      INTEGER :: irand,isheb,nbinx,nbiny,nbins(2),nin
      INTEGER :: indt
!$    INTEGER, EXTERNAL :: omp_get_max_threads


      TYPE(IOPLAN) :: planio
      CHARACTER(len=16)   :: suff
      CHARACTER(len=128)  :: pref
      CHARACTER(len=256)  :: odir,idir
      CHARACTER(len=1024) :: fnout
      CHARACTER(len=2048) :: fntmp
      CHARACTER(len=2048) :: fntmp1,fntmp2,fntmp3,fntmp4,fntmp5,fntmp6
      CHARACTER(len=4096) :: sstat
!
      NAMELIST / status / idir,odir,stat,mult,bench,outs,mean,trans,iswap

      NAMELIST / parameter / dt,step,tstep,sstep,cstep,rand,cort,seed
!
      NAMELIST / velocity / f0,u0,kdn,kup,nu,fparam0,fparam1,fparam2
      NAMELIST / velocity / fparam3,fparam4,fparam5,fparam6,fparam7
      NAMELIST / velocity / fparam8,fparam9,vparam0,vparam1,vparam2
      NAMELIST / velocity / vparam3,vparam4,vparam5,vparam6,vparam7
      NAMELIST / velocity / vparam8,vparam9
#ifdef SCALAR_
      NAMELIST / scalar / c0,s0,skdn,skup,kappa,cparam0,cparam1
      NAMELIST / scalar / cparam2,cparam3,cparam4,cparam5,cparam6
      NAMELIST / scalar / cparam7,cparam8,cparam9,sparam0,sparam1
      NAMELIST / scalar / sparam2,sparam3,sparam4,sparam5,sparam6
      NAMELIST / scalar / sparam7,sparam8,sparam9

#endif
#ifdef ROTATION_
      NAMELIST / rotation / omegax,omegay,omegaz
#endif
#ifdef BOUSSINESQ_
      NAMELIST / boussinesq / bvfreq,xmom,xtemp
#endif
      NAMELIST / shear / demean,ilamb,iobin,isheb,isolve,iswap
      NAMELIST / shear / dolog,oswap,idir,odir,pref,sstat
      NAMELIST / shear / dismin,dismax,ensmin,ensmax,jpdf,nbinx,nbiny
      NAMELIST / shear / irand,krmin,krmax
      NAMELIST / shear / btrunc,ktmin,ktmax

!
! Initializes the MPI and I/O libraries
      CALL MPI_INIT(ierr)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
      CALL range(1,nx/2+1,nprocs,myrank,ista,iend)
      CALL range(1,nz,nprocs,myrank,ksta,kend)
      CALL io_init(myrank,(/nx,ny,nz/),ksta,kend,planio)
      nth = 1
!$    nth = omp_get_max_threads()
!$    CALL fftp3d_init_threads(ierr)
!     tiny: minimum truncation for dealiasing
      tiny = 1e-5
     
      kmax   = real(nx/2+1,kind=GP)

      idir   = '.'
      odir   = '.'
      sstat  = '0'
      iswap  = 0
      oswap  = 0
      btrunc = 0
      demean = 0
      dolog  = 1
      iobin = 0
      ilamb  = 0
      irand  = 0
      isheb  = 0
      isolve = 0
      jpdf   = 3
      krmin  = tiny
      krmax  = kmax
      ktmin  = tiny
      ktmax  = kmax
      seed   = 1000
      pref   = 'ksplambda'

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
!     outs : = 0 writes velocity [and vector potential (MAGFIELD_)]
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
!
! Reads from the external file 'shear`.txt' the 
! parameters that will be used to compute the transfer
!     idir   : directory for unformatted input (field components)
!     odir   : directory for unformatted output (prolongated data)
!     sstat  : time index for which to compute SHEAR, or a ';--separated list
!     btrunc : if == 1, truncate spectral range to [ktmin,ktmax]
!     ktmin  : min wavenumber for truncation if btrunc=1
!     ktmax  : max wavenumber for truncation if btrunc=1
!     iswap  : do endian swap on input?
!     oswap  : do endian swap on output?
!     isolve : 0==>just max |eigenvale| field; 1==>max |eigenvalue| field plus
!              corresponding eigenvector inevx,evy,evz
!     ilamb  : 0: don't do evalue computation; 1: do it (is done for jpdf>0).
!     iobin : 1==>write eigenvalue field to disk; 0==>don't
!     irand  : randomize phases between [krmin,krmax] if 1; else, don't
!     isheb  : compute Shebaliln angles for total energy, and ke and pe?
!              ==0: don't compute them;
!              ==1: Compute angles
!              ==2: Compute angles and return
!     krmin  : min wavenumber for randomization if irand=1
!     krmax  : max wavenumber for randomization if irand=1
!     jpdf   : 0: do no pdfs
!              1: do 1d pdfs only
!              2: do joint pdfs only
!              3: do both 1d and joint pdfs
!     demean : demean the eigenvalue field?
!     dolog  : compute PDFs in log=space?

      IF (myrank.eq.0) THEN
write(*,*)'main: opening shear.inp...'
         OPEN(1,file='shear.inp',status='unknown',form="formatted")
         READ(1,NML=shear)
         CLOSE(1)
write(*,*)'main: shear.inp read.'
      ENDIF
      CALL MPI_BCAST(idir  ,256 ,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(odir  ,256 ,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(pref  ,128 ,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(sstat ,4096,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(btrunc,1   ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(demean,1   ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(dolog ,1   ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(isheb ,1   ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(iswap ,1   ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(isolve,1   ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(ilamb ,1   ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(iobin,1   ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(irand ,1   ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(jpdf  ,1   ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(ktmin ,1   ,GC_REAL      ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(ktmax ,1   ,GC_REAL      ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(nbinx ,1   ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(nbiny ,1   ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(oswap ,1   ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(krmin ,1   ,GC_REAL      ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(krmax ,1   ,GC_REAL      ,0,MPI_COMM_WORLD,ierr)
if (myrank.eq.0) write(*,*)'main: broadcast done.'
!
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
!            = 2 constant energy
!            = 3 slowly varying random phases (only for the velocity and
!                magnetic forcings)
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

!
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

! Reads parameters specifically for Boussinesq solver from the 
! namelist 'boussinesq' on the external file 'parameter.inp'
!     bvfreq: Brunt-Vaisala frequency (positive definite)
!     xmom  : multiplies bouyancy term in momentum equation
!     xtemp : multiplies temperature-current term in 
!             temperature/density equation
      xmom  = 1.0
      xtemp = 1.0
      IF (myrank.eq.0) THEN
         OPEN(1,file='parameter.inp',status='unknown',form="formatted")
         READ(1,NML=boussinesq)
         CLOSE(1)
      ENDIF
      CALL MPI_BCAST(bvfreq,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(xmom  ,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(xtemp ,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      xmom  = xmom  * bvfreq
      xtemp = xtemp * bvfreq

#ifdef ROTATION_
!
! Reads parameters for runs with rotation from the 
! namelist 'rotation' on the external file 'parameter.inp'
!     omegax, omegay, omegaz: amplitude of the uniform rotation
 
      omegax = 0.0_GP
      omegay = 0.0_GP
      omegaz = 0.0_GP
      IF (myrank.eq.0) THEN
         OPEN(1,file='parameter.inp',status='unknown',form="formatted")
         READ(1,NML=rotation)
         CLOSE(1)
      ENDIF
      CALL MPI_BCAST(omegax,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(omegay,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(omegaz,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
#endif

#ifdef SCALAR_
! Reads parameters for the passive/active scalar from the 
! namelist 'scalar' on the external file 'parameter.ip'
!     s0   : amplitude of the passive scalar source
!     c0   : amplitude of the initial concentration
!     skdn : minimum wave number in concentration/source
!     skup : maximum wave number in concentration/source
!     kappa: diffusivity
!     sparam0-9 : ten real numbers to control properties of 
!            the source
!     cparam0-9 : ten real numbers to control properties of
!            the initial concentration

      IF (myrank.eq.0) THEN
         OPEN(1,file='parameter.inp',status='unknown',form="formatted")
         READ(1,NML=scalar)
         CLOSE(1)
      ENDIF
      CALL MPI_BCAST(s0,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(c0,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(skdn,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(skup,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(kappa,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(sparam0,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(sparam1,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(sparam2,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(sparam3,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(sparam4,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(sparam5,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(sparam6,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(sparam7,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(sparam8,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(sparam9,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(cparam0,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(cparam1,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(cparam2,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(cparam3,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(cparam4,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(cparam5,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(cparam6,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(cparam7,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(cparam8,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(cparam9,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
#endif


      nbins(1) = nbinx
      nbins(2) = nbiny

      istrain = 0
      IF ( jpdf.GT.0    ) ilamb  = 1
      IF ( isolve .GT.0 ) ilamb  = 1
      IF ( jpdf.GT.0 .OR. isolve.GT.0 .OR. ilamb.GT.0 ) istrain = 1
      IF ( ilamb  .LE.0 ) iobin = 0

!
      ALLOCATE( ctmp(nz,ny,ista:iend) )
      IF ( istrain.GT.0 ) THEN
      ALLOCATE( sij (nz,ny,ista:iend) )
      ENDIF
      ALLOCATE( vx(nz,ny,ista:iend) )
      ALLOCATE( vy(nz,ny,ista:iend) )
      ALLOCATE( vz(nz,ny,ista:iend) )
#if defined(SCALAR_)
      ALLOCATE( th(nz,ny,ista:iend) )
#endif
      ALLOCATE( kx(nx), ky(ny), kz(nz) )
      ALLOCATE( kn2(nz,ny,ista:iend) )
#ifdef DEF_ARBSIZE_
!     anis = 1
!     ALLOCATE( kk2(nz,ny,ista:iend) )
#else
      IF ((nx.ne.ny).or.(ny.ne.nz)) THEN
!        anis = 1
         ALLOCATE( kk2(nz,ny,ista:iend) )
      ELSE
!        anis = 0
         kk2 => kn2
      ENDIF
#endif

      ALLOCATE( R1(nx,ny,ksta:kend) )
      ALLOCATE( R2(nx,ny,ksta:kend) )
      ALLOCATE( R3(nx,ny,ksta:kend) )
      ALLOCATE( R4(nx,ny,ksta:kend) )
      ALLOCATE( R5(nx,ny,ksta:kend) )
      ALLOCATE( R6(nx,ny,ksta:kend) )
      IF ( ilamb.GT.0 ) THEN
      ALLOCATE( lamb(nz,ny,ksta:kend) )
      ENDIF
      
      ALLOCATE( evx(nx,ny,ksta:kend) )
      ALLOCATE( evy(nx,ny,ksta:kend) )
      ALLOCATE( evz(nx,ny,ksta:kend) )
      ALLOCATE( beta(nx,ny,ksta:kend) )
      ALLOCATE( pv  (nx,ny,ksta:kend) )
      ALLOCATE( rig (nx,ny,ksta:kend) )
      ALLOCATE( rign(nx,ny,ksta:kend) )

if (myrank.eq.0) write(*,*)'main: creating plans...'
      CALL fftp3d_create_plan(planrc,(/nx,ny,nz/),FFTW_REAL_TO_COMPLEX, &
                             FFTW_ESTIMATE)
      CALL fftp3d_create_plan(plancr,(/nx,ny,nz/),FFTW_COMPLEX_TO_REAL, &
                             FFTW_ESTIMATE)

if (myrank.eq.0) write(*,*)'main: plans done.'
!
! Some constants for the FFT
!     kmax: maximum truncation for dealiasing
      IF ( irand.eq.0 ) THEN
         krmin  = tiny
         krmax  = kmax
      ENDIF
      IF ( btrunc.eq.0 ) THEN
         ktmin  = tiny
         ktmax  = kmax
      ENDIF
      krmin2 = krmin**2
      krmax2 = krmax**2
      ktmin2 = ktmin**2
      ktmax2 = ktmax**2
if (myrank.eq.0) write(*,*)'main: d-extrema done.'

!
! Builds the wave number and the square wave 
! number matrixes

      DO i = 1,nx/2
         kx(i) = real(i-1,kind=GP)
         kx(i+nx/2) = real(i-nx/2-1,kind=GP)
      END DO
      DO j = 1,ny/2
         ky(j) = real(j-1,kind=GP)
         ky(j+ny/2) = real(j-ny/2-1,kind=GP)
      END DO
      DO k = 1,nz/2
         kz(k) = real(k-1,kind=GP)
         kz(k+nz/2) = real(k-nz/2-1,kind=GP)
      END DO

if (myrank.eq.0) write(*,*)'main: calling parseind...'
      CALL parseind(sstat,';', istat , 4096, nstat) 
if (myrank.eq.0) write(*,*)'main: index parsing done: nstat=',nstat

      tmp = 1.0_GP/REAL(nx*ny*nz,KIND=GP)
      DO it = 1,nstat
if (myrank.eq.0) write(*,*)'main: writing ext: fmtext:', fmtext,' istat=',istat(it)
        WRITE(ext, fmtext) istat(it)
! read in appropriate file:
if (myrank.eq.0) write(*,*)'main: reading time index: ', ext, ' ...'
        CALL io_read(1,idir,'vx',ext,planio,R1)
        CALL io_read(1,idir,'vy',ext,planio,R2)
        CALL io_read(1,idir,'vz',ext,planio,R3)
#if defined(SCALAR_)
        CALL io_read(1,idir,'th',ext,planio,R4)
#endif
if (myrank.eq.0) write(*,*)'main: reading done.'

!
! Byte-swap on input:
        IF ( iswap .NE. 0 ) THEN
          CALL rarray_byte_swap(R1, nx*ny*(kend-ksta+1))
          CALL rarray_byte_swap(R2, nx*ny*(kend-ksta+1))
          CALL rarray_byte_swap(R3, nx*ny*(kend-ksta+1))
        ENDIF
!
! take FFT of component:
if (myrank.eq.0) write(*,*)'main: real 2 cmplex for vx...'
        CALL fftp3d_real_to_complex(planrc,R1,vx,MPI_COMM_WORLD)
if (myrank.eq.0) write(*,*)'main: real 2 cmplex for vy...'
        CALL fftp3d_real_to_complex(planrc,R2,vy,MPI_COMM_WORLD)
if (myrank.eq.0) write(*,*)'main: real 2 cmplex for vz...'
        CALL fftp3d_real_to_complex(planrc,R3,vz,MPI_COMM_WORLD)
#if defined(SCALAR_)
if (myrank.eq.0) write(*,*)'main: real 2 cmplex for th...'
        CALL fftp3d_real_to_complex(planrc,R4,th,MPI_COMM_WORLD)
#endif
if (myrank.eq.0) write(*,*)'main: real 2 cmplex done.'
        IF ( irand.GT.0 ) THEN
          IF (myrank.eq.0) phase = 2*pi*randu(seed)
          CALL MPI_BCAST(phase,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
          cdump = COS(phase)+im*SIN(phase)
          jdump = conjg(cdump)
          CALL Randomize(vx,vy,vz,krmin,krmax,phase)
        ENDIF

!$omp parallel do if (iend-ista.ge.nth) private (j,k)
        DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
          DO j = 1,ny
            DO k = 1,nz
              IF ((kk2(k,j,i).lt.ktmin2).or.(kk2(k,j,i).gt.ktmax2)) THEN
                vx(k,j,i) = 0.0_GP
                vy(k,j,i) = 0.0_GP
                vz(k,j,i) = 0.0_GP
#if defined(SCALAR_)
                th(k,j,i) = 0.0_GP
#endif
              ENDIF
            END DO
          END DO
        END DO

! Do new kspecpara:
       CALL specpara(vx,vy,vz,ext,1,1);
!!     cycle

! Do Shebalin angles:
#if defined(SCALAR_)
        IF ( isheb.GT.0 ) THEN
          CALL  Shebalin(sheb(1),0,vx,vy,vz,th)
          CALL  Shebalin(sheb(2),1,vx,vy,vz,th)
          CALL  Shebalin(sheb(3),2,vx,vy,vz,th)
          ! Print out Shebalin angles:
          IF ( myrank.EQ.0 ) THEN
            fnout = trim(odir) // '/' // 'shebalin.txt'
            OPEN(1,file=trim(fnout),position='append')
            WRITE(1,*)ext,sheb(1:3)
            CLOSE(1)
          ENDIF
        ENDIF
#else
        IF ( isheb.GT.0 ) THEN
          CALL  Shebalin(sheb(1),0,vx,vy,vz)
          ! Print out Shebalin angles:
          IF ( myrank.EQ.0 ) THEN
            fnout = trim(odir) // '/' // 'shebalin.txt'
            OPEN(1,file=trim(fnout),position='append')
            WRITE(1,*)ext,sheb(1)
            CLOSE(1)
          ENDIF
        ENDIF
#endif
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
      IF ( isheb.EQ.2 .AND. it.EQ.nstat ) THEN
        IF ( ALLOCATED(ctmp) ) DEALLOCATE(ctmp); IF ( ALLOCATED (sij) ) DEALLOCATE (sij)
        IF ( ALLOCATED  (vx) ) DEALLOCATE  (vx); IF ( ALLOCATED  (vy) ) DEALLOCATE  (vy)
        IF ( ALLOCATED  (vz) ) DEALLOCATE  (vz); IF ( ALLOCATED  (th) ) DEALLOCATE  (th)
        IF ( ALLOCATED(lamb) ) DEALLOCATE(lamb); IF ( ALLOCATED  (R1) ) DEALLOCATE  (R1)
        IF ( ALLOCATED  (R2) ) DEALLOCATE  (R2); IF ( ALLOCATED  (R3) ) DEALLOCATE  (R3)
        IF ( ALLOCATED  (R4) ) DEALLOCATE  (R4); IF ( ALLOCATED  (R5) ) DEALLOCATE  (R5)
        IF ( ALLOCATED  (R6) ) DEALLOCATE  (R6); IF ( ALLOCATED (evx) ) DEALLOCATE (evx)
        IF ( ALLOCATED (evy) ) DEALLOCATE (evy); IF ( ALLOCATED (evz) ) DEALLOCATE (evz)
        IF ( ALLOCATED(beta) ) DEALLOCATE(beta); 
        IF ( ALLOCATED  (pv) ) DEALLOCATE(pv); 
        IF ( ALLOCATED (rig) ) DEALLOCATE(rig); 
        IF ( ALLOCATED(rign) ) DEALLOCATE(rign); 
!       IF ( ALLOCATED (kk2) ) DEALLOCATE (kk2)
        STOP
      ENDIF

#if defined(SCALAR_)
!       indt = (istat(it)-1)*tstep
!       CALL tbouss(vx,vy,vz,th,indt,dt,omega,bvfreq)
!       CALL havgwrite(0,'shear'  ,ext,vx,vy,vz,th,omega,bvfreq) ! hor. avg shear
!       CALL havgwrite(1,'tgradz' ,ext,vx,vy,vz,th,omega,bvfreq) ! hor. avg dtheta/dz
!       CALL havgwrite(2,'hawdtdz',ext,vx,vy,vz,th,omega,bvfreq) ! hor. avg u_z * dtheta/dz
!       CALL havgwrite(3,'hahke'  ,ext,vx,vy,vz,th,omega,bvfreq) ! hor. avg hor. k.e.
!       CALL havgwrite(4,'havke'  ,ext,vx,vy,vz,th,omega,bvfreq) ! hor. avg vert. k.e.
!!      CALL havgwrite(5,'haphel' ,ext,vx,vy,vz,th,omega,bvfreq) ! hor. avg perp. helicity
!       CALL havgwrite(6,'haomzt' ,ext,vx,vy,vz,th,omega,bvfreq) ! hor. avg ometa_z * theta
!       CALL havgwrite(7,'hapv2'  ,ext,vx,vy,vz,th,omega,bvfreq) ! hor. avg pot'l vorticity^2
!       CALL havgwrite(8,'hasuph' ,ext,vx,vy,vz,th,omega,bvfreq) ! hor. avg super-helicity
!       CALL havgwrite(9,'hari'   ,ext,vx,vy,vz,th,omega,bvfreq) ! hor. avg Richardson no.
#endif

!write(*,*)'main: ktmin2=',ktmin2, ' ktmax2=',ktmax2

! Compute required strain rate components:
        inorm = 1

        IF ( istrain.GT.0 ) THEN
if (myrank.eq.0) write(*,*)'main: doing Strain components...'
        CALL Strain(vx,vy,vz,1,1,ktmin,ktmax,inorm,sij,ctmp)
        CALL fftp3d_complex_to_real(plancr,sij,R1,MPI_COMM_WORLD)
        CALL Strain(vx,vy,vz,1,2,ktmin,ktmax,inorm,sij,ctmp)
        CALL fftp3d_complex_to_real(plancr,sij,R2,MPI_COMM_WORLD)
        CALL Strain(vx,vy,vz,1,3,ktmin,ktmax,inorm,sij,ctmp)
        CALL fftp3d_complex_to_real(plancr,sij,R3,MPI_COMM_WORLD)
        CALL Strain(vx,vy,vz,2,2,ktmin,ktmax,inorm,sij,ctmp)
        CALL fftp3d_complex_to_real(plancr,sij,R4,MPI_COMM_WORLD)
        CALL Strain(vx,vy,vz,2,3,ktmin,ktmax,inorm,sij,ctmp)
        CALL fftp3d_complex_to_real(plancr,sij,R5,MPI_COMM_WORLD)
        ENDIF

! Compute required eigenvalue field of strain:
        IF ( ilamb.GT. 0 ) THEN
if (myrank.eq.0) write(*,*)'main: doing Eigen system...'
          IF ( isolve.GT.0 ) THEN
            CALL EigenSolveMax(R1,R2,R3,R4,R5,lamb,evx,evy,evz)
          ELSE 
            CALL EigenValMax(R1,R2,R3,R4,R5,lamb)
          ENDIF
!
! Remove mean, if desired:

        IF ( demean.GT.0 ) THEN
        sl = 0.0_GP
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
        DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
          DO j = 1,ny
            DO i = 1,nx
!$omp atomic
              sl = sl + lamb(i,j,k)
             ENDDO
           ENDDO
         ENDDO
         call MPI_ALLREDUCE(sl,sg,1,GC_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
         sg = sg / real(nx*ny*nz,kind=GP)

!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
        DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
          DO j = 1,ny
            DO i = 1,nx
!$omp atomic
              lamb(i,j,k) = lamb(i,j,k) - sg
           
             ENDDO
           ENDDO
         ENDDO

         ENDIF ! demean

        R6  = lamb
        CALL fftp3d_real_to_complex(planrc,R6,ctmp,MPI_COMM_WORLD)
!
! Compute power spectrum of e-value  and output it:
        IF ( btrunc.EQ.0 ) THEN
          fnout = trim(pref) // '.' // ext  // '.txt';
        ELSE
          WRITE(suff,'(a2,i5.5,a1,i5.5)') '_T', int(ktmin),'_',int(ktmax)
          fnout = trim(pref) // '.' // ext //  trim(suff) // '.txt'
        ENDIF
        fntmp = trim(odir) // '/' // trim(fnout)

        CALL pspectrum(ctmp,fntmp,int(kmax))

        ENDIF ! do lamb computation 
! 
! Prepare eignesystem for output if necessary
! (don't forget to undo swaps for later):
        IF ( iobin.GT.0 ) THEN
if (myrank.eq.0) write(*,*)'main: output lambda...'
          IF ( oswap .NE. 0 ) THEN
            CALL rarray_byte_swap(lamb, nx*ny*(kend-ksta+1))
          ENDIF
          CALL io_write(1,odir,'lmb',ext,planio,lamb)
        ENDIF
        IF ( isolve .GT. 0 ) THEN
          IF ( oswap .NE. 0 ) THEN
            CALL rarray_byte_swap(evx, nx*ny*(kend-ksta+1))
            CALL rarray_byte_swap(evy, nx*ny*(kend-ksta+1))
            CALL rarray_byte_swap(evz, nx*ny*(kend-ksta+1))
          ENDIF
          CALL io_write(1,odir,'evx',ext,planio,evx)
          CALL io_write(1,odir,'evy',ext,planio,evy)
          CALL io_write(1,odir,'evz',ext,planio,evz)
        ENDIF
!  
!  
! Compute and print joint and 1d pdfs for several quantities:
        IF ( jpdf.GT.0 ) THEN
          ! Note, the strain rate components RI & lambda will be modified on exit:
          IF ( oswap .NE. 0 ) THEN
            CALL rarray_byte_swap(lamb, nx*ny*(kend-ksta+1))
          ENDIF

if (myrank.eq.0) write(*,*)'main: call DoKPDF ...'
          CALL DoKPDF(R1,R2,R3,R4,R5,vx,vy,vz,th,lamb,bvfreq,(/omegax,omegay,omegaz/),nu,kappa,ext    , &
               odir,nbins,dolog,ctmp,sij,evx,evy,evz,beta,pv,rig,rign,jpdf,planio,iobin)
        ENDIF

        IF ( myrank.EQ. 0 ) THEN
          write(*,*)'main: fntmp=',trim(fntmp),' ktmin=',ktmin,' ktmax=',ktmax
          write(*,*)'main: time index ', trim(ext), ' done.'
        ENDIF

      ENDDO   ! time (stat) loop
!
      CALL fftp3d_destroy_plan(plancr)
      CALL fftp3d_destroy_plan(planrc)
      CALL MPI_FINALIZE(ierr)

      IF ( ALLOCATED(ctmp) ) DEALLOCATE(ctmp)
      IF ( ALLOCATED (sij) ) DEALLOCATE (sij)
      IF ( ALLOCATED  (vx) ) DEALLOCATE  (vx)
      IF ( ALLOCATED  (vy) ) DEALLOCATE  (vy)
      IF ( ALLOCATED  (vz) ) DEALLOCATE  (vz)
      IF ( ALLOCATED  (th) ) DEALLOCATE  (th)
      IF ( ALLOCATED(lamb) ) DEALLOCATE(lamb)
      IF ( ALLOCATED  (R1) ) DEALLOCATE  (R1)
      IF ( ALLOCATED  (R2) ) DEALLOCATE  (R2)
      IF ( ALLOCATED  (R3) ) DEALLOCATE  (R3)
      IF ( ALLOCATED  (R4) ) DEALLOCATE  (R4)
      IF ( ALLOCATED  (R5) ) DEALLOCATE  (R5)
      IF ( ALLOCATED  (R6) ) DEALLOCATE  (R6)
      IF ( ALLOCATED (evx) ) DEALLOCATE (evx)
      IF ( ALLOCATED (evy) ) DEALLOCATE (evy)
      IF ( ALLOCATED (evz) ) DEALLOCATE (evz)
!!    IF ( ALLOCATED (kk2) ) DEALLOCATE (kk2)

      END PROGRAM SHEAR3D
!-----------------------------------------------------------------
!-----------------------------------------------------------------


      SUBROUTINE EigenValMax(S11,S12,S13,S22,S23,lamb)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!
! Computes the max |eigenvalue| field for strain rate tensor, whos
! required components are specified
!
! Parameters
!     SIJ  : strain rate components (real)
!     lamb : returned eigenvalue field
!
      USE fprecision
      USE commtypes
      USE kes
      USE grid
      USE mpivars
!$    USE threads
      IMPLICIT NONE

      REAL   (KIND=GP), INTENT(INOUT), DIMENSION(nx,ny,ksta:kend) :: S11,S12,S13,S22,S23
      REAL   (KIND=GP), INTENT  (OUT), DIMENSION(nx,ny,ksta:kend) :: lamb
      REAL   (KIND=GP)                                          :: sa,sb,sc,sd,se,sf
      REAL   (KIND=GP)                                          :: a,b,c,d
      REAL   (KIND=GP)                                          :: del0,del1,ll(3),lmax
      COMPLEX(KIND=GP)                                          :: CC,u(3)
      COMPLEX(KIND=GP)                                          :: D0,D1

!
      INTEGER                                                   :: i,j,k,l

      u(1) = cmplx(1.0_GP , 0.0_GP)
      u(2) = cmplx(-0.5_GP, 0.5_GP*sqrt(3.0_GP))
      u(3) = cmplx(-0.5_GP,-0.5_GP*sqrt(3.0_GP))
!$omp parallel do if (kend-ksta.ge.nth) private (j,i,l,sa,sb,sc,sd,se,sf,a,b,c,d,del0,del1,D0,D1,CC,lmax)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i,l,sa,sb,sc,sd,se,sf,a,b,c,d,del0,del1,D0,D1,CC,lmax)
        DO j = 1,ny
          DO i = 1,nx
            sa = S11(i,j,k); sb = S12(i,j,k); sc = S13(i,j,k); 
            sd = S22(i,j,k); se = S23(i,j,k); sf = -(sa + sd);
            a    = 1.0_GP
            b    = 0.0_GP ! sa + sd + sf
            c    = sa*sd+sa*sf+sd*sf-sb**2 - sc**2 - se**2
            d    = sa*se**2 + sf*sb**2 + sd*sc**2 - sa*sd*sf - 2.0_GP*sb*sc*se
!           del0 = b**2 - 3.0_GP*a*c
!           del1 = 2.0_GP*b**3 - 9.0_GP*a*b*c + 27.0_GP*d*a**2
            del0 = -3.0_GP*a*c
            del1 =  27.0_GP*d*a**2
            D0   = cmplx(del0,0)
            D1   = cmplx(del1,0)
            CC   = (0.5*D1 + 0.5*sqrt(D1**2 -4.0_GP*D0**3) )**(1.0_GP/3.0_GP)
            lmax = 0.0_GP
            DO l = 1,3
!             ll(l) = real(-( b + u(l)*CC + D0/(u(l)*CC))/(3.0_GP*a),kind=GP)
              ll(l) = real(-( u(l)*CC + D0/(u(l)*CC))/(3.0_GP*a),kind=GP)
              lmax = max(lmax,abs(ll(l)))
            ENDDO
!if ( i.eq.10 .and. j.eq.10 .and. k.gt.10 .and. k.lt.15) then
!!write(*,*)'Eigen: sa=',sa,' sb=',sb,' sc=',sc,' sd=',sd, 'se=',se
!write(*,*)'Eigen: u=',u,' CC=', CC,' del0=',del0,' del1=',del1,' lmax=',lmax, ' ll=',ll
!endif

            lamb(i,j,k) = lmax
          END DO
        END DO
      END DO


      END SUBROUTINE EigenValMax
!-----------------------------------------------------------------
!-----------------------------------------------------------------


      SUBROUTINE EigenSolveMax(S11,S12,S13,S22,S23,lamb,evx,evy,evz)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!
! Computes the max |eigenvalue| field for strain rate tensor, whos
! required components are specified. Also computes the compnents of
! the eigenvector corresponding to this max eigenvalue at each point.
! If at a point the strain rate tensor is not of rank 2, meaning 3
! distinct e-values, and 3 e-vectors, then return (0,0,0) for the 
! e-vector of the max e-value.
!
! Parameters
!     SIJ  : strain rate components (real)
!     lamb : returned eigenvalue field
!     evx,
!     evy,
!     evz  : field of eigenvector components corresp. to lamb, returned
!
      USE fprecision
      USE commtypes
      USE kes
      USE grid
      USE mpivars
      USE threads
      IMPLICIT NONE

      REAL   (KIND=GP), INTENT   (IN), DIMENSION(nx,ny,ksta:kend) :: S11,S12,S13,S22,S23
      REAL   (KIND=GP), INTENT  (OUT), DIMENSION(nx,ny,ksta:kend) :: lamb,evx,evy,evz
      REAL   (KIND=GP)                                          :: s   (3,3),m1  (3,3),m2  (3,3),m3  (3,3)
      REAL   (KIND=GP)                                          :: eval(3)  ,evec(3)  ,row1(3)  ,row2(3)
      REAL   (KIND=GP)                                          :: ev1 (3)  ,ev2 (3)
      INTEGER                                                   :: i,j,k,l,m
      INTEGER                                                   :: GetRank
!$    INTEGER,EXTERNAL                                          :: omp_get_thread_num

!$omp parallel do if (kend-ksta.ge.nth) private (j,i,s,m1,m2,m3,evec,eval,ev1,ev2,row1,row2)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i,s,m1,m2,m3,evec,eval,ev1,ev2,row1,row2)
        DO j = 1,ny
          DO i = 1,nx
            s(1,1)=S11(i,j,k); s(1,2)=S12(i,j,k); s(1,3)=S13(i,j,k)
            s(2,1)=s(1,2); s(2,2)=S22(i,j,k); s(2,3)=S23(i,j,k)
            s(3,1)=s(1,3); s(3,2)=s(2,3); s(3,3)=-(s(1,1)+s(2,2))
            CALL GetRoots(s,eval)
            IF ( GetRank(s).LT.2 ) THEN
              evec = 0.0_GP
              CYCLE
            ENDIF
            m1 = s;
            m1(1,1) = s(1,1)-eval(1); m1(2,2) = s(2,2)-eval(1); m1(3,3) = s(3,3)-eval(1)
            IF ( GetRank(m1).LT.2 ) THEN
              evec = 0.0_GP
              CYCLE
            ENDIF
            row1(1:3) = m1(1,1:3)
            row2(1:3) = m1(2,1:3)
            CALL GetComplement1(row1,row2,ev1)

            m2 = s;
            m2(1,1) = s(1,1)-eval(2); m2(2,2) = s(2,2)-eval(2); m2(3,3) = s(3,3)-eval(2)
            IF ( GetRank(m2).LT.2 ) THEN
              evec(1:3) = 0.0_GP
              CYCLE
            ENDIF
            row1(1:3) = m2(1,1:3)
            row2(1:3) = m2(2,1:3)
            CALL GetComplement1(row1,row2,ev2)

            m3 = s;
            m3(1,1) = s(1,1)-eval(3); m3(2,2) = s(2,2)-eval(3); m3(3,3) = s(3,3)-eval(3)
            IF ( GetRank(m3).LT.2 ) THEN
              evec(1:3) = 0.0_GP
              CYCLE
            ENDIF
            CALL GetComplement1(ev1,ev2,evec)
            evx(i,j,k)=evec(1); evy(i,j,k)=evec(2); evz(i,j,k)=evec(3); lamb(i,j,k)=abs(eval(3))
          ENDDO
        ENDDO
      ENDDO


      END SUBROUTINE EigenSolveMax
!-----------------------------------------------------------------
!-----------------------------------------------------------------


      RECURSIVE SUBROUTINE GetComplement1(u,v,w)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!
! Computes the normalized cross product of 3-vectors u,v, and
! stores result in w
!
! Parameters
!     u,v  : input 3-vectors
!     w    : returned result
!
      USE fprecision
      USE ali
      IMPLICIT NONE

      REAL   (KIND=GP), INTENT(INOUT), DIMENSION(3) :: u,v
      REAL   (KIND=GP), INTENT(INOUT), DIMENSION(3) :: w
      REAL   (KIND=GP)                              :: cmag
      INTEGER                                       :: j  

      ! compute cross product:
      w(1) = u(2)*v(3) - u(3)*v(2)
      w(2) = u(3)*v(1) - u(1)*v(3)
      w(3) = u(1)*v(2) - u(2)*v(1)

      cmag = 0.0_GP
      DO j = 1,3
        cmag = cmag + w(j)*w(j)
      ENDDO
      cmag = sqrt(cmag) 
      IF ( cmag.GT.tiny ) cmag = 1.0_GP/cmag
      DO j = 1,3
        w(j) = w(j)*cmag
      ENDDO

      END SUBROUTINE GetComplement1
!-----------------------------------------------------------------
!-----------------------------------------------------------------

      RECURSIVE INTEGER FUNCTION GetRank(m)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!
! Computes the rank of the 3x3 matrix, m, and returns it
!
! Parameters
!     m  : input 3x3 matrix
!
      USE fprecision
      USE ali
      IMPLICIT NONE

      REAL   (KIND=GP), INTENT(INOUT), DIMENSION(3,3) :: m
      REAL   (KIND=GP)                                :: fabs,fmax,fmaxi,fsave
      INTEGER                                         :: i,j,maxrow,maxcol

      maxrow = 0
      maxcol = 0

      fmax = 0.0_GP
      DO i = 1, 3
        DO j = i, 3
          fabs = abs(m(i,j))
          IF ( fabs .GT. fmax ) THEN
            fmax = fabs
            maxrow = i
            maxcol = j
          ENDIF
        ENDDO
      ENDDO

      IF ( fmax .LT. tiny ) THEN
        GetRank = 0 ! e-value has multiplicity three
        RETURN
      ENDIF
     
      ! rank is >= 1. Swap row containing the max-magnitude entry
      ! with wow 1:

      IF ( maxrow .NE. 1 ) THEN
       DO j = 1, 3
          DO i = 1, 3
            fsave = m(1,j)
            m(1,j) = m(maxrow,j)
            m(maxrow,j) = fsave
          ENDDO
        ENDDO
      ENDIF

      ! row-reduce matrix:

      ! scale row 1 to generate a 1-valued pivot:
      fmaxi = 0.0_GP
      IF ( abs(m(1,maxcol)).GT.tiny ) fmaxi = 1.0_GP / m(1,maxcol)
      m(1,1) = m(1,1)*fmaxi
      m(1,2) = m(1,2)*fmaxi
      m(1,3) = m(1,3)*fmaxi
     
      ! eliminate maxcol column entries in rows 2 & 2:
      IF ( maxcol .EQ. 1 ) THEN 
        m(2,2) = m(2,2) - m(2,1)*m(1,2)
        m(2,3) = m(2,3) - m(2,1)*m(1,3)
        m(3,2) = m(3,2) - m(3,1)*m(1,2)
        m(3,3) = m(3,3) - m(3,1)*m(1,3)
        m(2,1) = 0.0_GP
        m(3,1) = 0.0_GP
      ELSE IF ( maxcol .EQ. 2 ) THEN
        m(2,1) = m(2,1) - m(2,2)*m(1,1)
        m(2,3) = m(2,3) - m(2,2)*m(1,3)
        m(3,1) = m(3,1) - m(3,2)*m(1,1)
        m(3,3) = m(3,3) - m(3,2)*m(1,3)
        m(2,2) = 0.0_GP
        m(3,2) = 0.0_GP
      ELSE
        m(2,1) = m(2,1) - m(2,3)*m(1,1)
        m(2,2) = m(2,2) - m(2,3)*m(1,2)
        m(3,1) = m(3,1) - m(3,3)*m(1,1)
        m(3,2) = m(3,2) - m(3,3)*m(1,2)
        m(2,3) = 0.0_GP
        m(3,3) = 0.0_GP
      ENDIF

      ! compute max-magnitude entry of the last 2 rows
      ! of row-reduced matrix:
      fmax = -1.0_GP
      maxrow = 0
      maxcol = 0

      DO j = 1, 3
        DO i = 2, 3
          fabs = abs(m(i,j))
          IF ( fabs .GT. fmax ) THEN
            fmax = fabs
            maxrow = i
            maxcol = j
          ENDIF
        ENDDO
      ENDDO

      IF ( fmax .lt. tiny ) THEN
        GetRank = 1 ! e-value has multiplicity 2
        RETURN
      ENDIF

      ! if row 2 has the max-magnitude entry, swap it with 
      ! row 1:
      IF ( maxrow .EQ. 3 ) THEN
        DO j = 1, 3
          fsave = m(2,j)
          m(2,j) = m(3,j)
          m(3,j) = fsave
        ENDDO
      ENDIF

      ! scale row 1 to generate a 1-vaued point:
      fmaxi = 0.0_GP
      IF ( abs(m(2,maxcol)).GT.tiny ) fmaxi = 1.0_GP / m(2,maxcol)
      m(2,1) = m(2,1)*fmaxi
      m(2,2) = m(2,2)*fmaxi
      m(2,3) = m(2,3)*fmaxi
 
      GetRank = 2 ! e-value has multiplicity 1

      END FUNCTION GetRank
!-----------------------------------------------------------------
!-----------------------------------------------------------------

      RECURSIVE SUBROUTINE GetRoots(m,eval)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!
! Computes the ordered eigen values of input 3x3 symmetric matrix, m,
! returns in eval vector from smallest to largest. It is assumed that
! the input matrix, m, is trace free, which applies to the straing
! rate matrix for incompressible flows only!
!
! Parameters
!     m    : input 3x3 matrix
!     eval : eigenvalues from smallest to largest returned 
!
      USE fprecision
      IMPLICIT NONE

      REAL   (KIND=GP), INTENT(INOUT), DIMENSION(3,3) :: m
      REAL   (KIND=GP), INTENT(INOUT), DIMENSION  (3) :: eval
      REAL   (KIND=GP)                                :: sa,sb,sc,sd,se,sf
      REAL   (KIND=GP)                                :: a,b,c,d
      REAL   (KIND=GP)                                :: del0,del1
      COMPLEX(KIND=GP)                                :: CC,u(3)
      COMPLEX(KIND=GP)                                :: D0,D1
      INTEGER                                         :: j,l

      u(1) = cmplx(1.0_GP , 0.0_GP)
      u(2) = cmplx(-0.5_GP, 0.5_GP*sqrt(3.0_GP))
      u(3) = cmplx(-0.5_GP,-0.5_GP*sqrt(3.0_GP))

      sa   = m(1,1); sb = m(1,2); sc = m(1,3); 
      sd   = m(2,2); se = m(2,3); sf = -(sa+sd)
      a    = 1.0_GP
      b    = 0.0_GP ! sa + sd + sf; trace-free
      c    = sa*sd+sa*sf+sd*sf-sb**2 - sc**2 - se**2
      d    = sa*se**2 + sf*sb**2 + sd*sc**2 - sa*sd*sf - 2.0_GP*sb*sc*se
!     del0 = b**2 - 3.0_GP*a*c
!     del1 = 2.0_GP*b**3 - 9.0_GP*a*b*c + 27.0_GP*d*a**2
      del0 = -3.0_GP*a*c
      del1 =  27.0_GP*d*a**2
      D0   = cmplx(del0,0)
      D1   = cmplx(del1,0)
      CC   = (0.5*D1 + 0.5*sqrt(D1**2 -4.0_GP*D0**3) )**(1.0_GP/3.0_GP)
      DO l = 1,3
!       eval(l) = real(-( b + u(l)*CC + D0/(u(l)*CC))/(3.0_GP*a),kind=GP)
        eval(l) = real(-( u(l)*CC + D0/(u(l)*CC))/(3.0_GP*a),kind=GP)
      ENDDO

      ! reorder from smallest to largest:in magnitude
      DO j = 1,2
        DO l = j+1,3
          IF ( abs(eval(l)).LT.abs(eval(j)) ) THEN
            a       = eval(j)
            eval(j) = eval(l)
            eval(l) = a
          ENDIF
        ENDDO
      ENDDO

      END SUBROUTINE GetRoots
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!
!
      SUBROUTINE DoKPDF(S11,S12,S13,S22,S23,vx,vy,vz,th,lambda,bvfreq,omega,nu,kappa,ext, &
                            odir,nbins,dolog,ctmp,vtmp,rtmp,rtmp1,rtmp2,beta,pv,rig,rign,kin,planio,iobin)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!
! Computes the joint PDF of energy dissipation/(2*nu) and ( enstrophy
! lambda, helicity, and other quantities), and output to files. 1d PDFs of dissipation, 
! enstrophy density, lambda, and helicity are also written.
!
! Note: after this call, the input data should be expected to be overwritten.
!
! Parameters
!     SIJ   : strain rate components (real)
!     vx,
!     vy,
!     vz    : complex velocities
!     lambda: e-value field
!     bvfreq: Brunt-Vaisalla frequency
!     omega : rotation freq vector
!     nu    : viscosity
!     kappa : scalar dissipation
!     ext   : string time index
!     odir  : output directory
!     nbins : 2-elem array with no. bins for (enstrophy,energy diss)
!     dolog : flag to do (1) or not (0) logs of ranges when computing bins
!     kin   : 0: compute no pdfs
!             1: compute 1d pdfs only
!             2: compute 2d pdfs only
!             3: compute 1d and 2d pdfs 
!     ctmp  : complex temp array of size vx,vy,vz
!     vtmp  : complex temp array of size vx,vy,vz
!     rtmp  : real temp array of size lambda
!     rtmp1 : real temp array of size lambda
!     rtmp2 : real temp array of size lambda
!     beta  : real temp array of size lambda
!
      USE fprecision
      USE commtypes
      USE kes
      USE grid
      USE mpivars
      USE threads
      USE fft
      USE var
      USE fftplans
      USE ali
      USE gutils
      USE iovar
      USE iompi
      IMPLICIT NONE

      REAL   (KIND=GP), INTENT(INOUT), DIMENSION(nx,ny,ksta:kend) :: lambda,rtmp,rtmp1,rtmp2,beta,S11,S12,S13,S22,S23,pv,rig, rign
      REAL   (KIND=GP), INTENT   (IN)                           :: bvfreq,omega(3),nu,kappa
      COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(nx,ny,ista:iend) :: ctmp,vtmp
      COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(nx,ny,ista:iend) :: vx,vy,vz,th
      REAL   (KIND=GP)                                          :: fact,fmin(2),fmax(2),xnorm,xnormi,xnormn
      REAL   (KIND=GP)                                          :: as11,as12,as13,as22,as23,as33
      REAL   (KIND=GP)                                          :: ss11,ss12,ss13,ss22,ss23,ss33
      REAL   (KIND=GP)                                          :: sf11,sf12,sf13,sf22,sf23,sf33
      REAL   (KIND=GP)                                          :: sg11,sg12,sg13,sg22,sg23,sg33
      REAL   (KIND=GP)                                          :: sw11,sw12,sw13,sw22,sw23,sw33
      DOUBLE PRECISION                                          :: s2,s3,s4,s5,s6
      REAL   (KIND=GP)                                          :: ria,dupa,upa,dthdxa,dthdya,dthdza,tha,ktra,ptra,wta,wupa,betaa,pva,wgta,wgt2pva
      REAL   (KIND=GP)                                          :: rif,dupf,upf,dthdxf,dthdyf,dthdzf,thf,ktrf,ptrf,wtf,wupf,betaf,pvf,wgtf,wgt2pvf
      REAL   (KIND=GP)                                          :: ris,dups,ups,dthdxs,dthdys,dthdzs,ths,ktrs,ptrs,wts,wups,betas,pvs,wgts,wgt2pvs
      REAL   (KIND=GP)                                          :: rid,dupg,upg,dthdxg,dthdyg,dthdzg,thg,ktrg,ptrg,wtg,wupg,betag,pvg,wgtg,wgt2pvg
      REAL   (KIND=GP)                                          :: riw,dupw,upw,dthdxw,dthdyw,dthdzw,thw,ktrw,ptrw,wtw,wupw,betaw,pvw,wgtw,wgt2pvw
      REAL   (KIND=GP)                                          :: omegax,omegay,omegaz
      REAL   (KIND=GP)                                          :: va1,va2,va3
      REAL   (KIND=GP)                                          :: vs1,vs2,vs3
      REAL   (KIND=GP)                                          :: vf1,vf2,vf3
      REAL   (KIND=GP)                                          :: vg1,vg2,vg3
      REAL   (KIND=GP)                                          :: wa1,wa2,wa3
      REAL   (KIND=GP)                                          :: vw1,vw2,vw3
      REAL   (KIND=GP)                                          :: ws1,ws2,ws3
      REAL   (KIND=GP)                                          :: wg1,wg2,wg3
      REAL   (KIND=GP)                                          :: wf1,wf2,wf3
      REAL   (KIND=GP)                                          :: ww1,ww2,ww3
      REAL   (KIND=GP)                                          :: adiss,fdiss,sdiss,gdiss,wdiss
      REAL   (KIND=GP)                                          :: fdisst,sdisst,gdisst,wdisst
      REAL   (KIND=GP)                                          :: aenst,fenst,senst,genst,wenst
      REAL   (KIND=GP)                                          :: aheli,fheli,sheli,gheli,wheli
      REAL   (KIND=GP)                                          :: alamb,flamb,slamb,glamb,wlamb
      REAL   (KIND=GP)                                          :: den  ,rmi
      REAL   (KIND=GP)                                          :: alpha
      INTEGER         , INTENT   (IN)                           :: dolog,iobin,kin,nbins(2)
      INTEGER                                                   :: i,j,k,nin,sr,indtime
      CHARACTER(len=*), INTENT   (IN)                           :: odir
      CHARACTER(len=*), INTENT   (IN)                           :: ext
      CHARACTER(len=1024)                                       :: fnout
      TYPE(IOPLAN)    , INTENT(INOUT)                           :: planio

      rmi = 1000. ! max |Richardson| to allow

      omegax = omega(1)
      omegay = omega(2)
      omegaz = omega(3)

if ( myrank.eq.0 ) write(*,*) 'DoDPDF: 0'
      read(ext,*) indtime ! convert time index string to integer

if ( myrank.eq.0 ) write(*,*) 'DoDPDF: 1'
      nin = nx*ny*(kend-ksta+1)
      xnormn = 1.0_GP/ real(nx*ny*nz,kind=GP)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
        DO j = 1,ny
          DO k = 1,nz
            vx(k,j,i) = vx(k,j,i)*xnormn
            vy(k,j,i) = vy(k,j,i)*xnormn
            vz(k,j,i) = vz(k,j,i)*xnormn
#if defined(SCALAR_)
            th(k,j,i) = th(k,j,i)*xnormn
#endif
          END DO
        END DO
      END DO

! At this point forward, complex velocity components are normalized!
if ( myrank.eq.0 ) write(*,*) 'DoDPDF: 2'

      ! compute S33 for flatness/skewness:
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
        DO j = 1,ny
          DO i = 1,nx
            rtmp(i,j,k) =  -(S11(i,j,k)+S22(i,j,k)) 
          ENDDO
        ENDDO
      ENDDO
if ( myrank.eq.0 ) write(*,*) 'DoDPDF: 3'

      CALL skewflat(S11 ,nx,ny,nz,as11,ss11,sf11,sg11,sw11,s2,s3,s4,s5,s6)
      CALL skewflat(S12 ,nx,ny,nz,as12,ss12,sf12,sg12,sw12,s2,s3,s4,s5,s6)
      CALL skewflat(S13 ,nx,ny,nz,as13,ss13,sf13,sg13,sw13,s2,s3,s4,s5,s6)
      CALL skewflat(S22 ,nx,ny,nz,as22,ss22,sf22,sg22,sw22,s2,s3,s4,s5,s6)
      CALL skewflat(S23 ,nx,ny,nz,as23,ss23,sf23,sg23,sw23,s2,s3,s4,s5,s6)
      CALL skewflat(rtmp,nx,ny,nz,as33,ss33,sf33,sg33,sw33,s2,s3,s4,s5,s6)
      ! Compute, write, 1d Sij pdfs:
   
if ( myrank.eq.0 ) write(*,*) 'DoDPDF: 4'
      If ( ibits(kin,0,1).EQ.1 ) THEN
      fnout = trim(odir) // '/' // 's11pdf.' // ext // '.txt'
      CALL dopdfr(S11 ,nx,ny,nz,fnout,nbins(1),0,fmin(1),fmax(1),0) 
write(*,*) 'S11 = ', S11(20:50,10,10)
stop
      fnout = trim(odir) // '/' // 's12pdf.' // ext // '.txt'
      CALL dopdfr(S12 ,nx,ny,nz,fnout,nbins(1),0,fmin(1),fmax(1),0) 
      fnout = trim(odir) // '/' // 's13pdf.' // ext // '.txt'
      CALL dopdfr(S13 ,nx,ny,nz,fnout,nbins(1),0,fmin(1),fmax(1),0) 
      fnout = trim(odir) // '/' // 's22pdf.' // ext // '.txt'
      CALL dopdfr(S22 ,nx,ny,nz,fnout,nbins(1),0,fmin(1),fmax(1),0) 
      fnout = trim(odir) // '/' // 's23pdf.' // ext // '.txt'
      CALL dopdfr(S23 ,nx,ny,nz,fnout,nbins(1),0,fmin(1),fmax(1),0) 
      fnout = trim(odir) // '/' // 's33pdf.' // ext // '.txt'
      CALL dopdfr(rtmp,nx,ny,nz,fnout,nbins(1),0,fmin(1),fmax(1),0) 
      ENDIF

      ! Compute normalized energy dissipation field/2nu, store in ! S11:
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
        DO j = 1,ny
          DO i = 1,nx
            S11(i,j,k) = 2.0_GP*( S11(i,j,k)*S11(i,j,k) &
                                 +S22(i,j,k)*S22(i,j,k) &
                                 +S22(i,j,k)*S11(i,j,k) &
                                 +S12(i,j,k)*S12(i,j,k) &
                                 +S13(i,j,k)*S13(i,j,k) &
                                 +S23(i,j,k)*S23(i,j,k) &
                                )*nu
          ENDDO
        ENDDO
      ENDDO

      CALL skewflat(S11   ,nx,ny,nz,adiss,sdiss,fdiss,gdiss,wdiss,s2,s3,s4,s5,s6) ! dissipation
!write(*,*)'DoKPDF: diss: s2=',s2,' s3=',s3,' s4,s5=',s4,s5,' sdiss=',sdiss,' fdiss=',fdiss
      CALL skewflat(lambda,nx,ny,nz,alamb,slamb,flamb,glamb,wlamb,s2,s3,s4,s5,s6) ! lambda
!write(*,*)'DoKPDF: lamb: s2=',s2,' s3=',s3,' s4,s5=',s4,s5,' slamb=',slamb,' flamb=',flamb
  
      ! Compute and write out dissipation spectrum:
      rtmp1 = S11
      CALL fftp3d_real_to_complex(planrc,rtmp1,ctmp,MPI_COMM_WORLD)
      fnout = trim(odir) // '/' // 'kdissv.' // ext // '.txt'
      CALL pspectrum(ctmp,fnout,int(kmax))

      ! Compute joint PDF for energy diss and lambda (order 
      ! switched, so that lambda is on x-axis, and energy diss on y axis:
      If ( ibits(kin,1,1).EQ.1 ) THEN
      fnout = trim(odir) // '/' // 'jpdf_lamb_diss_log00.' // ext // '.txt'
      CALL dojpdfr(lambda,'lambda',S11,'diss',nx,ny,nz,fnout,nbins,[0,0],fmin,fmax,[0,0])
      fnout = trim(odir) // '/' // 'jpdf_lamb_diss_log11.' // ext // '.txt'
      CALL dojpdfr(lambda,'lambda',S11,'diss',nx,ny,nz,fnout,nbins,[0,0],fmin,fmax,[1,1])
      ENDIF
      If ( ibits(kin,0,1).EQ.1 ) THEN
      ! Do 1d and jooint pdfs for diss and lambda:
      fnout = trim(odir) // '/' // 'disspdf_log.' // ext // '.txt'
      CALL dopdfr(S11   ,nx,ny,nz,fnout,nbins(1),0,fmin(2),fmax(2),1) 
      fnout = trim(odir) // '/' // 'disspdf.' // ext // '.txt'
      CALL dopdfr(S11   ,nx,ny,nz,fnout,nbins(1),0,fmin(2),fmax(2),0) 
      ! Compute, write, 1d lambda pdfs:
      fnout = trim(odir) // '/' // 'lambpdf_log.' // ext // '.txt'
      CALL dopdfr(lambda,nx,ny,nz,fnout,nbins(1),0,fmin(1),fmax(1),1) 
      fnout = trim(odir) // '/' // 'lambpdf.' // ext // '.txt'
      CALL dopdfr(lambda,nx,ny,nz,fnout,nbins(1),0,fmin(1),fmax(1),0) 
      ENDIF

if (myrank.eq.0) write(*,*)'DoKPDF: finding beta...'

      ! Compute statistics for beta = [epsilon_V/(2 nu)] / u^3
      ! Must restore dimensions by computing r * 2 nu * L_int
      vtmp = vx
      CALL fftp3d_complex_to_real(plancr,vtmp,rtmp,MPI_COMM_WORLD)
      vtmp = vy
      CALL fftp3d_complex_to_real(plancr,vtmp,rtmp1,MPI_COMM_WORLD)
      vtmp = vz
      CALL fftp3d_complex_to_real(plancr,vtmp,rtmp2,MPI_COMM_WORLD)
if (myrank.eq.0) write(*,*)'DoKPDF: beta FTs done...'

!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
        DO j = 1,ny
          DO i = 1,nx
            rtmp(i,j,k) = rtmp(i,j,k)**2 + rtmp1(i,j,k)**2 + rtmp2(i,j,k)**2
            beta(i,j,k) = S11(i,j,k) / rtmp(i,j,k)**(1.5) ! Sij Sij / u^3
          ENDDO
        ENDDO
      ENDDO
      If ( ibits(kin,0,1).EQ.1 ) THEN
if (myrank.eq.0) write(*,*)'DoKPDF: do beta pdf:'
      fnout = trim(odir) // '/' // 'betapdf.' // ext // '.txt'
      CALL dopdfr(beta,nx,ny,nz,fnout,nbins(1),0,fmin(1),fmax(1),0) 
      fnout = trim(odir) // '/' // 'betapdf_log.' // ext // '.txt'
      CALL dopdfr(beta,nx,ny,nz,fnout,nbins(1),0,fmin(1),fmax(1),1) 
      ENDIF
if (myrank.eq.0) write(*,*)'DoKPDF: do beta skewnews, flatness:'
      CALL skewflat(beta ,nx,ny,nz,betaa,betas,betaf,betag,betaw,s2,s3,s4,s5,s6)       ! beta = diss/u^3

      ! Compute enstrophy density,helicity field, v^2; store in S22, S13, S12:
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
        DO j = 1,ny
          DO i = 1,nx
            S22 (i,j,k) = 0.0_GP ! Enstrophy density: omega^2
            S12 (i,j,k) = 0.0_GP ! v^2
            S13 (i,j,k) = 0.0_GP ! v.omega
          ENDDO
        ENDDO
      ENDDO

!     fact  = 1.0_GP/real(nx*ny*nz,kind=GP)**2 ! no  longer needed; is in v_i now
      CALL rotor3(vy,vz,ctmp,1)
      CALL fftp3d_complex_to_real(plancr,ctmp,S23   ,MPI_COMM_WORLD)
      vtmp = vx
      CALL fftp3d_complex_to_real(plancr,vtmp,rtmp,MPI_COMM_WORLD)
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
        DO j = 1,ny
          DO i = 1,nx
            S22 (i,j,k) = S22 (i,j,k) + S23   (i,j,k)**2         ! Enstrophy dens
            S13 (i,j,k) = S13 (i,j,k) + rtmp  (i,j,k)*S23(i,j,k) ! v.omega
            S12 (i,j,k) = S12 (i,j,k) + rtmp  (i,j,k)**2         ! v^2
          ENDDO
        ENDDO
      ENDDO
      CALL skewflat(rtmp,nx,ny,nz,va1,vs1,vf1,vg1,vw1,s2,s3,s4,s5,s6) ! v1 
      CALL skewflat(S23 ,nx,ny,nz,wa1,ws1,wf1,wg1,ww1,s2,s3,s4,s5,s6) ! omega1 

      ! Compute, write, 1d vx, \omega_i pdfs:
      If ( ibits(kin,0,1).EQ.1 ) THEN
      fnout = trim(odir) // '/' // 'v1pdf.' // ext // '.txt'
      CALL dopdfr(rtmp,nx,ny,nz,fnout,nbins(1),0,fmin(1),fmax(1),0) 
      fnout = trim(odir) // '/' // 'w1pdf.' // ext // '.txt'
      CALL dopdfr(S23 ,nx,ny,nz,fnout,nbins(1),0,fmin(1),fmax(1),0) 
      ENDIF

      CALL rotor3(vz,vx,ctmp,2)
      CALL fftp3d_complex_to_real(plancr,ctmp,S23   ,MPI_COMM_WORLD)
      vtmp = vy
      CALL fftp3d_complex_to_real(plancr,vtmp,rtmp  ,MPI_COMM_WORLD)
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
        DO j = 1,ny
          DO i = 1,nx
            S22 (i,j,k) = S22 (i,j,k) + S23   (i,j,k)**2         ! Enstrophy dens
            S13 (i,j,k) = S13 (i,j,k) + rtmp  (i,j,k)*S23(i,j,k) ! v.omega
            S12 (i,j,k) = S12 (i,j,k) + rtmp  (i,j,k)**2         ! v^2
          ENDDO
        ENDDO
      ENDDO
      CALL skewflat(rtmp,nx,ny,nz,va2,vs2,vf2,vg2,vw2,s2,s3,s4,s5,s6) ! v2 
      CALL skewflat(S23 ,nx,ny,nz,wa2,ws2,wf2,wg2,ww2,s2,s3,s4,s5,s6) ! omega2 
      ! Compute, write, 1d vy, \omega_y pdfs:
      If ( ibits(kin,0,1).EQ.1 ) THEN
      fnout = trim(odir) // '/' // 'v2pdf.' // ext // '.txt'
      CALL dopdfr(rtmp,nx,ny,nz,fnout,nbins(1),0,fmin(1),fmax(1),0) 
      fnout = trim(odir) // '/' // 'w2pdf.' // ext // '.txt'
      CALL dopdfr(S23 ,nx,ny,nz,fnout,nbins(1),0,fmin(1),fmax(1),0) 
      ENDIF

      CALL rotor3(vx,vy,ctmp,3)
      CALL fftp3d_complex_to_real(plancr,ctmp,S23,MPI_COMM_WORLD)
      vtmp = vz
      CALL fftp3d_complex_to_real(plancr,vtmp,rtmp,MPI_COMM_WORLD)
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
        DO j = 1,ny
          DO i = 1,nx
            S22 (i,j,k) = S22 (i,j,k) + S23   (i,j,k)**2         ! Enstrophy dens
            S13 (i,j,k) = S13 (i,j,k) + rtmp  (i,j,k)*S23(i,j,k) ! v.omega
            S12 (i,j,k) = S12 (i,j,k) + rtmp  (i,j,k)**2         ! v^2
          ENDDO
        ENDDO
      ENDDO
      CALL skewflat(rtmp,nx,ny,nz,va3,vs3,vf3,vg3,vw3,s2,s3,s4,s5,s6)       ! v3 
      CALL skewflat(S23 ,nx,ny,nz,wa3,ws3,wf3,wg3,ww3,s2,s3,s4,s5,s6)       ! omega3 
      CALL skewflat(S22 ,nx,ny,nz,aenst,senst,fenst,genst,wenst,s2,s3,s4,s5,s6) ! enstrophy density
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
        DO j = 1,ny
          DO i = 1,nx
            S22 (i,j,k) = kappa*S22 (i,j,k) + S11(i,j,k)      ! total dissipation
          ENDDO
        ENDDO
      ENDDO
      CALL skewflat(S22 ,nx,ny,nz,adiss,sdiss,fdiss,gdiss,wdiss,s2,s3,s4,s5,s6) ! total dissipation
      If ( ibits(kin,0,1).EQ.1 ) THEN
      fnout = trim(odir) // '/' // 'disstpdf.' // ext // '.txt'
      CALL dopdfr(S22,nx,ny,nz,fnout,nbins(1),0,fmin(1),fmax(1),0) 
      fnout = trim(odir) // '/' // 'disstpdf_log.' // ext // '.txt'
      CALL dopdfr(S22,nx,ny,nz,fnout,nbins(1),0,fmin(1),fmax(1),1) 
      ENDIF

      ! Compute, write, 1d v_z, \omega_z pdfs:
      If ( ibits(kin,0,1).EQ.1 ) THEN
      fnout = trim(odir) // '/' // 'v3pdf.' // ext // '.txt'
      CALL dopdfr(rtmp,nx,ny,nz,fnout,nbins(1),0,fmin(1),fmax(1),0) 
      fnout = trim(odir) // '/' // 'w3pdf.' // ext // '.txt'
      CALL dopdfr(S23 ,nx,ny,nz,fnout,nbins(1),0,fmin(1),fmax(1),0) 
      ENDIF

!     Compute relative helicity:
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
        DO j = 1,ny
          DO i = 1,nx
            sr = sign(1.0,S13(i,j,k));
            xnorm = sqrt( S22(i,j,k)*S12(i,j,k) )
            xnormi = 0.0_GP
            IF ( xnorm .GT. tiny ) THEN
              xnormi = 1.0_GP/xnorm
            ENDIF
            S13 (i,j,k) = S13(i,j,k) * xnormi
            S13 (i,j,k) = sr*min(abs(S13(i,j,k)),1.0_GP)
          ENDDO
        ENDDO
      ENDDO
      CALL skewflat(S13,nx,ny,nz,aheli,sheli,fheli,gheli,wheli,s2,s3,s4,s5,s6)    ! v.omega/rel. helicity

! Stats for u_perp
! NOTE: S12 and S23 overwritten here:
      ctmp = vx;
      CALL fftp3d_complex_to_real(plancr,ctmp,S12,MPI_COMM_WORLD)
      ctmp = vy;
      CALL fftp3d_complex_to_real(plancr,ctmp,S23,MPI_COMM_WORLD)
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
        DO j = 1,ny
          DO i = 1,nx
            rtmp(i,j,k) =  sqrt(S12(i,j,k)**2 + S23(i,j,k)**2)
          ENDDO
        ENDDO
      ENDDO
      ! Compute dup/dz
      rtmp1 = rtmp;
      CALL fftp3d_real_to_complex(planrc,rtmp1,ctmp,MPI_COMM_WORLD)
      CALL derivk3(ctmp, vtmp, 3)
      vtmp = vtmp * xnormn
      CALL fftp3d_complex_to_real(plancr,vtmp,rtmp1,MPI_COMM_WORLD)
      If ( ibits(kin,0,1).EQ.1 ) THEN
      fnout = trim(odir) // '/' // 'uppdf.' // ext // '.txt'
      CALL dopdfr(rtmp ,nx,ny,nz,fnout,nbins(1),0,fmin(1),fmax(1),0) 
      fnout = trim(odir) // '/' // 'dupdzpdf.' // ext // '.txt'
      CALL dopdfr(rtmp1,nx,ny,nz,fnout,nbins(1),0,fmin(1),fmax(1),0) 
      ENDIF
      CALL skewflat(rtmp ,nx,ny,nz, upa, ups, upf,upg,upw,s2,s3,s4,s5,s6)       ! up
      CALL skewflat(rtmp1,nx,ny,nz,dupa,dups,dupf,dupg,dupw,s2,s3,s4,s5,s6)       ! dup/dz

      ! Compute stats for u_z u_perp = w up:
      ! NOTE: rtmp should still contain up:
      vtmp = vz
      CALL fftp3d_complex_to_real(plancr,vtmp,rtmp1,MPI_COMM_WORLD) !vz
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
        DO j = 1,ny
          DO i = 1,nx
            rtmp(i,j,k) =  rtmp1(i,j,k) * rtmp(i,j,k) ! vz * u_perp 
          ENDDO
        ENDDO
      ENDDO
      fnout = trim(odir) // '/' // 'wuppdf.' // ext // '.txt'
      CALL dopdfr(rtmp ,nx,ny,nz,fnout,nbins(1),0,fmin(1),fmax(1),0)  ! vz u_perp
      CALL skewflat(rtmp ,nx,ny,nz,wupa, wups, wupf,wupg,wupw,s2,s3,s4,s5,s6)       ! w up = vz u_perp

      ! Print out high order statistics  data:
      IF ( myrank.EQ.0 ) THEN
        fnout = trim(odir) // '/' // 'kavg.txt'
        OPEN(1,file=trim(fnout),position='append')
        WRITE(1,'(I,1x,20(g12.4,1x))',advance='no') indtime,as11,as12,as13,as22,as23,as33,va1,va2,va3,wa1,wa2,wa3,adiss,aenst,aheli,alamb,upa,dupa,wupa,betaa
        CLOSE(1)
        fnout = trim(odir) // '/' // 'skew.txt'
        OPEN(1,file=trim(fnout),position='append')
        WRITE(1,'(I,1x,20(g12.4,1x))',advance='no') indtime,ss11,ss12,ss13,ss22,ss23,ss33,vs1,vs2,vs3,ws1,ws2,ws3,sdiss,senst,sheli,slamb,ups,dups,wups,betas
        CLOSE(1)
        fnout = trim(odir) // '/' // 'flat.txt'
        OPEN(1,file=trim(fnout),position='append')
        WRITE(1,'(I,1x,20(g12.4,1x))',advance='no') indtime,sf11,sf12,sf13,sf22,sf23,sf33,vf1,vf2,vf3,wf1,wf2,wf3,fdiss,fenst,fheli,flamb,upf,dupf,wupf,betaf
        CLOSE(1)
        fnout = trim(odir) // '/' // 'glop.txt'
        OPEN(1,file=trim(fnout),position='append')
        WRITE(1,'(I,1x,20(g12.4,1x))',advance='no') indtime,sg11,sg12,sg13,sg22,sg23,sg33,vg1,vg2,vg3,wg1,wg2,wg3,gdiss,genst,gheli,glamb,upg,dupg,wupg,betag
        CLOSE(1)
        fnout = trim(odir) // '/' // 'whoa.txt'
        OPEN(1,file=trim(fnout),position='append')
        WRITE(1,'(I,1x,20(g12.4,1x))',advance='no') indtime,sw11,sw12,sw13,sw22,sw23,sw33,vw1,vw2,vw3,ww1,ww2,ww3,wdiss,wenst,wheli,wlamb,upw,dupw,wupw,betaw
        CLOSE(1)
      ENDIF


      If ( ibits(kin,0,1).EQ.1 ) THEN
      ! Compute, write, 1d enstrophy density pdf:
      fnout = trim(odir) // '/' // 'enstpdf_log.' // ext // '.txt'
      CALL dopdfr(S22,nx,ny,nz,fnout,nbins(1),0,fmin(1),fmax(1),1) 
      fnout = trim(odir) // '/' // 'enstpdf.' // ext // '.txt'
      CALL dopdfr(S22,nx,ny,nz,fnout,nbins(1),0,fmin(1),fmax(1),0) 
      ! Compute, write, 1d relative helicity pdf:
      fmin(1) = -1.0; fmax(1) = 1.0
      fnout = trim(odir) // '/' // 'helpdf.' // ext // '.txt'
      CALL dopdfr(S13,nx,ny,nz,fnout ,nbins(1),1,fmin(1),fmax(1),0) 
      ENDIF

      ! Compute joint PDF for various quantities 
      If ( ibits(kin,1,1).EQ.1 ) THEN
      fnout = trim(odir) // '/' // 'jpdf_enst_diss_log00.' // ext // '.txt'
      CALL dojpdfr(S22,'enst',S11,'diss',nx,ny,nz,fnout,nbins,[0,0],fmin,fmax,[0,0])
      fnout = trim(odir) // '/' // 'jpdf_enst_diss_log11.' // ext // '.txt'
      CALL dojpdfr(S22,'enst',S11,'diss',nx,ny,nz,fnout,nbins,[0,0],fmin,fmax,[1,1])
      ! Compute joint PDF for rel. helicity and diss :
      fnout = trim(odir) // '/' // 'jpdf_rhel_diss_log00.' // ext // '.txt'
      fmin(1) = -1.0; fmax(1) = 1.0
      CALL dojpdfr(S13,'rhel',S11,'diss',nx,ny,nz,fnout,nbins,[1,0],fmin,fmax,[0,0])
      fnout = trim(odir) // '/' // 'jpdf_hel_diss_log01.' // ext // '.txt'
      fmin(1) = -1.0; fmax(1) = 1.0
      CALL dojpdfr(S13,'rhel',S11,'diss',nx,ny,nz,fnout,nbins,[1,0],fmin,fmax,[0,1])
      ! Compute joint PDF for rel. helicity and enstroph. density:
      fnout = trim(odir) // '/' // 'jpdf_hel_enst_log01.' // ext // '.txt'
      fmin(1) = -1.0; fmax(1) = 1.0
      CALL dojpdfr(S13,'rhel',S22,'enst',nx,ny,nz,fnout,nbins,[1,0],fmin,fmax,[0,1])
      fnout = trim(odir) // '/' // 'jpdf_hel_enst_log01.' // ext // '.txt'
      fmin(1) = -1.0; fmax(1) = 1.0
      CALL dojpdfr(S13,'rhel',S22,'enst',nx,ny,nz,fnout,nbins,[1,0],fmin,fmax,[0,1])
      ! Compute joint PDF for rel. helicity and lambda:
      fnout = trim(odir) // '/' // 'jpdf_hel_lamb_log00.' // ext // '.txt'
      fmin(1) = -1.0; fmax(1) = 1.0
      CALL dojpdfr(S13,'rhel',lambda,'lamb',nx,ny,nz,fnout,nbins,[1,0],fmin,fmax,[0,0])
      fnout = trim(odir) // '/' // 'jpdf_hel_lamb_log01.' // ext // '.txt'
      fmin(1) = -1.0; fmax(1) = 1.0
      CALL dojpdfr(S13,'rhel',lambda,'lamb',nx,ny,nz,fnout,nbins,[1,0],fmin,fmax,[0,1])
      fnout = trim(odir) // '/' // 'jpdf_lamb_enst_log00.' // ext // '.txt'
      CALL dojpdfr(lambda,'lambda',S22,'enst',nx,ny,nz,fnout,nbins,[0,0],fmin,fmax,[0,0])
      fnout = trim(odir) // '/' // 'jpdf_lamb_enst_log11.' // ext // '.txt'
      CALL dojpdfr(lambda,'lambda',S22,'enst',nx,ny,nz,fnout,nbins,[0,0],fmin,fmax,[1,1])
      fnout = trim(odir) // '/' // 'jpdf_lamb_beta_log00.' // ext // '.txt'
      CALL dojpdfr(lambda,'lambda',beta,'beta',nx,ny,nz,fnout,nbins,[0,0],fmin,fmax,[0,0])
      fnout = trim(odir) // '/' // 'jpdf_lamb_beta_log11.' // ext // '.txt'
      CALL dojpdfr(lambda,'lambda',beta,'beta',nx,ny,nz,fnout,nbins,[0,0],fmin,fmax,[1,1])
      fnout = trim(odir) // '/' // 'jpdf_rhel_beta_log01.' // ext // '.txt'
      fmin(1) = -1.0; fmax(1) = 1.0
      CALL dojpdfr(S13,'rhel',beta,'beta',nx,ny,nz,fnout,nbins,[1,0],fmin,fmax,[0,1])
      fnout = trim(odir) // '/' // 'jpdf_enst_beta_log11.' // ext // '.txt'
      CALL dojpdfr(S22,'enst',beta,'beta',nx,ny,nz,fnout,nbins,[0,0],fmin,fmax,[1,1])
      fnout = trim(odir) // '/' // 'jpdf_diss_beta_log11.' // ext // '.txt'
      CALL dojpdfr(S11,'diss',beta,'beta',nx,ny,nz,fnout,nbins,[0,0],fmin,fmax,[1,1])
      ENDIF

      nin = nx*ny*(kend-ksta+1)

#if defined(SCALAR_)
! PDF for potl temperature:
! NOTE: S12 and S23 overwritten here
      ctmp = th;
      CALL fftp3d_complex_to_real(plancr,ctmp,rtmp1,MPI_COMM_WORLD)
      If ( ibits(kin,0,1).EQ.1 ) THEN
      fnout = trim(odir) // '/' // 'thpdf.' // ext // '.txt'
      CALL dopdfr(rtmp1,nx,ny,nz,fnout,nbins(1),0,fmin(1),fmax(1),0) 
      ENDIF
      CALL skewflat(rtmp1,nx,ny,nz,tha,ths,thf,thg,thw,s2,s3,s4,s5,s6)       ! theta

      ! Compute gradient Richardson no. and its pdf:
      CALL derivk3(vx, ctmp, 3)
      CALL fftp3d_complex_to_real(plancr,ctmp,S12,MPI_COMM_WORLD)
      CALL derivk3(vy, ctmp, 3)
      CALL fftp3d_complex_to_real(plancr,ctmp,S23,MPI_COMM_WORLD)
      CALL derivk3(th, ctmp, 1)
      CALL fftp3d_complex_to_real(plancr,ctmp,rtmp1,MPI_COMM_WORLD)
      If ( ibits(kin,0,1).EQ.1 ) THEN
      fnout = trim(odir) // '/' // 'dthdxpdf.' // ext // '.txt'
      CALL dopdfr(rtmp1,nx,ny,nz,fnout,nbins(1),0,fmin(1),fmax(1),0) 
      ENDIF
      CALL skewflat(rtmp1,nx,ny,nz,dthdxa,dthdxs,dthdxf,dthdxg,dthdxw,s2,s3,s4,s5,s6)       ! dtheta/dx
      CALL derivk3(th, ctmp, 2)
      CALL fftp3d_complex_to_real(plancr,ctmp,rtmp1,MPI_COMM_WORLD)
      If ( ibits(kin,0,1).EQ.1 ) THEN
      fnout = trim(odir) // '/' // 'dthdypdf.' // ext // '.txt'
      CALL dopdfr(rtmp1,nx,ny,nz,fnout,nbins(1),0,fmin(1),fmax(1),0) 
      ENDIF
      CALL skewflat(rtmp1,nx,ny,nz,dthdya,dthdys,dthdyf,dthdyg,dthdyw,s2,s3,s4,s5,s6)       ! dtheta/dy
      CALL derivk3(th, ctmp, 3)
      CALL fftp3d_complex_to_real(plancr,ctmp,rtmp1,MPI_COMM_WORLD)
      If ( ibits(kin,0,1).EQ.1 ) THEN
      fnout = trim(odir) // '/' // 'dthdzpdf.' // ext // '.txt'
      CALL dopdfr(rtmp1,nx,ny,nz,fnout,nbins(1),0,fmin(1),fmax(1),0) 
      ENDIF
      CALL skewflat(rtmp1,nx,ny,nz,dthdza,dthdzs,dthdzf,dthdzg,dthdzw,s2,s3,s4,s5,s6)       ! dtheta/dz
      IF ( iobin.GT.0 ) THEN
        rtmp2 = rtmp1
        IF ( oswap.NE.0 ) THEN
          CALL rarray_byte_swap(rtmp2, nx*ny*(kend-ksta+1))
        ENDIF
        CALL io_write(1,odir,'dthdz',ext,planio,rtmp2)
      ENDIF
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
        DO j = 1,ny
          DO i = 1,nx
            rtmp(i,j,k) = S12(i,j,k)**2 + S23(i,j,k)**2
          ENDDO
        ENDDO
      ENDDO
      IF ( iobin.GT.0 ) THEN
        IF ( oswap.NE.0 ) THEN
          rtmp2 = rtmp
          CALL rarray_byte_swap(rtmp2, nx*ny*(kend-ksta+1))
        ENDIF
        CALL io_write(1,odir,'vshw',ext,planio,rtmp2)  ! output (du_x/dz)^2 + (du_y/dz)^2
      ENDIF

      alpha = bvfreq/sqrt(4.0*omega(3)**2 + bvfreq**2)
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
        DO j = 1,ny
          DO i = 1,nx
            den = S12(i,j,k)**2 + S23(i,j,k)**2
            rig(i,j,k) = bvfreq*(bvfreq-rtmp1(i,j,k)) / den
            rign(i,j,k) = bvfreq*alpha*(bvfreq*alpha-rtmp1(i,j,k)) / den
          ENDDO
        ENDDO
      ENDDO
      IF ( iobin.GT.0 ) THEN
        rtmp2 = rig
        IF ( oswap.NE.0 ) THEN
          CALL rarray_byte_swap(rtmp2, nx*ny*(kend-ksta+1))
        ENDIF
        CALL io_write(1,odir,'rig',ext,planio,rtmp2)  ! output Ri_g
      ENDIF
      fnout = trim(odir) // '/' // 'riipdf.' // ext // '.txt'
      fmin(1) = -rmi;
      fmax(1) =  rmi;
      If ( ibits(kin,0,1).EQ.1 ) THEN
      CALL dopdfr(rig,nx,ny,nz,'riipdf.'//ext//'.txt',nbins(1),1,fmin(1),fmax(1),0) 
      ENDIF
      CALL skewflat(rig,nx,ny,nz,ria,ris,rif,rid,riw,s2,s3,s4,s5,s6)       ! Ri
      ! Compute joint PDF for Richardson no. and diss:
      If ( ibits(kin,1,1).EQ.1 ) THEN
      fmin  = [-rmi,0.0]; fmax  = [rmi,0.0]
      fnout = trim(odir) // '/' // 'jpdf_rig_beta_log00.' // ext // '.txt'
      CALL dojpdfr(rig,'Rig' ,beta,'beta',nx,ny,nz,fnout,nbins,[1,0],fmin,fmax,[0,0])
      fnout = trim(odir) // '/' // 'jpdf_rig_beta_log01.' // ext // '.txt'
      CALL dojpdfr(rig,'Rig' ,beta,'beta',nx,ny,nz,fnout,nbins,[1,0],fmin,fmax,[0,1])
      fnout = trim(odir) // '/' // 'jpdf_rign_beta_log00.' // ext // '.txt'
      CALL dojpdfr(rign,'Rign' ,beta,'beta',nx,ny,nz,fnout,nbins,[1,0],fmin,fmax,[0,0])
      fnout = trim(odir) // '/' // 'jpdf_rign_beta_log01.' // ext // '.txt'
      CALL dojpdfr(rign,'Rign' ,beta,'beta',nx,ny,nz,fnout,nbins,[1,0],fmin,fmax,[0,1])
      fnout = trim(odir) // '/' // 'jpdf_rig_diss_log00.' // ext // '.txt'
      fmin  = [-rmi,0.0]; fmax  = [rmi,0.0]
      CALL dojpdfr(rig,'Rig' ,S11 ,'diss',nx,ny,nz,fnout,nbins,[1,0],fmin,fmax,[0,0])
      fnout = trim(odir) // '/' // 'jpdf_rig_diss_log01.' // ext // '.txt'
      fmin  = [-rmi,0.0]; fmax  = [rmi,0.0]
      CALL dojpdfr(rig,'Rig' ,S11 ,'diss',nx,ny,nz,fnout,nbins,[1,0],fmin,fmax,[0,1])
      fnout = trim(odir) // '/' // 'jpdf_rig_enst_log00.' // ext // '.txt'
      fmin  = [-rmi,0.0]; fmax  = [rmi,0.0]
      CALL dojpdfr(rig,'Rig' ,S22 ,'enst',nx,ny,nz,fnout,nbins,[1,0],fmin,fmax,[0,0])
      fnout = trim(odir) // '/' // 'jpdf_rig_enst_log01.' // ext // '.txt'
      fmin  = [-rmi,0.0]; fmax  = [rmi,0.0]
      CALL dojpdfr(rig,'Rig' ,S22 ,'enst',nx,ny,nz,fnout,nbins,[1,0],fmin,fmax,[0,1])
      fnout = trim(odir) // '/' // 'jpdf_rhel_rig_log00.' // ext // '.txt'
      fmin  = [-1.0,-rmi]; fmax  = [1.0,rmi]
      CALL dojpdfr(S13,'rhel',rig,'Rig'  ,nx,ny,nz,fnout,nbins,[1,1],fmin,fmax,[0,0])
      fnout = trim(odir) // '/' // 'jpdf_rig_lamb_log00.' // ext // '.txt'
      fmin  = [-rmi,0.0]; fmax  = [rmi,0.0]
      CALL dojpdfr(rig,'Rgi',lambda,'lambda',nx,ny,nz,fnout,nbins,[1,0],fmin,fmax,[0,0])
      fnout = trim(odir) // '/' // 'jpdf_rig_lamb_log01.' // ext // '.txt'
      fmin  = [-rmi,0.0]; fmax  = [rmi,0.0]
      CALL dojpdfr(rig,'Rig',lambda,'lambda',nx,ny,nz,fnout,nbins,[1,0],fmin,fmax,[0,1])

      CALL derivk3(th, ctmp, 3)
      CALL fftp3d_complex_to_real(plancr,ctmp,rtmp1,MPI_COMM_WORLD)
      fnout = trim(odir) // '/' // 'jpdf_rig_dtdz_log00.' // ext // '.txt'
      fmin  = [-rmi,0.0]; fmax  = [rmi,0.0]
      CALL dojpdfr(rig,'Rig',rtmp1,'dtdz',nx,ny,nz,fnout,nbins,[1,0],fmin,fmax,[0,0])
      ENDIF

#if defined(SCALAR_)
      !! Compute stats for PV quantities:
      !! ...First, omega .Grad theta:
!!    xnormn = 1.0_GP/ real(nx*ny*nz,kind=GP)**2
if (myrank.eq.0) write(*,*)'DoKPDF: compute omega_x_Grad_x theta...'
      CALL derivk3(th, ctmp, 1)  ! x-deriv
      CALL fftp3d_complex_to_real(plancr,ctmp,rtmp1,MPI_COMM_WORLD)
      CALL rotor3(vy,vz,ctmp,1)
      CALL fftp3d_complex_to_real(plancr,ctmp,rtmp2,MPI_COMM_WORLD)
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
        DO j = 1,ny
          DO i = 1,nx
            rtmp(i,j,k) = rtmp1(i,j,k)*rtmp2(i,j,k) ! += omega_x Grad_x theta
          ENDDO
        ENDDO
      ENDDO

if (myrank.eq.0) write(*,*)'DoKPDF: compute omega_y_Grad_y theta...'
      CALL derivk3(th, ctmp, 2)  ! y-deriv
      CALL fftp3d_complex_to_real(plancr,ctmp,rtmp1,MPI_COMM_WORLD)
      CALL rotor3(vy,vz,ctmp,2)
      CALL fftp3d_complex_to_real(plancr,ctmp,rtmp2,MPI_COMM_WORLD)
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
        DO j = 1,ny
          DO i = 1,nx
            rtmp(i,j,k) = rtmp(i,j,k) + rtmp1(i,j,k)*rtmp2(i,j,k) ! += omega_y Grad_y theta
          ENDDO
        ENDDO
      ENDDO

if (myrank.eq.0) write(*,*)'DoKPDF: omega.Grad th = ', rtmp(1:10,1,ksta)
if (myrank.eq.0) write(*,*)'DoKPDF: compute omega_z_Grad_z theta...'
      CALL derivk3(th, ctmp, 3)  ! z-deriv
      CALL fftp3d_complex_to_real(plancr,ctmp,rtmp1,MPI_COMM_WORLD)
      CALL rotor3(vy,vz,ctmp,3)
      CALL fftp3d_complex_to_real(plancr,ctmp,rtmp2,MPI_COMM_WORLD)
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
        DO j = 1,ny
          DO i = 1,nx
            rtmp(i,j,k) = rtmp(i,j,k) + rtmp1(i,j,k)*rtmp2(i,j,k) ! += omega_z Grad_z theta
            !pv = fdtheta/dz -N \omega_z +  omega_z Grad_z theta:
            pv  (i,j,k) = 2.0*omegaz*rtmp1(i,j,k) - bvfreq*rtmp2(i,j,k) + rtmp(i,j,k) 
          ENDDO
        ENDDO
      ENDDO

      ! Here, rtmp contains omega. Graad theta (wgt)
      !       pc   contains full PV:
      If ( ibits(kin,0,1).EQ.1 ) THEN
if (myrank.eq.0) write(*,*)'DoKPDF: write PDFs of PV and omega Grad theta...'
      fnout = trim(odir) // '/' // 'pvpdf.' // ext // '.txt'
if (myrank.eq.0) write(*,*)'DoKPDF: pv=',pv(1:10,1,ksta);
      CALL dopdfr(pv  ,nx,ny,nz,fnout,nbins(1),0,fmin(1),fmax(1),0) 
      fnout = trim(odir) // '/' // 'wgradtpdf.' // ext // '.txt'
      CALL dopdfr(rtmp,nx,ny,nz,fnout,nbins(1),0,fmin(1),fmax(1),0) 
      ENDIF
      CALL skewflat(pv  ,nx,ny,nz,pva,pvs ,pvf ,pvg,pvw,s2,s3,s4,s5,s6)       ! PV
      CALL skewflat(rtmp,nx,ny,nz,wgta,wgts,wgtf,wgtg,wgtw,s2,s3,s4,s5,s6)       ! omega.Grad theta
      If ( ibits(kin,1,1).EQ.1 ) THEN
      fnout = trim(odir) // '/' // 'jpdf_wgt_diss_log00.' // ext // '.txt'
      CALL dojpdfr(rtmp,'wgt' ,S11 ,'diss',nx,ny,nz,fnout,nbins,[0,0],fmin,fmax,[0,0])
      fnout = trim(odir) // '/' // 'jpdf_wgt_diss_log01.' // ext // '.txt'
      CALL dojpdfr(rtmp,'wgt' ,S11 ,'diss',nx,ny,nz,fnout,nbins,[0,0],fmin,fmax,[0,1])
      fnout = trim(odir) // '/' // 'jpdf_pv_diss_log01.' // ext // '.txt'
      CALL dojpdfr(pv,'PV' ,S11,'diss',nx,ny,nz,fnout,nbins,[0,0],fmin,fmax,[0,1])
      fnout = trim(odir) // '/' // 'jpdf_pv_rhel_log00.' // ext // '.txt'
      fmin(2) = -1.0; fmax(2) = 1.0
      CALL dojpdfr(pv,'PV',S13,'rhel',nx,ny,nz,fnout,nbins,[0,1],fmin,fmax,[0,0])
      fnout = trim(odir) // '/' // 'jpdf_pv_rig_log00.' // ext // '.txt'
      CALL dojpdfr(pv,'PV',rig,'Rig',nx,ny,nz,fnout,nbins,[0,0],fmin,fmax,[0,0])
      fnout = trim(odir) // '/' // 'jpdf_pv_beta_log00.' // ext // '.txt'
      CALL dojpdfr(pv,'PV',beta,'beta',nx,ny,nz,fnout,nbins,[0,0],fmin,fmax,[0,0])
      fnout = trim(odir) // '/' // 'jpdf_pv_beta_log01.' // ext // '.txt'
      CALL dojpdfr(pv,'PV',beta,'beta',nx,ny,nz,fnout,nbins,[0,0],fmin,fmax,[0,1])
      ENDIF

      ! Ratio of NL PV to total PV
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
        DO j = 1,ny
          DO i = 1,nx
            rtmp(i,j,k) = abs(rtmp(i,j,k)/pv(i,j,k))
          ENDDO
        ENDDO
      ENDDO
      fnout = trim(odir) // '/' // 'wgt2pvpdf.' // ext // '.txt'
      CALL dopdfr(rtmp  ,nx,ny,nz,fnout,nbins(1),0,fmin(1),fmax(1),0) 
      CALL skewflat(rtmp,nx,ny,nz,wgt2pva,wgt2pvs,wgt2pvf,wgt2pvg,wgt2pvw,s2,s3,s4,s5,s6)       ! omega.Grad theta / PV
#endif

      !! Compute statistics for w\theta:
      ctmp = th;
      CALL fftp3d_complex_to_real(plancr,ctmp,rtmp1,MPI_COMM_WORLD)
      ctmp = vz;
      CALL fftp3d_complex_to_real(plancr,ctmp,S23  ,MPI_COMM_WORLD)
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
        DO j = 1,ny
          DO i = 1,nx
            rtmp1(i,j,k) = S23(i,j,k)*rtmp1(i,j,k)
          ENDDO
        ENDDO
      ENDDO
      IF ( iobin.GT.0 ) THEN
        rtmp2 = rtmp
        IF ( oswap.NE.0 ) THEN
          CALL rarray_byte_swap(rtmp2, nx*ny*(kend-ksta+1))
        ENDIF
        CALL io_write(1,odir,'wth',ext,planio,rtmp2)  ! output w\theta
      ENDIF
      If ( ibits(kin,0,1).EQ.1 ) THEN
      fnout = trim(odir) // '/' // 'wtpdf.' // ext // '.txt'
      CALL dopdfr(rtmp1,nx,ny,nz,fnout,nbins(1),0,fmin(1),fmax(1),0) 
      ENDIF
      CALL skewflat(rtmp1,nx,ny,nz,wta,wts,wtf,wtg,wtw,s2,s3,s4,s5,s6)       ! w theta
      If ( ibits(kin,1,1).EQ.1 ) THEN
      fnout = trim(odir) // '/' // 'jpdf_rig_wtheta_log00.' // ext // '.txt'
      fmin  = [-rmi,0.0]; fmax  = [rmi,0.0]
      CALL dojpdfr(rig,'Rig' ,rtmp1 ,'wtheta',nx,ny,nz,fnout,nbins,[1,0],fmin,fmax,[0,0])
      fnout = trim(odir) // '/' // 'jpdf_wtheta_beta_log00.' // ext // '.txt'
      CALL dojpdfr(rtmp1,'wtheta',beta,'beta',nx,ny,nz,fnout,nbins,[0,0],fmin,fmax,[0,0])
      fnout = trim(odir) // '/' // 'jpdf_wtheta_beta_log01.' // ext // '.txt'
      CALL dojpdfr(rtmp1,'wtheta',beta,'beta',nx,ny,nz,fnout,nbins,[0,0],fmin,fmax,[0,1])
      ENDIF

      ! Need xnormn bec. advect3 multiplies by -1/N^6:
      xnormn = -1.0  !!*real(nx*ny*nz,kind=GP)**2
      CALL advect3(vx,vy,vz,vx,ctmp)    ! v.Grad vx 
      ctmp = ctmp *xnorm
      CALL fftp3d_complex_to_real(plancr,ctmp,S11  ,MPI_COMM_WORLD)
      ctmp = vx
      CALL fftp3d_complex_to_real(plancr,ctmp,S23  ,MPI_COMM_WORLD)
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
        DO j = 1,ny
          DO i = 1,nx
            S11(i,j,k) = S11(i,j,k)*S23(i,j,k)  ! vx (v.Grad vx)
          ENDDO
        ENDDO
      ENDDO

      CALL advect3(vx,vy,vz,vy,ctmp)     !  v.Grad vy 
      ctmp = ctmp *xnorm
      CALL fftp3d_complex_to_real(plancr,ctmp,S22  ,MPI_COMM_WORLD)
      ctmp = vy
      CALL fftp3d_complex_to_real(plancr,ctmp,S23  ,MPI_COMM_WORLD)
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
        DO j = 1,ny
          DO i = 1,nx
            S11(i,j,k) = S11(i,j,k) + S22(i,j,k)*S23(i,j,k)  ! += vy (v.Grad vy)
          ENDDO
        ENDDO
      ENDDO
      
      CALL advect3(vx,vy,vz,vz,ctmp)    ! v.Grad vz 
      ctmp = ctmp *xnorm
      CALL fftp3d_complex_to_real(plancr,ctmp,S22  ,MPI_COMM_WORLD)
      ctmp = vy
      CALL fftp3d_complex_to_real(plancr,ctmp,S23  ,MPI_COMM_WORLD)
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
        DO j = 1,ny
          DO i = 1,nx
            S11(i,j,k) = S11(i,j,k) + S22(i,j,k)*S23(i,j,k)  ! += vz (v.Grad vz)
          ENDDO
        ENDDO
      ENDDO
      If ( ibits(kin,0,1).EQ.1 ) THEN
      fnout = trim(odir) // '/' // 'ktrpdf.' // ext // '.txt'
      CALL dopdfr(S11,nx,ny,nz,fnout,nbins(1),0,fmin(1),fmax(1),0) 
      ENDIF
      CALL skewflat(S11,nx,ny,nz,ktra,ktrs,ktrf,ktrg,ktrw,s2,s3,s4,s5,s6)       ! u.(u.Grad u)
      If ( ibits(kin,1,1).EQ.1 ) THEN
      fnout = trim(odir) // '/' // 'jpdf_ktr_wtheta_log00.' // ext // '.txt'
      CALL dojpdfr(S11,'ktrans' ,rtmp1 ,'wtheta',nx,ny,nz,fnout,nbins,[0,0],fmin,fmax,[0,0])
      fnout = trim(odir) // '/' // 'jpdf_ktr_beta_log00.' // ext // '.txt'
      CALL dojpdfr(S11,'ktrans' ,beta,'beta',nx,ny,nz,fnout,nbins,[0,0],fmin,fmax,[0,0])
      fnout = trim(odir) // '/' // 'jpdf_ktr_beta_log01.' // ext // '.txt'
      CALL dojpdfr(S11,'ktrans' ,beta,'beta',nx,ny,nz,fnout,nbins,[0,0],fmin,fmax,[0,1])
      fnout = trim(odir) // '/' // 'jpdf_ktr_pv_log00.' // ext // '.txt'
      CALL dojpdfr(S11,'ktrans' ,pv,'PV',nx,ny,nz,fnout,nbins,[0,0],fmin,fmax,[0,0])
      ENDIF

      CALL advect3(vx,vy,vz,th,ctmp)    ! v.Grad theta
      ctmp = ctmp *xnorm
      CALL fftp3d_complex_to_real(plancr,ctmp,S11  ,MPI_COMM_WORLD)
      ctmp = th
      CALL fftp3d_complex_to_real(plancr,ctmp,S22  ,MPI_COMM_WORLD)
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
        DO j = 1,ny
          DO i = 1,nx
            S11(i,j,k) = S11(i,j,k)*S22(i,j,k) ! += theta (v.Grad theta)
          ENDDO
        ENDDO
      ENDDO
      If ( ibits(kin,0,1).EQ.1 ) THEN
      fnout = trim(odir) // '/' // 'ptrpdf.' // ext // '.txt'
      CALL dopdfr(S11,nx,ny,nz,fnout,nbins(1),0,fmin(1),fmax(1),0) 
      ENDIF
      CALL skewflat(S11,nx,ny,nz,ptra,ptrs,ptrf,ptrg,ptrw,s2,s3,s4,s5,s6)       ! \theta.(u.Grad \theta)
      If ( ibits(kin,1,1).EQ.1 ) THEN
      fnout = trim(odir) // '/' // 'jpdf_ptr_wtheta_log00.' // ext // '.txt'
      CALL dojpdfr(S11,'ptrans' ,rtmp1 ,'wtheta',nx,ny,nz,fnout,nbins,[0,0],fmin,fmax,[0,0])
      fnout = trim(odir) // '/' // 'jpdf_ptr_beta_log00.' // ext // '.txt'
      CALL dojpdfr(S11,'ptrans' ,beta,'beta',nx,ny,nz,fnout,nbins,[0,0],fmin,fmax,[0,0])
      fnout = trim(odir) // '/' // 'jpdf_ptr_beta_log01.' // ext // '.txt'
      CALL dojpdfr(S11,'ptrans' ,beta,'beta',nx,ny,nz,fnout,nbins,[0,0],fmin,fmax,[0,1])
      fnout = trim(odir) // '/' // 'jpdf_ptr_PV_log00.' // ext // '.txt'
      CALL dojpdfr(S11,'ptrans' ,pv,'PV',nx,ny,nz,fnout,nbins,[0,0],fmin,fmax,[0,0])
      ENDIF

      ! Print out skewness and flatness SCALAR data:
      IF ( myrank.EQ.0 ) THEN
        fnout = trim(odir) // '/' // 'scavg.txt'
        OPEN(1,file=trim(fnout),position='append')
      ! Print out avg SCALAR data:
        WRITE(1,'(I,1x,11(g12.4,1x))',advance='no') indtime,tha,dthdxa,dthdya,dthdza,ria,wta,ktra,ptra,pva,wgta,wgt2pva
        CLOSE(1)
        fnout = trim(odir) // '/' // 'scskew.txt'
        OPEN(1,file=trim(fnout),position='append')
      ! Print out skewness and flatness SCALAR data:
        WRITE(1,'(I,1x,11(g12.4,1x))',advance='no') indtime,ths,dthdxs,dthdys,dthdzs,ris,wts,ktrs,ptrs,pvs,wgts,wgt2pvs
        CLOSE(1)
        fnout = trim(odir) // '/' // 'scflat.txt'
        OPEN(1,file=trim(fnout),position='append')
        WRITE(1,'(I,1x,11(g12.4,1x))',advance='no') indtime,thf,dthdxf,dthdyf,dthdzf,rif,wtf,ktrf,ptrf,pvf,wgtf,wgt2pvf
        CLOSE(1)
        fnout = trim(odir) // '/' // 'scglop.txt'
        OPEN(1,file=trim(fnout),position='append')
        WRITE(1,'(I,1x,12(g12.4,1x))',advance='no') indtime,thg,dthdxg,dthdyg,dthdzg,rid,wtg,ktrg,ptrg,pvg,wgtg,wgt2pvg
        CLOSE(1)
        fnout = trim(odir) // '/' // 'scwhoa.txt'
        OPEN(1,file=trim(fnout),position='append')
        WRITE(1,'(I,1x,12(g12.4,1x))',advance='no') indtime,thw,dthdxw,dthdyw,dthdzw,riw,wtw,ktrw,ptrw,pvw,wgtw,wgt2pvw
        CLOSE(1)
      ENDIF
#endif

      END SUBROUTINE DoKPDF
!
!
!
      SUBROUTINE skewflat(fx,nx,ny,nz,avg,skew,flat,glop,whoa,s2,s3,s4,s5,s6)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!
! Computes skewnewss and flatness of input random variable
! 
!
! Parameters
!     fx    : input real random variable. Must already be normalized!
!     nx,ny,
!     nz    : dimensions
!     avg   : average
!     skew  : skewness, valid only for MPI taks 0
!     flat  : flatness/kurtosis
!     glop  : Sum (x-<x>)^5/[Sum(x-<x>)^2]^5/2, 5th order moment ('glop')
!     whoa  : Sum (x-<x>)^6/[Sum(x-<x>)^2]^3, 6th order moment ('whoa')
!     s2-s6 : 2nd-6th order moments: sum( (x-<x>)^n ), where n=2-6.
!-----------------------------------------------------------------
      USE fprecision
      USE commtypes
      IMPLICIT NONE

      REAL(KIND=GP), INTENT (IN), DIMENSION(nx*ny*nz) :: fx
      REAL(KIND=GP), INTENT(OUT)                 :: skew,flat,glop,whoa
      REAL(KIND=GP), INTENT(OUT)                 :: avg,s2,s3,s4,s5,s6
      DOUBLE PRECISION                           :: davg,ds2,ds3,ds4,ds5,ds6
      DOUBLE PRECISION                           :: gs(5),s(5),xnorm
      INTEGER      , INTENT (IN)                 :: nx,ny,nz
      INTEGER                                    :: ierr,j

      INTEGER  nin

      nin = nx * ny * nz

      xnorm = 1.0_GP / dble(nx*ny*nz)
      ds2 = 0.0D0
!$omp parallel do default(shared) private(j) reduction(+:s2)
      DO j = 1, nin
        ds2 = ds2 + xnorm*dble(fx(j))
      ENDDO
!$omp end parallel do

      CALL MPI_ALLREDUCE(ds2, davg, 1, MPI_DOUBLE_PRECISION, &
                      MPI_SUM, MPI_COMM_WORLD,ierr)
      avg = davg

      ds2 = 0.0D0
!$omp parallel do default(shared) private(j) reduction(+:s2)
      DO j = 1, nin
        ds2 = ds2 + xnorm*(dble(fx(j))-davg)**2
      ENDDO
!$omp end parallel do

      ds3 = 0.0D0
!$omp parallel do default(shared) private(j) reduction(+:s3)
      DO j = 1, nin
        ds3 = ds3 + xnorm*(dble(fx(j))-davg)**3
      ENDDO
!$omp end parallel do

      ds4 = 0.0D0
!$omp parallel do default(shared) private(j) reduction(+:s4)
      DO j = 1, nin
        ds4 = ds4 + xnorm*(dble(fx(j))-davg)**4
      ENDDO
!$omp end parallel do

      ds5 = 0.0D0
!$omp parallel do default(shared) private(j) reduction(+:s5)
      DO j = 1, nin
        ds5 = ds5 + xnorm*(dble(fx(j))-davg)**5
      ENDDO
!$omp end parallel do

      ds6 = 0.0D0
!$omp parallel do default(shared) private(j) reduction(+:s6)
      DO j = 1, nin
        ds6 = ds6 + xnorm*(dble(fx(j))-davg)**6
      ENDDO
!$omp end parallel do

      s(1)=ds2; s(2)=ds3; s(3)=ds4; s(4) = ds5; s(5) = ds6
      CALL MPI_ALLREDUCE(s, gs, 5, MPI_DOUBLE_PRECISION, &
                      MPI_SUM, MPI_COMM_WORLD,ierr)
      if ( ierr.ne.MPI_SUCCESS ) then
        write(*,*)'skewflat: final allreduce failed'
        stop
      endif
      s2=gs(1); s3=gs(2); s4=gs(3); s5=gs(4); s6=gs(5)

      skew = real(s3 / (s2+1.0e-15)**(1.5),kind=GP)
      flat = real(s4 / (s2+1.0e-15)**(2.0),kind=GP)
      glop = real(s5 / (s2+1.0e-15)**(2.5),kind=GP)
      whoa = real(s6 / (s2+1.0e-15)**(3.0),kind=GP)

      END SUBROUTINE skewflat
!-----------------------------------------------------------------
!-----------------------------------------------------------------

!
!
      SUBROUTINE Randomize(fx,fy,fz,krmin,krmax,phase)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!
! Randomizes phases of vector field, (vx,vy,vz) between
! wavenumbers krmin,krmax
!
! Parameters
!     vx,
!     vy,
!     vz    : complex velocities
!     krmin : minimum wavenumber
!     krmax : maximum wavenumber
!     phase : phase 
!
      USE fprecision
      USE commtypes
      USE kes
      USE grid
      USE mpivars
      USE threads
      USE fft
      USE var
      USE fftplans
      USE gutils
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(nz,ny,ista:iend) :: fx,fy,fz
      REAL   (KIND=GP), INTENT   (IN)                           :: krmin,krmax,phase
      REAL   (KIND=GP)                                          :: krmin2,krmax2
      COMPLEX(KIND=GP)                                          :: cdump,jdump
      INTEGER                                                   :: i,j,k

      cdump = COS(phase)+im*SIN(phase)
      jdump = conjg(cdump)

      krmin2 = krmin**2
      krmax2 = krmax**2
      IF (ista.eq.1) THEN

      DO j = 2,ny/2+1
        IF ( kk2(1,j,1).gt.krmin2 .and. kk2(1,j,1).lt.krmax2 ) THEN
        fx     (1,j,1) = fx     (1,j,1)*cdump
        fx(1,ny-j+2,1) = fx(1,ny-j+2,1)*jdump
        fy     (1,j,1) = fy     (1,j,1)*cdump
        fy(1,ny-j+2,1) = fy(1,ny-j+2,1)*jdump
        fz     (1,j,1) = fz     (1,j,1)*cdump
        fz(1,ny-j+2,1) = fz(1,ny-j+2,1)*jdump
        ENDIF
      ENDDO

!$omp parallel do
      DO k = 2,nz/2+1
        IF ( kk2(k,1,1).gt.krmin2 .and. kk2(k,1,1).lt.krmax2 ) THEN
        fx     (k,1,1) = fx     (k,1,1)*cdump
        fx(nz-k+2,1,1) = fx(nz-k+2,1,1)*jdump
        fy    (kz,1,1) = fy     (k,1,1)*cdump
        fy(nz-k+2,1,1) = fy(nz-k+2,1,1)*jdump
        fz     (k,1,1) = fz     (k,1,1)*cdump
        fz(nz-k+2,1,1) = fz(nz-k+2,1,1)*jdump
        ENDIF
      END DO
!$omp end parallel do
      
!$omp parallel do private (k)
      DO j = 2,ny
        DO k = 2,nz/2+1
          IF ( kk2(k,j,1).gt.krmin2 .and. kk2(k,j,1).lt.krmax2 ) THEN
          fx          (k,j,1) = fx          (k,j,1)*cdump
          fx(nz-k+2,ny-j+2,1) = fx(nz-k+2,ny-j+2,1)*jdump
          fy          (k,j,1) = fy          (k,j,1)*cdump
          fy(nz-k+2,ny-j+2,1) = fy(nz-k+2,ny-j+2,1)*jdump
          fz          (k,j,1) = fz          (k,j,1)*cdump
          fz(nz-k+2,ny-j+2,1) = fz(nz-k+2,ny-j+2,1)*jdump
          ENDIF
        ENDDO
      ENDDO
!$omp end parallel do


!$omp parallel do if (iend-2.ge.nth) private (j,k)
      DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k)
        DO j = 1,ny
          DO k = 1,nz
            IF ( kk2(k,j,i).gt.krmin2 .and. kk2(k,j,i).lt.krmax2 ) THEN
            fx(k,j,i) = fx(k,j,i)*cdump
            fy(k,j,i) = fy(k,j,i)*cdump
            fz(k,j,i) = fz(k,j,i)*cdump
            ENDIF
          ENDDO
        ENDDO
      ENDDO

      ELSE
    
      DO  i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
           DO k = 1,nz
             IF ( kk2(k,j,i).gt.krmin2 .and. kk2(k,j,i).lt.krmax2 ) THEN
             fx(k,j,i) = fx(k,j,i)*cdump
             fy(k,j,i) = fy(k,j,i)*cdump
             fz(k,j,i) = fz(k,j,i)*cdump
            ENDIF
           ENDDO
         ENDDO
       ENDDO
       ENDIF
  
      END SUBROUTINE Randomize
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!
!
      SUBROUTINE Shebalin(asheb,kin,vx,vy,vz,th)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!
! Compute Shebalin angle measure of anisotrop.
!
! Parameters
!     asheb : returned Shebalin angle measurement
!     kin   : 
!       == 0: do total energy; 
!       == 1: do kinetic energy only; 
!       == 2: do potential energy only; (done only if th is specified)
!     vx,
!     vy,
!     vz    : complex velocities
!     th    : optional complex potl temperature
!
      USE fprecision
      USE commtypes
      USE kes
      USE grid
      USE mpivars
      USE threads
      USE fft
      USE var
      USE fftplans
      USE gutils
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(nz,ny,ista:iend)           :: vx,vy,vz
      COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(nz,ny,ista:iend), OPTIONAL :: th
      REAL   (KIND=GP), INTENT  (OUT)                           :: asheb
      REAL   (KIND=GP)                                          :: lasheb
      REAL   (KIND=GP)                                          :: tmp
      DOUBLE PRECISION                                          :: kprp,kpar,tmq,lprp,lpar,gprp,gpar,la(2),ga(2)
      INTEGER         , INTENT (IN)                             :: kin
      INTEGER                                                   :: i,j,k

      tmq  = 1.0D0/real(nx*ny*nz,kind=GP)**2

      lpar = 0.0D0
      lprp = 0.0D0
      IF (kin.eq.0 .OR. kin.eq.1) THEN ! total or kinetic energy
        IF (ista.eq.1) THEN
!$omp parallel do private (k,kprp,kpar,tmp) reduction(+:lprp,lpar)
           DO j = 1,ny
              DO k = 1,nz
                kprp = kx(1)**2 + ky(j)**2
                kpar = kz(k)**2
                tmp  = ( abs(vx(k,j,1))**2+abs(vy(k,j,1))**2+          &
                         abs(vz(k,j,1))**2 )*tmq
                lprp  = lprp + tmp*kprp
                lpar  = lpar + tmp*kpar
              END DO
           END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k,kprp,kpar,tmp) reduction(+:lprp,lpar)
           DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k,kprp,kpar,tmp) reduction(+:lprp,lpar)
              DO j = 1,ny
                 DO k = 1,nz
                    kprp = kx(i)**2 + ky(j)**2
                    kpar = kz(k)**2
                    tmp  = 2.0*( abs(vx(k,j,i))**2+abs(vy(k,j,i))**2+       &
                                 abs(vz(k,j,i))**2 )*tmq
                    lprp  = lprp + tmp*kprp
                    lpar  = lpar + tmp*kpar
                 END DO
              END DO
           END DO
        ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k,kprp,kpar,tmp) reduction(+:lprp,lpar)
           DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,kprp,kpar,tmp) reduction(+:lprp,lpar)
              DO j = 1,ny
                 DO k = 1,nz
                    kprp = kx(i)**2 + ky(j)**2
                    kpar = kz(k)**2
                    tmp  = 2.0*( abs(vx(k,j,i))**2+abs(vy(k,j,i))**2+       &
                                 abs(vz(k,j,i))**2 )*tmq
                    lprp  = lprp + tmp*kprp
                    lpar  = lpar + tmp*kpar
                 END DO
              END DO
           END DO
        ENDIF
      ENDIF

      IF (kin.EQ.1) THEN ! kinetic energy only
        la(1) = lprp
        la(2) = lpar
        CALL MPI_ALLREDUCE(la,ga,2,MPI_DOUBLE_PRECISION,MPI_SUM,   &
                        MPI_COMM_WORLD,ierr)
        asheb  = real(ga(1) / (ga(2) + tiny(1.0_GP)),kind=GP)
        RETURN
      ENDIF

      IF (.NOT.PRESENT(th) ) THEN
        WRITE(*,*) 'Shebalin: potential temperature not provided'
        STOP
      ENDIF

      IF (kin.EQ.2) THEN ! if potential energy only
        lprp = 0.0_GP
        lpar = 0.0_GP
      ENDIF

      IF (ista.eq.1) THEN ! add in the potential component
!$omp parallel do private (k,kprp,kpar,tmp) reduction(+:lprp,lpar)
         DO j = 1,ny
            DO k = 1,nz
              kprp = kx(1)**2 + ky(j)**2
              kpar = kz(k)**2
              tmp  = ( abs(th(k,j,1))**2 )*tmq
              lprp  = lprp + tmp*kprp
              lpar  = lpar + tmp*kpar
            END DO
         END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k,kprp,kpar,tmp) reduction(+:lprp,lpar)
         DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k,kprp,kpar,tmp) reduction(+:lprp,lpar)
            DO j = 1,ny
               DO k = 1,nz
                  kprp = kx(i)**2 + ky(j)**2
                  kpar = kz(k)**2
                  tmp  = 2.0*( abs(th(k,j,i))**2 )*tmq
                  lprp  = lprp + tmp*kprp
                  lpar  = lpar + tmp*kpar
               END DO
            END DO
         END DO
      ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k,kprp,kpar,tmp) reduction(+:lprp,lpar)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,kprp,kpar,tmp) reduction(+:lprp,lpar)
            DO j = 1,ny
               DO k = 1,nz
                  kprp = kx(i)**2 + ky(j)**2
                  kpar = kz(k)**2
                  tmp  = 2.0*( abs(vx(k,j,i))**2+abs(vy(k,j,i))**2+       &
                           abs(vz(k,j,i))**2 )*tmq
                  lprp  = lprp + tmp*kprp
                  lpar  = lpar + tmp*kpar
               END DO
            END DO
         END DO
      ENDIF
!
! Computes the reduction between nodes, store in return variable
      la(1) = lprp
      la(2) = lpar
      CALL MPI_ALLREDUCE(la,ga,2,MPI_DOUBLE_PRECISION,MPI_SUM,   &
                        MPI_COMM_WORLD,ierr)
      asheb  = real(ga(1) / (ga(2) + tiny(1.0_GP)),kind=GP)

      END SUBROUTINE Shebalin
!-----------------------------------------------------------------
!-----------------------------------------------------------------

      SUBROUTINE Strain(vx,vy,vz,ir,jc,ktmin,ktmax,inorm,sij,ctmp)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!
! Computes the complex strain rate component 
!
! Parameters
!     vi   : input velocities
!     sij  : returned complex component of strain rate tensor 
!     ir,jc: the row and col of sij
!     ktmin: truncaton min wavenumber for spherical truncation
!     ktmax: truncaton max wavenumber for spherical truncation
!     inorm : normalize (1), or not (0)
!
      USE fprecision
      USE commtypes
      USE kes
      USE grid
      USE mpivars
      USE ali
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT   (IN), DIMENSION(nz,ny,ista:iend) :: vx,vy,vz
      COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(nz,nz,ista:iend) :: ctmp,sij
      REAL   (KIND=GP), INTENT   (IN)                             :: ktmin,ktmax
      REAL   (KIND=GP)                                            :: ktmin2,ktmax2,tmp
!
      INTEGER         , INTENT   (IN)                           :: inorm,ir,jc
      INTEGER                                                   :: i,j,k

      IF ( ir.NE.1 .AND. ir.NE.2 .AND. ir.NE.3 &
      .AND.jc.NE.1 .AND. jc.NE.2 .AND. jc.NE.3 ) THEN
        WRITE(*,*)'Strain: invalid row/column specification: ', ir, jc
        STOP
      ENDIF

      ktmin2 = ktmin**2
      ktmax2 = ktmax**2

      IF ( ir.EQ.1 ) THEN
        CALL derivk3(vx, sij, jc)
        SELECT CASE (jc)
          CASE(1)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
            DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
              DO j = 1,ny
                DO k = 1,nz
                  ctmp(k,j,i) = sij(k,j,i)
                END DO
              END DO
            END DO
          CASE(2)
            CALL derivk3(vy, ctmp, 1)
          CASE(3)
            CALL derivk3(vz, ctmp, 1)
        END SELECT
      ELSE IF ( ir.EQ.2 ) THEN
        CALL derivk3(vy, sij, jc)
        SELECT CASE (jc)
          CASE(1)
            CALL derivk3(vx, ctmp, 2)
          CASE(2)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
            DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
              DO j = 1,ny
                DO k = 1,nz
                  ctmp(k,j,i) = sij(k,j,i)
                END DO
              END DO
            END DO
          CASE(3)
            CALL derivk3(vz, ctmp, 2)
        END SELECT
      ELSE IF ( ir.EQ.3 ) THEN
        CALL derivk3(vz, sij, jc)
        SELECT CASE (jc)
          CASE(1)
            CALL derivk3(vx, ctmp, 3)
          CASE(2)
            CALL derivk3(vy, ctmp, 3)
          CASE(3)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
            DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
              DO j = 1,ny
                DO k = 1,nz
                  ctmp(k,j,i) = sij(k,j,i)
                END DO
              END DO
            END DO
        END SELECT
      ENDIF

!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
        DO j = 1,ny
          DO k = 1,nz
            sij(k,j,i) = 0.50_GP*(sij(k,j,i)+ctmp(k,j,i)) 
          END DO
        END DO
      END DO


      ! truncate spherically:
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
        DO j = 1,ny
          DO k = 1,nz
            IF ((kk2(k,j,i).lt.ktmin2 ).or.(kk2(k,j,i).gt.ktmax2)) THEN
              sij(k,j,i) = 0.0_GP
            ENDIF
          END DO
        END DO
      END DO


      IF ( inorm.GT.0 ) THEN
        
        tmp = 1.0_GP/REAL(nx*ny*nz,KIND=GP)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
        DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
          DO j = 1,ny
            DO k = 1,nz
              sij(k,j,i) = sij(k,j,i)*tmp
            END DO
          END DO
        END DO

      ENDIF

      END SUBROUTINE Strain
!-----------------------------------------------------------------
!-----------------------------------------------------------------

