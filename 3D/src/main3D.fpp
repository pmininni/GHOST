!=================================================================
      PROGRAM MAIN3D
!=================================================================
! GHOST code: Geophysical High Order Suite for Turbulence
!
! Numerically integrates the incompressible HD/MHD/Hall-MHD 
! equations in 3 dimensions with periodic boundary conditions 
! and external forcing. A pseudo-spectral method is used to 
! compute spatial derivatives, while adjustable order 
! Runge-Kutta method is used to evolve the system in the time 
! domain. To compile, you need the FFTW library installed on 
! your system. You should link with the FFTP subroutines
! and use the FFTPLANS and MPIVARS modules (see the file 
! 'fftp_mod.f90').
!
! Notation: index 'i' is 'x' 
!           index 'j' is 'y'
!           index 'k' is 'z'
!
! Conditional compilation options:
!           HD_SOL        builds the hydrodynamic (HD) solver
!           PHD_SOL       builds the HD solver with passive scalar
!           MPHD_SOL      builds the HD solver with multiple passive scalars
!           MHD_SOL       builds the MHD solver
!           MHDB_SOL      builds the MHD solver with uniform B_0
!           HMHD_SOL      builds the Hall-MHD solver
!           HMHDB_SOL     builds the HMHD solver with uniform B_0
!           ROTH_SOL      builds the HD solver in a rotating frame
!           PROTH_SOL     builds the ROTH solver with passive scalar
!           BOUSS_SOL     builds the BOUSS solver 
!           ROTBOUSS_SOL  builds the BOUSS solver in a rotating frame 
!           LAHD_SOL      builds the Lagrangian-averaged HD solver
!           CAHD_SOL      builds the Clark-alpha HD solver
!           LHD_SOL       builds the Leray HD solver
!           LAMHD_SOL     builds the Lagrangian-averaged MHD solver
!           EDQNMHD_SOL   builds the EDQNM HD solver
!           EDQNMROTH_SOL builds the EDQNM ROTH solver
!
! 2003 Pablo D. Mininni.
!      Department of Physics, 
!      Facultad de Ciencias Exactas y Naturales.
!      Universidad de Buenos Aires.
!      e-mail: mininni@df.uba.ar
!
! 15 Feb 2007: Main program for all solvers (HD/MHD/HMHD)
! 21 Feb 2007: POSIX and MPI/IO support
! 10 Mar 2007: FFTW-2.x and FFTW-3.x support
! 25 Aug 2009: Hybrid MPI/OpenMP support (D. Rosenberg & P. Mininni)
! 30 Aug 2009: SINGLE/DOUBLE precision (D. Rosenberg & P. Mininni)
!
! References:
! Mininni PD, Rosenberg DL, Reddy R, Pouquet A.; P. Comp. 37, 123 (2011)
! Mininni PD, Gomez DO, Mahajan SM; Astrophys. J. 619, 1019 (2005)
! Gomez DO, Mininni PD, Dmitruk P; Phys. Scripta T116, 123 (2005)
! Gomez DO, Mininni PD, Dmitruk P; Adv. Sp. Res. 35, 899 (2005)
!=================================================================

!
! Definitions for conditional compilation

#ifdef HD_SOL
#define DNS_
#endif

#ifdef PHD_SOL
#define DNS_
#define SCALAR_
#endif

#ifdef MPHD_SOL
#define DNS_
#define MULTISCALAR_
#endif

#ifdef MHD_SOL
#define DNS_
#define MAGFIELD_
#endif

#ifdef MHDB_SOL
#define DNS_
#define MAGFIELD_
#define UNIFORMB_
#endif

#ifdef HMHD_SOL
#define DNS_
#define MAGFIELD_
#define HALLTERM_
#endif

#ifdef HMHDB_SOL
#define DNS_
#define MAGFIELD_
#define HALLTERM_
#define UNIFORMB_
#endif

#ifdef ROTH_SOL
#define DNS_
#define ROTATION_
#endif 

#ifdef PROTH_SOL
#define DNS_
#define SCALAR_
#define ROTATION_
#endif

#ifdef BOUSS_SOL
#define DNS_
#define SCALAR_
#define BOUSSINESQ_
#endif

#ifdef ROTBOUSS_SOL
#define DNS_
#define SCALAR_
#define ROTATION_
#define BOUSSINESQ_
#endif

#ifdef LAHD_SOL
#define ALPHAV_
#endif

#ifdef CAHD_SOL
#define ALPHAV_
#endif

#ifdef LHD_SOL
#define ALPHAV_
#endif

#ifdef LAMHD_SOL
#define MAGFIELD_
#define ALPHAV_
#define ALPHAB_
#endif

#ifdef EDQNMHD_SOL
#define DNS_
#define EDQNM_
#endif

#ifdef EDQNMROTH_SOL
#define DNS_
#define EDQNM_
#define ROTATION_
#endif

#ifdef DEF_GHOST_LAGP
#define LAGPART_           
#endif

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
      USE gtimer
#ifdef DNS_
      USE dns
#endif
#ifdef HALLTERM_
      USE hall
#endif
#ifdef ALPHAV_
      USE alpha
#endif
#ifdef EDQNM_
      USE edqnm
#endif
#if defined(DEF_GHOST_CUDA_)
      USE, INTRINSIC :: iso_c_binding
      USE cuda_bindings
      USE cutypes
#endif
#if defined(LAGPART_)
      USE class_GPart 
#endif

      IMPLICIT NONE

!
! Arrays for the fields and the external forcing

      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: vx,vy,vz
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: fx,fy,fz
#ifdef SCALAR_
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: th
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: fs
#endif
#ifdef MULTISCALAR_
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: th1,th2,th3
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: fs1,fs2,fs3
#endif
#ifdef MAGFIELD_
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: ax,ay,az
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: mx,my,mz
#endif

!
! Temporal data storage arrays

      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: C1,C2
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: C3,C4
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: C5,C6
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: C7,C8
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: M1,M2,M3
#ifdef SCALAR_
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: C20
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: M7
#endif
#ifdef MULTISCALAR_
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: C21,C22,C23,C24
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: M8,M9,M10
#endif
#ifdef MAGFIELD_
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: C9,C10,C11
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: C12,C13,C14
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: C15,C16,C17
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: M4,M5,M6
#endif
#ifdef HALLTERM_
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: C18
#endif
#ifdef EDQNM_ 
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: C19
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:)  :: tepq,thpq,tve,tvh
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:)  :: Eold,Hold
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:)  :: Eext,Hext
      REAL(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: Eden,Hden
#endif
      REAL(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: R1,R2,R3
      REAL(KIND=GP), ALLOCATABLE, DIMENSION (:)     :: Faux1,Faux2
#ifdef LAGPART_
      REAL(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: R4,R5
#endif

!
! Auxiliary variables

      COMPLEX(KIND=GP) :: cdump,jdump
      COMPLEX(KIND=GP) :: cdumq,jdumq
      COMPLEX(KIND=GP) :: cdumr,jdumr
      DOUBLE PRECISION :: tmp,tmq
      DOUBLE PRECISION :: eps,epm
!!    DOUBLE PRECISION :: cputime1,cputime2      ! rreddy@psc.edu
!!    DOUBLE PRECISION :: cputime3,cputime4      ! rreddy@psc.edu 
!!    DOUBLE PRECISION :: omptime1,omptime2
!!    DOUBLE PRECISION :: omptime3,omptime4
!$    DOUBLE PRECISION, EXTERNAL :: omp_get_wtime

      REAL(KIND=GP)    :: dt,nu,mu,kappa,bvfreq,xmom,xtemp
      REAL(KIND=GP)    :: kup,kdn
      REAL(KIND=GP)    :: rmp,rmq
      REAL(KIND=GP)    :: dump
      REAL(KIND=GP)    :: stat
      REAL(KIND=GP)    :: f0,u0
      REAL(KIND=GP)    :: phase,ampl,cort
      REAL(KIND=GP)    :: fparam0,fparam1,fparam2,fparam3,fparam4
      REAL(KIND=GP)    :: fparam5,fparam6,fparam7,fparam8,fparam9
      REAL(KIND=GP)    :: vparam0,vparam1,vparam2,vparam3,vparam4
      REAL(KIND=GP)    :: vparam5,vparam6,vparam7,vparam8,vparam9
#ifdef SCALAR_
      REAL(KIND=GP)    :: skup,skdn
      REAL(KIND=GP)    :: c0,s0
      REAL(KIND=GP)    :: cparam0,cparam1,cparam2,cparam3,cparam4
      REAL(KIND=GP)    :: cparam5,cparam6,cparam7,cparam8,cparam9
      REAL(KIND=GP)    :: sparam0,sparam1,sparam2,sparam3,sparam4
      REAL(KIND=GP)    :: sparam5,sparam6,sparam7,sparam8,sparam9
#endif
#ifdef MULTISCALAR_
      DOUBLE PRECISION :: tmp1,tmq1,tmp2,tmq2,tmp3,tmq3
      REAL(KIND=GP)    :: kappa1,kappa2,kappa3
      REAL(KIND=GP)    :: rmp1,rmq1,rmp2,rmq2,rmp3,rmq3
      REAL(KIND=GP)    :: skup,skdn
      REAL(KIND=GP)    :: c10,s10,c20,s20,c30,s30
      REAL(KIND=GP)    :: c1param0,c1param1,c1param2,c1param3,c1param4
      REAL(KIND=GP)    :: c1param5,c1param6,c1param7,c1param8,c1param9
      REAL(KIND=GP)    :: s1param0,s1param1,s1param2,s1param3,s1param4
      REAL(KIND=GP)    :: s1param5,s1param6,s1param7,s1param8,s1param9
      REAL(KIND=GP)    :: c2param0,c2param1,c2param2,c2param3,c2param4
      REAL(KIND=GP)    :: c2param5,c2param6,c2param7,c2param8,c2param9
      REAL(KIND=GP)    :: s2param0,s2param1,s2param2,s2param3,s2param4
      REAL(KIND=GP)    :: s2param5,s2param6,s2param7,s2param8,s2param9
      REAL(KIND=GP)    :: c3param0,c3param1,c3param2,c3param3,c3param4
      REAL(KIND=GP)    :: c3param5,c3param6,c3param7,c3param8,c3param9
      REAL(KIND=GP)    :: s3param0,s3param1,s3param2,s3param3,s3param4
      REAL(KIND=GP)    :: s3param5,s3param6,s3param7,s3param8,s3param9
#endif
#ifdef MAGFIELD_
      REAL(KIND=GP)    :: mkup,mkdn
      REAL(KIND=GP)    :: m0,a0
      REAL(KIND=GP)    :: mparam0,mparam1,mparam2,mparam3,mparam4
      REAL(KIND=GP)    :: mparam5,mparam6,mparam7,mparam8,mparam9
      REAL(KIND=GP)    :: aparam0,aparam1,aparam2,aparam3,aparam4
      REAL(KIND=GP)    :: aparam5,aparam6,aparam7,aparam8,aparam9
#endif
#ifdef UNIFORMB_
      REAL(KIND=GP)    :: bx0
      REAL(KIND=GP)    :: by0
      REAL(KIND=GP)    :: bz0
#endif
#ifdef ROTATION_
      REAL(KIND=GP)    :: omega
#endif

      INTEGER :: idevice, iret, ncuda, ngcuda, ppn
      INTEGER :: ini,step
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
#if defined(SCALAR_) || defined(MULTISCALAR_) 
      INTEGER :: injt
      INTEGER :: creset
#endif
#ifdef MAGFIELD_
      INTEGER :: dyna
      INTEGER :: corr
#endif
#ifdef LAGPART_
      REAL         :: rbal
      INTEGER      :: maxparts
      INTEGER      :: injtp
      INTEGER      :: cresetp
      INTEGER      :: ilginittype
      INTEGER      :: ilgintrptype
      INTEGER      :: ilgexchtype
      INTEGER      :: ilgouttype
      INTEGER      :: dolag
      TYPE (GPart) :: lagpart
#endif
!$    INTEGER, EXTERNAL :: omp_get_max_threads

#if defined(DEF_GHOST_CUDA_)
       TYPE(cudaDevicePropG) :: devprop
#endif
      TYPE(IOPLAN) :: planio
      CHARACTER (len=100) :: odir,idir

#ifdef LAGPART_
      CHARACTER(len=1024) :: lgseedfile
#endif

!
! Namelists for the input files

      NAMELIST / status / idir,odir,stat,mult,bench,outs,mean,trans
      NAMELIST / parameter / dt,step,tstep,sstep,cstep,rand,cort,seed
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
#ifdef MULTISCALAR_
      NAMELIST / mscalar / c10,s10,c20,s20,c30,s30
      NAMELIST / mscalar / skdn,skup
      NAMELIST / mscalar / kappa1,kappa2,kappa3
      NAMELIST / mscalar / c1param0,c1param1,c1param2,c1param3,c1param4 
      NAMELIST / mscalar / c1param5,c1param6,c1param7,c1param8,c1param9 
      NAMELIST / mscalar / s1param0,s1param1,s1param2,s1param3,s1param4 
      NAMELIST / mscalar / s1param5,s1param6,s1param7,s1param8,s1param9 
      NAMELIST / mscalar / c2param0,c2param1,c2param2,c2param3,c2param4 
      NAMELIST / mscalar / c2param5,c2param6,c2param7,c2param8,c2param9 
      NAMELIST / mscalar / s2param0,s2param1,s2param2,s2param3,s2param4 
      NAMELIST / mscalar / s2param5,s2param6,s2param7,s2param8,s2param9 
      NAMELIST / mscalar / c3param0,c3param1,c3param2,c3param3,c3param4 
      NAMELIST / mscalar / c3param5,c3param6,c3param7,c3param8,c3param9 
      NAMELIST / mscalar / s3param0,s3param1,s3param2,s3param3,s3param4 
      NAMELIST / mscalar / s3param5,s3param6,s3param7,s3param8,s3param9 
#endif
#if defined(SCALAR_) || defined(MULTISCALAR_)
      NAMELIST / inject / injt,creset
#endif
#ifdef MAGFIELD_
      NAMELIST / magfield / m0,a0,mkdn,mkup,mu,corr,mparam0,mparam1
      NAMELIST / magfield / mparam2,mparam3,mparam4,mparam5,mparam6
      NAMELIST / magfield / mparam7,mparam8,mparam9,aparam0,aparam1
      NAMELIST / magfield / aparam2,aparam3,aparam4,aparam5,aparam6
      NAMELIST / magfield / aparam7,aparam8,aparam9
      NAMELIST / dynamo / dyna
#endif
#ifdef UNIFORMB_
      NAMELIST / uniformb / bx0,by0,bz0
#endif
#ifdef HALLTERM_
      NAMELIST / hallparam / ep,gspe
#endif
#ifdef ROTATION_
      NAMELIST / rotation / omega
#endif
#ifdef BOUSSINESQ_
      NAMELIST / boussinesq / bvfreq,xmom,xtemp
#endif
#ifdef ALPHAV_
      NAMELIST / alphav / alpk
#endif
#ifdef ALPHAB_
      NAMELIST / alphab / alpm
#endif
#ifdef EDQNM_
      NAMELIST / edqnmles / kolmo,heli
#endif
#ifdef LAGPART_
      NAMELIST / plagpart / lgmult,maxparts,ilginittype,ilgintrptype
      NAMELIST / plagpart / ilgexchtype,ilgouttype,lgseedfile,injtp
      NAMELIST / plagpart / cresetp,dolag
#endif

! Initializes the MPI and I/O libraries
      CALL MPI_INIT_THREAD(MPI_THREAD_FUNNELED,provided,ierr)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)

! NOTE: On systems with a single GPU per node (e.g., titan)
!       we remove the following block. But on systems with 
!       multiple devices per node, this will have to be 
!       considered carefully, and possibly re-def'd:
#if 0   
  #if 0
! Initializes CUDA for Linux-based systems. This is a call to an
! NVIDIA-developed intermediary code that gets the GPU dev. no. 
! by looking in cpu_info and finding the device that resides on 
! its PCI bus:
     iret = cudaGetDeviceCount(ncuda)  ! diagnostic , for now
     CALL MPI_REDUCE(ncuda,ngcuda,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
     ppn     = G_PPN_
     idevice = -1
     iret    = setaffinity_for_nvidia(myrank,ppn,idevice)
     iret    = cudaSetDevice(idevice);
  #else
! Initializes CUDA by selecting device. The list of devices can
! be changed by the modifying the env. variable CUDA_VISIBLE_DEVICES:
     iret    = cudaGetDeviceCount(ncuda)
     idevice = mod(myrank,ncuda)
     iret    = cudaSetDevice(idevice);
     IF ( iret .EQ. cudaErrorInvalidDevice ) THEN
!    IF ( iret .NE. cudaSuccess ) THEN
       WRITE(*,*)'MAIN: Invalid CUDA device selected: ', &
       idevice, '; myrank=',myrank, '; NUM_CUDA_DEV=',ncuda
       STOP
     ENDIF
     CALL cudaGetDeviceProperties(devprop,idevice)
     IF ( devprop%major .GT. 999 ) THEN
       WRITE(*,*)'MAIN: CUDA device emulation not allowed!'
       STOP
     ENDIF
     iret = cudaGetDevice(idevice)
  #endif
#endif
      CALL range(1,n/2+1,nprocs,myrank,ista,iend)
      CALL range(1,n,nprocs,myrank,ksta,kend)
      CALL io_init(myrank,n,ksta,kend,planio)

!
! Initializes the FFT library
! Use FFTW_ESTIMATE or FFTW_MEASURE in short runs
! Use FFTW_PATIENT or FFTW_EXHAUSTIVE in long runs
! FFTW 2.x only supports FFTW_ESTIMATE or FFTW_MEASURE

      nth = 1
!$    nth = omp_get_max_threads()
#if !defined(DEF_GHOST_CUDA_)
!$    CALL fftp3d_init_threads(ierr)
#endif
      IF (bench.eq.2) THEN
         CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
         CALL GTStart(ihcpu2,GT_CPUTIME)
         CALL GTStart(ihomp2,GT_OMPTIME)
         CALL GTStart(ihwtm2,GT_WTIME)
      ENDIF
      CALL fftp3d_create_plan(planrc,n,FFTW_REAL_TO_COMPLEX, &
                             FFTW_ESTIMATE)
      CALL fftp3d_create_plan(plancr,n,FFTW_COMPLEX_TO_REAL, &
                             FFTW_ESTIMATE)
      IF (bench.eq.2) THEN
         CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
         CALL GTStop(ihcpu2)
         CALL GTStop(ihomp2)
         CALL GTStop(ihwtm2)
      ENDIF

!
! Allocates memory for distributed blocks

      ALLOCATE( C1(n,n,ista:iend),  C2(n,n,ista:iend) )
      ALLOCATE( C3(n,n,ista:iend),  C4(n,n,ista:iend) )
      ALLOCATE( C5(n,n,ista:iend),  C6(n,n,ista:iend) )
      ALLOCATE( C7(n,n,ista:iend),  C8(n,n,ista:iend) )
      ALLOCATE( vx(n,n,ista:iend) )
      ALLOCATE( vy(n,n,ista:iend) )
      ALLOCATE( vz(n,n,ista:iend) )
      ALLOCATE( fx(n,n,ista:iend) )
      ALLOCATE( fy(n,n,ista:iend) )
      ALLOCATE( fz(n,n,ista:iend) )
#ifdef SCALAR_
      ALLOCATE( C20(n,n,ista:iend) )
      ALLOCATE( th (n,n,ista:iend) )
      ALLOCATE( fs (n,n,ista:iend) )
#endif
#ifdef MULTISCALAR_
      ALLOCATE( C21(n,n,ista:iend) )
      ALLOCATE( C22(n,n,ista:iend) )
      ALLOCATE( C23(n,n,ista:iend) )
      ALLOCATE( C24(n,n,ista:iend) )
      ALLOCATE( th1(n,n,ista:iend) )
      ALLOCATE( th2(n,n,ista:iend) )
      ALLOCATE( th3(n,n,ista:iend) )
      ALLOCATE( fs1(n,n,ista:iend) )
      ALLOCATE( fs2(n,n,ista:iend) )
      ALLOCATE( fs3(n,n,ista:iend) )
#endif
#ifdef MAGFIELD_
      ALLOCATE( C9(n,n,ista:iend),  C10(n,n,ista:iend) )
      ALLOCATE( C11(n,n,ista:iend), C12(n,n,ista:iend) )
      ALLOCATE( C13(n,n,ista:iend), C14(n,n,ista:iend) )
      ALLOCATE( C15(n,n,ista:iend), C16(n,n,ista:iend) )
      ALLOCATE( C17(n,n,ista:iend) )
      ALLOCATE( ax(n,n,ista:iend) )
      ALLOCATE( ay(n,n,ista:iend) )
      ALLOCATE( az(n,n,ista:iend) )
      ALLOCATE( mx(n,n,ista:iend) )
      ALLOCATE( my(n,n,ista:iend) )
      ALLOCATE( mz(n,n,ista:iend) )
#endif
#ifdef HALLTERM_
      ALLOCATE( C18(n,n,ista:iend) )
#endif
      ALLOCATE( ka(n), ka2(n,n,ista:iend) )
      ALLOCATE( R1(n,n,ksta:kend) )
      ALLOCATE( R2(n,n,ksta:kend) )
      ALLOCATE( R3(n,n,ksta:kend) )

#ifdef LAGPART_
      ALLOCATE( R4(n,n,ksta:kend) )
      ALLOCATE( R5(n,n,ksta:kend) )
#endif

#ifdef EDQNM_
      ALLOCATE( C19(n,n,ista:iend) )
      ALLOCATE( Eden(n,n,ista:iend) )
      ALLOCATE( Hden(n,n,ista:iend) )
      ALLOCATE( tepq(n/2+1) )
      ALLOCATE( thpq(n/2+1) )
      ALLOCATE( tve (n/2+1) )
      ALLOCATE( tvh (n/2+1) )
      ALLOCATE( Eold(n/2+1) )
      ALLOCATE( Hold(n/2+1) )
      ALLOCATE( Eext(3*(n/2+1)) )
      ALLOCATE( Hext(3*(n/2+1)) )
#endif

!
! Some constants for the FFT
!     kmax: maximum truncation for dealiasing
!     tiny: minimum truncation for dealiasing

      kmax = (real(n,kind=GP)/3.0_GP)**2
#ifdef EDQNM_
      kmax = (real(n,kind=GP)/2.0_GP-0.5_GP)**2
#endif
      tiny  = 1e-5_GP
      tinyf = 1e-15_GP

!
! Builds arrays with the wavenumbers and the 
! square wavenumbers

      DO i = 1,n/2
         ka(i) = real(i-1,kind=GP)
         ka(i+n/2) = real(i-n/2-1,kind=GP)
      END DO
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,n
            DO k = 1,n
               ka2(k,j,i) = ka(i)**2+ka(j)**2+ka(k)**2
            END DO
         END DO
      END DO

! The following lines read the file 'parameter.txt'

!
! Reads general configuration flags from the namelist 
! 'status' on the external file 'parameter.txt'
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

      IF (myrank.eq.0) THEN
         OPEN(1,file='parameter.txt',status='unknown',form="formatted")
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

!
! Reads parameters that will be used to control the 
! time integration from the namelist 'parameter' on 
! the external file 'parameter.txt' 
!     dt   : time step size
!     step : total number of time steps to compute
!     tstep: number of steps between binary output
!     sstep: number of steps between power spectrum output
!     cstep: number of steps between output of global quantities
!     rand : = 0 constant force
!            = 1 random phases
!            = 2 constant energy
!     cort : time correlation of the external forcing
!     seed : seed for the random number generator

      IF (myrank.eq.0) THEN
         OPEN(1,file='parameter.txt',status='unknown',form="formatted")
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
! Reads parameters for the velocity field from the 
! namelist 'velocity' on the external file 'parameter.txt' 
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
         OPEN(1,file='parameter.txt',status='unknown',form="formatted")
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

#ifdef BOUSSINESQ_
!
! Reads parameters specifically for Boussinesq solver from the 
! namelist 'boussinesq' on the external file 'parameter.txt'
!     bvfreq: Brunt-Vaisala frequency (positive definite)
!     xmom  : multiplies bouyancy term in momentum equation
!     xtemp : multiplies temperature-current term in 
!             temperature/density equation
      xmom  = 1.0
      xtemp = 1.0
      IF (myrank.eq.0) THEN
         OPEN(1,file='parameter.txt',status='unknown',form="formatted")
         READ(1,NML=boussinesq)
         CLOSE(1)
      ENDIF
      CALL MPI_BCAST(bvfreq,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(xmom  ,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(xtemp ,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      xmom  = xmom  * bvfreq
      xtemp = xtemp * bvfreq
#endif


#ifdef SCALAR_
!
! Reads general configuration flags for runs with 
! a passive/active scalar from the namelist 'inject' 
! on the external file 'parameter.txt'
!     injt : = 0 when stat=0 generates initial v and th (SCALAR_)
!            = 1 when stat.ne.0 imports v and generates th (SCALAR_)
!     creset: = 0: don't reset counters; 1 = reset counters 

      injt   = 0
      creset = 1
      IF (myrank.eq.0) THEN
         OPEN(1,file='parameter.txt',status='unknown',form="formatted")
         READ(1,NML=inject)
         CLOSE(1)
      ENDIF
      CALL MPI_BCAST(injt  ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(creset,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

!
! Reads parameters for the passive/active scalar from the 
! namelist 'scalar' on the external file 'parameter.txt'
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
         OPEN(1,file='parameter.txt',status='unknown',form="formatted")
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

#ifdef MULTISCALAR_
!
! Reads general configuration flags for runs with 
! a passive/active scalar from the namelist 'inject' 
! on the external file 'parameter.txt'
!     injt : = 0 when stat=0 generates initial v and th (SCALAR_)
!            = 1 when stat.ne.0 imports v and generates th (SCALAR_)
!     creset: = 0: don't reset counters; 1 = reset counters 

      injt   = 0
      creset = 1
      IF (myrank.eq.0) THEN
         OPEN(1,file='parameter.txt',status='unknown',form="formatted")
         READ(1,NML=inject)
         CLOSE(1)
      ENDIF
      CALL MPI_BCAST(injt  ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(creset,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

!
! Reads parameters for the passive/active scalar from the 
! namelist 'mscalar' on the external file 'parameter.txt'
!     si0   : amplitude of the passive scalar source i
!     ci0   : amplitude of the initial concentration i
!     skdn  : minimum wave number in concentration/source
!     skup  : maximum wave number in concentration/source
!     kappa1: diffusivity for scalars 1
!     kappa2: diffusivity for scalars 2
!     kappa3: diffusivity for scalars 3
!     sparam0-9 : ten real numbers to control properties of 
!            the source
!     cparam0-9 : ten real numbers to control properties of
!            the initial concentration

      IF (myrank.eq.0) THEN
         OPEN(1,file='parameter.txt',status='unknown',form="formatted")
         READ(1,NML=mscalar)
         CLOSE(1)
      ENDIF
      CALL MPI_BCAST(c10     ,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(c20     ,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(c30     ,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(s10     ,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(s20     ,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(s30     ,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(skdn   ,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(skup   ,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(kappa1  ,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(kappa2  ,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(kappa3  ,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(s1param0,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(s1param1,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(s1param2,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(s1param3,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(s1param4,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(s1param5,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(s1param6,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(s1param7,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(s1param8,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(s1param9,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(c1param0,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(c1param1,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(c1param2,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(c1param3,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(c1param4,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(c1param5,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(c1param6,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(c1param7,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(c1param8,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(c1param9,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(s2param0,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(s2param1,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(s2param2,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(s2param3,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(s2param4,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(s2param5,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(s2param6,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(s2param7,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(s2param8,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(s2param9,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(c2param0,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(c2param1,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(c2param2,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(c2param3,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(c2param4,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(c2param5,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(c2param6,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(c2param7,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(c2param8,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(c2param9,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(s3param0,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(s3param1,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(s3param2,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(s3param3,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(s3param4,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(s3param5,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(s3param6,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(s3param7,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(s3param8,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(s3param9,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(c3param0,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(c3param1,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(c3param2,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(c3param3,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(c3param4,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(c3param5,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(c3param6,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(c3param7,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(c3param8,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(c3param9,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
#endif

#ifdef MAGFIELD_
!
! Reads general configuration flags for runs with 
! magnetic fields from the namelist 'dynamo' on 
! the external file 'parameter.txt'
!     dyna : = 0 when stat=0 generates initial v and B (MAGFIELD_)
!            = 1 when stat.ne.0 imports v and generates B (MAGFIELD_) 

      IF (myrank.eq.0) THEN
         OPEN(1,file='parameter.txt',status='unknown',form="formatted")
         READ(1,NML=dynamo)
         CLOSE(1)
      ENDIF
      CALL MPI_BCAST(dyna,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

!
! Reads parameters for the magnetic field from the 
! namelist 'magfield' on the external file 'parameter.txt' 
!     m0   : amplitude of the electromotive forcing
!     a0   : amplitude of the initial vector potential
!     mkdn : minimum wave number in B/electromotive forcing
!     mkup : maximum wave number in B/electromotive forcing
!     mu   : magnetic diffusivity
!     corr : = 0 no correlation between the random phases
!            = 1 correlation in the random phases generator
!     mparam0-9 : ten real numbers to control properties of 
!            the electromotive forcing
!     aparam0-9 : ten real numbers to control properties of
!            the initial conditions for the magnetic field

      IF (myrank.eq.0) THEN
         OPEN(1,file='parameter.txt',status='unknown',form="formatted")
         READ(1,NML=magfield)
         CLOSE(1)
      ENDIF
      CALL MPI_BCAST(m0,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(a0,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(mkdn,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(mkup,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(mu,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(corr,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(mparam0,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(mparam1,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(mparam2,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(mparam3,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(mparam4,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(mparam5,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(mparam6,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(mparam7,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(mparam8,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(mparam9,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(aparam0,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(aparam1,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(aparam2,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(aparam3,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(aparam4,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(aparam5,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(aparam6,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(aparam7,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(aparam8,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(aparam9,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
#endif

#ifdef UNIFORMB_
!
! Reads parameters for runs with a uniform magnetic 
! field from the namelist 'uniformb' on the external 
! file 'parameter.txt' 
!     bx0: uniform magnetic field in x
!     by0: uniform magnetic field in y
!     bz0: uniform magnetic field in z

      IF (myrank.eq.0) THEN
         OPEN(1,file='parameter.txt',status='unknown',form="formatted")
         READ(1,NML=uniformb)
         CLOSE(1)
      ENDIF
      CALL MPI_BCAST(bx0,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(by0,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(bz0,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
#endif

#ifdef HALLTERM_
!
! Reads parameters for runs with the Hall effect 
! from the namelist 'hallparam' on the external 
! file 'parameter.txt' 
!     ep  : amplitude of the Hall effect
!     gspe: = 0 skips generalized helicity spectrum computation
!           = 1 computes the spectrum of generalized helicity

      IF (myrank.eq.0) THEN
         OPEN(1,file='parameter.txt',status='unknown',form="formatted")
         READ(1,NML=hallparam)
         CLOSE(1)
      ENDIF
      CALL MPI_BCAST(ep,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(gspe,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#endif

#ifdef ROTATION_
!
! Reads parameters for runs with rotation from the 
! namelist 'rotation' on the external file 'parameter.txt'
!     omega: amplitude of the uniform rotation

      IF (myrank.eq.0) THEN
         OPEN(1,file='parameter.txt',status='unknown',form="formatted")
         READ(1,NML=rotation)
         CLOSE(1)
      ENDIF
      CALL MPI_BCAST(omega,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
#endif

#ifdef ALPHAV_
!
! Reads the value of alpha for the velocity field 
! in runs using Lagrangian averaged subgrid models
!     alpk: filter length for the velocity field

      IF (myrank.eq.0) THEN
         OPEN(1,file='parameter.txt',status='unknown',form="formatted")
         READ(1,NML=alphav)
         CLOSE(1)
      ENDIF
      CALL MPI_BCAST(alpk,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
#endif

#ifdef ALPHAB_
!
! Reads the value of alpha for the magnetic field 
! in runs using Lagrangian averaged subgrid models
!     alpm: filter length for the magnetic field

      IF (myrank.eq.0) THEN
         OPEN(1,file='parameter.txt',status='unknown',form="formatted")
         READ(1,NML=alphab)
         CLOSE(1)
      ENDIF
      CALL MPI_BCAST(alpm,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
#endif

#ifdef EDQNM_
!
! Reads the value of the Kolmogorov constant and a 
! flag for the helicity LES in runs using EDQNM-based 
! LES models
!     kolmo: Kolmogorov constant
!     heli:  = 0 helicity not taken into account
!            = 1 helicity taken into account

      kolmo = 1.4_GP !Default value
      heli = 0       !Default value
      IF (myrank.eq.0) THEN
         OPEN(1,file='parameter.txt',status='unknown',form="formatted")
         READ(1,NML=edqnmles)
         CLOSE(1)
      ENDIF
      CALL MPI_BCAST(kolmo,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(heli,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#endif

#ifdef LAGPART_
      maxparts     = 1000
!     injtp: = 0 when stat=0 generates initial v and initial part seeds 
!            = 1 when stat.ne.0 imports v and generates initial 
!              part seeds. If injtp=0 when stat.ne.0, then the 
!              particle restart file is read that corresp. to 
!              stat.
!     cresetp= 0: don't reset counters when injtp=1;
!            = 1: _do_ reset counters when injtp=1.
!     ilginittype : Inititialization type. either GPINIT_RANDLOC or GPINIT_USERLOC
!     ilgintrptype: Interpolation type: only GPINTRP_CSPLINE currently
!     ilgexchtype : Boundary exchange type: GPEXCHTYPE_VDB (voxel db) or GPEXCHTYPE_NN (nearest-neighbor)
!     ilgouttype  : Particle output type: 0 = binary; 1= ASCII
!     lgseedfile  : Name of seed file if ilginittype=GPINIT_USERLOC
!     dolag       : 1 = run with particles; 0 = don't 
      injtp        = 0
      cresetp      = 0
      ilginittype  = GPINIT_RANDLOC
      ilgintrptype = GPINTRP_CSPLINE
      ilgexchtype  = GPEXCHTYPE_VDB
      ilgouttype   = 0
      lgmult       = 1
      lgseedfile   = 'gplag.dat'
      rbal         = 0.0
      dolag        = 1
!
      IF (myrank.eq.0) THEN
         OPEN(1,file='parameter.txt',status='unknown',form="formatted")
         READ(1,NML=plagpart)
         CLOSE(1)
      ENDIF
      CALL MPI_BCAST(maxparts    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(injtp       ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(cresetp     ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(lgmult      ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(ilginittype ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(ilgintrptype,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(ilgexchtype ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(ilgouttype  ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(dolag       ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(lgseedfile,1024,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      IF ( mod(tstep,lgmult).NE.0 ) THEN
        WRITE(*,*)'main: lgmult must divide tstep evenly'
        STOP
      ENDIF  
      pstep = tstep/lgmult
#endif

!
! Initializes arrays to keep track of the forcing 
! amplitude if constant energy is used, and 
! allocates arrays to compute mean flows if needed.

      ampl = 1.0_GP
      timef = fstep
      IF (rand.eq.2) THEN
         ALLOCATE( Faux1(10) )
         ALLOCATE( Faux2(10) )
         DO i = 1,10
            Faux1(i) = ampl
            Faux2(i) = 0.0_GP
         END DO
      ENDIF
      IF (mean.eq.1) THEN
         ALLOCATE( M1(n,n,ista:iend) )
         ALLOCATE( M2(n,n,ista:iend) )
         ALLOCATE( M3(n,n,ista:iend) )
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
            DO j = 1,n
               DO k = 1,n
                  M1(k,j,i) = 0.0_GP
                  M2(k,j,i) = 0.0_GP
                  M3(k,j,i) = 0.0_GP
               END DO
            END DO
         END DO
#ifdef SCALAR_
         ALLOCATE( M7(n,n,ista:iend) )
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
            DO j = 1,n
               DO k = 1,n
                  M7(k,j,i) = 0.0_GP
               END DO
            END DO
         END DO
#endif
#ifdef MULTISCALAR_
         ALLOCATE( M8 (n,n,ista:iend) )
         ALLOCATE( M9 (n,n,ista:iend) )
         ALLOCATE( M10(n,n,ista:iend) )
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
            DO j = 1,n
               DO k = 1,n
                  M8 (k,j,i) = 0.0_GP
                  M9 (k,j,i) = 0.0_GP
                  M10(k,j,i) = 0.0_GP
               END DO
            END DO
         END DO
#endif
#ifdef MAGFIELD_
         ALLOCATE( M4(n,n,ista:iend) )
         ALLOCATE( M5(n,n,ista:iend) )
         ALLOCATE( M6(n,n,ista:iend) )
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
            DO j = 1,n
               DO k = 1,n
                  M4(k,j,i) = 0.0_GP
                  M5(k,j,i) = 0.0_GP
                  M6(k,j,i) = 0.0_GP
               END DO
            END DO
         END DO
#endif
      ENDIF

#ifdef LAGPART_
      IF ( dolag.GT.0 ) THEN
        CALL lagpart%GPart_ctor(MPI_COMM_WORLD,maxparts,ilginittype,&
             ilgintrptype,3,ilgexchtype,ilgouttype,csize,nstrip)
        CALL lagpart%SetRandSeed(seed)
        CALL lagpart%SetSeedFile(trim(lgseedfile))
      ENDIF
#endif

!
! Sets the external forcing

      INCLUDE 'initialfv.f90'           ! mechanical forcing
#ifdef SCALAR_
      INCLUDE 'initialfs.f90'           ! scalar source
#endif
#ifdef MUJLTISCALAR_
      INCLUDE 'initialfms.f90'          ! multiscalar sources
#endif
#ifdef MAGFIELD_
      INCLUDE 'initialfb.f90'           ! electromotive forcing
#endif

! If stat=0 we start a new run.
! Generates initial conditions for the fields.
 IC : IF (stat.eq.0) THEN

      ini = 1
      sind = 0                          ! index for the spectrum
      tind = 0                          ! index for the binaries
      pind = 0                          ! index for the particles
      timet = tstep
      timec = cstep
      times = sstep
      timep = pstep
      INCLUDE 'initialv.f90'            ! initial velocity
#ifdef SCALAR_
      INCLUDE 'initials.f90'            ! initial concentration
#endif
#ifdef MULTISCALAR_
      INCLUDE 'initialms.f90'           ! multiscalar initial concentration
#endif
#ifdef MAGFIELD_
      INCLUDE 'initialb.f90'            ! initial vector potential
#endif
#ifdef LAGPART_
      IF ( dolag.GT.0 ) THEN
        CALL lagpart%Init()
      ENDIF
#endif

      ELSE

! If stat.ne.0 a previous run is continued

      ini = int((stat-1)*tstep)
      tind = int(stat)
      sind = int(real(ini,kind=GP)/real(sstep,kind=GP)+1)
!     pind = int(real(ini,kind=GP)/real(pstep,kind=GP)+1)
      pind = int((stat-1)*lgmult+1)
      WRITE(ext, fmtext) tind
      times = 0
      timet = 0
      timec = 0
      timep = 0

      CALL io_read(1,idir,'vx',ext,planio,R1)
      CALL io_read(1,idir,'vy',ext,planio,R2)
      CALL io_read(1,idir,'vz',ext,planio,R3)
      CALL fftp3d_real_to_complex(planrc,R1,vx,MPI_COMM_WORLD)
      CALL fftp3d_real_to_complex(planrc,R2,vy,MPI_COMM_WORLD)
      CALL fftp3d_real_to_complex(planrc,R3,vz,MPI_COMM_WORLD)

      IF (mean.eq.1) THEN
         CALL io_read(1,idir,'mean_vx',ext,planio,R1)
         CALL io_read(1,idir,'mean_vy',ext,planio,R2)
         CALL io_read(1,idir,'mean_vz',ext,planio,R3)
         CALL fftp3d_real_to_complex(planrc,R1,M1,MPI_COMM_WORLD)
         CALL fftp3d_real_to_complex(planrc,R2,M2,MPI_COMM_WORLD)
         CALL fftp3d_real_to_complex(planrc,R3,M3,MPI_COMM_WORLD)
         dump = real(ini,kind=GP)/cstep
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
            DO j = 1,n
               DO k = 1,n
                  M1(k,j,i) = dump*M1(k,j,i)
                  M2(k,j,i) = dump*M2(k,j,i)
                  M3(k,j,i) = dump*M3(k,j,i)
               END DO
            END DO
         END DO
      ENDIF

#ifdef SCALAR_
 INJ: IF (injt.eq.0) THEN
         CALL io_read(1,idir,'th',ext,planio,R1)
         CALL fftp3d_real_to_complex(planrc,R1,th,MPI_COMM_WORLD)
         IF (mean.eq.1) THEN
            CALL io_read(1,idir,'mean_th',ext,planio,R1)
            CALL fftp3d_real_to_complex(planrc,R1,M7,MPI_COMM_WORLD)
            dump = real(ini,kind=GP)/cstep
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
            DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
               DO j = 1,n
                  DO k = 1,n
                     M7(k,j,i) = dump*M7(k,j,i)
                 END DO
               END DO
            END DO
         ENDIF
      ELSE
         INCLUDE 'initials.f90'      ! initial concentration
         IF (creset.ne.0) THEN
         ini = 1                     ! resets all counters (the
         sind = 0                    ! run starts at t=0)
         tind = 0
         pind = 0
         timet = tstep
         timec = cstep
         times = sstep
         timep = pstep
         ENDIF
      ENDIF INJ
#endif

#ifdef MULTISCALAR_
 INJ: IF (injt.eq.0) THEN
         CALL io_read(1,idir,'th1',ext,planio,R1)
         CALL fftp3d_real_to_complex(planrc,R1,th1,MPI_COMM_WORLD)
         CALL io_read(1,idir,'th2',ext,planio,R1)
         CALL fftp3d_real_to_complex(planrc,R1,th1,MPI_COMM_WORLD)
         CALL io_read(1,idir,'th3',ext,planio,R1)
         CALL fftp3d_real_to_complex(planrc,R1,th1,MPI_COMM_WORLD)
         IF (mean.eq.1) THEN
            CALL io_read(1,idir,'mean_th1',ext,planio,R1)
            CALL fftp3d_real_to_complex(planrc,R1,M8 ,MPI_COMM_WORLD)
            CALL io_read(1,idir,'mean_th2',ext,planio,R1)
            CALL fftp3d_real_to_complex(planrc,R1,M9 ,MPI_COMM_WORLD)
            CALL io_read(1,idir,'mean_th3',ext,planio,R1)
            CALL fftp3d_real_to_complex(planrc,R1,M10,MPI_COMM_WORLD)
            dump = real(ini,kind=GP)/cstep
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
            DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
               DO j = 1,n
                  DO k = 1,n
                     M8 (k,j,i) = dump*M8 (k,j,i)
                     M9 (k,j,i) = dump*M9 (k,j,i)
                     M10(k,j,i) = dump*M10(k,j,i)
                 END DO
               END DO
            END DO
         ENDIF
      ELSE
         INCLUDE 'initialms.f90'     ! initial concentrations
         IF (creset.ne.0) THEN
         ini = 1                     ! resets all counters (the
         sind = 0                    ! run starts at t=0)
         tind = 0
         pind = 0
         timet = tstep
         timec = cstep
         times = sstep
         timep = pstep
         ENDIF
      ENDIF INJ
#endif

#ifdef MAGFIELD_
 DYN: IF (dyna.eq.0) THEN
         CALL io_read(1,idir,'ax',ext,planio,R1)
         CALL io_read(1,idir,'ay',ext,planio,R2)
         CALL io_read(1,idir,'az',ext,planio,R3)
         CALL fftp3d_real_to_complex(planrc,R1,ax,MPI_COMM_WORLD)
         CALL fftp3d_real_to_complex(planrc,R2,ay,MPI_COMM_WORLD)
         CALL fftp3d_real_to_complex(planrc,R3,az,MPI_COMM_WORLD)
         IF (mean.eq.1) THEN
            CALL io_read(1,idir,'mean_bx',ext,planio,R1)
            CALL io_read(1,idir,'mean_by',ext,planio,R2)
            CALL io_read(1,idir,'mean_bz',ext,planio,R3)
            CALL fftp3d_real_to_complex(planrc,R1,M4,MPI_COMM_WORLD)
            CALL fftp3d_real_to_complex(planrc,R2,M5,MPI_COMM_WORLD)
            CALL fftp3d_real_to_complex(planrc,R3,M6,MPI_COMM_WORLD)
            dump = real(ini,kind=GP)/cstep
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
            DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
               DO j = 1,n
                  DO k = 1,n
                     M4(k,j,i) = dump*M4(k,j,i)
                     M5(k,j,i) = dump*M5(k,j,i)
                     M6(k,j,i) = dump*M6(k,j,i)
                 END DO
               END DO
            END DO
         ENDIF
      ELSE
         INCLUDE 'initialb.f90'      ! initial vector potential
         ini = 1                     ! resets all counters (the
         sind = 0                    ! dynamo run starts at t=0)
         tind = 0
         pind = 0
         timet = tstep
         timec = cstep
         times = sstep
         timep = pstep
      ENDIF DYN
#endif

#ifdef LAGPART_
      IF ( dolag.GT.0 ) THEN
        IF (injtp.eq.0) THEN
          WRITE(ext, fmtext) pind
          CALL lagpart%io_read(1,idir,'xlg',ext)
        ELSE
          CALL lagpart%Init()
          IF (cresetp.ne.0) THEN
            ini = 1                   ! resets all counters (the
            sind = 0                  ! particle run starts at t=0)
            tind = 0
            pind = 0
            timet = tstep
            timec = cstep
            times = sstep
            timep = pstep
          ENDIF
        ENDIF
      ENDIF
#endif

      ENDIF IC

!
! Time integration scheme starts here.
! Does ord iterations of Runge-Kutta. If 
! we are doing a benchmark, we measure 
! cputime before starting.

      IF (bench.eq.1) THEN
         CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
         CALL GTStart(ihcpu1,GT_CPUTIME)
         CALL GTStart(ihomp1,GT_OMPTIME)
         CALL GTStart(ihwtm1,GT_WTIME)
! Must reset FFT internal timers after initialization:
         comtime = 0.0D0
         ffttime = 0.0D0
         tratime = 0.0D0
#if defined(DEF_GHOST_CUDA_)
         memtime = 0.0D0
#endif

      ENDIF

 RK : DO t = ini,step

! Updates the external forcing. Every 'fsteps'
! the phase or amplitude is changed according 
! to the value of 'rand'.
         IF (timef.eq.fstep) THEN
            timef = 0

            IF (rand.eq.1) THEN      ! randomizes phases

               IF (myrank.eq.0) phase = 2*pi*randu(seed)
               CALL MPI_BCAST(phase,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
               cdump = COS(phase)+im*SIN(phase)
               jdump = conjg(cdump)
#if defined(SCALAR_) || defined(MULTISCALAR_)
               IF (myrank.eq.0) phase = 2*pi*randu(seed)
               CALL MPI_BCAST(phase,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
               cdumr = COS(phase)+im*SIN(phase)
               jdumr = conjg(cdumr)
#endif
#ifdef MAGFIELD_
               IF (myrank.eq.0) phase = 2*pi*randu(seed)
               CALL MPI_BCAST(phase,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
               cdumq = corr*cdump+(1-corr)*(COS(phase)+im*SIN(phase))
               jdumq = corr*jdump+(1-corr)*(COS(phase)-im*SIN(phase))
#endif
               IF (ista.eq.1) THEN
!$omp parallel do
                  DO j = 2,n/2+1
                     fx(1,j,1) = fx(1,j,1)*cdump
                     fx(1,n-j+2,1) = fx(1,n-j+2,1)*jdump
                     fy(1,j,1) = fy(1,j,1)*cdump
                     fy(1,n-j+2,1) = fy(1,n-j+2,1)*jdump
                     fz(1,j,1) = fz(1,j,1)*cdump
                     fz(1,n-j+2,1) = fz(1,n-j+2,1)*jdump
#ifdef SCALAR_
                     fs(1,j,1) = fs(1,j,1)*cdumr
                     fs(1,n-j+2,1) = fs(1,n-j+2,1)*jdumr
#endif
#ifdef MULTISCALAR_
                     fs1(1,j,1) = fs1(1,j,1)*cdumr
                     fs1(1,n-j+2,1) = fs1(1,n-j+2,1)*jdumr
                     fs2(1,j,1) = fs2(1,j,1)*cdumr
                     fs2(1,n-j+2,1) = fs2(1,n-j+2,1)*jdumr
                     fs3(1,j,1) = fs3(1,j,1)*cdumr
                     fs3(1,n-j+2,1) = fs3(1,n-j+2,1)*jdumr
#endif
#ifdef MAGFIELD_
                     mx(1,j,1) = mx(1,j,1)*cdumq
                     mx(1,n-j+2,1) = mx(1,n-j+2,1)*jdumq
                     my(1,j,1) = my(1,j,1)*cdumq
                     my(1,n-j+2,1) = my(1,n-j+2,1)*jdumq
                     mz(1,j,1) = mz(1,j,1)*cdumq
                     mz(1,n-j+2,1) = mz(1,n-j+2,1)*jdumq
#endif
                  END DO
!$omp parallel do
                  DO k = 2,n/2+1
                     fx(k,1,1) = fx(k,1,1)*cdump
                     fx(n-k+2,1,1) = fx(n-k+2,1,1)*jdump
                     fy(k,1,1) = fy(k,1,1)*cdump
                     fy(n-k+2,1,1) = fy(n-k+2,1,1)*jdump
                     fz(k,1,1) = fz(k,1,1)*cdump
                     fz(n-k+2,1,1) = fz(n-k+2,1,1)*jdump
#ifdef SCALAR_
                     fs(k,1,1) = fs(k,1,1)*cdumr
                     fs(n-k+2,1,1) = fs(n-k+2,1,1)*jdumr
#endif
#ifdef MULTISCALAR_
                     fs1(k,1,1) = fs1(k,1,1)*cdumr
                     fs1(n-k+2,1,1) = fs1(n-k+2,1,1)*jdumr
                     fs2(k,1,1) = fs2(k,1,1)*cdumr
                     fs2(n-k+2,1,1) = fs2(n-k+2,1,1)*jdumr
                     fs3(k,1,1) = fs3(k,1,1)*cdumr
                     fs3(n-k+2,1,1) = fs3(n-k+2,1,1)*jdumr
#endif
#ifdef MAGFIELD_
                     mx(k,1,1) = mx(k,1,1)*cdumq
                     mx(n-k+2,1,1) = mx(n-k+2,1,1)*jdumq
                     my(k,1,1) = my(k,1,1)*cdumq
                     my(n-k+2,1,1) = my(n-k+2,1,1)*jdumq
                     mz(k,1,1) = mz(k,1,1)*cdumq
                     mz(n-k+2,1,1) = mz(n-k+2,1,1)*jdumq
#endif
                  END DO
!$omp parallel do private (k)
                  DO j = 2,n
                     DO k = 2,n/2+1
                        fx(k,j,1) = fx(k,j,1)*cdump
                        fx(n-k+2,n-j+2,1) = fx(n-k+2,n-j+2,1)*jdump
                        fy(k,j,1) = fy(k,j,1)*cdump
                        fy(n-k+2,n-j+2,1) = fy(n-k+2,n-j+2,1)*jdump
                        fz(k,j,1) = fz(k,j,1)*cdump
                        fz(n-k+2,n-j+2,1) = fz(n-k+2,n-j+2,1)*jdump
#ifdef SCALAR_
                        fs(k,j,1) = fs(k,j,1)*cdumr
                        fs(n-k+2,n-j+2,1) = fs(n-k+2,n-j+2,1)*jdumr
#endif
#ifdef MULTISCALAR_
                        fs1(k,j,1) = fs1(k,j,1)*cdumr
                        fs1(n-k+2,n-j+2,1) = fs1(n-k+2,n-j+2,1)*jdumr
                        fs2(k,j,1) = fs2(k,j,1)*cdumr
                        fs2(n-k+2,n-j+2,1) = fs2(n-k+2,n-j+2,1)*jdumr
                        fs3(k,j,1) = fs3(k,j,1)*cdumr
                        fs3(n-k+2,n-j+2,1) = fs3(n-k+2,n-j+2,1)*jdumr
#endif
#ifdef MAGFIELD_
                        mx(k,j,1) = mx(k,j,1)*cdumq
                        mx(n-k+2,n-j+2,1) = mx(n-k+2,n-j+2,1)*jdumq
                        my(k,j,1) = my(k,j,1)*cdumq
                        my(n-k+2,n-j+2,1) = my(n-k+2,n-j+2,1)*jdumq
                        mz(k,j,1) = mz(k,j,1)*cdumq
                        mz(n-k+2,n-j+2,1) = mz(n-k+2,n-j+2,1)*jdumq
#endif
                     END DO
                  END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k)
                  DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k)
                     DO j = 1,n
                        DO k = 1,n
                           fx(k,j,i) = fx(k,j,i)*cdump
                           fy(k,j,i) = fy(k,j,i)*cdump
                           fz(k,j,i) = fz(k,j,i)*cdump
#ifdef SCALAR_
                           fs(k,j,i) = fs(k,j,i)*cdumr
#endif
#ifdef MULTISCALAR_
                           fs1(k,j,i) = fs1(k,j,i)*cdumr
                           fs2(k,j,i) = fs2(k,j,i)*cdumr
                           fs3(k,j,i) = fs3(k,j,i)*cdumr
#endif
#ifdef MAGFIELD_
                           mx(k,j,i) = mx(k,j,i)*cdumq
                           my(k,j,i) = my(k,j,i)*cdumq
                           mz(k,j,i) = mz(k,j,i)*cdumq
#endif
                        END DO
                     END DO
                  END DO
               ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
                  DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
                     DO j = 1,n
                        DO k = 1,n
                           fx(k,j,i) = fx(k,j,i)*cdump
                           fy(k,j,i) = fy(k,j,i)*cdump
                           fz(k,j,i) = fz(k,j,i)*cdump
#ifdef SCALAR_
                           fs(k,j,i) = fs(k,j,i)*cdumr
#endif
#ifdef MULTISCALAR_
                           fs1(k,j,i) = fs1(k,j,i)*cdumr
                           fs2(k,j,i) = fs2(k,j,i)*cdumr
                           fs3(k,j,i) = fs3(k,j,i)*cdumr
#endif
#ifdef MAGFIELD_
                           mx(k,j,i) = mx(k,j,i)*cdumq
                           my(k,j,i) = my(k,j,i)*cdumq
                           mz(k,j,i) = mz(k,j,i)*cdumq
#endif
                        END DO
                     END DO
                  END DO
               ENDIF

            ELSE IF (rand.eq.2) THEN ! constant energy

#ifdef HD_SOL
              INCLUDE 'hd_adjustfv.f90'
#endif
#ifdef PHD_SOL
              INCLUDE 'hd_adjustfv.f90'
#endif
#ifdef MHD_SOL
              INCLUDE 'mhd_adjustfv.f90'
#endif
#ifdef MHDB_SOL
              INCLUDE 'mhd_adjustfv.f90'
#endif
#ifdef HMHD_SOL
              INCLUDE 'mhd_adjustfv.f90'
#endif
#ifdef HMHDB_SOL
              INCLUDE 'mhd_adjustfv.f90'
#endif
#ifdef ROTH_SOL
              INCLUDE 'hd_adjustfv.f90'
#endif
#ifdef PROTH_SOL
              INCLUDE 'hd_adjustfv.f90'
#endif
#if defined(BOUSS_SOL) || defined(ROTBOUSS_SOL)
              INCLUDE 'hd_adjustfv.f90'
#endif
#ifdef LAHD_SOL
              INCLUDE 'lahd_adjustfv.f90'
#endif
#ifdef CAHD_SOL
              INCLUDE 'lahd_adjustfv.f90'
#endif
#ifdef LHD_SOL
              INCLUDE 'hd_adjustfv.f90'
#endif
#ifdef LAMHD_SOL
              INCLUDE 'lahd_adjustfv.f90'
#endif
#ifdef EDQNMHD_SOL
              INCLUDE 'edqnmhd_adjustfv.f90'
#endif
#ifdef EDQNMROTH_SOL
              INCLUDE 'edqnmhd_adjustfv.f90'
#endif

            ENDIF

         ENDIF

! Every 'tstep' steps, stores the fields 
! in binary files

         IF ((timet.eq.tstep).and.(bench.eq.0)) THEN
            timet = 0
            tind = tind+1
            IF (rand.eq.2) THEN
               CALL energy(fx,fy,fz,tmp,1)
               IF (myrank.eq.0) THEN
                  OPEN(1,file='force.txt',position='append')
                  WRITE(1,*) (t-1)*dt,sqrt(tmp)
                  CLOSE(1)
               ENDIF
            ENDIF
            WRITE(ext, fmtext) tind
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
            DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
               DO j = 1,n
                  DO k = 1,n
                     C1(k,j,i) = vx(k,j,i)/real(n,kind=GP)**3
                     C2(k,j,i) = vy(k,j,i)/real(n,kind=GP)**3
                     C3(k,j,i) = vz(k,j,i)/real(n,kind=GP)**3
                  END DO
               END DO
            END DO
            IF (outs.ge.1) THEN
               CALL rotor3(C2,C3,C4,1)
               CALL rotor3(C1,C3,C5,2)
               CALL rotor3(C1,C2,C6,3)
               CALL fftp3d_complex_to_real(plancr,C4,R1,MPI_COMM_WORLD)
               CALL fftp3d_complex_to_real(plancr,C5,R2,MPI_COMM_WORLD)
               CALL fftp3d_complex_to_real(plancr,C6,R3,MPI_COMM_WORLD)
               CALL io_write(1,odir,'wx',ext,planio,R1)
               CALL io_write(1,odir,'wy',ext,planio,R2)
               CALL io_write(1,odir,'wz',ext,planio,R3)
            ENDIF
            CALL fftp3d_complex_to_real(plancr,C1,R1,MPI_COMM_WORLD)
            CALL fftp3d_complex_to_real(plancr,C2,R2,MPI_COMM_WORLD)
            CALL fftp3d_complex_to_real(plancr,C3,R3,MPI_COMM_WORLD)
            CALL io_write(1,odir,'vx',ext,planio,R1)
            CALL io_write(1,odir,'vy',ext,planio,R2)
            CALL io_write(1,odir,'vz',ext,planio,R3)
            IF (mean.eq.1) THEN
               dump = real(cstep,kind=GP)/t
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
               DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
                  DO j = 1,n
                     DO k = 1,n
                        C1(k,j,i) = dump*M1(k,j,i)/real(n,kind=GP)**3
                        C2(k,j,i) = dump*M2(k,j,i)/real(n,kind=GP)**3
                        C3(k,j,i) = dump*M3(k,j,i)/real(n,kind=GP)**3
                     END DO
                  END DO
               END DO
               CALL fftp3d_complex_to_real(plancr,C1,R1,MPI_COMM_WORLD)
               CALL fftp3d_complex_to_real(plancr,C2,R2,MPI_COMM_WORLD)
               CALL fftp3d_complex_to_real(plancr,C3,R3,MPI_COMM_WORLD)
               CALL io_write(1,odir,'mean_vx',ext,planio,R1)
               CALL io_write(1,odir,'mean_vy',ext,planio,R2)
               CALL io_write(1,odir,'mean_vz',ext,planio,R3)
            ENDIF
#ifdef SCALAR_
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
            DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
               DO j = 1,n
                  DO k = 1,n
                     C1(k,j,i) = th(k,j,i)/real(n,kind=GP)**3
                  END DO
               END DO
            END DO
            CALL fftp3d_complex_to_real(plancr,C1,R1,MPI_COMM_WORLD)
            CALL io_write(1,odir,'th',ext,planio,R1)
            IF (mean.eq.1) THEN
               dump = real(cstep,kind=GP)/t
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
               DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
                  DO j = 1,n
                     DO k = 1,n
                        C1(k,j,i) = dump*M7(k,j,i)/real(n,kind=GP)**3
                     END DO
                  END DO
               END DO
               CALL fftp3d_complex_to_real(plancr,C1,R1,MPI_COMM_WORLD)
               CALL io_write(1,odir,'mean_th',ext,planio,R1)
            ENDIF
#endif
#ifdef MULTISCALAR_
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
            DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
               DO j = 1,n
                  DO k = 1,n
                     C1(k,j,i) = th1(k,j,i)/real(n,kind=GP)**3
                     C2(k,j,i) = th2(k,j,i)/real(n,kind=GP)**3
                     C3(k,j,i) = th3(k,j,i)/real(n,kind=GP)**3
                  END DO
               END DO
            END DO
            CALL fftp3d_complex_to_real(plancr,C1,R1,MPI_COMM_WORLD)
            CALL io_write(1,odir,'th1',ext,planio,R1)
            CALL fftp3d_complex_to_real(plancr,C2,R1,MPI_COMM_WORLD)
            CALL io_write(1,odir,'th2',ext,planio,R1)
            CALL fftp3d_complex_to_real(plancr,C3,R1,MPI_COMM_WORLD)
            CALL io_write(1,odir,'th3',ext,planio,R1)
            IF (mean.eq.1) THEN
               dump = real(cstep,kind=GP)/t
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
               DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
                  DO j = 1,n
                     DO k = 1,n
                        C1(k,j,i) = dump*M8 (k,j,i)/real(n,kind=GP)**3
                        C2(k,j,i) = dump*M9 (k,j,i)/real(n,kind=GP)**3
                        C3(k,j,i) = dump*M10(k,j,i)/real(n,kind=GP)**3
                     END DO
                  END DO
               END DO
               CALL fftp3d_complex_to_real(plancr,C1,R1,MPI_COMM_WORLD)
               CALL io_write(1,odir,'mean_th1',ext,planio,R1)
               CALL fftp3d_complex_to_real(plancr,C2,R1,MPI_COMM_WORLD)
               CALL io_write(1,odir,'mean_th2',ext,planio,R1)
               CALL fftp3d_complex_to_real(plancr,C3,R1,MPI_COMM_WORLD)
               CALL io_write(1,odir,'mean_th3',ext,planio,R1)
            ENDIF
#endif
#ifdef MAGFIELD_
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
            DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
               DO j = 1,n
                  DO k = 1,n
                     C1(k,j,i) = ax(k,j,i)/real(n,kind=GP)**3
                     C2(k,j,i) = ay(k,j,i)/real(n,kind=GP)**3
                     C3(k,j,i) = az(k,j,i)/real(n,kind=GP)**3
                  END DO
               END DO
            END DO
            IF (outs.ge.1) THEN
               CALL rotor3(C2,C3,C4,1)
               CALL rotor3(C1,C3,C5,2)
               CALL rotor3(C1,C2,C6,3)
               CALL fftp3d_complex_to_real(plancr,C4,R1,MPI_COMM_WORLD)
               CALL fftp3d_complex_to_real(plancr,C5,R2,MPI_COMM_WORLD)
               CALL fftp3d_complex_to_real(plancr,C6,R3,MPI_COMM_WORLD)
               CALL io_write(1,odir,'bx',ext,planio,R1)
               CALL io_write(1,odir,'by',ext,planio,R2)
               CALL io_write(1,odir,'bz',ext,planio,R3)
            ENDIF
            IF (outs.eq.2) THEN
               CALL laplak3(C1,C4)
               CALL laplak3(C2,C5)
               CALL laplak3(C3,C6)
               CALL fftp3d_complex_to_real(plancr,C4,R1,MPI_COMM_WORLD)
               CALL fftp3d_complex_to_real(plancr,C5,R2,MPI_COMM_WORLD)
               CALL fftp3d_complex_to_real(plancr,C6,R3,MPI_COMM_WORLD)
               CALL io_write(1,odir,'jx',ext,planio,R1)
               CALL io_write(1,odir,'jy',ext,planio,R2)
               CALL io_write(1,odir,'jz',ext,planio,R3)
            ENDIF
            CALL fftp3d_complex_to_real(plancr,C1,R1,MPI_COMM_WORLD)
            CALL fftp3d_complex_to_real(plancr,C2,R2,MPI_COMM_WORLD)
            CALL fftp3d_complex_to_real(plancr,C3,R3,MPI_COMM_WORLD)
            CALL io_write(1,odir,'ax',ext,planio,R1)
            CALL io_write(1,odir,'ay',ext,planio,R2)
            CALL io_write(1,odir,'az',ext,planio,R3)
            IF (mean.eq.1) THEN
               dump = real(cstep,kind=GP)/t
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
               DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
                  DO j = 1,n
                     DO k = 1,n
                        C1(k,j,i) = dump*M4(k,j,i)/real(n,kind=GP)**3
                        C2(k,j,i) = dump*M5(k,j,i)/real(n,kind=GP)**3
                        C3(k,j,i) = dump*M6(k,j,i)/real(n,kind=GP)**3
                     END DO
                  END DO
               END DO
               CALL fftp3d_complex_to_real(plancr,C1,R1,MPI_COMM_WORLD)
               CALL fftp3d_complex_to_real(plancr,C2,R2,MPI_COMM_WORLD)
               CALL fftp3d_complex_to_real(plancr,C3,R3,MPI_COMM_WORLD)
               CALL io_write(1,odir,'mean_bx',ext,planio,R1)
               CALL io_write(1,odir,'mean_by',ext,planio,R2)
               CALL io_write(1,odir,'mean_bz',ext,planio,R3)
            ENDIF
#endif
         ENDIF

#ifdef LAGPART_
         IF ( dolag.GT.0 ) THEN
           IF ((timep.eq.pstep).and.(bench.eq.0)) THEN
             INCLUDE 'hd_lagpartout.f90'
#if defined(SCALAR_)
             INCLUDE 'scalar_lagpartout.f90'
#endif
#if defined(MULTISCALAR_)
             INCLUDE 'mscalar_lagpartout.f90'
#endif
#if defined(MAGFIELD_)
             INCLUDE 'mhd_lagpartout.f90'
#endif
           ENDIF
         ENDIF
#endif

! Every 'cstep' steps, generates external files 
! with global quantities. If mean=1 also updates 
! the mean fields.

         IF ((timec.eq.cstep).and.(bench.eq.0)) THEN
            timec = 0
#ifdef HD_SOL
            INCLUDE 'hd_global.f90'
#endif
#ifdef PHD_SOL
            INCLUDE 'phd_global.f90'
#endif
#ifdef MPHD_SOL
            INCLUDE 'mphd_global.f90'
#endif
#ifdef MHD_SOL
            INCLUDE 'mhd_global.f90'
#endif
#ifdef MHDB_SOL
            INCLUDE 'mhd_global.f90'
#endif
#ifdef HMHD_SOL
            INCLUDE 'hmhd_global.f90'
#endif
#ifdef HMHDB_SOL
            INCLUDE 'hmhd_global.f90'
#endif
#ifdef ROTH_SOL
            INCLUDE 'hd_global.f90'
#endif
#ifdef PROTH_SOL
            INCLUDE 'phd_global.f90'
#endif
#if defined(BOUSS_SOL) 
            INCLUDE 'bouss_global.f90'
#endif
#if defined(ROTBOUSS_SOL)
            INCLUDE 'rotbouss_global.f90'
#endif
#ifdef LAHD_SOL
            INCLUDE 'lahd_global.f90'
#endif
#ifdef CAHD_SOL
            INCLUDE 'lahd_global.f90'
#endif
#ifdef LHD_SOL
            INCLUDE 'hd_global.f90'
#endif
#ifdef LAMHD_SOL
            INCLUDE 'lamhd_global.f90'
#endif
#ifdef EDQNMHD_SOL
            INCLUDE 'hd_global.f90'
#endif
#ifdef EDQNMROTH_SOL
            INCLUDE 'hd_global.f90'
#endif
            IF (mean.eq.1) THEN
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
               DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
                  DO j = 1,n
                     DO k = 1,n
                        M1(k,j,i) = M1(k,j,i)+vx(k,j,i)
                        M2(k,j,i) = M2(k,j,i)+vy(k,j,i)
                        M3(k,j,i) = M3(k,j,i)+vz(k,j,i)
                     END DO
                  END DO
               END DO
#ifdef SCALAR_
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
               DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
                  DO j = 1,n
                     DO k = 1,n
                        M7(k,j,i) = M7(k,j,i)+th(k,j,i)
                     END DO
                  END DO
               END DO
#endif
#ifdef MULTISCALAR_
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
               DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
                  DO j = 1,n
                     DO k = 1,n
                        M8 (k,j,i) = M8 (k,j,i)+th1(k,j,i)
                        M9 (k,j,i) = M9 (k,j,i)+th2(k,j,i)
                        M10(k,j,i) = M10(k,j,i)+th3(k,j,i)
                     END DO
                  END DO
               END DO
#endif
#ifdef MAGFIELD_
               CALL rotor3(ay,az,C1,1)
               CALL rotor3(ax,az,C2,2)
               CALL rotor3(ax,ay,C3,3)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
               DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
                  DO j = 1,n
                     DO k = 1,n
                        M4(k,j,i) = M4(k,j,i)+vx(k,j,i)
                        M5(k,j,i) = M5(k,j,i)+vy(k,j,i)
                        M6(k,j,i) = M6(k,j,i)+vz(k,j,i)
                     END DO
                  END DO
               END DO
#endif
            ENDIF
         ENDIF

! Every 'sstep' steps, generates external files 
! with the power spectrum.

         IF ((times.eq.sstep).and.(bench.eq.0)) THEN
            times = 0
            sind = sind+1
            WRITE(ext, fmtext) sind
#ifdef HD_SOL
            INCLUDE 'hd_spectrum.f90'
#endif
#ifdef PHD_SOL
            INCLUDE 'phd_spectrum.f90'
#endif
#ifdef MPHD_SOL
            INCLUDE 'mphd_spectrum.f90'
#endif
#ifdef MHD_SOL
            INCLUDE 'mhd_spectrum.f90'
#endif
#ifdef MHDB_SOL
            INCLUDE 'mhdb_spectrum.f90'
#endif
#ifdef HMHD_SOL
            INCLUDE 'hmhd_spectrum.f90'
#endif
#ifdef HMHDB_SOL
            INCLUDE 'hmhdb_spectrum.f90'
#endif
#ifdef ROTH_SOL
            INCLUDE 'roth_spectrum.f90'
#endif
#ifdef PROTH_SOL
            INCLUDE 'proth_spectrum.f90'
#endif
#ifdef BOUSS_SOL
            INCLUDE 'bouss_spectrum.f90'
#endif
#ifdef ROTBOUSS_SOL
            INCLUDE 'rotbouss_spectrum.f90'
#endif
#ifdef LAHD_SOL
            INCLUDE 'lahd_spectrum.f90'
#endif
#ifdef CAHD_SOL
            INCLUDE 'lahd_spectrum.f90'
#endif
#ifdef LHD_SOL
            INCLUDE 'hd_spectrum.f90'
#endif
#ifdef LAMHD_SOL
            INCLUDE 'lamhd_spectrum.f90'
#endif
#ifdef EDQNMHD_SOL
            INCLUDE 'edqnmhd_spectrum.f90'
#endif
#ifdef EDQNMROTH_SOL
            INCLUDE 'edqnmroth_spectrum.f90'
#endif
         ENDIF

! Runge-Kutta step 1
! Copies the fields into auxiliary arrays

!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,n
         DO k = 1,n

#ifdef HD_SOL
         INCLUDE 'hd_rkstep1.f90'
#endif
#ifdef PHD_SOL
         INCLUDE 'phd_rkstep1.f90'
#endif
#ifdef MPHD_SOL
         INCLUDE 'mphd_rkstep1.f90'
#endif
#ifdef MHD_SOL
         INCLUDE 'mhd_rkstep1.f90'
#endif
#ifdef MHDB_SOL
         INCLUDE 'mhd_rkstep1.f90'
#endif
#ifdef HMHD_SOL
         INCLUDE 'mhd_rkstep1.f90'
#endif
#ifdef HMHDB_SOL
         INCLUDE 'mhd_rkstep1.f90'
#endif
#ifdef ROTH_SOL
         INCLUDE 'hd_rkstep1.f90'
#endif
#ifdef PROTH_SOL
         INCLUDE 'phd_rkstep1.f90'
#endif
#if defined(BOUSS_SOL) || defined(ROTBOUSS_SOL)
         INCLUDE 'phd_rkstep1.f90'
#endif
#ifdef LAHD_SOL
         INCLUDE 'hd_rkstep1.f90'
#endif
#ifdef CAHD_SOL
         INCLUDE 'hd_rkstep1.f90'
#endif
#ifdef LHD_SOL
         INCLUDE 'hd_rkstep1.f90'
#endif
#ifdef LAMHD_SOL
         INCLUDE 'mhd_rkstep1.f90'
#endif
#ifdef EDQNMHD_SOL
         INCLUDE 'hd_rkstep1.f90'
#endif
#ifdef EDQNMROTH_SOL
         INCLUDE 'hd_rkstep1.f90'
#endif

         END DO
         END DO
         END DO

#ifdef LAGPART_
         IF ( dolag.GT.0 ) THEN
           CALL lagpart%SetStep()
         ENDIF
#endif

! Runge-Kutta step 2
! Evolves the system in time

         DO o = ord,1,-1
#ifdef HD_SOL
         INCLUDE 'hd_rkstep2.f90'
#endif
#ifdef PHD_SOL
         INCLUDE 'phd_rkstep2.f90'
#endif
#ifdef MPHD_SOL
         INCLUDE 'mphd_rkstep2.f90'
#endif
#ifdef MHD_SOL
         INCLUDE 'mhd_rkstep2.f90'
#endif
#ifdef MHDB_SOL
         INCLUDE 'mhdb_rkstep2.f90'
#endif
#ifdef HMHD_SOL
         INCLUDE 'hmhd_rkstep2.f90'
#endif
#ifdef HMHDB_SOL
         INCLUDE 'hmhdb_rkstep2.f90'
#endif
#ifdef ROTH_SOL
         INCLUDE 'roth_rkstep2.f90'
#endif
#ifdef PROTH_SOL
         INCLUDE 'proth_rkstep2.f90'
#endif
#ifdef BOUSS_SOL
         INCLUDE 'bouss_rkstep2.f90'
#endif
#ifdef ROTBOUSS_SOL
         INCLUDE 'rotbouss_rkstep2.f90'
#endif
#ifdef LAHD_SOL
         INCLUDE 'lahd_rkstep2.f90'
#endif
#ifdef CAHD_SOL
         INCLUDE 'lahd_rkstep2.f90'
#endif
#ifdef LHD_SOL
         INCLUDE 'lhd_rkstep2.f90'
#endif
#ifdef LAMHD_SOL
         INCLUDE 'lamhd_rkstep2.f90'
#endif
#ifdef EDQNMHD_SOL
         INCLUDE 'edqnmhd_rkstep2.f90'
#endif
#ifdef EDQNMROTH_SOL
         INCLUDE 'edqnmroth_rkstep2.f90'
#endif

#ifdef LAGPART_
      IF ( dolag.GT.0 ) THEN
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
        DO j = 1,n
          DO k = 1,n
            C7(k,j,i) = vx(k,j,i)/real(n,kind=GP)**3
          END DO
        END DO
      END DO
      CALL fftp3d_complex_to_real(plancr,C7,R1,MPI_COMM_WORLD)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
        DO j = 1,n
          DO k = 1,n
            C7(k,j,i) = vy(k,j,i)/real(n,kind=GP)**3
          END DO
        END DO
      END DO
      CALL fftp3d_complex_to_real(plancr,C7,R2,MPI_COMM_WORLD)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
        DO j = 1,n
          DO k = 1,n
            C7(k,j,i) = vz(k,j,i)/real(n,kind=GP)**3
          END DO
        END DO
      END DO
      CALL fftp3d_complex_to_real(plancr,C7,R3,MPI_COMM_WORLD)
      CALL lagpart%Step(R1,R2,R3,dt,1.0_GP/real(o,kind=GP),R4,R5)
      CALL lagpart%FinalStep()
      ENDIF
#endif

         END DO

         timet = timet+1
         times = times+1
         timec = timec+1
         timef = timef+1
         timep = timep+1

      END DO RK


!
! End of Runge-Kutta

! Computes the benchmark

      IF (bench.gt.0) THEN
         CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
         CALL GTStop(ihcpu1)
         CALL GTStop(ihomp1)
         CALL GTStop(ihwtm1)
#ifdef LAGPART_
        IF ( dolag.GT.0 ) THEN
          rbal = rbal + lagpart%GetLoadBal()
        ENDIF
#endif
         IF (myrank.eq.0) THEN
            OPEN(1,file='benchmark.txt',position='append')
#if defined(DEF_GHOST_CUDA_)
            WRITE(1,*) n,(step-ini+1),nprocs,nth, &
                       GTGetTime(ihcpu1)/(step-ini+1), &
                       GTGetTime(ihomp1)/(step-ini+1), &
                       GTGetTime(ihwtm1)/(step-ini+1), &
                       ffttime/(step-ini+1), tratime/(step-ini+1), &
                       comtime/(step-ini+1), memtime/(step-ini+1)
#else
            WRITE(1,*) n,(step-ini+1),nprocs,nth, &
                       GTGetTime(ihcpu1)/(step-ini+1), &
                       GTGetTime(ihomp1)/(step-ini+1), &
                       GTGetTime(ihwtm1)/(step-ini+1), &
                       ffttime/(step-ini+1), tratime/(step-ini+1), &
                       comtime/(step-ini+1)

#endif
            IF (bench.eq.2) THEN
               WRITE(1,*) 'FFTW: Create_plan = ', &
                       GTGetTime(ihcpu2)/(step-ini+1), &
                       GTGetTime(ihomp2)/(step-ini+1), &
                       GTGetTime(ihwtm2)/(step-ini+1)
            ENDIF
            CLOSE(1)
#if defined(LAGPART_)
            IF ( dolag.GT.0 ) THEN
              OPEN(2,file='gpbenchmark.txt',position='append')
                WRITE(2,*) n,maxparts,rbal/(step-ini+1),(step-ini+1),nprocs,nth, &
                           lagpart%GetTime   (GPTIME_STEP)/(step-ini+1), &
                           lagpart%GetTime   (GPTIME_COMM)/(step-ini+1), &
                           lagpart%GetTime (GPTIME_SPLINE)/(step-ini+1), &
                           lagpart%GetTime (GPTIME_TRANSP)/(step-ini+1), &
                           lagpart%GetTime (GPTIME_DATAEX)/(step-ini+1), &
                           lagpart%GetTime (GPTIME_INTERP)/(step-ini+1), &
                           lagpart%GetTime(GPTIME_PUPDATE)/(step-ini+1)
              CLOSE(2)
            ENDIF
#endif
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
      DEALLOCATE( vx,vy,vz,fx,fy,fz )
      DEALLOCATE( C1,C2,C3,C4,C5,C6,C7,C8 )
      DEALLOCATE( ka,ka2 )
      IF (mean.eq.1) DEALLOCATE( M1,M2,M3 )
      IF (rand.eq.2) DEALLOCATE( Faux1, Faux2 )
#ifdef SCALAR_
      DEALLOCATE( th,fs )
      DEALLOCATE( C20 )
      IF (mean.eq.1) DEALLOCATE( M7 )
#endif
#ifdef MULTISCALAR_
      DEALLOCATE( th1,fs1,th2,fs2,th3,fs3 )
      DEALLOCATE( C21,C22,C23,C24 )
      IF (mean.eq.1) DEALLOCATE( M8,M9,M10 )
#endif
#ifdef MAGFIELD_
      DEALLOCATE( ax,ay,az,mx,my,mz )
      DEALLOCATE( C9,C10,C11,C12,C13,C14,C15,C16,C17 )
      IF (mean.eq.1) DEALLOCATE( M4,M5,M6 )
#endif
#ifdef HALLTERM_
      DEALLOCATE( C18 )
#endif
#ifdef EDQNM_
      DEALLOCATE( C19 )
      DEALLOCATE( tepq,thpq,tve,tvh,Eext,Hext )
#endif 

      END PROGRAM MAIN3D
