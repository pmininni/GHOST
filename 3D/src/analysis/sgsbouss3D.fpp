!=================================================================
      PROGRAM MAIN3D
!=================================================================
! SGSBOUSS3D: Reads binary GHOST data at full resolution,
! truncates, and computes the LES version of the RHS of the 
! equations: the filtered operators acting on the
! filtered data. This should be built with the truncated
! grid resolution. The resolution provided in the
! input file should be the full grid resolution of the
! input data. Thus, the 'truncated' grid here is actually
! the full-res grid.

! This utility is meant to be used togther with 
! LESBOUS3D, and the results will be subtracted from the 
! time derivative estimates from that utility to compute 
! the full truth solution for the SGS.
!
! Notation: index 'i' is 'x'
!           index 'j' is 'y'
!           index 'k' is 'z'
!
! Conditional compilation options:
!           HD_SOL         builds the hydrodynamic (HD) solver
!           PHD_SOL        builds the HD solver with passive scalar
!           MPHD_SOL       builds the HD solver with multiple scalars
!           MHD_SOL        builds the MHD solver
!           MHDB_SOL       builds the MHD solver with uniform B_0
!           RMHDB_SOL      builds the MHD solver with B_0 and rotation
!           HMHD_SOL       builds the Hall-MHD solver
!           HMHDB_SOL      builds the HMHD solver with uniform B_0
!           COMPRHD_SOL    builds the compressible HD solver
!           CMHD_SOL       builds the compressible MHD solver
!           CMHDB_SOL      builds the compressible MHD solver with B_0
!           CHMHD_SOL      builds the compressible Hall-MHD solver
!           CHMHDB_SOL     builds the compressible HMHD solver with B_0
!           ROTH_SOL       builds the HD solver in a rotating frame
!           PROTH_SOL      builds the ROTH solver with passive scalar
!           MPROTH_SOL     builds the ROTH solver with multi-scalar
!           BOUSS_SOL      builds the BOUSS solver
!           ROTBOUSS_SOL   builds the BOUSS solver in a rotating frame
!           ROTBOUMHDB_SOL builds BOUSS and MHD with rotation and B_0
!           MPROTBOUSS_SOL builds the BOUSS eq, rotating, multi-scalar
!           GPE_SOL        builds the Gross-Pitaevskii Equation solver
!           ARGL_SOL       builds the Advective Real Ginzburg Landau
!           RGPE_SOL       builds the rotating GPE solver
!           RARGL_SOL      builds the rotating ARGL solver
!           LAHD_SOL       builds the Lagrangian-averaged HD solver
!           CAHD_SOL       builds the Clark-alpha HD solver
!           LHD_SOL        builds the Leray HD solver
!           LAMHD_SOL      builds the Lagrangian-averaged MHD solver
!           EDQNMHD_SOL    builds the EDQNM HD solver
!           EDQNMROTH_SOL  builds the EDQNM ROTH solver
!
! 2003 Pablo D. Mininni.
!      Department of Physics,
!      Facultad de Ciencias Exactas y Naturales.
!      Universidad de Buenos Aires.
!      e-mail: mininni@df.uba.ar
!
! 15 Feb 2007: Main program for all solvers (HD/MHD/HMHD)
! 21 Feb 2007: POSIX and MPI-IO support
! 10 Mar 2007: FFTW-2.x and FFTW-3.x support
! 25 Aug 2009: Hybrid MPI/OpenMP support (D. Rosenberg & P. Mininni)
! 30 Aug 2009: SINGLE/DOUBLE precision (D. Rosenberg & P. Mininni)
! 10 Feb 2011: Hybrid MPI/OpenMP/CUDA support (D. Rosenberg)
! 21 Nov 2016: Anisotropic boxes (A. Alexakis & P. Mininni)
!
! References:
! Mininni PD, Rosenberg DL, Reddy R, Pouquet A.; P.Comp.37, 123 (2011)
! Mininni PD, Gomez DO, Mahajan SM; Astrophys. J. 619, 1019 (2005)
! Gomez DO, Mininni PD, Dmitruk P; Phys. Scripta T116, 123 (2005)
! Gomez DO, Mininni PD, Dmitruk P; Adv. Sp. Res. 35, 899 (2005)
!=================================================================

!
! Definitions for conditional compilation
#include "lesbouss3D.h"

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
      USE gutils
      USE gtrunc

#ifdef DNS_
      USE dns
#endif
#ifdef HALLTERM_
      USE hall
#endif
#ifdef WAVEFUNCTION_
      USE newtmod
      USE hbar
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
      IMPLICIT NONE

!
! Arrays for the fields and external forcings

#if defined(VELOC_) || defined(ADVECT_)
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: vx,vy,vz
#endif
#ifdef VELOC_
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: fx,fy,fz
#endif
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
#ifdef WAVEFUNCTION_
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: zre,zim
#endif
#ifdef QFORCE_
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: fre,fim
#endif

!
! Temporal data storage arrays

      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:), TARGET :: C1,C2
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:), TARGET :: C3,C4
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:), TARGET :: C5,C6
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:), TARGET :: C7,C8
#ifdef VELOC_
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: M1,M2,M3
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
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:)     :: tepq,thpq,tve,tvh
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:)     :: Eold,Hold
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:)     :: Eext,Hext
      REAL(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:)    :: Eden,Hden
#endif
#ifdef SCALAR_
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: C20
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: M7
#endif
#ifdef MULTISCALAR_
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: C21,C22,C23,C24
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: M8,M9,M10
#endif
#ifdef COMPR_AUX_ARR_
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: C25,C26,C27
#endif
#ifdef TRAP_
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: C28,C29
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: C30
#endif

#ifdef WAVEFUNCTION_
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:)     :: iold,qold,kold,cold
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:)     :: inew,qnew,knew,cnew
#endif

      REAL(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:), TARGET    :: R1,R2,R3
#ifdef VELOC_
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: fxold,fyold,fzold
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: fxnew,fynew,fznew
      REAL(KIND=GP), ALLOCATABLE, DIMENSION (:)        :: Faux1,Faux2
#endif
#ifdef MAGFIELD_
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: mxold,myold,mzold
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: mxnew,mynew,mznew
#endif
#ifdef ADVECT_
      REAL(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:)    :: vsq
#endif
#ifdef TRAP_
      REAL(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:)    :: Vtrap
      REAL(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:)    :: Vlinx,Vliny
#endif

      REAL(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:)    :: R4,R5,R6
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
      REAL(KIND=GP)    :: phase,ampl,cort
      REAL(KIND=GP)    :: fparam0,fparam1,fparam2,fparam3,fparam4
      REAL(KIND=GP)    :: fparam5,fparam6,fparam7,fparam8,fparam9
      REAL(KIND=GP)    :: vparam0,vparam1,vparam2,vparam3,vparam4
      REAL(KIND=GP)    :: vparam5,vparam6,vparam7,vparam8,vparam9
#ifdef VELOC_
      REAL(KIND=GP)    :: bvfreq,xmom,xtemp
#endif
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
#ifdef MULTISCALAR_
      DOUBLE PRECISION :: tmp1,tmq1,tmp2,tmq2,tmp3,tmq3
      REAL(KIND=GP)    :: kappa1,kappa2,kappa3
      REAL(KIND=GP)    :: rmp1,rmq1,rmp2,rmq2,rmp3,rmq3
      REAL(KIND=GP)    :: cc10,ss10,cc20,ss20,cc30,ss30
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
#ifdef COMPRESSIBLE_
      REAL(KIND=GP)    :: smach, gam1, cp1, nu2
#endif
#ifdef CMHD_
      REAL(KIND=GP)    :: amach, cp2
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
      REAL(KIND=GP)    :: omegax,omegay,omegaz
#endif
#ifdef WAVEFUNCTION_
      REAL(KIND=GP)    :: cspeed,lambda,rho0,kttherm,V0
      REAL(KIND=GP)    :: zparam0,zparam1,zparam2,zparam3,zparam4
      REAL(KIND=GP)    :: zparam5,zparam6,zparam7,zparam8,zparam9
#endif

      INTEGER :: idevice, iret, ncuda, ngcuda, ppn
      INTEGER :: ini,step
      INTEGER :: tstep,cstep
      INTEGER :: sstep,fstep
      INTEGER :: bench,trans
      INTEGER :: outs,mean
      INTEGER :: seed,rand
      INTEGER :: anis
      INTEGER :: mult=1
      INTEGER :: t,o
      INTEGER :: i,j,k
      INTEGER :: ki,kj,kk
      INTEGER :: pind,tind,sind
      INTEGER :: timet,timec
      INTEGER :: times,timef
      INTEGER :: timep,pstep,lgmult=1
      INTEGER :: ihcpu1,ihcpu2
      INTEGER :: ihomp1,ihomp2
      INTEGER :: ihwtm1,ihwtm2
#if defined(SCALAR_) || defined(MULTISCALAR_)
      INTEGER :: injt,injtm
      INTEGER :: creset
#endif
#ifdef MAGFIELD_
      INTEGER :: dyna
      INTEGER :: corr
#endif
#ifdef WAVEFUNCTION_
      INTEGER :: cflow
#endif

#if defined(DEF_GHOST_CUDA_)
      TYPE(cudaDevicePropG) :: devprop
#endif
      TYPE(IOPLAN)          :: planio
      CHARACTER(len=1024)   :: odir,idir
      LOGICAL               :: bbenchexist

      ! App data:
      COMPLEX(KIND=GP), ALLOCATABLE, TARGET, DIMENSION (:,:,:) :: CT1, CT2
      REAL(KIND=GP)   , ALLOCATABLE, TARGET, DIMENSION (:,:,:) :: RT1
      INTEGER             :: istat(4096), npkeep, nstat
      INTEGER             :: commtrunc, grouptrunc, n(3), nt(3)
      CHARACTER(len=1024) :: iidir, sparam
      CHARACTER(len=64)   :: ext1
      CHARACTER(len=1024) :: sstat
      TYPE(IOPLAN)        :: planiot
      TYPE(FFTPLAN)       :: planrct

      commtrunc  = MPI_COMM_NULL
      grouptrunc = MPI_GROUP_NULL

!
! Namelists for the input files

      NAMELIST / status / idir,odir,stat,mult,bench,outs,mean,trans,iswap
      NAMELIST / parameter / dt,step,tstep,sstep,cstep,rand,cort,seed
#ifdef DEF_ARBSIZE_
      NAMELIST / boxparams / Lx,Ly,Lz,Dkk
#endif
#if defined(VELOC_) || defined(ADVECT_)
      NAMELIST / velocity / f0,u0,kdn,kup,nu,fparam0,fparam1,fparam2
      NAMELIST / velocity / fparam3,fparam4,fparam5,fparam6,fparam7
      NAMELIST / velocity / fparam8,fparam9,vparam0,vparam1,vparam2
      NAMELIST / velocity / vparam3,vparam4,vparam5,vparam6,vparam7
      NAMELIST / velocity / vparam8,vparam9
#endif
#ifdef SCALAR_
      NAMELIST / scalar / c0,s0,skdn,skup,kappa,cparam0,cparam1
      NAMELIST / scalar / cparam2,cparam3,cparam4,cparam5,cparam6
      NAMELIST / scalar / cparam7,cparam8,cparam9,sparam0,sparam1
      NAMELIST / scalar / sparam2,sparam3,sparam4,sparam5,sparam6
      NAMELIST / scalar / sparam7,sparam8,sparam9
#endif
#ifdef MULTISCALAR_
      NAMELIST / mscalar / cc10,ss10,cc20,ss20,cc30,ss30
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
      NAMELIST / inject / injt,injtm,creset
#endif
#ifdef COMPRESSIBLE_
      NAMELIST / compressible / smach, gam1, nu2
#endif
#ifdef CMHD_
      NAMELIST / cmhdb / amach
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
      NAMELIST / rotation / omegax,omegay,omegaz
#endif
#ifdef BOUSSINESQ_
      NAMELIST / boussinesq / bvfreq,xmom,xtemp
#endif
#ifdef WAVEFUNCTION_
      NAMELIST / wavefunction / cspeed,lambda,rho0,kttherm,V0
      NAMELIST / wavefunction / cflow,iter_max_newt,iter_max_bicg
      NAMELIST / wavefunction / cflow_newt,dt_newt,tol_newt,tolbicg_rel
      NAMELIST / wavefunction / zparam0,zparam1,zparam2,zparam3,zparam4
      NAMELIST / wavefunction / zparam5,zparam6,zparam7,zparam8,zparam9
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

      ! App NAMELIST
      NAMELIST / regrid / idir, odir, sstat, iswap, nxt, nyt, nzt


!
! Initialization

! Initializes the MPI and I/O libraries
      CALL MPI_INIT_THREAD(MPI_THREAD_FUNNELED,provided,ierr)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)

! Initialization of offloading to GPUs using OpenMP (this is independent
! of CUDA initialization in systems with NVIDIA GPUs). GHOST
! assumes the number of MPI jobs in each node is equal to the
! number of GPUs available in the node. The user must ensure this
! condition is fulfilled.
#if defined(DO_HYBRIDoffl)
      CALL init_offload(myrank,numdev,hostdev,targetdev)
#endif

! NOTE: On systems with a single GPU per node (e.g., Titan)
!       we remove the following block. But on systems with 
!       multiple devices per node, this will have to be 
!       considered carefully, and possibly re-def'd:
#if defined(DEF_GHOST_CUDA_)
#if defined(CUDA_BIND_LINUX_)
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
#endif
#if defined(CUDA_BIND_DEVICE_)
! Initializes CUDA by selecting device. The list of devices can
! be changed by modifying the env. variable CUDA_VISIBLE_DEVICES:
     iret    = cudaGetDeviceCount(ncuda)
     idevice = mod(myrank,ncuda)
     iret    = cudaSetDevice(idevice);
     IF ( iret .EQ. cudaErrorInvalidDevice ) THEN
       WRITE(*,*)'MAIN: Invalid CUDA device selected: ', &
       idevice, '; myrank=',myrank, '; NUM_CUDA_DEV=',ncuda
       STOP
     ENDIF
     CALL cudaGetDeviceProperties(devprop,idevice)
     IF ( nstreams .GT. 1 .AND. devprop%deviceOverlap .EQ. 0 ) THEN
       WRITE(*,*)'MAIN: Async transfer and computation overlap not supported!'
!      STOP
     ENDIF
     iret = cudaGetDevice(idevice)
#endif
#endif

     CALL range(1,nx/2+1,nprocs,myrank,ista,iend)
     CALL range(1,nz,nprocs,myrank,ksta,kend)
     CALL io_init_comm(MPI_COMM_WORLD,(/nx,ny,nz/),ksta,kend,planio)

!
! Allocates memory for distributed blocks

      ALLOCATE( C1(nz,ny,ista:iend),  C2(nz,ny,ista:iend) )
      ALLOCATE( C3(nz,ny,ista:iend),  C4(nz,ny,ista:iend) )
      ALLOCATE( C5(nz,ny,ista:iend),  C6(nz,ny,ista:iend) )
      ALLOCATE( C7(nz,ny,ista:iend),  C8(nz,ny,ista:iend) )
#if defined(VELOC_) || defined(ADVECT_)
      ALLOCATE( vx(nz,ny,ista:iend) )
      ALLOCATE( vy(nz,ny,ista:iend) )
      ALLOCATE( vz(nz,ny,ista:iend) )
#endif
#ifdef VELOC_
      ALLOCATE( fx(nz,ny,ista:iend) )
      ALLOCATE( fy(nz,ny,ista:iend) )
      ALLOCATE( fz(nz,ny,ista:iend) )
#endif
#ifdef MAGFIELD_
      ALLOCATE( C9 (nz,ny,ista:iend), C10(nz,ny,ista:iend) )
      ALLOCATE( C11(nz,ny,ista:iend), C12(nz,ny,ista:iend) )
      ALLOCATE( C13(nz,ny,ista:iend), C14(nz,ny,ista:iend) )
      ALLOCATE( C15(nz,ny,ista:iend), C16(nz,ny,ista:iend) )
      ALLOCATE( C17(nz,ny,ista:iend) )
      ALLOCATE( ax(nz,ny,ista:iend) )
      ALLOCATE( ay(nz,ny,ista:iend) )
      ALLOCATE( az(nz,ny,ista:iend) )
      ALLOCATE( mx(nz,ny,ista:iend) )
      ALLOCATE( my(nz,ny,ista:iend) )
      ALLOCATE( mz(nz,ny,ista:iend) )
#endif
#ifdef HALLTERM_
      ALLOCATE( C18(nz,ny,ista:iend) )
#endif
#ifdef EDQNM_
      n = nx ! EDQNM solvers only work in cubic boxes
      ALLOCATE( C19(nz,ny,ista:iend) )
#endif
#ifdef SCALAR_
      ALLOCATE( C20(nz,ny,ista:iend) )
      ALLOCATE( th (nz,ny,ista:iend) )
      ALLOCATE( fs (nz,ny,ista:iend) )
#endif
#ifdef MULTISCALAR_
      ALLOCATE( C21(nz,ny,ista:iend), C22(nz,ny,ista:iend) )
      ALLOCATE( C23(nz,ny,ista:iend), C24(nz,ny,ista:iend) )
      ALLOCATE( th1(nz,ny,ista:iend) )
      ALLOCATE( th2(nz,ny,ista:iend) )
      ALLOCATE( th3(nz,ny,ista:iend) )
      ALLOCATE( fs1(nz,ny,ista:iend) )
      ALLOCATE( fs2(nz,ny,ista:iend) )
      ALLOCATE( fs3(nz,ny,ista:iend) )
#endif
#ifdef COMPR_AUX_ARR_
      ALLOCATE( C25(nz,ny,ista:iend) )
      ALLOCATE( C26(nz,ny,ista:iend) )
      ALLOCATE( C27(nz,ny,ista:iend) )
#endif
#ifdef TRAP_
      ALLOCATE( C28(nz,ny,ista:iend), C29(nz,ny,ista:iend) )
      ALLOCATE( C30(nz,ny,ista:iend) )
#endif
#ifdef WAVEFUNCTION_
      ALLOCATE( zre(nz,ny,ista:iend), zim(nz,ny,ista:iend) )
#endif
#ifdef QFORCE_
      ALLOCATE( fre(nz,ny,ista:iend), fim(nz,ny,ista:iend) )
#endif

      ALLOCATE( kx(nx), ky(ny), kz(nz) )
      ALLOCATE( kn2(nz,ny,ista:iend) )
#ifdef DEF_ARBSIZE_
      anis = 1
      ALLOCATE( kk2(nz,ny,ista:iend) )
#else
      IF ((nx.ne.ny).or.(ny.ne.nz)) THEN
         anis = 1
         ALLOCATE( kk2(nz,ny,ista:iend) )
      ELSE
         anis = 0
         kk2 => kn2
      ENDIF
#endif

      ALLOCATE( R1(nx,ny,ksta:kend) )
      ALLOCATE( R2(nx,ny,ksta:kend) )
      ALLOCATE( R3(nx,ny,ksta:kend) )
#ifdef ADVECT_
      ALLOCATE( vsq(nx,ny,ksta:kend) )
#endif
#ifdef TRAP_
      ALLOCATE( Vtrap(nx,ny,ksta:kend) )
      ALLOCATE( Vlinx(nx,ny,ksta:kend), Vliny(nx,ny,ksta:kend) )
#endif

      ALLOCATE( R4(nx,ny,ksta:kend) )
      ALLOCATE( R5(nx,ny,ksta:kend) )
      ALLOCATE( R6(nx,ny,ksta:kend) )

#ifdef EDQNM_
      ALLOCATE( Eden(nz,ny,ista:iend) )
      ALLOCATE( Hden(nz,ny,ista:iend) )
      ALLOCATE( tepq(n/2+1) )
      ALLOCATE( thpq(n/2+1) )
      ALLOCATE( tve (n/2+1) )
      ALLOCATE( tvh (n/2+1) )
      ALLOCATE( Eold(n/2+1) )
      ALLOCATE( Hold(n/2+1) )
      ALLOCATE( Eext(3*(n/2+1)) )
      ALLOCATE( Hext(3*(n/2+1)) )
#endif

! Reads from the external app file 'sgsml.inp' the 
! parameters that will be used to compute the transfer
!     idir   : directory for unformatted input (field components)
!     odir   : directory for unformatted output (truncated data)
!     sstat  : file list of time indices, separated by ';'
!     iswap  : do endian swap on input?
!     nxt    : truncated linear size in x direction
!     nyt    : truncated linear size in y direction
!     nzt    : truncated linear size in z direction
      idir   = '.'
      odir   = '.'
      sstat  = ''
      iswap  = 0
      nxt    = 0
      nyt    = 0
      nzt    = 0

      IF (myrank.eq.0) THEN
         OPEN(1,file='sgsml.inp',status='unknown')
         READ(1,NML=regrid)
         CLOSE(1)
      ENDIF
      CALL MPI_BCAST(idir  ,1024,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(odir  ,1024,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(sstat ,1024,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(iswap ,1   ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(nxt   ,1   ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(nyt   ,1   ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(nzt   ,1   ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)


      ! Check input quantities:
      IF ( nxt .LT. nx .OR. nxt .LT. 1 ) THEN
        IF ( myrank .eq. 0) &
          PRINT*, 'MAIN: truncation specification incorrect; input nxt must be less than Nx'
        CALL MPI_Finalize(ierr)
        STOP
      ENDIF
      IF ( nyt .LT. ny .OR. nyt .LT. 1 ) THEN
        IF ( myrank .eq. 0) &
          PRINT*, 'MAIN: truncation specification incorrect; input nyt must be less than Ny'
        CALL MPI_Finalize(ierr)
        STOP
      ENDIF
      IF ( nzt .LT. nz .OR. nzt .LT. 1 ) THEN
        IF ( myrank .eq. 0) &
          PRINT*, 'MAIN: truncation specification incorrect; input nzt must be less than Nz'
        CALL MPI_Finalize(ierr)
        STOP
      ENDIF


!
! The following lines read the file 'input_dir/parameter.inp':

  sparam = trim(idir) // '/' // 'parameter.inp' 


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
         OPEN(1,file=sparam,status='unknown',form="formatted")
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

      write(*,*) 'main: sparam=', sparam, '  dt_0=', dt

      Lx  = 1.0_GP
      Ly  = 1.0_GP
      Lz  = 1.0_GP
      Dkx = 1.0_GP
      Dky = 1.0_GP
      Dkz = 1.0_GP
      Dkk = 0.0_GP
#ifdef DEF_ARBSIZE_
! Reads parameters to set the box size
!     Lx  : Length in x (in units of 2.pi, =1 gives a side of length 2.pi)
!     Ly  : Length in y
!     Lz  : Length in z
!     Dkk : Width of Fourier shells for 2D and 3D spectra
!           Default = min(1/Lx, 1/Ly, 1/Lz)

      IF (myrank.eq.0) THEN
         OPEN(1,file=sparam,status='unknown',form="formatted")
         READ(1,NML=boxparams)
         CLOSE(1)
      ENDIF
      CALL MPI_BCAST(Lx,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(Ly,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(Lz,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(Dkk,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      Dkx = 1.0_GP/Lx
      Dky = 1.0_GP/Ly
      Dkz = 1.0_GP/Lz
#endif
      IF (Dkk.lt.1e-5) Dkk = min(Dkx,Dky,Dkz)

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
         OPEN(1,file=sparam,status='unknown',form="formatted")
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

#ifdef BOUSSINESQ_
! Reads parameters specifically for Boussinesq solver from the
! namelist 'boussinesq' on the external file 'parameter.inp'
!     bvfreq: Brunt-Vaisala frequency (positive definite)
!     xmom  : multiplies bouyancy term in momentum equation
!     xtemp : multiplies temperature-current term in
!             temperature/density equation
      xmom  = 1.0
      xtemp = 1.0
      IF (myrank.eq.0) THEN
         OPEN(1,file=sparam,status='unknown',form="formatted")
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
! Reads general configuration flags for runs with
! a passive/active scalar from the namelist 'inject'
! on the external file 'parameter.inp'
!     injt : = 0 when stat=0 generates initial v and th (SCALAR_)
!            = 1 when stat.ne.0 imports v and generates th (SCALAR_)
!     creset: = 0: don't reset counters; 1 = reset counters

      injt   = 0
      creset = 1
      IF (myrank.eq.0) THEN
         OPEN(1,file=sparam,status='unknown',form="formatted")
         READ(1,NML=inject)
         CLOSE(1)
      ENDIF
      CALL MPI_BCAST(injt  ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(creset,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

! Reads parameters for the passive/active scalar from the
! namelist 'scalar' on the external file 'parameter.inp'
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
         OPEN(1,file=sparam,status='unknown',form="formatted")
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
! Reads general configuration flags for runs with
! a passive/active scalar from the namelist 'inject'
! on the external file 'parameter.inp'
!     injtm : = 0 when stat=0 generates initial v,th (SCALAR_), th[1-3]
!             = 1 when stat.ne.0 imports v,th and generates th[1-3]
!     creset: = 0: don't reset counters; 1 = reset counters

      injtm  = 0
      creset = 1
      IF (myrank.eq.0) THEN
         OPEN(1,file=sparam,status='unknown',form="formatted")
         READ(1,NML=inject)
         CLOSE(1)
      ENDIF
      CALL MPI_BCAST(injt  ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(injtm ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(creset,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

! Reads parameters for the passive/active scalar from the
! namelist 'mscalar' on the external file 'parameter.inp'
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
         OPEN(1,file=sparam,status='unknown',form="formatted")
         READ(1,NML=mscalar)
         CLOSE(1)
      ENDIF
      CALL MPI_BCAST(cc10   ,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(cc20   ,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(cc30   ,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(ss10   ,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(ss20   ,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(ss30   ,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
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

#ifdef COMPRESSIBLE_
! Reads parameters for the compressible runs from the
! namelist 'compressible' on the external file 'parameter.inp'
!     smach : sound Mach number
!     gam1  : gamma parameter for polytropic eq. of state
!     nu2   : second viscosity for divergence (velocity) term

      IF (myrank.eq.0) THEN
         OPEN(1,file=sparam,status='unknown',form="formatted")
         READ(1,NML=compressible)
         CLOSE(1)
      ENDIF
      CALL MPI_BCAST(smach,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(gam1,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(nu2,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      gam1 = gam1 - 1.0_GP
      cp1  = 2.0_GP / (gam1*smach*smach)
      nu2  = nu2 + nu/3.0_GP
#endif

#ifdef CMHD_
! Reads parameters for the compressible MHD runs from the
! namelist 'cmhdb' on the external file 'parameter.inp'
!     amach : Alfvenic Mach number

      IF (myrank.eq.0) THEN
         OPEN(1,file=sparam,status='unknown',form="formatted")
         READ(1,NML=cmhdb)
         CLOSE(1)
      ENDIF
      CALL MPI_BCAST(amach,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      cp2 = 1.0_GP / (amach*amach)
#endif

#ifdef MAGFIELD_
! Reads general configuration flags for runs with
! magnetic fields from the namelist 'dynamo' on
! the external file 'parameter.inp'
!     dyna : = 0 when stat=0 generates initial v and B (MAGFIELD_)
!            = 1 when stat.ne.0 imports v and generates B (MAGFIELD_)

      IF (myrank.eq.0) THEN
         OPEN(1,file=sparam,status='unknown',form="formatted")
         READ(1,NML=dynamo)
         CLOSE(1)
      ENDIF
      CALL MPI_BCAST(dyna,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

! Reads parameters for the magnetic field from the
! namelist 'magfield' on the external file 'parameter.inp'
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
         OPEN(1,file=sparam,status='unknown',form="formatted")
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
! Reads parameters for runs with a uniform magnetic
! field from the namelist 'uniformb' on the external
! file 'parameter.inp'
!     bx0: uniform magnetic field in x
!     by0: uniform magnetic field in y
!     bz0: uniform magnetic field in z

      bx0 = 0.0_GP
      by0 = 0.0_GP
      bz0 = 0.0_GP
      IF (myrank.eq.0) THEN
         OPEN(1,file=sparam,status='unknown',form="formatted")
         READ(1,NML=uniformb)
         CLOSE(1)
      ENDIF
      CALL MPI_BCAST(bx0,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(by0,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(bz0,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
#endif

#ifdef HALLTERM_
! Reads parameters for runs with the Hall effect
! from the namelist 'hallparam' on the external
! file 'parameter.inp'
!     ep  : amplitude of the Hall effect
!     gspe: = 0 skips generalized helicity spectrum computation
!           = 1 computes the spectrum of generalized helicity

      IF (myrank.eq.0) THEN
         OPEN(1,file=sparam,status='unknown',form="formatted")
         READ(1,NML=hallparam)
         CLOSE(1)
      ENDIF
      CALL MPI_BCAST(ep,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(gspe,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#endif

#ifdef ROTATION_
! Reads parameters for runs with rotation from the
! namelist 'rotation' on the external file 'parameter.inp'
!     omegax: amplitude of the uniform rotation along x
!     omegay: amplitude of the uniform rotation along y
!     omegaz: amplitude of the uniform rotation along z

      omegax = 0.0_GP
      omegay = 0.0_GP
      omegaz = 0.0_GP
      IF (myrank.eq.0) THEN
         OPEN(1,file=sparam,status='unknown',form="formatted")
         READ(1,NML=rotation)
         CLOSE(1)
      ENDIF
      CALL MPI_BCAST(omegax,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(omegay,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(omegaz,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
#endif

#ifdef WAVEFUNCTION_
! Reads parameters specifically for the GPE and ARGL solvers
! from the namelist 'wavefunction' on the external file
! 'parameter.inp'
!     cspeed : speed of sound
!     lambda : coherence length
!     rho0   : density at infinity (or at the zero of the potential)
!     kttherm: KT with T=thermalization temperature (for ARGL)
!     V0     : potential amplitude (for solvers with trapping potentials)
!     cflow  : =1 if generating counterflow (ARGL)
!     cflow_newt   : =1 if mean flow is needed for Newton method (ARGL)
!     dt_newt      : time step (preconditioner) for Newton method (ARGL)
!     iter_max_newt: max number of iterations for Newton method (ARGL)
!     iter_max_bicg: max number of iterations for biconjugate gradient (ARGL)
!     tol_newt     : tolerance for the Newton method (ARGL)
!     tolbicg_rel  : relarive tolerance for biconjugate gradient (ARGL)
!     zparam0-9    : ten real numbers to control properties of
!              the wavefunction

      rho0 = 1.0_GP        !Default value
      kttherm = 0.0_GP     !Default value
      V0 = 0.0_GP          !Default value
      cflow = 0            !Default value
      cflow_newt = 0       !Default value
      dt_newt = dt         !Default value
      iter_max_newt = 0    !Default value (no Newton done after ARGL)
      iter_max_bicg = 0    !Default value
      tol_newt = 0.0_GP    !Default value
      tolbicg_rel = 0.0_GP !Default value
      IF (myrank.eq.0) THEN
         OPEN(1,file=sparam,status='unknown',form="formatted")
         READ(1,NML=wavefunction)
         CLOSE(1)
         alpha = cspeed*lambda/sqrt(2.0_GP)
         omegag = cspeed/(lambda*sqrt(2.0_GP))
         beta= omegag/rho0
      ENDIF
      CALL MPI_BCAST(cspeed,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(lambda,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(rho0,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(alpha,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(beta,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(omegag,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(kttherm,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(V0,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(cflow,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(cflow_newt,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(dt_newt,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(iter_max_newt,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(iter_max_bicg,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(tol_newt,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(tolbicg_rel,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(zparam0,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(zparam1,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(zparam2,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(zparam3,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(zparam4,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(zparam5,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(zparam6,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(zparam7,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(zparam8,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(zparam9,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
#endif

#ifdef ALPHAV_
! Reads the value of alpha for the velocity field
! in runs using Lagrangian averaged subgrid models
!     alpk: filter length for the velocity field

      IF (myrank.eq.0) THEN
         OPEN(1,file=sparam,status='unknown',form="formatted")
         READ(1,NML=alphav)
         CLOSE(1)
      ENDIF
      CALL MPI_BCAST(alpk,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
#endif

#ifdef ALPHAB_
! Reads the value of alpha for the magnetic field
! in runs using Lagrangian averaged subgrid models
!     alpm: filter length for the magnetic field

      IF (myrank.eq.0) THEN
         OPEN(1,file=sparam,status='unknown',form="formatted")
         READ(1,NML=alphab)
         CLOSE(1)
      ENDIF
      CALL MPI_BCAST(alpm,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
#endif

#ifdef EDQNM_
! Reads the value of the Kolmogorov constant and a
! flag for the helicity LES in runs using EDQNM-based
! LES models
!     kolmo: Kolmogorov constant
!     heli:  = 0 helicity not taken into account
!            = 1 helicity taken into account

      kolmo = 1.4_GP !Default value
      heli = 0       !Default value
      IF (myrank.eq.0) THEN
         OPEN(1,file=sparam,status='unknown',form="formatted")
         READ(1,NML=edqnmles)
         CLOSE(1)
      ENDIF
      CALL MPI_BCAST(kolmo,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(heli,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#endif


! Before continuing, we verify that all parameters and compilation
! options are compatible with the SOLVER being used

      INCLUDE SOLVERCHECK_

! Initializes arrays and constants for the pseudospectral method

! Some constants for the FFT
!     kmax: maximum truncation for dealiasing
!     tiny: minimum truncation for dealiasing

      kmax =     1.0_GP/9.0_GP
      nmax =     int(max(nx*Dkx,ny*Dky,nz*Dkz)/Dkk)
      nmaxperp = int(max(nx*Dkx,ny*Dky)/Dkk)
#ifndef DEF_ARBSIZE_
      IF (anis.eq.0)  kmax = kmax*real(nx,kind=GP)**2
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
      IF (anis.eq.1) THEN
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
      IF (anis.eq.1) THEN
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
! FFTW 2.x only supports FFTW_ESTIMATE or FFTW_MEASURE

      nth = 1
!$    nth = omp_get_max_threads()
#if !defined(DEF_GHOST_CUDA_)
!$    CALL fftp3d_init_threads(ierr)
#endif
      n (1) = nx ; n (2) = ny ; n (3) = nz
      nt(1) = nxt; nt(2) = nyt; nt(3) = nzt
      CALL fftp3d_create_plan_comm(planrc,n,FFTW_REAL_TO_COMPLEX, &
                             FFTW_ESTIMATE, MPI_COMM_WORLD)
      CALL fftp3d_create_plan_comm(plancr,n,FFTW_COMPLEX_TO_REAL, &
                             FFTW_ESTIMATE, MPI_COMM_WORLD)

      ! Parse input index set, store in istat:
      CALL parseind(sstat,';', istat , 4096, nstat)

      ! NOTE: nt refers to full grid here, and n is 'truncated' grid:

      ! Plans for 'truncated' grid. Remember, nt > n here:
!     CALL trrange(1,nzt    ,nz    ,nprocs,myrank,ktsta,ktend)
!     CALL trrange(1,nzt/2+1,nz/2+1,nprocs,myrank,itsta,itend)
      ktsta = ksta; ktend = kend;
      itsta = ista; itend = iend;
      CALL fftp3d_create_trplan_comm(planrct,nt,n,FFTW_REAL_TO_COMPLEX,FFTW_MEASURE,MPI_COMM_WORLD)

      ! Reset indices:
      CALL range(1,nx/2+1,nprocs,myrank,ista,iend)
      CALL range(1,nz    ,nprocs,myrank,ksta,kend)

      ! Get new communicator, props  for truncated grid: 
!!    CALL create_trcomm(nt, n, MPI_COMM_WORLD, commtrunc, grouptrunc)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD, ntprocs, ierr)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD, mytrank, ierr)
      CALL io_init_comm(MPI_COMM_WORLD,(/nxt,nyt,nzt/),ktsta,ktend,planiot)

       
! Now that we have nmax we can allocate some temporary
! arrays needed to compute transfer functions in quantum
! flows

#ifdef WAVEFUNCTION_
      ALLOCATE( iold(nmax/2+1), qold(nmax/2+1) )
      ALLOCATE( kold(nmax/2+1), cold(nmax/2+1) )
      ALLOCATE( inew(nmax/2+1), qnew(nmax/2+1) )
      ALLOCATE( knew(nmax/2+1), cnew(nmax/2+1) )
#endif

! Initializes arrays to keep track of the forcing 
! if slowly evolving phases are used, and allocates
! arrays to compute mean flows if needed.

#ifdef VELOC_
      ampl = 1.0_GP
      timef = fstep
      IF (rand.eq.2) THEN
         ALLOCATE( fxold(nz,ny,ista:iend) )
         ALLOCATE( fyold(nz,ny,ista:iend) )
         ALLOCATE( fzold(nz,ny,ista:iend) )
         ALLOCATE( fxnew(nz,ny,ista:iend) )
         ALLOCATE( fynew(nz,ny,ista:iend) )
         ALLOCATE( fznew(nz,ny,ista:iend) )
      ENDIF
      IF (rand.eq.3) THEN
         ALLOCATE( Faux1(10) )
         ALLOCATE( Faux2(10) )
         DO i = 1,10
            Faux1(i) = ampl
            Faux2(i) = 0.0_GP
         END DO
      ENDIF
      IF (mean.eq.1) THEN
         ALLOCATE( M1(nz,ny,ista:iend) )
         ALLOCATE( M2(nz,ny,ista:iend) )
         ALLOCATE( M3(nz,ny,ista:iend) )
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
            DO j = 1,ny
               DO k = 1,nz
                  M1(k,j,i) = 0.0_GP
                  M2(k,j,i) = 0.0_GP
                  M3(k,j,i) = 0.0_GP
               END DO
            END DO
         END DO
#ifdef SCALAR_
         ALLOCATE( M7(nz,ny,ista:iend) )
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
            DO j = 1,ny
               DO k = 1,nz
                  M7(k,j,i) = 0.0_GP
               END DO
            END DO
         END DO
#endif
#ifdef MULTISCALAR_
         ALLOCATE( M8 (nz,ny,ista:iend) )
         ALLOCATE( M9 (nz,ny,ista:iend) )
         ALLOCATE( M10(nz,ny,ista:iend) )
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
            DO j = 1,ny
               DO k = 1,nz
                  M8 (k,j,i) = 0.0_GP
                  M9 (k,j,i) = 0.0_GP
                  M10(k,j,i) = 0.0_GP
               END DO
            END DO
         END DO
#endif
#ifdef MAGFIELD_
         ALLOCATE( M4(nz,ny,ista:iend) )
         ALLOCATE( M5(nz,ny,ista:iend) )
         ALLOCATE( M6(nz,ny,ista:iend) )
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
            DO j = 1,ny
               DO k = 1,nz
                  M4(k,j,i) = 0.0_GP
                  M5(k,j,i) = 0.0_GP
                  M6(k,j,i) = 0.0_GP
               END DO
            END DO
         END DO
#endif
      ENDIF
#endif
#ifdef MAGFIELD_
      IF (rand.eq.2) THEN
         ALLOCATE( mxold(nz,ny,ista:iend) )
         ALLOCATE( myold(nz,ny,ista:iend) )
         ALLOCATE( mzold(nz,ny,ista:iend) )
         ALLOCATE( mxnew(nz,ny,ista:iend) )
         ALLOCATE( mynew(nz,ny,ista:iend) )
         ALLOCATE( mznew(nz,ny,ista:iend) )
      ENDIF
#endif


#ifdef ADVECT_
!
! If doing a simulation with advection (for wavefunctions),
! we compute the square of the velocity field in real space
! and dealias it in Fourier space

!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
            DO k = 1,nz
               C1(k,j,i) = vx(k,j,i)
               C2(k,j,i) = vy(k,j,i)
               C3(k,j,i) = vz(k,j,i)
            END DO
         END DO
      END DO
      CALL fftp3d_complex_to_real(plancr,C1,R1,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,C2,R2,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,C3,R3,MPI_COMM_WORLD)
      rmp = 1.0_GP/(4*alpha*beta*(real(nx,kind=GP)* &
            real(ny,kind=GP)*real(nz,kind=GP))**2)
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
         DO j = 1,ny
            DO i = 1,nx
               vsq(i,j,k) = (R1(i,j,k)**2+R2(i,j,k)**2+ &
	                     R3(i,j,k)**2)*rmp
            END DO
         END DO
      END DO
      CALL fftp3d_real_to_complex(planrc,vsq,C1,MPI_COMM_WORLD)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
            DO k = 1,nz
               IF (kn2(k,j,i).gt.kmax) THEN
                  C1(k,j,i) = 0.0_GP
               ENDIF
            END DO
         END DO
      END DO
      CALL fftp3d_complex_to_real(plancr,C1,vsq,MPI_COMM_WORLD)
#endif
  
      ! Allocate 'trunc' arrays
      ALLOCATE( CT1(nzt,nyt,itsta:itend)) ! full res
      ALLOCATE( CT2(nzt,nyt,itsta:itend)) ! full res
      ALLOCATE( RT1(nxt,nyt,ktsta:ktend)) ! full res

! Cycle over all input times, and do analysis:
 LSTAT : DO t = 1, nstat

!       Read binary data
      WRITE(ext, fmtext) istat(t)
#ifdef VELOC_
        iidir = trim(idir) // '/outs'
        write(*,*) ' ext=', ext, ' iidir=', trim(iidir)
     write(*,*) 'main: input vx ...'
        CALL io_read(1,iidir,'vx',ext,planiot,RT1)
     write(*,*) 'main: vxR=', RT1(1:10,10,10)
        CALL fftp3d_real_to_complex(planrct,RT1,vx,MPI_COMM_WORLD)
     write(*,*) 'main: vxC=', CT1(1:10,10,10)
!       CALL trunc(CT1, nt, n, kmax, 1, CT2, vx)
!    write(*,*) 'main: vxT =', vx(1:10,10,10)

     write(*,*) 'main: input vy ...'
        CALL io_read(1,iidir,'vy',ext,planiot,RT1)
        CALL fftp3d_real_to_complex(planrct,RT1,vy,MPI_COMM_WORLD)
!       CALL trunc(CT1, nt, n, kmax, 1, CT2, vy)

     write(*,*) 'main: input vz ...'
        CALL io_read(1,iidir,'vz',ext,planiot,RT1)
        CALL fftp3d_real_to_complex(planrct,RT1,vz,MPI_COMM_WORLD)
!       CALL trunc(CT1, nt, n, kmax, 1, CT2, vz)
#endif
#ifdef BOUSSINESQ_
     write(*,*) 'main: input th ...'
        CALL io_read(1,iidir,'th',ext,planiot,RT1)
        CALL fftp3d_real_to_complex(planrct,RT1,th,MPI_COMM_WORLD)
!       CALL trunc(CT1, nt, n, kmax, 1, CT2, th)
#endif
        write(*,*) ' data loaded: index=', ext

        ! Compute time derivative estimate at this time step:
        INCLUDE RKSTEP1_
        DO o = ord,1,-1
           INCLUDE RKSTEP2_
        END DO 

#ifdef BOUSSINESQ_
        rmp = 1.0_GP/ &
              (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))
	! Compute RHS of dq/dt = RHS for each primitive.
        ! These compute the RHS of the filtered spatial
        ! operators acting on the filitered state:
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
        DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
          DO j = 1,ny
            DO k = 1,nz
              C1 (k,j,i) = rmp * ( vx(k,j,i) - C1 (k,j,i) )/dt
              C2 (k,j,i) = rmp * ( vy(k,j,i) - C2 (k,j,i) )/dt
              C3 (k,j,i) = rmp * ( vz(k,j,i) - C3 (k,j,i) )/dt
              C20(k,j,i) = rmp * ( th(k,j,i) - C20(k,j,i) )/dt
            ENDDO
          ENDDO
        ENDDO
        CALL fftp3d_complex_to_real(plancr,C1 ,R1,MPI_COMM_WORLD)
        CALL io_write(1,odir,'rhsvx_T',ext,planio,R1)
        CALL fftp3d_complex_to_real(plancr,C2 ,R1,MPI_COMM_WORLD)
        CALL io_write(1,odir,'rhsvy_T',ext,planio,R1)
        CALL fftp3d_complex_to_real(plancr,C3 ,R1,MPI_COMM_WORLD)
        CALL io_write(1,odir,'rhsvz_T',ext,planio,R1)
        CALL fftp3d_complex_to_real(plancr,C20,R1,MPI_COMM_WORLD)
        CALL io_write(1,odir,'rhsth_T',ext,planio,R1)
#endif
        write(*,*) ' done: index=', ext
      
      END DO LSTAT

!
! End of STAT loop


      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

!
! End of MAIN3D

      CALL GTFree(ihcpu1)
      CALL GTFree(ihomp1)
      CALL GTFree(ihwtm1)
      CALL GTFree(ihcpu2)
      CALL GTFree(ihomp2)
      CALL GTFree(ihwtm2)

      CALL fftp3d_destroy_plan(plancr)
      CALL fftp3d_destroy_plan(planrc)

      IF ( commtrunc .NE. MPI_COMM_NULL ) THEN
         CALL MPI_COMM_FREE(commtrunc, ierr)
      endif
      IF ( grouptrunc .NE. MPI_GROUP_NULL ) THEN
        CALL MPI_GROUP_FREE(grouptrunc, ierr)
      ENDIF

      CALL MPI_FINALIZE(ierr)

      DEALLOCATE( R1,R2,R3 )
      DEALLOCATE( RT1 )

      DEALLOCATE( C1,C2,C3,C4,C5,C6,C7,C8 )
      DEALLOCATE( CT1,CT2 )
      DEALLOCATE( kx,ky,kz )
      IF (anis.eq.1) THEN
         DEALLOCATE( kk2 )
      ELSE
         NULLIFY( kk2 )
      ENDIF
      DEALLOCATE( kn2 )
#ifdef VELOC_
      DEALLOCATE( fx,fy,fz )
      IF (mean.eq.1) DEALLOCATE( M1,M2,M3 )
      IF (rand.eq.2) DEALLOCATE( fxold, fyold, fzold )
      IF (rand.eq.2) DEALLOCATE( fxnew, fynew, fznew )
      IF (rand.eq.3) DEALLOCATE( Faux1, Faux2 )
#endif
#if defined(VELOC_) || defined (ADVECT_)
      DEALLOCATE( vx,vy,vz )
#endif
#ifdef ADVECT_
      DEALLOCATE( vsq )
#endif
#ifdef MAGFIELD_
      DEALLOCATE( ax,ay,az,mx,my,mz )
      DEALLOCATE( C9,C10,C11,C12,C13,C14,C15,C16,C17 )
      IF (mean.eq.1) DEALLOCATE( M4,M5,M6 )
      IF (rand.eq.2) DEALLOCATE( mxold, myold, mzold )
      IF (rand.eq.2) DEALLOCATE( mxnew, mynew, mznew )
#endif
#ifdef HALLTERM_
      DEALLOCATE( C18 )
#endif
#ifdef EDQNM_
      DEALLOCATE( C19 )
      DEALLOCATE( tepq,thpq,tve,tvh,Eext,Hext )
#endif
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
#ifdef COMPR_AUX_ARR_
      DEALLOCATE( C25,C26,C27 )
#endif
#ifdef WAVEFUNCTION_
      DEALLOCATE( zre,zim )
#endif
#ifdef QFORCE_
      DEALLOCATE( fre,fim )
#endif
#ifdef TRAP_
      DEALLOCATE( C28,C29,C30 )
      DEALLOCATE( Vtrap,Vlinx,Vliny )
#endif

      DEALLOCATE( R4,R5,R6 )

      END PROGRAM MAIN3D


