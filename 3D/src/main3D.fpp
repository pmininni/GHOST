!=================================================================
      PROGRAM MAIN3D
!=================================================================
! GHOST code: Geophysical High Order Suite for Turbulence
!
! Numerically integrates several fluid dynamics equations
! in 2 and 3 dimensions with periodic boundary conditions 
! and external forcing. A pseudo-spectral method is used to 
! compute spatial derivatives, while an adjustable order 
! Runge-Kutta method is used to evolve the system in the time 
! domain. To compile, you need the FFTW library installed on 
! your system. The parallel FFT is in the FFTP subroutines
! and uses the FFTPLANS and MPIVARS modules (see the file 
! 'fftp_mod.f90' for details).
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
      USE boxsize
      USE gtimer
      USE fftplans
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
#if defined(PART_)
      USE class_GPart
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

      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: C1,C2
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: C3,C4
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: C5,C6
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: C7,C8
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
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: C30,C31
#endif

#ifdef WAVEFUNCTION_
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:)     :: iold,qold,kold,cold
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:)     :: inew,qnew,knew,cnew
#endif

      REAL(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:)    :: R1,R2,R3
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
#ifdef PART_
      REAL(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:)    :: R4,R5,R6
#endif
#if defined(INERPART_)
      REAL(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:)    :: Rv1,Rv2
#endif
#if defined(TESTPART_) && defined(MAGFIELD_)
      REAL(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:)    :: Rb1,Rb2,Rb3
      REAL(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:)    :: Rj1,Rj2,Rj3
#endif

!
! Auxiliary variables

      COMPLEX(KIND=GP) :: cdump,jdump
      COMPLEX(KIND=GP) :: cdumq,jdumq
      COMPLEX(KIND=GP) :: cdumr,jdumr
      DOUBLE PRECISION :: tmp,tmq,tmr
      DOUBLE PRECISION :: eps,epm
!$    DOUBLE PRECISION, EXTERNAL :: omp_get_wtime

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
#ifdef PART_
      REAL         :: rbal
      REAL(KIND=GP):: tbeta(3)
      INTEGER      :: maxparts
      INTEGER      :: injtp
      INTEGER      :: cresetp
      INTEGER      :: ilginittype
      INTEGER      :: ilgintrptype
      INTEGER      :: ilgexchtype
      INTEGER      :: ilgouttype
      INTEGER      :: ilgwrtunit
      INTEGER      :: ilgcoll
      INTEGER      :: ilgfpfiletype
      INTEGER      :: blgdofp
      INTEGER      :: blgfpfilecoll
      INTEGER      :: dolag
      INTEGER      :: dopacc
      INTEGER      :: nwpart
#endif
#ifdef LAGPART_
      TYPE (GPart) :: lagpart,lagfp
#endif
#if defined(INERPART_)
      INTEGER               :: dolightp, donldrag
      REAL(KIND=GP)         :: tau, grav, gamma
      TYPE (InerGPart)      :: lagpart
#endif
#if defined(TESTPART_) && defined(MAGFIELD_)
      INTEGER               :: dokinelv, dokinelp
      REAL(KIND=GP)         :: gyrof, vtherm, dii      
      TYPE (TestGPart)      :: lagpart
#endif
!$    INTEGER, EXTERNAL     :: omp_get_max_threads

#if defined(DEF_GHOST_CUDA_)
      TYPE(cudaDevicePropG) :: devprop
#endif
      TYPE(IOPLAN)          :: planio
      CHARACTER(len=100)    :: odir,idir
#ifdef PART_
      CHARACTER(len=1024)   :: lgseedfile,slgfpfile
#endif
      LOGICAL               :: bbenchexist

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
#ifdef PART_
      NAMELIST / plagpart / lgmult,maxparts,ilginittype,ilgintrptype
      NAMELIST / plagpart / ilgexchtype,ilgouttype,ilgwrtunit,lgseedfile
      NAMELIST / plagpart / injtp,ilgcoll,cresetp,dolag,dopacc
      NAMELIST / plagpart / blgdofp,ilgfpfiletype,blgfpfilecoll,slgfpfile
#endif
#if defined(INERPART_)
      NAMELIST / pinerpart / tau,grav,gamma,dolightp,donldrag
#endif
#if defined(TESTPART_) && defined(MAGFIELD_)
      NAMELIST / ptestpart / gyrof,vtherm,dii,dokinelv,dokinelp
#endif

!
! Initialization

! Initializes the MPI and I/O libraries
      CALL MPI_INIT_THREAD(MPI_THREAD_FUNNELED,provided,ierr)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)

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
     CALL io_init(myrank,(/nx,ny,nz/),ksta,kend,planio)

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
      ALLOCATE( C30(nz,ny,ista:iend), C31(nz,ny,ista:iend) )
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
#ifdef PART_
      ALLOCATE( R4(nx,ny,ksta:kend) )
      ALLOCATE( R5(nx,ny,ksta:kend) )
      ALLOCATE( R6(nx,ny,ksta:kend) )
#endif
#if defined (INERPART_)
      ALLOCATE( Rv1(nx,ny,ksta:kend) )
      ALLOCATE( Rv2(nx,ny,ksta:kend) )
#endif
#if defined (TESTPART_) && defined(MAGFIELD_)
      ALLOCATE( Rb1(nx,ny,ksta:kend) )        
      ALLOCATE( Rb2(nx,ny,ksta:kend) )
      ALLOCATE( Rb3(nx,ny,ksta:kend) )
      ALLOCATE( Rj1(nx,ny,ksta:kend) )        
      ALLOCATE( Rj2(nx,ny,ksta:kend) )
      ALLOCATE( Rj3(nx,ny,ksta:kend) )
#endif
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

!
! The following lines read the file 'parameter.inp'

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
#ifdef DEF_ARBSIZE_
! Reads parameters to set the box size 
!     Lx  : Length in x (in units of 2.pi, =1 gives a side of length 2.pi)
!     Ly  : Length in y
!     Lz  : Length in z
!     Dkk : Width of Fourier shells for 2D and 3D spectra
!           Default = min(1/Lx, 1/Ly, 1/Lz)

      IF (myrank.eq.0) THEN
         OPEN(1,file='parameter.inp',status='unknown',form="formatted")
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
         OPEN(1,file='parameter.inp',status='unknown',form="formatted")
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
         OPEN(1,file='parameter.inp',status='unknown',form="formatted")
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
         OPEN(1,file='parameter.inp',status='unknown',form="formatted")
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
         OPEN(1,file='parameter.inp',status='unknown',form="formatted")
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
         OPEN(1,file='parameter.inp',status='unknown',form="formatted")
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
         OPEN(1,file='parameter.inp',status='unknown',form="formatted")
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
         OPEN(1,file='parameter.inp',status='unknown',form="formatted")
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
         OPEN(1,file='parameter.inp',status='unknown',form="formatted")
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
         OPEN(1,file='parameter.inp',status='unknown',form="formatted")
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
         OPEN(1,file='parameter.inp',status='unknown',form="formatted")
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
         OPEN(1,file='parameter.inp',status='unknown',form="formatted")
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
!     rho0   : density at infinity
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
         OPEN(1,file='parameter.inp',status='unknown',form="formatted")
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
         OPEN(1,file='parameter.inp',status='unknown',form="formatted")
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
         OPEN(1,file='parameter.inp',status='unknown',form="formatted")
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
         OPEN(1,file='parameter.inp',status='unknown',form="formatted")
         READ(1,NML=edqnmles)
         CLOSE(1)
      ENDIF
      CALL MPI_BCAST(kolmo,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(heli,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#endif

#ifdef PART_
! Reads parameters for runs with particles
!     maxparts    : Maximum number of particles
!     injtp       : = 0 when stat=0 generates initial v and part seeds 
!                   = 1 when stat.ne.0 imports v and generates initial 
!                     part seeds. If injtp=0 when stat.ne.0, then the 
!                     particle restart file is read that corresp. to 
!                     stat.
!     cresetp     : = 0 don't reset counters when injtp=1;
!                   = 1: _do_ reset counters when injtp=1.
!     lgmult      : Multiplier for part output (must divide tstep evenly)
!     ilginittype : Inititialization type: either GPINIT_RANDLOC or
!                   GPINIT_USERLOC
!     ilgintrptype: Interpolation type: only GPINTRP_CSPLINE currently
!     ilgexchtype : Boundary exchange type: GPEXCHTYPE_VDB (voxel db)
!                   or GPEXCHTYPE_NN (nearest-neighbor)
!     ilgouttype  : Particle output type: 0=binary; 1=ASCII
!     ilgwrtunit  : Units for part position write: 0=box units; 1=grid units
!     lgseedfile  : Name of seed file if ilginittype=GPINIT_USERLOC
!     ilgcoll     : 1=binary collective I/O; 0=task 0 binary (posix) I/O
!     dolag       : 1=run with particles; 0=don't 
!     dopacc      : 1=compute acceleration internally; 0=don't
!     XlgfpX      : Parameters for 'fixed point' particles (test of
!                   frozen-in, only for LAGPART)

      maxparts     = 1000
      injtp        = 0
      cresetp      = 0
      ilginittype  = GPINIT_RANDLOC
      ilgintrptype = GPINTRP_CSPLINE
      ilgexchtype  = GPEXCHTYPE_VDB
      ilgouttype   = 0
      ilgwrtunit   = 0
      lgmult       = 1
      lgseedfile   = 'gplag.dat'
      ilgcoll      = 1
      rbal         = 0.0
      dolag        = 1
      dopacc       = 1
      nwpart       = 0
      blgdofp      = 0
      ilgfpfiletype= 0
      blgfpfilecoll= 1
      slgfpfile    = 'xlgInitRndSeed.000.txt'
      IF (myrank.eq.0) THEN
         OPEN(1,file='parameter.inp',status='unknown',form="formatted")
         READ(1,NML=plagpart)
         CLOSE(1)
      ENDIF
#if defined(VPART_)
      dopacc       = 0 ! Lag. accel. not supported for particles w/velocity
#endif
      CALL MPI_BCAST(maxparts     ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(injtp        ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(cresetp      ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(lgmult       ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(ilginittype  ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(ilgintrptype ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(ilgexchtype  ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(ilgouttype   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(ilgwrtunit   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(dolag        ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(dopacc       ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(ilgcoll      ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(blgdofp      ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(ilgfpfiletype,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(blgfpfilecoll,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(slgfpfile ,1024,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(lgseedfile,1024,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      IF ( mod(tstep,lgmult).NE.0 ) THEN
        WRITE(*,*)'main: lgmult must divide tstep evenly'
        STOP
      ENDIF  
      pstep = tstep/lgmult
#endif

#if defined(INERPART_)
! Reads parameters for runs with inertial particles
!     tau      : Stokes time
!     grav     : Effective gravity acceleration
!     gamma    : Mass ratio (= m_f/m_p, f:fluid, p:particle)
!     dolightp : = 0 do not compute mass ratio terms (heavy particles);
!                = 1 compute mass ratio terms (light particles);
!     donldrag : = 0 compute linear drag;
!                = 1 compute non-linear drag (only for light particles);

      gamma    = 0.0_GP
      dolightp = 0
      donldrag = 0
      IF (myrank.eq.0) THEN
         OPEN(1,file='parameter.inp',status='unknown',form="formatted")
         READ(1,NML=pinerpart)
      CLOSE(1)
      ENDIF
      CALL MPI_BCAST(tau     ,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(grav    ,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(gamma   ,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(dolightp,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(donldrag,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#endif

#if defined(TESTPART_) && defined(MAGFIELD_)
! Reads parameters for runs with test particles
!     gyrof    : Gyrofrequency
!     vtherm   : Thermal velocity of the test particles
!     dii      : Ion inertial length scale (=epsilon in Hall-MHD)
!     dokinelv : = 0 Compute the uxB term in the Lorentz force
!                = 1 Compute also electron velocity corrections (Hall-MHD)
!     dokinelp : = 0 Compute Ohmic currents in the Lorentz force
!                = 1 Compute also the electron pressure (ambipolar diff.)

      dii = 0.0_GP
      dokinelv = 0
      dokinelv = 0
      IF (myrank.eq.0) THEN
         OPEN(1,file='parameter.inp',status='unknown',form="formatted")
         READ(1,NML=ptestpart)
         CLOSE(1)
      ENDIF
      CALL MPI_BCAST(gyrof    ,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(vtherm   ,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(dii      ,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(dokinelv,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(dokinelp,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
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
      tiny  = min(1e-5_GP,.1_GP/real(nmax,kind=GP))
      tinyf = min(1e-15_GP,.1_GP/real(nmax,kind=GP))

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

#ifdef PART_
      rmp   = 0.5/dt
      tbeta = (/-rmp, 0.0, rmp/)
      IF ( dolag.GT.0 ) THEN
        CALL lagpart%GPart_ctor(MPI_COMM_WORLD,maxparts,ilginittype,&
             ilgintrptype,3,ilgexchtype,ilgouttype,ilgcoll,csize, &
             nstrip,dopacc,ilgwrtunit)
        CALL lagpart%SetRandSeed(seed)
        CALL lagpart%SetSeedFile(trim(lgseedfile))
#if defined(INERPART_) && !defined(ROTATION_)
        CALL lagpart%InerGPart_ctor(tau,grav,gamma,nu,donldrag)
#endif
#if defined(INERPART_) && defined(ROTATION_)
        CALL lagpart%InerGPart_ctor(tau,grav,gamma,nu,donldrag, &
	om=(/omegax,omegay,omegaz/),                            &
	x0=(/(real(nx+1,kind=GP)/2),(real(ny+1,kind=GP)/2),(real(nz+1,kind=GP)/2)/))
#endif
#if defined(TESTPART_) && defined(MAGFIELD_)
        CALL lagpart%TestGPart_ctor()
#endif
      ENDIF
#endif
#ifdef LAGPART_
      IF ( blgdofp.GT.0 ) THEN ! 'Fixed point' particles: interpolates
                               ! Eulerian fields at given points, and
                               ! outputs them without evolving them.
        CALL lagfp%GPart_ctor(MPI_COMM_WORLD,maxparts,GPINIT_USERLOC,&
             ilgintrptype,3,ilgexchtype,ilgfpfiletype,blgfpfilecoll,&
             csize,nstrip,0)
        CALL lagfp%SetRandSeed(seed)
        CALL lagfp%SetSeedFile(trim(slgfpfile))
      ENDIF
#endif

!
! Sets the external forcing

#ifdef VELOC_
      INCLUDE 'initialfv.f90'           ! mechanical forcing
#endif
#ifdef SCALAR_
      INCLUDE 'initialfs.f90'           ! scalar source
#endif
#ifdef MULTISCALAR_
      INCLUDE 'initialfms.f90'          ! multiscalar sources
#endif
#ifdef MAGFIELD_
      INCLUDE 'initialfb.f90'           ! electromotive forcing
#endif
#ifdef ADVECT_
      timef = fstep
      INCLUDE 'initialfq.f90'           ! quantum thermal forcing
#endif
#ifdef TRAP_
      INCLUDE 'potentialq.f90'          ! quantum trapping potential
#endif

! If stat=0 we start a new run.
! Generates initial conditions for the fields and particles.
 IC : IF (stat.eq.0) THEN

      ini = 1
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
#ifdef SCALAR_
      INCLUDE 'initials.f90'            ! initial scalar density
#endif
#ifdef MULTISCALAR_
      INCLUDE 'initialms.f90'           ! initial multiscalar conc.
#endif
#ifdef MAGFIELD_
      INCLUDE 'initialb.f90'            ! initial vector potential
#endif
#ifdef WAVEFUNCTION_
      INCLUDE 'initialz.f90'            ! initial wave function
#endif
#ifdef PART_
      IF ( dolag.GT.0 ) THEN
        CALL lagpart%Init()
#if defined(INERPART_)
         rmp = 1.0_GP/(real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
           DO j = 1,ny
             DO k = 1,nz
               C7(k,j,i) = vx(k,j,i)*rmp
             END DO
           END DO
         END DO
         CALL fftp3d_complex_to_real(plancr,C7,R1,MPI_COMM_WORLD)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
           DO j = 1,ny
             DO k = 1,nz
               C7(k,j,i) = vy(k,j,i)*rmp
             END DO
           END DO
         END DO
         CALL fftp3d_complex_to_real(plancr,C7,R2,MPI_COMM_WORLD)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
           DO j = 1,ny
             DO k = 1,nz
               C7(k,j,i) = vz(k,j,i)*rmp
             END DO
           END DO
         END DO
         CALL fftp3d_complex_to_real(plancr,C7,R3,MPI_COMM_WORLD)
         CALL lagpart%InitVel(R1,R2,R3,Rv1,Rv2)
#endif
#if defined(TESTPART_) && defined(MAGFIELD_)
        CALL lagpart%InitVel(vtherm)
#endif
      ENDIF
#ifdef LAGPART_
      IF ( blgdofp.GT.0 ) THEN
        CALL lagfp%Init()
      ENDIF
#endif
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

#ifdef VELOC_
      CALL io_read(1,idir,'vx',ext,planio,R1)
      CALL io_read(1,idir,'vy',ext,planio,R2)
      CALL io_read(1,idir,'vz',ext,planio,R3)
      CALL fftp3d_real_to_complex(planrc,R1,vx,MPI_COMM_WORLD)
      CALL fftp3d_real_to_complex(planrc,R2,vy,MPI_COMM_WORLD)
      CALL fftp3d_real_to_complex(planrc,R3,vz,MPI_COMM_WORLD)

      IF (rand.eq.2) THEN
         CALL io_read(1,idir,'fxnew',ext,planio,R1)
         CALL io_read(1,idir,'fynew',ext,planio,R2)
         CALL io_read(1,idir,'fznew',ext,planio,R3)
         CALL fftp3d_real_to_complex(planrc,R1,fxnew,MPI_COMM_WORLD)
         CALL fftp3d_real_to_complex(planrc,R2,fynew,MPI_COMM_WORLD)
         CALL fftp3d_real_to_complex(planrc,R3,fznew,MPI_COMM_WORLD)

         CALL io_read(1,idir,'fxold',ext,planio,R1)
         CALL io_read(1,idir,'fyold',ext,planio,R2)
         CALL io_read(1,idir,'fzold',ext,planio,R3)
         CALL fftp3d_real_to_complex(planrc,R1,fxold,MPI_COMM_WORLD)
         CALL fftp3d_real_to_complex(planrc,R2,fyold,MPI_COMM_WORLD)
         CALL fftp3d_real_to_complex(planrc,R3,fzold,MPI_COMM_WORLD)
      ENDIF

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
            DO j = 1,ny
               DO k = 1,nz
                  M1(k,j,i) = dump*M1(k,j,i)
                  M2(k,j,i) = dump*M2(k,j,i)
                  M3(k,j,i) = dump*M3(k,j,i)
               END DO
            END DO
         END DO
      ENDIF
#endif

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
               DO j = 1,ny
                  DO k = 1,nz
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
 INJM: IF (injtm.eq.0) THEN
         CALL io_read(1,idir,'th1',ext,planio,R1)
         CALL fftp3d_real_to_complex(planrc,R1,th1,MPI_COMM_WORLD)
         CALL io_read(1,idir,'th2',ext,planio,R1)
         CALL fftp3d_real_to_complex(planrc,R1,th2,MPI_COMM_WORLD)
         CALL io_read(1,idir,'th3',ext,planio,R1)
         CALL fftp3d_real_to_complex(planrc,R1,th3,MPI_COMM_WORLD)
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
               DO j = 1,ny
                  DO k = 1,nz
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
      ENDIF INJM
#endif

#ifdef MAGFIELD_
      IF (rand.eq.2) THEN
         CALL io_read(1,idir,'mxnew',ext,planio,R1)
         CALL io_read(1,idir,'mynew',ext,planio,R2)
         CALL io_read(1,idir,'mznew',ext,planio,R3)
         CALL fftp3d_real_to_complex(planrc,R1,mxnew,MPI_COMM_WORLD)
         CALL fftp3d_real_to_complex(planrc,R2,mynew,MPI_COMM_WORLD)
         CALL fftp3d_real_to_complex(planrc,R3,mznew,MPI_COMM_WORLD)

         CALL io_read(1,idir,'mxold',ext,planio,R1)
         CALL io_read(1,idir,'myold',ext,planio,R2)
         CALL io_read(1,idir,'mzold',ext,planio,R3)
         CALL fftp3d_real_to_complex(planrc,R1,mxold,MPI_COMM_WORLD)
         CALL fftp3d_real_to_complex(planrc,R2,myold,MPI_COMM_WORLD)
         CALL fftp3d_real_to_complex(planrc,R3,mzold,MPI_COMM_WORLD)
      ENDIF
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
               DO j = 1,ny
                  DO k = 1,nz
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

#ifdef WAVEFUNCTION_
      CALL io_read(1,idir,'phi_re',ext,planio,R1)
      CALL io_read(1,idir,'phi_im',ext,planio,R2)
      CALL fftp3d_real_to_complex(planrc,R1,zre,MPI_COMM_WORLD)
      CALL fftp3d_real_to_complex(planrc,R2,zim,MPI_COMM_WORLD)
#endif

#ifdef ADVECT_
      INCLUDE 'initialv.f90'         ! Recreate velocity for advection
#endif

#ifdef PART_
      IF ( dolag.GT.0 ) THEN
        IF (injtp.eq.0) THEN
          WRITE(lgext, lgfmtext) pind
          CALL lagpart%io_read (1,idir,'xlg',lgext)
#if defined(INERPART_)
          CALL lagpart%io_readv(1,idir,'vip',lgext)
#endif
#if defined(TESTPART_) && defined(MAGFIELD_)
          CALL lagpart%io_readv(1,idir,'vtp',lgext)
#endif
        ELSE
          CALL lagpart%Init()
#if defined(INERPART_)
          rmp = 1.0_GP/(real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
          DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
            DO j = 1,ny
              DO k = 1,nz
                C7(k,j,i) = vx(k,j,i)*rmp
              END DO
            END DO
          END DO
          CALL fftp3d_complex_to_real(plancr,C7,R1,MPI_COMM_WORLD)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
          DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
            DO j = 1,ny
              DO k = 1,nz
                C7(k,j,i) = vy(k,j,i)*rmp
              END DO
            END DO
          END DO
          CALL fftp3d_complex_to_real(plancr,C7,R2,MPI_COMM_WORLD)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
          DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
            DO j = 1,ny
              DO k = 1,nz
                C7(k,j,i) = vz(k,j,i)*rmp
              END DO
            END DO
          END DO
          CALL fftp3d_complex_to_real(plancr,C7,R3,MPI_COMM_WORLD)
          CALL lagpart%InitVel(R1,R2,R3,Rv1,Rv2)
#endif
#if defined(TESTPART_) && defined(MAGFIELD_)
          CALL lagpart%InitVel(vtherm)
#endif
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
#ifdef LAGPART_
      IF ( blgdofp.GT.0 ) THEN
        CALL lagfp%Init()
      ENDIF
#endif
#endif

      ENDIF IC

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

! Updates the external forcing. Every 'fsteps'
! the phase or amplitude is changed according 
! to the value of 'rand'.
 TF :    IF (timef.eq.fstep) THEN
            timef = 0

            IF (rand.eq.1) THEN      ! randomizes phases

#ifdef VELOC_
               IF (myrank.eq.0) phase = 2*pi*randu(seed)
               CALL MPI_BCAST(phase,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
               cdump = COS(phase)+im*SIN(phase)
	       CALL phaseshift(fx,cdump)
	       CALL phaseshift(fy,cdump)
	       CALL phaseshift(fz,cdump)
#endif
#if defined(SCALAR_) || defined(MULTISCALAR_)
               IF (myrank.eq.0) phase = 2*pi*randu(seed)
               CALL MPI_BCAST(phase,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
               cdumr = COS(phase)+im*SIN(phase)
#endif
#if defined(SCALAR_)
	       CALL phaseshift(fs,cdumr)
#endif
#if defined(MULTISCALAR_)
	       CALL phaseshift(fs1,cdumr)
	       CALL phaseshift(fs2,cdumr)
	       CALL phaseshift(fs3,cdumr)
#endif
#ifdef MAGFIELD_
               IF (myrank.eq.0) phase = 2*pi*randu(seed)
               CALL MPI_BCAST(phase,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
               cdumq = corr*cdump+(1-corr)*(COS(phase)+im*SIN(phase))
	       CALL phaseshift(mx,cdump)
	       CALL phaseshift(my,cdump)
	       CALL phaseshift(mz,cdump)
#endif
#ifdef QFORCE_
               INCLUDE 'initialfq.f90'  ! generates a new forcing function
#endif

            ELSE IF (rand.eq.2) THEN ! slowly varying phases

#ifdef VELOC_
! Keeps a copy of the old forcing
               DO i = ista,iend
                  DO j = 1,ny
                     DO k = 1,nz
                        fxold(k,j,i) = fx(k,j,i) 
                        fyold(k,j,i) = fy(k,j,i) 
                        fzold(k,j,i) = fz(k,j,i) 
                     END DO
                  END DO
               END DO
! Creates a new random forcing. 'initialfv.f90'
! should be chosen to generate a forcing with random
! phases ('initialfv.f90_patterson' is recommended).
               INCLUDE 'initialfv.f90'
! Copies the new forcing to arrays for the target forcing
               DO i = ista,iend
                  DO j = 1,ny
                     DO k = 1,nz
                        fxnew(k,j,i) = fx(k,j,i) 
                        fynew(k,j,i) = fy(k,j,i) 
                        fznew(k,j,i) = fz(k,j,i) 
                     END DO
                  END DO
               END DO
#endif

#ifdef MAGFIELD_
! Keeps a copy of the old forcing
               DO i = ista,iend
                  DO j = 1,ny
                     DO k = 1,nz
                        mxold(k,j,i) = mx(k,j,i) 
                        myold(k,j,i) = my(k,j,i) 
                        mzold(k,j,i) = mz(k,j,i) 
                     END DO
                  END DO
               END DO
! Creates a new random forcing. 'initialfb.f90'
! should be chosen to generate forcing with random
! phases ('initialfb.f90_patterson' is recommended).
               INCLUDE 'initialfb.f90'
! Copies the new forcing to arrays for the target forcing
               DO i = ista,iend
                  DO j = 1,ny
                     DO k = 1,nz
                        mxnew(k,j,i) = mx(k,j,i) 
                        mynew(k,j,i) = my(k,j,i) 
                        mznew(k,j,i) = mz(k,j,i) 
                     END DO
                  END DO
               END DO
#endif

            ELSE IF (rand.eq.3) THEN ! user-defined forcing scheme

               INCLUDE 'incscheme.f90'

            END IF

         END IF TF

         IF (rand.eq.2) THEN ! Updates forcing if slowly varying
            rmp = FLOAT(timef+1)/float(fstep)
#ifdef VELOC_
            DO i = ista,iend
               DO j = 1,ny
                  DO k = 1,nz
                     fx(k,j,i) = (1-rmp)*fxold(k,j,i)+rmp*fxnew(k,j,i)
                     fy(k,j,i) = (1-rmp)*fyold(k,j,i)+rmp*fynew(k,j,i)
                     fz(k,j,i) = (1-rmp)*fzold(k,j,i)+rmp*fznew(k,j,i)
                  END DO
               END DO
            END DO
#endif
#ifdef MAGFIELD_
            DO i = ista,iend
               DO j = 1,ny
                  DO k = 1,nz
                     mx(k,j,i) = (1-rmp)*mxold(k,j,i)+rmp*mxnew(k,j,i)
                     my(k,j,i) = (1-rmp)*myold(k,j,i)+rmp*mynew(k,j,i)
                     mz(k,j,i) = (1-rmp)*mzold(k,j,i)+rmp*mznew(k,j,i)
                  END DO
               END DO
            END DO
#endif
         END IF

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
            IF (rand.eq.2) THEN
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
               DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
                  DO j = 1,ny
                     DO k = 1,nz
                        C1(k,j,i) = fxold(k,j,i)*rmp
                        C2(k,j,i) = fyold(k,j,i)*rmp
                        C3(k,j,i) = fzold(k,j,i)*rmp
                     END DO
                  END DO
               END DO
               CALL fftp3d_complex_to_real(plancr,C1,R1,MPI_COMM_WORLD)
               CALL fftp3d_complex_to_real(plancr,C2,R2,MPI_COMM_WORLD)
               CALL fftp3d_complex_to_real(plancr,C3,R3,MPI_COMM_WORLD)
               CALL io_write(1,odir,'fxold',ext,planio,R1)
               CALL io_write(1,odir,'fyold',ext,planio,R2)
               CALL io_write(1,odir,'fzold',ext,planio,R3)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
               DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
                  DO j = 1,ny
                     DO k = 1,nz
                        C1(k,j,i) = fxnew(k,j,i)*rmp
                        C2(k,j,i) = fynew(k,j,i)*rmp
                        C3(k,j,i) = fznew(k,j,i)*rmp
                     END DO
                  END DO
               END DO
               CALL fftp3d_complex_to_real(plancr,C1,R1,MPI_COMM_WORLD)
               CALL fftp3d_complex_to_real(plancr,C2,R2,MPI_COMM_WORLD)
               CALL fftp3d_complex_to_real(plancr,C3,R3,MPI_COMM_WORLD)
               CALL io_write(1,odir,'fxnew',ext,planio,R1)
               CALL io_write(1,odir,'fynew',ext,planio,R2)
               CALL io_write(1,odir,'fznew',ext,planio,R3)
            END IF
            IF (mean.eq.1) THEN
               dump = real(cstep,kind=GP)/t
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
               DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
                  DO j = 1,ny
                     DO k = 1,nz
                        C1(k,j,i) = dump*M1(k,j,i)*rmp
                        C2(k,j,i) = dump*M2(k,j,i)*rmp
                        C3(k,j,i) = dump*M3(k,j,i)*rmp
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
#endif
#ifdef SCALAR_
            rmp = 1.0_GP/ &
	          (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
            DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
               DO j = 1,ny
                  DO k = 1,nz
                     C1(k,j,i) = th(k,j,i)*rmp
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
                  DO j = 1,ny
                     DO k = 1,nz
                        C1(k,j,i) = dump*M7(k,j,i)*rmp
                     END DO
                  END DO
               END DO
               CALL fftp3d_complex_to_real(plancr,C1,R1,MPI_COMM_WORLD)
               CALL io_write(1,odir,'mean_th',ext,planio,R1)
            ENDIF
#endif
#ifdef MULTISCALAR_
            rmp = 1.0_GP/(real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
            DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
               DO j = 1,ny
                  DO k = 1,nz
                     C1(k,j,i) = th1(k,j,i)*rmp
                     C2(k,j,i) = th2(k,j,i)*rmp
                     C3(k,j,i) = th3(k,j,i)*rmp
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
                  DO j = 1,ny
                     DO k = 1,nz
                        C1(k,j,i) = dump*M8 (k,j,i)*rmp
                        C2(k,j,i) = dump*M9 (k,j,i)*rmp
                        C3(k,j,i) = dump*M10(k,j,i)*rmp
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
            rmp = 1.0_GP/(real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
            DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
               DO j = 1,ny
                  DO k = 1,nz
                     C1(k,j,i) = ax(k,j,i)*rmp
                     C2(k,j,i) = ay(k,j,i)*rmp
                     C3(k,j,i) = az(k,j,i)*rmp
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
            IF (rand.eq.2) THEN
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
               DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
                  DO j = 1,ny
                     DO k = 1,nz
                        C1(k,j,i) = mxold(k,j,i)*rmp
                        C2(k,j,i) = myold(k,j,i)*rmp
                        C3(k,j,i) = mzold(k,j,i)*rmp
                     END DO
                  END DO
               END DO
               CALL fftp3d_complex_to_real(plancr,C1,R1,MPI_COMM_WORLD)
               CALL fftp3d_complex_to_real(plancr,C2,R2,MPI_COMM_WORLD)
               CALL fftp3d_complex_to_real(plancr,C3,R3,MPI_COMM_WORLD)
               CALL io_write(1,odir,'mxold',ext,planio,R1)
               CALL io_write(1,odir,'myold',ext,planio,R2)
               CALL io_write(1,odir,'mzold',ext,planio,R3)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
               DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
                  DO j = 1,ny
                     DO k = 1,nz
                        C1(k,j,i) = mxnew(k,j,i)*rmp
                        C2(k,j,i) = mynew(k,j,i)*rmp
                        C3(k,j,i) = mznew(k,j,i)*rmp
                     END DO
                  END DO
               END DO
               CALL fftp3d_complex_to_real(plancr,C1,R1,MPI_COMM_WORLD)
               CALL fftp3d_complex_to_real(plancr,C2,R2,MPI_COMM_WORLD)
               CALL fftp3d_complex_to_real(plancr,C3,R3,MPI_COMM_WORLD)
               CALL io_write(1,odir,'mxnew',ext,planio,R1)
               CALL io_write(1,odir,'mynew',ext,planio,R2)
               CALL io_write(1,odir,'mznew',ext,planio,R3)
            END IF
            IF (mean.eq.1) THEN
               dump = real(cstep,kind=GP)/t
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
               DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
                  DO j = 1,ny
                     DO k = 1,nz
                        C1(k,j,i) = dump*M4(k,j,i)*rmp
                        C2(k,j,i) = dump*M5(k,j,i)*rmp
                        C3(k,j,i) = dump*M6(k,j,i)*rmp
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
#ifdef WAVEFUNCTION_
            rmp = 1.0_GP/(real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
            DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
               DO j = 1,ny
                  DO k = 1,nz
                     C1(k,j,i) = zre(k,j,i)*rmp
                     C2(k,j,i) = zim(k,j,i)*rmp
                  END DO
               END DO
            END DO
            CALL fftp3d_complex_to_real(plancr,C1,R1,MPI_COMM_WORLD)
            CALL fftp3d_complex_to_real(plancr,C2,R2,MPI_COMM_WORLD)
            CALL io_write(1,odir,'phi_re',ext,planio,R1)
            CALL io_write(1,odir,'phi_im',ext,planio,R2)
            IF (outs.ge.1) THEN
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
               DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
                  DO j = 1,ny
                     DO i = 1,nx
                        R3(i,j,k) = R1(i,j,k)**2+R2(i,j,k)**2
                     END DO
                  END DO
               END DO
               CALL io_write(1,odir,'rho',ext,planio,R3)
            ENDIF
#endif
         ENDIF

#ifdef PART_
         IF ( dolag.GT.0 ) THEN
           IF ((timep.eq.pstep).and.(bench.eq.0)) THEN
#if defined(LAGPART_)
             INCLUDE 'lagpartout_velfield.f90'
#endif
#if defined(MAGFIELD_)
             INCLUDE 'lagpartout_magfield.f90'
#endif
#if defined(SCALAR_)
             INCLUDE 'lagpartout_scalar.f90'
#endif
#if defined(MULTISCALAR_)
             INCLUDE 'lagpartout_mscalar.f90'
#endif
#if defined(INERPART_)
             INCLUDE 'inerpartout_velfield.f90'
#endif
#if defined(TESTPART_) && defined(MAGFIELD_)
             INCLUDE 'testpartout_magfield.f90'
#endif
           ENDIF
         ENDIF
#endif

! Every 'cstep' steps, generates external files 
! with global quantities. If mean=1 also updates 
! the mean fields.

         IF ((timec.eq.cstep).and.(bench.eq.0)) THEN
            timec = 0
            INCLUDE GLOBALOUTPUT_

            IF (mean.eq.1) THEN ! Update mean fields
#ifdef VELOC_
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
               DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
                  DO j = 1,ny
                     DO k = 1,nz
                        M1(k,j,i) = M1(k,j,i)+vx(k,j,i)
                        M2(k,j,i) = M2(k,j,i)+vy(k,j,i)
                        M3(k,j,i) = M3(k,j,i)+vz(k,j,i)
                     END DO
                  END DO
               END DO
#endif
#ifdef SCALAR_
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
               DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
                  DO j = 1,ny
                     DO k = 1,nz
                        M7(k,j,i) = M7(k,j,i)+th(k,j,i)
                     END DO
                  END DO
               END DO
#endif
#ifdef MULTISCALAR_
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
               DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
                  DO j = 1,ny
                     DO k = 1,nz
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
                  DO j = 1,ny
                     DO k = 1,nz
                        M4(k,j,i) = M4(k,j,i)+C1(k,j,i)
                        M5(k,j,i) = M5(k,j,i)+C2(k,j,i)
                        M6(k,j,i) = M6(k,j,i)+C3(k,j,i)
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
            INCLUDE SPECTROUTPUT_
         ENDIF

! Runge-Kutta step 1
! Copies the fields into auxiliary arrays

!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
         DO k = 1,nz

            INCLUDE RKSTEP1_

         END DO
         END DO
         END DO

#ifdef PART_
         IF ( dolag.GT.0 ) THEN
           CALL lagpart%SetStep()
#if defined(INERPART_)
           CALL lagpart%SetStepVel()
#endif
#if defined(TESTPART_) && defined(MAGFIELD_)
           CALL lagpart%SetStepVel()
#endif
         ENDIF
#endif

! Runge-Kutta step 2
! Evolves the system in time

         DO o = ord,1,-1

         INCLUDE RKSTEP2_

#ifdef PART_
         IF ( dolag.GT.0 ) THEN
         rmp = 1.0_GP/(real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
           DO j = 1,ny
             DO k = 1,nz
               C7(k,j,i) = vx(k,j,i)*rmp
             END DO
           END DO
         END DO
         CALL fftp3d_complex_to_real(plancr,C7,R1,MPI_COMM_WORLD)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
           DO j = 1,ny
             DO k = 1,nz
               C7(k,j,i) = vy(k,j,i)*rmp
             END DO
           END DO
         END DO
         CALL fftp3d_complex_to_real(plancr,C7,R2,MPI_COMM_WORLD)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
           DO j = 1,ny
             DO k = 1,nz
               C7(k,j,i) = vz(k,j,i)*rmp
             END DO
           END DO
         END DO
         CALL fftp3d_complex_to_real(plancr,C7,R3,MPI_COMM_WORLD)
#ifdef LAGPART_
         CALL lagpart%Step(R1,R2,R3,dt,1.0_GP/real(o,kind=GP),R4,R5,R6)
#endif

#if defined(INERPART_)
         IF ( dolightp.EQ.0 ) THEN  ! Heavy particles
            CALL lagpart%StepInerp(R1,R2,R3,dt,1.0_GP/real(o,kind=GP),Rv1,Rv2)
         ELSE                       ! Light particles
            CALL gradre3(vx,vy,vz,C4,C5,C6)
            rmp = 1.0_GP/(real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))
            rmq = real(o,kind=GP)/dt
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
            DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
              DO j = 1,ny
                DO k = 1,nz   ! Du_x/Dt = (du/dt + u.grad(u))_x
                  C7(k,j,i) = ((vx(k,j,i)-C1(k,j,i))*rmq + C4(k,j,i))*rmp
                END DO
              END DO
            END DO
            CALL fftp3d_complex_to_real(plancr,C7,R4,MPI_COMM_WORLD)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
            DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
              DO j = 1,ny
                DO k = 1,nz   ! Du_y/Dt = (du/dt + u.grad(u))_y
                  C7(k,j,i) = ((vy(k,j,i)-C2(k,j,i))*rmq + C5(k,j,i))*rmp
                END DO
              END DO
            END DO
            CALL fftp3d_complex_to_real(plancr,C7,R5,MPI_COMM_WORLD)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
            DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
              DO j = 1,ny
                DO k = 1,nz   ! Du_z/Dt = (du/dt + u.grad(u))_z
                  C7(k,j,i) = ((vz(k,j,i)-C3(k,j,i))*rmq + C6(k,j,i))*rmp
                END DO
              END DO
            END DO
            CALL fftp3d_complex_to_real(plancr,C7,R6,MPI_COMM_WORLD)
            CALL lagpart%StepLitep(R1,R2,R3,R4,R5,R6,dt, &
                                   1.0_GP/real(o,kind=GP),Rv1,Rv2)
         ENDIF
#endif

#if defined(TESTPART_) && defined(MAGFIELD_)
         rmp = 1.0_GP/(real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
           DO j = 1,ny
             DO k = 1,nz
               C11(k,j,i) = ax(k,j,i)*rmp
               C12(k,j,i) = ay(k,j,i)*rmp
               C13(k,j,i) = az(k,j,i)*rmp
             END DO
           END DO
         END DO
         CALL rotor3(C12,C13,C14,1) ! bx
         CALL rotor3(C11,C13,C15,2) ! by
         CALL rotor3(C11,C12,C16,3) ! bz
#ifdef UNIFORMB_
         IF (myrank.eq.0) THEN      ! b = b + B_0
            C14(1,1,1) = bx0
            C15(1,1,1) = by0
            C16(1,1,1) = bz0
         ENDIF
#endif
         CALL fftp3d_complex_to_real(plancr,C14,Rb1,MPI_COMM_WORLD)
         CALL fftp3d_complex_to_real(plancr,C15,Rb2,MPI_COMM_WORLD)
         CALL fftp3d_complex_to_real(plancr,C16,Rb3,MPI_COMM_WORLD)
         Rb1 = gyrof*Rb1            ! gyrofreq*bx
         Rb2 = gyrof*Rb2            ! gyrofreq*by
         Rb3 = gyrof*Rb3            ! gyrofreq*bz
         IF ( dokinelv.EQ.1 ) THEN  ! Compute electron velocity corrections
           CALL laplak3(C11,C14)    ! -jx
           CALL laplak3(C12,C15)    ! -jy
           CALL laplak3(C13,C16)    ! -jz
#ifdef COMPRESSIBLE_
           CALL divide(th,C14,C15,C16) ! -j/rho
#endif
           CALL fftp3d_complex_to_real(plancr,C14,Rj1,MPI_COMM_WORLD)
           CALL fftp3d_complex_to_real(plancr,C15,Rj2,MPI_COMM_WORLD)
           CALL fftp3d_complex_to_real(plancr,C16,Rj3,MPI_COMM_WORLD)
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
           DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
             DO j = 1,ny
               DO i = 1,nx
                 R1(i,j,k) = R1(i,j,k) + dii*Rj1(i,j,k) ! u - epsilon.j/rho
                 R2(i,j,k) = R2(i,j,k) + dii*Rj2(i,j,k) ! u - epsilon.j/rho
                 R3(i,j,k) = R3(i,j,k) + dii*Rj3(i,j,k) ! u - epsilon.j/rho
               END DO
             END DO
           END DO
	 ENDIF
         CALL laplak3(C11,C14)      ! -jx (already normalized)
         CALL laplak3(C12,C15)      ! -jy
         CALL laplak3(C13,C16)      ! -jz
         IF ( dokinelp.EQ.0 ) THEN  ! Compute only Ohmic current
           rmq = -gyrof*mu
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
           DO i = ista,iend         ! gyrofreq*(j/Rm-.5*epsilon*grad(p)/rho)
!$omp parallel do if (iend-ista.lt.nth) private (k)
             DO j = 1,ny
               DO k = 1,nz
                 C14(k,j,i) = rmq*C14(k,j,i) ! gyrofreq*j/Rm
                 C15(k,j,i) = rmq*C15(k,j,i)
                 C16(k,j,i) = rmq*C16(k,j,i)
               END DO
             END DO
           END DO
#ifdef COMPRESSIBLE_
         ELSEIF ( dokinelp.EQ.1 ) THEN ! Electron pressure correction
           CALL gradpstate(cp1,gam1,th,C11,C12,C13) ! grad(p_gas)
           CALL divide(th,C11,C12,C13)              ! grad(p_gas)/rho     
           rmq = -gyrof*mu
           rms = .5_GP*gyrof*dii*rmp
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
           DO i = ista,iend         ! gyrofreq*(j/Rm-.5*epsilon*grad(p)/rho)
!$omp parallel do if (iend-ista.lt.nth) private (k)
             DO j = 1,ny
               DO k = 1,nz
                 C14(k,j,i) = rmq*C14(k,j,i) - rms*C11(k,j,i)
                 C15(k,j,i) = rmq*C15(k,j,i) - rms*C12(k,j,i)
                 C16(k,j,i) = rmq*C16(k,j,i) - rms*C13(k,j,i)
               END DO
             END DO
           END DO
#endif
	 ENDIF
         CALL fftp3d_complex_to_real(plancr,C14,Rj1,MPI_COMM_WORLD)
         CALL fftp3d_complex_to_real(plancr,C15,Rj2,MPI_COMM_WORLD)
         CALL fftp3d_complex_to_real(plancr,C16,Rj3,MPI_COMM_WORLD)
         CALL lagpart%StepTestp(R1,R2,R3,Rb1,Rb2,Rb3,Rj1,Rj2,Rj3,dt &
                               ,1.0_GP/real(o,kind=GP),R4,R5,R6)
#endif
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
#ifdef PART_
        IF ( dolag.GT.0 ) THEN
          rbal = rbal + lagpart%GetLoadBal()
        ENDIF
#endif
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
#if defined(PART_)
            IF ( dolag.GT.0 ) THEN
              OPEN(1,file='gpbenchmark.txt',position='append')
              WRITE(1,*) nx,ny,nz,maxparts,rbal/(step-ini+1),            &
                           (step-ini+1),nprocs,nth,                      &
                           lagpart%GetTime   (GPTIME_STEP)/(step-ini+1), &
                           lagpart%GetTime   (GPTIME_COMM)/(step-ini+1), &
                           lagpart%GetTime (GPTIME_SPLINE)/(step-ini+1), &
                           lagpart%GetTime (GPTIME_TRANSP)/(step-ini+1), &
                           lagpart%GetTime (GPTIME_DATAEX)/(step-ini+1), &
                           lagpart%GetTime (GPTIME_INTERP)/(step-ini+1), &
                           lagpart%GetTime(GPTIME_PUPDATE)/(step-ini+1)
              CLOSE(1)
            ENDIF
#endif
         ENDIF
      ENDIF

#if defined(PART_)
! Write particle write times to screen:
      IF ( myrank .EQ. 0 ) THEN
        WRITE(*,*) 'main: <GPTIME_GPWRITE>=', &
	           lagpart%GetTime(GPTIME_GPWRITE)/nwpart
      ENDIF
#endif
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

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

      DEALLOCATE( C1,C2,C3,C4,C5,C6,C7,C8 )
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
#ifdef WAVEFUNCTION_
      DEALLOCATE( zre, zim )
#endif
#ifdef QFORCE_
      DEALLOCATE( fre, fim )
#endif
#ifdef EDQNM_
      DEALLOCATE( C19 )
      DEALLOCATE( tepq,thpq,tve,tvh,Eext,Hext )
#endif
#ifdef PART_
      DEALLOCATE( R4,R5,R6 )
#endif
#if defined(INERPART_)
      DEALLOCATE( Rv1,Rv2 )
#endif
#if defined(TESTPART_) && defined(MAGFIELD_)
      DEALLOCATE( Rb1,Rb2,Rb3 )
      DEALLOCATE( Rj1,Rj2,Rj3 )
#endif
      END PROGRAM MAIN3D
