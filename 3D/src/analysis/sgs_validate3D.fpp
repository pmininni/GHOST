!=================================================================
      PROGRAM MAIN3D
!=================================================================
! GHOST code: Geophysical High Order Suite for Turbulence
!
! Validates SGS data model for ROTBOUSS_SGS solver
!
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
#if defined(PART_) || defined(PIC_)
      USE class_GPart
#endif
#if defined(VELOCSGS_) && defined(SCALARSGS_)
      USE class_GSGSmodel
#endif

      IMPLICIT NONE

!
! Arrays for the fields and external forcings

#if defined(VELOC_) || defined(ADVECT_)
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: vx,vy,vz
#endif
#if defined(MOM_) 
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: sx,sy,sz
#endif
#ifdef VELOC_
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: fx,fy,fz
#endif
#ifdef DENSITY_
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: rho
      CHARACTER                                        :: srho           
#endif
#ifdef SCALAR_
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: th
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: fs
#endif
#ifdef MULTISCALAR_
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: th1,th2,th3
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: fs1,fs2,fs3
#endif
#ifdef CPIC_
      TYPE(ChargPIC)                                   :: picpart
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: rhoc,Temp
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: ux,uy,uz
#endif
#ifdef ELECFIELD_
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: phi
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
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: C1SGS,C2SGS,C3SGS
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: SGS1,SGS2,SGS3
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: LSGS1,LSGS2,LSGS3
#ifdef VELOC_
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: M1,M2,M3
#endif
#ifdef CPIC_
      REAL(KIND=GP)   , ALLOCATABLE, DIMENSION (:,:,:) :: Re1,Re2,Re3
      REAL(KIND=GP)   , ALLOCATABLE, DIMENSION (:,:,:) :: Rb1,Rb2,Rb3
      REAL(KIND=GP)   , ALLOCATABLE, DIMENSION (:,:,:) :: Rj1,Rj2,Rj3
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
#ifdef SCALARSGS_
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: SGSth,LSGSth
      CHARACTER(len=len(ext)) sgsext
#endif
#ifdef MULTISCALAR_
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: C21,C22,C23,C24
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: M8,M9,M10
#endif
#ifdef DENSITY_
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: M11
#endif
#ifdef COMPR_AUX_ARR_
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: C25,C26,C27
#endif
#ifdef TRAP_
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: C28,C29
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: C30
#endif
#ifdef COMPI_AUX_ARR_
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: C31,C32,C33
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: C34,C35,C36
#endif
#ifdef PENALTY_
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: C37,C38
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: C39,C40
#endif

#ifdef WAVEFUNCTION_
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:)     :: iold,qold,kold,cold
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:)     :: inew,qnew,knew,cnew
#endif

#ifdef VELOC_
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: fxold,fyold,fzold
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: fxnew,fynew,fznew
      REAL(KIND=GP), ALLOCATABLE, DIMENSION (:)        :: Faux1,Faux2
#endif
#ifdef VELOCSGS_
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: C1SGS,C2SGS,C3SGS,SGS1,SGS2,SGS3
#endif
#ifdef MAGFIELD_
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: mxold,myold,mzold
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: mxnew,mynew,mznew
#endif

#ifdef PENALTY_
      REAL(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:)    :: chi
#endif
#ifdef ADVECT_
      REAL(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:)    :: vsq
#endif
#ifdef TRAP_
      REAL(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:)    :: Vtrap
      REAL(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:)    :: Vlinx,Vliny
#endif

      REAL(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:)    :: R1,R2,R3
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
      INTEGER          :: Stokeshyp
      REAL(KIND=GP)    :: smach, gam1, cp1, nu2, rho0
#endif
#ifdef CMHD_
      REAL(KIND=GP)    :: amach, cp2
#endif
#ifdef ELECSTAT_
      REAL(KIND=GP)    :: kde,kde2
#endif
#ifdef PIC_
      REAL         :: rbal
      INTEGER      :: splord, picdiv, partpcell
      INTEGER      :: picinittype
      INTEGER      :: picexchtype,picouttype
      INTEGER      :: picwrtunit,piccoll
#endif
#ifdef CPIC_
      REAL(KIND=GP)    :: r0,T0,delT
      REAL(KIND=GP)    :: krd,kru,kud,kuu,ktd,ktu
      REAL(KIND=GP)    :: rparam0,rparam1,rparam2,rparam3,rparam4
      REAL(KIND=GP)    :: rparam5,rparam6,rparam7,rparam8,rparam9
      REAL(KIND=GP)    :: uparam0,uparam1,uparam2,uparam3,uparam4
      REAL(KIND=GP)    :: uparam5,uparam6,uparam7,uparam8,uparam9
      REAL(KIND=GP)    :: tparam0,tparam1,tparam2,tparam3,tparam4
      REAL(KIND=GP)    :: tparam5,tparam6,tparam7,tparam8,tparam9
#endif
#ifdef HYBPIC_
      INTEGER          :: Bmult = 10, iB
      REAL(KIND=GP)    :: ekin,dii,betae,gammae,gam1,cp1,filstr
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
      INTEGER :: i,j,k, it
      INTEGER :: ki,kj,kk
      INTEGER :: pind,tind,sind
      INTEGER :: sindsgs
      INTEGER :: timet,timec
      INTEGER :: times,timef
      INTEGER :: timessgs
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
#ifdef PENALTY_
      INTEGER      :: shape, ischi
      REAL(KIND=GP):: inveta, x0, y0, z0, radius
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

#if defined(DEF_GHOST_CUDA_)
      TYPE(cudaDevicePropG) :: devprop
#endif
#if defined(VELOCSGS_) || defined(SCALARSGS_)
      TYPE (GSGSmodel)      :: mlsgs
      TYPE(GSGSmodelTraits) :: mlsgstraits
      LOGICAL               :: sgs_doproj, sgs_dodealias
      INTEGER               :: sgs_nx, sgs_ny, sgs_nz
      INTEGER               :: sgs_nichannel, sgs_nochannel
      REAL(KIND=GP)         :: sgs_vfactor, sgs_thfactor
      CHARACTER(len=1024)   :: sgs_model_path, sgs_model_type
      CHARACTER(len=1024)   :: sgs_in_name, sgs_out_name
#endif
      TYPE(IOPLAN),TARGET   :: planio

      ! Data specific to app:
      CHARACTER(len=100)    :: odir, idir
      REAL(kind=GP) omega(3),xnormn
      REAL(kind=GP) fmin,fmax
      INTEGER :: ic,ir,it,jc
      INTEGER :: istat(4096),nstat,prtbin,doSGSinj
      INTEGER :: nbinx,nbiny,nbins(2)
      CHARACTER(len=64) :: ext1
      CHARACTER(len=4096) :: sstat
      CHARACTER(len=256) :: fnout


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
#if defined(VELOCSGS_) && defined(SCALARSGS_)
      NAMELIST / mlsgsnml / sgs_doproj, sgs_dodealias, sgs_nx, sgs_ny, &
                            sgs_nz, sgs_vfactor, sgs_thfactor
      NAMELIST / mlsgsnml / sgs_nichannel, sgs_nochannel 
      NAMELIST / mlsgsnml / sgs_model_path, sgs_model_type
      NAMELIST / mlsgsnml / sgs_in_name, sgs_out_name
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
      NAMELIST / compressible / Stokeshyp, smach, gam1, nu2, rho0
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
#ifdef PENALTY_ 
      NAMELIST / penalty / shape,inveta,x0,y0,z0,radius
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

      NAMELIST / validate_mlsgs / iswap,oswap
      NAMELIST / validate_mlsgs/ idir,odir,sstat
      NAMELIST / validate_mlsgs/ nbinx,nbiny,prtbin,doSGSinj

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

     use_mlsgs = .FALSE. ! ML-SGS terms

     CALL range(1,nx/2+1,nprocs,myrank,ista,iend)
     CALL range(1,nz,nprocs,myrank,ksta,kend)
     CALL io_init(myrank,(/nx,ny,nz/),ksta,kend,planio)

!
! Allocates memory for distributed blocks

      ALLOCATE( C1(nz,ny,ista:iend),  C2(nz,ny,ista:iend) )
      ALLOCATE( C3(nz,ny,ista:iend),  C4(nz,ny,ista:iend) )
      ALLOCATE( C5(nz,ny,ista:iend),  C6(nz,ny,ista:iend) )
      ALLOCATE( C7(nz,ny,ista:iend),  C8(nz,ny,ista:iend) )
      ALLOCATE( C1SGS(nz,ny,ista:iend), C2SGS(nz,ny,ista:iend) )
      ALLOCATE( C3SGS(nz,ny,ista:iend), SGS1(nz,ny,ista:iend) )
      ALLOCATE( SGS2(nz,ny,ista:iend), SGS3(nz,ny,ista:iend) )
      ALLOCATE( SGSth(nz,ny,ista:iend) )
      ALLOCATE( LSGS1(nz,ny,ista:iend), LSGS2(nz,ny,ista:iend) )
      ALLOCATE( LSGS3(nz,ny,ista:iend), LSGSth(nz,ny,ista:iend) )
#if defined(VELOC_) || defined(ADVECT_)
      ALLOCATE( vx(nz,ny,ista:iend) )
      ALLOCATE( vy(nz,ny,ista:iend) )
      ALLOCATE( vz(nz,ny,ista:iend) )
#endif
#ifdef VELOCSGS_
      ALLOCATE( SGS1(nz,ny,ista:iend) )
      ALLOCATE( SGS2(nz,ny,ista:iend) )
      ALLOCATE( SGS3(nz,ny,ista:iend) )
      ALLOCATE( C1SGS(nz,ny,ista:iend) )
      ALLOCATE( C2SGS(nz,ny,ista:iend) )
      ALLOCATE( C3SGS(nz,ny,ista:iend) )
#endif

#if defined(MOM_) 
      ALLOCATE( sx(nz,ny,ista:iend) )
      ALLOCATE( sy(nz,ny,ista:iend) )
      ALLOCATE( sz(nz,ny,ista:iend) )
#endif
#ifdef VELOC_
      ALLOCATE( fx(nz,ny,ista:iend) )
      ALLOCATE( fy(nz,ny,ista:iend) )
      ALLOCATE( fz(nz,ny,ista:iend) )
#endif
#ifdef CPIC_
      ALLOCATE(rhoc(nz,ny,ista:iend))
      ALLOCATE(  ux(nz,ny,ista:iend))
      ALLOCATE(  uy(nz,ny,ista:iend))
      ALLOCATE(  uz(nz,ny,ista:iend))
      ALLOCATE(Temp(nz,ny,ista:iend))
#endif
#ifdef ELECFIELD_
      ALLOCATE( phi(nz,ny,ista:iend))
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
#ifdef DENSITY_
      ALLOCATE( rho(nz,ny,ista:iend) )
      srho = 'rhospect'
#endif
#ifdef SCALAR_
      ALLOCATE( C20(nz,ny,ista:iend) )
      ALLOCATE( th (nz,ny,ista:iend) )
      ALLOCATE( fs (nz,ny,ista:iend) )
#endif
#ifdef SCALARSGS_
      ALLOCATE( SGSth(nz,ny,ista:iend) )
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
#ifdef COMPI_AUX_ARR_
      ALLOCATE( C31(nz,ny,ista:iend) )
      ALLOCATE( C32(nz,ny,ista:iend) )
      ALLOCATE( C33(nz,ny,ista:iend) )
      ALLOCATE( C34(nz,ny,ista:iend) )
      ALLOCATE( C35(nz,ny,ista:iend) )
      ALLOCATE( C36(nz,ny,ista:iend) )
#endif
#ifdef PENALTY_
      ALLOCATE( C37(nz,ny,ista:iend),  C38(nz,ny,ista:iend) )
      ALLOCATE( C39(nz,ny,ista:iend),  C40(nz,ny,ista:iend) )
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
      ALLOCATE( Hinv(nz,ny,ista:iend) ) ! Voigt inv. Helmholtz operator

      ALLOCATE( R1(nx,ny,ksta:kend) )
      ALLOCATE( R2(nx,ny,ksta:kend) )
      ALLOCATE( R3(nx,ny,ksta:kend) )
#ifdef PENALTY_
      ALLOCATE( chi(nx,ny,ksta:kend) )
#endif
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
#if defined (CPIC_)
      ALLOCATE( Rb1(nx,ny,ksta:kend) )
      ALLOCATE( Rb2(nx,ny,ksta:kend) )
      ALLOCATE( Rb3(nx,ny,ksta:kend) )
      ALLOCATE( Re1(nx,ny,ksta:kend) )
      ALLOCATE( Re2(nx,ny,ksta:kend) )
      ALLOCATE( Re3(nx,ny,ksta:kend) )
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

#if defined(VELOC_) || defined(ADVECT_) ||  defined(VELOCSGS_) 
! Reads parameters for the velocity field from the
! namelist 'velocity' on the external file 'parameter.inp'
!     f0   : amplitude of the mechanical forcing
!     u0   : amplitude of the initial velocity field
!     kdn  : minimum wave number in v/mechanical forcing
!     kup  : maximum wave number in v/mechanical forcing
!     nu   : kinematic viscosity
!     fparam0-9   : ten real numbers to control properties of
!                   the mechanical forcing
!     vparam0-9   : ten real numbers to control properties of
!                   the initial conditions for the velocity field

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

#if defined(VELOCSGS_) && defined(SCALARSGS_) 
! Reads parameters for the SGS ML model from the 
! namelist 'mlsgs' on the external file 'parameter.inp'
!     do_projection: do SGS projection?
!     do_dealias : do explicit dealiasing of SGS terms?
!     nx,ny,nz   : sizes that model thinks it has
!     nchannel   : no. channels/features
!     model_path : path to model (+name)
!     model_type : model type (required by ONNx)
!     in_name    : input tensor name
!     out_name   : output tensor name
       
      sgs_doproj     = .true.
      sgs_dodealias   = .true.
      sgs_nx         = nx
      sgs_ny         = ny
      sgs_nz         = nz
      sgs_vfactor    = 1.0
      sgs_thfactor   = 1.0
      sgs_nichannel  = 4
      sgs_nochannel  = 4
      sgs_model_path = './mymodel'
      sgs_model_type = 'CNN'
      sgs_in_name    = 'state'
      sgs_out_name   = 'SGS'

      IF (myrank.eq.0) THEN
         OPEN(1,file='parameter.inp',status='unknown',form="formatted")
         READ(1,NML=mlsgsnml)
         CLOSE(1)
      ENDIF
      CALL MPI_BCAST(sgs_doproj    ,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(sgs_dodealias ,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(sgs_nx        ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(sgs_ny        ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(sgs_nz        ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(sgs_nichannel ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(sgs_nochannel ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(sgs_vfactor   ,1,GC_REAL    ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(sgs_thfactor  ,1,GC_REAL    ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(sgs_model_path,1024,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(sgs_model_type,1024,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(sgs_in_name   ,1024,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(sgs_out_name  ,1024,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)

      mlsgstraits%do_projection = sgs_doproj
      mlsgstraits%do_dealias = sgs_dodealias
      mlsgstraits%nx         = sgs_nx
      mlsgstraits%ny         = sgs_ny
      mlsgstraits%nz         = sgs_nz
      mlsgstraits%vfactor    = sgs_vfactor
      mlsgstraits%thfactor   = sgs_thfactor
      mlsgstraits%nichannel  = sgs_nichannel
      mlsgstraits%nochannel  = sgs_nochannel
      mlsgstraits%model_path = sgs_model_path
      mlsgstraits%model_type = sgs_model_type
      mlsgstraits%in_name    = sgs_in_name
      mlsgstraits%out_name   = sgs_out_name
      mlsgstraits%odir       = odir
      mlsgstraits%idir       = idir
      mlsgstraits%planio     => planio
      write(*,*)'            main: sgs_nx=',sgs_nx,' sgs_ny=',sgs_ny, ' sgs_nz=', sgs_nz
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

#ifdef PENALTY_
! Reads parameters for runs with the penalty method from the
! namelist 'penalty' on the external file 'parameter.inp'
!     shape : shape of the obstacle
!     inveta: inverse of the eta penalty parameter
!     x0    : x coordinate of the center of the obstacle
!     y0    : y coordinate of the center of the obstacle
!     z0    : z coordinate of the center of the obstacle
!     radius: radius of the obstacle

      inveta = 1000.0_GP
      ischi  = 0 ! The penalty function must be computed after a restart
      IF (myrank.eq.0) THEN
         OPEN(1,file='parameter.inp',status='unknown',form="formatted")
         READ(1,NML=penalty)
         CLOSE(1)
      ENDIF
      CALL MPI_BCAST(shape,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(inveta,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(x0    ,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(y0    ,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(z0    ,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(radius,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
#endif

#ifdef COMPRESSIBLE_
! Reads parameters for the compressible runs from the
! namelist 'compressible' on the external file 'parameter.inp'
!     smach    : sound Mach number
!     gam1     : gamma parameter for polytropic eq. of state
!     nu2      : second (bulk) viscosity for divergence (velocity) term
!     rho0     : reference density
!     Stokeshyp: if 1, then nu2 = -2/3 * nu; else set by user
      Stokeshyp = 0
      rho0      = 1.0

      IF (myrank.eq.0) THEN
         OPEN(1,file='parameter.inp',status='unknown',form="formatted")
         READ(1,NML=compressible)
         CLOSE(1)
      ENDIF
      CALL MPI_BCAST(smach    ,1,GC_REAL    ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(gam1     ,1,GC_REAL    ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(Stokeshyp,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(nu2      ,1,GC_REAL    ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(rho0     ,1,GC_REAL    ,0,MPI_COMM_WORLD,ierr)
      gam1 = gam1 - 1.0_GP
      cp1  = 2.0_GP / (gam1*smach*smach)
      IF ( Stokeshyp .GE. 1 ) THEN
        ! s.t. nu2 + 2 nu/d = 0 ; d == dimensionality
        nu2  = -2.0_GP * nu / 3.0_GP
      ELSE IF (Stokeshyp .EQ. 0)  THEN
        nu2  = nu2 + nu/3.0
      ! IF < 0, then accept nu2 as set
      ENDIF
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

#ifdef HYBPIC_
! namelist 'phybrid' on the external file 'parameter.inp'
!     gammae: barotropic exponent for fluid electrons
!     betae : electronic plasma beta
!     dii   : ion inertial length scale
!     Bmult : magnetic field steps per ion step

      gammae = 0.0_GP
      dii    = 1.0_GP
      IF (myrank.eq.0) THEN
         OPEN(1,file='parameter.inp',status='unknown',form="formatted")
         READ(1,NML=phybrid)
         CLOSE(1)
      ENDIF
      CALL MPI_BCAST(gammae ,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(betae  ,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(dii    ,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(filstr,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(Bmult,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      gam1 = gammae-1.0_GP
      cp1 = dii*betae/2
      IF (gammae.GT.1) THEN
         cp1 = cp1*gammae/gam1
      END IF
#endif

#ifdef ELECSTAT_
! namelist 'elecstat' on the external file 'parameter.inp'
!     kde  : inverse debye length
      
      kde  = 0.0_GP
      IF (myrank.eq.0) THEN
         OPEN(1,file='parameter.inp',status='unknown',form="formatted")
         READ(1,NML=elecstat)
         CLOSE(1)
      ENDIF
      CALL MPI_BCAST(kde ,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      kde2 = kde*kde
#endif

#ifdef PIC_
! namelist 'ppic' on the external file 'parameter.inp'
!     splord     : spline order (0, 1, 2 or 3)
!     partpcell  : number of particles per cell
!     picdiv     : divisor for pic particle output
!     picinittype: initialization type for locations (0=random,1=lattice, 2=user)
!     picexchtype: boundary exchange type (0=nearest neighbor, 1=voxeldb)
!     picouttype : output type (0=binary, 1=ASCII)
!     piccoll    : I/O method when using binary output (0=task 0,1=collective)
!     picwrtunit :   ! Write part. positions in box units (=1) or grid units(=0)
!     picseedfile: filename for picinittype=2 case

      splord    = 0
      picdiv    = 1
      partpcell = 1
      rbal      = 0.0
      spicfpfile= 'xpicInitRndSeed.000.txt'
      IF (myrank.eq.0) THEN
         OPEN(1,file='parameter.inp',status='unknown',form="formatted")
         READ(1,NML=ppic)
         CLOSE(1)
      ENDIF
      CALL MPI_BCAST(splord,        1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(picdiv,        1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(partpcell,     1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(picinittype   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(picexchtype   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(picouttype    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(picwrtunit    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(piccoll       ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(spicfpfile ,1024,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(picseedfile,1024,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
#endif

#ifdef CPIC_
! namelist 'picinitflds' on the external file 'parameter.inp'
!     r0        : amplitude of the density fluctuations (mean density=1)
!     krd       : minimum wavenumber for density fluctuations
!     kru       : maximum wavenumber for density fluctuations
!     rparam0   : rparam0-9 can be used to control the ion density
!     rparam0-9 : ten real numbers to control properties of
!                 the initial ion density fluctuations
!     u0        : amplitude of the initial ion mean velocity field
!     kud       : minimum wavenumber for ion mean velocity
!     kuu       : maximum wavenumber for ion mean velocity
!     uparam0-9 : ten real numbers to control properties of
!                 the initial ion mean velocity fluctuations
!     T0        : mean temperature of the ions
!     delT      : amplitude of the temperature fluctuations
!     ktd       : minimum wavenumber for temperature fluctuations
!     ktu       : maximum wavenumber for temperature fluctuations
!     tparam0-9 : ten real numbers to control properties of
!                 the initial ion temperature fluctuations

      IF (myrank.eq.0) THEN
         OPEN(1,file='parameter.inp',status='unknown',form="formatted")
         READ(1,NML=picinitflds)
         CLOSE(1)
      ENDIF
      CALL MPI_BCAST(r0  ,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(krd ,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(kru ,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(u0  ,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(kud ,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(kuu ,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(T0  ,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(delT,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(ktd ,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(ktu ,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(rparam0,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(rparam1,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(rparam2,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(rparam3,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(rparam4,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(rparam5,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(rparam6,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(rparam7,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(rparam8,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(rparam9,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(uparam0,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(uparam1,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(uparam2,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(uparam3,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(uparam4,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(uparam5,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(uparam6,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(uparam7,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(uparam8,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(uparam9,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(tparam0,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(tparam1,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(tparam2,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(tparam3,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(tparam4,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(tparam5,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(tparam6,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(tparam7,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(tparam8,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(tparam9,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
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

! Reads from the external file 'validate_mlsgs.inp' the 
! parameters that will be used to compute the transfer
!     idir   : directory for unformatted input (field components)
!     odir   : directory for unformatted output (prolongated data)
!     sstat  : time index for which to compute VOIGT, or a
!     ';--separated list
!     iswap  : do endian swap on input?
!     oswap  : do endian swap on output? Not used.
!
!     Defaults:
      idir   = '.' ! location of label data: vx,vy,vz, SGS*
      odir   = '.' ! location of all output from this app
      sstat  = '0'
      iswap  = 0
      oswap  = 0
      nbinx  = 100
      nbiny  = 100
      prtbin = 0   ! don't print binary data
      doSGSinj = 0 ! don't examine SGSinj terms
      doSGSinj = 0 ! don't examine SGSinj terms

      IF (myrank.eq.0) THEN
         OPEN(1,file='validate_mlsgs.inp',status='unknown',form="formatted")
         READ(1,NML=validate_mlsgs)
         CLOSE(1)
      ENDIF
      CALL MPI_BCAST(idir     ,1024,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(odir     ,1024,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(sstat    ,4096,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(oswap    ,1   ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(iswap    ,1   ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(nbinx    ,1   ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(nbiny    ,1   ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(prtbin   ,1   ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(doSGSinj ,1   ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)


! Before continuing, we verify that all parameters and compilation
! options are compatible with the SOLVER being used

#if !(defined(ROTBOUSS_SGS_SOL))
     write(*,*) 'main: This app must use the ROTBOUSS_SGS_SOL solver'
     STOP
#endif

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
      IF (anis.eq.1) THEN
         rmp = 1.0_GP/real(nx,kind=GP)**2
         rmq = 1.0_GP/real(ny,kind=GP)**2
         rms = 1.0_GP/real(nz,kind=GP)**2
      ELSE
         rmp = 1.0_GP
	 rmq = 1.0_GP
	 rms = 1.0_GP
      ENDIF

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

     ! Create ML-SGS interface:
#if defined(ROTBOUSS_SGS_SOL)
      CALL mlsgs%GSGS_ctor(MPI_COMM_WORLD, (/nx ,ny ,nz /),&
           (/ista,iend,ksta,kend/), anis, (/Dkx,Dky,Dkz/), &
           plancr, planrc , mlsgstraits)
#endif

      ! Parse input index set, store in istat:
      CALL parseind(sstat,';', istat , 4096, nstat)

      rmp = 1.0_GP/ &
      (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))

      DO it = 1,nstat
        WRITE(ext1, fmtext) istat(it)
        ext = trim(ext1)
        nbins(1) = nbinx ; nbins(2) = nbiny

if (myrank.eq.0) write(*,*)'main: Reading time index: ', ext, '...'
#ifdef MOM_
if (myrank.eq.0) write(*,*)'main: Reading vx_T...'
        CALL io_read(1,idir,'vx_T',ext,planio,R1)
        CALL fftp3d_real_to_complex(planrc,R1,vx,MPI_COMM_WORLD)
if (myrank.eq.0) write(*,*)'main: Reading vy_T...'
        CALL io_read(1,idir,'vy_T',ext,planio,R1)
        CALL fftp3d_real_to_complex(planrc,R1,vy,MPI_COMM_WORLD)
if (myrank.eq.0) write(*,*)'main: Reading vz_T...'
        CALL io_read(1,idir,'vz_T',ext,planio,R1)
        CALL fftp3d_real_to_complex(planrc,R1,vz,MPI_COMM_WORLD)
if (myrank.eq.0) write(*,*)'main: Reading th_T...'
        CALL io_read(1,idir,'th_T',ext,planio,R1)
        CALL fftp3d_real_to_complex(planrc,R1,th,MPI_COMM_WORLD)

        ! Get label SGS terms:
if (myrank.eq.0) write(*,*)'main: Reading LSGS1_T...'
        CALL io_read(1,idir,'SGS1_T',ext,planio,R1)
        CALL fftp3d_real_to_complex(planrc,R1,LSGS1,MPI_COMM_WORLD)
if (myrank.eq.0) write(*,*)'main: Reading LSGS2_T...'
        CALL io_read(1,idir,'SGS2_T',ext,planio,R1)
        CALL fftp3d_real_to_complex(planrc,R1,LSGS2,MPI_COMM_WORLD)
if (myrank.eq.0) write(*,*)'main: Reading LSGS3_T...'
        CALL io_read(1,idir,'SGS3_T',ext,planio,R1)
        CALL fftp3d_real_to_complex(planrc,R1,LSGS3,MPI_COMM_WORLD)
if (myrank.eq.0) write(*,*)'main: Reading LSGSth_T...'
        CALL io_read(1,idir,'SGSth_T',ext,planio,R1)
        CALL fftp3d_real_to_complex(planrc,R1,LSGSth,MPI_COMM_WORLD)

        ! Do analysis: Compute model SGS terms:
        CALL mlsgs%sgs_model(vx, vy, vz, th, &
                             C1SGS, C2SGS, C3SGS, R1, &
                             SGS1 , SGS2 , SGS3 , SGSth)

        CALL DOCOMPARE(SGS1 , SGS2 ,SGS3 ,SGSth , &
                       LSGS1, LSGS2,LSGS4,LSGSth, &
                       istat(it), odir, nbins,   &
                       C1, C2, R1, R2)


        if (myrank.eq.0) write(*,*)'main: Time index ', ext, ' done.'
      ENDDO ! end, it-loop

!
!       CALL mlsgs%prtbininject(vx, vy, vz, th, &
!                               SGS1, SGS2, SGS3, SGSth, ext, &
!                               C1SGS , C2SGS , R1, R2, R3)


      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

!
! End of MAIN3D

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
      DEALLOCATE( Hinv)
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
#ifdef VELOCSGS_
      DEALLOCATE( SGS1 ,SGS2 ,SGS3 ) 
      DEALLOCATE( C1SGS,C2SGS,C3SGS) 
#endif
#ifdef MOM_ 
      DEALLOCATE( sx,sy,sz )
#endif
#ifdef ADVECT_
      DEALLOCATE( vsq )
#endif
#ifdef CPIC_
      DEALLOCATE(rhoc)
      DEALLOCATE( ux )
      DEALLOCATE( uy )
      DEALLOCATE( uz )
      DEALLOCATE(Temp)
#endif
#ifdef ELECFIELD_
      DEALLOCATE(phi)
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
#ifdef DENSITY_
      DEALLOCATE( rho )
      DEALLOCATE( M11 )
#endif
#ifdef SCALAR_
      DEALLOCATE( th,fs )
      DEALLOCATE( C20 )
      IF (mean.eq.1) DEALLOCATE( M7 )
#endif
#ifdef SCALARSGS_
      DEALLOCATE( SGSth )
#endif
#ifdef MULTISCALAR_
      DEALLOCATE( th1,fs1,th2,fs2,th3,fs3 )
      DEALLOCATE( C21,C22,C23,C24 )
      IF (mean.eq.1) DEALLOCATE( M8,M9,M10 )
#endif
#ifdef COMPR_AUX_ARR_
      DEALLOCATE( C25,C26,C27 )
#endif
#ifdef COMPI_AUX_ARR_
      DEALLOCATE( C31,C32,C33,C34,C35,C36 )
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
#ifdef PENALTY_
      DEALLOCATE( C31, C32, C33, C34, chi )
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
      DEALLOCATE( C1SGS,C2SGS,C3SGS,SGS1,SGS2,SGS3,SGSth)
      DEALLOCATE( LSGS1,LSGS2,LSGS3,LSGSth)
      END PROGRAM MAIN3D


      SUBROUTINE DOCOMPARE(SGS1 ,SGS2 ,SGS3 ,SGSth , &
                           LSGS1,LSGS2,LSGS4,LSGSth, &
                           indtime, odir, nbins, C1, C2, R1, R2) 
!-----------------------------------------------------------------
!-----------------------------------------------------------------
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
      USE iovar
      USE filefmt
      USE boxsize
 
        IMPLICIT NONE

        COMPLEX(KIND=GP), INTENT   (IN), DIMENSION(nz,ny,ista:iend):: SGS1,SGS2,SGS3,SGSth
        COMPLEX(KIND=GP), INTENT   (IN), DIMENSION(nz,ny,ista:iend):: LSGS1,LSGS2,LSGS3,LSGSth
        COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(nx,ny,ista:iend):: C1,C2
        REAL   (KIND=GP), INTENT(INOUT), DIMENSION(nx,ny,ksta:kend):: R1,R2

        INTEGER         , INTENT   (IN)                            :: indtime
        INTEGER         , INTENT   (IN)                            :: nbins(2)

        LOGICAL              :: bexist
        INTEGER              :: ierr, i, kz, n
        REAL   (KIND=GP)     :: g5,s2,s3,s4,s5,s6,tmp,w6
        REAL   (KIND=GP)     :: av(10),sk(10),ku(10),var(10)
        REAL   (KIND=GP)     :: avL(10),skL(10),kuL(10),varL(10)
        REAL   (KIND=GP)     :: sgs1_corr(10), sgs2_corr(10), &
                                sgs3_corr(10), sgsth_corr(10)
        REAL   (KIND=GP)     :: fmin(2)

        CHARACTER(len=1024)  :: fnout
        CHARACTER(len=1024), INTENT(IN)  :: odir
        CHARACTER(len=128)   :: sfld(10)
        CHARACTER(len=128)   :: hdrfmt, rowfmt



        knz    = kend - ksta + 1
        tmp    = 1.0_GP/ ( real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP) )

        C1 = SGS1; ! SGS1
        C1 = C1 * tmp 
        CALL fftp3d_complex_to_real(plancr,C1,R1,MPI_COMM_WORLD)
        n     = n + 1; sfld(n) = 'sgs1';
        CALL skewflat(R1,nx,ny,knz,av(n),sk(n),ku(n),g5,w6,var(n),s3,s4,s5,s6)
        fnout = trim(odir) // '/' // 'SGS1pdf.' // ext // '.txt'
        CALL dopdfr(R1,nx,ny,knz,fnout,nbins(1),0,fmin(1),fmax(1),0)

        C1 = LSGS1;
        C1 = C1 * tmp 
        CALL fftp3d_complex_to_real(plancr,C1,R2,MPI_COMM_WORLD)
        fnout = trim(odir) // '/' // 'SGS1Lpdf.' // ext // '.txt'
        CALL skewflat(R2,nx,ny,knz,avL(n),skL(n),kuL(n),g5,w6,varL(n),s3,s4,s5,s6)
        CALL dopdfr(R2,nx,ny,knz,fnout,nbins(1),0,fmin(1),fmax(1),0)
        sgs1_corr(n) = compute_corr(R1,R2)

        C1 = SGS2;    ! SGS2
        C1 = C1 * tmp 
        CALL fftp3d_complex_to_real(plancr,C1,R1,MPI_COMM_WORLD)
        n     = n + 1; sfld(n) = 'sgs2';
        CALL skewflat(R1,nx,ny,knz,avL(n),skL(n),kuL(n),g5,w6,var(n),s3,s4,s5,s6)
        fnout = trim(odir) // '/' // 'SGS2pdf.' // ext // '.txt'
        CALL dopdfr(R1,nx,ny,knz,fnout,nbins(1),0,fmin(1),fmax(1),0)

        C1 = LSGS2;
        C1 = C1 * tmp 
        CALL fftp3d_complex_to_real(plancr,C1,R2,MPI_COMM_WORLD)
        fnout = trim(odir) // '/' // 'SGS2Lpdf.' // ext // '.txt'
        CALL skewflat(R2,nx,ny,knz,avL(n),skL(n),kuL(n),g5,w6,varL(n),s3,s4,s5,s6)
        CALL dopdfr(R2,nx,ny,knz,fnout,nbins(1),0,fmin(1),fmax(1),0)
        sgs2_corr(n) = compute_corr(R1,R2)

        C1 = SGS3;    ! SGS3
        C1 = C1 * tmp 
        CALL fftp3d_complex_to_real(plancr,C1,R1,MPI_COMM_WORLD)
        n  = n + 1; sfld(n) = 'sgs3';
        CALL skewflat(R1,nx,ny,knz,av(n),sk(n),ku(n),g5,w6,var(n),s3,s4,s5,s6)
        fnout = trim(odir) // '/' // 'SGS3pdf.' // ext // '.txt'
        CALL dopdfr(R1,nx,ny,knz,fnout,nbins(1),0,fmin(1),fmax(1),0)

        C1 = LSGS3;
        C1 = C1 * tmp 
        CALL fftp3d_complex_to_real(plancr,C1,R2,MPI_COMM_WORLD)
        fnout = trim(odir) // '/' // 'SGS3Lpdf.' // ext // '.txt'
        CALL skewflat(R2,nx,ny,knz,avL(n),skL(n),kuL(n),g5,w6,varL(n),s3,s4,s5,s6)
        CALL dopdfr(R2,nx,ny,knz,fnout,nbins(1),0,fmin(1),fmax(1),0)
        sgs3_corr(n) = compute_corr(R1,R2)

        C1 = SGSth;    ! SGSth
        C1 = C1 * tmp 
        CALL fftp3d_complex_to_real(plancr,C1,R1,MPI_COMM_WORLD)
        n     = n + 1; sfld(n) = 'sgsth';
        CALL skewflat(R1,nx,ny,knz,av(n),sk(n),ku(n),g5,w6,var(n),s3,s4,s5,s6)
        fnout = trim(odir) // '/' // 'SGSthpdf.' // ext // '.txt'
        CALL dopdfr(R1,nx,ny,knz,fnout,nbins(1),0,fmin(1),fmax(1),0)

        C1 = LSGSth;
        C1 = C1 * tmp 
        CALL fftp3d_complex_to_real(plancr,C1,R2,MPI_COMM_WORLD)
        fnout = trim(odir) // '/' // 'SGSthLpdf.' // ext // '.txt'
        CALL skewflat(R2,nx,ny,knz,avL(n),skL(n),kuL(n),g5,w6,varL(n),s3,s4,s5,s6)
        CALL dopdfr(R2,nx,ny,knz,fnout,nbins(1),0,fmin(1),fmax(1),0)
        sgsth_corr(n) = compute_corr(R1,R2)

        ! Write data to files:
        IF ( myrank.EQ.0 ) THEN
          inquire( file='var.txt', exist=bexist )

          ! Create format for statistical data:
          WRITE(rowfmt,'(A, I4, A)') '(I4,',n,'(2X,E14.6))'
          WRITE(hdrfmt,'(A, I4, A)') '(A,',n,'(2X,A))'
  
          fnout = trim(odir) // '/' // 'var_inference.txt'
          OPEN(2,file=trim(fnout),position='append')
          if ( .NOT. bexist ) THEN
          WRITE(2,hdrfmt,advance='yes') '#itime', (sfld(j), j=1,n)
          ENDIF
          WRITE(2,rowfmt,advance='no') indtime, (var(j), j=1,n)
          CLOSE(2)
  
          fnout = trim(odir) // '/' // 'var_label.txt'
          OPEN(2,file=trim(fnout),position='append')
          if ( .NOT. bexist ) THEN
          WRITE(2,hdrfmt,advance='yes') '#itime', (sfld(j), j=1,n)
          ENDIF
          WRITE(2,rowfmt,advance='no') indtime, (varL(j), j=1,n)
          CLOSE(2)
  
          fnout = trim(odir) // '/' // 'skew_inference.txt'
          OPEN(2,file=trim(fnout),position='append')
          if ( .NOT. bexist ) THEN
          WRITE(2,hdrfmt,advance='yes') '#itime', (sfld(j), j=1,n)
          ENDIF
          WRITE(2,rowfmt,advance='no') indtime, (sk(j), j=1,n)
          CLOSE(2)

          fnout = trim(odir) // '/' // 'skew_label.txt'
          OPEN(2,file=trim(fnout),position='append')
          if ( .NOT. bexist ) THEN
          WRITE(2,hdrfmt,advance='yes') '#itime', (sfld(j), j=1,n)
          ENDIF
          WRITE(2,rowfmt,advance='no') indtime, (skL(j), j=1,n)
          CLOSE(2)
  
          fnout = trim(odir) // '/' // 'flat_inference.txt'
          OPEN(2,file=trim(fnout),position='append')
          if ( .NOT. bexist ) THEN
          WRITE(2,hdrfmt,advance='yes') '#itime', (sfld(j), j=1,n)
          ENDIF
  
          fnout = trim(odir) // '/' // 'flat_label.txt'
          OPEN(2,file=trim(fnout),position='append')
          if ( .NOT. bexist ) THEN
          WRITE(2,hdrfmt,advance='yes') '#itime', (sfld(j), j=1,n)
          ENDIF
          WRITE(2,rowfmt,advance='no') indtime, (kuL(j), j=1,n)
          CLOSE(2)
  
          WRITE(rowfmt,'(A, I4, A)') '(I4,',4,'(2X,E14.6))'
          WRITE(hdrfmt,'(A, I4, A)') '(A,',4,'(2X,A))'

          fnout = trim(odir) // '/' // 'sgs_corr.txt'
          OPEN(2,file=trim(fnout),position='append')
          if ( .NOT. bexist ) THEN
          WRITE(2,hdrfmt,advance='yes') '#itime', 'sgs1', 'sgs2', 'sgs3', 'sgsth'
          ENDIF
          DO j = 1, n
          WRITE(2,rowfmt,advance='no') indtime, sgs1_corr(j),sgs2_corr(j), sgs3_corr(j), sgsth_corr(j)
          ENDDO
          CLOSE(2)
      ENDIF

      END SUBROUTINE 


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
!     nz    : local dimensions
!     avg   : average
!     skew  : skewness
!     flat  : flatness/kurtosis
!     glop  : Sum (x-<x>)^5/[Sum(x-<x>)^2]^5/2, 5th order moment ('glop')
!     whoa  : Sum (x-<x>)^6/[Sum(x-<x>)^2]^3, 6th order moment ('whoa')
!     s2-s6 : 2nd-6th order moments: sum( (x-<x>)^n ), where n=2-6.
!-----------------------------------------------------------------
      USE fprecision
      USE commtypes
      IMPLICIT NONE

      REAL(KIND=GP), INTENT(INOUT), DIMENSION(*)   :: fx
      REAL(KIND=GP), INTENT  (OUT)                 :: skew,flat,glop,whoa
      REAL(KIND=GP), INTENT  (OUT)                 :: avg,s2,s3,s4,s5,s6
      DOUBLE PRECISION                             :: davg,ds2,ds3,ds4,ds5,ds6
      DOUBLE PRECISION                             :: gs(5),s(5),xnorm
      INTEGER      , INTENT   (IN)                 :: nx,ny,nz
      INTEGER                                      :: gnz,ierr,j

      INTEGER  nin

      nin = nx * ny * nz
      CALL MPI_ALLREDUCE(nz, gnz, 1, MPI_INTEGER, &
                      MPI_SUM, MPI_COMM_WORLD,ierr)

      xnorm = 1.0_GP / (dble(nx)*dble(ny)*dble(gnz))
      ds2 = 0.0D0
!$omp parallel do default(shared) private(j) reduction(+:s2)
      DO j = 1, nin
        ds2 = ds2 + dble(fx(j))
      ENDDO
!$omp end parallel do

      CALL MPI_ALLREDUCE(ds2, davg, 1, MPI_DOUBLE_PRECISION, &
                      MPI_SUM, MPI_COMM_WORLD,ierr)
      davg = davg * xnorm
      avg  = davg

      ds2 = 0.0D0
!$omp parallel do default(shared) private(j) reduction(+:s2)
      DO j = 1, nin
        ds2 = ds2 + (dble(fx(j))-davg)**2
      ENDDO
!$omp end parallel do

      ds3 = 0.0D0
!$omp parallel do default(shared) private(j) reduction(+:s3)
      DO j = 1, nin
        ds3 = ds3 + (dble(fx(j))-davg)**3
      ENDDO
!$omp end parallel do

      ds4 = 0.0D0
!$omp parallel do default(shared) private(j) reduction(+:s4)
      DO j = 1, nin
        ds4 = ds4 + (dble(fx(j))-davg)**4
      ENDDO
!$omp end parallel do

      ds5 = 0.0D0
!$omp parallel do default(shared) private(j) reduction(+:s5)
      DO j = 1, nin
        ds5 = ds5 + (dble(fx(j))-davg)**5
      ENDDO
!$omp end parallel do

      ds6 = 0.0D0
!$omp parallel do default(shared) private(j) reduction(+:s6)
      DO j = 1, nin
        ds6 = ds6 + (dble(fx(j))-davg)**6
      ENDDO
!$omp end parallel do

      s(1)=ds2; s(2)=ds3; s(3)=ds4; s(4) = ds5; s(5) = ds6
      CALL MPI_ALLREDUCE(s, gs, 5, MPI_DOUBLE_PRECISION, &
                      MPI_SUM, MPI_COMM_WORLD,ierr)
      if ( ierr.ne.MPI_SUCCESS ) then
        write(*,*)'skewflat: final allreduce failed'
        stop
      endif
      s2=gs(1)*xnorm; s3=gs(2)*xnorm; s4=gs(3)*xnorm; s5=gs(4)*xnorm; s6=gs(5)*xnorm
!     s2=gs(1); s3=gs(2); s4=gs(3); s5=gs(4); s6=gs(5)

      skew = real( s3 / ( s2**1.5 + 1.0e-15 ), kind=GP )
      flat = real( s4 / ( s2**2.0 + 1.0e-15 ), kind=GP )
      glop = real( s5 / ( s2**2.5 + 1.0e-15 ), kind=GP )
      whoa = real( s6 / ( s2**3.0 + 1.0e-15 ), kind=GP )

      RETURN
      END SUBROUTINE skewflat
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!
!
      FUNCTION compute_corr(R1 ,R2) result(gsum)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
! Compute correlation <R1 R2>. Normally the correlation
! would be computed <R1 R1> / ( <R1> <R2> ), 
! but <R1> = <R2> = 0, typically, so we may want
! to normalize by
!    <R1 R2> / ( max(R1) max(R2) )
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
      USE iovar
      USE filefmt
      USE boxsize
 
        IMPLICIT NONE

        REAL   (KIND=GP), INTENT(INOUT), DIMENSION(nx,ny,ksta:kend):: R1,R2

        INTEGER              :: ierr, kz, n
        REAL   (KIND=GP)     :: lsum, gsum, tmp


        CHARACTER(len=1024)  :: fnout
        CHARACTER(len=1024), INTENT(IN)  :: odir
        CHARACTER(len=128)   :: sfld(10), sfldL(10)

        tmp    = 1.0_GP/ ( real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP) )


        lsum = 0.0_GP;
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
        DO k = ksta,kend                                            
!$omp parallel do if (kend-ksta.lt.nth) private (i)               
           DO j = 1,ny  
              DO i = 1,nx                                           
                 lsum = lsum + r1(i,j,k)*r2(i,j,k) 
              END DO
           END DO
        END DO

        CALL MPI_ALLREDUCE(fmin1,gmin,1, GC_REAL,      &
                           MPI_SUM,MPI_COMM_WORLD,ierr)
        gsum = gsum * tmp

      END SUBROUTINE compute_corr
