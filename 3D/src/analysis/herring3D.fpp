!=================================================================
      PROGRAM HERRING3D
!=================================================================
! GHOST code: Geophysical High Order Suite for Turbulence
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
      USE boxsize
      USE gutils
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
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: th
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: fs

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
#ifdef VELOC_
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: M1,M2,M3
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
#ifdef MAGFIELD_
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: C9,C10,C11
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: C12,C13,C14
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: C15,C16,C17
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: M4,M5,M6
#endif
#ifdef HALLTERM_
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: C18
#endif
#ifdef WAVEFUNCTION_
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:)     :: iold,qold,kold,cold
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:)     :: inew,qnew,knew,cnew
#endif
#ifdef EDQNM_
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: C19
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:)     :: tepq,thpq,tve,tvh
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:)     :: Eold,Hold
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:)     :: Eext,Hext
      REAL(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:)    :: Eden,Hden
#endif

      REAL(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:)    :: R1,R2,R3,R4,R5
      REAL(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:)    :: R6,R7,R8,R9,R10
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

      REAL(KIND=GP)    :: bvfreq,xmom,xtemp


      REAL(KIND=GP)    :: kappa
      REAL(KIND=GP)    :: c0,s0
      REAL(KIND=GP)    :: cparam0,cparam1,cparam2,cparam3,cparam4
      REAL(KIND=GP)    :: cparam5,cparam6,cparam7,cparam8,cparam9
      REAL(KIND=GP)    :: sparam0,sparam1,sparam2,sparam3,sparam4
      REAL(KIND=GP)    :: sparam5,sparam6,sparam7,sparam8,sparam9

      REAL(KIND=GP)    :: skup,skdn

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

      REAL(KIND=GP)    :: omegax,omegay,omegaz

#ifdef WAVEFUNCTION_
      REAL(KIND=GP)    :: cspeed,lambda,rho0,kttherm
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
      REAL(KIND=GP)         :: tau, grav
      TYPE (InerGPart)      :: lagpart
#endif
#if defined(TESTPART_) && defined(MAGFIELD_)
      REAL(KIND=GP)         :: gyrof, vtherm
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

      ! Data specific to HERRING3D:
      DOUBLE PRECISION, DIMENSION(3,3) :: bij,dij,gij,vij
      DOUBLE PRECISION                 :: bden,dden,gden,vden
      REAL(kind=GP) sav,ssk,sku,sg5,sw6,ss2,ss3,ss4,ss5,ss6
      REAL(kind=GP) ktmin,ktmax,omega(3),xnormn
      INTEGER :: ic,ir,it,jc
      INTEGER :: bAniso,bHPDF,dolog,inorm,istat(4096),jpdf,nstat
      INTEGER :: nbinx,nbiny,nbins(2)
      INTEGER :: btrunc,useaccum
      LOGICAL :: accum
      CHARACTER(len=64) :: ext1
      CHARACTER(len=4096) :: sstat


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

      NAMELIST / scalar / c0,s0,skdn,skup,kappa,cparam0,cparam1
      NAMELIST / scalar / cparam2,cparam3,cparam4,cparam5,cparam6
      NAMELIST / scalar / cparam7,cparam8,cparam9,sparam0,sparam1
      NAMELIST / scalar / sparam2,sparam3,sparam4,sparam5,sparam6
      NAMELIST / scalar / sparam7,sparam8,sparam9

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
#ifdef WAVEFUNCTION_
      NAMELIST / wavefunction / cspeed,lambda,rho0,kttherm
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
      NAMELIST / pinerpart / tau,grav
#endif
#if defined(TESTPART_) && defined(MAGFIELD_)
      NAMELIST / ptestpart / gyrof,vtherm
#endif
      NAMELIST / shear / iswap,jpdf
      NAMELIST / shear / dolog,useaccum,bAniso,bHPDF,oswap,idir,odir,sstat
      NAMELIST / shear / btrunc,ktmin,ktmax,nbinx,nbiny

!
! Initialization

! Initializes the MPI and I/O libraries
!     CALL MPI_INIT_THREAD(MPI_THREAD_FUNNELED,provided,ierr)
      CALL MPI_INIT(ierr)
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
      ALLOCATE( th (nz,ny,ista:iend) )
#if defined(VELOC_) || defined(ADVECT_)
      ALLOCATE( vx(nz,ny,ista:iend) )
      ALLOCATE( vy(nz,ny,ista:iend) )
      ALLOCATE( vz(nz,ny,ista:iend) )
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
#ifdef DENSITY_
      ALLOCATE( rho(nz,ny,ista:iend) )
      srho = 'rhospect'
#endif
#ifdef SCALAR_
      ALLOCATE( C20(nz,ny,ista:iend) )
      ALLOCATE( fs (nz,ny,ista:iend) )
#endif
#ifdef MULTISCALAR_
      ALLOCATE( C21(nz,ny,ista:iend) )
      ALLOCATE( C22(nz,ny,ista:iend) )
      ALLOCATE( C23(nz,ny,ista:iend) )
      ALLOCATE( C24(nz,ny,ista:iend) )
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
#ifdef WAVEFUNCTION_
      ALLOCATE ( zre(nz,ny,ista:iend), zim(nz,ny,ista:iend) )
#endif
#ifdef QFORCE_
      ALLOCATE ( fre(nz,ny,ista:iend), fim(nz,ny,ista:iend) )
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

      ALLOCATE( R1 (nx,ny,ksta:kend) )
      ALLOCATE( R2 (nx,ny,ksta:kend) )
      ALLOCATE( R3 (nx,ny,ksta:kend) )
      ALLOCATE( R4 (nx,ny,ksta:kend) )
      ALLOCATE( R5 (nx,ny,ksta:kend) )
      ALLOCATE( R6 (nx,ny,ksta:kend) )
      ALLOCATE( R7 (nx,ny,ksta:kend) )
      ALLOCATE( R8 (nx,ny,ksta:kend) )
      ALLOCATE( R9 (nx,ny,ksta:kend) )
      ALLOCATE( R10(nx,ny,ksta:kend) )
#ifdef ADVECT_
      ALLOCATE( vsq(nx,ny,ksta:kend) )
#endif
#ifdef WAVEFUNCTION_
      ALLOCATE( iold(nmax/2+1), qold(nmax/2+1) )
      ALLOCATE( kold(nmax/2+1), cold(nmax/2+1) )
      ALLOCATE( inew(nmax/2+1), qnew(nmax/2+1) )
      ALLOCATE( knew(nmax/2+1), cnew(nmax/2+1) )
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
      ALLOCATE( C19(nz,ny,ista:iend) )
      ALLOCATE( Eden(nz,ny,ista:iend) )
      ALLOCATE( Hden(nz,ny,ista:iend) )
      ALLOCATE( tepq(n/2+1) )      !!!!!!! CHECK LATER !!!!!!!!
      ALLOCATE( thpq(n/2+1) )      ! Here we should have nmax !
      ALLOCATE( tve (n/2+1) )
      ALLOCATE( tvh (n/2+1) )
      ALLOCATE( Eold(n/2+1) )
      ALLOCATE( Hold(n/2+1) )
      ALLOCATE( Eext(3*(n/2+1)) )
      ALLOCATE( Hext(3*(n/2+1)) )
#endif
      th = 0.0
      vx = 0.0
      vy = 0.0
      vz = 0.0
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
!     grav     : Gravitational acceleration

      IF (myrank.eq.0) THEN
         OPEN(1,file='parameter.inp',status='unknown',form="formatted")
         READ(1,NML=pinerpart)
      CLOSE(1)
      ENDIF
      CALL MPI_BCAST(tau    ,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(grav   ,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
#endif

#if defined(TESTPART_) && defined(MAGFIELD_)
! Reads parameters for runs with test particles
!     gyrof    : Gyrofrequency
!     vtherm   : Thermal velocity of the test particles

      IF (myrank.eq.0) THEN
         OPEN(1,file='parameter.inp',status='unknown',form="formatted")
         READ(1,NML=ptestpart)
         CLOSE(1)
      ENDIF
      CALL MPI_BCAST(gyrof    ,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(vtherm   ,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
#endif

! Reads from the external file 'shear`.txt' the 
! parameters that will be used to compute the transfer
!     idir   : directory for unformatted input (field components)
!     odir   : directory for unformatted output (prolongated data)
!     sstat  : time index for which to compute HERRING, or a
!     ';--separated list
!     btrunc : if == 1, truncate spectral range to [ktmin,ktmax]. Not used.
!     ktmin  : min wavenumber for truncation if btrunc=1. Not used.
!     ktmax  : max wavenumber for truncation if btrunc=1. Not used.
!     iswap  : do endian swap on input?
!     oswap  : do endian swap on output? Not used.
!     jpdf   : 0: do no pdfs
!              1: do 1d pdfs only
!              2: do joint pdfs only
!              3: do both 1d and joint pdfs
!     dolog  : compute PDFs in log=space?
!     useaccum: do accumulation over all specified time steps to 
!               compute aniso tensors?
!     bHPDF  : Do PDFs as in Herring m.s.
!     bAniso : Do anisotropy tensor calculations
!
!     Defaults:
      idir   = '.'
      odir   = '.'
      sstat  = '0'
      iswap  = 0
      oswap  = 0
      btrunc = 0
      dolog  = 1
      useaccum = 0
      bHPDF  = 1
      bAniso = 1
      jpdf   = 3
      ktmin  = tiny
      ktmax  = kmax


      IF (myrank.eq.0) THEN
write(*,*)'main: opening herring.inp...'
         OPEN(1,file='herring.inp',status='unknown',form="formatted")
         READ(1,NML=shear)
         CLOSE(1)
write(*,*)'main: herring.inp read.'
      ENDIF
      CALL MPI_BCAST(idir   ,256 ,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(odir   ,256 ,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(sstat  ,4096,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(btrunc ,1   ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(dolog  ,1   ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(useaccum,1   ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(bHPDF  ,1   ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(bAniso ,1   ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(iswap  ,1   ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(jpdf   ,1   ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(ktmin  ,1   ,GC_REAL      ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(ktmax  ,1   ,GC_REAL      ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(nbinx  ,1   ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(nbiny  ,1   ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(oswap  ,1   ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
if (myrank.eq.0) write(*,*)'main: broadcast done.'
! Befor
! options are compatible with the SOLVER being used

      ! Initialize anisotropy tensor data:
      IF ( bAniso .gt. 0 ) THEN
        accum = .FALSE.
        IF ( useaccum .gt. 0 ) THEN
          accum = .TRUE.
        ENDIF
        bij = 0.0; dij = 0.0; gij = 0.0; vij = 0.0;
        bden= 0.0; dden= 0.0; gden= 0.0; vden= 0.0;
if (myrank.eq.0) write(*,*)'main: accum = ', accum
      ENDIF

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
      kmax = (real(n,kind=GP)/2.0_GP-0.5_GP)**2 !!!!!! CHECK !!!!!!!
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

      ! Parse input index set, store in istat:
      CALL parseind(sstat,';', istat , 4096, nstat)
if (myrank.eq.0) write(*,*)'main: index parsing done: nstat=',nstat

      tmp = 1.0_GP/REAL(nx*ny*nz,KIND=GP)
      DO it = 1,nstat
        WRITE(ext1, fmtext) istat(it)
        ext = trim(ext1)
if (myrank.eq.0) write(*,*)'main: Reading time index: ', ext, '...' 
#ifdef MOM_
if (myrank.eq.0) write(*,*)'main: Reading sx...'
        CALL io_read(1,idir,'sx',ext,planio,R1)
if (myrank.eq.0) write(*,*)'main: Reading sy...'
        CALL io_read(1,idir,'sy',ext,planio,R2)
if (myrank.eq.0) write(*,*)'main: Reading sz...'
        CALL io_read(1,idir,'sz',ext,planio,R3)
        CALL fftp3d_real_to_complex(planrc,R1,sx,MPI_COMM_WORLD)
        CALL fftp3d_real_to_complex(planrc,R2,sy,MPI_COMM_WORLD)
        CALL fftp3d_real_to_complex(planrc,R3,sz,MPI_COMM_WORLD)
# ifdef DENSITY_
if (myrank.eq.0) write(*,*)'main: Reading rho...'
        CALL io_read(1,idir,'rho',ext,planio,R1)
if (myrank.eq.0) write(*,*)'main: FFT rho...'
        CALL fftp3d_real_to_complex(planrc,R1,rho,MPI_COMM_WORLD)
if (myrank.eq.0) write(*,*)'main: call mom2vel...'
        CALL mom2vel(rho,sx,sy,sz,0,vx,vy,vz)
#  endif
#else
        CALL io_read(1,idir,'vx',ext,planio,R1)
        CALL io_read(1,idir,'vy',ext,planio,R2)
        CALL io_read(1,idir,'vz',ext,planio,R3)
        CALL fftp3d_real_to_complex(planrc,R1,vx,MPI_COMM_WORLD)
        CALL fftp3d_real_to_complex(planrc,R2,vy,MPI_COMM_WORLD)
        CALL fftp3d_real_to_complex(planrc,R3,vz,MPI_COMM_WORLD)
#endif
#ifdef SCALAR_
        CALL io_read(5,idir,'th',ext,planio,R1)
        CALL fftp3d_real_to_complex(planrc,R1,th,MPI_COMM_WORLD)
#endif
if (myrank.eq.0) write(*,*)'main: Time index ', ext, ' read.' 

!if ( myrank.eq.0 ) write(*,*) 'main: real(vz)=',R1(16,1:16,kend)
if ( myrank.eq.0 ) then
write(*,*) 'main: max(th)=',maxval(R4), ' min(th)=',minval(R4)
endif

if (myrank.eq.0) write(*,*)'main: Do fftr2c... ' 
      xnormn = 1.0_GP/ ( real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP) )
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

        IF ( bHPDF .gt. 0 ) THEN
!if ( myrank.eq.0 ) write(*,*) 'main: vz=',vz(16,1:16,iend)
if (myrank.eq.0) write(*,*)'main: call DoHPDF ...'
          omega(1) = omegax; omega(2) = omegay; omega(3) = omegaz 
          nbins(1) = nbinx ; nbins(2) = nbiny
          CALL DoHPDF(vx,vy,vz,th,bvfreq,omega,nu,kappa,istat(it), &
                      odir,jpdf,nbins,dolog,planio,C1,C2,R1,R2,R3,R4,&
                      R5,R6,R7,R8,R9,R10)
          IF ( myrank.EQ. 0 ) THEN
            write(*,*)'main: time index ', ext, ' done.'
          ENDIF
        ENDIF

        IF ( bAniso .gt. 0 ) THEN
          IF ( useaccum .EQ. 0 ) THEN
            bij = 0.0; dij = 0.0; gij = 0.0; vij = 0.0;
            bden= 0.0; dden= 0.0; gden= 0.0; vden= 0.0;
          ENDIF
          IF ( useaccum .GT. 0 .AND. it .EQ. nstat ) accum = .FALSE.
write(*,*)'main: nstat=', nstat, ' it=', it, ' accum=', accum
          CALL DoAniso(vx,vy,vz,th,istat(it),odir,C1,C2,R1,R2,R3,accum,bden,dden,gden,vden, bij,dij,gij,vij)
        ENDIF

      ENDDO ! end, it-loop

      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
!
! End of HERRING3D

      CALL GTFree(ihcpu1)
      CALL GTFree(ihomp1)
      CALL GTFree(ihwtm1)
      CALL GTFree(ihcpu2)
      CALL GTFree(ihomp2)
      CALL GTFree(ihwtm2)

      CALL MPI_FINALIZE(ierr)
      CALL fftp3d_destroy_plan(plancr)
      CALL fftp3d_destroy_plan(planrc)
      DEALLOCATE( R1,R2,R3,R4,R5,R6,R7,R8,R9,R10 )

      DEALLOCATE( C1,C2 )
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
#ifdef MOM_ 
      DEALLOCATE( sx,sy,sz )
#endif
#ifdef ADVECT_
      DEALLOCATE( vsq )
#endif
#ifdef DENSITY_
      DEALLOCATE( rho )
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
#if defined(INERPART_)
      DEALLOCATE( Rv1,Rv2 )
#endif
#if defined(TESTPART_) && defined(MAGFIELD_)
      DEALLOCATE( Rb1,Rb2,Rb3 )
      DEALLOCATE( Rj1,Rj2,Rj3 )
#endif
      END PROGRAM HERRING3D


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
      SUBROUTINE DoHPDF(vx,vy,vz,th,bvfreq,omega,nu,kappa,indtime, &
                        odir,kin,nbins,dolog,planio,          &
                        ctmp,vtmp,R1,R2,R3,R4,R5,R6,pv,Ri,dissv,dissp)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!
! Computes the 1d and joint PDFs, misc statistics 
! Note: after this call, the input data should be expected to be overwritten.
!
! Parameters
!     vx,
!     vy,
!     vz    : complex velocities
!     th    : pot. temp
!     bvfreq: Brunt-Vaisalla frequency
!     f     : 2x rotation rate
!     omega : rotation freq vector (rank=3)
!     nu    : viscosity
!     kappa : scalar dissipation
!     indtime: integter time index
!     odir  : output directory
!     kin   : 0: compute no pdfs
!             1: compute 1d pdfs only
!             2: compute 2d pdfs only
!             3: compute 1d and 2d pdfs 
!     nbins : 2-elem array with no. bins for (enstrophy,energy diss)
!     dolog : flag to do (1) or not (0) logs of ranges when computing bins
!     planio: io plan, in case we want to output binary data
!     c/vtmp: complex temp array of size vx,vy,vz
!     R1-R6 : real temp arrays
!     pv, Ri, 
!     dissv,
!     dissp : persistent storage for these quantities
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
      USE filefmt
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT   (IN), DIMENSION(nz,ny,ista:iend):: vx,vy,vz,th
      COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(nz,ny,ista:iend):: ctmp,vtmp

      REAL   (KIND=GP), INTENT(INOUT), DIMENSION(nx,ny,ksta:kend):: R1,R2,R3,R4,R5,R6
      REAL   (KIND=GP), INTENT(INOUT), DIMENSION(nx,ny,ksta:kend):: pv,Ri,dissv,dissp
      REAL   (KIND=GP), INTENT   (IN)                            :: bvfreq,omega(3),nu,kappa
      REAL   (KIND=GP)                                           :: fact,fmin(2),fmax(2),xnorm,xnormi,xnormn
      REAL   (KIND=GP)                                           :: omegax,omegay,omegaz
      REAL   (KIND=GP)                                           :: av(100),sk(100),ku(100),g5(100),w6(100)
      REAL   (KIND=GP)                                           :: s2,s3,s4,s5,s6
      REAL   (KIND=GP)                                           :: den  ,rmi
      REAL   (KIND=GP)                                           :: alpha,f
      REAL   (KIND=GP)                                           :: ktmin,ktmax 
      INTEGER         , INTENT   (IN)                            :: dolog,kin,indtime,nbins(2)
      INTEGER                                                    :: i,inorm,j,k,knz,n,nin,sr
      INTEGER                                                    :: btrunc
      LOGICAL                                                    :: bexist
      CHARACTER(len=*), INTENT   (IN)                            :: odir
      CHARACTER(len=1024)                                        :: fnout
      CHARACTER(len=16)                                          :: sfld(100)
      CHARACTER(len=128)                                         :: hdrfmt, rowfmt
      TYPE(IOPLAN)    , INTENT(INOUT)                            :: planio

      rmi = 1000. ! max |Richardson| to allow
      inorm = 0;

      btrunc = 0;
      ktmin  = 0.0
      ktmax  = 0.0

    
      WRITE(ext, fmtext) indtime

      omegax = omega(1)
      omegay = omega(2)
      omegaz = omega(3)
      f      = 2 * omegaz

if ( myrank.eq.0 ) write(*,*) 'DoHPDF: vz=',vz(16,1:16,iend)
      knz    = kend - ksta + 1

if ( myrank.eq.0 ) write(*,*) 'DoHPDF: entering... ext=', trim(ext)
!!    read((ext,*) indtime ! convert time index string to integer

if ( myrank.eq.0 ) write(*,*) 'DoHPDF: normalize input data...'
      nin = nx*ny*(kend-ksta+1)

      av = 0.0; sk = 0.0; ku = 0.0; g5 = 0.0; w6 = 0.0;

      n = 0

      ! At this point forward, complex velocity components are normalized!

      ! Compute 'persistent' fields...
      ! ... First, dissipation:
if ( myrank.eq.0 ) write(*,*) 'DoHPDF: call compute_dissv...'
      CALL compute_dissv(vx,vy,vz,btrunc,ktmin,ktmax,inorm,ctmp,vtmp, &
                  R1,R2,R3,R4,R5,R6,dissv)
#ifdef SCALAR_
if ( myrank.eq.0 ) write(*,*) 'DoHPDF: call compute_dissp...'
      CALL compute_dissp(th,ctmp,R1,R2,R3,dissp)
#endif
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
        DO j = 1,ny
          DO i = 1,nx
!           dissv(i,j,k) =  0.01*real(2*j+k+i)
            dissv(i,j,k) =  dissv(i,j,k) * 2.0*nu
            dissp(i,j,k) =  dissp(i,j,k) *     kappa
          ENDDO
        ENDDO
      ENDDO

      ! Compute for flatness/skewness:
#ifdef SCALAR_
if ( myrank.eq.0 ) write(*,*) 'DoHPDF: call compute_PV...'
      ! ... PV :
      CALL compute_PV(vx,vy,vz,th,bvfreq,f,inorm,ctmp,R1,R2,R6,pv) 

      ! ... gradient Richardson :
if ( myrank.eq.0 ) write(*,*) 'DoHPDF: call compute_Rig...'
      CALL compute_Rig(vx,vy,vz,th,bvfreq,f,1,ctmp,R1,R2,R3,Ri)
#endif

      n = n + 1; sfld(n) = 'dissv' 
if ( myrank.eq.0 ) write(*,*) 'DoHPDF: call stat ', sfld(n)
      CALL skewflat(dissv,nx,ny,knz,av(n),sk(n),ku(n),g5(n),w6(n),s2,s3,s4,s5,s6)
      If ( ibits(kin,0,1).EQ.1 ) THEN
        fnout = trim(odir) // '/' // 'dissvpdf.' // ext // '.txt'
        CALL dopdfr(dissv ,nx,ny,knz,fnout,nbins(1),0,fmin(1),fmax(1),0) 
      ENDIF

#ifdef SCALAR_
      n = n + 1; sfld(n) = 'dissp' 
if ( myrank.eq.0 ) write(*,*) 'DoHPDF: call stat ', sfld(n)
      CALL skewflat(dissp,nx,ny,knz,av(n),sk(n),ku(n),g5(n),w6(n),s2,s3,s4,s5,s6)
      If ( ibits(kin,0,1).EQ.1 ) THEN
        fnout = trim(odir) // '/' // 'dissppdf.' // ext // '.txt'
        CALL dopdfr(dissp ,nx,ny,knz,fnout,nbins(1),0,fmin(1),fmax(1),0) 
      ENDIF

      n = n + 1; sfld(n) = 'PV' 
if ( myrank.eq.0 ) write(*,*) 'DoHPDF: call stat ', sfld(n)
      CALL skewflat(pv,nx,ny,knz,av(n),sk(n),ku(n),g5(n),w6(n),s2,s3,s4,s5,s6)
      If ( ibits(kin,0,1).EQ.1 ) THEN
        fnout = trim(odir) // '/' // 'pvpdf.' // ext // '.txt'
        CALL dopdfr(pv,nx,ny,knz,fnout,nbins(1),0,fmin(1),fmax(1),0) 
      ENDIF

      If ( ibits(kin,1,1).EQ.1 ) THEN ! joint pdfs for pv:
        fnout = trim(odir) // '/' // 'jpdf_pv_dissv_log00.' // ext // '.txt'
        CALL dojpdfr(pv,'PV',dissv,'dissv',nx,ny,knz,fnout,nbins,[0,0],fmin,fmax,[0,0])
        fnout = trim(odir) // '/' // 'jpdf_pv_Ri_log00.' // ext // '.txt'
        CALL dojpdfr(pv,'PV',dissv,'Ri',nx,ny,knz,fnout,nbins,[0,0],fmin,fmax,[0,0])
        fnout = trim(odir) // '/' // 'jpdf_Ri_dissv_log00.' // ext // '.txt'
        CALL dojpdfr(Ri,'Ri',pv,'PV',nx,ny,knz,fnout,nbins,[0,0],fmin,fmax,[0,0])
      ENDIF

#endif

      ! Compute statistics for omega_z:
      CALL rotor3(vy,vz,ctmp,3)
      CALL fftp3d_complex_to_real(plancr,ctmp,R1,MPI_COMM_WORLD)
      n = n + 1; sfld(n) = 'z-curl v' 
if ( myrank.eq.0 ) write(*,*) 'DoHPDF: call stat ', sfld(n)
      CALL skewflat(R1,nx,ny,knz,av(n),sk(n),ku(n),g5(n),w6(n),s2,s3,s4,s5,s6)
      If ( ibits(kin,0,1).EQ.1 ) THEN
        fnout = trim(odir) // '/' // 'omzpdf.' // ext // '.txt'
        CALL dopdfr(R1,nx,ny,knz,fnout,nbins(1),0,fmin(1),fmax(1),0) 
      ENDIF
      If ( ibits(kin,1,1).EQ.1 ) THEN ! joint pdfs for omegaz:
        fnout = trim(odir) // '/' // 'jpdf_omz_dissv.' // ext // '.txt'
        CALL dojpdfr(R1,'omz',dissv,'dissv',nx,ny,knz,fnout,nbins,[0,0],fmin,fmax,[0,0])
#ifdef SCALAR_
        fnout = trim(odir) // '/' // 'jpdf_omz_Ri.' // ext // '.txt'
        CALL dojpdfr(R1,'omz',Ri,'Ri',nx,ny,knz,fnout,nbins,[0,0],fmin,fmax,[0,0])
        fnout = trim(odir) // '/' // 'jpdf_omz_PV.' // ext // '.txt'
        CALL dojpdfr(R1,'omz',pv,'PV',nx,ny,knz,fnout,nbins,[0,0],fmin,fmax,[0,0])
#endif
      ENDIF

      ! Compute statistics for v_z:
      ctmp = vz;
      CALL fftp3d_complex_to_real(plancr,ctmp,R1,MPI_COMM_WORLD)
      n = n + 1; sfld(n) = 'v_z' 
if ( myrank.eq.0 ) write(*,*) 'DoHPDF: call stat ', sfld(n)
      CALL skewflat(R1,nx,ny,knz,av(n),sk(n),ku(n),g5(n),w6(n),s2,s3,s4,s5,s6)
      If ( ibits(kin,0,1).EQ.1 ) THEN
        fnout = trim(odir) // '/' // 'vzpdf.' // ext // '.txt'
        CALL dopdfr(R1,nx,ny,knz,fnout,nbins(1),0,fmin(1),fmax(1),0) 
      ENDIF
      If ( ibits(kin,1,1).EQ.1 ) THEN ! joint pdfs for pv:
        fnout = trim(odir) // '/' // 'jpdf_vz_dissv.' // ext // '.txt'
        CALL dojpdfr(R1,'vz',dissv,'dissv',nx,ny,knz,fnout,nbins,[0,0],fmin,fmax,[0,0])
#ifdef SCALAR_
        fnout = trim(odir) // '/' // 'jpdf_vz_Ri.' // ext // '.txt'
        CALL dojpdfr(R1,'vz',Ri,'Ri',nx,ny,knz,fnout,nbins,[0,0],fmin,fmax,[0,0])
        fnout = trim(odir) // '/' // 'jpdf_vz_PV.' // ext // '.txt'
        CALL dojpdfr(R1,'vz',pv,'PV',nx,ny,knz,fnout,nbins,[0,0],fmin,fmax,[0,0])
#endif
      ENDIF

      ! Compute statistics for v_\perp:
      CALL compute_vperp(vx,vy,vz, ctmp,R2,R1)
      n = n + 1; sfld(n) = 'v_perp' 
if ( myrank.eq.0 ) write(*,*) 'DoHPDF: call stat ', sfld(n)
      CALL skewflat(R1,nx,ny,knz,av(n),sk(n),ku(n),g5(n),w6(n),s2,s3,s4,s5,s6)
      If ( ibits(kin,0,1).EQ.1 ) THEN
        fnout = trim(odir) // '/' // 'vperppdf.' // ext // '.txt'
        CALL dopdfr(R1,nx,ny,knz,fnout,nbins(1),0,fmin(1),fmax(1),0) 
      ENDIF
      If ( ibits(kin,1,1).EQ.1 ) THEN ! joint pdfs for pv:
        fnout = trim(odir) // '/' // 'jpdf_vperp_dissv.' // ext // '.txt'
        CALL dojpdfr(R1,'vz',dissv,'dissv',nx,ny,knz,fnout,nbins,[0,0],fmin,fmax,[0,0])
#ifdef SCALAR_
        fnout = trim(odir) // '/' // 'jpdf_vperp_Ri.' // ext // '.txt'
        CALL dojpdfr(R1,'vz',Ri,'Ri',nx,ny,knz,fnout,nbins,[0,0],fmin,fmax,[0,0])
        fnout = trim(odir) // '/' // 'jpdf_vperp_PV.' // ext // '.txt'
        CALL dojpdfr(R1,'vz',pv,'PV',nx,ny,knz,fnout,nbins,[0,0],fmin,fmax,[0,0])
#endif
      ENDIF

#ifdef SCALAR_
      ! Compute statistics for theta:
      ctmp = th;
      CALL fftp3d_complex_to_real(plancr,ctmp,R2,MPI_COMM_WORLD)
      n = n + 1; sfld(n) = 'th' 
if ( myrank.eq.0 ) write(*,*) 'DoHPDF: call stat ', sfld(n)
      CALL skewflat(R2,nx,ny,knz,av(n),sk(n),ku(n),g5(n),w6(n),s2,s3,s4,s5,s6)
      If ( ibits(kin,0,1).EQ.1 ) THEN
        fnout = trim(odir) // '/' // 'thpdf.' // ext // '.txt'
        CALL dopdfr(R2,nx,ny,knz,fnout,nbins(1),0,fmin(1),fmax(1),0) 
      ENDIF
      If ( ibits(kin,1,1).EQ.1 ) THEN ! joint pdfs for theta:
        fnout = trim(odir) // '/' // 'jpdf_th_Ri.' // ext // '.txt'
        CALL dojpdfr(R2,'vz',Ri,'Ri',nx,ny,knz,fnout,nbins,[0,0],fmin,fmax,[0,0])
        fnout = trim(odir) // '/' // 'jpdf_th_dissv.' // ext // '.txt'
        CALL dojpdfr(R2,'vz',dissv,'dissv',nx,ny,knz,fnout,nbins,[0,0],fmin,fmax,[0,0])
        fnout = trim(odir) // '/' // 'jpdf_th_PV.' // ext // '.txt'
        CALL dojpdfr(R2,'vz',pv,'PV',nx,ny,knz,fnout,nbins,[0,0],fmin,fmax,[0,0])
      ENDIF
#endif

      CALL compute_vort(vx,vy,vz,ctmp,R1,R2,R3,R4)
      n = n + 1; sfld(n) = 'vort' 
      CALL skewflat(R4,nx,ny,knz,av(n),sk(n),ku(n),g5(n),w6(n),s2,s3,s4,s5,s6)
      If ( ibits(kin,0,1).EQ.1 ) THEN
        fnout = trim(odir) // '/' // 'vortpdf.' // ext // '.txt'
        CALL dopdfr(R4,nx,ny,knz,fnout,nbins(1),0,fmin(1),fmax(1),0) 
      ENDIF

#ifdef SCALAR_
      ! Compute statistics for buoyancy flux: N theta v_z:
      ctmp = th;
      CALL fftp3d_complex_to_real(plancr,ctmp,R1,MPI_COMM_WORLD)
      ctmp = vz;
      CALL fftp3d_complex_to_real(plancr,ctmp,R2,MPI_COMM_WORLD)
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
        DO j = 1,ny
          DO i = 1,nx
            R2(i,j,k) =  bvfreq * R1(i,j,k) * R2(i,j,k)
          ENDDO
        ENDDO
      ENDDO
      n = n + 1; sfld(n) = 'Bf' 
if ( myrank.eq.0 ) write(*,*) 'DoHPDF: call stat ', sfld(n)
      CALL skewflat(R2,nx,ny,knz,av(n),sk(n),ku(n),g5(n),w6(n),s2,s3,s4,s5,s6)
      If ( ibits(kin,0,1).EQ.1 ) THEN
        fnout = trim(odir) // '/' // 'Bfpdf.' // ext // '.txt'
        CALL dopdfr(R2,nx,ny,knz,fnout,nbins(1),0,fmin(1),fmax(1),0) 
      ENDIF
      If ( ibits(kin,1,1).EQ.1 ) THEN ! joint pdfs for theta:
        fnout = trim(odir) // '/' // 'jpdf_Bf_Ri.' // ext // '.txt'
        CALL dojpdfr(R2,'Bf',Ri,'Ri',nx,ny,knz,fnout,nbins,[0,0],fmin,fmax,[0,0])
        fnout = trim(odir) // '/' // 'jpdf_Bf_dissv.' // ext // '.txt'
        CALL dojpdfr(R2,'Bf',dissv,'dissv',nx,ny,knz,fnout,nbins,[0,0],fmin,fmax,[0,0])
        fnout = trim(odir) // '/' // 'jpdf_Bf_PV.' // ext // '.txt'
        CALL dojpdfr(R2,'Bf',pv,'PV',nx,ny,knz,fnout,nbins,[0,0],fmin,fmax,[0,0])
      ENDIF
#endif

#ifdef SCALAR_
      ! Compute statistics for mixing ratio, Gamf = Bf/dissv:
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
        DO j = 1,ny
          DO i = 1,nx
            R3(i,j,k) =  R2(i,j,k) / dissv(i,j,k)
          ENDDO
        ENDDO
      ENDDO
      n = n + 1; sfld(n) = 'Gamf' 
if ( myrank.eq.0 ) write(*,*) 'DoHPDF: call stat ', sfld(n)
      CALL skewflat(R3,nx,ny,knz,av(n),sk(n),ku(n),g5(n),w6(n),s2,s3,s4,s5,s6)
      If ( ibits(kin,0,1).EQ.1 ) THEN
        fnout = trim(odir) // '/' // 'Gamfpdf.' // ext // '.txt'
        CALL dopdfr(R3,nx,ny,knz,fnout,nbins(1),0,fmin(1),fmax(1),0) 
      ENDIF
      If ( ibits(kin,1,1).EQ.1 ) THEN ! joint pdfs for theta:
        fnout = trim(odir) // '/' // 'jpdf_Gamf_Ri.' // ext // '.txt'
        CALL dojpdfr(R3,'Gamf',Ri,'Ri',nx,ny,knz,fnout,nbins,[0,0],fmin,fmax,[0,0])
        fnout = trim(odir) // '/' // 'jpdf_Gamf_dissv.' // ext // '.txt'
        CALL dojpdfr(R3,'Gamf',dissv,'dissv',nx,ny,knz,fnout,nbins,[0,0],fmin,fmax,[0,0])
        fnout = trim(odir) // '/' // 'jpdf_Gamf_PV.' // ext // '.txt'
        CALL dojpdfr(R3,'Gamf',pv,'PV',nx,ny,knz,fnout,nbins,[0,0],fmin,fmax,[0,0])
      ENDIF

      ! Compute statistics for mixing ratio, Gamr = dissp/dissv:
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
        DO j = 1,ny
          DO i = 1,nx
            R4(i,j,k) =  dissp(i,j,k) / dissv(i,j,k)
          ENDDO
        ENDDO
      ENDDO
      n = n + 1; sfld(n) = 'Gamr' 
if ( myrank.eq.0 ) write(*,*) 'DoHPDF: call stat ', sfld(n)
      CALL skewflat(R4,nx,ny,knz,av(n),sk(n),ku(n),g5(n),w6(n),s2,s3,s4,s5,s6)
      If ( ibits(kin,0,1).EQ.1 ) THEN
        fnout = trim(odir) // '/' // 'Gamrpdf.' // ext // '.txt'
        CALL dopdfr(R4,nx,ny,knz,fnout,nbins(1),0,fmin(1),fmax(1),0) 
      ENDIF
      If ( ibits(kin,1,1).EQ.1 ) THEN ! joint pdfs for theta:
        fnout = trim(odir) // '/' // 'jpdf_Gamr_Ri.' // ext // '.txt'
        CALL dojpdfr(R4,'Gamr',Ri,'Ri',nx,ny,knz,fnout,nbins,[0,0],fmin,fmax,[0,0])
        fnout = trim(odir) // '/' // 'jpdf_Gamr_dissv.' // ext // '.txt'
        CALL dojpdfr(R4,'Gamr',dissv,'dissv',nx,ny,knz,fnout,nbins,[0,0],fmin,fmax,[0,0])
        fnout = trim(odir) // '/' // 'jpdf_Gamr_PV.' // ext // '.txt'
        CALL dojpdfr(R4,'Gamr',pv,'PV',nx,ny,knz,fnout,nbins,[0,0],fmin,fmax,[0,0])
      ENDIF
#endif

#ifdef SCALAR_
      ! Compute statistics for d theta/dz:
      CALL derivk3(th, ctmp, 3)
      CALL fftp3d_complex_to_real(plancr,ctmp,R1,MPI_COMM_WORLD)
      n = n + 1; sfld(n) = 'dthdz' 
if ( myrank.eq.0 ) write(*,*) 'DoHPDF: call stat ', sfld(n)
      CALL skewflat(R1,nx,ny,knz,av(n),sk(n),ku(n),g5(n),w6(n),s2,s3,s4,s5,s6)
      If ( ibits(kin,0,1).EQ.1 ) THEN
        fnout = trim(odir) // '/' // 'dthdzpdf.' // ext // '.txt'
        CALL dopdfr(R1,nx,ny,knz,fnout,nbins(1),0,fmin(1),fmax(1),0) 
      ENDIF
      If ( ibits(kin,1,1).EQ.1 ) THEN ! joint pdfs for theta:
        fnout = trim(odir) // '/' // 'jpdf_dthdz_Ri.' // ext // '.txt'
        CALL dojpdfr(R1,'vz',Ri,'Ri',nx,ny,knz,fnout,nbins,[0,0],fmin,fmax,[0,0])
        fnout = trim(odir) // '/' // 'jpdf_dthdz_dissv.' // ext // '.txt'
        CALL dojpdfr(R1,'vz',dissv,'dissv',nx,ny,knz,fnout,nbins,[0,0],fmin,fmax,[0,0])
        fnout = trim(odir) // '/' // 'jpdf_dthdz_PV.' // ext // '.txt'
        CALL dojpdfr(R1,'vz',pv,'PV',nx,ny,knz,fnout,nbins,[0,0],fmin,fmax,[0,0])
      ENDIF
#endif

      ! Compute statistics for gradient Ri:
      n = n + 1; sfld(n) = 'Rig' 
if ( myrank.eq.0 ) write(*,*) 'DoHPDF: call stat ', sfld(n)
      fmin(1) = -10.0; fmax(1) = 10.0;
      CALL skewflatb(Ri,nx,ny,knz,1,fmin(1),fmax(1),av(n),sk(n),ku(n),g5(n),w6(n),s2,s3,s4,s5,s6)
if ( abs(av(n)) .gt. fmax(1) ) then
  write(*,*) 'DoHPDF: BAD ', sfld(n), ' avg=', av(n)
  stop
endif
      If ( ibits(kin,0,1).EQ.1 ) THEN
        fnout = trim(odir) // '/' // 'Rigpdf.' // ext // '.txt'
        CALL dopdfr(Ri,nx,ny,knz,fnout,nbins(1),1,fmin(1),fmax(1),0) 
      ENDIF

      ! Compute statistics for Ri = N tau_shear: 
      ! NOTE: Ri overwritten!
      CALL compute_Rig(vx,vy,vz,th,bvfreq,f,2,ctmp,R1,R2,R3,Ri)

      n = n + 1; sfld(n) = 'Ri' 
if ( myrank.eq.0 ) write(*,*) 'DoHPDF: call stat ', sfld(n)
      fmin(1) = -10.0; fmax(1) = 10.0;
      CALL skewflatb(Ri,nx,ny,knz,1,fmin(1),fmax(1),av(n),sk(n),ku(n),g5(n),w6(n),s2,s3,s4,s5,s6)
if ( abs(av(n)) .gt. fmax(1) ) then
  write(*,*) 'DoHPDF: BAD ', sfld(n), ' avg=', av(n)
  stop
endif
      If ( ibits(kin,0,1).EQ.1 ) THEN
        fnout = trim(odir) // '/' // 'Ripdf.' // ext // '.txt'
        CALL dopdfr(Ri,nx,ny,knz,fnout,nbins(1),1,fmin(1),fmax(1),0) 
      ENDIF

      ! Create format for statistical data:
      WRITE(rowfmt,'(A, I4, A)') '(I4,',n,'(2X,E14.6))'
      WRITE(hdrfmt,'(A, I4, A)') '(A,',n,'(2X,A))'

!if ( myrank.eq.0 ) write(*,*) 'DoHPDF: do file write: rowfmt', rowfmt, ' hdrfmt=',hdrfmt
      ! Print out high order statistics  data:
      IF ( myrank.EQ.0 ) THEN
        inquire( file='kavg.txt', exist=bexist )
        
        fnout = trim(odir) // '/' // 'kavg.txt'
        OPEN(2,file=trim(fnout),position='append')
        if ( .NOT. bexist ) THEN
        WRITE(2,hdrfmt,advance='yes') '#itime', (sfld(j), j=1,n)
        ENDIF
        WRITE(2,rowfmt,advance='no') indtime, (av(j), j=1,n)
        CLOSE(2)

        fnout = trim(odir) // '/' // 'skew.txt'
        OPEN(2,file=trim(fnout),position='append')
        if ( .NOT. bexist ) THEN
        WRITE(2,hdrfmt,advance='yes') '#itime', (sfld(j), j=1,n)
        ENDIF
        WRITE(2,rowfmt,advance='no') indtime, (sk(j), j=1,n)
        CLOSE(2)

        fnout = trim(odir) // '/' // 'flat.txt'
        OPEN(2,file=trim(fnout),position='append')
        if ( .NOT. bexist ) THEN
        WRITE(2,hdrfmt,advance='yes') '#itime', (sfld(j), j=1,n)
        ENDIF
        WRITE(2,rowfmt,advance='no') indtime, (ku(j), j=1,n)
        CLOSE(2)

        fnout = trim(odir) // '/' // 'glop.txt'
        OPEN(2,file=trim(fnout),position='append')
        if ( .NOT. bexist ) THEN
        WRITE(2,hdrfmt,advance='yes') '#itime', (sfld(j), j=1,n)
        ENDIF
        WRITE(2,rowfmt,advance='no') indtime, (g5(j), j=1,n)
        CLOSE(2)

        fnout = trim(odir) // '/' // 'whoa.txt'
        OPEN(2,file=trim(fnout),position='append')
        if ( .NOT. bexist ) THEN
        WRITE(2,hdrfmt,advance='yes') '#itime', (sfld(j), j=1,n)
        ENDIF
        WRITE(2,rowfmt,advance='no') indtime, (w6(j), j=1,n)
        CLOSE(2)
      ENDIF

#if 0
write(*,*) 'DoHPDF: hdrfmt= ', trim(hdrfmt)
if ( indtime .eq. 233 ) then
  if ( myrank.eq.0 ) then
    write(*,*) '                          n=', n
    do j=1,n
      write(*,*) '       .................j=',j, ' name=', sfld(j), ' av=',av(j), ' sk=', sk(j), ' ku=', ku(j)
    enddo
  endif
  stop
endif
#endif
      
if ( myrank.eq.0 ) write(*,*) 'DoHPDF: done.'

      END SUBROUTINE DoHPDF
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

      END SUBROUTINE skewflat
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!
!
      SUBROUTINE skewflatb(fx,nx,ny,nz,ifixdr,fmin,fmax,avg,skew,flat,glop,whoa,s2,s3,s4,s5,s6)
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
      USE gutils
      IMPLICIT NONE

      REAL(KIND=GP), INTENT(INOUT), DIMENSION(*)   :: fx
      REAL(KIND=GP), INTENT   (IN)                 :: fmin,fmax
      REAL(KIND=GP), INTENT  (OUT)                 :: skew,flat,glop,whoa
      REAL(KIND=GP), INTENT  (OUT)                 :: avg,s2,s3,s4,s5,s6
      DOUBLE PRECISION                             :: davg,ds2,ds3,ds4,ds5,ds6
      DOUBLE PRECISION                             :: dtot,ltot, gs(5),s(5),xnorm
      DOUBLE PRECISION                             :: dkeep
      INTEGER      , INTENT   (IN)                 :: ifixdr,nx,ny,nz
      INTEGER                                      :: gnz,ierr,j

      INTEGER  nin

      nin = nx * ny * nz
      IF ( .NOT. ALLOCATED(ikeep_) .OR. nin.GT.nikeep_ ) THEN
        IF ( ALLOCATED(ikeep_ ) ) DEALLOCATE(ikeep_)
        ALLOCATE(ikeep_(nin))
      ENDIF
      nikeep_ = nin


      IF ( ifixdr .EQ. 0 ) THEN
!$omp parallel do default(shared) private(j) reduction(+:s2)
        DO j = 1, nin
          ikeep_(j) = j
        ENDDO
      ELSE
!$omp parallel do default(shared) private(j) reduction(+:s2)
        nikeep_ = 0
        DO j = 1, nin
          IF ( fx(j) .GE. fmin .AND. fx(j) .LE. fmax ) THEN
            nikeep_ = nikeep_ + 1
            ikeep_(nikeep_) = j
          ENDIF
        ENDDO
      ENDIF

      ltot = dble(nikeep_)
      CALL MPI_ALLREDUCE(ltot, dtot, 1, MPI_DOUBLE_PRECISION, &
                      MPI_SUM, MPI_COMM_WORLD,ierr)
      
      xnorm = 1.0D0 / dtot
      ds2 = 0.0D0
!$omp parallel do default(shared) private(j) reduction(+:s2)
      DO j = 1, nikeep_
        ds2 = ds2 + dble(fx(ikeep_(j)))
      ENDDO
!$omp end parallel do

      CALL MPI_ALLREDUCE(ds2, davg, 1, MPI_DOUBLE_PRECISION, &
                      MPI_SUM, MPI_COMM_WORLD,ierr)
dkeep = davg
      davg = davg * xnorm
      avg  = davg

if ( davg .gt. fmax .or. davg .lt. fmin ) then
       write(*,*)'skewflatb: fmin=',fmin, ' fmax=', fmax, ' nkeep=', nikeep_, ' nin=', nin, ' xnorm=', xnorm, ' davg=', davg, ' dkeep=', dkeep
endif
!      write(*,*)'skewflatb: nkeep=', nikeep_, ' nin=', nin,' davg=', davg

      ds2 = 0.0D0
!$omp parallel do default(shared) private(j) reduction(+:s2)
      DO j = 1, nikeep_
        ds2 = ds2 + (dble(fx(ikeep_(j)))-davg)**2
      ENDDO
!$omp end parallel do

      ds3 = 0.0D0
!$omp parallel do default(shared) private(j) reduction(+:s3)
      DO j = 1, nikeep_
        ds3 = ds3 + (dble(fx(ikeep_(j)))-davg)**3
      ENDDO
!$omp end parallel do

      ds4 = 0.0D0
!$omp parallel do default(shared) private(j) reduction(+:s4)
      DO j = 1, nikeep_
        ds4 = ds4 + (dble(fx(ikeep_(j)))-davg)**4
      ENDDO
!$omp end parallel do

      ds5 = 0.0D0
!$omp parallel do default(shared) private(j) reduction(+:s5)
      DO j = 1, nikeep_
        ds5 = ds5 + (dble(fx(ikeep_(j)))-davg)**5
      ENDDO
!$omp end parallel do

      ds6 = 0.0D0
!$omp parallel do default(shared) private(j) reduction(+:s6)
      DO j = 1, nikeep_
        ds6 = ds6 + (dble(fx(ikeep_(j)))-davg)**6
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

      END SUBROUTINE skewflatb
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


      SUBROUTINE compute_dissv(vx,vy,vz,btrunc,ktmin,ktmax,inorm,ctmp1,ctmp2, &
                                S11,S12,S13,S22,S23,S33,diss)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!
! Computes the real energy dissipation field
!
! Parameters
!     vi   : input velocities
!     ktmin: truncaton min wavenumber for spherical truncation
!     ktmax: truncaton max wavenumber for spherical truncation
!     inorm: normalize (1), or not (0)
!     ctmp*: complex temp array
!     SIJ  : real temp arrays
!     diss : real dissipation rate field, returned
!
      USE fprecision
      USE commtypes
      USE kes
      USE grid
      USE mpivars
      USE ali
      USE fft
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT   (IN), DIMENSION(nz,ny,ista:iend) :: vx,vy,vz
      COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(nz,ny,ista:iend) :: ctmp1,ctmp2
      REAL   (KIND=GP), INTENT(INOUT), DIMENSION(nx,ny,ksta:kend) :: S11,S12,S13,S22,S23,S33
      REAL   (KIND=GP), INTENT(INOUT), DIMENSION(nx,ny,ksta:kend) :: diss
      REAL   (KIND=GP), INTENT   (IN)                             :: ktmin,ktmax
      REAL   (KIND=GP)                                            :: ktmin2,ktmax2,tmp
!
      INTEGER         , INTENT   (IN)                             :: btrunc,inorm
      INTEGER                                                     :: i,j,k

#if 1
S11 = 0.; S12 = 0.; S13=0.; S22 = 0.; S23 = 0.; S33 = 0.
      CALL Strain(vx,vy,vz,1,1,btrunc,ktmin,ktmax,inorm,ctmp1,ctmp2)
      CALL fftp3d_complex_to_real(plancr,ctmp2,diss,MPI_COMM_WORLD)
      CALL Strain(vx,vy,vz,1,2,btrunc,ktmin,ktmax,inorm,ctmp1,ctmp2)
      CALL fftp3d_complex_to_real(plancr,ctmp2,S12,MPI_COMM_WORLD)
      CALL Strain(vx,vy,vz,1,3,btrunc,ktmin,ktmax,inorm,ctmp1,ctmp2)
      CALL fftp3d_complex_to_real(plancr,ctmp2,S13,MPI_COMM_WORLD)
      CALL Strain(vx,vy,vz,2,2,btrunc,ktmin,ktmax,inorm,ctmp1,ctmp2)
      CALL fftp3d_complex_to_real(plancr,ctmp2,S22,MPI_COMM_WORLD)

      CALL Strain(vx,vy,vz,2,3,btrunc,ktmin,ktmax,inorm,ctmp1,ctmp2)
      CALL fftp3d_complex_to_real(plancr,ctmp2,S23,MPI_COMM_WORLD)
      CALL Strain(vx,vy,vz,3,3,btrunc,ktmin,ktmax,inorm,ctmp1,ctmp2)
      CALL fftp3d_complex_to_real(plancr,ctmp2,S33,MPI_COMM_WORLD)
#endif

#if 0
      ! Compute S33 for flatness/skewness:
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
        DO j = 1,ny
          DO i = 1,nx
            S33(i,j,k) =  -(S11(i,j,k)+S22(i,j,k)) 
          ENDDO
        ENDDO
      ENDDO
#endif


      ! Compute normalized energy dissipation field/2nu, store in ! S11:
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
        DO j = 1,ny
          DO i = 1,nx
#if 1
            diss(i,j,k) =          S11(i,j,k)*S11(i,j,k) &
                            +      S22(i,j,k)*S22(i,j,k) &
                            +      S33(i,j,k)*S33(i,j,k) &
                            +  2.0*S12(i,j,k)*S12(i,j,k) &
                            +  2.0*S13(i,j,k)*S13(i,j,k) &
                            +  2.0*S23(i,j,k)*S23(i,j,k)  
#else
            diss(i,j,k) =   0.01*(2*j+k+i)
#endif
          ENDDO
        ENDDO
      ENDDO

      END SUBROUTINE compute_dissv
!-----------------------------------------------------------------
!-----------------------------------------------------------------

      SUBROUTINE compute_dissp(th,ctmp1,R1,R2,R3,diss)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!
! Computes the real potential energy dissipation field
!
! Parameters
!     th   : input potential temp
!     ctmp*: complex temp array
!     R*   : real temp arrays
!     diss : real dissipation rate field, returned
!
      USE fprecision
      USE commtypes
      USE kes
      USE grid
      USE mpivars
      USE ali
      USE fft
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT   (IN), DIMENSION(nz,ny,ista:iend) :: th
      COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(nz,ny,ista:iend) :: ctmp1
      REAL   (KIND=GP), INTENT(INOUT), DIMENSION(nx,ny,ksta:kend) :: R1,R2,R3
      REAL   (KIND=GP), INTENT(INOUT), DIMENSION(nx,ny,ksta:kend) :: diss
      INTEGER                                                     :: i,j,k


      ! Compute gradient:
      CALL derivk3(th, ctmp1, 1)  ! x-deriv
      CALL fftp3d_complex_to_real(plancr,ctmp1,R1,MPI_COMM_WORLD)
      CALL derivk3(th, ctmp1, 2)  ! y-deriv
      CALL fftp3d_complex_to_real(plancr,ctmp1,R2,MPI_COMM_WORLD)
      CALL derivk3(th, ctmp1, 3)  ! z-deriv
      CALL fftp3d_complex_to_real(plancr,ctmp1,R3,MPI_COMM_WORLD)


      ! Compute dissipation field/2nu:
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
        DO j = 1,ny
          DO i = 1,nx
            diss(i,j,k) =          R1(i,j,k)*R1(i,j,k) &
                            +      R2(i,j,k)*R2(i,j,k) &
                            +      R3(i,j,k)*R3(i,j,k) 
          ENDDO
        ENDDO
      ENDDO

      END SUBROUTINE compute_dissp
!-----------------------------------------------------------------
!-----------------------------------------------------------------

      SUBROUTINE compute_vort(vx,vy,vz,ctmp1,R1,R2,R3,vort)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!
! Computes the real |curl v|^2
!
! Parameters
!     v*   : input velocity components
!     ctmp*: complex temp array
!     R*   : real temp arrays
!     vort : real vorticity field
!
      USE fprecision
      USE commtypes
      USE kes
      USE grid
      USE mpivars
      USE ali
      USE fft
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT   (IN), DIMENSION(nz,ny,ista:iend) :: vx,vy,vz
      COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(nz,ny,ista:iend) :: ctmp1
      REAL   (KIND=GP), INTENT(INOUT), DIMENSION(nx,ny,ksta:kend) :: R1,R2,R3
      REAL   (KIND=GP), INTENT(INOUT), DIMENSION(nx,ny,ksta:kend) :: vort
      INTEGER                                                     :: i,j,k


      ! Compute gradient:
      CALL rotor3(vy,vz,ctmp1,1)
      CALL fftp3d_complex_to_real(plancr,ctmp1,R1,MPI_COMM_WORLD)

      CALL rotor3(vz,vx,ctmp1,2)
      CALL fftp3d_complex_to_real(plancr,ctmp1,R2,MPI_COMM_WORLD)

      CALL rotor3(vx,vy,ctmp1,3)
      CALL fftp3d_complex_to_real(plancr,ctmp1,R3,MPI_COMM_WORLD)

      ! Compute |curl v|^2
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
        DO j = 1,ny
          DO i = 1,nx
            vort (i,j,k) =         R1(i,j,k)*R1(i,j,k) &
                            +      R2(i,j,k)*R2(i,j,k) &
                            +      R3(i,j,k)*R3(i,j,k) 
          ENDDO
        ENDDO
      ENDDO

      END SUBROUTINE compute_vort
!-----------------------------------------------------------------
!-----------------------------------------------------------------

      SUBROUTINE compute_PV(vx,vy,vz,th,bvfreq,f,inorm,ctmp1,R1,R2,omgth,pv)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!
! Computes the real PV field:
            !pv = f dtheta/dz -N \omega_z +  curl v . Grad theta 
!
! Parameters
!     vx,vy,
!     vz,th : compled input field
!     bvfreq: Brunt-Vaisalla frequency, N
!     f     : 2x rotation rate in z
!     inorm : normalize (1), or not (0)
!     ctmp* : complex temp array
!     R1-R5 : real temp arrays
!     omgth : curl v . Grad theta, returned
!     pv    : real PV field, returned
!
      USE fprecision
      USE commtypes
      USE kes
      USE grid
      USE mpivars
      USE ali
      USE fft
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT   (IN), DIMENSION(nz,ny,ista:iend) :: vx,vy,vz,th
      COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(nz,ny,ista:iend) :: ctmp1
      REAL   (KIND=GP), INTENT   (IN)                             :: bvfreq,f
      REAL   (KIND=GP), INTENT(INOUT), DIMENSION(nx,ny,ksta:kend) :: R1,R2
      REAL   (KIND=GP), INTENT(INOUT), DIMENSION(nx,ny,ksta:kend) :: omgth,pv
!
      INTEGER         , INTENT   (IN)                             :: inorm
      INTEGER                                                     :: i,j,k

      !! Compute stats for PV quantities:
      !! ...First, omega .Grad theta:
!!    xnormn = 1.0_GP/ real(nx*ny*nz,kind=GP)**2
      CALL derivk3(th, ctmp1, 1)  ! x-deriv
      CALL fftp3d_complex_to_real(plancr,ctmp1,R1,MPI_COMM_WORLD)
      CALL rotor3(vy,vz,ctmp1,1)
      CALL fftp3d_complex_to_real(plancr,ctmp1,R2,MPI_COMM_WORLD)
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
        DO j = 1,ny
          DO i = 1,nx
            omgth(i,j,k) = R1(i,j,k)*R2(i,j,k) ! += omega_x Grad_x theta
          ENDDO
        ENDDO
      ENDDO

      CALL derivk3(th, ctmp1, 2)  ! y-deriv
      CALL fftp3d_complex_to_real(plancr,ctmp1,R1,MPI_COMM_WORLD)
      CALL rotor3(vz,vx,ctmp1,2)
      CALL fftp3d_complex_to_real(plancr,ctmp1,R2,MPI_COMM_WORLD)
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
        DO j = 1,ny
          DO i = 1,nx
            omgth(i,j,k) = omgth(i,j,k) + R1(i,j,k)*R2(i,j,k) ! += omega_y Grad_y theta
          ENDDO
        ENDDO
      ENDDO

      CALL derivk3(th, ctmp1, 3)  ! z-deriv
      CALL fftp3d_complex_to_real(plancr,ctmp1,R1,MPI_COMM_WORLD) ! dtheta/dz
      CALL rotor3(vx,vy,ctmp1,3)
      CALL fftp3d_complex_to_real(plancr,ctmp1,R2,MPI_COMM_WORLD) ! omega_z
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
        DO j = 1,ny
          DO i = 1,nx
            omgth(i,j,k) = omgth(i,j,k) + R1(i,j,k)*R2(i,j,k) ! += omega_z Grad_z theta
            !pv = f dtheta/dz -N \omega_z +  curl v . Grad theta:
            pv  (i,j,k) = f*R1(i,j,k) - bvfreq*R2(i,j,k) + omgth(i,j,k) 
          ENDDO
        ENDDO
      ENDDO

      END SUBROUTINE compute_PV
!-----------------------------------------------------------------
!-----------------------------------------------------------------

      SUBROUTINE compute_Rig(vx,vy,vz,th,bvfreq,f,itype,ctmp1,R1,R2,R3,ri)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!
! Computes the real gradient Richardsond number field:
!   type      Form
!    1      Rig = N (N - d\theta/dz) / ( (dvx/dz)^2 + (dvy/dz)^2 ) (itype 1)
!
!    2      Ri  = N / sqrt( (dvx/dz)^2 + (dvy/dz)^2 ) (itype 2) == N tau_shear
!
!    3      Rig = N*a (N*a - d\theta/dz) / ( (dvx/dz)^2 + (dvy/dz)^2 ) (itype 3)
!          
!           where
! 
!           a = N / sqrt(f^2 + N^2)
!
! Parameters
!     vx,vy,
!     vz,th: compled input field
!     bvfreq: Brunt-Vaisalla frequency
!     f    : 2x rotation rate in z
!     itype: 1, 2, 3
!     ctmp*: complex temp array
!     R1-R3: real temp arrays
!     ri   : real Rig field, returned
!
      USE fprecision
      USE commtypes
      USE kes
      USE grid
      USE mpivars
      USE ali
      USE fft
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT   (IN), DIMENSION(nz,ny,ista:iend) :: vx,vy,vz,th
      COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(nz,ny,ista:iend) :: ctmp1
      REAL   (KIND=GP), INTENT   (IN)                             :: bvfreq,f
      REAL   (KIND=GP), INTENT(INOUT), DIMENSION(nx,ny,ksta:kend) :: R1,R2,R3
      REAL   (KIND=GP), INTENT(INOUT), DIMENSION(nx,ny,ksta:kend) :: ri
      REAL   (KIND=GP)                                            :: alpha, den, xnormn
!
      INTEGER         , INTENT   (IN)                             :: itype
      INTEGER                                                     :: i,j,k

      ! Compute gradient Richardson no. and its pdf:
      CALL derivk3(vx, ctmp1, 3)
      CALL fftp3d_complex_to_real(plancr,ctmp1,R1,MPI_COMM_WORLD)
      CALL derivk3(vy, ctmp1, 3)
      CALL fftp3d_complex_to_real(plancr,ctmp1,R2,MPI_COMM_WORLD)
      CALL derivk3(th, ctmp1, 3)
      CALL fftp3d_complex_to_real(plancr,ctmp1,R3,MPI_COMM_WORLD)

      alpha  = bvfreq/sqrt(f**2 + bvfreq**2)
      xnormn = 1.0_GP/ ( real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP) )

      IF ( itype .EQ. 1 ) THEN
        ! NOTE: no normalization required as vx, vy, vz, and theta were
        ! normalized on read-in
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
        DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
          DO j = 1,ny
            DO i = 1,nx
              den = R1(i,j,k)**2 + R2(i,j,k)**2
              den = 1.0 / (den + 1.0e-15)
              ri(i,j,k) = bvfreq*(bvfreq-R3(i,j,k)) * den
            ENDDO
          ENDDO
        ENDDO
      ELSE IF ( itype .EQ. 2 ) THEN
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
        DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
          DO j = 1,ny
            DO i = 1,nx
              den = sqrt( R1(i,j,k)**2 + R2(i,j,k)**2 ) 
              den = 1.0 / (den + 1.0e-15)
              ri(i,j,k) = bvfreq * den
            ENDDO
          ENDDO
        ENDDO
      ELSE IF ( itype .EQ. 3 ) THEN
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
        DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
          DO j = 1,ny
            DO i = 1,nx
              den = R1(i,j,k)**2 + R2(i,j,k)**2 
              den = 1.0 / (den + 1.0e-15)
              ri(i,j,k) = bvfreq*alpha*(bvfreq*alpha-R3(i,j,k)) * den
            ENDDO
          ENDDO
        ENDDO
      ELSE
        stop 'Subroutine Rig: invalid itype'
      ENDIF

      END SUBROUTINE compute_Rig
!-----------------------------------------------------------------
!-----------------------------------------------------------------

      SUBROUTINE compute_vperp(vx,vy,vz, ctmp1,R1,vp)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!
! Computes real v_perp
!
! Parameters
!     vx,vy,
!     vz     input field
!     ctmp*: complex temp array
!     R1-R3: real temp arrays
!     vp   : real vperp field, returned
!
      USE fprecision
      USE commtypes
      USE kes
      USE grid
      USE mpivars
      USE ali
      USE fft
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT   (IN), DIMENSION(nz,ny,ista:iend) :: vx,vy,vz
      COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(nz,ny,ista:iend) :: ctmp1
      REAL   (KIND=GP), INTENT(INOUT), DIMENSION(nx,ny,ksta:kend) :: R1
      REAL   (KIND=GP), INTENT(INOUT), DIMENSION(nx,ny,ksta:kend) :: vp
!
      INTEGER                                                     :: i,j,k

      ! Compute gradient Richardson no. and its pdf:
      ctmp1 = vx
      CALL fftp3d_complex_to_real(plancr,ctmp1,vp,MPI_COMM_WORLD)
      ctmp1 = vy
      CALL fftp3d_complex_to_real(plancr,ctmp1,R1,MPI_COMM_WORLD)

!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
        DO j = 1,ny
          DO i = 1,nx
            vp(i,j,k) = sqrt(vp(i,j,k)*vp(i,j,k) + R1(i,j,k)*R1(i,j,k))
          ENDDO
        ENDDO
      ENDDO

      END SUBROUTINE compute_vperp
!-----------------------------------------------------------------
!-----------------------------------------------------------------

      SUBROUTINE DoAniso(vx,vy,vz,th,indtime,odir,C1,C2,R1,R2,R3,accum,bdenom,ddenom,gdenom,vdenom, bij,dij,gij,vij)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!
! Computes the 1d and joint PDFs, misc statistics 
! Note: after this call, the input data should be expected to be overwritten.
!
! Parameters
!     vx,
!     vy,
!     vz     : complex velocities
!     th     : pot. temp
!     indtime: integter time index
!     odir   : output directory
!     accum  : if TRUE, continues to accumulate the aniso tensors and normalizations.
!              If FALSE, final accumulation is done, and global sums are done to
!              compute tensors
!     bij,gij,
!     vij,dij: aniso tensors, returned. First time in, should be initialized to 0
!     bdenom,
!     ...   ,
!     ddenom: Tensor normalizations, returned. First time in, should be initialized to 0
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
      USE filefmt
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT   (IN), DIMENSION(nz,ny,ista:iend):: vx,vy,vz,th
      COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(nz,ny,ista:iend):: C1,C2

      REAL   (KIND=GP), INTENT(INOUT), DIMENSION(nx,ny,ksta:kend):: R1,R2,R3
      DOUBLE PRECISION,                DIMENSION(4,3)            :: invar
      DOUBLE PRECISION, INTENT(INOUT), DIMENSION(3,3), TARGET    :: bij,vij,gij,dij
      DOUBLE PRECISION, INTENT(INOUT)                            :: bdenom,vdenom,gdenom,ddenom
      TYPE(PMAT)                                                 :: pm(4)
      LOGICAL                                                    :: bexist
      LOGICAL         , INTENT   (IN)                            :: accum
      INTEGER         , INTENT   (IN)                            :: indtime
      INTEGER                                                    :: i,j
      CHARACTER(len=*), INTENT   (IN)                            :: odir
      CHARACTER(len=1024)                                        :: fnout
      CHARACTER(len=128)                                         :: rowfmt

      WRITE(rowfmt,'(A, I4, A)') '(I4,',12,'(2X,E14.6))'

      pm(1).mat => bij
      pm(2).mat => dij
      pm(3).mat => gij
      pm(4).mat => vij
      CALL anisobij(vx,vy,vz,C1,R1,R2,R3,accum,bdenom,bij)
      CALL anisodij(vx,vy,vz,C1,C2,R1,R2,accum,ddenom,dij)
      CALL anisogij(th,C1,R1,R2,R3,accum,gdenom,gij)
      CALL anisovij(vx,vy,vz,C1,C2,R1,R2,R3,accum,vdenom,vij)

      IF (.not. accum) THEN
        DO i = 1, 3 ! which invariant, I, II, III
          DO j = 1, 4 ! which tensor
            CALL invariant(pm(j).mat, i, invar(j,i))
          ENDDO
        ENDDO
      ENDIF


      IF (myrank.eq.0 .AND. .not. accum) THEN
      fnout = trim(odir) // '/' // 'invar.txt'
      inquire( file=fnout, exist=bexist )
      OPEN(2,file=trim(fnout),position='append')
      if ( .NOT. bexist ) THEN
      WRITE(2,'(A, 4x, 12(A, 3x))') '#itime', 'bI', 'bII', 'bIII', 'dI', 'dII', 'dIII', 'gI', 'gII', 'gIII', 'vI', 'vII', 'vIII'
      ENDIF
      WRITE(2,rowfmt,advance='no') &
                          indtime, invar(1,1), invar(1,2), invar(1,3), &
                          invar(2,1), invar(2,2), invar(2,3), &
                          invar(3,1), invar(3,2), invar(3,3), &
                          invar(4,1), invar(4,2), invar(4,3) 
      CLOSE(2)
      ENDIF


      END SUBROUTINE DoAniso
