!=================================================================
! GHOST GPSplineInt interpolation class
! The spline used (that on the 'moments' or 2nd derivatives)
! as well as the tridiagonal solve for the moments, are
! derived from several sources including the presentation 
! given in the book 'The theory of splines and their application'
! 1967, by Ahlberg, Nilson, Walsh.
!
! There are 3 main public calls to this class:
!
!   CALL this%intop_%PartUpdate3D(this%px_,this%py_,this%pz_,this%nparts_)
!   CALL this%intop_%CompSpline3D(evar,tmp1,tmp2)
!   CALL this%intop_%DoInterp3D(lag,nl)
! 
! The first updates the interpolation points; the second computes the
! control point function, which, here, is global; the third ccarries out
! the actual interpolation to the interpolation points. 
!
! Note that this class assumes a regularly spaced grid.
!
! 2013 A. Pumir (ENS, Lyon)
!      D. Rosenberg (NCCS: ORNL)
!
! 15 Jan 2013: Initial version
! 19 Jul 2016: CompSpline3D loops optimized for speed (P. Mininni)
! 10 May 2023: Added support for GPU offloading (F. Pugliese)
!=================================================================
MODULE class_GPSplineInt
      USE mpivars
      USE fprecision
      USE class_GPartComm
      USE gtimer
!$    USE offloading
!$    USE, INTRINSIC :: iso_c_binding
      IMPLICIT NONE
  
      PRIVATE
      TYPE, PUBLIC :: GPSplineInt
        PRIVATE
        ! Member data:
        REAL(KIND=GP),POINTER,DIMENSION  (:,:)   :: esplfld2_
        REAL(KIND=GP),POINTER,DIMENSION(:,:,:)   :: esplfld_
        REAL(KIND=GP),POINTER,DIMENSION    (:)   :: ax_,bx_,betx_,cx_,gamx_,px_,xxx_
        REAL(KIND=GP),POINTER,DIMENSION    (:)   :: ay_,by_,bety_,cy_,gamy_,py_,xxy_
        REAL(KIND=GP),POINTER,DIMENSION    (:)   :: az_,bz_,betz_,cz_,gamz_,pz_,xxz_
        REAL(KIND=GP),POINTER,DIMENSION    (:)   :: xrk_,yrk_,zrk_
        REAL(KIND=GP),POINTER,DIMENSION  (:,:)   :: wrkl_
        REAL(KIND=GP)                            :: dxi_(3),xbnds_(3,2),zetax_,zetay_,zetaz_
        TYPE(GPartComm),POINTER                  :: gpcomm_
        INTEGER      ,POINTER,DIMENSION  (:,:)   :: ilg_,jlg_,klg_
        INTEGER                                  :: maxint_
        INTEGER                                  :: ierr_,ider_(3),nd_(3)
        INTEGER                                  :: ibnds_(3,2),obnds_(3,2)
        INTEGER                                  :: ldims_(3),odims_(3)
        INTEGER                                  :: hdataex_,htransp_
        INTEGER                                  :: ntot_,ttot_
        INTEGER                                  :: rank_
        CHARACTER(len=1024)                      :: serr_
#if defined(DO_HYBRIDoffl)
!$      TYPE(c_ptr)                              :: dptrax,dptrpx,dptrbetx,dptrgamx,dptrxxx
!$      TYPE(c_ptr)                              :: dptray,dptrpy,dptrbety,dptrgamy,dptrxxy
!$      TYPE(c_ptr)                              :: dptraz,dptrpz,dptrbetz,dptrgamz,dptrxxz
!$      TYPE(c_ptr)                              :: dptrxrk,dptrilg,dptrwrkl
!$      TYPE(c_ptr)                              :: dptryrk,dptrjlg,dptresplfld
!$      TYPE(c_ptr)                              :: dptrzrk,dptrklg,dptresplfld2
#endif
      CONTAINS
        ! Public methods:
        PROCEDURE,PUBLIC :: GPSplineInt_ctor
        FINAL            :: GPSplineInt_dtor

        ! ...Aliases:
        PROCEDURE,PUBLIC :: Init         => GPSplineInt_Init
        PROCEDURE,PUBLIC :: DoInterp2D   => GPSplineInt_Interp2D
        PROCEDURE,PUBLIC :: DoInterp3D   => GPSplineInt_Interp3D
        PROCEDURE,PUBLIC :: SetDeriv     => GPSplineInt_SetDeriv
        PROCEDURE,PUBLIC :: PartUpdate2D => GPSplineInt_PartUpdate2D
        PROCEDURE,PUBLIC :: PartUpdate3D => GPSplineInt_PartUpdate3D
        PROCEDURE,PUBLIC :: CompSpline2D => GPSplineInt_CompSpline2D
        PROCEDURE,PUBLIC :: CompSpline3D => GPSplineInt_CompSpline3D
        PROCEDURE,PUBLIC :: ResizeArrays => GPSplineInt_ResizeArrays

!!      GENERIC  ,PUBLIC :: PartUpdate   => PartUpdate2D,PartUpdate3D
      END TYPE GPSplineInt

      ! Private methods:
      PRIVATE :: GPSplineInt_Init          , GPSplineInt_MatInvQ
      PRIVATE :: GPSplineInt_Interp2D      , GPSplineInt_Interp3D
      PRIVATE :: GPSplineInt_DoAlloc       , GPSplineInt_DoDealloc
      PRIVATE :: GPSplineInt_CompSpline2D  , GPSplineInt_CompSpline3D
      PRIVATE :: GPSplineInt_PartUpdate2D  , GPSplineInt_PartUpdate3D

! Methods:
  CONTAINS

  SUBROUTINE GPSplineInt_ctor(this,rank,nd,ibnds,xbnds,obnds,nzghost,maxpart,gpcomm, &
                              hdataex,htransp)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  Main explicit constructor
!  ARGUMENTS:
!    this    : 'this' class instance
!    rank    : rank of solution (1, 2, or 3)
!    nd      : global grid size in each direction
!    ibnds   : starting and ending indices for each direction for this MPI task.
!              Integer array of size (rank,2).
!    xbnds   : task's slab bounds
!    obnds   : starting and ending indices for each transposed array for this task
!              Integer array of size (rank,2).
!    nzghost : no. z-ghost zones required for interpolation
!    maxpart : max no. interpolation points/Lag. particles
!    gpcomm  : GHOST particle communicator object
!    hdataex : handle to data exchange. Must be valid.
!    htransp : handle to timer for transpose. Must be valid.
!-----------------------------------------------------------------
    IMPLICIT NONE
    CLASS(GPSplineInt)                           :: this
    TYPE(GPartComm),TARGET                       :: gpcomm
    INTEGER        ,INTENT(IN)                   :: hdataex,htransp,maxpart,nzghost,rank
    INTEGER        ,INTENT(IN),DIMENSION  (rank) :: nd
    INTEGER        ,INTENT(IN),DIMENSION(rank,2) :: ibnds,obnds
    INTEGER                                      :: j,k
    REAL(KIND=GP)  ,INTENT(IN),DIMENSION(rank,2) :: xbnds

    this%gpcomm_ => gpcomm
    this%maxint_  = maxpart
    this%rank_    = rank
    this%ider_    = 0
    this%ldims_   = 0
    this%odims_   = 0
    this%ntot_    = 1
    this%ttot_    = 1
    j = GTValidHandle(htransp)
    IF ( j.NE.GTERR_GOOD_HANDLE ) THEN
      WRITE(*,*) 'GPSplineInt_ctor: invalid transpose timer handle: ',j
      STOP
    ENDIF
    this%htransp_  = htransp
    j = GTValidHandle(hdataex)
    IF ( j.NE.GTERR_GOOD_HANDLE ) THEN
      WRITE(*,*) 'GPSplineInt_ctor: invalid data exch. timer handle: ',j
      STOP
    ENDIF
    this%hdataex_  = hdataex

    DO j = 1, this%rank_
      DO k = 1,2
        this%ibnds_(j,k)  = ibnds(j,k)
        this%obnds_(j,k)  = obnds(j,k)
        this%xbnds_(j,k)  = real(ibnds(j,k),kind=GP)-1.0_GP
      ENDDO
      this%ldims_(j)  = ibnds(j,2) - ibnds(j,1) + 1
      this%odims_(j)  = obnds(j,2) - obnds(j,1) + 1
      this%nd_   (j)  = nd   (j)
      this%ntot_ = this%ntot_*this%ldims_(j)
      this%ttot_ = this%ttot_*this%odims_(j)
    ENDDO
    ! Note: Expand by no. ghost zones:
    this%xbnds_(3,1)  = this%xbnds_(3,1)-real(nzghost,kind=GP)
    this%xbnds_(3,2)  = this%xbnds_(3,2)+real(nzghost,kind=GP)

    IF ( this%rank_.LT.2 .OR. this%rank_.GT.3 ) THEN
      WRITE(*,*)'GPSplineInt::ctor: Invalid rank'
      STOP
    ENDIF

    IF ( this%rank_.EQ.2 .AND. this%ldims_(2).NE.(this%ibnds_(2,2)-this%ibnds_(2,1)+1) ) THEN
      WRITE(*,*) 'GPSplineInt::ctor: Inconsistent 2-indices: rank=2'
      STOP
    ENDIF

    IF ( this%rank_.EQ.3 .AND. this%ldims_(3).NE.(this%ibnds_(3,2)-this%ibnds_(3,1)+1) ) THEN
      WRITE(*,*) 'GPSplineInt::ctor: Inconsitent 3-indices; rank=3'
      STOP
    ENDIF

    CALL GPSplineInt_Init(this)

  END SUBROUTINE GPSplineInt_ctor
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  
  SUBROUTINE GPSplineInt_dtor(this)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  Main explicit destructor
!  ARGUMENTS:
!    this    : 'this' class instance
!-----------------------------------------------------------------

    IMPLICIT NONE
    TYPE(GPSplineInt),INTENT(INOUT)        :: this

    CALL GPSplineInt_DoDealloc(this)

  END SUBROUTINE GPSplineInt_dtor
!-----------------------------------------------------------------
!-----------------------------------------------------------------


  SUBROUTINE GPSplineInt_Interp2D(this,fp,np)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  DESCRIPTION:  Interp2D
!               Carry out 2d interpolation. Spline control points
!               must have been computed on field to be interpolated
!               by call to GPSplineInt_CompSpline, before entering.
!  ARGUMENTS  : 
!    this     : 'this' class instance
!    fp       : field interpolated at the points xp, yp, returned.
!    xp,yp    : x,y locations of interpolation points
!    np       : no. interploation points (lenghts of fp,xp,yp).
!-----------------------------------------------------------------
    IMPLICIT NONE
    CLASS(GPSplineInt)                        :: this
    INTEGER      ,INTENT   (IN)               :: np
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(np) :: fp

!

  END SUBROUTINE GPSplineInt_Interp2D
!-----------------------------------------------------------------
!-----------------------------------------------------------------


  SUBROUTINE GPSplineInt_Interp3D(this,fp,np)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  DESCRIPTION:  Interp3D
!               Carry out 3d interpolation. Spline control points
!               must have been computed on field to be interpolated
!               by call to GPSplineInt_CompSpline, before entering.
!               In addition, the normalized positions of the particles
!               in each interval must be computed prior to entry,
!               as must be the indices into the full grid that
!               define the range of control points. 
!  ARGUMENTS  : 
!    this     : 'this' class instance
!    fp       : field interpolated at the points xp, yp, returned.
!    np       : no. interploation points (lenghts of fp,xp,yp,zp).
!-----------------------------------------------------------------
    IMPLICIT NONE
    CLASS(GPSplineInt)                       :: this
    INTEGER      ,INTENT   (IN)              :: np
    INTEGER                                  :: lag,nx,ny,nxy,nz
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(np):: fp
    REAL(KIND=GP)                            :: four,half,halfm,sixth,six,three,threeh,two
    REAL(KIND=GP)                            :: xx,xsm,xxm,yy,ysm,yym,zz,zsm,zzm
    REAL(KIND=GP)                            :: xx1,xx2,xx3,xx4,yy1,yy2,yy3,yy4,zz1,zz2,zz3,zz4
!
    nx = this%ldims_(1)
    ny = this%ldims_(2)
    nz = this%ldims_(3)
    nxy = nx*ny
!
! Do the interpolation stricto sensu :
    sixth  = 1.0_GP/6.0_GP
    four   = 4.0_GP
    three  = 3.0_GP
    six    = 6.0_GP
    half   = 0.5_GP
    halfm  = -half
    threeh = 3.0_GP/2.0_GP
    two    = 2.0_GP
!
#if defined(DO_HYBRIDoffl)
!$omp target update to(this%ilg_,this%jlg_,this%klg_,this%esplfld_,this%wrkl_,      &
!$omp                  this%xrk_,this%yrk_,this%zrk_) device(targetdev)
!$omp target data map(to:np,nx,ny,nz,nxy,sixth,four,three,six,half,halfm,threeh,    &
!$omp                    two,this%dxi_) map(alloc:fp(1:np)) device(targetdev)
#endif
    IF ( this%ider_(1).EQ.0 ) THEN
    xsm=1.0_GP
#if defined(DO_HYBRIDoffl)
!$omp target teams distribute parallel do private(xx,xxm) map(this%wrkl_,this%xrk_)&
!$omp                                         device(targetdev)
#else
!$omp parallel do private(xx,xxm)
#endif
    DO lag=1,np
      xx = this%xrk_(lag)
      xxm = (1.0_GP-xx)
      this%wrkl_(1,lag) = sixth*xxm*xxm*xxm
      this%wrkl_(2,lag) = sixth*(four+xx*xx*(three*xx-six))
      this%wrkl_(3,lag) = sixth*(four+xxm*xxm*(three*xxm-six))
    ENDDO
    ENDIF
!
    IF ( this%ider_(1).EQ.1 ) THEN
    xsm=0.0_GP
#if defined(DO_HYBRIDoffl)
!$omp target teams distribute parallel do private(xx,xxm) map(this%wrkl_,this%xrk_)&
!$omp                                                     device(targetdev)
#else
!$omp parallel do private(xx,xxm)
#endif
    DO lag=1,np
      xx = this%xrk_(lag)
      xxm = (1.0_GP-xx)
      this%wrkl_(1,lag) = halfm*xxm*xxm*this%dxi_(1)
      this%wrkl_(2,lag) = xx*(threeh*xx-two)*this%dxi_(1)
      this%wrkl_(3,lag) = -xxm*(threeh*xxm-two)*this%dxi_(1)
    ENDDO
    ENDIF
!
    IF ( this%ider_(2).EQ.0 ) THEN
    ysm=1.0_GP
#if defined(DO_HYBRIDoffl)
!$omp target teams distribute parallel do private(yy,yym) map(this%wrkl_,this%yrk_)&
!$omp                                                     device(targetdev)
#else
!$omp parallel do private(yy,yym)
#endif
    DO lag=1,np
      yy = this%yrk_(lag)
      yym = (1.0_GP-yy)
      this%wrkl_(4,lag) = sixth*yym*yym*yym
      this%wrkl_(5,lag) = sixth*(four+yy*yy*(three*yy-six))
      this%wrkl_(6,lag) = sixth*(four+yym*yym*(three*yym-six))
    ENDDO
    ENDIF
!
    IF( this%ider_(2).EQ.1 ) THEN
    ysm=0.0_GP
#if defined(DO_HYBRIDoffl)
!$omp target teams distribute parallel do private(yy,yym) map(this%wrkl_,this%yrk_)&
!$omp                                                     device(targetdev)
#else
!$omp parallel do private(yy,yym)
#endif
    DO lag=1,np
      yy = this%yrk_(lag)
      yym = (1.0_GP-yy)
      this%wrkl_(4,lag) = halfm*yym*yym*this%dxi_(2)
      this%wrkl_(5,lag) = yy*(threeh*yy-two)*this%dxi_(2)
      this%wrkl_(6,lag) = -yym*(threeh*yym-two)*this%dxi_(2)
    ENDDO
    ENDIF
!
    IF ( this%ider_(3).EQ.0 ) THEN
    zsm=1.0_GP
#if defined(DO_HYBRIDoffl)
!$omp target teams distribute parallel do private(zz,zzm) map(this%wrkl_,this%zrk_)&
!$omp                                                     device(targetdev)
#else
!$omp parallel do private(zz,zzm)
#endif
    DO lag=1,np
      zz = this%zrk_(lag)
      zzm = (1.0_GP-zz)
      this%wrkl_(7,lag) = sixth*zzm*zzm*zzm
      this%wrkl_(8,lag) = sixth*(four+zz*zz*(three*zz-six))
      this%wrkl_(9,lag) = sixth*(four+zzm*zzm*(three*zzm-six))
    ENDDO
    ENDIF
!
    IF( this%ider_(3).EQ.1 ) THEN
    zsm=0.0
#if defined(DO_HYBRIDoffl)
!$omp target teams distribute parallel do private(zz,zzm) map(this%wrkl_,this%zrk_)&
!$omp                                                     device(targetdev)
#else
!$omp parallel do private(zz,zzm)
#endif
    DO lag=1,np
      zz = this%zrk_(lag)
      zzm = (1.0_GP-zz)
      this%wrkl_(7,lag) = halfm*zzm*zzm*this%dxi_(3)
      this%wrkl_(8,lag) = zz*(threeh*zz-two)*this%dxi_(3)
      this%wrkl_(9,lag) = -zzm*(threeh*zzm-two)*this%dxi_(3)
    ENDDO
    ENDIF
!
    IF(this%ider_(1).eq.0) xsm = 1.0_GP
    IF(this%ider_(1).eq.1) xsm = 0.0_GP
    IF(this%ider_(2).eq.0) ysm = 1.0_GP
    IF(this%ider_(2).eq.1) ysm = 0.0_GP
    IF(this%ider_(3).eq.0) zsm = 1.0_GP
    IF(this%ider_(3).eq.1) zsm = 0.0_GP
!
#if defined(DO_HYBRIDoffl)
!$omp target teams distribute parallel do private(xx1,xx2,xx3,xx4,yy1,yy2,yy3,yy4,    &
!$omp                                             zz1,zz2,zz3,zz4) map(to:xsm,ysm,zsm)&
!$omp     map(this%wrkl_,this%ilg_,this%jlg_,this%klg_,this%esplfld_) device(targetdev)
#else
!$omp parallel do private(xx1,xx2,xx3,xx4,yy1,yy2,yy3,yy4,zz1,zz2,zz3,zz4)
#endif
    DO lag=1,np
      xx1 = this%wrkl_(1,lag)
      xx2 = this%wrkl_(2,lag)
      xx3 = this%wrkl_(3,lag)
      xx4 = xsm - xx1 - xx2 - xx3
      yy1 = this%wrkl_(4,lag)
      yy2 = this%wrkl_(5,lag)
      yy3 = this%wrkl_(6,lag)
      yy4 = ysm - yy1 - yy2 - yy3
      zz1 = this%wrkl_(7,lag)
      zz2 = this%wrkl_(8,lag)
      zz3 = this%wrkl_(9,lag)
      zz4 = zsm - zz1 - zz2 - zz3
!write(*,*)'part=',lag,' xx1=',xx1,' xx2=',xx2,' xx3=',xx3,' xx4=',xx4
!write(*,*)'part=',lag,' yy1=',yy1,' yy2=',yy2,' yy3=',yy3,' ty4=',yy4
!write(*,*)'part=',lag,' zz1=',zz1,' zz2=',zz2,' zz3=',zz3,' zz4=',zz4

      fp(lag) =  &
        this%esplfld_(this%ilg_(1,lag),this%jlg_(1,lag),this%klg_(1,lag))*xx1*yy1*zz1 &
      + this%esplfld_(this%ilg_(2,lag),this%jlg_(1,lag),this%klg_(1,lag))*xx2*yy1*zz1 &
      + this%esplfld_(this%ilg_(3,lag),this%jlg_(1,lag),this%klg_(1,lag))*xx3*yy1*zz1 &
      + this%esplfld_(this%ilg_(4,lag),this%jlg_(1,lag),this%klg_(1,lag))*xx4*yy1*zz1 &
      + this%esplfld_(this%ilg_(1,lag),this%jlg_(2,lag),this%klg_(1,lag))*xx1*yy2*zz1 &
      + this%esplfld_(this%ilg_(2,lag),this%jlg_(2,lag),this%klg_(1,lag))*xx2*yy2*zz1 &
      + this%esplfld_(this%ilg_(3,lag),this%jlg_(2,lag),this%klg_(1,lag))*xx3*yy2*zz1 &
      + this%esplfld_(this%ilg_(4,lag),this%jlg_(2,lag),this%klg_(1,lag))*xx4*yy2*zz1 &
      + this%esplfld_(this%ilg_(1,lag),this%jlg_(3,lag),this%klg_(1,lag))*xx1*yy3*zz1 &
      + this%esplfld_(this%ilg_(2,lag),this%jlg_(3,lag),this%klg_(1,lag))*xx2*yy3*zz1 &
      + this%esplfld_(this%ilg_(3,lag),this%jlg_(3,lag),this%klg_(1,lag))*xx3*yy3*zz1 &
      + this%esplfld_(this%ilg_(4,lag),this%jlg_(3,lag),this%klg_(1,lag))*xx4*yy3*zz1 &
      + this%esplfld_(this%ilg_(1,lag),this%jlg_(4,lag),this%klg_(1,lag))*xx1*yy4*zz1 &
      + this%esplfld_(this%ilg_(2,lag),this%jlg_(4,lag),this%klg_(1,lag))*xx2*yy4*zz1 &
      + this%esplfld_(this%ilg_(3,lag),this%jlg_(4,lag),this%klg_(1,lag))*xx3*yy4*zz1 &
      + this%esplfld_(this%ilg_(4,lag),this%jlg_(4,lag),this%klg_(1,lag))*xx4*yy4*zz1 & 
!  k = 2
      + this%esplfld_(this%ilg_(1,lag),this%jlg_(1,lag),this%klg_(2,lag))*xx1*yy1*zz2 &
      + this%esplfld_(this%ilg_(2,lag),this%jlg_(1,lag),this%klg_(2,lag))*xx2*yy1*zz2 &
      + this%esplfld_(this%ilg_(3,lag),this%jlg_(1,lag),this%klg_(2,lag))*xx3*yy1*zz2 &
      + this%esplfld_(this%ilg_(4,lag),this%jlg_(1,lag),this%klg_(2,lag))*xx4*yy1*zz2 &
      + this%esplfld_(this%ilg_(1,lag),this%jlg_(2,lag),this%klg_(2,lag))*xx1*yy2*zz2 &
      + this%esplfld_(this%ilg_(2,lag),this%jlg_(2,lag),this%klg_(2,lag))*xx2*yy2*zz2 &
      + this%esplfld_(this%ilg_(3,lag),this%jlg_(2,lag),this%klg_(2,lag))*xx3*yy2*zz2 &
      + this%esplfld_(this%ilg_(4,lag),this%jlg_(2,lag),this%klg_(2,lag))*xx4*yy2*zz2 &
      + this%esplfld_(this%ilg_(1,lag),this%jlg_(3,lag),this%klg_(2,lag))*xx1*yy3*zz2 &
      + this%esplfld_(this%ilg_(2,lag),this%jlg_(3,lag),this%klg_(2,lag))*xx2*yy3*zz2 &
      + this%esplfld_(this%ilg_(3,lag),this%jlg_(3,lag),this%klg_(2,lag))*xx3*yy3*zz2 &
      + this%esplfld_(this%ilg_(4,lag),this%jlg_(3,lag),this%klg_(2,lag))*xx4*yy3*zz2 &
      + this%esplfld_(this%ilg_(1,lag),this%jlg_(4,lag),this%klg_(2,lag))*xx1*yy4*zz2 &
      + this%esplfld_(this%ilg_(2,lag),this%jlg_(4,lag),this%klg_(2,lag))*xx2*yy4*zz2 &
      + this%esplfld_(this%ilg_(3,lag),this%jlg_(4,lag),this%klg_(2,lag))*xx3*yy4*zz2 &
      + this%esplfld_(this%ilg_(4,lag),this%jlg_(4,lag),this%klg_(2,lag))*xx4*yy4*zz2 &
!  k = 3
      + this%esplfld_(this%ilg_(1,lag),this%jlg_(1,lag),this%klg_(3,lag))*xx1*yy1*zz3 &
      + this%esplfld_(this%ilg_(2,lag),this%jlg_(1,lag),this%klg_(3,lag))*xx2*yy1*zz3 &
      + this%esplfld_(this%ilg_(3,lag),this%jlg_(1,lag),this%klg_(3,lag))*xx3*yy1*zz3 &
      + this%esplfld_(this%ilg_(4,lag),this%jlg_(1,lag),this%klg_(3,lag))*xx4*yy1*zz3 &
      + this%esplfld_(this%ilg_(1,lag),this%jlg_(2,lag),this%klg_(3,lag))*xx1*yy2*zz3 &
      + this%esplfld_(this%ilg_(2,lag),this%jlg_(2,lag),this%klg_(3,lag))*xx2*yy2*zz3 &
      + this%esplfld_(this%ilg_(3,lag),this%jlg_(2,lag),this%klg_(3,lag))*xx3*yy2*zz3 &
      + this%esplfld_(this%ilg_(4,lag),this%jlg_(2,lag),this%klg_(3,lag))*xx4*yy2*zz3 &
      + this%esplfld_(this%ilg_(1,lag),this%jlg_(3,lag),this%klg_(3,lag))*xx1*yy3*zz3 &
      + this%esplfld_(this%ilg_(2,lag),this%jlg_(3,lag),this%klg_(3,lag))*xx2*yy3*zz3 &
      + this%esplfld_(this%ilg_(3,lag),this%jlg_(3,lag),this%klg_(3,lag))*xx3*yy3*zz3 &
      + this%esplfld_(this%ilg_(4,lag),this%jlg_(3,lag),this%klg_(3,lag))*xx4*yy3*zz3 &
      + this%esplfld_(this%ilg_(1,lag),this%jlg_(4,lag),this%klg_(3,lag))*xx1*yy4*zz3 &
      + this%esplfld_(this%ilg_(2,lag),this%jlg_(4,lag),this%klg_(3,lag))*xx2*yy4*zz3 &
      + this%esplfld_(this%ilg_(3,lag),this%jlg_(4,lag),this%klg_(3,lag))*xx3*yy4*zz3 &
      + this%esplfld_(this%ilg_(4,lag),this%jlg_(4,lag),this%klg_(3,lag))*xx4*yy4*zz3 &
!  k = 4
      + this%esplfld_(this%ilg_(1,lag),this%jlg_(1,lag),this%klg_(4,lag))*xx1*yy1*zz4 &
      + this%esplfld_(this%ilg_(2,lag),this%jlg_(1,lag),this%klg_(4,lag))*xx2*yy1*zz4 &
      + this%esplfld_(this%ilg_(3,lag),this%jlg_(1,lag),this%klg_(4,lag))*xx3*yy1*zz4 &
      + this%esplfld_(this%ilg_(4,lag),this%jlg_(1,lag),this%klg_(4,lag))*xx4*yy1*zz4 &
      + this%esplfld_(this%ilg_(1,lag),this%jlg_(2,lag),this%klg_(4,lag))*xx1*yy2*zz4 &
      + this%esplfld_(this%ilg_(2,lag),this%jlg_(2,lag),this%klg_(4,lag))*xx2*yy2*zz4 &
      + this%esplfld_(this%ilg_(3,lag),this%jlg_(2,lag),this%klg_(4,lag))*xx3*yy2*zz4 &
      + this%esplfld_(this%ilg_(4,lag),this%jlg_(2,lag),this%klg_(4,lag))*xx4*yy2*zz4 &
      + this%esplfld_(this%ilg_(1,lag),this%jlg_(3,lag),this%klg_(4,lag))*xx1*yy3*zz4 &
      + this%esplfld_(this%ilg_(2,lag),this%jlg_(3,lag),this%klg_(4,lag))*xx2*yy3*zz4 &
      + this%esplfld_(this%ilg_(3,lag),this%jlg_(3,lag),this%klg_(4,lag))*xx3*yy3*zz4 &
      + this%esplfld_(this%ilg_(4,lag),this%jlg_(3,lag),this%klg_(4,lag))*xx4*yy3*zz4 &
      + this%esplfld_(this%ilg_(1,lag),this%jlg_(4,lag),this%klg_(4,lag))*xx1*yy4*zz4 &
      + this%esplfld_(this%ilg_(2,lag),this%jlg_(4,lag),this%klg_(4,lag))*xx2*yy4*zz4 &
      + this%esplfld_(this%ilg_(3,lag),this%jlg_(4,lag),this%klg_(4,lag))*xx3*yy4*zz4 &
      + this%esplfld_(this%ilg_(4,lag),this%jlg_(4,lag),this%klg_(4,lag))*xx4*yy4*zz4
   ENDDO
#if defined(DO_HYBRIDoffl)
!$omp target update from(fp(1:np))
!$omp end target data
#endif

  END SUBROUTINE GPSplineInt_Interp3D
!-----------------------------------------------------------------
!-----------------------------------------------------------------


  SUBROUTINE GPSplineInt_PartUpdate2D(this,xp,yp,np)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : PartUpdate2D
!  DESCRIPTION: Carries update update of interpolation positions.
!               Should be called before any interpolations are
!               done. 
!
!               NOTE: All interp points _must_ reside in the sub-
!               domain represented by this task. No checking is
!               done to insure this, and bad indices will be computed
!               if this isn't done.
!  ARGUMENTS  : 
!    this     : 'this' class instance
!    xp,yp    : particle x, y positions, or interpolation points
!               relative to grid set in constructor call 
!    np       : no. interp points or particles.
!-----------------------------------------------------------------
    IMPLICIT NONE
    CLASS(GPSplineInt)                     :: this
    REAL(KIND=GP),INTENT(IN),DIMENSION(*) :: xp,yp
    INTEGER      ,INTENT(IN)              :: np
    INTEGER                               :: j,nx,ny
    LOGICAL                               :: bok,btmp

    ! Compute interval-normalized positions and
    ! indices into control point array in x, y directions:
    nx = this%ldims_(1)
    ny = this%ldims_(2)
 
    ! x-coords:
!$omp parallel do 
    DO j = 1, np
       this%ilg_(1,j) = xp(j)*this%dxi_(1)
       this%xrk_  (j) = xp(j)*this%dxi_(1) - real(this%ilg_(1,j),kind=GP)
       this%ilg_(2,j) = this%ilg_(1,j) + 1
       this%ilg_(3,j) = modulo(this%ilg_(2,j),nx) + 1
       this%ilg_(4,j) = modulo(this%ilg_(3,j),nx) + 1
       this%ilg_(1,j) = modulo(nx+this%ilg_(2,j)-2,nx) + 1
    ENDDO
      
    ! y-coords:
    ! Note that the yp are on the _global_ grid, so
    ! we have to restrict to this MPI task's 'slab' 
    ! before computing local indices and spline coordinates....
    ! 
    ! First, do a check:
    bok = .true.
!!!$omp parallel do private(btmp)
    DO j = 1, np
      btmp = (yp(j).GE.this%xbnds_(2,1).AND.yp(j).LT.this%xbnds_(2,2))
!!!$omp critical
      bok = bok .AND. btmp
!!!$omp end critical
    ENDDO
    IF ( .NOT. bok ) THEN
      WRITE(*,*) 'GPSplineInt::PartUpdate2D: Invalid particle y-range'
      STOP
    ENDIF
!$omp parallel do 
    DO j = 1, np
       this%jlg_(1,j) = (yp(j)-this%xbnds_(2,1))*this%dxi_(2)
       this%yrk_  (j) = (yp(j)-this%xbnds_(2,1))*this%dxi_(2) &
                      - real(this%jlg_(1,j),kind=GP)
       this%jlg_(2,j) = this%jlg_(1,j) + 1
       this%jlg_(3,j) = this%jlg_(2,j) + 1
       this%jlg_(4,j) = this%jlg_(3,j) + 1
       this%jlg_(1,j) = this%jlg_(2,j) - 1
    ENDDO

  END SUBROUTINE GPSplineInt_PartUpdate2D
!-----------------------------------------------------------------
!-----------------------------------------------------------------


  SUBROUTINE GPSplineInt_PartUpdate3D(this,xp,yp,zp,np)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : PartUpdate3D
!  DESCRIPTION: Carries update update of interpolation positions.
!               Should be called before any interpolations are
!               done. 
!
!               NOTE: All interp points _must_ reside in the sub-
!               domain represented by this task. No checking is
!               done to insure this, and bad indices will be computed
!               if this isn't done.
!  ARGUMENTS  : 
!    this     : 'this' class instance
!    xp,yp,zp : particle x,y,z positions, or interpolation points
!               relative to grid set in constructor call 
!    np       : no. interp points or particles.
!-----------------------------------------------------------------
    IMPLICIT NONE
    CLASS(GPSplineInt)                    :: this
    REAL(KIND=GP),INTENT(IN),DIMENSION(*) :: xp,yp,zp
    INTEGER      ,INTENT(IN)              :: np
    INTEGER                               :: j,kmax,kmin,nx,ny,nz
    LOGICAL                               :: bok,btmp

    ! Compute interval-normalized positions and
    ! indices into control point array in x, y directions:
    nx = this%ldims_(1)
    ny = this%ldims_(2)
    nz = this%ldims_(3)
    kmax = nz+this%gpcomm_%GetNumGhost()-1
    kmin = this%gpcomm_%GetNumGhost()-1
#if defined(DO_HYBRIDoffl)
!$omp target update to(this%ilg_,this%jlg_,this%klg_,                 &
!$omp                  this%xrk_,this%yrk_,this%zrk_) device(targetdev)
#endif
    ! x-coords:
#if defined(DO_HYBRIDoffl)
!$omp target teams distribute parallel do map(to:np,nx,this%dxi_,this%xbnds_,xp(1:np))&
!$omp                                     map(this%ilg_,this%xrk_) device(targetdev)
#else
!$omp parallel do 
#endif
    DO j = 1, np
      this%ilg_(1,j) = (xp(j)-this%xbnds_(1,1))*this%dxi_(1)
      this%xrk_  (j) = (xp(j)-this%xbnds_(1,1))*this%dxi_(1) - real(this%ilg_(1,j),kind=GP)
      this%ilg_(2,j) = modulo(this%ilg_(1,j),nx) + 1
      this%ilg_(3,j) = modulo(this%ilg_(2,j),nx) + 1
      this%ilg_(4,j) = modulo(this%ilg_(3,j),nx) + 1
      this%ilg_(1,j) = modulo(nx+this%ilg_(2,j)-2,nx) + 1 
    ENDDO
      
    ! y-coords:
#if defined(DO_HYBRIDoffl)
!$omp target teams distribute parallel do map(to:np,ny,this%dxi_,this%xbnds_,yp(1:np))&
!$omp                                     map(this%jlg_,this%yrk_) device(targetdev)
#else
!$omp parallel do 
#endif
    DO j = 1, np
      this%jlg_(1,j) = (yp(j)-this%xbnds_(2,1))*this%dxi_(2)
      this%yrk_  (j) = (yp(j)-this%xbnds_(2,1))*this%dxi_(2) - real(this%jlg_(1,j),kind=GP)
      this%jlg_(2,j) = modulo(this%jlg_(1,j),ny) + 1 
      this%jlg_(3,j) = modulo(this%jlg_(2,j),ny) + 1
      this%jlg_(4,j) = modulo(this%jlg_(3,j),ny) + 1
      this%jlg_(1,j) = modulo(ny+this%jlg_(2,j)-2,ny) + 1 
    ENDDO

    ! z-coords:
    ! Note that the zp are on the _global_ grid, so
    ! we have to restrict to this MPI task's slab 
    ! before computing local indices and spline coordinates....
    ! 
    ! First, do a check:
    bok = .true.
!!!$omp parallel do private(btmp)
    DO j = 1, np
      btmp = (zp(j).GE.this%xbnds_(3,1).AND.zp(j).LT.this%xbnds_(3,2))
!!!$omp critical
      bok = bok .AND. btmp
!!!$omp end critical
      IF ( .NOT. bok ) THEN
        WRITE(*,*) this%rank_, ' GPSplineInt::PartUpdate3D: Invalid particle z-range'
        WRITE(*,*) this%rank_, ' GPSplineInt::zbnd_0=',this%xbnds_(3,1),';  zbnd_1=',this%xbnds_(3,2), 'zp=',zp(j)
        STOP
      ENDIF
    ENDDO
#if defined(DO_HYBRIDoffl)
!$omp target teams distribute parallel do map(to:np,km,this%dxi_,this%xbnds_,zp(1:np))&
!$omp                                     map(this%klg_,this%zrk_) device(targetdev)
#else
!$omp parallel do 
#endif
    DO j = 1, np
      this%klg_(1,j) = (zp(j)-this%xbnds_(3,1))*this%dxi_(3)
      this%klg_(1,j) = max(min(this%klg_(1,j),kmax),kmin)
      this%zrk_  (j) = (zp(j)-this%xbnds_(3,1))*this%dxi_(3) &
                     - real(this%klg_(1,j),kind=GP)
      this%klg_(2,j) = this%klg_(1,j) + 1
      this%klg_(3,j) = this%klg_(2,j) + 1
      this%klg_(4,j) = this%klg_(3,j) + 1
    ENDDO
#if defined(DO_HYBRIDoffl)
!$omp target update from(this%ilg_,this%jlg_,this%klg_,this%xrk_,this%yrk_,this%zrk_)&
!$omp               device(targetdev)
#endif

  END SUBROUTINE GPSplineInt_PartUpdate3D
!-----------------------------------------------------------------
!-----------------------------------------------------------------


  SUBROUTINE GPSplineInt_SetDeriv(this,idir,ido)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : SetDeriv
!  DESCRIPTION: Sets flag to interpolate derivative of the field
!               in a given direction. Must be called before call to
!               Interp2(3)D is made.
!  ARGUMENTS  : 
!    this     : 'this' class instance
!    idir     : coord. direction (1,2, ... rank)
!    ido      : >0: do derivative; <= 0: don't do derivative
!-----------------------------------------------------------------

    IMPLICIT NONE
    CLASS(GPSplineInt)                      :: this
    INTEGER,INTENT(IN)                      :: idir,ido


    IF ( idir .LT. 1 .OR. idir .GT. this%rank_ ) THEN
      WRITE(*,*) 'GPSplineInt::SetDeriv: Invalid coordinate direction'
      STOP
    ENDIF
    IF ( ido .GT. 0 ) THEN
      this%ider_(idir) = 1
    ELSE
      this%ider_(idir) = 0
    ENDIF

  END SUBROUTINE GPSplineInt_SetDeriv
!-----------------------------------------------------------------
!-----------------------------------------------------------------


  SUBROUTINE GPSplineInt_DoAlloc(this)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : DoAlloc
!  DESCRIPTION: performs allocations based on member data
!  ARGUMENTS  : 
!    this     : 'this' class instance
!-----------------------------------------------------------------
    IMPLICIT NONE
    CLASS(GPSplineInt)                     :: this
    INTEGER                                :: nzg
!$  INTEGER(c_int)                         :: rc
!$  TYPE(c_ptr)                            :: hptr
  
    CALL GPSplineInt_DoDealloc(this)
    nzg = this%gpcomm_%GetNumGhost()

    IF ( this%rank_ .EQ. 2 ) THEN
      ALLOCATE( this%esplfld2_(this%ldims_(1),(this%ldims_(2)+2*nzg)) )
    ELSE IF ( this%rank_ .EQ. 3 ) THEN
      ALLOCATE(this%esplfld_(this%ldims_(1),this%ldims_(2),(this%ldims_(3)+2*nzg)) )
#if defined(DO_HYBRIDoffl)
!$    this%dptresplfld = omp_target_alloc(this%ldims_(1)*this%ldims_(2)*(this%ldims_(3) &
!$                                      +2*nzg)*C_SIZEOF(this%esplfld_(1,1,1)),targetdev)
!$    rc = omp_target_associate_ptr(C_LOC(this%esplfld_),this%dptresplfld,   &
!$                      this%ldims_(1)*this%ldims_(2)*(this%ldims_(3)+2*nzg) &
!$                      *C_SIZEOF(this%esplfld_(1,1,1)),0_c_size_t,targetdev)
#endif
    ENDIF

    ALLOCATE(this%ax_   (this%nd_(1)) )
    ALLOCATE(this%bx_   (this%nd_(1)) )
    ALLOCATE(this%cx_   (this%nd_(1)) )
    ALLOCATE(this%betx_ (this%nd_(1)) )
    ALLOCATE(this%gamx_ (this%nd_(1)) )
    ALLOCATE(this%px_   (this%nd_(1)) )
    ALLOCATE(this%xxx_  (this%nd_(1)) )
#if defined(DO_HYBRIDoffl)
!$  this%dptrax = omp_target_alloc(this%nd_(1)*C_SIZEOF(this%ax_(1)),targetdev)
!$  rc = omp_target_associate_ptr(C_LOC(this%ax_),this%dptrax,          &
!$                      this%nd_(1)*C_SIZEOF(this%ax_(1)),0_c_size_t,targetdev)
!$  this%dptrpx = omp_target_alloc(this%nd_(1)*C_SIZEOF(this%px_(1)),targetdev)
!$  rc = omp_target_associate_ptr(C_LOC(this%px_),this%dptrpx,          &
!$                      this%nd_(1)*C_SIZEOF(this%px_(1)),0_c_size_t,targetdev)
!$  this%dptrxxx = omp_target_alloc(this%nd_(1)*C_SIZEOF(this%xxx_(1)),targetdev)
!$  rc = omp_target_associate_ptr(C_LOC(this%xxx_),this%dptrxxx,          &
!$                       this%nd_(1)*C_SIZEOF(this%xxx_(1)),0_c_size_t,targetdev)
!$  this%dptrbetx = omp_target_alloc(this%nd_(1)*C_SIZEOF(this%betx_(1)),targetdev)
!$  rc = omp_target_associate_ptr(C_LOC(this%betx_),this%dptrbetx,          &
!$                        this%nd_(1)*C_SIZEOF(this%betx_(1)),0_c_size_t,targetdev)
!$  this%dptrgamx = omp_target_alloc(this%nd_(1)*C_SIZEOF(this%gamx_(1)),targetdev)
!$  rc = omp_target_associate_ptr(C_LOC(this%gamx_),this%dptrgamx,          &
!$                        this%nd_(1)*C_SIZEOF(this%gamx_(1)),0_c_size_t,targetdev)
#endif

    ALLOCATE(this%ay_   (this%nd_(2)) )
    ALLOCATE(this%by_   (this%nd_(2)) )
    ALLOCATE(this%cy_   (this%nd_(2)) )
    ALLOCATE(this%bety_ (this%nd_(2)) )
    ALLOCATE(this%gamy_ (this%nd_(2)) )
    ALLOCATE(this%py_   (this%nd_(2)) )
    ALLOCATE(this%xxy_  (this%nd_(2)) )
#if defined(DO_HYBRIDoffl)
!$  this%dptray = omp_target_alloc(this%nd_(2)*C_SIZEOF(this%ay_(1)),targetdev)
!$  rc = omp_target_associate_ptr(C_LOC(this%ay_),this%dptray,          &
!$                      this%nd_(2)*C_SIZEOF(this%ay_(1)),0_c_size_t,targetdev)
!$  this%dptrpy = omp_target_alloc(this%nd_(2)*C_SIZEOF(this%py_(1)),targetdev)
!$  rc = omp_target_associate_ptr(C_LOC(this%py_),this%dptrpy,          &
!$                      this%nd_(2)*C_SIZEOF(this%py_(1)),0_c_size_t,targetdev)
!$  this%dptrxxy = omp_target_alloc(this%nd_(2)*C_SIZEOF(this%xxy_(1)),targetdev)
!$  rc = omp_target_associate_ptr(C_LOC(this%xxy_),this%dptrxxy,          &
!$                       this%nd_(2)*C_SIZEOF(this%xxy_(1)),0_c_size_t,targetdev)
!$  this%dptrbety = omp_target_alloc(this%nd_(2)*C_SIZEOF(this%bety_(1)),targetdev)
!$  rc = omp_target_associate_ptr(C_LOC(this%bety_),this%dptrbety,          &
!$                        this%nd_(2)*C_SIZEOF(this%bety_(1)),0_c_size_t,targetdev)
!$  this%dptrgamy = omp_target_alloc(this%nd_(2)*C_SIZEOF(this%gamy_(1)),targetdev)
!$  rc = omp_target_associate_ptr(C_LOC(this%gamy_),this%dptrgamy,          &
!$                        this%nd_(2)*C_SIZEOF(this%gamy_(1)),0_c_size_t,targetdev)
#endif

    ALLOCATE(this%wrkl_(9,this%maxint_))
    ALLOCATE(this%ilg_(4,this%maxint_))
    ALLOCATE(this%jlg_(4,this%maxint_))
    ALLOCATE(this%klg_(4,this%maxint_))
    ALLOCATE(this%xrk_(this%maxint_))
    ALLOCATE(this%yrk_(this%maxint_))
    ALLOCATE(this%zrk_(this%maxint_))
#if defined(DO_HYBRIDoffl)
!$  this%dptrwrkl = omp_target_alloc(9*this%maxint_*C_SIZEOF(this%wrkl_(1,1)),    &
!$                                   targetdev)
!$  rc = omp_target_associate_ptr(C_LOC(this%wrkl_),this%dptrwrkl,9*this%maxint_* &
!$                                C_SIZEOF(this%wrkl_(1,1)),0_c_size_t,targetdev)
!$  this%dptrilg = omp_target_alloc(4*this%maxint_*C_SIZEOF(this%ilg_(1,1)),targetdev)
!$  rc = omp_target_associate_ptr(C_LOC(this%ilg_),this%dptrilg,4*this%maxint_*   &
!$                                C_SIZEOF(this%ilg_(1,1)),0_c_size_t,targetdev)
!$  this%dptrjlg = omp_target_alloc(4*this%maxint_*C_SIZEOF(this%jlg_(1,1)),targetdev)
!$  rc = omp_target_associate_ptr(C_LOC(this%jlg_),this%dptrjlg,4*this%maxint_*   &
!$                                C_SIZEOF(this%jlg_(1,1)),0_c_size_t,targetdev)
!$  this%dptrklg = omp_target_alloc(4*this%maxint_*C_SIZEOF(this%klg_(1,1)),targetdev)
!$  rc = omp_target_associate_ptr(C_LOC(this%klg_),this%dptrklg,4*this%maxint_*   &
!$                                C_SIZEOF(this%klg_(1,1)),0_c_size_t,targetdev)
!$  this%dptrxrk = omp_target_alloc(this%maxint_*C_SIZEOF(this%xrk_(1)),targetdev)
!$  rc = omp_target_associate_ptr(C_LOC(this%xrk_),this%dptrxrk,this%maxint_*   &
!$                                C_SIZEOF(this%xrk_(1)),0_c_size_t,targetdev)
!$  this%dptryrk = omp_target_alloc(this%maxint_*C_SIZEOF(this%yrk_(1)),targetdev)
!$  rc = omp_target_associate_ptr(C_LOC(this%yrk_),this%dptryrk,this%maxint_*   &
!$                                C_SIZEOF(this%yrk_(1)),0_c_size_t,targetdev)
!$  this%dptrzrk = omp_target_alloc(this%maxint_*C_SIZEOF(this%zrk_(1)),targetdev)
!$  rc = omp_target_associate_ptr(C_LOC(this%zrk_),this%dptrzrk,this%maxint_*   &
!$                                C_SIZEOF(this%zrk_(1)),0_c_size_t,targetdev)
#endif

    IF ( this%rank_ .GT. 2 ) THEN
    ALLOCATE(this%az_   (this%nd_(3)) )
    ALLOCATE(this%bz_   (this%nd_(3)) )
    ALLOCATE(this%cz_   (this%nd_(3)) )
    ALLOCATE(this%betz_ (this%nd_(3)) )
    ALLOCATE(this%gamz_ (this%nd_(3)) )
    ALLOCATE(this%pz_   (this%nd_(3)) )
    ALLOCATE(this%xxz_  (this%nd_(3)) )
#if defined(DO_HYBRIDoffl)
!$  this%dptraz = omp_target_alloc(this%nd_(3)*C_SIZEOF(this%az_(1)),targetdev)
!$  rc = omp_target_associate_ptr(C_LOC(this%az_),this%dptraz,          &
!$                      this%nd_(3)*C_SIZEOF(this%az_(1)),0_c_size_t,targetdev)
!$  this%dptrpz = omp_target_alloc(this%nd_(3)*C_SIZEOF(this%pz_(1)),targetdev)
!$  rc = omp_target_associate_ptr(C_LOC(this%pz_),this%dptrpz,          &
!$                      this%nd_(3)*C_SIZEOF(this%pz_(1)),0_c_size_t,targetdev)
!$  this%dptrxxz = omp_target_alloc(this%nd_(3)*C_SIZEOF(this%xxz_(1)),targetdev)
!$  rc = omp_target_associate_ptr(C_LOC(this%xxz_),this%dptrxxz,          &
!$                       this%nd_(3)*C_SIZEOF(this%xxz_(1)),0_c_size_t,targetdev)
!$  this%dptrbetz = omp_target_alloc(this%nd_(3)*C_SIZEOF(this%betz_(1)),targetdev)
!$  rc = omp_target_associate_ptr(C_LOC(this%betz_),this%dptrbetz,          &
!$                        this%nd_(3)*C_SIZEOF(this%betz_(1)),0_c_size_t,targetdev)
!$  this%dptrgamz = omp_target_alloc(this%nd_(3)*C_SIZEOF(this%gamz_(1)),targetdev)
!$  rc = omp_target_associate_ptr(C_LOC(this%gamz_),this%dptrgamz,          &
!$                        this%nd_(3)*C_SIZEOF(this%gamz_(1)),0_c_size_t,targetdev)
#endif
    ENDIF

  END SUBROUTINE GPSplineInt_DoAlloc
!-----------------------------------------------------------------
!-----------------------------------------------------------------


  SUBROUTINE GPSplineInt_DoDealloc(this)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : DoDealloc
!  DESCRIPTION: performs de-allocation of all class quantities
!  ARGUMENTS  : 
!    this     : 'this' class instance
!-----------------------------------------------------------------

    IMPLICIT NONE
    CLASS(GPSplineInt)                      :: this
!$  INTEGER(c_int)                             :: rc

    IF ( ASSOCIATED (this%esplfld_  ) ) THEN
          DEALLOCATE (this%esplfld_)
#if defined(DO_HYBRIDoffl)
!$        call omp_target_free(this%dptresplfld,targetdev)
!$        rc = omp_target_disassociate_ptr(C_LOC(this%esplfld_),targetdev)
#endif
    ENDIF
    IF ( ASSOCIATED (this%esplfld2_ ) ) DEALLOCATE(this%esplfld2_)

    IF ( ASSOCIATED    (this%ax_   ) ) THEN
          DEALLOCATE      (this%ax_)
#if defined(DO_HYBRIDoffl)
!$        call omp_target_free(this%dptrax,targetdev)
!$        rc = omp_target_disassociate_ptr(C_LOC(this%ax_),targetdev)
#endif
    ENDIF
    IF ( ASSOCIATED    (this%bx_   ) ) DEALLOCATE      (this%bx_)
    IF ( ASSOCIATED    (this%cx_   ) ) DEALLOCATE      (this%cx_)
    IF ( ASSOCIATED    (this%px_   ) ) THEN
          DEALLOCATE      (this%px_)
#if defined(DO_HYBRIDoffl)
!$        call omp_target_free(this%dptrpx,targetdev)
!$        rc = omp_target_disassociate_ptr(C_LOC(this%px_),targetdev)
#endif
    ENDIF
    IF ( ASSOCIATED    (this%gamx_ ) ) THEN
          DEALLOCATE    (this%gamx_)
#if defined(DO_HYBRIDoffl)
!$        call omp_target_free(this%dptrgamx,targetdev)
!$        rc = omp_target_disassociate_ptr(C_LOC(this%gamx_),targetdev)
#endif
    ENDIF
    IF ( ASSOCIATED    (this%betx_ ) ) THEN
          DEALLOCATE    (this%betx_)
#if defined(DO_HYBRIDoffl)
!$        call omp_target_free(this%dptrbetx,targetdev)
!$        rc = omp_target_disassociate_ptr(C_LOC(this%betx_),targetdev)
#endif
    ENDIF
    IF ( ASSOCIATED    (this%xxx_  ) ) THEN
          DEALLOCATE     (this%xxx_)
#if defined(DO_HYBRIDoffl)
!$        call omp_target_free(this%dptrxxx,targetdev)
!$        rc = omp_target_disassociate_ptr(C_LOC(this%xxx_),targetdev)
#endif
    ENDIF

    IF ( ASSOCIATED    (this%ay_   ) ) THEN
          DEALLOCATE      (this%ay_)
#if defined(DO_HYBRIDoffl)
!$        call omp_target_free(this%dptray,targetdev)
!$        rc = omp_target_disassociate_ptr(C_LOC(this%ay_),targetdev)
#endif
    ENDIF
    IF ( ASSOCIATED    (this%by_   ) ) DEALLOCATE      (this%by_)
    IF ( ASSOCIATED    (this%cy_   ) ) DEALLOCATE      (this%cy_)
    IF ( ASSOCIATED    (this%py_   ) ) THEN
          DEALLOCATE      (this%py_)
#if defined(DO_HYBRIDoffl)
!$        call omp_target_free(this%dptrpy,targetdev)
!$        rc = omp_target_disassociate_ptr(C_LOC(this%py_),targetdev)
#endif
    ENDIF
    IF ( ASSOCIATED    (this%gamy_ ) ) THEN
          DEALLOCATE    (this%gamy_)
#if defined(DO_HYBRIDoffl)
!$        call omp_target_free(this%dptrgamy,targetdev)
!$        rc = omp_target_disassociate_ptr(C_LOC(this%gamy_),targetdev)
#endif
    ENDIF
    IF ( ASSOCIATED    (this%bety_ ) ) THEN
          DEALLOCATE    (this%bety_)
#if defined(DO_HYBRIDoffl)
!$        call omp_target_free(this%dptrbety,targetdev)
!$        rc = omp_target_disassociate_ptr(C_LOC(this%bety_),targetdev)
#endif
    ENDIF
    IF ( ASSOCIATED    (this%xxy_  ) ) THEN
          DEALLOCATE     (this%xxy_)
#if defined(DO_HYBRIDoffl)
!$        call omp_target_free(this%dptrxxy,targetdev)
!$        rc = omp_target_disassociate_ptr(C_LOC(this%xxy_),targetdev)
#endif
    ENDIF

    IF ( ASSOCIATED    (this%az_   ) ) THEN
          DEALLOCATE      (this%az_)
#if defined(DO_HYBRIDoffl)
!$        call omp_target_free(this%dptraz,targetdev)
!$        rc = omp_target_disassociate_ptr(C_LOC(this%az_),targetdev)
#endif
    ENDIF
    IF ( ASSOCIATED    (this%bz_   ) ) DEALLOCATE      (this%bz_)
    IF ( ASSOCIATED    (this%cz_   ) ) DEALLOCATE      (this%cz_)
    IF ( ASSOCIATED    (this%pz_   ) ) THEN
          DEALLOCATE      (this%pz_)
#if defined(DO_HYBRIDoffl)
!$        call omp_target_free(this%dptrpz,targetdev)
!$        rc = omp_target_disassociate_ptr(C_LOC(this%pz_),targetdev)
#endif
    ENDIF
    IF ( ASSOCIATED    (this%gamz_ ) ) THEN 
          DEALLOCATE    (this%gamz_)
#if defined(DO_HYBRIDoffl)
!$        call omp_target_free(this%dptrgamz,targetdev)
!$        rc = omp_target_disassociate_ptr(C_LOC(this%gamz_),targetdev)
#endif
    ENDIF
    IF ( ASSOCIATED    (this%betz_ ) ) THEN 
          DEALLOCATE    (this%betz_)
#if defined(DO_HYBRIDoffl)
!$        call omp_target_free(this%dptrbetz,targetdev)
!$        rc = omp_target_disassociate_ptr(C_LOC(this%betz_),targetdev)
#endif
    ENDIF
    IF ( ASSOCIATED    (this%xxz_  ) ) THEN
          DEALLOCATE     (this%xxz_)
#if defined(DO_HYBRIDoffl)
!$        call omp_target_free(this%dptrxxz,targetdev)
!$        rc = omp_target_disassociate_ptr(C_LOC(this%xxz_),targetdev)
#endif
    ENDIF

    IF ( ASSOCIATED      (this%wrkl_) ) THEN
          DEALLOCATE    (this%wrkl_)
#if defined(DO_HYBRIDoffl)
!$        call omp_target_free(this%dptrwrkl,targetdev)
!$        rc = omp_target_disassociate_ptr(C_LOC(this%wrkl_),targetdev)
#endif
    ENDIF
    IF ( ASSOCIATED       (this%ilg_) ) THEN
          DEALLOCATE     (this%ilg_)
#if defined(DO_HYBRIDoffl)
!$        call omp_target_free(this%dptrilg,targetdev)
!$        rc = omp_target_disassociate_ptr(C_LOC(this%ilg_),targetdev)
#endif
    ENDIF
    IF ( ASSOCIATED       (this%jlg_) ) THEN
          DEALLOCATE     (this%jlg_)
#if defined(DO_HYBRIDoffl)
!$        call omp_target_free(this%dptrjlg,targetdev)
!$        rc = omp_target_disassociate_ptr(C_LOC(this%jlg_),targetdev)
#endif
    ENDIF
    IF ( ASSOCIATED       (this%klg_) ) THEN
          DEALLOCATE     (this%klg_)
#if defined(DO_HYBRIDoffl)
!$        call omp_target_free(this%dptrklg,targetdev)
!$        rc = omp_target_disassociate_ptr(C_LOC(this%klg_),targetdev)
#endif
    ENDIF
    IF ( ASSOCIATED       (this%xrk_) ) THEN
          DEALLOCATE     (this%xrk_)
#if defined(DO_HYBRIDoffl)
!$        call omp_target_free(this%dptrxrk,targetdev)
!$        rc = omp_target_disassociate_ptr(C_LOC(this%xrk_),targetdev)
#endif
    ENDIF
    IF ( ASSOCIATED       (this%yrk_) ) THEN
            DEALLOCATE     (this%yrk_)
#if defined(DO_HYBRIDoffl)
!$        call omp_target_free(this%dptryrk,targetdev)
!$        rc = omp_target_disassociate_ptr(C_LOC(this%yrk_),targetdev)
#endif
    ENDIF
    IF ( ASSOCIATED       (this%zrk_) ) THEN
          DEALLOCATE     (this%zrk_)
#if defined(DO_HYBRIDoffl)
!$        call omp_target_free(this%dptrzrk,targetdev)
!$        rc = omp_target_disassociate_ptr(C_LOC(this%zrk_),targetdev)
#endif
    ENDIF

  END SUBROUTINE GPSplineInt_DoDealloc
!-----------------------------------------------------------------
!-----------------------------------------------------------------


  SUBROUTINE GPSplineInt_Init(this)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : Init
!  DESCRIPTION: Does initialization of object. Called by constructor
!  ARGUMENTS  : 
!    this     : 'this' class instance
!-----------------------------------------------------------------
    IMPLICIT NONE
    CLASS(GPSplineInt)        :: this
    INTEGER                   :: j

    CALL GPSplineInt_DoAlloc(this)

    ! Compute interval widths:
    DO j = 1, this%rank_
      this%dxi_(j) = 1.0_GP  !real(this%ldims_(j)-1,kind=GP)/ ( this%xbnds_(j,2) - this%xbnds_(j,1) )
    ENDDO

    IF (this%nd_(1).GT.2) THEN
    CALL GPSplineInt_MatInvQ(this,this%nd_(1),this%ax_,this%bx_,this%cx_,&
    this%px_,this%gamx_,this%betx_,this%xxx_,this%zetax_)
    END IF
    IF (this%nd_(2).GT.2) THEN
    CALL GPSplineInt_MatInvQ(this,this%nd_(2),this%ay_,this%by_,this%cy_,&
    this%py_,this%gamy_,this%bety_,this%xxy_,this%zetay_)
    END IF
    IF ( this%rank_ .GT. 2 ) THEN
    CALL GPSplineInt_MatInvQ(this,this%nd_(3),this%az_,this%bz_,this%cz_,&
    this%pz_,this%gamz_,this%betz_,this%xxz_,this%zetaz_)
    ENDIF

  END SUBROUTINE GPSplineInt_Init
!-----------------------------------------------------------------
!-----------------------------------------------------------------


  SUBROUTINE GPSplineInt_MatInvQ(this,n,a,b,c,p,gam,bet,xx,zeta)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : MatInvQ
!  DESCRIPTION: Computes quantities for matrix inversion
!  ARGUMENTS  : 
!    this     : 'this' class instance
!-----------------------------------------------------------------

    IMPLICIT NONE
    CLASS(GPSplineInt)                       :: this
    INTEGER      ,INTENT   (IN)              :: n
    INTEGER                                  :: i
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(n) :: a,b,c,p,gam,bet,xx
    REAL(KIND=GP),INTENT(INOUT)              :: zeta
    REAL(KIND=GP)                            :: beta,sixth,twotrd

!  Setup the arrays for the inversion:
        sixth = 1.0_GP/6.0_GP
        twotrd= 2.0_GP/3.0_GP
        DO i = 1, n
          a(i) = sixth
          c(i) = sixth
          b(i) = twotrd
        ENDDO
!
!  Initialize the other arrays :
        bet (1) = 1./b(1)
        p   (1) = a(1)*bet(1)
        xx  (1) = c(n)
        beta    = b(n)
!
        DO i= 2, n-2
          gam (i) = c(i-1)*bet(i-1)
          bet (i) = 1./(b(i)-a(i)*gam(i))
          p   (i) = -p(i-1)*a(i)*bet(i)
          beta    = beta - xx(i-1)*p(i-1)
          xx  (i) = -xx(i-1)*gam(i)
        ENDDO
!  ** n-1 **
        gam (n-1) = c(n-2)*bet(n-2)
        bet (n-1) = 1./(b(n-1)-a(n-1)*gam(n-1))
        gam   (n) = (c(n-1)-p(n-2)*a(n-1))*bet(n-1)
        zeta      = a(n) - xx(n-2)*gam(n-1)
        beta      = beta - xx(n-2)*p(n-2)
!  ** n  **
        bet   (n) = 1./(beta - zeta*gam(n))

  END SUBROUTINE GPSplineInt_MatInvQ
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GPSplineInt_CompSpline2D(this,field,tmp1,tmp2)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : GPSplineInt_CompSpline2D
!  DESCRIPTION: Computes the spline coefficients, by soving tridiagonal
!            system. Field must have dimensions set in contructor.
!            Spline coeffs are stored in member data array.
!            
!  ARGUMENTS  :
!    this    : 'this' class instance
!    field   : Field to be interpolated. Must have dimensions set in 
!              constructor. This field is overwirtten with tensor 
!              product of full set of spline coeffs, on exit
!    tmp1/2  : temp arrays of size n X (jend-jsta+1) as in
!              instantiating code.
!-----------------------------------------------------------------
    IMPLICIT NONE
    CLASS(GPSplineInt)                         :: this
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(*)  :: field,tmp1,tmp2
    INTEGER                                   :: nx,ny

    WRITE(*,*) 'GPSplineInt::CompSpline2D: Under construction!'
    STOP

  END SUBROUTINE GPSplineInt_CompSpline2D
!-----------------------------------------------------------------
!-----------------------------------------------------------------


  SUBROUTINE GPSplineInt_CompSpline3D(this,field,tmp1,tmp2)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : GPSplineInt_CompSpline3D
!  DESCRIPTION: Computes the spline coefficients, by soving tridiagonal
!               system. Field must have dimensions set in contructor.
!               Spline coeffs are stored in member data array.
!            
!  ARGUMENTS  :
!    this    : 'this' class instance
!    field   : Field to be interpolated. Must have dimensions set in 
!              constructor. On exit this field is overwirtten with 
!              (partial) tensor product of full set of spline coeffs
!    tmp1/2:   temp arrays of size n X n X (kend-ksta+1) as in
!              instantiating code.
!-----------------------------------------------------------------
    USE mpivars

    IMPLICIT NONE
    CLASS(GPSplineInt)                                :: this
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(this%ntot_) :: field
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(this%ntot_) :: tmp1,tmp2
    REAL(KIND=GP),DIMENSION(this%ttot_)               :: tmptr,tmpt2
    INTEGER                                           :: i,j,k
    INTEGER                                           :: jm,km
    INTEGER                                           :: nx,nxy,ny,nz
    INTEGER                                           :: ibnds(3,2)
    INTEGER                                           :: obnds(3,2)

    nx = this%ldims_(1)
    ny = this%ldims_(2)
    nz = this%ldims_(3)
    nxy = nx*ny

#if defined(DO_HYBRIDoffl)
!General offloading of class properties
!$omp target update to(this%ax_,this%ay_,this%az_,this%px_,this%py_,this%pz_,  &
!$omp                  this%betx_,this%bety_,this%betz_,this%gamx_,this%gamy_, &
!$omp                  this%gamz_,this%xxx_,this%xxy_,this%xxz_) device(targetdev)

!Offloading array necessary for x-y computation (field) and transposing (tmp2)
!$omp target data map(to:nx,ny,nz,nxy,field,this%zetax_,this%zetay_)           &
!$omp             map(alloc:tmp2) device(targetdev)
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!! x-field computation !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#if defined(DO_HYBRIDoffl)
!$omp target teams loop collapse(1) private(j,km,jm) map(this%betx_) device(targetdev)
#else
!$omp parallel do private(j,km,jm) 
#endif
    DO k=1,nz
      km = k-1
      DO j=1,ny
        jm = j-1
        tmp2 (1+jm*nx+km*nxy) = field(1+jm*nx+km*nxy)*this%betx_(1)
      ENDDO
    ENDDO
#if defined(DO_HYBRIDoffl)
!$omp target teams loop collapse(1) private(j,km,jm) device(targetdev)
#else
!$omp parallel do private(j,km,jm)
#endif
     DO k=1,nz
      km = k-1
      DO j=1,ny
        jm = j-1
        tmp2(nx+jm*nx+km*nxy) = field(nx+jm*nx+km*nxy)
      ENDDO
    ENDDO
#if defined(DO_HYBRIDoffl)
!$omp target teams loop collapse(2) private(i,jm) map(this%ax_,this%betx_)& 
!$omp                                             device(targetdev)
    DO km=0,nz-1
#else
    DO k=1,nz
      km = k-1
#endif
      DO j=1,ny
        jm = j-1
        DO i=2,nx-2
           tmp2(i+jm*nx+km*nxy) =  &
           ( field(i+jm*nx+km*nxy) - this%ax_(i)*tmp2(i-1+jm*nx+km*nxy) )*this%betx_(i)  
        ENDDO
      ENDDO
    ENDDO
#if defined(DO_HYBRIDoffl)
!$omp target teams loop collapse(2) private(i,jm) map(this%xxx_) device(targetdev)
    DO km=0,nz-1
#else
    DO k=1,nz
      km = k-1
#endif
      DO j=1,ny
        jm = j-1
        DO i=2,nx-2
           tmp2(nx+jm*nx+km*nxy) = tmp2(nx+jm*nx+km*nxy) - this%xxx_(i-1)*tmp2(i-1+jm*nx+km*nxy)
        ENDDO
      ENDDO
    ENDDO
#if defined(DO_HYBRIDoffl)
!$omp target teams loop collapse(1) private(j,km,jm) map(this%ax_,this%zetax_,&
!$omp                 this%betx_,this%gamx_,this%xxx_) device(targetdev)
#endif
    DO k=1,nz
      km = k-1
      DO j=1,ny
        jm = j-1
!  ** n-1 **
        tmp2   (nx-1+jm*nx+km*nxy) = &
        (field(nx-1+jm*nx+km*nxy) - this%ax_(nx-1)*tmp2(nx-2+jm*nx+km*nxy))*this%betx_(nx-1)
        tmp2     (nx+jm*nx+km*nxy) = tmp2(nx+jm*nx+km*nxy) - this%xxx_(nx-2)*tmp2(nx-1+jm*nx+km*nxy)
!  ** n  **
        tmp2(nx+jm*nx+km*nxy) = (tmp2(nx+jm*nx+km*nxy) - tmp2(nx-1+jm*nx+km*nxy)*this%zetax_) &
                     * this%betx_(nx)
!  Backsubstitution phase :
        tmp2(nx-1+jm*nx+km*nxy) = tmp2(nx-1+jm*nx+km*nxy) - this%gamx_(nx)*tmp2(nx+jm*nx+km*nxy)
      ENDDO
    ENDDO
!
#if defined(DO_HYBRIDoffl)
!$omp target teams loop collapse(2) private(i,jm)                      &
!$omp                   map(this%px_,this%gamx_) device(targetdev)
    DO km=0,nz-1
#else
    DO k=1,nz
    km = k-1
#endif
      DO j=1,ny
        jm = j-1
        DO  i=nx-2,1,-1
        tmp2(i+jm*nx+km*nxy) = tmp2(i+jm*nx+km*nxy) &
                      - this%gamx_(i+1)*tmp2(i+1+jm*nx+km*nxy) - this%px_(i)*tmp2(nx+jm*nx+km*nxy)
        ENDDO
      ENDDO
    ENDDO

!  Copy splfld -> field :
#if defined(DO_HYBRIDoffl)
!$omp target teams loop collapse(2) private(i,jm) device(targetdev)
    DO km=0,nz-1
#else
!$omp parallel do private(i,j,km,jm)
    DO k=1,nz
      km = k-1
#endif
      DO j=1,ny
        jm = j-1
        DO i=1,nx
          field(i+jm*nx+km*nxy) = tmp2(i+jm*nx+km*nxy)
        ENDDO
      ENDDO
    ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!! y-field computation !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Compute the j-equation (coefficients in y) :
#if defined(DO_HYBRIDoffl)
!$omp target teams loop collapse(1) private(i,km) map(this%bety_) device(targetdev)
#else
!$omp parallel do private(i,km)
#endif
    DO k=1,nz
      km = k-1
      DO i=1,nx
        tmp2(i+   km*nxy) = field(i+km*nxy)*this%bety_(1)
      ENDDO
    ENDDO
#if defined(DO_HYBRIDoffl)
!$omp target teams loop collapse(1) private(i,km) device(targetdev)
#else
!$omp parallel do private(i,km)
#endif
    DO k=1,nz
      km = k-1
      DO i=1,nx
        tmp2(i+(ny-1)*nx+km*nxy) = field(i+(ny-1)*nx+km*nxy)
      ENDDO
    ENDDO
!!
#if defined(DO_HYBRIDoffl)
!$omp target teams loop collapse(2) private(j,jm) map(this%ay_,this%bety_)&
!$omp                                             device(targetdev)
    DO km=0,nz-1
      DO i=1,nx
        DO j=2,ny-2
          jm = j-1
#else
    DO k=1,nz
       km = k-1
       DO j=2,ny-2
         jm = j-1
         DO i=1,nx
#endif
           tmp2    (i+jm*nx+km*nxy) = &
           ( field(i+jm*nx+km*nxy) - this%ay_(j)*tmp2(i+(jm-1)*nx+km*nxy) )*this%bety_(j)
         ENDDO
       ENDDO
    ENDDO
#if defined(DO_HYBRIDoffl)
!$omp target teams loop collapse(2) private(j,jm) map(this%xxy_) device(targetdev)
    DO km=0,nz-1
      DO i=1,nx
         DO j=2,ny-2
           jm = j-1
#else
    DO k=1,nz
       km = k-1
       DO j=2,ny-2
         jm = j-1
         DO i=1,nx
#endif
           tmp2   (i+(ny-1)*nx+km*nxy) = tmp2(i+(ny-1)*nx+km*nxy) - this%xxy_(j-1)*tmp2(i+(jm-1)*nx+km*nxy)
         ENDDO
       ENDDO
    ENDDO
!  ** n-1 **
#if defined(DO_HYBRIDoffl)
!$omp target teams loop collapse(1) private(i,km) map(this%ay_,this%zetay_,&
!$omp                 this%bety_,this%gamy_,this%xxy_) device(targetdev)
#endif
    DO k=1,nz
      km = k-1
      DO i=1,nx
        tmp2    (i+(ny-2)*nx+km*nxy) = &
        (field(i+(ny-2)*nx+km*nxy) - this%ay_(ny-1)*tmp2(i+(ny-3)*nx+km*nxy))*this%bety_(ny-1)
        tmp2      (i+(ny-1)*nx+km*nxy) = tmp2(i+(ny-1)*nx+km*nxy) - this%xxy_(ny-2)*tmp2(i+(ny-2)*nx+km*nxy)
!  ** n  **
        tmp2(i+(ny-1)*nx+km*nxy) = (tmp2(i+(ny-1)*nx+km*nxy) - tmp2(i+(ny-2)*nx+km*nxy)*this%zetay_) &
                       * this%bety_(ny)
!  Backsubstitution phase :
        tmp2(i+(ny-2)*nx+km*nxy) = tmp2(i+(ny-2)*nx+km*nxy) - this%gamy_(ny)*tmp2(i+(ny-1)*nx+km*nxy)
      ENDDO
    ENDDO
!
#if defined(DO_HYBRIDoffl)
!$omp target teams loop collapse(2) private(j,jm)                      &
!$omp                   map(this%py_,this%gamy_) device(targetdev)
    DO km=0,nz-1
      DO i=1,nx
        DO j=ny-2,1,-1
          jm = j-1
#else
    DO k=1,nz
      km = k-1
      DO j=ny-2,1,-1
        jm = j-1
        DO i=1,nx
#endif
          tmp2(i+jm*nx+km*nxy) = tmp2(i+jm*nx+km*nxy)  &
                         - this%gamy_(j+1)*tmp2(i+(jm+1)*nx+km*nxy) - this%py_(j)*tmp2(i+(ny-1)*nx+km*nxy)
        ENDDO
      ENDDO
    ENDDO
#if defined(DO_HYBRIDoffl)
!$omp target update from(tmp2) device(targetdev)
!$omp end target data 
#endif

! Note tmp2 now contains the contributions for full tensor product
! field from x, y only....

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!! z-field computation !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Must transpose data so that it is z-y complete, in order to
! carry out computation in z:
!!    call gpcomm%GTranspose (vt,obnds,v ,ibnds,3) ! make z-y complete
!!    call gpcomm%GITranspose(vb,ibnds,vt,obnds,3) ! make x-y complete again

    ibnds(1,1) = 1
    ibnds(1,2) = this%nd_(1)
    ibnds(2,1) = 1
    ibnds(2,2) = this%nd_(2)
    ibnds(3,1) = 1 
    ibnds(3,2) = this%ldims_(3)
    obnds(1,1) = 1
    obnds(1,2) = this%nd_(3)
    obnds(2,1) = 1
    obnds(2,2) = this%nd_(2)
    obnds(3,1) = 1  
    obnds(3,2) = this%odims_(3)

    CALL GTStart(this%htransp_)
    CALL this%gpcomm_%GTranspose(tmptr,obnds,tmp2,ibnds,3,tmpt2)
    CALL GTAcc(this%htransp_)

    nx = this%odims_(1)
    ny = this%odims_(2)
    nz = this%odims_(3)
    nxy = nx*ny
#if defined(DO_HYBRIDoffl)
!Offload transposed array for z-computation (tmptr) and transposing (tmpt2)
!$omp target data map(to:nx,ny,nz,nxy,tmptr,this%zetaz_) map(alloc:tmpt2)      &
!$omp                                                    device(targetdev)
#endif
!  Compute the k-equation (coefficient in z) :
!  Note that the we work on the transposed system here...;
!  
!  tmpt2  <== spline coeff. field;
!  tmptr <== 'field', semi-'tensored' with spl. coeffs
#if defined(DO_HYBRIDoffl)
!$omp target teams loop collapse(1) private(j,km,jm) map(this%betz_) device(targetdev)
#else
!$omp parallel do private(j,jm,km)
#endif
    DO k=1,nz
      km = k-1
      DO j=1,ny
        jm = j-1
        tmpt2 (1+jm*nx+km*nxy) = tmptr(1+jm*nx+km*nxy)*this%betz_(1)
      ENDDO
    ENDDO
#if defined(DO_HYBRIDoffl)
!$omp target teams loop collapse(1) private(j,jm) device(targetdev)
#endif
    DO k=1,nz
      km = k-1
      DO j=1,ny
        jm = j-1
        tmpt2(nx+jm*nx+km*nxy) = tmptr(nx+jm*nx+km*nxy)
      ENDDO
    ENDDO
!
#if defined(DO_HYBRIDoffl)
!$omp target teams loop collapse(2) private(i,jm) map(this%az_,this%betz_)&
!$omp                                             device(targetdev)
    DO km=0,nz-1
#else
    DO k=1,nz
      km = k-1
#endif
      DO j=1,ny
        jm = j-1
        DO i=2,nx-2
          tmpt2 (i+jm*nx+km*nxy) =  &
          (tmptr(i+jm*nx+km*nxy) - this%az_(i)*tmpt2(i-1+jm*nx+km*nxy) )*this%betz_(i) 
        ENDDO
      ENDDO
    ENDDO
#if defined(DO_HYBRIDoffl)
!$omp target teams loop collapse(2) private(i,jm) map(this%xxz_) device(targetdev)
    DO km=0,nz-1
#else
    DO k=1,nz
      km = k-1
#endif
      DO j=1,ny
        jm = j-1
        DO i=2,nx-2
          tmpt2(nx+jm*nx+km*nxy) = &
          tmpt2(nx+jm*nx+km*nxy) - this%xxz_(i-1)*tmpt2(i-1+jm*nx+km*nxy)
        ENDDO
      ENDDO
    ENDDO
!
#if defined(DO_HYBRIDoffl)
!$omp target teams loop collapse(1) private(j,jm) map(this%az_,this%zetaz_,&
!$omp                   this%betz_,this%gamz_,this%xxz_) device(targetdev)
#endif
    DO k=1,nz
      km = k-1
      DO j=1,ny
        jm = j-1
!  ** n-1 **
        tmpt2 (nx-1+jm*nx+km*nxy) = &
        (tmptr(nx-1+jm*nx+km*nxy) - this%az_(nx-1)*tmpt2(nx-2+jm*nx+km*nxy))*this%betz_(nx-1)
        tmpt2 (nx+jm*nx+km*nxy)   = &
        tmpt2(nx+jm*nx+km*nxy) - this%xxz_(nx-2)*tmpt2(nx-1+jm*nx+km*nxy)
!  ** n  **
        tmpt2 (nx+jm*nx+km*nxy)   = &
        (tmpt2(nx+jm*nx+km*nxy) - tmpt2(nx-1+jm*nx+km*nxy)*this%zetaz_)*this%betz_(nx)
!  Backsubstitution phase :
        tmpt2(nx-1+jm*nx+km*nxy)  = &
        tmpt2(nx-1+jm*nx+km*nxy) - this%gamz_(nx)*tmpt2(nx+jm*nx+km*nxy)
      ENDDO
    ENDDO
!
#if defined(DO_HYBRIDoffl)
!$omp target teams loop collapse(2) private(i,km)                      &
!$omp                   map(this%pz_,this%gamz_) device(targetdev)
    DO jm=0,ny-1
#else
    DO j=1,ny
      jm = j-1
#endif
      DO k=1,nz
        km = k-1
          DO i=nx-2,1,-1
          tmpt2(i+jm*nx+km*nxy) = tmpt2(i+jm*nx+km*nxy) &
          - this%gamz_(i+1)*tmpt2(i+1+jm*nx+km*nxy) - this%pz_(i)*tmpt2(nx+jm*nx+km*nxy)
        ENDDO
      ENDDO
    ENDDO
#if defined(DO_HYBRIDoffl)
!$omp target update from(tmpt2) device(targetdev)
!$omp end target data
#endif
!
!  Transpose back to standard orientation:
    CALL GTStart(this%htransp_)
    CALL this%gpcomm_%GITranspose(field,ibnds,tmpt2,obnds,3,tmptr)
    CALL GTAcc(this%htransp_)
! 
! Spline coeff field is now slab-decomposed to be xy complete, and 
! distributed in z. To do interpolation, we must set up 
! spline coeffs in array extended in z by 2 control points
! so that spline interpolation can be done:
    CALL GTStart(this%hdataex_)
    CALL this%gpcomm_%SlabDataExchangeSF(this%esplfld_,field)
    CALL GTAcc(this%hdataex_)
!!write(*,*)'esplfld1=',this%esplfld_(1:10,1:10,1)
!!write(*,*)'esplfld2=',this%esplfld_(1:10,1:10,2)
!!write(*,*)'v1     =',field(1:100)
!!write(*,*)'v2     =',field(1:100+nxy)

   RETURN

  END SUBROUTINE GPSplineInt_CompSpline3D
!-----------------------------------------------------------------
!-----------------------------------------------------------------

 SUBROUTINE GPSplineInt_ResizeArrays(this,newmparts,onlyinc)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : Resize_Arrays
!  DESCRIPTION: Resize all arrays in the GPSplineInt class.
!               Communicator _must_ be resized separately (i.e., 
!               from GPart%ResizeArrays).
!               Not compatible with offloading.
!
!  ARGUMENTS  :
!    this     : 'this' class instance
!    newmparts: new number of particles
!    onlyinc  : if true, will only resize to increase array size
!-----------------------------------------------------------------
    IMPLICIT NONE
    CLASS(GPSplineInt),INTENT(INOUT)   :: this
    INTEGER           ,INTENT(IN)      :: newmparts
    LOGICAL           ,INTENT(IN)      :: onlyinc
    INTEGER                            :: n

    n = SIZE(this%wrkl_, 2)
    IF ((n.lt.newmparts).OR.((n.gt.newmparts).AND..NOT.onlyinc)) THEN
      CALL Resize_PIntArrayRank2(this%ilg_,newmparts,.false.)
      CALL Resize_PIntArrayRank2(this%jlg_,newmparts,.false.)
      CALL Resize_PIntArrayRank2(this%klg_,newmparts,.false.)
      CALL Resize_PArrayRank2(this%wrkl_  ,newmparts,.false.)
      CALL Resize_PArrayRank1(this%xrk_   ,newmparts,.false.)
      CALL Resize_PArrayRank1(this%yrk_   ,newmparts,.false.)
      CALL Resize_PArrayRank1(this%zrk_   ,newmparts,.false.)
      this%maxint_ = newmparts
    END IF
! Communicator resized from GPart class
    RETURN

 END SUBROUTINE GPSplineInt_ResizeArrays
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE Resize_PIntArrayRank2(a,new_size,docpy)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : Resize_PIntArrayRank2
!  DESCRIPTION: Resize input integer 2D pointer array to new size
!  ARGUMENTS  :
!    a       : array to be resized
!    new_size: new size for array (second axis)
!    docpy   : if true, keeps previous data in resized array
!-----------------------------------------------------------------
!$  USE threads

    IMPLICIT NONE
    INTEGER,INTENT(INOUT),POINTER,DIMENSION(:,:) :: a
    INTEGER      ,INTENT(IN)                     :: new_size
    LOGICAL      ,INTENT(IN)                     :: docpy
    INTEGER              ,POINTER,DIMENSION(:,:) :: temp
    INTEGER                                      :: i,n,shp(2)

    shp = SHAPE(a)
    ALLOCATE ( temp(shp(1),new_size) )

    IF (docpy) THEN
      n = shp(2)
      n = MIN(n,new_size)
!$omp parallel do if(n.gt.NMIN_OMP)
      DO i = 1,n
        temp(:,i) = a(:,i)
      END DO
    END IF
    DEALLOCATE(a)
    a => temp
!    CALL move_alloc(temp, a)

    RETURN
  END SUBROUTINE Resize_PIntArrayRank2
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE Resize_PArrayRank2(a,new_size,docpy)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : Resize_ArrayRank2
!  DESCRIPTION: Resize input 2D array to new size
!  ARGUMENTS  :
!    a       : array to be resized
!    new_size: new size for array (second axis)
!    docpy   : if true, keeps previous data in resized array
!-----------------------------------------------------------------
!$  USE threads

    IMPLICIT NONE
    REAL(KIND=GP),INTENT(INOUT),POINTER,DIMENSION(:,:) :: a
    INTEGER      ,INTENT(IN)                           :: new_size
    LOGICAL      ,INTENT(IN)                           :: docpy
    REAL(KIND=GP)              ,POINTER,DIMENSION(:,:) :: temp
    INTEGER                                            :: i,n,shp(2)

    shp = SHAPE(a)
    ALLOCATE ( temp(shp(1),new_size) )

    IF (docpy) THEN
      n = shp(2)
      n = MIN(n,new_size)
!$omp parallel do if(n.gt.NMIN_OMP)
      DO i = 1,n
        temp(:,i) = a(:,i)
      END DO
    END IF

    DEALLOCATE(a)
    a => temp
!    CALL move_alloc(temp, a)

    RETURN

  END SUBROUTINE Resize_PArrayRank2
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE Resize_PArrayRank1(a,new_size,docpy)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : Resize_ArrayRank1
!  DESCRIPTION: Resize input 1D array to new size
!  ARGUMENTS  :
!    a       : array to be resized
!    new_size: new size for array
!    docpy   : if true, keeps previous data in resized array
!-----------------------------------------------------------------
!$  USE threads

    IMPLICIT NONE
    REAL(KIND=GP),INTENT(INOUT),POINTER,DIMENSION(:) :: a
    INTEGER      ,INTENT(IN)                         :: new_size
    LOGICAL      ,INTENT(IN)                         :: docpy
    REAL(KIND=GP)              ,POINTER,DIMENSION(:) :: temp
    INTEGER                                          :: i,n

    ALLOCATE ( temp(new_size) )
    IF (docpy) THEN
      n = SIZE(a)
      n = MIN(n,new_size)
!$omp parallel do if(n.gt.NMIN_OMP)
      DO i = 1,n
        temp(i) = a(i)
      END DO
    END IF

    DEALLOCATE(a)
    a => temp

    RETURN
  END SUBROUTINE Resize_PArrayRank1
!-----------------------------------------------------------------
!-----------------------------------------------------------------

END MODULE class_GPSplineInt
