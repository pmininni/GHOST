!=================================================================
! GHOST GPSplineInt interpolation class
! The spline used (that on the 'momoents' or 2nd derivatives)
! as well as the tridiagonal solve for the moments, are
! derived from several sources including the presentation 
! given in the book 'The theory of splines and their application'
! 1967, by Ahlberg, Nilson, Walsh.
!
! Note that this class assumes a regularly spaced grid.
!
! 2013 A. Pumir (ENS, Lyon)
!      D. Rosenberg (NCCS: ORNL)
!
! 15 Jan 2013: Initial version
!=================================================================
MODULE class_GPSplineInt
      USE mpivars
      USE fprecision
      USE class_GPartComm
      IMPLICIT NONE

      TYPE, PUBLIC :: GPSplineInt
        PRIVATE
        ! Member data:
        REAL(KIND=GP),ALLOCATABLE,DIMENSION  (:,:)   :: esplfld2_
        REAL(KIND=GP),ALLOCATABLE,DIMENSION(:,:,:)   :: esplfld_
        REAL(KIND=GP),ALLOCATABLE,DIMENSION    (:)   :: ax_,bx_,betx_,cx_,gamx_,px_,xxx_
        REAL(KIND=GP),ALLOCATABLE,DIMENSION    (:)   :: ay_,by_,bety_,cy_,gamy_,py_,xxy_
        REAL(KIND=GP),ALLOCATABLE,DIMENSION    (:)   :: az_,bz_,betz_,cz_,gamz_,pz_,xxz_
        REAL(KIND=GP),ALLOCATABLE,DIMENSION    (:)   :: xrk_,yrk_,zrk_
        REAL(KIND=GP),ALLOCATABLE,DIMENSION  (:,:)   :: wrkl_
        REAL(KIND=GP)                                :: dxi_(3),xbnds_(3,2),zetax_,zetay_,zetaz_
        TYPE(GPartComm),POINTER                      :: gpcomm_
        INTEGER      ,ALLOCATABLE,DIMENSION  (:,:)   :: ilg_,jlg_,klg_
        INTEGER                                      :: maxint_
        INTEGER                                      :: ierr_,ibnds_(3,2),ldims_(3),ider_(3),nd_(3)
        INTEGER                                      :: ntot_
        INTEGER                                      :: rank_
        CHARACTER(len=1024)                          :: serr_
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

  SUBROUTINE GPSplineInt_ctor(this,rank,nd,ibnds,xbnds,maxpart,gpcomm)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  Main explicit constructor
!  ARGUMENTS:
!    this    : 'this' class instance
!    rank    : rank of solution (1, 2, or 3)
!    nd      : global grid size in each direction
!    ibnds   : starting and ending indices for each direction for this MPI task.
!              Integer array of size (rank,2).
!    xbnds   : starting and stopping coordinates in each direction for this
!              task. Real array of size (rank,2)
!    maxpart : max no. interpolation points/Lag. particles
!    gpcomm  : GHOST particle communicator object
!-----------------------------------------------------------------
    IMPLICIT NONE
    CLASS(GPSplineInt)                         :: this
    TYPE(GPartComm), TARGET                    :: gpcomm
    INTEGER      ,INTENT(IN)                   :: rank,maxpart
    INTEGER      ,INTENT(IN),DIMENSION  (rank) :: nd
    INTEGER      ,INTENT(IN),DIMENSION(rank,2) :: ibnds
    INTEGER                                    :: j,k
    REAL(KIND=GP),INTENT(IN),DIMENSION(rank,2) :: xbnds


    
    this%gpcomm_ => gpcomm
    this%maxint_  = maxpart
    this%rank_    = rank
    this%ider_    = 0
    this%ldims_   = 0
    this%ntot_    = 1
    DO j = 1, this%rank_
      DO k = 1,2
        this%ibnds_(j,k)  = ibnds(j,k)
        this%xbnds_(j,k)  = xbnds(j,k)
      ENDDO
      this%ldims_(j)  = ibnds(j,2) - ibnds(j,1) + 1
      this%nd_   (j)  = nd   (j)
      this%ntot_ = this%ntot_*this%ldims_(j)
    ENDDO

    IF ( this%rank_.LT.2 .OR. this%rank_.GT.3 ) THEN
      WRITE(*,*)'GPSplineInt::ctor: Invalid rank'
      STOP
    ENDIF

    IF ( this%rank_.EQ.2 .AND. this%ldims_(2).NE.(this%ibnds_(2,2)-this%ibnds_(2,1)+1) ) THEN
      WRITE(*,*) 'GPSplineInt::ctor: Inconsistent 2-indices: rank=2'
      STOP
    ENDIF

    IF ( this%rank_ .EQ. 3 .AND. this%ldims_(3) .NE.  (this%ibnds_(3,2)-this%ibnds_(3,1)+1) ) THEN
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
    USE fprecision
    IMPLICIT NONE
    CLASS(GPSplineInt)                        :: this
    INTEGER      ,INTENT   (IN)               :: np
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(np) :: fp

!

  END SUBROUTINE GPSplineInt_Interp2D
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
    USE fprecision
    IMPLICIT NONE
    CLASS(GPSplineInt)                        :: this
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
    sixth  = 1.0/6.0
    four   = 4.0
    three  = 3.0
    six    = 6.0
    halfm  = -half
    threeh = 3.0/2.0
    two    = 2.0
!
    IF ( this%ider_(1).EQ.0 ) THEN
    xsm=1.0
    DO lag=1,np
      xx = this%xrk_(lag)
      xxm = (1.0-xx)
       this%wrkl_(1,lag) = sixth*xxm*xxm*xxm
       this%wrkl_(2,lag) = sixth*(four+xx*xx*(three*xx-six))
       this%wrkl_(3,lag) = sixth*(four+xxm*xxm*(three*xxm-six))
    ENDDO
    ENDIF
!
    IF ( this%ider_(1).EQ.1 ) THEN
    xsm=0.0
    DO lag=1,np
      xx = this%xrk_(lag)
      xxm = (1.0-xx)
       this%wrkl_(1,lag) = halfm*xxm*xxm*this%dxi_(1)
       this%wrkl_(2,lag) = xx*(threeh*xx-two)*this%dxi_(1)
       this%wrkl_(3,lag) = -xxm*(threeh*xxm-two)*this%dxi_(1)
    ENDDO
    ENDIF
!
    IF ( this%ider_(2).EQ.0 ) THEN
    ysm=1.0
    DO lag=1,np
      yy = this%yrk_(lag)
      yym = (1.0-yy)
       this%wrkl_(4,lag) = sixth*yym*yym*yym
       this%wrkl_(5,lag) = sixth*(four+yy*yy*(three*yy-six))
       this%wrkl_(6,lag) = sixth*(four+yym*yym*(three*yym-six))
    ENDDO
    ENDIF
!
    IF( this%ider_(2).EQ.1 ) THEN
    ysm=0.0
    DO lag=1,np
      yy = this%yrk_(lag)
      yym = (1.0-yy)
       this%wrkl_(4,lag) = halfm*yym*yym*this%dxi_(2)
       this%wrkl_(5,lag) = yy*(threeh*yy-two)*this%dxi_(2)
       this%wrkl_(6,lag) = -yym*(threeh*yym-two)*this%dxi_(2)
    ENDDO
    ENDIF
!
    IF ( this%ider_(3).EQ.0 ) THEN
    zsm=1.0
    DO lag=1,np
      zz = this%zrk_(lag)
      zzm = (1.0-zz)
       this%wrkl_(7,lag) = sixth*zzm*zzm*zzm
       this%wrkl_(8,lag) = sixth*(four+zz*zz*(three*zz-six))
       this%wrkl_(9,lag) = sixth*(four+zzm*zzm*(three*zzm-six))
    ENDDO
    ENDIF
!
    IF( this%ider_(3).EQ.1 ) THEN
    zsm=0.0
    DO lag=1,np
      zz = this%zrk_(lag)
      zzm = (1.0-zz)
       this%wrkl_(7,lag) = halfm*zzm*zzm*this%dxi_(3)
       this%wrkl_(8,lag) = zz*(threeh*zz-two)*this%dxi_(3)
       this%wrkl_(9,lag) = -zzm*(threeh*zzm-two)*this%dxi_(3)
    ENDDO
    ENDIF
!
    IF(this%ider_(1).eq.0) xsm = 1.
    IF(this%ider_(1).eq.1) xsm = 0.
    IF(this%ider_(1).eq.0) ysm = 1.
    IF(this%ider_(1).eq.1) ysm = 0.
    IF(this%ider_(1).eq.0) ysm = 1.
    IF(this%ider_(1).eq.1) ysm = 0.

!
    DO lag=1,np
      xx1 = this%wrkl_(1,lag)
      xx2 = this%wrkl_(2,lag)
      xx3 = this%wrkl_(3,lag)
      xx4 = xsm - xx1 -xx2 - xx3
      yy1 = this%wrkl_(4,lag)
      yy2 = this%wrkl_(5,lag)
      yy3 = this%wrkl_(6,lag)
      yy4 = ysm - yy1 -yy2 - yy3
      zz1 = this%wrkl_(7,lag)
      zz2 = this%wrkl_(8,lag)
      zz3 = this%wrkl_(9,lag)
      zz4 = zsm - zz1 -zz2 - zz3

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
    USE fprecision

    IMPLICIT NONE
    CLASS(GPSplineInt)                     :: this
    REAL(KIND=GP),INTENT(IN),DIMENSION(*) :: xp,yp
    INTEGER      ,INTENT(IN)              :: np
    INTEGER                               :: j,nx,ny
    LOGICAL                               :: bok

    ! Compute interval-normalized positions and
    ! indices into control point array in x, y directions:
    nx = this%ldims_(1)
    ny = this%ldims_(2)
 
    ! x-coords:
    DO j = 1, np
       this%ilg_(1,j) = xp(j)*this%dxi_(1)
       this%xrk_  (j) = xp(j)*this%dxi_(1) - real(this%ilg_(1,j),kind=GP)
       this%ilg_(2,j) = this%ilg_(1,j) + 1
       this%ilg_(3,j) = mod(this%ilg_(2,j),nx) + 1
       this%ilg_(4,j) = mod(this%ilg_(3,j),nx) + 1
       this%ilg_(1,j) = mod(nx+this%ilg_(2,j)-2,nx) + 1
    ENDDO
      
    ! y-coords:
    ! Note that the yp are on the _global_ grid, so
    ! we have to restrict to this MPI task's 'slab' 
    ! before computing local indices and spline coordinates....
    ! 
    ! First, do a check:
    bok = .true.
    DO j = 1, np
      bok = bok .AND. (yp(j).GE.this%xbnds_(2,1).AND.yp(j).LE.this%xbnds_(2,2))
    ENDDO
    IF ( .NOT. bok ) THEN
      WRITE(*,*) 'GPSplineInt::PartUpdate2D: Invalid particle z-range'
      STOP
    ENDIF
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
    USE fprecision

    IMPLICIT NONE
    CLASS(GPSplineInt)                     :: this
    REAL(KIND=GP),INTENT(IN),DIMENSION(*) :: xp,yp,zp
    INTEGER      ,INTENT(IN)              :: np
    INTEGER                               :: j,nx,ny,nz
    LOGICAL                               :: bok

    ! Compute interval-normalized positions and
    ! indices into control point array in x, y directions:
    nx = this%ldims_(1)
    ny = this%ldims_(2)
    nz = this%ldims_(3)
 
    ! x-coords:
    DO j = 1, np
      this%ilg_(1,j) = xp(j)*this%dxi_(1)
      this%xrk_  (j) = xp(j)*this%dxi_(1) - real(this%ilg_(1,j),kind=GP)
      this%ilg_(2,j) = this%ilg_(1,j) + 1
      this%ilg_(3,j) = mod(this%ilg_(2,j),nx) + 1
      this%ilg_(4,j) = mod(this%ilg_(3,j),nx) + 1
      this%ilg_(1,j) = mod(nx+this%ilg_(2,j)-2,nx) + 1
    ENDDO
      
    ! y-coords:
    DO j = 1, np
      this%jlg_(1,j) = yp(j)*this%dxi_(2)
      this%yrk_  (j) = yp(j)*this%dxi_(2) - real(this%jlg_(1,j),kind=GP)
      this%jlg_(2,j) = this%jlg_(1,j) + 1
      this%jlg_(3,j) = mod(this%jlg_(2,j),ny) + 1
      this%jlg_(4,j) = mod(this%jlg_(3,j),ny) + 1
      this%jlg_(1,j) = mod(ny+this%jlg_(2,j)-2,ny) + 1
    ENDDO

    ! z-coords:
    ! Note that the zp are on the _global_ grid, so
    ! we have to restrict to this M{PI task's slab 
    ! before computing local indices and spline coordinates....
    ! 
    ! First, do a check:
    bok = .true.
    DO j = 1, np
      bok = bok .AND. (zp(j).GE.this%xbnds_(3,1).AND.zp(j).LE.this%xbnds_(3,2))
    ENDDO
    IF ( .NOT. bok ) THEN
      WRITE(*,*) 'GPSplineInt::PartUpdate3D: Invalid particle z-range'
      STOP
    ENDIF
    DO j = 1, np
      this%klg_(1,j) = (zp(j)-this%xbnds_(3,1))*this%dxi_(3)
      this%zrk_  (j) = (zp(j)-this%xbnds_(3,1))*this%dxi_(3) &
                     - real(this%klg_(1,j),kind=GP)
      this%klg_(2,j) = this%klg_(1,j) + 1
      this%klg_(3,j) = this%klg_(2,j) + 1
      this%klg_(4,j) = this%klg_(3,j) + 1
      this%klg_(1,j) = this%klg_(2,j) - 1
    ENDDO

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
  
    CALL GPSplineInt_DoDealloc(this)
    nzg = this%gpcomm_%GetNumGhost()

    IF ( this%rank_ .EQ. 2 ) THEN
      ALLOCATE( this%esplfld2_(this%ldims_(1),(this%ldims_(2)+2*nzg)) )
    ELSE IF ( this%rank_ .EQ. 3 ) THEN
      ALLOCATE(this%esplfld_(this%ldims_(1),this%ldims_(2),(this%ldims_(3)+2*nzg)) )
    ENDIF

    ALLOCATE(this%ax_   (this%ldims_(1)) )
    ALLOCATE(this%bx_   (this%ldims_(1)) )
    ALLOCATE(this%cx_   (this%ldims_(1)) )
    ALLOCATE(this%gamx_ (this%ldims_(1)) )
    ALLOCATE(this%xxx_  (this%ldims_(1)) )

    ALLOCATE(this%ay_   (this%ldims_(2)) )
    ALLOCATE(this%by_   (this%ldims_(2)) )
    ALLOCATE(this%cy_   (this%ldims_(2)) )
    ALLOCATE(this%gamy_ (this%ldims_(2)) )
    ALLOCATE(this%xxy_  (this%ldims_(2)) )

    ALLOCATE(this%wrkl_(9,this%maxint_))
    ALLOCATE(this%ilg_(4,this%maxint_))
    ALLOCATE(this%jlg_(4,this%maxint_))
    ALLOCATE(this%klg_(4,this%maxint_))
    ALLOCATE(this%xrk_(this%maxint_))
    ALLOCATE(this%yrk_(this%maxint_))
    ALLOCATE(this%zrk_(this%maxint_))

    IF ( this%rank_ .GT. 2 ) THEN
    ALLOCATE(this%az_   (this%ldims_(3)) )
    ALLOCATE(this%bz_   (this%ldims_(3)) )
    ALLOCATE(this%cz_   (this%ldims_(3)) )
    ALLOCATE(this%gamz_ (this%ldims_(3)) )
    ALLOCATE(this%xxz_  (this%ldims_(3)) )
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


    IF ( ALLOCATED (this%esplfld_  ) ) DEALLOCATE (this%esplfld_)
    IF ( ALLOCATED (this%esplfld2_ ) ) DEALLOCATE(this%esplfld2_)

    IF ( ALLOCATED     (this%ax_   ) ) DEALLOCATE      (this%ax_)
    IF ( ALLOCATED     (this%bx_   ) ) DEALLOCATE      (this%bx_)
    IF ( ALLOCATED     (this%cx_   ) ) DEALLOCATE      (this%cx_)
    IF ( ALLOCATED     (this%px_   ) ) DEALLOCATE      (this%px_)
    IF ( ALLOCATED     (this%gamx_ ) ) DEALLOCATE    (this%gamx_)
    IF ( ALLOCATED     (this%xxx_  ) ) DEALLOCATE     (this%xxx_)

    IF ( ALLOCATED     (this%ay_   ) ) DEALLOCATE      (this%ay_)
    IF ( ALLOCATED     (this%by_   ) ) DEALLOCATE      (this%by_)
    IF ( ALLOCATED     (this%cy_   ) ) DEALLOCATE      (this%cy_)
    IF ( ALLOCATED     (this%py_   ) ) DEALLOCATE      (this%py_)
    IF ( ALLOCATED     (this%gamy_ ) ) DEALLOCATE    (this%gamy_)
    IF ( ALLOCATED     (this%xxy_  ) ) DEALLOCATE     (this%xxy_)

    IF ( ALLOCATED     (this%az_   ) ) DEALLOCATE      (this%az_)
    IF ( ALLOCATED     (this%bz_   ) ) DEALLOCATE      (this%bz_)
    IF ( ALLOCATED     (this%cz_   ) ) DEALLOCATE      (this%cz_)
    IF ( ALLOCATED     (this%pz_   ) ) DEALLOCATE      (this%pz_)
    IF ( ALLOCATED     (this%gamz_ ) ) DEALLOCATE    (this%gamz_)
    IF ( ALLOCATED     (this%xxz_  ) ) DEALLOCATE     (this%xxz_)

    IF ( ALLOCATED      (this%wrkl_) ) DEALLOCATE    (this%wrkl_)
    IF ( ALLOCATED       (this%ilg_) ) DEALLOCATE     (this%ilg_)
    IF ( ALLOCATED       (this%jlg_) ) DEALLOCATE     (this%jlg_)
    IF ( ALLOCATED       (this%klg_) ) DEALLOCATE     (this%klg_)
    IF ( ALLOCATED       (this%xrk_) ) DEALLOCATE     (this%xrk_)
    IF ( ALLOCATED       (this%yrk_) ) DEALLOCATE     (this%yrk_)
    IF ( ALLOCATED       (this%zrk_) ) DEALLOCATE     (this%zrk_)

  END SUBROUTINE GPSplineInt_DoDealloc
!-----------------------------------------------------------------
!-----------------------------------------------------------------


  SUBROUTINE GPSplineInt_Init(this)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : Init
!  DESCRIPTION: Does initialzation of object. Called by constructor
!  ARGUMENTS  : 
!    this     : 'this' class instance
!-----------------------------------------------------------------
    USE fprecision

    IMPLICIT NONE
    CLASS(GPSplineInt)        :: this
    REAL(KIND=GP)            :: beta,sixth,twotrd
    INTEGER                  :: j


    CALL GPSplineInt_DoAlloc(this)

    ! Compute interval widths:
    DO j = 1, this%rank_
      this%dxi_(j) = real(this%ldims_(j)-1,kind=GP)/ ( this%xbnds_(j,2) - this%xbnds_(j,1) )
    ENDDO

    CALL GPSplineInt_MatInvQ(this,this%ldims_(1),this%ax_,this%bx_,this%cx_,&
         this%px_,this%gamx_,this%betx_,this%xxx_,this%zetax_)
    CALL GPSplineInt_MatInvQ(this,this%ldims_(2),this%ay_,this%by_,this%cy_,&
    this%py_,this%gamy_,this%bety_,this%xxy_,this%zetay_)
    IF ( this%rank_ .GT. 2 ) THEN
    CALL GPSplineInt_MatInvQ(this,this%ldims_(3),this%az_,this%bz_,this%cz_,&
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
    USE fprecision
    USE commtypes

    IMPLICIT NONE
    CLASS(GPSplineInt)                       :: this
    INTEGER      ,INTENT   (IN)              :: n
    INTEGER                                  :: i
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(n) :: a,b,c,p,gam,bet,xx
    REAL(KIND=GP),INTENT(INOUT)              :: zeta
    REAL(KIND=GP)                            :: beta,sixth,twotrd


!  Setup the arrays for the inversion:
        sixth = 1.0/6.0
        twotrd=2.0/3.0
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
    USE fprecision

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
    USE fprecision

    IMPLICIT NONE
    CLASS(GPSplineInt)                                :: this
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(this%ntot_) :: field
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(this%ntot_) :: tmp1,tmp2
    INTEGER                                           :: i,j,k
    INTEGER                                           :: jm,km
    INTEGER                                           :: nx,nxy,ny,nz,idims(3)

    nx = this%ldims_(1)
    ny = this%ldims_(2)
    nz = this%ldims_(3)
    nxy = nx*ny

!   tmp2 = reshape(tmp2,(/nx,ny,nz/))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!! x-field computation !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DO k=1,nz
      km = j-1
      DO j=1,ny
        jm = j-1
        tmp2 (1+jm*nx+km*nxy) = field(1+jm*nx+km*nxy)*this%betx_(1)
        tmp2(nx+jm*nx+km*nxy) = field(nx+jm*nx+km*nxy)
      ENDDO
    ENDDO
!
    DO i=2,nx-2
      DO k=1,nz
        km = k-1
        DO j=1,ny
           jm = j-1
           tmp2(i+jm*nx+km*nxy) =  &
           ( field(i+jm*nx+km*nxy) - this%ax_(i)*tmp2(i-1+jm*nx+km*nxy) )*this%betx_(i)  
           tmp2(nx+jm*nx+km*nxy) = tmp2(nx+jm*nx+km*nxy) - this%xxx_(i-1)*tmp2(i-1+jm*nx+km*nxy)
        ENDDO
      ENDDO
    ENDDO
!
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
    DO  i=nx-2,1,-1
      DO j=1,ny
        DO k=1,nz
        tmp2(i+jm*nx+km*nxy) = tmp2(i+jm*nx+km*nxy) &
                      - this%gamx_(i+1)*tmp2(i+1+jm*nx+km*nxy) - this%px_(i)*tmp2(nx+jm*nx+km*nxy)
        ENDDO
      ENDDO
    ENDDO

!  Copy splfld -> field :
    DO k=1,nz
      km = k-1
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
    DO k=1,nz
      km = k-1
      DO i=1,nx
        tmp2(i+   km*nxy) = field(i+km*nxy)*this%bety_(1)
        tmp2(i+(ny-1)*nx+km*nxy) = field(i+(ny-1)*nx+km*nxy)
      ENDDO
    ENDDO
!
    DO j=2,ny-2
      jm = j-1
      DO k=1,nz
         km = k-1
         DO i=1,nx
           tmp2    (i+jm*nx+km*nxy) = &
           ( field(i+jm*nx+km*nxy) - this%ay_(j)*tmp2(i+(jm-1)*nx+km*nxy) )*this%bety_(j)
           tmp2   (i+(ny-1)*nx+km*nxy) = tmp2(i+(ny-1)*jm+km*nxy) - this%xxy_(j-1)*tmp2(i+(jm-1)*nx+km*nxy)
         ENDDO
       ENDDO
    ENDDO
!  ** n-1 **
    DO k=1,nz
      km = k-1
      DO i=1,nx
        tmp2    (i+(ny-2)*nx+km*nxy) = &
        (field(i+(ny-1)*nx+km*nxy) - this%ay_(ny-1)*tmp2(i+(ny-3)*nx+km*nxy))*this%bety_(ny-1)
        tmp2      (i+(ny-1)*nx+km*nxy) = tmp2(i+(ny-1)*nx+km*nxy) - this%xxy_(ny-2)*tmp2(i+(ny-2)*nx+km*nxy)
!  ** n  **
        tmp2(i+(ny-1)*nx+km*nxy) = (tmp2(i+(ny-1)*nx+km*nxy) - tmp2(i+(ny-2)*nx+km*nxy)*this%zetay_) &
                       * this%bety_(ny)
!  Backsubstitution phase :
      tmp2(i+(ny-2)*nx+km*nxy) = tmp2(i+(ny-2)*nx+km*nxy) - this%gamy_(ny)*tmp2(i+(ny-1)*nx+km*nxy)
      ENDDO
    ENDDO
!
    DO  j=ny-2,1,-1
      DO i=1,nx
        DO k=1,nz
          tmp2(i+jm*nx+km*nxy) = tmp2(i+jm*nx+km*nxy)  &
                         - this%gamy_(j+1)*tmp2(i+(jm+1)*nx+km*nxy) - this%py_(j)*tmp2(i+(ny-1)*nx+km*nxy)
        ENDDO
      ENDDO
    ENDDO
!
!  Copy splfld -> field :
!!  DO k=1,nz
!!    km = k-1
!!    DO j=1,ny
!!      jm = j-1
!!      DO i=1,nx
!!        field(i+jm*nx+km*nxy) = tmp2(i,j,k)
!!      ENDDO
!!    ENDDO
!!  ENDDO

! Note tmp2 now contains the contributions for full tensor product
! field from x, y only....

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!! z-field computation !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Must transpose data so that it is z-y complete, in order to
! carry out computation in z:
    CALL this%gpcomm_%GTranspose(field,tmp2,this%ldims_,3,tmp1)
    nx = this%nd_(3)
    ny = this%nd_(2)
    nz = this%ldims_(3)
    nxy = nx*ny
!
!  Compute the k-equation (coefficient in z) :
!  Note that the we work on the transposed system here...;
!  
!  tmp2  <== spline coeff. field;
!  field <== 'field', semi-'tensored' with spl. coeffs
    DO k=1,nz
      km = k-1
      DO j=1,ny
        jm = j-1
        tmp2 (1+jm*nx+km*nxy) = field(1+jm*nx+km*nxy)*this%betx_(1)
        tmp2(nx+jm*nx+km*nxy) = field(nx+jm*nx+km*nxy)
      ENDDO
    ENDDO
!
    DO i=2,nx-2
      DO k=1,nz
        km = k-1
        DO j=1,ny
          jm = j-1
          tmp2 (i+jm*nx+km*nxy) =  &
          (field(i+jm*nx+km*nxy) - this%ax_(i)*tmp2(i-1+jm*nx+km*nxy) )*this%betx_(i) 
          tmp2(nx+jm*nx+km*nxy) = tmp2(nx+jm*nx+km*nxy) - this%xxx_(i-1)*tmp2(i-1+jm*nx+km*nxy)
        ENDDO
      ENDDO
    ENDDO
!
    DO k=1,nz
      km = k-1
      DO j=1,ny
        jm = j-1
!  ** n-1 **
        tmp2 (nx-1+jm*nx+km*nxy) = &
        (field(nx-1+jm*nx+km*nxy) - this%ax_(nx-1)*tmp2(nx-2+jm*nx+km*nxy))*this%betx_(nx-1)
        tmp2 (nx+jm*nx+km*nxy) = tmp2(nx+jm*nx+km*nxy) - this%xxx_(nx-2)*tmp2(nx-1+jm*nx+km*nxy)
!  ** n  **
        tmp2 (nx+jm*nx+km*nxy) = (tmp2(nx+jm*nx+km*nxy) - tmp2(nx-1+jm*nx+km*nxy)*this%zetax_) &
                                * this%betx_(nx)
!  Backsubstitution phase :
        tmp2(nx-1+jm*nx+km*nxy) = tmp2(nx-1+jm*nx+km*nxy) - this%gamx_(nx)*tmp2(nx+jm*nx+km*nxy)
      ENDDO
    ENDDO
!
    DO  i=nx-2,1,-1
      DO j=1,ny
        jm = j-1
        DO k=1,nz
          km = k-1
          tmp2(i+jm*nx+km*nxy) = tmp2(i+jm*nx+km*nxy) &
          - this%gamx_(i+1)*tmp2(i+1+jm*nx+km*nxy) - this%px_(i)*tmp2(nx+jm*nx+km*nxy)
        ENDDO
      ENDDO
    ENDDO
!
!  Transpose back to standard coords:
    idims(1) = this%nd_   (1)
    idims(2) = this%nd_   (2)
    idims(3) = this%ldims_(3)
    CALL this%gpcomm_%GTranspose(field,tmp2,idims,3,tmp1)
 
!! DO k=1,nz
!!   DO j=1,ny
!!     DO i=1,nx
!!         field(i+j*nx+k*nxy) = this%splfld_(i+j*nx+k*nxy)
!!     ENDDO
!!   ENDDO
!! ENDDO
!
! 
! Spline coeff field is now slab-decomposed to be xy complete, and 
! distributed in z. To do interpolation, we must set up 
! spline coeffs in array extended in z by 2 control points
! so that spline interpolation can be done:
     CALL this%gpcomm_%SlabDataExchangeSF(this%esplfld_,field)

   RETURN

  END SUBROUTINE GPSplineInt_CompSpline3D
!-----------------------------------------------------------------
!-----------------------------------------------------------------


END MODULE class_GPSplineInt
