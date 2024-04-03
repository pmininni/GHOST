!=================================================================
! GHOST GPICSplineInt interpolation and deposition class
! The splines used are derived from 'The Hybrid Multiscale
! Simulation Technology: An Introduction to Astrophysical and
! Laboratory Plasmas' 2010, by A. Lipatov.
! The structure is derived from the GPSplineInt interpolation
! class by A. Pumir and D. Rosenberg and retains the basic
! structure for interpolation:
!
!   CALL this%intop_%PartUpdate3D(this%px_,this%py_,this%pz_,this%nparts_)
!   CALL this%intop_%CompSpline3D(evar)
!   CALL this%intop_%DoInterp3D(lag,nl)
!
! The first updates the interpolation points; the second computes the
! control point function, which, here, is _not global_ but _local_
! (thus, global continuity is not ensured); the third ccarries out
! the actual interpolation to the interpolation points.
! A similar structure is used for deposition of particle quantities
! into the grid:
!
!   CALL this%intop_%PartUpdate(this%px_,this%py_,this%pz_,this%nparts_)
!   CALL this%intop_%CompSpline(evar)
!   CALL this%intop_%DoDeposit(lag,nl,field)
!
! Note that this class assumes a regularly spaced grid.
!
! 2023 F. Pugliese (DF, FCEN, UBA)
!=================================================================
MODULE class_GPICSplineInt
      USE mpivars
      USE fprecision
      USE class_GPartComm
      USE gtimer
      IMPLICIT NONE

      PRIVATE
      TYPE, PUBLIC :: GPICSplineInt
        PRIVATE
        ! Member data:
        REAL(KIND=GP),ALLOCATABLE,DIMENSION(:,:,:)   :: esplfld_
        DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:):: esplflddp_
        REAL(KIND=GP),ALLOCATABLE,DIMENSION  (:,:)   :: esplfld2_
        REAL(KIND=GP),ALLOCATABLE,DIMENSION    (:)   :: xrk_,yrk_,zrk_
        REAL(KIND=GP)                                :: dxi_(3),xbnds_(3,2)
        TYPE(GFieldComm),POINTER                     :: gfcomm_
        INTEGER      ,ALLOCATABLE,DIMENSION  (:,:)   :: ilg_,jlg_,klg_
        INTEGER                                      :: maxint_,splord_
        INTEGER                                      :: ierr_,ider_(3),nd_(3)
        INTEGER                                      :: ibnds_(3,2),obnds_(3,2)
        INTEGER                                      :: ldims_(3),odims_(3)
        INTEGER                                      :: hdataex_,htransp_
        INTEGER                                      :: ntot_,ttot_
        INTEGER                                      :: rank_
        CHARACTER(len=1024)                          :: serr_
      CONTAINS
        ! Public methods:
        PROCEDURE,PUBLIC :: GPICSplineInt_ctor
        FINAL            :: GPICSplineInt_dtor

        ! ...Aliases:
        PROCEDURE,PUBLIC :: Init         => GPICSplineInt_Init
        PROCEDURE,PUBLIC :: DoInterp     => GPICSplineInt_Interp
        PROCEDURE,PUBLIC :: DoDeposit    => GPICSplineInt_Deposit
        PROCEDURE,PUBLIC :: PartUpdate   => GPICSplineInt_PartUpdate
        PROCEDURE,PUBLIC :: CompSpline   => GPICSplineInt_CompSpline
        PROCEDURE,PUBLIC :: ResizeArrays => GPICSplineInt_ResizeArrays

      END TYPE GPICSplineInt

      ! Private methods:
      PRIVATE :: GPICSplineInt_Init         , GPICSplineInt_CompSpline
      PRIVATE :: GPICSplineInt_PartUpdate   , GPICSplineInt_ResizeArrays
      PRIVATE :: GPICSplineInt_Deposit3D    , GPICSplineInt_Interp3D
      PRIVATE :: GPICSplineInt_Deposit2D    , GPICSplineInt_Interp2D
      PRIVATE :: GPICSplineInt_Deposit1D    , GPICSplineInt_Interp1D
      PRIVATE :: GPICSplineInt_DoAlloc      , GPICSplineInt_DoDealloc

! Methods:
  CONTAINS

  SUBROUTINE GPICSplineInt_ctor(this,rank,nd,ibnds,xbnds,obnds,splord,nzghost, &
                                maxpart,gfcomm,hdataex,htransp)
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
!    splord  : order used for interpolation and deposition
!    nzghost : no. z-ghost zones required for interpolation
!    maxpart : max no. interpolation points/Lag. particles
!    gfcomm  : GHOST particle and field communicator object
!    hdataex : handle to data exchange. Must be valid.
!    htransp : handle to timer for transpose. Must be valid.
!-----------------------------------------------------------------
    IMPLICIT NONE
    CLASS(GPICSplineInt)                         :: this
    TYPE(GFieldComm),TARGET                      :: gfcomm
    INTEGER        ,INTENT(IN)                   :: hdataex,htransp,maxpart,nzghost,rank,splord
    INTEGER        ,INTENT(IN),DIMENSION(rank)   :: nd
    INTEGER        ,INTENT(IN),DIMENSION(rank,2) :: ibnds,obnds
    INTEGER                                      :: j,k
    REAL(KIND=GP)  ,INTENT(IN),DIMENSION(rank,2) :: xbnds

    this%gfcomm_ => gfcomm
    this%maxint_  = maxpart
    this%rank_    = rank
    this%ider_    = 0
    this%ldims_   = 0
    this%odims_   = 0
    this%ntot_    = 1
    this%ttot_    = 1
    this%splord_  = splord
    j = GTValidHandle(htransp)
    IF ( j.NE.GTERR_GOOD_HANDLE ) THEN
      WRITE(*,*) 'GPICSplineInt_ctor: invalid transpose timer handle: ',j
      STOP
    ENDIF
    this%htransp_  = htransp
    j = GTValidHandle(hdataex)
    IF ( j.NE.GTERR_GOOD_HANDLE ) THEN
      WRITE(*,*) 'GPICSplineInt_ctor: invalid data exch. timer handle: ',j
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
      WRITE(*,*)'GPICSplineInt::ctor: Invalid rank'
      STOP
    ENDIF

    IF ( this%rank_.EQ.2 .AND. this%ldims_(2).NE.(this%ibnds_(2,2)-this%ibnds_(2,1)+1) ) THEN
      WRITE(*,*) 'GPICSplineInt::ctor: Inconsistent 2-indices: rank=2'
      STOP
    ENDIF

    IF ( this%rank_.EQ.3 .AND. this%ldims_(3).NE.(this%ibnds_(3,2)-this%ibnds_(3,1)+1) ) THEN
      WRITE(*,*) 'GPICSplineInt::ctor: Inconsitent 3-indices; rank=3'
      STOP
    ENDIF

    CALL GPICSplineInt_Init(this)

  END SUBROUTINE GPICSplineInt_ctor
!-----------------------------------------------------------------
!-----------------------------------------------------------------


  SUBROUTINE GPICSplineInt_dtor(this)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  Main explicit destructor
!  ARGUMENTS:
!    this    : 'this' class instance
!-----------------------------------------------------------------

    IMPLICIT NONE
    TYPE(GPICSplineInt),INTENT(INOUT)        :: this

    CALL GPICSplineInt_DoDealloc(this)

  END SUBROUTINE GPICSplineInt_dtor
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GPICSplineInt_Interp(this,fp,np)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  DESCRIPTION:  Interp
!               Carry out interpolation. In addition, the normalized 
!               positions of the particles in each interval must be 
!               computed prior to entry, as must be the indices into 
!               the full grid that define the range of control points.
!  ARGUMENTS  :
!    this     : 'this' class instance
!    fp       : field interpolated at the points xp, yp, returned.
!    np       : no. interploation points (lenghts of fp,xp,yp,zp).
!-----------------------------------------------------------------
    IMPLICIT NONE
    CLASS(GPICSplineInt)                      :: this
    INTEGER      ,INTENT   (IN)               :: np
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(np) :: fp

    IF (this%ldims_(1).EQ.1) THEN
      CALL GPICSplineInt_Interp1D(this,fp,np)
    ELSE IF (this%ldims_(2).EQ.1) THEN
      CALL GPICSplineInt_Interp2D(this,fp,np)
    ELSE
      CALL GPICSplineInt_Interp3D(this,fp,np)
    END IF

  END SUBROUTINE GPICSplineInt_Interp
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GPICSplineInt_Interp1D(this,fp,np)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  DESCRIPTION:  Interp1D
!               Carry out 2d interpolation.
!  ARGUMENTS  :
!    this     : 'this' class instance
!    fp       : field interpolated at the points xp, yp, returned.
!    np       : no. interploation points (lenghts of fp,xp,yp,zp).
!-----------------------------------------------------------------
    IMPLICIT NONE
    CLASS(GPICSplineInt)                      :: this
    INTEGER      ,INTENT   (IN)               :: np
    INTEGER                                   :: lag,nx,ny,nxy,nz
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(np) :: fp
    REAL(KIND=GP)                             :: zz,zzm
    REAL(KIND=GP)                             :: zz1,zz2,zz3,zz4
    REAL(KIND=GP)                             :: two,three,four,six,sixth,eighth
!
    IF (this%splord_.eq.0) THEN
      DO lag=1, np
        fp(lag) = this%esplfld_(1,1,this%klg_(1,lag))
      ENDDO
    ELSE IF (this%splord_.eq.1) THEN
!$omp parallel do private(zz,zzm)
      DO lag=1, np
        zz  = this%zrk_(lag)
        zzm = 1.0_GP - zz
        fp(lag) = this%esplfld_(1,1,this%klg_(1,lag))*zzm     &
                + this%esplfld_(1,1,this%klg_(2,lag))*zz
      ENDDO
    ELSE IF (this%splord_.eq.2) THEN
      two    = 2.0_GP
      eighth = 1.0_GP/8.0_GP
!$omp parallel do private(zz,zz1,zz2,zz3)
      DO lag=1, np
        zz  = this%zrk_(lag)
        zz1 = eighth*(1-two*zz)*(1-two*zz)
        zz3 = eighth*(1+two*zz)*(1+two*zz)
        zz2 = 1.0_GP - zz1 - zz3
        fp(lag) = this%esplfld_(1,1,this%klg_(1,lag))*zz1    &
                + this%esplfld_(1,1,this%klg_(2,lag))*zz2    &
                + this%esplfld_(1,1,this%klg_(3,lag))*zz3
      ENDDO
    ELSE IF (this%splord_.eq.3) THEN
      sixth  = 1.0_GP/6.0_GP
      six    = 6.0_GP
      four   = 4.0_GP
      three  = 3.0_GP
!$omp parallel do private(zz,zzm,zz1,zz2,zz3,zz4)
      DO lag=1,np
        zz = this%zrk_(lag)
        zzm = (1.0_GP-zz)
        zz1 = sixth*zzm*zzm*zzm
        zz2 = sixth*(four+zz *zz *(three*zz -six))
        zz3 = sixth*(four+zzm*zzm*(three*zzm-six))
        zz4 = 1.0_GP - zz1 - zz2 - zz3
        fp(lag) = this%esplfld_(1,1,this%klg_(1,lag))*zz1 &
                + this%esplfld_(1,1,this%klg_(2,lag))*zz2 &
                + this%esplfld_(1,1,this%klg_(3,lag))*zz3 &
                + this%esplfld_(1,1,this%klg_(4,lag))*zz4
      ENDDO
    ENDIF

  END SUBROUTINE GPICSplineInt_Interp1D
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GPICSplineInt_Interp2D(this,fp,np)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  DESCRIPTION:  Interp2D
!               Carry out 2d interpolation.
!  ARGUMENTS  :
!    this     : 'this' class instance
!    fp       : field interpolated at the points xp, yp, returned.
!    np       : no. interploation points (lenghts of fp,xp,yp,zp).
!-----------------------------------------------------------------
    IMPLICIT NONE
    CLASS(GPICSplineInt)                      :: this
    INTEGER      ,INTENT   (IN)               :: np
    INTEGER                                   :: lag,nx,ny,nxy,nz
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(np) :: fp
    REAL(KIND=GP)                             :: xx,xxm,zz,zzm
    REAL(KIND=GP)                             :: xx1,xx2,xx3,xx4
    REAL(KIND=GP)                             :: zz1,zz2,zz3,zz4
    REAL(KIND=GP)                             :: two,three,four,six,sixth,eighth
!
    IF (this%splord_.eq.0) THEN
      DO lag=1, np
        fp(lag) = this%esplfld_(this%ilg_(1,lag),1,this%klg_(1,lag))
      ENDDO
    ELSE IF (this%splord_.eq.1) THEN
!$omp parallel do private(xx,xxm,zz,zzm)
      DO lag=1, np
        xx  = this%xrk_(lag)
        xxm = 1.0_GP - xx
        zz  = this%zrk_(lag)
        zzm = 1.0_GP - zz
        fp(lag) =  &
         (this%esplfld_(this%ilg_(1,lag),1,this%klg_(1,lag))*xxm         &
        + this%esplfld_(this%ilg_(2,lag),1,this%klg_(1,lag))*xx)*zzm     &
        +(this%esplfld_(this%ilg_(1,lag),1,this%klg_(2,lag))*xxm         &
        + this%esplfld_(this%ilg_(2,lag),1,this%klg_(2,lag))*xx)*zz
      ENDDO
    ELSE IF (this%splord_.eq.2) THEN
      two    = 2.0_GP
      eighth = 1.0_GP/8.0_GP
!$omp parallel do private(xx,xx1,xx2,xx3,zz,zz1,zz2,zz3)
      DO lag=1, np
        xx  = this%xrk_(lag)
        xx1 = eighth*(1-two*xx)*(1-two*xx)
        xx3 = eighth*(1+two*xx)*(1+two*xx)
        xx2 = 1.0_GP - xx1 - xx3
        zz  = this%zrk_(lag)
        zz1 = eighth*(1-two*zz)*(1-two*zz)
        zz3 = eighth*(1+two*zz)*(1+two*zz)
        zz2 = 1.0_GP - zz1 - zz3
        fp(lag) =  &
         (this%esplfld_(this%ilg_(1,lag),1,this%klg_(1,lag))*xx1         &
        + this%esplfld_(this%ilg_(2,lag),1,this%klg_(1,lag))*xx2         &
        + this%esplfld_(this%ilg_(3,lag),1,this%klg_(1,lag))*xx3)*zz1    &
  ! k = 2
        +(this%esplfld_(this%ilg_(1,lag),1,this%klg_(2,lag))*xx1         &
        + this%esplfld_(this%ilg_(2,lag),1,this%klg_(2,lag))*xx2         &
        + this%esplfld_(this%ilg_(3,lag),1,this%klg_(2,lag))*xx3)*zz2    &
  ! k = 3
        +(this%esplfld_(this%ilg_(1,lag),1,this%klg_(3,lag))*xx1         &
        + this%esplfld_(this%ilg_(2,lag),1,this%klg_(3,lag))*xx2         &
        + this%esplfld_(this%ilg_(3,lag),1,this%klg_(3,lag))*xx3)*zz3
      ENDDO
    ELSE IF (this%splord_.eq.3) THEN
      sixth  = 1.0_GP/6.0_GP
      six    = 6.0_GP
      four   = 4.0_GP
      three  = 3.0_GP
!$omp parallel do private(xx,xxm,xx1,xx2,xx3,xx4,zz,zzm,zz1,zz2,zz3,zz4)
      DO lag=1,np
        xx = this%xrk_(lag)
        xxm = (1.0_GP-xx)
        xx1 = sixth*xxm*xxm*xxm
        xx2 = sixth*(four+xx *xx *(three*xx -six))
        xx3 = sixth*(four+xxm*xxm*(three*xxm-six))
        xx4 = 1.0_GP - xx1 - xx2 - xx3

        zz = this%zrk_(lag)
        zzm = (1.0_GP-zz)
        zz1 = sixth*zzm*zzm*zzm
        zz2 = sixth*(four+zz *zz *(three*zz -six))
        zz3 = sixth*(four+zzm*zzm*(three*zzm-six))
        zz4 = 1.0_GP - zz1 - zz2 - zz3
  !write(*,*)'part=',lag,' xx1=',xx1,' xx2=',xx2,' xx3=',xx3,' xx4=',xx4
  !write(*,*)'part=',lag,' yy1=',yy1,' yy2=',yy2,' yy3=',yy3,' ty4=',yy4
  !write(*,*)'part=',lag,' zz1=',zz1,' zz2=',zz2,' zz3=',zz3,' zz4=',zz4
        fp(lag) =  &
         (this%esplfld_(this%ilg_(1,lag),1,this%klg_(1,lag))*xx1         &
        + this%esplfld_(this%ilg_(2,lag),1,this%klg_(1,lag))*xx2         &
        + this%esplfld_(this%ilg_(3,lag),1,this%klg_(1,lag))*xx3         &
        + this%esplfld_(this%ilg_(4,lag),1,this%klg_(1,lag))*xx4)*zz1    &
  !  k = 2
        +(this%esplfld_(this%ilg_(1,lag),1,this%klg_(2,lag))*xx1         &
        + this%esplfld_(this%ilg_(2,lag),1,this%klg_(2,lag))*xx2         &
        + this%esplfld_(this%ilg_(3,lag),1,this%klg_(2,lag))*xx3         &
        + this%esplfld_(this%ilg_(4,lag),1,this%klg_(2,lag))*xx4)*zz2    &
  !  k = 3
        +(this%esplfld_(this%ilg_(1,lag),1,this%klg_(3,lag))*xx1         &
        + this%esplfld_(this%ilg_(2,lag),1,this%klg_(3,lag))*xx2         &
        + this%esplfld_(this%ilg_(3,lag),1,this%klg_(3,lag))*xx3         &
        + this%esplfld_(this%ilg_(4,lag),1,this%klg_(3,lag))*xx4)*zz3 &
  !  k = 4
        +(this%esplfld_(this%ilg_(1,lag),1,this%klg_(4,lag))*xx1         &
        + this%esplfld_(this%ilg_(2,lag),1,this%klg_(4,lag))*xx2         &
        + this%esplfld_(this%ilg_(3,lag),1,this%klg_(4,lag))*xx3         &
        + this%esplfld_(this%ilg_(4,lag),1,this%klg_(4,lag))*xx4)*zz4
      ENDDO
    ENDIF

  END SUBROUTINE GPICSplineInt_Interp2D
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GPICSplineInt_Interp3D(this,fp,np)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  DESCRIPTION:  Interp3D
!               Carry out 3d interpolation. ¿¿Spline control points
!               must have been computed on field to be interpolated
!               by call to GPICSplineInt_CompSpline, before entering.??
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
    CLASS(GPICSplineInt)                      :: this
    INTEGER      ,INTENT   (IN)               :: np
    INTEGER                                   :: lag,nx,ny,nxy,nz
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(np) :: fp
    REAL(KIND=GP)                             :: xx,xxm,yy,yym,zz,zzm
    REAL(KIND=GP)                             :: xx1,xx2,xx3,xx4
    REAL(KIND=GP)                             :: yy1,yy2,yy3,yy4
    REAL(KIND=GP)                             :: zz1,zz2,zz3,zz4
    REAL(KIND=GP)                             :: two,three,four,six,sixth,eighth
!
    IF (this%splord_.eq.0) THEN
      DO lag=1, np
        fp(lag) = this%esplfld_(this%ilg_(1,lag),this%jlg_(1,lag),this%klg_(1,lag))
      ENDDO
    ELSE IF (this%splord_.eq.1) THEN
!$omp parallel do private(xx,xxm,yy,yym,zz,zzm)
      DO lag=1, np
        xx  = this%xrk_(lag)
        xxm = 1.0_GP - xx
        yy  = this%yrk_(lag)
        yym = 1.0_GP - yy
        zz  = this%zrk_(lag)
        zzm = 1.0_GP - zz
        fp(lag) =  &
        ((this%esplfld_(this%ilg_(1,lag),this%jlg_(1,lag),this%klg_(1,lag))*xxm         &
        + this%esplfld_(this%ilg_(2,lag),this%jlg_(1,lag),this%klg_(1,lag))*xx)*yym     &
        +(this%esplfld_(this%ilg_(1,lag),this%jlg_(2,lag),this%klg_(1,lag))*xxm         &
        + this%esplfld_(this%ilg_(2,lag),this%jlg_(2,lag),this%klg_(1,lag))*xx)*yy)*zzm &
      + ((this%esplfld_(this%ilg_(1,lag),this%jlg_(1,lag),this%klg_(2,lag))*xxm         &
        + this%esplfld_(this%ilg_(2,lag),this%jlg_(1,lag),this%klg_(2,lag))*xx)*yym     &
        +(this%esplfld_(this%ilg_(1,lag),this%jlg_(2,lag),this%klg_(2,lag))*xxm         &
        + this%esplfld_(this%ilg_(2,lag),this%jlg_(2,lag),this%klg_(2,lag))*xx)*yy)*zz
      ENDDO
    ELSE IF (this%splord_.eq.2) THEN
      two    = 2.0_GP
      eighth = 1.0_GP/8.0_GP
!$omp parallel do private(xx,xx1,xx2,xx3,yy,yy1,yy2,yy3,zz,zz1,zz2,zz3)
      DO lag=1, np
        xx  = this%xrk_(lag)
        xx1 = eighth*(1-two*xx)*(1-two*xx)
        xx3 = eighth*(1+two*xx)*(1+two*xx)
        xx2 = 1.0_GP - xx1 - xx3
        yy  = this%yrk_(lag)
        yy1 = eighth*(1-two*yy)*(1-two*yy)
        yy3 = eighth*(1+two*yy)*(1+two*yy)
        yy2 = 1.0_GP - yy1 - yy3
        zz  = this%zrk_(lag)
        zz1 = eighth*(1-two*zz)*(1-two*zz)
        zz3 = eighth*(1+two*zz)*(1+two*zz)
        zz2 = 1.0_GP - zz1 - zz3
        fp(lag) =  &
        ((this%esplfld_(this%ilg_(1,lag),this%jlg_(1,lag),this%klg_(1,lag))*xx1         &
        + this%esplfld_(this%ilg_(2,lag),this%jlg_(1,lag),this%klg_(1,lag))*xx2         &
        + this%esplfld_(this%ilg_(3,lag),this%jlg_(1,lag),this%klg_(1,lag))*xx3)*yy1    &
        +(this%esplfld_(this%ilg_(1,lag),this%jlg_(2,lag),this%klg_(1,lag))*xx1         &
        + this%esplfld_(this%ilg_(2,lag),this%jlg_(2,lag),this%klg_(1,lag))*xx2         &
        + this%esplfld_(this%ilg_(3,lag),this%jlg_(2,lag),this%klg_(1,lag))*xx3)*yy2    &
        +(this%esplfld_(this%ilg_(1,lag),this%jlg_(3,lag),this%klg_(1,lag))*xx1         &
        + this%esplfld_(this%ilg_(2,lag),this%jlg_(3,lag),this%klg_(1,lag))*xx2         &
        + this%esplfld_(this%ilg_(3,lag),this%jlg_(3,lag),this%klg_(1,lag))*xx3)*yy3)*zz1 &
  ! k = 2
      + ((this%esplfld_(this%ilg_(1,lag),this%jlg_(1,lag),this%klg_(2,lag))*xx1         &
        + this%esplfld_(this%ilg_(2,lag),this%jlg_(1,lag),this%klg_(2,lag))*xx2         &
        + this%esplfld_(this%ilg_(3,lag),this%jlg_(1,lag),this%klg_(2,lag))*xx3)*yy1    &
        +(this%esplfld_(this%ilg_(1,lag),this%jlg_(2,lag),this%klg_(2,lag))*xx1         &
        + this%esplfld_(this%ilg_(2,lag),this%jlg_(2,lag),this%klg_(2,lag))*xx2         &
        + this%esplfld_(this%ilg_(3,lag),this%jlg_(2,lag),this%klg_(2,lag))*xx3)*yy2    &
        +(this%esplfld_(this%ilg_(1,lag),this%jlg_(3,lag),this%klg_(2,lag))*xx1         &
        + this%esplfld_(this%ilg_(2,lag),this%jlg_(3,lag),this%klg_(2,lag))*xx2         &
        + this%esplfld_(this%ilg_(3,lag),this%jlg_(3,lag),this%klg_(2,lag))*xx3)*yy3)*zz2 &
  ! k = 3
      + ((this%esplfld_(this%ilg_(1,lag),this%jlg_(1,lag),this%klg_(3,lag))*xx1         &
        + this%esplfld_(this%ilg_(2,lag),this%jlg_(1,lag),this%klg_(3,lag))*xx2         &
        + this%esplfld_(this%ilg_(3,lag),this%jlg_(1,lag),this%klg_(3,lag))*xx3)*yy1    &
        +(this%esplfld_(this%ilg_(1,lag),this%jlg_(2,lag),this%klg_(3,lag))*xx1         &
        + this%esplfld_(this%ilg_(2,lag),this%jlg_(2,lag),this%klg_(3,lag))*xx2         &
        + this%esplfld_(this%ilg_(3,lag),this%jlg_(2,lag),this%klg_(3,lag))*xx3)*yy2    &
        +(this%esplfld_(this%ilg_(1,lag),this%jlg_(3,lag),this%klg_(3,lag))*xx1         &
        + this%esplfld_(this%ilg_(2,lag),this%jlg_(3,lag),this%klg_(3,lag))*xx2         &
        + this%esplfld_(this%ilg_(3,lag),this%jlg_(3,lag),this%klg_(3,lag))*xx3)*yy3)*zz3
      ENDDO
    ELSE IF (this%splord_.eq.3) THEN
      sixth  = 1.0_GP/6.0_GP
      six    = 6.0_GP
      four   = 4.0_GP
      three  = 3.0_GP
!$omp parallel do private(xx,xxm,xx1,xx2,xx3,xx4,yy,yym,yy1,yy2,yy3,yy4,zz,zzm,zz1,zz2,zz3,zz4)
      DO lag=1,np
        xx = this%xrk_(lag)
        xxm = (1.0_GP-xx)
        xx1 = sixth*xxm*xxm*xxm
        xx2 = sixth*(four+xx *xx *(three*xx -six))
        xx3 = sixth*(four+xxm*xxm*(three*xxm-six))
        xx4 = 1.0_GP - xx1 - xx2 - xx3

        yy = this%yrk_(lag)
        yym = (1.0_GP-yy)
        yy1 = sixth*yym*yym*yym
        yy2 = sixth*(four+yy *yy *(three*yy -six))
        yy3 = sixth*(four+yym*yym*(three*yym-six))
        yy4 = 1.0_GP - yy1 - yy2 - yy3

        zz = this%zrk_(lag)
        zzm = (1.0_GP-zz)
        zz1 = sixth*zzm*zzm*zzm
        zz2 = sixth*(four+zz *zz *(three*zz -six))
        zz3 = sixth*(four+zzm*zzm*(three*zzm-six))
        zz4 = 1.0_GP - zz1 - zz2 - zz3
  !write(*,*)'part=',lag,' xx1=',xx1,' xx2=',xx2,' xx3=',xx3,' xx4=',xx4
  !write(*,*)'part=',lag,' yy1=',yy1,' yy2=',yy2,' yy3=',yy3,' ty4=',yy4
  !write(*,*)'part=',lag,' zz1=',zz1,' zz2=',zz2,' zz3=',zz3,' zz4=',zz4
        fp(lag) =  &
        ((this%esplfld_(this%ilg_(1,lag),this%jlg_(1,lag),this%klg_(1,lag))*xx1         &
        + this%esplfld_(this%ilg_(2,lag),this%jlg_(1,lag),this%klg_(1,lag))*xx2         &
        + this%esplfld_(this%ilg_(3,lag),this%jlg_(1,lag),this%klg_(1,lag))*xx3         &
        + this%esplfld_(this%ilg_(4,lag),this%jlg_(1,lag),this%klg_(1,lag))*xx4)*yy1    &
        +(this%esplfld_(this%ilg_(1,lag),this%jlg_(2,lag),this%klg_(1,lag))*xx1         &
        + this%esplfld_(this%ilg_(2,lag),this%jlg_(2,lag),this%klg_(1,lag))*xx2         &
        + this%esplfld_(this%ilg_(3,lag),this%jlg_(2,lag),this%klg_(1,lag))*xx3         &
        + this%esplfld_(this%ilg_(4,lag),this%jlg_(2,lag),this%klg_(1,lag))*xx4)*yy2    &
        +(this%esplfld_(this%ilg_(1,lag),this%jlg_(3,lag),this%klg_(1,lag))*xx1         &
        + this%esplfld_(this%ilg_(2,lag),this%jlg_(3,lag),this%klg_(1,lag))*xx2         &
        + this%esplfld_(this%ilg_(3,lag),this%jlg_(3,lag),this%klg_(1,lag))*xx3         &
        + this%esplfld_(this%ilg_(4,lag),this%jlg_(3,lag),this%klg_(1,lag))*xx4)*yy3    &
        +(this%esplfld_(this%ilg_(1,lag),this%jlg_(4,lag),this%klg_(1,lag))*xx1         &
        + this%esplfld_(this%ilg_(2,lag),this%jlg_(4,lag),this%klg_(1,lag))*xx2         &
        + this%esplfld_(this%ilg_(3,lag),this%jlg_(4,lag),this%klg_(1,lag))*xx3         &
        + this%esplfld_(this%ilg_(4,lag),this%jlg_(4,lag),this%klg_(1,lag))*xx4)*yy4)*zz1 &
  !  k = 2
      + ((this%esplfld_(this%ilg_(1,lag),this%jlg_(1,lag),this%klg_(2,lag))*xx1         &
        + this%esplfld_(this%ilg_(2,lag),this%jlg_(1,lag),this%klg_(2,lag))*xx2         &
        + this%esplfld_(this%ilg_(3,lag),this%jlg_(1,lag),this%klg_(2,lag))*xx3         &
        + this%esplfld_(this%ilg_(4,lag),this%jlg_(1,lag),this%klg_(2,lag))*xx4)*yy1    &
        +(this%esplfld_(this%ilg_(1,lag),this%jlg_(2,lag),this%klg_(2,lag))*xx1         &
        + this%esplfld_(this%ilg_(2,lag),this%jlg_(2,lag),this%klg_(2,lag))*xx2         &
        + this%esplfld_(this%ilg_(3,lag),this%jlg_(2,lag),this%klg_(2,lag))*xx3         &
        + this%esplfld_(this%ilg_(4,lag),this%jlg_(2,lag),this%klg_(2,lag))*xx4)*yy2    &
        +(this%esplfld_(this%ilg_(1,lag),this%jlg_(3,lag),this%klg_(2,lag))*xx1         &
        + this%esplfld_(this%ilg_(2,lag),this%jlg_(3,lag),this%klg_(2,lag))*xx2         &
        + this%esplfld_(this%ilg_(3,lag),this%jlg_(3,lag),this%klg_(2,lag))*xx3         &
        + this%esplfld_(this%ilg_(4,lag),this%jlg_(3,lag),this%klg_(2,lag))*xx4)*yy3    &
        +(this%esplfld_(this%ilg_(1,lag),this%jlg_(4,lag),this%klg_(2,lag))*xx1         &
        + this%esplfld_(this%ilg_(2,lag),this%jlg_(4,lag),this%klg_(2,lag))*xx2         &
        + this%esplfld_(this%ilg_(3,lag),this%jlg_(4,lag),this%klg_(2,lag))*xx3         &
        + this%esplfld_(this%ilg_(4,lag),this%jlg_(4,lag),this%klg_(2,lag))*xx4)*yy4)*zz2 &
  !  k = 3
      + ((this%esplfld_(this%ilg_(1,lag),this%jlg_(1,lag),this%klg_(3,lag))*xx1         &
        + this%esplfld_(this%ilg_(2,lag),this%jlg_(1,lag),this%klg_(3,lag))*xx2         &
        + this%esplfld_(this%ilg_(3,lag),this%jlg_(1,lag),this%klg_(3,lag))*xx3         &
        + this%esplfld_(this%ilg_(4,lag),this%jlg_(1,lag),this%klg_(3,lag))*xx4)*yy1    &
        +(this%esplfld_(this%ilg_(1,lag),this%jlg_(2,lag),this%klg_(3,lag))*xx1         &
        + this%esplfld_(this%ilg_(2,lag),this%jlg_(2,lag),this%klg_(3,lag))*xx2         &
        + this%esplfld_(this%ilg_(3,lag),this%jlg_(2,lag),this%klg_(3,lag))*xx3         &
        + this%esplfld_(this%ilg_(4,lag),this%jlg_(2,lag),this%klg_(3,lag))*xx4)*yy2    &
        +(this%esplfld_(this%ilg_(1,lag),this%jlg_(3,lag),this%klg_(3,lag))*xx1         &
        + this%esplfld_(this%ilg_(2,lag),this%jlg_(3,lag),this%klg_(3,lag))*xx2         &
        + this%esplfld_(this%ilg_(3,lag),this%jlg_(3,lag),this%klg_(3,lag))*xx3         &
        + this%esplfld_(this%ilg_(4,lag),this%jlg_(3,lag),this%klg_(3,lag))*xx4)*yy3    &
        +(this%esplfld_(this%ilg_(1,lag),this%jlg_(4,lag),this%klg_(3,lag))*xx1         &
        + this%esplfld_(this%ilg_(2,lag),this%jlg_(4,lag),this%klg_(3,lag))*xx2         &
        + this%esplfld_(this%ilg_(3,lag),this%jlg_(4,lag),this%klg_(3,lag))*xx3         &
        + this%esplfld_(this%ilg_(4,lag),this%jlg_(4,lag),this%klg_(3,lag))*xx4)*yy4)*zz3 &
  !  k = 4
      + ((this%esplfld_(this%ilg_(1,lag),this%jlg_(1,lag),this%klg_(4,lag))*xx1         &
        + this%esplfld_(this%ilg_(2,lag),this%jlg_(1,lag),this%klg_(4,lag))*xx2         &
        + this%esplfld_(this%ilg_(3,lag),this%jlg_(1,lag),this%klg_(4,lag))*xx3         &
        + this%esplfld_(this%ilg_(4,lag),this%jlg_(1,lag),this%klg_(4,lag))*xx4)*yy1    &
        +(this%esplfld_(this%ilg_(1,lag),this%jlg_(2,lag),this%klg_(4,lag))*xx1         &
        + this%esplfld_(this%ilg_(2,lag),this%jlg_(2,lag),this%klg_(4,lag))*xx2         &
        + this%esplfld_(this%ilg_(3,lag),this%jlg_(2,lag),this%klg_(4,lag))*xx3         &
        + this%esplfld_(this%ilg_(4,lag),this%jlg_(2,lag),this%klg_(4,lag))*xx4)*yy2    &
        +(this%esplfld_(this%ilg_(1,lag),this%jlg_(3,lag),this%klg_(4,lag))*xx1         &
        + this%esplfld_(this%ilg_(2,lag),this%jlg_(3,lag),this%klg_(4,lag))*xx2         &
        + this%esplfld_(this%ilg_(3,lag),this%jlg_(3,lag),this%klg_(4,lag))*xx3         &
        + this%esplfld_(this%ilg_(4,lag),this%jlg_(3,lag),this%klg_(4,lag))*xx4)*yy3    &
        +(this%esplfld_(this%ilg_(1,lag),this%jlg_(4,lag),this%klg_(4,lag))*xx1         &
        + this%esplfld_(this%ilg_(2,lag),this%jlg_(4,lag),this%klg_(4,lag))*xx2         &
        + this%esplfld_(this%ilg_(3,lag),this%jlg_(4,lag),this%klg_(4,lag))*xx3         &
        + this%esplfld_(this%ilg_(4,lag),this%jlg_(4,lag),this%klg_(4,lag))*xx4)*yy4)*zz4
      ENDDO
    ENDIF

  END SUBROUTINE GPICSplineInt_Interp3D
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GPICSplineInt_Deposit(this,prop,np,proj)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  DESCRIPTION: Deposit
!               Carry out deposition into grid points.
!               In addition, the normalized positions of the particles
!               in each interval must be computed prior to entry, as
!               must be the indices into the full grid that define the
!               range of control points. Temporary deposition are stored
!               in member data array.
!
!  ARGUMENTS  :
!    this     : 'this' class instance
!    prop     : particle property to be interpolated
!    np       : no. particles (lenght of prop).
!    proj     : returned deposition of the property into regular grid
!               points. Must have dimensions set in constructor.
!-----------------------------------------------------------------
!$  USE threads
    IMPLICIT NONE
    CLASS(GPICSplineInt)                         :: this
    INTEGER      ,INTENT   (IN)                  :: np
    REAL(KIND=GP),INTENT   (IN),DIMENSION(np)    :: prop
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(:,:,:) :: proj
 
   IF (this%ldims_(1).EQ.1) THEN
      CALL GPICSplineInt_Deposit1D(this,prop,np,proj)
    ELSE IF (this%ldims_(2).EQ.1) THEN
      CALL GPICSplineInt_Deposit2D(this,prop,np,proj)
    ELSE
      CALL GPICSplineInt_Deposit3D(this,prop,np,proj)
    END IF

  END SUBROUTINE GPICSplineInt_Deposit
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GPICSplineInt_Deposit1D(this,prop,np,proj)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  DESCRIPTION: Deposit1D
!               Carry out 1d deposition into grid points.
!
!  ARGUMENTS  :
!    this     : 'this' class instance
!    prop     : particle property to be interpolated
!    np       : no. particles (lenght of prop).
!    proj     : returned deposition of the property into regular grid
!               points. Must have dimensions set in constructor.
!-----------------------------------------------------------------
!$  USE threads
    IMPLICIT NONE
    CLASS(GPICSplineInt)                         :: this
    INTEGER      ,INTENT   (IN)                  :: np
    INTEGER                                      :: lag,k,nx,ny,nz
    REAL(KIND=GP),INTENT   (IN),DIMENSION(np)    :: prop
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(:,:,:) :: proj
    REAL(KIND=GP)                                :: zz(4),zzz,zzm
    REAL(KIND=GP)                                :: two,three,four,six,sixth,eighth
  
    nx = 1
    ny = 1
    nz = this%ldims_(3) + 2*this%gfcomm_%GetNumGhost()
  
!$omp parallel do if((ny*nz).gt.nth) private(i)
    DO k = 1,nz
      this%esplflddp_(1,1,k) = 0.0D0
    END DO
  
    IF (this%splord_.eq.0) THEN
!$omp parallel do
      DO lag=1, np
        this%esplflddp_(1,1,this%klg_(1,lag)) =  &
        this%esplflddp_(1,1,this%klg_(1,lag)) + prop(lag)
      ENDDO
    ELSE IF (this%splord_.eq.1) THEN
!$omp parallel do private(zz)
      DO lag=1, np
        zz(2) = this%zrk_(lag)
        zz(1) = 1.0_GP - zz(2)
        DO k=1,2
          this%esplflddp_(1,1,this%klg_(k,lag)) = &
          this%esplflddp_(1,1,this%klg_(k,lag)) + prop(lag)*zz(k)
        ENDDO
      ENDDO
    ELSE IF (this%splord_ .EQ. 2) THEN
      two    = 2.0_GP
      eighth = 1.0_GP/8.0_GP
!$omp parallel do private(zzz,zz)
      DO lag=1, np
        zzz   = this%zrk_(lag)
        zz(1) = eighth*(1-two*zzz)*(1-two*zzz)
        zz(3) = eighth*(1+two*zzz)*(1+two*zzz)
        zz(2) = 1.0_GP - zz(1) - zz(3)
        DO k=1,3
          this%esplflddp_(1,1,this%klg_(k,lag)) = &
          this%esplflddp_(1,1,this%klg_(k,lag)) + prop(lag)*zz(k)
        ENDDO
      ENDDO
    ELSE IF (this%splord_ .EQ. 3) THEN
      sixth  = 1.0_GP/6.0_GP
      six    = 6.0_GP
      four   = 4.0_GP
      three  = 3.0_GP
!$omp parallel do private(xxx,xxm,xx,zzz,zzm,zz)
      DO lag=1,np
        zzz = this%zrk_(lag)
        zzm = (1.0_GP-zzz)
        zz(1) = sixth*zzm*zzm*zzm
        zz(2) = sixth*(four+zzz*zzz*(three*zzz-six))
        zz(3) = sixth*(four+zzm*zzm*(three*zzm-six))
        zz(4) = 1.0_GP - zz(1) - zz(2) - zz(3)
        DO k=1,4
          this%esplflddp_(1,1,this%klg_(k,lag)) = &
          this%esplflddp_(1,1,this%klg_(k,lag)) + prop(lag)*zz(k)
        ENDDO
      ENDDO
    ENDIF
  
!$omp parallel do if((ny*nz).gt.nth) private(i)
    DO k = 1,nz
      this%esplfld_(1,1,k) = REAL(this%esplflddp_(1,1,k),kind=GP)
    END DO
  
    CALL GTStart(this%hdataex_)
    CALL this%gfcomm_%SlabDataReturnSF(proj,this%esplfld_,1)!UNPACK_SUM)
    CALL GTAcc(this%hdataex_)
  
    RETURN

  END SUBROUTINE GPICSplineInt_Deposit1D
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GPICSplineInt_Deposit2D(this,prop,np,proj)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  DESCRIPTION: Deposit2D
!               Carry out 3d deposition into grid points.
!
!  ARGUMENTS  :
!    this     : 'this' class instance
!    prop     : particle property to be interpolated
!    np       : no. particles (lenght of prop).
!    proj     : returned deposition of the property into regular grid
!               points. Must have dimensions set in constructor.
!-----------------------------------------------------------------
!$  USE threads
    IMPLICIT NONE
    CLASS(GPICSplineInt)                         :: this
    INTEGER      ,INTENT   (IN)                  :: np
    INTEGER                                      :: lag,i,k,nx,ny,nz
    REAL(KIND=GP),INTENT   (IN),DIMENSION(np)    :: prop
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(:,:,:) :: proj
    REAL(KIND=GP)                                :: xx(4),xxx,xxm
    REAL(KIND=GP)                                :: zz(4),zzz,zzm
    REAL(KIND=GP)                                :: two,three,four,six,sixth,eighth
  
    nx = this%ldims_(1)
    ny = this%ldims_(2)
    nz = this%ldims_(3) + 2*this%gfcomm_%GetNumGhost()
  
!$omp parallel do if((ny*nz).gt.nth) private(i)
    DO k = 1,nz
      DO i = 1,nx
        this%esplflddp_(i,1,k) = 0.0D0
      END DO
    END DO
  
    IF (this%splord_.eq.0) THEN
!$omp parallel do
      DO lag=1, np
        this%esplflddp_(this%ilg_(1,lag),1,this%klg_(1,lag)) =  &
        this%esplflddp_(this%ilg_(1,lag),1,this%klg_(1,lag)) + prop(lag)
      ENDDO
    ELSE IF (this%splord_.eq.1) THEN
!$omp parallel do private(xx,zz)
      DO lag=1, np
        xx(2) = this%xrk_(lag)
        xx(1) = 1.0_GP - xx(2)
        zz(2) = this%zrk_(lag)
        zz(1) = 1.0_GP - zz(2)
        DO k=1,2
          DO i=1,2
            this%esplflddp_(this%ilg_(i,lag),1,this%klg_(k,lag)) = &
            this%esplflddp_(this%ilg_(i,lag),1,this%klg_(k,lag)) &
                          + prop(lag)*xx(i)*zz(k)
          ENDDO
        ENDDO
      ENDDO
    ELSE IF (this%splord_ .EQ. 2) THEN
      two    = 2.0_GP
      eighth = 1.0_GP/8.0_GP
!$omp parallel do private(xxx,xx,zzz,zz)
      DO lag=1, np
        xxx   = this%xrk_(lag)
        xx(1) = eighth*(1-two*xxx)*(1-two*xxx)
        xx(3) = eighth*(1+two*xxx)*(1+two*xxx)
        xx(2) = 1.0_GP - xx(1) - xx(3)
   
        zzz   = this%zrk_(lag)
        zz(1) = eighth*(1-two*zzz)*(1-two*zzz)
        zz(3) = eighth*(1+two*zzz)*(1+two*zzz)
        zz(2) = 1.0_GP - zz(1) - zz(3)
        DO k=1,3
          DO i=1,3
            this%esplflddp_(this%ilg_(i,lag),1,this%klg_(k,lag)) = &
              this%esplflddp_(this%ilg_(i,lag),1,this%klg_(k,lag)) &
                          + prop(lag)*xx(i)*zz(k)
          ENDDO
        ENDDO
      ENDDO
    ELSE IF (this%splord_ .EQ. 3) THEN
      sixth  = 1.0_GP/6.0_GP
      six    = 6.0_GP
      four   = 4.0_GP
      three  = 3.0_GP
!$omp parallel do private(xxx,xxm,xx,zzz,zzm,zz)
      DO lag=1,np
        xxx = this%xrk_(lag)
        xxm = (1.0_GP-xxx)
        xx(1) = sixth*xxm*xxm*xxm
        xx(2) = sixth*(four+xxx*xxx*(three*xxx-six))
        xx(3) = sixth*(four+xxm*xxm*(three*xxm-six))
        xx(4) = 1.0_GP - xx(1) - xx(2) - xx(3)
  
        zzz = this%zrk_(lag)
        zzm = (1.0_GP-zzz)
        zz(1) = sixth*zzm*zzm*zzm
        zz(2) = sixth*(four+zzz*zzz*(three*zzz-six))
        zz(3) = sixth*(four+zzm*zzm*(three*zzm-six))
        zz(4) = 1.0_GP - zz(1) - zz(2) - zz(3)
        DO k=1,4
          DO i=1,4
            this%esplflddp_(this%ilg_(i,lag),1,this%klg_(k,lag)) = &
              this%esplflddp_(this%ilg_(i,lag),1,this%klg_(k,lag)) &
                          + prop(lag)*xx(i)*zz(k)
          ENDDO
        ENDDO
      ENDDO
    ENDIF
  
!$omp parallel do if((ny*nz).gt.nth) private(i)
    DO k = 1,nz
      DO i = 1,nx
        this%esplfld_(i,1,k) = REAL(this%esplflddp_(i,1,k),kind=GP)
      END DO
    END DO
  
    CALL GTStart(this%hdataex_)
    CALL this%gfcomm_%SlabDataReturnSF(proj,this%esplfld_,1)!UNPACK_SUM)
    CALL GTAcc(this%hdataex_)
  
    RETURN

  END SUBROUTINE GPICSplineInt_Deposit2D
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GPICSplineInt_Deposit3D(this,prop,np,proj)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  DESCRIPTION: Deposit3D
!               Carry out 3d deposition into grid points.
!               In addition, the normalized positions of the particles
!               in each interval must be computed prior to entry, as
!               must be the indices into the full grid that define the
!               range of control points. Temporary deposition are stored
!               in member data array.
!
!  ARGUMENTS  :
!    this     : 'this' class instance
!    prop     : particle property to be interpolated
!    np       : no. particles (lenght of prop).
!    proj     : returned deposition of the property into regular grid
!               points. Must have dimensions set in constructor.
!
!-----------------------------------------------------------------
!$  USE threads
    IMPLICIT NONE
    CLASS(GPICSplineInt)                         :: this
    INTEGER      ,INTENT   (IN)                  :: np
    INTEGER                                      :: lag,i,j,k,nx,ny,nz
    REAL(KIND=GP),INTENT   (IN),DIMENSION(np)    :: prop
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(:,:,:) :: proj
    REAL(KIND=GP)                                :: xx(4),xxx,xxm
    REAL(KIND=GP)                                :: yy(4),yyy,yym
    REAL(KIND=GP)                                :: zz(4),zzz,zzm
    REAL(KIND=GP)                                :: two,three,four,six,sixth,eighth
  
    nx = this%ldims_(1)
    ny = this%ldims_(2)
    nz = this%ldims_(3) + 2*this%gfcomm_%GetNumGhost()
  
  !$omp parallel do if((ny*nz).gt.nth) private(i) collapse(2)
    DO k = 1,nz
      DO j = 1,ny
        DO i = 1,nx
          this%esplflddp_(i,j,k) = 0.0D0
        END DO
      END DO
    END DO
  
    IF (this%splord_.eq.0) THEN
  !$omp parallel do
      DO lag=1, np
        this%esplflddp_(this%ilg_(1,lag),this%jlg_(1,lag),this%klg_(1,lag)) =  &
        this%esplflddp_(this%ilg_(1,lag),this%jlg_(1,lag),this%klg_(1,lag)) + prop(lag)
      ENDDO
    ELSE IF (this%splord_.eq.1) THEN
  !$omp parallel do private(xx,yy,zz)
      DO lag=1, np
        xx(2) = this%xrk_(lag)
        xx(1) = 1.0_GP - xx(2)
        yy(2) = this%yrk_(lag)
        yy(1) = 1.0_GP - yy(2)
        zz(2) = this%zrk_(lag)
        zz(1) = 1.0_GP - zz(2)
        DO k=1,2
          DO j=1,2
            DO i=1,2
              this%esplflddp_(this%ilg_(i,lag),this%jlg_(j,lag),this%klg_(k,lag)) = &
                this%esplflddp_(this%ilg_(i,lag),this%jlg_(j,lag),this%klg_(k,lag)) &
                            + prop(lag)*xx(i)*yy(j)*zz(k)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ELSE IF (this%splord_ .EQ. 2) THEN
      two    = 2.0_GP
      eighth = 1.0_GP/8.0_GP
  !$omp parallel do private(xxx,xx,yyy,yy,zzz,zz)
      DO lag=1, np
        xxx   = this%xrk_(lag)
        xx(1) = eighth*(1-two*xxx)*(1-two*xxx)
        xx(3) = eighth*(1+two*xxx)*(1+two*xxx)
        xx(2) = 1.0_GP - xx(1) - xx(3)
  
        yyy   = this%yrk_(lag)
        yy(1) = eighth*(1-two*yyy)*(1-two*yyy)
        yy(3) = eighth*(1+two*yyy)*(1+two*yyy)
        yy(2) = 1.0_GP - yy(1) - yy(3)
  
        zzz   = this%zrk_(lag)
        zz(1) = eighth*(1-two*zzz)*(1-two*zzz)
        zz(3) = eighth*(1+two*zzz)*(1+two*zzz)
        zz(2) = 1.0_GP - zz(1) - zz(3)
        DO k=1,3
          DO j=1,3
            DO i=1,3
              this%esplflddp_(this%ilg_(i,lag),this%jlg_(j,lag),this%klg_(k,lag)) = &
                this%esplflddp_(this%ilg_(i,lag),this%jlg_(j,lag),this%klg_(k,lag)) &
                            + prop(lag)*xx(i)*yy(j)*zz(k)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ELSE IF (this%splord_ .EQ. 3) THEN
      sixth  = 1.0_GP/6.0_GP
      six    = 6.0_GP
      four   = 4.0_GP
      three  = 3.0_GP
  !$omp parallel do private(xxx,xxm,xx,yyy,yym,yy,zzz,zzm,zz)
      DO lag=1,np
        xxx = this%xrk_(lag)
        xxm = (1.0_GP-xxx)
        xx(1) = sixth*xxm*xxm*xxm
        xx(2) = sixth*(four+xxx*xxx*(three*xxx-six))
        xx(3) = sixth*(four+xxm*xxm*(three*xxm-six))
        xx(4) = 1.0_GP - xx(1) - xx(2) - xx(3)
  
        yyy = this%yrk_(lag)
        yym = (1.0_GP-yyy)
        yy(1) = sixth*yym*yym*yym
        yy(2) = sixth*(four+yyy*yyy*(three*yyy-six))
        yy(3) = sixth*(four+yym*yym*(three*yym-six))
        yy(4) = 1.0_GP - yy(1) - yy(2) - yy(3)
  
        zzz = this%zrk_(lag)
        zzm = (1.0_GP-zzz)
        zz(1) = sixth*zzm*zzm*zzm
        zz(2) = sixth*(four+zzz*zzz*(three*zzz-six))
        zz(3) = sixth*(four+zzm*zzm*(three*zzm-six))
        zz(4) = 1.0_GP - zz(1) - zz(2) - zz(3)
        DO k=1,4
          DO j=1,4
            DO i=1,4
              this%esplflddp_(this%ilg_(i,lag),this%jlg_(j,lag),this%klg_(k,lag)) = &
                this%esplflddp_(this%ilg_(i,lag),this%jlg_(j,lag),this%klg_(k,lag)) &
                            + prop(lag)*xx(i)*yy(j)*zz(k)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDIF
  
  !$omp parallel do if((ny*nz).gt.nth) private(i) collapse(2)
    DO k = 1,nz
      DO j = 1,ny
        DO i = 1,nx
          this%esplfld_(i,j,k) = REAL(this%esplflddp_(i,j,k),kind=GP)
        END DO
      END DO
    END DO
  
  !  this%esplfld_ = this%esplflddp_
  
  !  proj = 0
    CALL GTStart(this%hdataex_)
    CALL this%gfcomm_%SlabDataReturnSF(proj,this%esplfld_,1)!UNPACK_SUM)
    CALL GTAcc(this%hdataex_)
  
    RETURN

  END SUBROUTINE GPICSplineInt_Deposit3D
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GPICSplineInt_PartUpdate(this,xp,yp,zp,np)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : PartUpdate
!  DESCRIPTION: Carries update update of interpolation positions.
!               Should be called before any interpolations are
!               done.
!
!               NOTE: All interp points _must_ reside in the sub-
!               domain represented by this task. No checking is
!               done to insure this, and bad indices will be
!               computed if this isn't done.
!  ARGUMENTS  :
!    this     : 'this' class instance
!    xp,yp,zp : particle x,y,z positions, or interpolation points
!               relative to grid set in constructor call
!    np       : no. interp points or particles.
!-----------------------------------------------------------------
    IMPLICIT NONE
    CLASS(GPICSplineInt)                  :: this
    REAL(KIND=GP),INTENT(IN),DIMENSION(*) :: xp,yp,zp
    REAL(KIND=GP)                         :: half
    INTEGER      ,INTENT(IN)              :: np
    INTEGER                               :: j,kmax,kmin,nx,ny,nz
    LOGICAL                               :: btmp

    ! Compute interval-normalized positions and
    ! indices into control point array in x, y directions:
    nx = this%ldims_(1)
    ny = this%ldims_(2)
    nz = this%ldims_(3)
    kmax = nz+this%gfcomm_%GetNumGhost()-1
    kmin =    this%gfcomm_%GetNumGhost()-1

    ! Note that the zp are on the _global_ grid, so
    ! we have to restrict to this MPI task's slab
    ! before computing local indices and spline coordinates....
    !
    ! First, do a check:
!$omp parallel do if(np.ge.NMIN_OMP) private(btmp)
    DO j = 1, np
      btmp = (zp(j).GE.this%xbnds_(3,1).AND.zp(j).LT.this%xbnds_(3,2))
      IF ( .NOT. btmp ) THEN
        WRITE(*,*) this%rank_, ' GPICSplineInt::PartUpdate3D: Invalid particle z-range'
        WRITE(*,*) this%rank_, ' GPICSplineInt::zbnd_0=',this%xbnds_(3,1),';  zbnd_1=',this%xbnds_(3,2), 'zp=',zp(j)
        STOP
      ENDIF
    ENDDO

    IF (this%splord_.eq.0) THEN
      half = 1.0_GP/2.0_GP
!$omp parallel do if(np.ge.NMIN_OMP)
      DO j = 1, np
        ! x-coord
        IF (nx.GT.1) THEN
          this%ilg_(1,j) = (xp(j)-this%xbnds_(1,1))*this%dxi_(1)
          this%xrk_  (j) = (xp(j)-this%xbnds_(1,1))*this%dxi_(1) - real(this%ilg_(1,j),kind=GP)
          IF (this%xrk_(j).gt.half) THEN
            this%ilg_(1,j) = this%ilg_(1,j) + 1
          ENDIF
          this%ilg_(1,j) = modulo(this%ilg_(1,j), nx) + 1
        END IF
        ! y-coord
        IF (ny.GT.1) THEN
          this%jlg_(1,j) = (yp(j)-this%xbnds_(2,1))*this%dxi_(2)
          this%yrk_  (j) = (yp(j)-this%xbnds_(2,1))*this%dxi_(2) - real(this%jlg_(1,j),kind=GP)
          IF (this%yrk_(j).gt.half) THEN
            this%jlg_(1,j) = this%jlg_(1,j) + 1
          ENDIF
          this%jlg_(1,j) = modulo(this%jlg_(1,j), ny) + 1
        END IF
        ! z-coord
        this%klg_(1,j) = (zp(j)-this%xbnds_(3,1))*this%dxi_(3)
        this%klg_(1,j) = max(min(this%klg_(1,j),kmax),kmin)
        this%zrk_  (j) = (zp(j)-this%xbnds_(3,1))*this%dxi_(3) - real(this%klg_(1,j),kind=GP)
        IF (this%zrk_(j).gt.half) THEN
          this%klg_(1,j) = this%klg_(1,j) + 1
        ENDIF
      ENDDO
    ELSE IF (this%splord_.eq.1) THEN
!$omp parallel do if(np.ge.NMIN_OMP)
      DO j = 1, np
        ! x-coord
        IF (nx.GT.1) THEN
          this%ilg_(1,j) = (xp(j)-this%xbnds_(1,1))*this%dxi_(1)
          this%xrk_  (j) = (xp(j)-this%xbnds_(1,1))*this%dxi_(1) - real(this%ilg_(1,j),kind=GP)
          this%ilg_(1,j) = modulo(this%ilg_(1,j), nx) + 1
          this%ilg_(2,j) = modulo(this%ilg_(1,j), nx) + 1
        END IF
        ! y-coord
        IF (ny.GT.1) THEN
          this%jlg_(1,j) = (yp(j)-this%xbnds_(2,1))*this%dxi_(2)
          this%yrk_  (j) = (yp(j)-this%xbnds_(2,1))*this%dxi_(2) - real(this%jlg_(1,j),kind=GP)
          this%jlg_(1,j) = modulo(this%jlg_(1,j), ny) + 1
          this%jlg_(2,j) = modulo(this%jlg_(1,j), ny) + 1
        END IF
        ! z-coord
        this%klg_(1,j) = (zp(j)-this%xbnds_(3,1))*this%dxi_(3)
        this%klg_(1,j) = max(min(this%klg_(1,j),kmax),kmin)
        this%zrk_  (j) = (zp(j)-this%xbnds_(3,1))*this%dxi_(3) - real(this%klg_(1,j),kind=GP)
        this%klg_(1,j) = this%klg_(1,j) + 1
        this%klg_(2,j) = this%klg_(1,j) + 1
      ENDDO
    ELSE IF (this%splord_.eq.2) THEN
      half = 1.0_GP/2.0_GP
!$omp parallel do if(np.ge.NMIN_OMP)
      DO j = 1, np
        ! x-coord
        IF (nx.GT.1) THEN
          this%ilg_(1,j) = (xp(j)-this%xbnds_(1,1))*this%dxi_(1)
          this%xrk_  (j) = (xp(j)-this%xbnds_(1,1))*this%dxi_(1) - real(this%ilg_(1,j),kind=GP)
          IF (this%xrk_(j).gt.half) THEN
            this%ilg_(1,j) = this%ilg_(1,j) + 1
            this%xrk_  (j) = this%xrk_  (j) - 1.0_GP
          ENDIF
          this%ilg_(2,j) = modulo(this%ilg_(1,j),nx) + 1
          this%ilg_(3,j) = modulo(this%ilg_(2,j),nx) + 1
          this%ilg_(1,j) = modulo(nx+this%ilg_(2,j)-2,nx) + 1
        END IF
        ! y-coord
        IF (ny.GT.1) THEN
          this%jlg_(1,j) = (yp(j)-this%xbnds_(2,1))*this%dxi_(2)
          this%yrk_  (j) = (yp(j)-this%xbnds_(2,1))*this%dxi_(2) - real(this%jlg_(1,j),kind=GP)
          IF (this%yrk_(j).gt.half) THEN
            this%jlg_(1,j) = this%jlg_(1,j) + 1
            this%yrk_  (j) = this%yrk_  (j) - 1.0_GP
          ENDIF
          this%jlg_(2,j) = modulo(this%jlg_(1,j),ny) + 1
          this%jlg_(3,j) = modulo(this%jlg_(2,j),ny) + 1
          this%jlg_(1,j) = modulo(nx+this%jlg_(2,j)-2,ny) + 1
        END IF
        ! z-coord
        this%klg_(1,j) = (zp(j)-this%xbnds_(3,1))*this%dxi_(3)
        this%klg_(1,j) = max(min(this%klg_(1,j),kmax),kmin)
        this%zrk_  (j) = (zp(j)-this%xbnds_(3,1))*this%dxi_(3) - real(this%klg_(1,j),kind=GP)
        IF (this%zrk_(j).gt.half) THEN
          this%klg_(1,j) = this%klg_(1,j) + 1
          this%zrk_  (j) = this%zrk_  (j) - 1.0_GP
        ENDIF
        this%klg_(2,j) = this%klg_(1,j) + 1
        this%klg_(3,j) = this%klg_(2,j) + 1
      ENDDO
    ELSE IF (this%splord_.eq.3) THEN
    ! x-coords:
      IF (nx.GT.1) THEN
!$omp parallel do 
        DO j = 1, np
          this%ilg_(1,j) = (xp(j)-this%xbnds_(1,1))*this%dxi_(1)
          this%xrk_  (j) = (xp(j)-this%xbnds_(1,1))*this%dxi_(1) - real(this%ilg_(1,j),kind=GP)
          this%ilg_(2,j) = modulo(this%ilg_(1,j),nx) + 1
          this%ilg_(3,j) = modulo(this%ilg_(2,j),nx) + 1
          this%ilg_(4,j) = modulo(this%ilg_(3,j),nx) + 1
          this%ilg_(1,j) = modulo(nx+this%ilg_(2,j)-2,nx) + 1
        ENDDO
      END IF
    ! y-coords:
      IF (ny.GT.1) THEN
!$omp parallel do 
        DO j = 1, np
          this%jlg_(1,j) = (yp(j)-this%xbnds_(2,1))*this%dxi_(2)
          this%yrk_  (j) = (yp(j)-this%xbnds_(2,1))*this%dxi_(2) - real(this%jlg_(1,j),kind=GP)
          this%jlg_(2,j) = modulo(this%jlg_(1,j),ny) + 1
          this%jlg_(3,j) = modulo(this%jlg_(2,j),ny) + 1
          this%jlg_(4,j) = modulo(this%jlg_(3,j),ny) + 1
          this%jlg_(1,j) = modulo(ny+this%jlg_(2,j)-2,ny) + 1
        ENDDO
      END IF
      ! z-coords:
!$omp parallel do
      DO j = 1, np
        this%klg_(1,j) = (zp(j)-this%xbnds_(3,1))*this%dxi_(3)
        this%klg_(1,j) = max(min(this%klg_(1,j),kmax),kmin)
        this%zrk_  (j) = (zp(j)-this%xbnds_(3,1))*this%dxi_(3) - real(this%klg_(1,j),kind=GP)
        this%klg_(2,j) = this%klg_(1,j) + 1
        this%klg_(3,j) = this%klg_(2,j) + 1
        this%klg_(4,j) = this%klg_(3,j) + 1
      ENDDO
    ELSE
      PRINT *, "GPICSPlineInt_PartUpdate: Spline order must be 0, 1, 2 or 3."
      STOP
    ENDIF
  END SUBROUTINE GPICSplineInt_PartUpdate
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GPICSplineInt_DoAlloc(this)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : DoAlloc
!  DESCRIPTION: performs allocations based on member data
!  ARGUMENTS  :
!    this     : 'this' class instance
!-----------------------------------------------------------------
    USE grid

    IMPLICIT NONE
    CLASS(GPICSplineInt)                   :: this
    INTEGER                                :: nzg

    CALL GPICSplineInt_DoDealloc(this)
    nzg = this%gfcomm_%GetNumGhost()

    IF ( this%rank_ .EQ. 2 ) THEN
      ALLOCATE( this%esplfld2_(this%ldims_(1),(this%ldims_(2)+2*nzg)) )
    ELSE IF ( this%rank_ .EQ. 3 ) THEN
      ALLOCATE(this%esplfld_(this%ldims_(1),this%ldims_(2),(this%ldims_(3)+2*nzg)) )
      ALLOCATE(this%esplflddp_(this%ldims_(1),this%ldims_(2),(this%ldims_(3)+2*nzg)) )
    ENDIF

    IF (ny.NE.1) THEN
      ALLOCATE(this%jlg_(this%splord_+1,this%maxint_))
      ALLOCATE(this%yrk_(this%maxint_))
    END IF
    IF (nx.NE.1) THEN
      ALLOCATE(this%ilg_(this%splord_+1,this%maxint_))    
      ALLOCATE(this%xrk_(this%maxint_))
    END IF
    ALLOCATE(this%klg_(this%splord_+1,this%maxint_))
    ALLOCATE(this%zrk_(this%maxint_))

  END SUBROUTINE GPICSplineInt_DoAlloc
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GPICSplineInt_DoDealloc(this)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : DoDealloc
!  DESCRIPTION: performs de-allocation of all class quantities
!  ARGUMENTS  :
!    this     : 'this' class instance
!-----------------------------------------------------------------

    IMPLICIT NONE
    CLASS(GPICSplineInt)                      :: this

    IF ( ALLOCATED (this%esplfld_  ) ) DEALLOCATE (this%esplfld_)
    IF ( ALLOCATED (this%esplfld2_ ) ) DEALLOCATE(this%esplfld2_)

    IF ( ALLOCATED       (this%ilg_) ) DEALLOCATE     (this%ilg_)
    IF ( ALLOCATED       (this%jlg_) ) DEALLOCATE     (this%jlg_)
    IF ( ALLOCATED       (this%klg_) ) DEALLOCATE     (this%klg_)
    IF ( ALLOCATED       (this%xrk_) ) DEALLOCATE     (this%xrk_)
    IF ( ALLOCATED       (this%yrk_) ) DEALLOCATE     (this%yrk_)
    IF ( ALLOCATED       (this%zrk_) ) DEALLOCATE     (this%zrk_)

  END SUBROUTINE GPICSplineInt_DoDealloc
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GPICSplineInt_Init(this)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : Init
!  DESCRIPTION: Does initialization of object. Called by constructor
!  ARGUMENTS  :
!    this     : 'this' class instance
!-----------------------------------------------------------------
    IMPLICIT NONE
    CLASS(GPICSplineInt)      :: this
    INTEGER                   :: j

    CALL GPICSplineInt_DoAlloc(this)

    ! Compute interval widths:
    DO j = 1, this%rank_
      this%dxi_(j) = 1.0_GP  !real(this%ldims_(j)-1,kind=GP)/ ( this%xbnds_(j,2) - this%xbnds_(j,1) )
    ENDDO

  END SUBROUTINE GPICSplineInt_Init
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GPICSplineInt_CompSpline(this,field)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : GPSplineInt_CompSpline
!  DESCRIPTION: Computes the spline coefficients, by trivially copying the
!               field to extended grid. Field must have dimensions set in
!               contructor. Spline coeffs are stored in member data array.
!
!  ARGUMENTS  :
!    this    : 'this' class instance
!    field   : Field to be interpolated. Must have dimensions set in
!              constructor.
!-----------------------------------------------------------------
!$  USE threads

    IMPLICIT NONE
    CLASS(GPICSplineInt)                              :: this
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(this%ntot_) :: field
    INTEGER                                           :: i,j,k,nx,ny,nz

    nx = this%ldims_(1)
    ny = this%ldims_(2)
    nz = this%ldims_(3) + 2*this%gfcomm_%GetNumGhost()
!$omp parallel do if(ny*nz.gt.nth) private(i) collapse(2)
    DO k = 1,nz
      DO j = 1,ny
        DO i = 1,nx
          this%esplfld_(i,j,k) = 0.0_GP
        END DO
      END DO
    END DO

!    this%esplfld_ = 0

    CALL GTStart(this%hdataex_)
    CALL this%gfcomm_%SlabDataExchangeSF(this%esplfld_,field)
    CALL GTAcc(this%hdataex_)

   RETURN

 END SUBROUTINE GPICSplineInt_CompSpline
!-----------------------------------------------------------------
!-----------------------------------------------------------------

 SUBROUTINE GPICSplineInt_ResizeArrays(this,newmparts,onlyinc)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : Resize_Arrays
!  DESCRIPTION: Resize all arrays in the GPSplineInt class.
!               Not compatible with offloading (for now).
!
!  ARGUMENTS  :
!    this     : 'this' class instance
!    newmparts: new number of particles
!    onlyinc  : if true, will only resize to increase array size
!-----------------------------------------------------------------
    IMPLICIT NONE
    CLASS(GPICSplineInt),INTENT(INOUT)   :: this
    INTEGER             ,INTENT(IN)      :: newmparts
    LOGICAL             ,INTENT(IN)      :: onlyinc
    INTEGER                              :: n

    n = SIZE(this%xrk_)
    IF ((n.lt.newmparts).OR.((n.gt.newmparts).AND..NOT.onlyinc)) THEN
      IF (this%ldims_(1).GT.1) CALL Resize_IntArrayRank2(this%ilg_,newmparts,.false.)
      IF (this%ldims_(2).GT.1) CALL Resize_IntArrayRank2(this%jlg_,newmparts,.false.)
      CALL Resize_IntArrayRank2(this%klg_,newmparts,.false.)
      IF (this%ldims_(1).GT.1) CALL Resize_ArrayRank1(this%xrk_   ,newmparts,.false.)
      IF (this%ldims_(2).GT.1) CALL Resize_ArrayRank1(this%yrk_   ,newmparts,.false.)
      CALL Resize_ArrayRank1(this%zrk_   ,newmparts,.false.)
      this%maxint_ = newmparts
    END IF

    CALL this%gfcomm_%ResizeArrays(newmparts,onlyinc)

    RETURN

 END SUBROUTINE GPICSplineInt_ResizeArrays
!-----------------------------------------------------------------
!-----------------------------------------------------------------

END MODULE class_GPICSplineInt
