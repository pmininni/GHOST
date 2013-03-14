!=================================================================
! GHOST GSplineInt interpolation class
! The spline used (that on the 'momoents' or 2nd derivatives)
! as well as the tridiagonal solve for the moments, are
! derived from the presentation given in the book
! 'The theory of splines and their application'
! 1967, by Ahlberg, Nilson, Walsh
!
! 2012 D. Rosenberg
!      NCAR
!
! 15 May 2012: Initial version
!=================================================================
MODULE class_GSplineInt
      USE mpivars
      USE fprecision
      USE gmxm
      IMPLICIT NONE

      ENUM, BIND(C) :: GPINITTYPE
        ENUMERATOR :: GSPLINEINIT_NOFFT=0
        ENUMERATOR :: GSPLINEINIT_FFT   
      ENDENUM GPINITTYPE
      
      TYPE, PUBLIC :: GSplineInt
        PRIVATE
        ! Member data:
        REAL(KIND=GP),ALLOCATABLE,DIMENSION  (:,:)   :: 
        REAL(KIND=GP),ALLOCATABLE,DIMENSION    (:)   :: d_,q_,s_,t_,u_,v_ 
        REAL(KIND=GP)                                :: 
        TYPE(GPINIT)                                 :: inttype_
        INTEGER                                      :: maxint_
        INTEGER                                      :: ierr_,idims_(3),ntot_
        INTEGER                                      :: nrank_
        CHARACTER(len=1024)                          :: serr_
      CONTAINS
        ! Public methods:
        PROCEDURE,PUBLIC :: GSplineInt_ctor
        FINAL            :: GSplineInt_dtor

        ! ...Aliases:
        PROCEDURE,PUBLIC :: init        => GSplineInt_Init
        PROCEDURE,PUBLIC :: DoInterp    => GSplineInt_Interp3D 
        PROCEDURE,PUBLIC :: SetType     => GSplineInt_SetType
        PROCEDURE,PUBLIC :: GetType     => GSplineInt_GetType
        PROCEDURE,PUBLIC :: SetRank     => GSplineInt_SetRank
        PROCEDURE,PUBLIC :: GetRank     => GSplineInt_GetRank
      END TYPE GSplineInt

      ! Private methods:
      PRIVATE :: GSplineInt_Init          , GSplineInt_Interp3D
      PRIVATE :: GSplineInt_SetType       , GSplineInt_GetType
      PRIVATE :: GSplineInt_SetRank       , GSplineInt_GetRank
      PRIVATE :: GSplineInt_DoAlloc       , GSplineInt_DoDealloc
      PRIVATE :: GSplineInt_DoTridiag     ,GSplineInt_FindMoments1d
      PRIVATE :: GSplineInt_FindMoments2d , GSplineInt_FindMoments3d

! Methods:
  CONTAINS

  SUBROUTINE GSplineInt_ctor(this,rank,idims,inittype,maxpart)
!-----------------------------------------------------------------
!  Main explicit constructor
!  ARGUMENTS:
!    this    : 'this' class instance
!    rank    : rank of solution (1, 2, or 3)
!    idims   : array of spatial dimensions of each rank direction, x, y, or z
!    inttype : GPINITTYPE of polynomial interpolation to set up for
!    maxpart : max no. interpolation points/Lag. particles
!-----------------------------------------------------------------
    USE grid
    USE mpivars

    IMPLICIT NONE
    CLASS(GSplineInt)                 :: this
    TYPE(GPINITTYPE), INTENT(IN)      :: inttype
    INTEGER,INTENT(IN)                :: rank,idims(*),maxpart

    INTEGER                           :: j

    this%maxint_  = maxpart
    this%inttype_ = inttype
    this%nrank_   = rank

    IF ( this%nrank_.LT.1 .OR. this%nrank_.GT.3 ) THEN
      WRITE(*,*)'GSplineInt::ctor: Invalid rank'
      STOP
    ENDIF

    ntot_  = 1
    idims_ = 0
    DO j = 1, rank
      idims_(j) = idims(j)
      ntot_ = ntot_ * idims(j)
    ENDDO

    CALL GSplineInt_DoAlloc()

  END SUBROUTINE GSplineInt_ctor
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GSplineInt_dtor(this)
!-----------------------------------------------------------------
!  Main explicit destructor
!  ARGUMENTS:
!    this    : 'this' class instance
!-----------------------------------------------------------------

    IMPLICIT NONE
    CLASS(GSplineInt)                      :: this

    CALL GSplineInt_DoDealloc()

  END SUBROUTINE GSplineInt_dtor
!-----------------------------------------------------------------
!-----------------------------------------------------------------


  SUBROUTINE GSplineInt_DoAlloc(this)
!-----------------------------------------------------------------
!  METHOD     : DoAlloc
!  DESCRIPTION: performs allocations based on member data
!  ARGUMENTS  : 
!    this     : 'this' class instance
!-----------------------------------------------------------------

    IMPLICIT NONE
    CLASS(GSplineInt)                      :: this
  
    CALL GSplineInt_DoDealloc(this)

    ALLOCATE(MX_(this%idims_(1),this%maxint_))
    IF ( this%nrank_.GT. 1 ) THEN
      ALLOCATE(MY_(this%idims_(1),this%maxint_))
    ENDIF
    IF ( this%nrank_.GT. 2 ) THEN
      ALLOCATE(MZ_(this%idims_(1),this%maxint_))
    ENDIF
    ALLOCATE(icell_  (this%maxint_,this%nrank_))
    ALLOCATE(d_(0:this%maxint_) )
    ALLOCATE(q_(0:this%maxint_) )
    ALLOCATE(s_(0:this%maxint_) )
    ALLOCATE(t_(0:this%maxint_) )
    ALLOCATE(u_(0:this%maxint_) )
    ALLOCATE(v_(0:this%maxint_) )

  END SUBROUTINE GSplineInt_DoAlloc
!-----------------------------------------------------------------
!-----------------------------------------------------------------


  SUBROUTINE GSplineInt_DoDealloc(this)
!-----------------------------------------------------------------
!  METHOD     : DoDealloc
!  DESCRIPTION: performs de-allocation of all class quantities
!  ARGUMENTS  : 
!    this     : 'this' class instance
!-----------------------------------------------------------------

    IMPLICIT NONE
    CLASS(GSplineInt)                      :: this


    IF ( ALLOCATED(MX_ ) ) DEALLOCATE(MX_)
    IF ( ALLOCATED(MY_ ) ) DEALLOCATE(MY_)
    IF ( ALLOCATED(MZ_ ) ) DEALLOCATE(MZ_)
    IF ( ALLOCATED(d_  ) ) DEALLOCATE(d_)
    IF ( ALLOCATED(q_  ) ) DEALLOCATE(q_)
    IF ( ALLOCATED(s_  ) ) DEALLOCATE(s_)
    IF ( ALLOCATED(t_  ) ) DEALLOCATE(t_)
    IF ( ALLOCATED(u_  ) ) DEALLOCATE(u_)
    IF ( ALLOCATED(v_  ) ) DEALLOCATE(v_)

  END SUBROUTINE GSplineInt_DoDealloc
!-----------------------------------------------------------------
!-----------------------------------------------------------------


  SUBROUTINE GSplineInt_SetRank(this,rank)
!-----------------------------------------------------------------
!  METHOD     : SetRank
!  DESCRIPTION: Set rank (dimension) of interpolation 
!  ARGUMENTS  : 
!    this    : 'this' class instance
!    rank    : integer rank
!-----------------------------------------------------------------
    IMPLICIT NONE
    CLASS(GSplineInt)                      :: this
    INTEGER,INTENT(IN)                   :: rank

    this%nrank_ = rank
    IF ( this%nrank_.LT.1 .OR. this%nrank_.GT.3 ) THEN
      WRITE(*,*)'GSplineInt::SetRank: Invalid rank'
      STOP
    ENDIF
    CALL GSplineInt_DoDealloc()
    CALL GSplineInt_DoAlloc()

  END SUBROUTINE GSplineInt_SetRank
!-----------------------------------------------------------------
!-----------------------------------------------------------------


  INTEGER FUNCTION GSplineInt_GetRank(this)
!-----------------------------------------------------------------
!  METHOD     : GetRank
!  DESCRIPTION: Set rank (dimension) of interpolation 
!  ARGUMENTS  : 
!    this    : 'this' class instance
!    rank    : integer rank
!-----------------------------------------------------------------
    IMPLICIT NONE
    CLASS(GSplineInt)                      :: this

    GSplineInt_GetRank = this%nrank_

  END SUBROUTINE GSplineInt_GetRank
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GSplineInt_SetType(this,itype)
!-----------------------------------------------------------------
!  METHOD     : SetType
!  DESCRIPTION: Set init type of interp.
!  ARGUMENTS  : 
!    this    : 'this' class instance
!    itype   : GPINITTYPE type 
!-----------------------------------------------------------------
    IMPLICIT NONE
    CLASS(GSplineInt)                    :: this
    TYPE(GPINITTYPE),INTENT(IN)          :: itype

    this%inttype_ = itype

  END SUBROUTINE GSplineInt_SetType
!-----------------------------------------------------------------
!-----------------------------------------------------------------


  GPINITTYPE FUNCTION GSplineInt_GetType(this)
!-----------------------------------------------------------------
!  METHOD     : GetType
!  DESCRIPTION: Get init type of interp.
!  ARGUMENTS  : 
!    this    : 'this' class instance
!-----------------------------------------------------------------
    IMPLICIT NONE
    CLASS(GSplineInt)                      :: this

    GSplineInt_GetType = this%inttype_

  END SUBROUTINE GSplineInt_GetType
!-----------------------------------------------------------------
!-----------------------------------------------------------------


  SUBROUTINE GSplineInt_FindMoments1d(this,v,n)
!-----------------------------------------------------------------
!  METHOD     : FindMoments1d
!  DESCRIPTION: Computes the spline: finds 'mooments' of spline, 
!            defined as the second derivatives at the spatial
!            Lagrange points. The moment equation is
!        0.5 M_(j-1) + 2M_j + 0.5 M_(j+1) = d_j
! where
!
!        d_j = 6*{ [(y_(j+1)-y_j)/h_(j+1)] - [(y_j-y_(j-1))/h_j] }/(h_j + j_(j+1))
! and
!        h_j = 1 for all j on our grids.
!            
!  ARGUMENTS  :
!    this    : 'this' class instance
!    v       : 1-d field to be interpolated
!    n       : no. interp. points
!-----------------------------------------------------------------
    IMPLICIT NONE
    CLASS(GSplineInt)                      :: this
    REAL(KIND=GP),INTENT(IN)               :: v(*)

    INTEGER                                :: k

    d_(1) = 3.0_GP*(v(2)-2.0_GP*v(1)+v(n))
    d_(n) = 3.0_GP*(v(1)-2.0_GP*v(n)+v(n-1))
    DO k = 2, n-1
      d_(k) = 3.0_GP*(v(k+1)-2.0_GP*v(k)+v(k-1))
    ENDDO
   
    CALL GSplineInt_DoTridiag(this,d_,n,Mx_)

  END SUBROUTINE GSplineInt_FindMoments1d
!-----------------------------------------------------------------
!-----------------------------------------------------------------


  SUBROUTINE GSplineInt_FindMoments2d(this,idir,v,jc,n,jsta,jend)
!-----------------------------------------------------------------
!  METHOD     : FindMoments2d
!  DESCRIPTION: Computes the spline: finds 'mooments' of spline, 
!            defined as the second derivatives at the spatial
!            Lagrange points. The moment equation is
!        0.5 M_(j-1) + 2M_j + 0.5 M_(j+1) = d_j
! where
!
!        d_j = 6*{ [(y_(j+1)-y_j)/h_(j+1)] - [(y_j-y_(j-1))/h_j] }/(h_j + j_(j+1))
! and
!        h_j = 1 for all j on our grids.
!            
!  ARGUMENTS  :
!    this    : 'this' class instance
!    v       : 2-d field to be interpolated, treated as a 1d array
!    idir    : coordinate direction to consider
!    jc      : direction held constant (1, ..., nx) or ( 1, ..., ny)
!    n       : no. interp. points in x, y; also size of v(nx,ny)
!    jsta,jend: 2-indices of y-decomposition
!-----------------------------------------------------------------
    IMPLICIT NONE
    CLASS(GSplineInt)                      :: this
    REAL(KIND=GP),INTENT(IN)               :: v(*)
    INTEGER,INTENT(IN)                     :: idir,jc,nx,ny

    INTEGER                                :: jm,k,nj

    STOP 'Not implemented yet'

#if 0
    IF ( idir.EQ.1 ) THEN
      !  jc==y-index
      jm = jc - 1
      nj = nx*jm
      d_(1) = 3.0_GP*(v(2+nx*jm)-2.0_GP*v(1+nx*jm)+v(n+nx*jm))
      d_(n) = 3.0_GP*(v(1+nx*jm)-2.0_GP*v(n+nx*jm)+v(n-1+nx*jm))
      DO k = 2, nx-1
        d_(k) = 3.0_GP*(v(k+1+nj)-2.0_GP*v(k+nj)+v(k-1+nj))
      ENDDO
      CALL GSplineInt_DoTridiag(this,d_,nx,Mx_)
    ELSE IF ( idir.EQ.2 ) THEN
      !  jc==x-index
      d_(1) = 3.0_GP*(v(jc+nx)-2.0_GP*v(jc)+v(jc+nx*(n-1)))
      d_(n) = 3.0_GP*(v(jc)-2.0_GP*v(jc+nx*(n-1))+v(jc+nx*(n-2)))
      DO k = 2, ny-1
        d_(k) = 3.0_GP*(v(jc+nx*(k))-2.0_GP*v(jc+nx*(k-1))+v(jc+nx*(k-2)))
      ENDDO
      CALL GSplineInt_DoTridiag(this,d_,ny,My_)
    ENDIF
#endif


  END SUBROUTINE GSplineInt_FindMoments2d
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GSplineInt_FindMoments3d(this,idir,v,jc,kc,n,ksta,kend,kghost)
!-----------------------------------------------------------------
!  METHOD     : FindMoments3d
!  DESCRIPTION: Computes the spline: finds 'moments' of spline, 
!            defined as the second derivatives at the spatial
!            Lagrange points. The moment equation is
!        0.5 M_(j-1) + 2M_j + 0.5 M_(j+1) = d_j
! where
!
!        d_j = 6*{ [(y_(j+1)-y_j)/h_(j+1)] - [(y_j-y_(j-1))/h_j] }/(h_j + j_(j+1))
! and
!        h_j = 1 for all j on our grids.
!
! Data, v, must enter method as x-y complete, and distributed
! in z. If idir==3, this data is transposed, and the z-moments
! are computed from the transposed data.
!            
!  ARGUMENTS  :
!    this    : 'this' class instance
!    v       : 3-d field to be interpolated, treated as a 1d array
!    idir    : coordinate direction to consider
!    jc,kc   : directions held constant, in cyclic order; refer to
!              fundamental grid of GHOST, not the one padded with
!              kghost zones in z-direction
!    n       : no. interp. points in x, y; also size of v(n,n,ksta:kend)
!    ksta,kend: k-range of z-decomposition
!    kghost  : number of 'ghost' in z-direction for top and bottom of
!              slab.
!-----------------------------------------------------------------
    IMPLICIT NONE
    CLASS(GSplineInt)                      :: this
    REAL(KIND=GP),INTENT(IN)               :: v(*)
    INTEGER,INTENT(IN)                     :: idir,jc,kc,ksta,kend,n

    INTEGER                                :: im,jm,k,km,nj,nk,nn

! Using fact that h_j = 1 for our regular grids, we
! have that
!      d_j = 3 * ( (y_(j+1)-2*y_j+y_(j-1) )
! for the periodic values in x-y (complete, so of size n).
!
! For z-direction, we assume kghost extra zones pad top and bottom of
! of slab, and use non-periodic splines in that direction 
    nn = n*n
    IF ( idir.EQ.1 ) THEN
      !  jc==y-index, kc=z-index
      jm = jc-1
      km = kc+kghost-1
      nj = n*jm
      nk = nn*km
      d_(1) = 3.0_GP*(v(2+nj+nk)-2.0_GP*v(1+nj+nk)+v(n+nj+nk))
      d_(n) = 3.0_GP*(v(1+nj+nk)-2.0_GP*v(n+nj+nk)+v(n-1+nj+nk))
      DO k = 2, n-1
        d_(k) = 3.0_GP*(v(k+1+nj+nk)-2.0_GP*v(k+nj+nk)+v(k-1+nj+nk))
      ENDDO
      CALL GSplineInt_DoTridiag(this,d_,n,Mx_)
    ELSE IF ( idir.EQ.2 ) THEN
      jm = jc-1
      km = kc+kghost-1
      nk = nn*km
      !  jc==z-index, kc=x-index
      d_(1) = 3.0_GP*(v(kc+n +nk)-2.0_GP*v(kc+nn*km)+v(kc+n *(n-1)+nk))
      d_(n) = 3.0_GP*(v(kc+nk)-2.0_GP*v(kc+n *(n-1)+nk)+v(kc+n *(n-2)+nk))
      DO k = 2, n-1
        d_(k) = 3.0_GP*(v(kc+n*(k)+nk)-2.0_GP*v(jc+n*(k-1)+nk)+v(jc+n*(k-2)+nk))
      ENDDO
      CALL GSplineInt_DoTridiag(this,d_,n,My_)
    ELSE IF ( idir.EQ.3 ) THEN
!     Transpose data, so that it is z-y complete, then
!     strip out the z-dependence:
      CALL GSplineInt_Transpose(v,vt_,n,ksta,kend)
      !  jc==x-index, kc=y-index
      jm = jc-1
      km = kc-1
      nj = n *km
      nk = nn*jm
      d_(1) = 3.0_GP*(vt_(2+n*jm+nn*km)-2.0_GP*vt_(1+n*jm+nn*km)+vt_(n+n*jm+nn*km))
      d_(n) = 3.0_GP*(vt_(1+nj+nk)-2.0_GP*vt_(n+nj+nk)+vt_(n-1+nj+nk))
      DO k = 2, kend-ksta
        d_(k) = 3.0_GP*(vt_(k+1+nj+nk)-2.0_GP*vt_(k+nj+nk)+vt_(k-1+nj+nk))
      ENDDO
      CALL GSplineInt_DoTridiag(this,d_,n,Mz_)
    ENDIF


  END SUBROUTINE GSplineInt_FindMoments3d
!-----------------------------------------------------------------
!-----------------------------------------------------------------


  SUBROUTINE GSplineInt_DoTridiag(this,d,n,x)
!-----------------------------------------------------------------
!  METHOD     : DoTridiag
!  DESCRIPTION: Carries out solve of (periodic) tridiagonal system
!            b1 x1     + c1 x2         + ... a1 xn = d1
!            a2 x1     + b2 x2         + ... c2 xn = d2
!                       ...
!            a(n-1) x1 + b(n-1) x(n-1) + ... bn xn = dn
!             
!            Matrix given by coeffieient arrays a(*), b(*), c(*), and
!            rhs given by d(*), where dimension is the size
!            of the coordinate system, n. solution, x, returned
!            in array x
!            
!            Method uses the fact that we are on an evenly 
!            spaced grid of interval width 1, so, a(k)=c(k) = 0.5, 
!            b(k) = 2. The algorithm used is that from the Ahlberg
!            text, pp. 14-15.
!            
!  ARGUMENTS  :
!    this    : 'this' class instance
!-----------------------------------------------------------------
    IMPLICIT NONE
    CLASS(GSplineInt)                      :: this

    REAL(KIND=GP),INTENT(IN)               :: d
    REAL(KIND=GP),INTENT(OUT)              :: x
    INTEGER,INTENT(IN)                     :: n
    REAL(KIND=GP)                          :: xnc,p,pinv
    INTEGER                                :: i,j,k

    q_(0) = 0.0_GP
    u_(0) = 0.0_GP
    s_(0) = 0.0_GP
    DO k = 1, n
      pinv  = 1.0_GP/(0.5_GP*q_(k-1) + 2.0_GP)
      q_(k) = -0.5_GP*pinv
      u_(k) = (d(k) - 0.5_GP*u_(k-1)) * pinv
      s_(k) = -0.5_GP*s_(k-1)*pinv
    ENDDO

    t_(n) = 0.0_GP
    v_(n) = 0.0_GP
    DO k = n-1, 1, -1
      t_(k) = q_(k)*t_(k+1) + s_(k)
      v_(k) = q_(k)*v_(k+1) + u_(k)
    ENDDO
    xnc  = 0.5*(t_(1) + t_(n-1)) + 2.0_GP
    x(n) = (d(n) - 0.5_GP*(v_(1) + v_(n-1)))/xnc
    xnc  = x(n)

    DO k = n-1, 1, -1
      x(k) = t_(k)*xnc + v_(k)
    ENDDO

  END SUBROUTINE GSplineInt_DoTridiag
!-----------------------------------------------------------------
!-----------------------------------------------------------------


END MODULE class_GSplineInt
