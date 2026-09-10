!=================================================================
! GPSplineInt: cubic spline interpolation of a slab-decomposed
! field at the particle positions.
!
! PartUpdate3D updates the interpolation points.
! Then, two ways of computing the spline coefficients are available:
!   CompSpline3D: solves the periodic tridiagonal systems in the
!                 three directions on the host (transposes for z)
!                 and fills the extended field;
!   SetCoeffs3D : takes a field that already holds the spline
!                 coefficients (computed by the caller) and only
!                 fills the extended field.
! Finally, DoInterp3D carries out the actual interpolation.
!
! The spline and the tridiagonal solves are derived from several
! sources including the presentation given in the book 'The theory
! of splines and their application', 1967, Ahlberg, Nilson, Walsh.
!
! The interpolation reads the spline coefficients on the extended
! grid esplfld_(nx,ny,nzl+2*nzghost), with nzghost planes of the
! neighbor slabs below and above the local slab, filled by the
! ghost-plane exchange of GPartComm. The cubic B-spline stencil of
! a point between the planes k and k+1 spans the planes k-1 to k+2,
! so a particle in the first cell above the slab needs three ghost
! planes above it (and one in the first cell below needs two):
! GPSI_NZGHOST = 3 planes are exchanged on each side, and particles
! may be interpolated up to one cell above and two cells below the
! slab (this happens in the stages of multi-stage steppers, which
! exchange particles between tasks only at the end of the step).
! PartUpdate3D stops if a particle lies outside that range (the
! time step is too long for the stepper). The particle-sized arrays
! (cell indices ilg_,jlg_,klg_, fractional positions xrk_,yrk_,
! zrk_ and weights wrkl_) and the coefficient field have device
! copies in offload builds; PartUpdate3D and DoInterp3D run their
! kernels on the device while gdev_active is set. The kernels are
! module procedures taking explicit-shape arrays.
!
! 2013: A. Pumir (ENS, Lyon)
!       D. Rosenberg (NCCS: ORNL) - Initial version
! 2016: CompSpline3D loops optimized for speed (P. Mininni)
! 2026: GPU offloading (P. Mininni)
!=================================================================
MODULE class_GPSplineInt
      USE mpivars
      USE fprecision
      USE class_GPartComm
      USE gtimer
      USE gmem
      USE gdevice, ONLY: gdev_active
      IMPLICIT NONE
      ! Ghost planes on each side of the slab (see the header comment)
      INTEGER, PARAMETER, PUBLIC :: GPSI_NZGHOST = 3

      PRIVATE
      TYPE, PUBLIC :: GPSplineInt
        PRIVATE
        ! Device-resident arrays
        REAL(KIND=GP),ALLOCATABLE,DIMENSION(:,:,:) :: esplfld_
        REAL(KIND=GP),ALLOCATABLE,DIMENSION    (:) :: xrk_,yrk_,zrk_
        REAL(KIND=GP),ALLOCATABLE,DIMENSION  (:,:) :: wrkl_
        INTEGER      ,ALLOCATABLE,DIMENSION  (:,:) :: ilg_,jlg_,klg_
        REAL(KIND=GP),ALLOCATABLE,DIMENSION    (:) :: tmptr_,tmpt2_  ! z-complete layout (il,ny,nz)
        ! Factorization of the tridiagonal systems (computed on the
        ! host; the five arrays the solves use have device copies)
        REAL(KIND=GP),ALLOCATABLE,DIMENSION    (:) :: ax_,bx_,betx_,cx_,gamx_,px_,xxx_
        REAL(KIND=GP),ALLOCATABLE,DIMENSION    (:) :: ay_,by_,bety_,cy_,gamy_,py_,xxy_
        REAL(KIND=GP),ALLOCATABLE,DIMENSION    (:) :: az_,bz_,betz_,cz_,gamz_,pz_,xxz_
        REAL(KIND=GP)                              :: dxi_(3),xbnds_(3,2),zetax_,zetay_,zetaz_
        REAL(KIND=GP)                              :: zlo_ ! lowest valid particle z
        TYPE(GPartComm),POINTER                    :: gpcomm_
        INTEGER                                    :: maxint_
        INTEGER                                    :: ierr_,ider_(3),nd_(3)
        INTEGER                                    :: ibnds_(3,2),obnds_(3,2)
        INTEGER                                    :: ldims_(3),odims_(3)
        INTEGER                                    :: hdataex_,htransp_
        INTEGER                                    :: ntot_,ttot_
        INTEGER                                    :: rank_
      CONTAINS
        PROCEDURE,PUBLIC :: GPSplineInt_ctor
        FINAL            :: GPSplineInt_dtor
        PROCEDURE,PUBLIC :: Init         => GPSplineInt_Init
        PROCEDURE,PUBLIC :: DoInterp3D   => GPSplineInt_Interp3D
        PROCEDURE,PUBLIC :: SetDeriv     => GPSplineInt_SetDeriv
        PROCEDURE,PUBLIC :: PartUpdate3D => GPSplineInt_PartUpdate3D
        PROCEDURE,PUBLIC :: CompSpline3D => GPSplineInt_CompSpline3D
        PROCEDURE,PUBLIC :: SetCoeffs3D  => GPSplineInt_SetCoeffs3D
        PROCEDURE,PUBLIC :: ResizeArrays => GPSplineInt_ResizeArrays
      END TYPE GPSplineInt

  CONTAINS

!=================================================================
! Constructor, destructor, allocator, PartUpdate3D and helpers
!=================================================================

!-----------------------------------------------------------------
!  METHOD     : GPSplineInt_ctor
!  DESCRIPTION: Constructor. Stores the grid and slab bounds, the
!               transposed bounds used by CompSpline3D, the
!               particle buffer size and the timer handles. In z,
!               xbnds_(3,1) is the origin of the extended grid
!               (the plane below the first ghost plane, so that
!               INT(z-xbnds_(3,1)) is the extended index of the
!               plane at or below z), xbnds_(3,2) and zlo_ are the
!               bounds of the z range that can be interpolated.
!-----------------------------------------------------------------
  SUBROUTINE GPSplineInt_ctor(this,rank,nd,ibnds,xbnds,obnds,maxpart,gpcomm, &
                              hdataex,htransp)
    IMPLICIT NONE
    CLASS(GPSplineInt)                           :: this
    TYPE(GPartComm),TARGET                       :: gpcomm
    INTEGER        ,INTENT(IN)                   :: hdataex,htransp,maxpart,rank
    INTEGER                                      :: nzg
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
    IF ( this%rank_.NE.3 ) THEN
      WRITE(*,*)'GPSplineInt::ctor: only rank 3 is supported'
      STOP
    ENDIF
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
    ! Extended grid in z: nzg ghost planes on each side. A stencil
    ! (4 planes from the one below the point) fits in the extended
    ! grid for points from the second ghost plane below the slab to
    ! the end of the first ghost cell above it.
    nzg = gpcomm%GetNumGhost()
    IF ( nzg.LT.3 ) THEN
      WRITE(*,*)'GPSplineInt::ctor: at least 3 ghost planes are needed'
      STOP
    ENDIF
    ! (here xbnds_(3,1) = ksta-1 and xbnds_(3,2) = kend-1 are the
    ! first and last planes of the slab)
    this%zlo_         = this%xbnds_(3,1)-real(nzg-1,kind=GP) ! ksta-3
    this%xbnds_(3,2)  = this%xbnds_(3,2)+real(nzg-1,kind=GP) ! kend+1
    this%xbnds_(3,1)  = this%xbnds_(3,1)-real(nzg+1,kind=GP) ! ksta-5
    CALL GPSplineInt_Init(this)
  END SUBROUTINE GPSplineInt_ctor


!-----------------------------------------------------------------
! Main destructor
!-----------------------------------------------------------------
  SUBROUTINE GPSplineInt_dtor(this)
    IMPLICIT NONE
    TYPE(GPSplineInt),INTENT(INOUT) :: this
    CALL GPSplineInt_DoDealloc(this)
  END SUBROUTINE GPSplineInt_dtor


!-----------------------------------------------------------------
!  METHOD     : Init
!  DESCRIPTION: Allocates the arrays and factorizes the periodic
!               tridiagonal spline matrices of the three directions
!-----------------------------------------------------------------
  SUBROUTINE GPSplineInt_Init(this)
    IMPLICIT NONE
    CLASS(GPSplineInt)        :: this
    INTEGER                   :: j

    CALL GPSplineInt_DoAlloc(this)
    DO j = 1, this%rank_
      this%dxi_(j) = 1.0_GP
    ENDDO
    CALL GPSplineInt_MatInvQ(this,this%nd_(1),this%ax_,this%bx_,this%cx_,&
         this%px_,this%gamx_,this%betx_,this%xxx_,this%zetax_)
    CALL GPSplineInt_MatInvQ(this,this%nd_(2),this%ay_,this%by_,this%cy_,&
         this%py_,this%gamy_,this%bety_,this%xxy_,this%zetay_)
    CALL GPSplineInt_MatInvQ(this,this%nd_(3),this%az_,this%bz_,this%cz_,&
         this%pz_,this%gamz_,this%betz_,this%xxz_,this%zetaz_)
    CALL gupdate_to(this%ax_); CALL gupdate_to(this%px_); CALL gupdate_to(this%gamx_)
    CALL gupdate_to(this%betx_); CALL gupdate_to(this%xxx_)
    CALL gupdate_to(this%ay_); CALL gupdate_to(this%py_); CALL gupdate_to(this%gamy_)
    CALL gupdate_to(this%bety_); CALL gupdate_to(this%xxy_)
    CALL gupdate_to(this%az_); CALL gupdate_to(this%pz_); CALL gupdate_to(this%gamz_)
    CALL gupdate_to(this%betz_); CALL gupdate_to(this%xxz_)
  END SUBROUTINE GPSplineInt_Init


!-----------------------------------------------------------------
! Allocator
!-----------------------------------------------------------------
  SUBROUTINE GPSplineInt_DoAlloc(this)
    IMPLICIT NONE
    CLASS(GPSplineInt)        :: this
    INTEGER                   :: nzg

    CALL GPSplineInt_DoDealloc(this)
    nzg = this%gpcomm_%GetNumGhost()
    CALL galloc(this%esplfld_,this%ldims_(1),this%ldims_(2),1,this%ldims_(3)+2*nzg)
    CALL galloc(this%xrk_ ,this%maxint_)
    CALL galloc(this%yrk_ ,this%maxint_)
    CALL galloc(this%zrk_ ,this%maxint_)
    CALL galloc(this%wrkl_,9,this%maxint_)
    CALL galloc(this%ilg_ ,4,this%maxint_)
    CALL galloc(this%jlg_ ,4,this%maxint_)
    CALL galloc(this%klg_ ,4,this%maxint_)
    CALL galloc(this%tmptr_,this%ttot_)
    CALL galloc(this%tmpt2_,this%ttot_)
    ALLOCATE(this%bx_(this%nd_(1)),this%cx_(this%nd_(1)))
    ALLOCATE(this%by_(this%nd_(2)),this%cy_(this%nd_(2)))
    ALLOCATE(this%bz_(this%nd_(3)),this%cz_(this%nd_(3)))
    CALL galloc(this%ax_,this%nd_(1)); CALL galloc(this%px_  ,this%nd_(1)); CALL galloc(this%gamx_,this%nd_(1))
    CALL galloc(this%betx_,this%nd_(1)); CALL galloc(this%xxx_,this%nd_(1))
    CALL galloc(this%ay_,this%nd_(2)); CALL galloc(this%py_  ,this%nd_(2)); CALL galloc(this%gamy_,this%nd_(2))
    CALL galloc(this%bety_,this%nd_(2)); CALL galloc(this%xxy_,this%nd_(2))
    CALL galloc(this%az_,this%nd_(3)); CALL galloc(this%pz_  ,this%nd_(3)); CALL galloc(this%gamz_,this%nd_(3))
    CALL galloc(this%betz_,this%nd_(3)); CALL galloc(this%xxz_,this%nd_(3))
  END SUBROUTINE GPSplineInt_DoAlloc


!-----------------------------------------------------------------
! Deallocator
!-----------------------------------------------------------------
  SUBROUTINE GPSplineInt_DoDealloc(this)
    IMPLICIT NONE
    CLASS(GPSplineInt)        :: this
    CALL gfree(this%esplfld_)
    CALL gfree(this%xrk_) ; CALL gfree(this%yrk_) ; CALL gfree(this%zrk_)
    CALL gfree(this%wrkl_)
    CALL gfree(this%ilg_) ; CALL gfree(this%jlg_) ; CALL gfree(this%klg_)
    CALL gfree(this%tmptr_); CALL gfree(this%tmpt2_)
    IF ( ALLOCATED(this%bx_) ) DEALLOCATE(this%bx_,this%cx_,this%by_,this%cy_,this%bz_,this%cz_)
    CALL gfree(this%ax_); CALL gfree(this%px_); CALL gfree(this%gamx_); CALL gfree(this%betx_); CALL gfree(this%xxx_)
    CALL gfree(this%ay_); CALL gfree(this%py_); CALL gfree(this%gamy_); CALL gfree(this%bety_); CALL gfree(this%xxy_)
    CALL gfree(this%az_); CALL gfree(this%pz_); CALL gfree(this%gamz_); CALL gfree(this%betz_); CALL gfree(this%xxz_)
  END SUBROUTINE GPSplineInt_DoDealloc


!-----------------------------------------------------------------
!  METHOD     : SetDeriv
!  DESCRIPTION: Selects whether DoInterp3D returns the field
!               (ido=0) or its derivative (ido=1) along idir
!-----------------------------------------------------------------
  SUBROUTINE GPSplineInt_SetDeriv(this,idir,ido)
    IMPLICIT NONE
    CLASS(GPSplineInt)        :: this
    INTEGER      ,INTENT(IN)  :: idir,ido
    IF ( idir.LT.1 .OR. idir.GT.this%rank_ ) THEN
      WRITE(*,*) 'GPSplineInt::SetDeriv: invalid coordinate direction'
      STOP
    ENDIF
    IF ( ido.LT.0 .OR. ido.GT.1 ) THEN
      WRITE(*,*) 'GPSplineInt::SetDeriv: invalid derivative order'
      STOP
    ENDIF
    this%ider_(idir) = ido
  END SUBROUTINE GPSplineInt_SetDeriv


!-----------------------------------------------------------------
!  METHOD     : PartUpdate3D
!  DESCRIPTION: Computes, for the np particles at xp,yp,zp (grid
!               units, global in z), the indices of the 4 control
!               points of the stencil in each direction and the
!               fractional positions in the cell. x and y are
!               periodic; in z the indices refer to the extended
!               local grid, and every particle must lie between
!               zlo_ and xbnds_(3,2) (from two cells below to one
!               cell above the slab of this task).
!-----------------------------------------------------------------
  SUBROUTINE GPSplineInt_PartUpdate3D(this,xp,yp,zp,np)
    IMPLICIT NONE
    CLASS(GPSplineInt)                    :: this
    REAL(KIND=GP),INTENT(IN),DIMENSION(*) :: xp,yp,zp
    INTEGER      ,INTENT(IN)              :: np
    INTEGER                               :: j,kmax,kmin,nx,ny,nz,nbad

    IF ( np.LE.0 ) RETURN
    nx = this%ldims_(1)
    ny = this%ldims_(2)
    nz = this%ldims_(3)
    ! Bounds of the first stencil index in the extended grid (a
    ! safety clamp: particles in the valid z range never reach them)
    kmax = nz+2*this%gpcomm_%GetNumGhost()-3
    kmin = 1

    CALL gpsi_update_xy(np,xp,this%xbnds_(1,1),this%dxi_(1),nx,this%ilg_,this%xrk_)
    CALL gpsi_update_xy(np,yp,this%xbnds_(2,1),this%dxi_(2),ny,this%jlg_,this%yrk_)
    CALL gpsi_check_z(np,zp,this%zlo_,this%xbnds_(3,2),nbad)
    IF ( nbad.GT.0 ) THEN
      WRITE(*,*) myrank, ' GPSplineInt::PartUpdate3D: ',nbad,' particles out of the z-range'
      WRITE(*,*) myrank, ' GPSplineInt::PartUpdate3D: the particles moved more than one', &
                         ' grid cell out of the slab during the time step (dt too long)'
      DO j = 1, np
        IF ( .NOT.(zp(j).GE.this%zlo_.AND.zp(j).LT.this%xbnds_(3,2)) ) THEN
          WRITE(*,*) myrank, ' GPSplineInt::zbnd_0=',this%zlo_,';  zbnd_1=',this%xbnds_(3,2), 'zp=',zp(j)
        ENDIF
      ENDDO
      STOP
    ENDIF
    CALL gpsi_update_z(np,zp,this%xbnds_(3,1),this%dxi_(3),kmin,kmax,this%klg_,this%zrk_)
  END SUBROUTINE GPSplineInt_PartUpdate3D


!-----------------------------------------------------------------
! Stencil indices (periodic) and fractional position in one of
! the two periodic directions
!-----------------------------------------------------------------
  SUBROUTINE gpsi_update_xy(n,xp,xb,dxi,nx,ilg,xrk)
    IMPLICIT NONE
    INTEGER      ,INTENT(IN)    :: n,nx
    REAL(KIND=GP),INTENT(IN)    :: xp(n),xb,dxi
    INTEGER      ,INTENT(INOUT) :: ilg(4,n)
    REAL(KIND=GP),INTENT(INOUT) :: xrk(n)
    INTEGER                     :: j,i1
#if defined(GHOST_GPU)
!$omp target teams distribute parallel do if(target: gdev_active) private(i1)
#else
!$omp parallel do private(i1)
#endif
    DO j = 1, n
      i1       = (xp(j)-xb)*dxi
      xrk  (j) = (xp(j)-xb)*dxi - real(i1,kind=GP)
      ilg(2,j) = modulo(i1,nx) + 1
      ilg(3,j) = modulo(ilg(2,j),nx) + 1
      ilg(4,j) = modulo(ilg(3,j),nx) + 1
      ilg(1,j) = modulo(nx+ilg(2,j)-2,nx) + 1
    ENDDO
  END SUBROUTINE gpsi_update_xy


!-----------------------------------------------------------------
! Number of particles outside the ghost-extended slab
!-----------------------------------------------------------------
  SUBROUTINE gpsi_check_z(n,zp,zlo,zhi,nbad)
    IMPLICIT NONE
    INTEGER      ,INTENT(IN)    :: n
    REAL(KIND=GP),INTENT(IN)    :: zp(n),zlo,zhi
    INTEGER      ,INTENT(OUT)   :: nbad
    INTEGER                     :: j
    nbad = 0
#if defined(GHOST_GPU)
!$omp target teams distribute parallel do if(target: gdev_active) reduction(+:nbad)
#else
!$omp parallel do reduction(+:nbad)
#endif
    DO j = 1, n
      IF ( .NOT.(zp(j).GE.zlo.AND.zp(j).LT.zhi) ) nbad = nbad + 1
    ENDDO
  END SUBROUTINE gpsi_check_z


!-----------------------------------------------------------------
! Stencil indices (clamped to the extended local grid) and
! fractional position in z
!-----------------------------------------------------------------
  SUBROUTINE gpsi_update_z(n,zp,zb,dxi,kmin,kmax,klg,zrk)
    IMPLICIT NONE
    INTEGER      ,INTENT(IN)    :: n,kmin,kmax
    REAL(KIND=GP),INTENT(IN)    :: zp(n),zb,dxi
    INTEGER      ,INTENT(INOUT) :: klg(4,n)
    REAL(KIND=GP),INTENT(INOUT) :: zrk(n)
    INTEGER                     :: j,k1
#if defined(GHOST_GPU)
!$omp target teams distribute parallel do if(target: gdev_active) private(k1)
#else
!$omp parallel do private(k1)
#endif
    DO j = 1, n
      k1       = (zp(j)-zb)*dxi
      zrk  (j) = (zp(j)-zb)*dxi - real(k1,kind=GP)
      klg(1,j) = max(min(k1-1,kmax),kmin)
      klg(2,j) = klg(1,j) + 1
      klg(3,j) = klg(2,j) + 1
      klg(4,j) = klg(3,j) + 1
    ENDDO
  END SUBROUTINE gpsi_update_z


!=================================================================
! Interpolator (DoInterp3D) and helpers
!=================================================================

!-----------------------------------------------------------------
!  METHOD     : DoInterp3D
!  DESCRIPTION: Interpolates the field whose spline coefficients
!               are in esplfld_ at the np particles, into fp
!               (or its derivative, see SetDeriv). PartUpdate3D
!               and CompSpline3D/SetCoeffs3D must be called first.
!-----------------------------------------------------------------
  SUBROUTINE GPSplineInt_Interp3D(this,fp,np)
    IMPLICIT NONE
    CLASS(GPSplineInt)                       :: this
    INTEGER      ,INTENT   (IN)              :: np
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(*) :: fp
    REAL(KIND=GP)                            :: xsm,ysm,zsm
    INTEGER                                  :: nzg

    IF ( np.LE.0 ) RETURN
    nzg = this%gpcomm_%GetNumGhost()
    CALL gpsi_weights(np,this%xrk_,this%ider_(1),this%dxi_(1),0,this%wrkl_)
    CALL gpsi_weights(np,this%yrk_,this%ider_(2),this%dxi_(2),3,this%wrkl_)
    CALL gpsi_weights(np,this%zrk_,this%ider_(3),this%dxi_(3),6,this%wrkl_)
    xsm = 1.0_GP; IF ( this%ider_(1).EQ.1 ) xsm = 0.0_GP
    ysm = 1.0_GP; IF ( this%ider_(2).EQ.1 ) ysm = 0.0_GP
    zsm = 1.0_GP; IF ( this%ider_(3).EQ.1 ) zsm = 0.0_GP
    CALL gpsi_interp(np,this%ldims_(1),this%ldims_(2),this%ldims_(3)+2*nzg,xsm,ysm,zsm, &
                     this%wrkl_,this%ilg_,this%jlg_,this%klg_,this%esplfld_,fp)
  END SUBROUTINE GPSplineInt_Interp3D


!-----------------------------------------------------------------
! Basis weights of one direction from the fractional positions
! (slots ioff+1..ioff+3 of wrkl; the 4th weight is completed in
! the interpolation kernel as xsm-w1-w2-w3)
!-----------------------------------------------------------------
  SUBROUTINE gpsi_weights(n,xrk,ider,dxi,ioff,wrkl)
    IMPLICIT NONE
    INTEGER      ,INTENT(IN)    :: n,ider,ioff
    REAL(KIND=GP),INTENT(IN)    :: xrk(n),dxi
    REAL(KIND=GP),INTENT(INOUT) :: wrkl(9,n)
    REAL(KIND=GP)               :: xx,xxm,sixth,four,three,six,halfm,threeh,two
    INTEGER                     :: lag
    sixth  = 1.0_GP/6.0_GP
    four   = 4.0_GP
    three  = 3.0_GP
    six    = 6.0_GP
    halfm  = -0.5_GP
    threeh = 3.0_GP/2.0_GP
    two    = 2.0_GP
    IF ( ider.EQ.0 ) THEN
#if defined(GHOST_GPU)
!$omp target teams distribute parallel do if(target: gdev_active) private(xx,xxm)
#else
!$omp parallel do private(xx,xxm)
#endif
      DO lag=1,n
        xx = xrk(lag)
        xxm = (1.0_GP-xx)
        wrkl(ioff+1,lag) = sixth*xxm*xxm*xxm
        wrkl(ioff+2,lag) = sixth*(four+xx*xx*(three*xx-six))
        wrkl(ioff+3,lag) = sixth*(four+xxm*xxm*(three*xxm-six))
      ENDDO
    ELSE
#if defined(GHOST_GPU)
!$omp target teams distribute parallel do if(target: gdev_active) private(xx,xxm)
#else
!$omp parallel do private(xx,xxm)
#endif
      DO lag=1,n
        xx = xrk(lag)
        xxm = (1.0_GP-xx)
        wrkl(ioff+1,lag) = halfm*xxm*xxm*dxi
        wrkl(ioff+2,lag) = xx*(threeh*xx-two)*dxi
        wrkl(ioff+3,lag) = -xxm*(threeh*xxm-two)*dxi
      ENDDO
    ENDIF
  END SUBROUTINE gpsi_weights


!-----------------------------------------------------------------
! The 4x4x4 tensor-product sum over the control points, in the
! same term order as the original expression (same round-off)
!-----------------------------------------------------------------
  SUBROUTINE gpsi_interp(n,nx,ny,nez,xsm,ysm,zsm,wrkl,ilg,jlg,klg,e,fp)
    IMPLICIT NONE
    INTEGER      ,INTENT(IN)    :: n,nx,ny,nez
    REAL(KIND=GP),INTENT(IN)    :: xsm,ysm,zsm
    REAL(KIND=GP),INTENT(IN)    :: wrkl(9,n)
    INTEGER      ,INTENT(IN)    :: ilg(4,n),jlg(4,n),klg(4,n)
    REAL(KIND=GP),INTENT(IN)    :: e(nx,ny,nez)
    REAL(KIND=GP),INTENT(INOUT) :: fp(n)
    REAL(KIND=GP)               :: xx1,xx2,xx3,xx4,yy1,yy2,yy3,yy4,zz1,zz2,zz3,zz4
    INTEGER                     :: lag,i1,i2,i3,i4,j1,j2,j3,j4,k1,k2,k3,k4
#if defined(GHOST_GPU)
!$omp target teams distribute parallel do if(target: gdev_active) &
!$omp   private(xx1,xx2,xx3,xx4,yy1,yy2,yy3,yy4,zz1,zz2,zz3,zz4,i1,i2,i3,i4,j1,j2,j3,j4,k1,k2,k3,k4)
#else
!$omp parallel do private(xx1,xx2,xx3,xx4,yy1,yy2,yy3,yy4,zz1,zz2,zz3,zz4,i1,i2,i3,i4,j1,j2,j3,j4,k1,k2,k3,k4)
#endif
    DO lag=1,n
      xx1 = wrkl(1,lag)
      xx2 = wrkl(2,lag)
      xx3 = wrkl(3,lag)
      xx4 = xsm - xx1 - xx2 - xx3
      yy1 = wrkl(4,lag)
      yy2 = wrkl(5,lag)
      yy3 = wrkl(6,lag)
      yy4 = ysm - yy1 - yy2 - yy3
      zz1 = wrkl(7,lag)
      zz2 = wrkl(8,lag)
      zz3 = wrkl(9,lag)
      zz4 = zsm - zz1 - zz2 - zz3
      i1 = ilg(1,lag); i2 = ilg(2,lag); i3 = ilg(3,lag); i4 = ilg(4,lag)
      j1 = jlg(1,lag); j2 = jlg(2,lag); j3 = jlg(3,lag); j4 = jlg(4,lag)
      k1 = klg(1,lag); k2 = klg(2,lag); k3 = klg(3,lag); k4 = klg(4,lag)
      fp(lag) = e(i1,j1,k1)*xx1*yy1*zz1 &
      + e(i2,j1,k1)*xx2*yy1*zz1 &
      + e(i3,j1,k1)*xx3*yy1*zz1 &
      + e(i4,j1,k1)*xx4*yy1*zz1 &
      + e(i1,j2,k1)*xx1*yy2*zz1 &
      + e(i2,j2,k1)*xx2*yy2*zz1 &
      + e(i3,j2,k1)*xx3*yy2*zz1 &
      + e(i4,j2,k1)*xx4*yy2*zz1 &
      + e(i1,j3,k1)*xx1*yy3*zz1 &
      + e(i2,j3,k1)*xx2*yy3*zz1 &
      + e(i3,j3,k1)*xx3*yy3*zz1 &
      + e(i4,j3,k1)*xx4*yy3*zz1 &
      + e(i1,j4,k1)*xx1*yy4*zz1 &
      + e(i2,j4,k1)*xx2*yy4*zz1 &
      + e(i3,j4,k1)*xx3*yy4*zz1 &
      + e(i4,j4,k1)*xx4*yy4*zz1 &
      + e(i1,j1,k2)*xx1*yy1*zz2 &
      + e(i2,j1,k2)*xx2*yy1*zz2 &
      + e(i3,j1,k2)*xx3*yy1*zz2 &
      + e(i4,j1,k2)*xx4*yy1*zz2 &
      + e(i1,j2,k2)*xx1*yy2*zz2 &
      + e(i2,j2,k2)*xx2*yy2*zz2 &
      + e(i3,j2,k2)*xx3*yy2*zz2 &
      + e(i4,j2,k2)*xx4*yy2*zz2 &
      + e(i1,j3,k2)*xx1*yy3*zz2 &
      + e(i2,j3,k2)*xx2*yy3*zz2 &
      + e(i3,j3,k2)*xx3*yy3*zz2 &
      + e(i4,j3,k2)*xx4*yy3*zz2 &
      + e(i1,j4,k2)*xx1*yy4*zz2 &
      + e(i2,j4,k2)*xx2*yy4*zz2 &
      + e(i3,j4,k2)*xx3*yy4*zz2 &
      + e(i4,j4,k2)*xx4*yy4*zz2 &
      + e(i1,j1,k3)*xx1*yy1*zz3 &
      + e(i2,j1,k3)*xx2*yy1*zz3 &
      + e(i3,j1,k3)*xx3*yy1*zz3 &
      + e(i4,j1,k3)*xx4*yy1*zz3 &
      + e(i1,j2,k3)*xx1*yy2*zz3 &
      + e(i2,j2,k3)*xx2*yy2*zz3 &
      + e(i3,j2,k3)*xx3*yy2*zz3 &
      + e(i4,j2,k3)*xx4*yy2*zz3 &
      + e(i1,j3,k3)*xx1*yy3*zz3 &
      + e(i2,j3,k3)*xx2*yy3*zz3 &
      + e(i3,j3,k3)*xx3*yy3*zz3 &
      + e(i4,j3,k3)*xx4*yy3*zz3 &
      + e(i1,j4,k3)*xx1*yy4*zz3 &
      + e(i2,j4,k3)*xx2*yy4*zz3 &
      + e(i3,j4,k3)*xx3*yy4*zz3 &
      + e(i4,j4,k3)*xx4*yy4*zz3 &
      + e(i1,j1,k4)*xx1*yy1*zz4 &
      + e(i2,j1,k4)*xx2*yy1*zz4 &
      + e(i3,j1,k4)*xx3*yy1*zz4 &
      + e(i4,j1,k4)*xx4*yy1*zz4 &
      + e(i1,j2,k4)*xx1*yy2*zz4 &
      + e(i2,j2,k4)*xx2*yy2*zz4 &
      + e(i3,j2,k4)*xx3*yy2*zz4 &
      + e(i4,j2,k4)*xx4*yy2*zz4 &
      + e(i1,j3,k4)*xx1*yy3*zz4 &
      + e(i2,j3,k4)*xx2*yy3*zz4 &
      + e(i3,j3,k4)*xx3*yy3*zz4 &
      + e(i4,j3,k4)*xx4*yy3*zz4 &
      + e(i1,j4,k4)*xx1*yy4*zz4 &
      + e(i2,j4,k4)*xx2*yy4*zz4 &
      + e(i3,j4,k4)*xx3*yy4*zz4 &
      + e(i4,j4,k4)*xx4*yy4*zz4
    ENDDO
  END SUBROUTINE gpsi_interp


!=================================================================
! Computation of coefficients and resizing, spline coefficients
! by the periodic tridiagonal solves, and helpers 
!=================================================================

!-----------------------------------------------------------------
!  METHOD     : SetCoeffs3D
!  DESCRIPTION: Fills the extended coefficient field from a
!               field that already holds the spline coefficients
!               (ghost-plane exchange only)
!-----------------------------------------------------------------
  SUBROUTINE GPSplineInt_SetCoeffs3D(this,field)
    IMPLICIT NONE
    CLASS(GPSplineInt)                                :: this
    REAL(KIND=GP),INTENT(IN),DIMENSION(this%ntot_)    :: field
    CALL GTStart(this%hdataex_)
    CALL this%gpcomm_%SlabDataExchangeSF(this%esplfld_,field)
    CALL GTAcc(this%hdataex_)
  END SUBROUTINE GPSplineInt_SetCoeffs3D


!-----------------------------------------------------------------
!  METHOD     : ResizeArrays
!  DESCRIPTION: Resizes the particle-sized arrays (no data kept)
!-----------------------------------------------------------------
  SUBROUTINE GPSplineInt_ResizeArrays(this,newmparts,onlyinc)
    IMPLICIT NONE
    CLASS(GPSplineInt),INTENT(INOUT)   :: this
    INTEGER           ,INTENT(IN)      :: newmparts
    LOGICAL           ,INTENT(IN)      :: onlyinc
    INTEGER                            :: n
    n = SIZE(this%wrkl_, 2)
    IF ((n.lt.newmparts).OR.((n.gt.newmparts).AND..NOT.onlyinc)) THEN
      CALL gresize(this%ilg_ ,4,newmparts,.false.)
      CALL gresize(this%jlg_ ,4,newmparts,.false.)
      CALL gresize(this%klg_ ,4,newmparts,.false.)
      CALL gresize(this%wrkl_,9,newmparts,.false.)
      CALL gresize(this%xrk_ ,newmparts,.false.)
      CALL gresize(this%yrk_ ,newmparts,.false.)
      CALL gresize(this%zrk_ ,newmparts,.false.)
      this%maxint_ = newmparts
    END IF
  END SUBROUTINE GPSplineInt_ResizeArrays


!-----------------------------------------------------------------
!  METHOD     : MatInvQ
!  DESCRIPTION: Computes quantities for matrix inversion
!  ARGUMENTS  : 
!    this     : 'this' class instance
!-----------------------------------------------------------------
  SUBROUTINE GPSplineInt_MatInvQ(this,n,a,b,c,p,gam,bet,xx,zeta)
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
!  METHOD     : CompSpline3D
!  DESCRIPTION: Computes the spline coefficients of the slab
!               field(nx,ny,kl) (periodic cubic spline) and fills
!               the extended coefficient field esplfld_. The
!               periodic tridiagonal systems are solved pencil by
!               pencil, in x, then y, then z; every pencil is
!               independent and the arithmetic in each pencil is
!               that of the original sweeps, so host and device
!               builds give the same results while their loops are
!               organized differently:
!               - Device (GHOST_GPU): one thread per pencil, running
!                 values in registers. A recurrence must never run
!                 along the contiguous index (each step of every
!                 thread would touch a different cache line, 16x
!                 slower), so x is solved on the slab transposed to
!                 (ny,nx,kl) with the register-tile transpose
!                 gpsi_xytr, y on (nx,ny,kl) along the second index
!                 (gpsi_solve2), and z on the z-complete layout
!                 (il,ny,nz) along the third index (gpsi_solve3);
!                 that layout keeps x contiguous, so the MPI
!                 transpose is a pure block copy (gpartcomm).
!               - Host: the opposite. x is solved along the
!                 contiguous index (gpsi_solve1, host only, streams
!                 through the caches); y and z advance plane by
!                 plane with the contiguous index innermost, so the
!                 sweeps vectorize and stream (host branches of
!                 gpsi_solve2/3). No local transposes are needed.
!               tmp2 is a temporary of the size of the field.
!-----------------------------------------------------------------
  SUBROUTINE GPSplineInt_CompSpline3D(this,field,tmp2)
    USE mpivars
    IMPLICIT NONE
    CLASS(GPSplineInt)                                :: this
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(this%ntot_) :: field
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(this%ntot_) :: tmp2
    INTEGER                                           :: nx,ny,nz,il

    nx = this%ldims_(1)
    ny = this%ldims_(2)
    nz = this%ldims_(3)
    il = this%odims_(3)

    ! Solves in x and y on the slab. On the device a recurrence must
    ! not run along the fastest index (each step of every thread would
    ! touch a different cache line, 16 times slower than a solve along
    ! a slower index), so the x solve is done on the slab transposed to
    ! (ny,nx,kl). On the host the opposite holds: the solve along the
    ! contiguous index streams through the caches and the transposes
    ! would only add memory traffic (25% slower at 256^3 with 8 MPI
    ! tasks), so the host solves x in place. The order of the
    ! operations in each pencil is that of the original sweeps in
    ! both cases, and the results are identical
#if defined(GHOST_GPU)
    CALL GTStart(this%htransp_)
    CALL gpsi_xytr(nx,ny,nz,field,tmp2)                          ! tmp2(ny,nx,kl)
    CALL GTAcc(this%htransp_)
    CALL gpsi_solve2(ny,nx,nz,this%ax_,this%betx_,this%gamx_,this%px_,this%xxx_,this%zetax_,tmp2,field)
    CALL GTStart(this%htransp_)
    CALL gpsi_xytr(ny,nx,nz,field,tmp2)                          ! tmp2(nx,ny,kl)
    CALL GTAcc(this%htransp_)
#else
    CALL gpsi_solve1(nx,ny,nz,this%ax_,this%betx_,this%gamx_,this%px_,this%xxx_,this%zetax_,field,tmp2)
#endif
    CALL gpsi_solve2(nx,ny,nz,this%ay_,this%bety_,this%gamy_,this%py_,this%xxy_,this%zetay_,tmp2,field)

    ! Solve in z on the z-complete layout (il,ny,nz), along the third
    ! index. With a single task that layout is the slab itself
    IF ( nprocs .EQ. 1 ) THEN
      CALL gpsi_solve3(nx,ny,nz,this%az_,this%betz_,this%gamz_,this%pz_,this%xxz_,this%zetaz_,field,tmp2)
      CALL GTStart(this%hdataex_)
      CALL this%gpcomm_%SlabDataExchangeSF(this%esplfld_,tmp2)
      CALL GTAcc(this%hdataex_)
      RETURN
    ENDIF
    CALL GTStart(this%htransp_)
    CALL this%gpcomm_%GTranspose(this%tmptr_,field)
    CALL GTAcc(this%htransp_)
    CALL gpsi_solve3(il,ny,this%nd_(3),this%az_,this%betz_,this%gamz_,this%pz_,this%xxz_,this%zetaz_, &
                     this%tmptr_,this%tmpt2_)
    CALL GTStart(this%htransp_)
    CALL this%gpcomm_%GITranspose(field,this%tmpt2_)
    CALL GTAcc(this%htransp_)

    ! Ghost planes of the coefficients
    CALL GTStart(this%hdataex_)
    CALL this%gpcomm_%SlabDataExchangeSF(this%esplfld_,field)
    CALL GTAcc(this%hdataex_)
  END SUBROUTINE GPSplineInt_CompSpline3D


!-----------------------------------------------------------------
! Local transpose of the two fastest indices: t(j,i,k) = f(i,j,k).
! Each thread transposes a TS x TS block through registers, so that
! it reads and writes TS consecutive elements (32 bytes) at a time;
! consecutive threads take consecutive j blocks, which makes the
! writes of a wavefront contiguous. Measured on an MI210 at 256^3
! this is 5.6 times faster than one element per thread (the partial
! line writes of the latter multiply the traffic by 16). Used only
! in GHOST_GPU builds (the host solves x in place); in those builds
! the host fallback (gdev_active unset) transposes the blocks of one
! plane per OpenMP thread.
!-----------------------------------------------------------------
  SUBROUTINE gpsi_xytr(n1,n2,n3,f,t)
    IMPLICIT NONE
    INTEGER      ,INTENT(IN)    :: n1,n2,n3
    REAL(KIND=GP),INTENT(IN)    :: f(n1,n2,n3)
    REAL(KIND=GP),INTENT(INOUT) :: t(n2,n1,n3)
    INTEGER      ,PARAMETER     :: TS = 8
    REAL(KIND=GP)               :: b(TS,TS)
    INTEGER                     :: i,j,k,ib,jb,ii,jj,nb1,nb2

    nb1 = (n1+TS-1)/TS
    nb2 = (n2+TS-1)/TS
    IF ( MOD(n1,TS).EQ.0 .AND. MOD(n2,TS).EQ.0 ) THEN
      ! Whole blocks only: no bounds tests in the inner loops, so that
      ! the compiler issues vector loads and stores of a block row
#if defined(GHOST_GPU)
!$omp target teams distribute parallel do collapse(3) if(target: gdev_active) private(b,ii,jj)
#else
!$omp parallel do collapse(2) private(b,ii,jj,jb)
#endif
      DO k = 1,n3
        DO ib = 0,nb1-1
          DO jb = 0,nb2-1
            DO jj = 1,TS
              DO ii = 1,TS
                b(ii,jj) = f(ib*TS+ii,jb*TS+jj,k)
              ENDDO
            ENDDO
            DO ii = 1,TS
              DO jj = 1,TS
                t(jb*TS+jj,ib*TS+ii,k) = b(ii,jj)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ELSE
#if defined(GHOST_GPU)
!$omp target teams distribute parallel do collapse(3) if(target: gdev_active) private(b,i,j,ii,jj)
#else
!$omp parallel do collapse(2) private(b,i,j,ii,jj,jb)
#endif
      DO k = 1,n3
        DO ib = 0,nb1-1
          DO jb = 0,nb2-1
            DO jj = 1,TS
              j = jb*TS+jj
              DO ii = 1,TS
                i = ib*TS+ii
                IF ( i.LE.n1 .AND. j.LE.n2 ) b(ii,jj) = f(i,j,k)
              ENDDO
            ENDDO
            DO ii = 1,TS
              i = ib*TS+ii
              DO jj = 1,TS
                j = jb*TS+jj
                IF ( i.LE.n1 .AND. j.LE.n2 ) t(j,i,k) = b(ii,jj)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDIF
  END SUBROUTINE gpsi_xytr


!-----------------------------------------------------------------
! Periodic tridiagonal solve along the second index of f(n1,n2,n3),
! result in t (a,bet,gam,p,xx,zeta: factorization of MatInvQ).
! - Device: one pencil (i,k) per thread; consecutive threads take
!   consecutive i, so the loads and stores of a wavefront are
!   coalesced. The running values of the recurrences are carried in
!   registers (tp, tn), so no step waits for the value the thread
!   has just stored. Do not move the recurrence to the first index.
! - Host: the recurrence advances plane by plane; each sweep runs
!   over the contiguous index i innermost, so it vectorizes and
!   streams through memory (one OpenMP thread per k plane). The
!   accumulators of the last row are kept in the vector tnv(i).
! Both forms do the same operations on each pencil, in the same
! order, and give identical results.
!-----------------------------------------------------------------
  SUBROUTINE gpsi_solve2(n1,n2,n3,a,bet,gam,p,xx,zeta,f,t)
    IMPLICIT NONE
    INTEGER      ,INTENT(IN)    :: n1,n2,n3
    REAL(KIND=GP),INTENT(IN)    :: a(n2),bet(n2),gam(n2),p(n2),xx(n2),zeta
    REAL(KIND=GP),INTENT(IN)    :: f(n1,n2,n3)
    REAL(KIND=GP),INTENT(INOUT) :: t(n1,n2,n3)
    REAL(KIND=GP)               :: tp,tn
    REAL(KIND=GP),ALLOCATABLE   :: tnv(:)
    INTEGER                     :: i,j,k
#if defined(GHOST_GPU)
    ! Device: one pencil per thread, (k,i) collapsed with i fastest
!$omp target teams distribute parallel do collapse(2) if(target: gdev_active) private(j,tp,tn)
    DO k = 1,n3
      DO i = 1,n1
        tp = f(i,1,k)*bet(1)
        t(i,1,k) = tp
        tn = f(i,n2,k)
        DO j = 2,n2-2
          tp = ( f(i,j,k) - a(j)*tp )*bet(j)
          t(i,j,k) = tp
        ENDDO
        DO j = 2,n2-2
          tn = tn - xx(j-1)*t(i,j-1,k)
        ENDDO
        tp = (f(i,n2-1,k) - a(n2-1)*t(i,n2-2,k))*bet(n2-1)
        tn = tn - xx(n2-2)*tp
        tn = (tn - tp*zeta)*bet(n2)
        tp = tp - gam(n2)*tn
        t(i,n2,k)   = tn
        t(i,n2-1,k) = tp
        DO j = n2-2,1,-1
          tp = t(i,j,k) - gam(j+1)*tp - p(j)*tn
          t(i,j,k) = tp
        ENDDO
      ENDDO
    ENDDO
#else
    ! Host: the recurrence runs over planes, with the contiguous index
    ! innermost (vectorized, streaming), and the same operations per
    ! pencil as above; tnv holds the last-row accumulators of one plane
    ALLOCATE(tnv(n1))
!$omp parallel do private(i,j,tp,tn,tnv)
    DO k = 1,n3
      DO i = 1,n1
        t(i,1,k) = f(i,1,k)*bet(1)
        tnv(i)   = f(i,n2,k)
      ENDDO
      DO j = 2,n2-2
        DO i = 1,n1
          t(i,j,k) = ( f(i,j,k) - a(j)*t(i,j-1,k) )*bet(j)
          tnv(i)   = tnv(i) - xx(j-1)*t(i,j-1,k)
        ENDDO
      ENDDO
      DO i = 1,n1
        tp = (f(i,n2-1,k) - a(n2-1)*t(i,n2-2,k))*bet(n2-1)
        tn = tnv(i) - xx(n2-2)*tp
        tn = (tn - tp*zeta)*bet(n2)
        tp = tp - gam(n2)*tn
        t(i,n2,k)   = tn
        t(i,n2-1,k) = tp
      ENDDO
      DO j = n2-2,1,-1
        DO i = 1,n1
          t(i,j,k) = t(i,j,k) - gam(j+1)*t(i,j+1,k) - p(j)*t(i,n2,k)
        ENDDO
      ENDDO
    ENDDO
    DEALLOCATE(tnv)
#endif
  END SUBROUTINE gpsi_solve2


!-----------------------------------------------------------------
! Periodic tridiagonal solve along the third index of f(n1,n2,n3),
! result in t.
! - Device: one pencil (i,j) per thread with the running values in
!   registers; consecutive threads take consecutive i (coalesced
!   accesses). This is the fastest of the three device solves: a
!   whole (i,j) plane of threads is in flight at each step.
! - Host: the recurrence advances plane by plane in k; each sweep
!   over a plane (i,j) has the contiguous index innermost and is
!   split over the OpenMP threads by rows j, so it vectorizes and
!   streams. The last-plane accumulators are kept in tnv(i,j).
! Same operations per pencil in both forms, identical results.
!-----------------------------------------------------------------
  SUBROUTINE gpsi_solve3(n1,n2,n3,a,bet,gam,p,xx,zeta,f,t)
    IMPLICIT NONE
    INTEGER      ,INTENT(IN)    :: n1,n2,n3
    REAL(KIND=GP),INTENT(IN)    :: a(n3),bet(n3),gam(n3),p(n3),xx(n3),zeta
    REAL(KIND=GP),INTENT(IN)    :: f(n1,n2,n3)
    REAL(KIND=GP),INTENT(INOUT) :: t(n1,n2,n3)
    REAL(KIND=GP)               :: tp,tn
    REAL(KIND=GP),ALLOCATABLE   :: tnv(:,:)
    INTEGER                     :: i,j,k
#if defined(GHOST_GPU)
    ! Device: one pencil per thread, (j,i) collapsed with i fastest
!$omp target teams distribute parallel do collapse(2) if(target: gdev_active) private(k,tp,tn)
    DO j = 1,n2
      DO i = 1,n1
        tp = f(i,j,1)*bet(1)
        t(i,j,1) = tp
        tn = f(i,j,n3)
        DO k = 2,n3-2
          tp = ( f(i,j,k) - a(k)*tp )*bet(k)
          t(i,j,k) = tp
        ENDDO
        DO k = 2,n3-2
          tn = tn - xx(k-1)*t(i,j,k-1)
        ENDDO
        tp = (f(i,j,n3-1) - a(n3-1)*t(i,j,n3-2))*bet(n3-1)
        tn = tn - xx(n3-2)*tp
        tn = (tn - tp*zeta)*bet(n3)
        tp = tp - gam(n3)*tn
        t(i,j,n3)   = tn
        t(i,j,n3-1) = tp
        DO k = n3-2,1,-1
          tp = t(i,j,k) - gam(k+1)*tp - p(k)*tn
          t(i,j,k) = tp
        ENDDO
      ENDDO
    ENDDO
#else
    ! Host: the recurrence runs over planes (i,j), each sweep streams
    ! through a contiguous plane; same operations per pencil as above
    ALLOCATE(tnv(n1,n2))
!$omp parallel do private(i)
    DO j = 1,n2
      DO i = 1,n1
        t(i,j,1) = f(i,j,1)*bet(1)
        tnv(i,j) = f(i,j,n3)
      ENDDO
    ENDDO
    DO k = 2,n3-2
!$omp parallel do private(i)
      DO j = 1,n2
        DO i = 1,n1
          t(i,j,k)  = ( f(i,j,k) - a(k)*t(i,j,k-1) )*bet(k)
          tnv(i,j)  = tnv(i,j) - xx(k-1)*t(i,j,k-1)
        ENDDO
      ENDDO
    ENDDO
!$omp parallel do private(i,tp,tn)
    DO j = 1,n2
      DO i = 1,n1
        tp = (f(i,j,n3-1) - a(n3-1)*t(i,j,n3-2))*bet(n3-1)
        tn = tnv(i,j) - xx(n3-2)*tp
        tn = (tn - tp*zeta)*bet(n3)
        tp = tp - gam(n3)*tn
        t(i,j,n3)   = tn
        t(i,j,n3-1) = tp
      ENDDO
    ENDDO
    DO k = n3-2,1,-1
!$omp parallel do private(i)
      DO j = 1,n2
        DO i = 1,n1
          t(i,j,k) = t(i,j,k) - gam(k+1)*t(i,j,k+1) - p(k)*t(i,j,n3)
        ENDDO
      ENDDO
    ENDDO
    DEALLOCATE(tnv)
#endif
  END SUBROUTINE gpsi_solve3


#if !defined(GHOST_GPU)
!-----------------------------------------------------------------
! Periodic tridiagonal solve along the first (contiguous) index of
! f(n1,n2,n3), result in t. Host builds only: each OpenMP thread
! walks one pencil (j,k) along contiguous memory with the running
! values in registers (tp, tn), the cache-friendly form on a CPU
! (sequential access, hardware prefetch). On the device this access
! pattern is 16 times slower than the solves along the other indices
! (one cache line per thread and step), so GHOST_GPU builds solve x
! on the slab transposed by gpsi_xytr with gpsi_solve2 instead.
!-----------------------------------------------------------------
  SUBROUTINE gpsi_solve1(n1,n2,n3,a,bet,gam,p,xx,zeta,f,t)
    IMPLICIT NONE
    INTEGER      ,INTENT(IN)    :: n1,n2,n3
    REAL(KIND=GP),INTENT(IN)    :: a(n1),bet(n1),gam(n1),p(n1),xx(n1),zeta
    REAL(KIND=GP),INTENT(IN)    :: f(n1,n2,n3)
    REAL(KIND=GP),INTENT(INOUT) :: t(n1,n2,n3)
    REAL(KIND=GP)               :: tp,tn
    INTEGER                     :: i,j,k
!$omp parallel do collapse(2) private(i,tp,tn)
    DO k = 1,n3
      DO j = 1,n2
        tp = f(1,j,k)*bet(1)
        t(1,j,k) = tp
        tn = f(n1,j,k)
        DO i = 2,n1-2
          tp = ( f(i,j,k) - a(i)*tp )*bet(i)
          t(i,j,k) = tp
        ENDDO
        DO i = 2,n1-2
          tn = tn - xx(i-1)*t(i-1,j,k)
        ENDDO
        tp = (f(n1-1,j,k) - a(n1-1)*t(n1-2,j,k))*bet(n1-1)
        tn = tn - xx(n1-2)*tp
        tn = (tn - tp*zeta)*bet(n1)
        tp = tp - gam(n1)*tn
        t(n1,j,k)   = tn
        t(n1-1,j,k) = tp
        DO i = n1-2,1,-1
          tp = t(i,j,k) - gam(i+1)*tp - p(i)*tn
          t(i,j,k) = tp
        ENDDO
      ENDDO
    ENDDO
  END SUBROUTINE gpsi_solve1
#endif

END MODULE class_GPSplineInt
