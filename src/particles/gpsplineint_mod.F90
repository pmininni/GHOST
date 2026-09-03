!=================================================================
! GPSplineInt: cubic spline interpolation of a slab-decomposed
! field at the particle positions.
!
! The interpolation reads the spline coefficients on the extended
! grid esplfld_(nx,ny,nzl+2*nzghost), with nzghost planes of the
! neighbor slabs below and above the local slab, filled by the
! ghost-plane exchange of GPartComm. The particle-sized arrays
! (cell indices ilg_,jlg_,klg_, fractional positions xrk_,yrk_,
! zrk_ and weights wrkl_) and the coefficient field have device
! copies in offload builds; PartUpdate3D and DoInterp3D run their
! kernels on the device while gdev_active is set. The kernels are
! module procedures taking explicit-shape arrays.
!
! Two ways of computing the spline coefficients are available:
!   CompSpline3D: solves the periodic tridiagonal systems in the
!                 three directions on the host (transposes for z)
!                 and fills the extended field;
!   SetCoeffs3D : takes a field that already holds the spline
!                 coefficients (computed in Fourier space by the
!                 caller) and only fills the extended field.
!=================================================================
MODULE class_GPSplineInt
      USE mpivars
      USE fprecision
      USE class_GPartComm
      USE gtimer
      USE gmem
      USE gdevice, ONLY: gdev_active
      IMPLICIT NONE

      PRIVATE
      TYPE, PUBLIC :: GPSplineInt
        PRIVATE
        ! Device-resident arrays
        REAL(KIND=GP),ALLOCATABLE,DIMENSION(:,:,:) :: esplfld_
        REAL(KIND=GP),ALLOCATABLE,DIMENSION    (:) :: xrk_,yrk_,zrk_
        REAL(KIND=GP),ALLOCATABLE,DIMENSION  (:,:) :: wrkl_
        INTEGER      ,ALLOCATABLE,DIMENSION  (:,:) :: ilg_,jlg_,klg_
        REAL(KIND=GP),ALLOCATABLE,DIMENSION    (:) :: tmptr_,tmpt2_  ! z-complete layout
        ! Factorization of the tridiagonal systems (computed on the
        ! host; the five arrays the solves use have device copies)
        REAL(KIND=GP),ALLOCATABLE,DIMENSION    (:) :: ax_,bx_,betx_,cx_,gamx_,px_,xxx_
        REAL(KIND=GP),ALLOCATABLE,DIMENSION    (:) :: ay_,by_,bety_,cy_,gamy_,py_,xxy_
        REAL(KIND=GP),ALLOCATABLE,DIMENSION    (:) :: az_,bz_,betz_,cz_,gamz_,pz_,xxz_
        REAL(KIND=GP)                              :: dxi_(3),xbnds_(3,2),zetax_,zetay_,zetaz_
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

!-----------------------------------------------------------------
!  METHOD     : GPSplineInt_ctor
!  DESCRIPTION: Constructor. Stores the grid and slab bounds
!               (xbnds widened by nzghost planes in z), the
!               transposed bounds used by CompSpline3D, the
!               particle buffer size and the timer handles.
!-----------------------------------------------------------------
  SUBROUTINE GPSplineInt_ctor(this,rank,nd,ibnds,xbnds,obnds,nzghost,maxpart,gpcomm, &
                              hdataex,htransp)
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
    this%xbnds_(3,1)  = this%xbnds_(3,1)-real(nzghost,kind=GP)
    this%xbnds_(3,2)  = this%xbnds_(3,2)+real(nzghost,kind=GP)
    CALL GPSplineInt_Init(this)
  END SUBROUTINE GPSplineInt_ctor


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
!               local grid, and every particle must lie in the
!               ghost-extended slab of this task.
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
    kmax = nz+this%gpcomm_%GetNumGhost()-1
    kmin = this%gpcomm_%GetNumGhost()-1

    CALL gpsi_update_xy(np,xp,this%xbnds_(1,1),this%dxi_(1),nx,this%ilg_,this%xrk_)
    CALL gpsi_update_xy(np,yp,this%xbnds_(2,1),this%dxi_(2),ny,this%jlg_,this%yrk_)
    CALL gpsi_check_z(np,zp,this%xbnds_(3,1),this%xbnds_(3,2),nbad)
    IF ( nbad.GT.0 ) THEN
      DO j = 1, np
        IF ( .NOT.(zp(j).GE.this%xbnds_(3,1).AND.zp(j).LT.this%xbnds_(3,2)) ) THEN
          WRITE(*,*) myrank, ' GPSplineInt::PartUpdate3D: Invalid particle z-range'
          WRITE(*,*) myrank, ' GPSplineInt::zbnd_0=',this%xbnds_(3,1),';  zbnd_1=',this%xbnds_(3,2), 'zp=',zp(j)
          STOP
        ENDIF
      ENDDO
      WRITE(*,*) myrank, ' GPSplineInt::PartUpdate3D: ',nbad,' particles out of the z-range'
      STOP
    ENDIF
    CALL gpsi_update_z(np,zp,this%xbnds_(3,1),this%dxi_(3),kmin,kmax,this%klg_,this%zrk_)
  END SUBROUTINE GPSplineInt_PartUpdate3D


! Stencil indices (periodic) and fractional position in one of
! the two periodic directions
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


! Number of particles outside the ghost-extended slab
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


! Stencil indices (clamped to the extended local grid) and
! fractional position in z
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


! Basis weights of one direction from the fractional positions
! (slots ioff+1..ioff+3 of wrkl; the 4th weight is completed in
! the interpolation kernel as xsm-w1-w2-w3)
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


! The 4x4x4 tensor-product sum over the control points, in the
! same term order as the original expression (same round-off)
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


!=================================================================
! Spline coefficients by the periodic tridiagonal solves
!=================================================================

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
!  METHOD     : CompSpline3D
!  DESCRIPTION: Computes the spline coefficients of the slab
!               field(nx,ny,kl) (periodic cubic spline) and fills
!               the extended coefficient field esplfld_. The
!               periodic tridiagonal systems are solved pencil by
!               pencil: in x (field -> tmp2), in y (tmp2 -> field),
!               then the field is transposed to the z-complete
!               layout, solved in z (tmptr_ -> tmpt2_) and
!               transposed back. Every pencil is independent, so
!               the kernels run one thread per pencil on the device
!               (one OpenMP thread per pencil on the host) with the
!               same arithmetic as the original sweeps. tmp2 is a
!               temporary of the size of the field.
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

    ! Solves in x and y on the slab
    CALL gpsi_solve1(nx,ny,nz,this%ax_,this%betx_,this%gamx_,this%px_,this%xxx_,this%zetax_,field,tmp2)
    CALL gpsi_solve2(nx,ny,nz,this%ay_,this%bety_,this%gamy_,this%py_,this%xxy_,this%zetay_,tmp2,field)

    ! Solve in z on the z-complete layout
    CALL GTStart(this%htransp_)
    CALL this%gpcomm_%GTranspose(this%tmptr_,field)
    CALL GTAcc(this%htransp_)
    CALL gpsi_solve1(this%nd_(3),ny,il,this%az_,this%betz_,this%gamz_,this%pz_,this%xxz_,this%zetaz_, &
                     this%tmptr_,this%tmpt2_)
    CALL GTStart(this%htransp_)
    CALL this%gpcomm_%GITranspose(field,this%tmpt2_)
    CALL GTAcc(this%htransp_)

    ! Ghost planes of the coefficients
    CALL GTStart(this%hdataex_)
    CALL this%gpcomm_%SlabDataExchangeSF(this%esplfld_,field)
    CALL GTAcc(this%hdataex_)
  END SUBROUTINE GPSplineInt_CompSpline3D


! Periodic tridiagonal solve along the first index of f(n1,n2,n3),
! one pencil (j,k) per thread, result in t (a,bet,gam,p,xx,zeta:
! factorization of MatInvQ). The running values of the recurrences
! are carried in registers (tp, tn) so that no step waits for the
! value the thread has just stored; the arithmetic is the same as
! the original sweeps.
  SUBROUTINE gpsi_solve1(n1,n2,n3,a,bet,gam,p,xx,zeta,f,t)
    IMPLICIT NONE
    INTEGER      ,INTENT(IN)    :: n1,n2,n3
    REAL(KIND=GP),INTENT(IN)    :: a(n1),bet(n1),gam(n1),p(n1),xx(n1),zeta
    REAL(KIND=GP),INTENT(IN)    :: f(n1,n2,n3)
    REAL(KIND=GP),INTENT(INOUT) :: t(n1,n2,n3)
    REAL(KIND=GP)               :: tp,tn
    INTEGER                     :: i,j,k
#if defined(GHOST_GPU)
!$omp target teams distribute parallel do collapse(2) if(target: gdev_active) private(i,tp,tn)
#else
!$omp parallel do collapse(2) private(i,tp,tn)
#endif
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


! Periodic tridiagonal solve along the second index of f(n1,n2,n3),
! one pencil (i,k) per thread (consecutive threads read consecutive
! i), result in t
  SUBROUTINE gpsi_solve2(n1,n2,n3,a,bet,gam,p,xx,zeta,f,t)
    IMPLICIT NONE
    INTEGER      ,INTENT(IN)    :: n1,n2,n3
    REAL(KIND=GP),INTENT(IN)    :: a(n2),bet(n2),gam(n2),p(n2),xx(n2),zeta
    REAL(KIND=GP),INTENT(IN)    :: f(n1,n2,n3)
    REAL(KIND=GP),INTENT(INOUT) :: t(n1,n2,n3)
    REAL(KIND=GP)               :: tp,tn
    INTEGER                     :: i,j,k
#if defined(GHOST_GPU)
!$omp target teams distribute parallel do collapse(2) if(target: gdev_active) private(j,tp,tn)
#else
!$omp parallel do collapse(2) private(j,tp,tn)
#endif
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
  END SUBROUTINE gpsi_solve2

END MODULE class_GPSplineInt
