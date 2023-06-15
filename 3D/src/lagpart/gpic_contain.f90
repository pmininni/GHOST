!=================================================================
! GPIC SUBROUTINES
!=================================================================

  SUBROUTINE GPIC_ctor(this,comm,ppc,inittype,iinterp,intorder,iexchtyp, &
                        iouttyp,bcoll,csize,nstrip,intacc,wrtunit)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  Explicit constructor for particle-in-cell. Should be called
!  after calling GPart_ctor.
!
!  ARGUMENTS:
!    this    : 'this' class instance
!    comm    : MPI communicator
!    ppc     : number of particle per cell
!    inittype: GPINIT-typed quantity to give the type of particle
!              initialization
!    iinterp : GPINTRP-type quantity providing the interpolation
!    intorder: for variable-order (e.g., Lagrange) interpolation,
!              the order (1, 2, 3...). Sets the number of 'ghost' zones
!              of data transferred between MPI tasks.
!              scheme
!    iexchtyp: format for exchanging particle data: GPEXCHTYPE_NN=nearest
!    neighbor,
!              suggesting that data is passed between neighboring
!              tasks only; GPEXCHTYPE_VDB form suggests that all tasks retain a
!              copy
!              of the 'master' particle d.b., so that passing between
!              neighbors is not necessary. But the VDB form does require
!              a global reduction, and may be more expensive.
!    iouttup : output type: 0==binary, 1=ASCII
!    bcoll   : if doing binary I/O, do collective (==1); or not (==0)
!    csize   : cache size param for local transposes
!    nstrip  : 'strip-mining' size for local transposes
!    intacc  : compute acceleration internally to class (==1); or not (==0).
!    Storage
!              allocated only if intacc==1.
!    wrtunit : (optional) write particle positions in box units (==1) (i.e.,
!              x,y,z in [0,2.pi]), or in grid units (==0) (x,y,z in [0,N]).
!-----------------------------------------------------------------
    USE var
    USE grid
    USE boxsize
    USE commtypes

    IMPLICIT NONE
    CLASS(GPIC)      ,INTENT(INOUT)     :: this
    INTEGER          ,INTENT   (IN)     :: bcoll,comm,ppc
    INTEGER          ,INTENT   (IN)     :: csize,nstrip,intacc
    INTEGER, OPTIONAL,INTENT   (IN)     :: wrtunit
    INTEGER                             :: disp(3),lens(3),types(3),szreal
    INTEGER          ,INTENT   (IN)     :: iexchtyp,iinterp,inittype
    INTEGER          ,INTENT   (IN)     :: intorder,iouttyp
    INTEGER                             :: maxparts

!    this%icv_      = real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP)  &
!                     *Dkx*Dky*Dkz/(2*pi)**3
!    this%icv_      = real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP)maxparts
    this%icv_      = 1.0_GP/real(ppc,kind=GP)
    maxparts       = ppc*nx*ny*nz

    CALL this%GPart_ctor(comm,maxparts,inittype,iinterp,3,iexchtyp, &
                        iouttyp,bcoll,csize,nstrip,intacc,wrtunit)
    this%intorder_ = intorder
    CALL this%gfcomm_%GPartComm_ctor(GPCOMM_INTRFC_SF,maxparts,    &
         this%nd_,this%intorder_/2+1,this%comm_,this%htimers_(GPTIME_COMM))
    CALL this%gfcomm_%SetCacheParam(csize,nstrip)
    CALL this%gfcomm_%Init()
    CALL this%gfcomm_%GFieldComm_ctor()

    CALL this%picspl_%GPICSplineInt_ctor(3,this%nd_,this%libnds_,this%lxbnds_, &
         this%tibnds_,this%intorder_,this%intorder_/2+1,this%maxparts_,        &
         this%gfcomm_,this%htimers_(GPTIME_DATAEX),this%htimers_(GPTIME_TRANSP))

    ALLOCATE ( this%prop_(this%maxparts_) )

  END SUBROUTINE GPIC_ctor
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GPIC_Init(this)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : Init
!  DESCRIPTION: Initializes particle locations before integration.
!               Call after construction.
!  ARGUMENTS  :
!    this    : 'this' class instance
!-----------------------------------------------------------------
    IMPLICIT NONE
    CLASS(GPIC)   ,INTENT(INOUT)    :: this
    INTEGER                         :: j

    IF      ( this%inittype_ .EQ. GPINIT_RANDLOC ) THEN
      CALL GPIC_InitRandSeed (this)
    ELSE IF ( this%inittype_ .EQ. GPINIT_USERLOC ) THEN
      CALL GPart_InitUserSeed (this)
    ENDIF

  END SUBROUTINE GPIC_Init
!-----------------------------------------------------------------
!-----------------------------------------------------------------
  SUBROUTINE GPIC_InitRandSeed(this)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : InitRandSeed
!  DESCRIPTION: Initializes particle locations by dividing
!               maxparts evenly among cells, and randomly
!               selecting positions within each cell.
!  ARGUMENTS  :
!    this    : 'this' class instance
!-----------------------------------------------------------------
    USE random
    USE grid
    USE mpivars

    IMPLICIT NONE
    CLASS(GPIC)   ,INTENT(INOUT)      :: this
    REAL(KIND=GP)                     :: r
    INTEGER                           :: ppc,ib,lag,i,j,k,l

    ppc = this%maxparts_/(nx*ny*nz)
    this%nparts_ = nx*ny*(kend - ksta + 1)*ppc
    ib = nx*ny*(ksta-1)*ppc
    lag = 1
    DO i = 1,nx
      DO j = 1,ny
        DO k = ksta,kend
          DO l = 1,ppc
            this%id_(lag) = lag + ib
            CALL prandom_number(r)
            this%px_(lag) = i - 1.0_GP  + r
            CALL prandom_number(r)
            this%py_(lag) = j - 1.0_GP  + r
            CALL prandom_number(r)
            this%pz_(lag) = k - 1.50_GP + r
            lag = lag + 1
          END DO
        END DO
      END DO
    END DO

    PRINT *, myrank, this%nparts_, lag, ib

    CALL this%gpcomm_%VDBSynch(this%vdb_,this%maxparts_,this%id_, &
                          this%px_,this%py_,this%pz_,this%nparts_,this%ptmp1_)
    CALL this%gpcomm_%VDBSynch(this%gptmp0_,this%maxparts_,this%id_, &
                          this%px_,this%py_,this%pz_,this%nparts_,this%ptmp1_)
    CALL GPart_GetLocalWrk(this,this%id_,this%px_,this%py_,this%pz_, &
                           this%nparts_,this%vdb_,this%maxparts_)

    IF ( this%wrtunit_ .EQ. 1 ) THEN ! rescale coordinates to box units
       this%ptmp0_(1,:) = this%vdb_(1,:)*this%delta_(1)
       this%ptmp0_(2,:) = this%vdb_(2,:)*this%delta_(2)
       this%ptmp0_(3,:) = this%vdb_(3,:)*this%delta_(3)
       CALL GPart_ascii_write_lag(this,1,'.','xlgInitRndSeed','000',0.0_GP,&
            this%maxparts_,this%ptmp0_(1,:),this%ptmp0_(2,:),this%ptmp0_(3,:))
    ELSE
       CALL GPart_ascii_write_lag(this,1,'.','xlgInitRndSeed','000',0.0_GP,&
            this%maxparts_,this%vdb_(1,:),this%vdb_(2,:),this%vdb_(3,:))
    ENDIF

    PRINT *, myrank, this%maxparts_, this%nparts_
    IF ( .NOT.GPart_PartNumConsistent(this,this%nparts_) ) THEN
      IF ( this%myrank_.eq.0 ) THEN
        WRITE(*,*) 'GPIC_InitRandSeed: Invalid particle after GetLocalWrk call'
        STOP
      ENDIF
    ENDIF

  END SUBROUTINE GPIC_InitRandSeed
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GPIC_dtor(this)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  Explicit destructor for particle-in-cell. Should be called
!  after calling GPart_dtor.
!
!  ARGUMENTS:
!    this    : 'this' class instance
!-----------------------------------------------------------------
    USE var
    USE grid
    USE boxsize
    USE mpivars
    USE commtypes
    USE random

    IMPLICIT NONE
    TYPE(GPIC)       ,INTENT(INOUT)     :: this

    IF ( ALLOCATED ( this%prop_ ) ) DEALLOCATE ( this%prop_ )

  END SUBROUTINE GPIC_dtor
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GPIC_GetDensity(this,dens)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : GetDensity
!  DESCRIPTION: Deposits particle density into grid.
!
!  ARGUMENTS  :
!    this     : 'this' class instance
!    dens     : Eulerian field containing particle density (OUT)
!-----------------------------------------------------------------
    USE grid
    USE mpivars

    IMPLICIT NONE
    CLASS(GPIC)  ,INTENT(INOUT)                            :: this
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(nx,ny,ksta:kend) :: dens
    INTEGER                                                :: lag

    DO lag=1,this%nparts_
      this%prop_(lag) = this%icv_
    END DO

    CALL GPIC_LagToEuler(this,this%prop_,this%nparts_,dens,.true.)

    RETURN

  END SUBROUTINE GPIC_GetDensity
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GPIC_PerturbPositions(this,d,k,drp)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : PerturbPositions
!  DESCRIPTION: Add small perturbation to particle positions in
!               the x-direction  dx = d*cos(k*x)

!  ARGUMENTS  :
!    this    : 'this' class instance
!    vz,vy,vz: compoments of velocity field, in real space, partially
!              updated, possibly. These will be overwritten!
!    xk      : multiplicative RK time stage factor
!-----------------------------------------------------------------
    USE var
    USE grid
    USE boxsize
    USE mpivars

    IMPLICIT NONE
    CLASS(GPIC)  ,INTENT(INOUT)                            :: this
    REAL(KIND=GP),INTENT(IN)                               :: d
    INTEGER      ,INTENT(IN)                               :: k,drp
    INTEGER                                                :: lag,ng
    REAL(KIND=GP)                                          :: kk,dd
   
    IF (drp.EQ.0) THEN
      dd = REAL(nx,kind=GP)*d/(2*pi*Lx)
      kk = 2*pi*k/REAL(nx,kind=GP)
      DO lag=1,this%nparts_
        PRINT *, this%px_(lag), d*COS(kk*this%px_(lag))
        this%px_(lag) = this%px_(lag) + dd*COS(kk*this%px_(lag))
      END DO
    ELSE IF (drp.EQ.1) THEN
      dd = REAL(ny,kind=GP)*d/(2*pi*Ly)
      kk = 2*pi*k/REAL(ny,kind=GP)
      DO lag=1,this%nparts_
        this%py_(lag) = this%py_(lag) + dd*COS(kk*this%py_(lag))
      END DO
    ELSE IF (drp.EQ.2) THEN
      dd = REAL(nz,kind=GP)*d/(2*pi*Lz)
      kk = 2*pi*k/REAL(nz,kind=GP)
      DO lag=1,this%nparts_
        this%pz_(lag) = this%pz_(lag) + dd*COS(kk*this%pz_(lag))
      END DO
    END IF

!    CALL GPart_MakePeriodicP(this,this%px_,this%py_,this%pz_,this%nparts_,7)

    ! If using nearest-neighbor interface, do particle exchange
    ! between nearest-neighbor tasks BEFORE z-PERIODIZING particle coordinates:
    IF ( this%iexchtype_.EQ.GPEXCHTYPE_NN ) THEN
      CALL GTStart(this%htimers_(GPTIME_COMM))
      CALL this%gpcomm_%PartExchangeV(this%id_,this%px_,this%py_,this%pz_, &
           this%nparts_,this%lxbnds_(3,1),this%lxbnds_(3,2))
      CALL GTAcc(this%htimers_(GPTIME_COMM))
    ENDIF

    ! Enforce periodicity in x, y, & z:
    CALL GPart_MakePeriodicP(this,this%px_,this%py_,this%pz_,this%nparts_,7)

    ! If using VDB interface, do synch-up, and get local work:
    IF ( this%iexchtype_.EQ.GPEXCHTYPE_VDB ) THEN

      IF ( .NOT.GPart_PartNumConsistent(this,this%nparts_) ) THEN
        IF ( this%myrank_.eq.0 ) THEN
          WRITE(*,*) 'GPIC_EndStageRKK: Inconsistent particle count'
        ENDIF
      ENDIF
      ! Synch up VDB, if necessary:
      CALL GTStart(this%htimers_(GPTIME_COMM))
      CALL this%gpcomm_%VDBSynch(this%vdb_,this%maxparts_,this%id_, &
                     this%px_,this%py_,this%pz_,this%nparts_,this%ptmp1_)
      CALL this%gpcomm_%VDBSynch(this%gptmp0_,this%maxparts_,this%id_, &
                     this%ptmp0_(1,:),this%ptmp0_(2,:),this%ptmp0_(3,:),&
                     this%nparts_,this%ptmp1_)
      CALL GTAcc(this%htimers_(GPTIME_COMM))

      ! If using VDB, get local particles to work on:
      ! GPart_GetLocalWrk_aux also synchronizes auxiliary RK arrays,
      ! and is needed if the call is done in the middle of a RK step.
      CALL GPart_GetLocalWrk_aux(this,this%id_,this%px_,this%py_,this%pz_,&
                       this%ptmp0_(1,:),this%ptmp0_(2,:),this%ptmp0_(3,:),&
                       this%nparts_,this%vdb_,this%gptmp0_,this%maxparts_)
      CALL MPI_ALLREDUCE(this%nparts_,ng,1,MPI_INTEGER,   &
                         MPI_SUM,this%comm_,this%ierr_)

      IF ( this%myrank_.EQ.0 .AND. ng.NE.this%maxparts_) THEN
        WRITE(*,*)'GPIC_EndStageRKK: inconsistent d.b.: expected: ', &
                 this%maxparts_, '; found: ',ng
        CALL GPART_ascii_write_lag(this,1,'.','xlgerr','000',0.0_GP, &
                                   this%maxparts_,this%vdb_)
        STOP
      ENDIF

    ENDIF

    IF ( this%intacc_.EQ.0 ) RETURN

    RETURN

  END SUBROUTINE GPIC_PerturbPositions

  SUBROUTINE GPIC_EndStageRKK(this,vx,vy,vz,xk)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : EndStageRKK
!  DESCRIPTION: Called at the end of all RK-like stages to
!               complete Lagrangian particle update.

!  ARGUMENTS  :
!    this    : 'this' class instance
!    vz,vy,vz: compoments of velocity field, in real space, partially
!              updated, possibly. These will be overwritten!
!    xk      : multiplicative RK time stage factor
!-----------------------------------------------------------------
    USE fprecision
    USE commtypes
    USE mpivars
    USE grid

    IMPLICIT NONE
    CLASS(GPIC)  ,INTENT(INOUT)                            :: this
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(nx,ny,ksta:kend) :: vx,vy,vz
    REAL(KIND=GP),INTENT   (IN)                            :: xk
    INTEGER                                                :: j,ng

    ! u(t+dt) = u*: done already

    ! If using nearest-neighbor interface, do particle exchange
    ! between nearest-neighbor tasks BEFORE z-PERIODIZING particle coordinates:
    IF ( this%iexchtype_.EQ.GPEXCHTYPE_NN ) THEN
      CALL GTStart(this%htimers_(GPTIME_COMM))
      CALL this%gpcomm_%PartExchangeV(this%id_,this%px_,this%py_,this%pz_, &
           this%nparts_,this%lxbnds_(3,1),this%lxbnds_(3,2))
      CALL GTAcc(this%htimers_(GPTIME_COMM))
    ENDIF

    ! Enforce periodicity in x, y, & z:
    CALL GPart_MakePeriodicP(this,this%px_,this%py_,this%pz_,this%nparts_,7)

    ! If using VDB interface, do synch-up, and get local work:
    IF ( this%iexchtype_.EQ.GPEXCHTYPE_VDB ) THEN

      IF ( .NOT.GPart_PartNumConsistent(this,this%nparts_) ) THEN
        IF ( this%myrank_.eq.0 ) THEN
          WRITE(*,*) 'GPIC_EndStageRKK: Inconsistent particle count'
        ENDIF
      ENDIF
      ! Synch up VDB, if necessary:
      CALL GTStart(this%htimers_(GPTIME_COMM))
      CALL this%gpcomm_%VDBSynch(this%vdb_,this%maxparts_,this%id_, &
                     this%px_,this%py_,this%pz_,this%nparts_,this%ptmp1_)
      CALL this%gpcomm_%VDBSynch(this%gptmp0_,this%maxparts_,this%id_, &
                     this%ptmp0_(1,:),this%ptmp0_(2,:),this%ptmp0_(3,:),&
                     this%nparts_,this%ptmp1_)
      CALL GTAcc(this%htimers_(GPTIME_COMM))

      ! If using VDB, get local particles to work on:
      ! GPart_GetLocalWrk_aux also synchronizes auxiliary RK arrays,
      ! and is needed if the call is done in the middle of a RK step.
      CALL GPart_GetLocalWrk_aux(this,this%id_,this%px_,this%py_,this%pz_,&
                       this%ptmp0_(1,:),this%ptmp0_(2,:),this%ptmp0_(3,:),&
                       this%nparts_,this%vdb_,this%gptmp0_,this%maxparts_)
      CALL MPI_ALLREDUCE(this%nparts_,ng,1,MPI_INTEGER,   &
                         MPI_SUM,this%comm_,this%ierr_)

      IF ( this%myrank_.EQ.0 .AND. ng.NE.this%maxparts_) THEN
        WRITE(*,*)'GPIC_EndStageRKK: inconsistent d.b.: expected: ', &
                 this%maxparts_, '; found: ',ng
        CALL GPART_ascii_write_lag(this,1,'.','xlgerr','000',0.0_GP, &
                                   this%maxparts_,this%vdb_)
        STOP
      ENDIF

    ENDIF

    IF ( this%intacc_.EQ.0 ) RETURN

    ! If doing internal acceleration, synch up past time levels:
    CALL GPart_synch_acc(this)

    ! Set t^n+1 velocity based on most recent Lag.particle positions:
    ! NOTE: vx, vy, vz are overwirtten on exit:
    CALL GPIC_EulerToLag(this,this%lvx_,this%nparts_,vx,.true. )
    CALL GPIC_EulerToLag(this,this%lvy_,this%nparts_,vy,.false.)
    CALL GPIC_EulerToLag(this,this%lvz_,this%nparts_,vz,.false.)

!$omp parallel do
    DO j = 1, this%nparts_
      this%vk2_(1,j) = this%lvx_(j)
      this%vk2_(2,j) = this%lvy_(j)
      this%vk2_(3,j) = this%lvz_(j)
    ENDDO

    ! set particle ids owned by task at lag times:
    this%idm_ = 0
    this%npartsm_ = this%nparts_
    this%idm_(1:this%nparts_) = this%id_(1:this%nparts_)

    RETURN

  END SUBROUTINE GPIC_EndStageRKK
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GPIC_StepRKK(this, vx, vy, vz, dt, xk)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : StepRKK
!  DESCRIPTION: Carries out one stage of explicit RK-like time
!               inegration step.  Intended for explicit step within
!               an outer stepper method of the form:
!
!               X = X_0 + dt * V(X(t),t) * xk,
!
!               Note that the vx, vy, vz, will be overwritten here.
!               Uses picspline method for interpolation.
!  ARGUMENTS  :
!    this    : 'this' class instance
!    vz,vy,vz: compoments of velocity field, in real space, partially
!              updated, possibly. These will be overwritten!
!    dt      : integration timestep
!    xk      : multiplicative RK time stage factor
!-----------------------------------------------------------------
    USE fprecision
    USE commtypes
    USE mpivars
    USE grid

    IMPLICIT NONE
    CLASS(GPIC)  ,INTENT(INOUT)                            :: this
    INTEGER                                                :: i,j
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(nx,ny,ksta:kend) :: vx,vy,vz
    REAL(KIND=GP),INTENT   (IN)                            :: dt,xk
    REAL(KIND=GP)                                          :: dtfact
    REAL(KIND=GP),ALLOCATABLE  ,DIMENSION              (:) :: lid,gid

    CALL GTStart(this%htimers_(GPTIME_STEP))

    ! Find F(u*):
    ! ... x:
    CALL GPIC_EulerToLag(this,this%lvx_,this%nparts_,vx,.true.)
    ! ux* <-- ux + dt * F(U*)*xk:
    dtfact = dt*xk*this%invdel_(1)
  !$omp parallel do
    DO j = 1, this%nparts_
      this%px_(j) = this%ptmp0_(1,j) + dtfact*this%lvx_(j)
    ENDDO

    !
    ! ... y:
    ! Exchange bdy data for velocities, so that we
    ! can perform local interpolations:
    CALL GPIC_EulerToLag(this,this%lvy_,this%nparts_,vy,.false.)
    ! uy* <-- uy + dt * F(U*)*xk:
    dtfact = dt*xk*this%invdel_(2)
  !$omp parallel do
    DO j = 1, this%nparts_
      this%py_(j) = this%ptmp0_(2,j) + dtfact*this%lvy_(j)
    ENDDO

    ! ... z:
    ! Exchange bdy data for velocities, so that we
    ! can perform local interpolations:
    CALL GPIC_EulerToLag(this,this%lvz_,this%nparts_,vz,.false.)
    ! uz* <-- uz + dt * F(U*)*xk:
    dtfact = dt*xk*this%invdel_(3)
  !$omp parallel do
    DO j = 1, this%nparts_
      this%pz_(j) = this%ptmp0_(3,j) + dtfact*this%lvz_(j)
    ENDDO

    ! Enforce periodicity in x-y only:
    CALL GPart_MakePeriodicP(this,this%px_,this%py_,this%pz_,this%nparts_,3)

    CALL GTAcc(this%htimers_(GPTIME_STEP))

  !   ALLOCATE  (lid(this%maxparts_))
  !   ALLOCATE  (gid(this%maxparts_))
  !   lid = 0
  !   gid = 0
  !   do j=1,this%maxparts_
  !     lid(this%id_(j)+1) = 1
  !   enddo
  !   call mpi_allreduce(lid,gid,this%maxparts_,MPI_INTEGER,MPI_SUM,this%comm_,this%ierr_)
  !   if ( this%myrank_.eq. 0 ) then
  !     do j=1,this%maxparts_
  !       if ( gid(j) .gt.1 ) then
  !         write(*,*)this%myrank_,': StepRKK: multiplicity > 1: id=',j-1,' p=', &
  !                   this%vdb_(1,j),this%vdb_(2,j),this%vdb_(3,j)
  !       else if ( gid(j).eq.0 ) then
  !         write(*,*)this%myrank_,': StepRKK: particle missing: id=',j-1
  !       endif
  !     enddo
  !   endif
  !
  !   DEALLOCATE(lid,gid)

  ! At this point, the vx, vy, vz should be intact:
    CALL GPIC_EndStageRKK(this,vx,vy,vz,xk)

  END SUBROUTINE GPIC_StepRKK
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GPIC_EulerToLag(this,lag,nl,evar,doupdate)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : EulerToLag
!  DESCRIPTION: Computes the Eulerian to Lagrangian
!               transformation by interpolating Eulerian field
!               evar to current position of Lagrangian paricles in
!               d.b. Array lag must be large enough to accommodate
!               max. no. particles; no checking is done. Note
!               that 'evar' array must have local dimensions
!               for a real array in GHOST (nx X ny X (kend-ksta+1)).
!               Global computation of spline or other interpolation
!               operator will be done here.
!
!  ARGUMENTS  :
!    this    : 'this' class instance (IN)
!    lag     : real array with the containing Lagrangian value
!              of evar field on output (OUT)
!    nl      : no. Lag. points in lag
!    evar    : Eulerian variable (IN)
!    doupdate: if true, do interp point update in interpolator; else don't
!
!-----------------------------------------------------------------
    USE grid
    USE fprecision
    USE commtypes
    USE mpivars

    IMPLICIT NONE
    CLASS(GPIC)  ,INTENT(INOUT)                            :: this
    INTEGER      ,INTENT   (IN)                            :: nl
    LOGICAL      ,INTENT   (IN)                            :: doupdate
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(nx,ny,ksta:kend) :: evar
    REAL(KIND=GP),INTENT(INOUT),DIMENSION             (nl) :: lag
    INTEGER                                                :: j

    IF ( doupdate ) THEN
      CALL GTStart(this%htimers_(GPTIME_PUPDATE))
      CALL this%picspl_%PartUpdate3D(this%px_,this%py_,this%pz_,this%nparts_)
      CALL GTAcc(this%htimers_(GPTIME_PUPDATE))
    ENDIF
    CALL GTStart(this%htimers_(GPTIME_SPLINE))
    CALL this%picspl_%CompSpline3D(evar)
    CALL GTAcc(this%htimers_(GPTIME_SPLINE))

    CALL GTStart(this%htimers_(GPTIME_INTERP))
    CALL this%picspl_%DoInterp3D(lag,nl)
    CALL GTAcc(this%htimers_(GPTIME_INTERP))

  END SUBROUTINE GPIC_EulerToLag
!-----------------------------------------------------------------
!-----------------------------------------------------------------

SUBROUTINE GPIC_LagToEuler(this,lag,nl,evar,doupdate)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : EulerToLag
!  DESCRIPTION: Computes the Lagrangian to Eulerian
!               transformation by depositing Lagrangian particles
!               property in their current position (in d.b.) to
!               Eulerian field evar interpolating Eulerian field
!               Array lag must be large enough to accommodate
!               max. no. particles; no checking is done. Note
!               that 'evar' array must have local dimensions
!               for a real array in GHOST (nx X ny X (kend-ksta+1)).
!
!  ARGUMENTS  :
!    this    : 'this' class instance (IN)
!    lag     : real array with the containing Lagrangian property
!              to be deposited to the grid (IN)
!    nl      : no. Lag. points in lag
!    evar    : Eulerian variable (OUT)
!    doupdate: if true, do interp point update in interpolator; else don't
!
!-----------------------------------------------------------------
    USE grid
    USE fprecision
    USE commtypes
    USE mpivars

    IMPLICIT NONE
    CLASS(GPIC)  ,INTENT(INOUT)                            :: this
    INTEGER      ,INTENT   (IN)                            :: nl
    LOGICAL      ,INTENT   (IN)                            :: doupdate
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(nx,ny,ksta:kend) :: evar
    REAL(KIND=GP),INTENT(INOUT),DIMENSION             (nl) :: lag
    INTEGER                                                :: j

    IF ( doupdate ) THEN
      CALL GTStart(this%htimers_(GPTIME_PUPDATE))
      CALL this%picspl_%PartUpdate3D(this%px_,this%py_,this%pz_,this%nparts_)
      CALL GTAcc(this%htimers_(GPTIME_PUPDATE))
    ENDIF

    CALL GTStart(this%htimers_(GPTIME_INTERP))
    CALL this%picspl_%DoDeposit3D(lag,nl,evar)
    CALL GTAcc(this%htimers_(GPTIME_INTERP))

  END SUBROUTINE GPIC_LagToEuler
!-----------------------------------------------------------------
!-----------------------------------------------------------------


!=================================================================
! VGPIC SUBROUTINES
!=================================================================

  SUBROUTINE VGPIC_io_write_pdbv(this, iunit, dir, spref, nmb, time)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : io_write_pdbv
!  DESCRIPTION: Does write of particle velocity d.b. to file.
!               Position of the particle structure in file is the
!               particle's id. Main entry point for both ASCII and
!               binary writes.
!  ARGUMENTS  :
!    this    : 'this' class instance
!    iunit   : unit number
!    dir     : output directory
!    spref   : filename prefix
!    nmd     : time index
!    time    : real time
!
!-----------------------------------------------------------------
    USE fprecision
    USE commtypes
    USE mpivars

    IMPLICIT NONE
    CLASS(VGPIC),INTENT(INOUT)        :: this
    REAL(KIND=GP),INTENT   (IN)       :: time
    REAL(KIND=GP)                     :: prec(3)
    INTEGER,INTENT(IN)                :: iunit
    INTEGER                           :: fh,j,nt
    INTEGER(kind=MPI_OFFSET_KIND)     :: offset
    CHARACTER(len=*),INTENT(IN)       :: dir
    CHARACTER(len=*),INTENT(IN)       :: nmb
    CHARACTER(len=*),INTENT(IN)       :: spref

    CALL this%gpcomm_%VDBSynch(this%ptmp0_,this%maxparts_,this%id_, &
         this%pvx_,this%pvy_,this%pvz_,this%nparts_,this%ptmp1_)

    IF ( this%iouttype_ .EQ. 0 ) THEN
      IF ( this%bcollective_.EQ. 1 ) THEN
        ! pass in the current linear _local_ particle velocity arrays
        CALL GPart_binary_write_lag_co(this,iunit,dir,spref,nmb,time,this%nparts_, &
             this%pvx_,this%pvy_,this%pvz_)
      ELSE
        ! pass in the synched-up VDB (copied to ptmp0_):
        CALL GPart_binary_write_lag_t0(this,iunit,dir,spref,nmb,time,this%maxparts_, &
             this%ptmp0_(1,:),this%ptmp0_(2,:),this%ptmp0_(3,:))
      ENDIF
    ELSE
      ! pass in the synched-up VDB (copied to ptmp0_):
      CALL GPart_ascii_write_lag(this,iunit,dir,spref,nmb,time,this%maxparts_, &
           this%ptmp0_(1,:),this%ptmp0_(2,:),this%ptmp0_(3,:))
    ENDIF

  END SUBROUTINE VGPIC_io_write_pdbv
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE VGPIC_io_read_pdbv(this, iunit, dir, spref, nmb)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : io_read_pdbv
!  DESCRIPTION: Does read of test particle velocity from file.
!               This is the main entry point for both binary and
!               ASCII reads.
!  ARGUMENTS  :
!    this    : 'this' class instance
!    iunit   : unit number
!    dir     : input directory
!    spref   : filename prefix
!    nmb     : time index
!-----------------------------------------------------------------
    USE fprecision
    USE commtypes
    USE mpivars

    IMPLICIT NONE
    CLASS(VGPIC)    ,INTENT(INOUT)            :: this
    REAL(KIND=GP)                             :: time
    INTEGER,INTENT(IN)                        :: iunit
    INTEGER                                   :: ng
    CHARACTER(len=*),INTENT   (IN)            :: dir
    CHARACTER(len=*),INTENT   (IN)            :: nmb
    CHARACTER(len=*),INTENT   (IN)            :: spref

    CALL GTStart(this%htimers_(GPTIME_GPREAD))
    IF ( this%iouttype_ .EQ. 0 ) THEN        ! Binary files
      IF ( this%bcollective_ .EQ. 1 ) THEN   ! collective binary
        IF (len_trim(nmb).gt.0 ) THEN
        CALL GPart_binary_read_pdb_co(this,iunit, &
        trim(dir) // '/' // trim(spref) // '.' // nmb // '.lag',time,this%ptmp0_)
        ELSE
        CALL GPart_binary_read_pdb_co(this,iunit, trim(spref),time,this%ptmp0_)
        ENDIF
      ELSE                      ! master thread binary
        IF (len_trim(nmb).gt.0 ) THEN
        CALL GPart_binary_read_pdb_t0(this,iunit,&
         trim(dir) // '/' // trim(spref) // '.' // nmb // '.lag',time,this%ptmp0_)
        ELSE
        CALL GPart_binary_read_pdb_t0(this,iunit, trim(spref),time,this%ptmp0_)
        ENDIF
      ENDIF
    ELSE                         ! ASCII files
      IF (len_trim(nmb).gt.0 ) THEN
      CALL GPart_ascii_read_pdb (this,iunit,&
            trim(dir) // '/' // trim(spref) // '.' // nmb // '.txt',time,this%ptmp0_)
      ELSE
      CALL GPart_ascii_read_pdb (this,iunit,trim(spref),time,this%ptmp0_)
      ENDIF
    ENDIF
    CALL GTAcc(this%htimers_(GPTIME_GPREAD))

    ! Store in particle velocity arrays
    CALL GPart_CopyLocalWrk(this,this%pvx_,this%pvy_,this%pvz_, &
                            this%vdb_,this%ptmp0_,this%maxparts_)

    CALL MPI_ALLREDUCE(this%nparts_,ng,1,MPI_INTEGER,   &
                       MPI_SUM,this%comm_,this%ierr_)
    IF ( this%myrank_.EQ.0 .AND. ng.NE.this%maxparts_ ) THEN
      WRITE(*,*)'InerGPart_io_read_pdbv: inconsistent d.b.: expected: ', &
                 this%maxparts_, '; found: ',ng
      STOP
    ENDIF

  END SUBROUTINE VGPIC_io_read_pdbv
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE VGPIC_GetFlux(this,jx,jy,jz)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : GetFlux
!  DESCRIPTION: Deposits particle current density into grid.
!
!  ARGUMENTS  :
!    this     : 'this' class instance
!    jx,jy,jz : Eulerian field containing particle current
!               density (OUT)
!-----------------------------------------------------------------
    USE grid
    USE mpivars

    IMPLICIT NONE
    CLASS(VGPIC)  ,INTENT(INOUT)                           :: this
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(nx,ny,ksta:kend) :: jx,jy,jz
    INTEGER                                                :: lag

! x-coord
    DO lag=1,this%nparts_
      this%prop_(lag) = this%icv_*this%pvx_(lag)
    END DO

    CALL GPIC_LagToEuler(this,this%prop_,this%nparts_,jx,.true.)

! y-coord
    DO lag=1,this%nparts_
      this%prop_(lag) = this%icv_*this%pvy_(lag)
    END DO
    CALL GPIC_LagToEuler(this,this%prop_,this%nparts_,jy,.false.)

! z-coord
    DO lag=1,this%nparts_
      this%prop_(lag) = this%icv_*this%pvz_(lag)
    END DO
    CALL GPIC_LagToEuler(this,this%prop_,this%nparts_,jz,.false.)

    RETURN

  END SUBROUTINE VGPIC_GetFlux
!-----------------------------------------------------------------
!-----------------------------------------------------------------


!=================================================================
! ChargPIC SUBROUTINES
!=================================================================

  SUBROUTINE ChargPIC_ctor(this)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  Explicit constructor for test particles. Should be called
!  after calling GPart_ctor. Calling GPIC_ctor is not required.
!
!  ARGUMENTS:
!    this    : 'this' class instance
!-----------------------------------------------------------------
    USE var
    USE grid
    USE boxsize
    USE commtypes

    IMPLICIT NONE
    CLASS(ChargPIC)   , INTENT(INOUT)   :: this

    ALLOCATE(this%pvx_     (this%maxparts_))
    ALLOCATE(this%pvy_     (this%maxparts_))
    ALLOCATE(this%pvz_     (this%maxparts_))
    ALLOCATE(this%lbx_     (this%maxparts_))
    ALLOCATE(this%lby_     (this%maxparts_))
    ALLOCATE(this%lbz_     (this%maxparts_))
    ALLOCATE(this%lfx_     (this%maxparts_))
    ALLOCATE(this%lfy_     (this%maxparts_))
    ALLOCATE(this%lfz_     (this%maxparts_))
    ALLOCATE(this%ttmp0_ (3,this%maxparts_))

  END SUBROUTINE ChargPIC_ctor
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE ChargPIC_dtor(this)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  Explicit destructor for test particles. Should be called
!  after GPart_dtor.
!
!  ARGUMENTS:
!    this    : 'this' class instance
!-----------------------------------------------------------------

    IMPLICIT NONE
    TYPE(ChargPIC)    ,INTENT(INOUT)    :: this

    IF ( ALLOCATED  (this%prop_) ) DEALLOCATE  (this%prop_)
    IF ( ALLOCATED   (this%pvx_) ) DEALLOCATE   (this%pvx_)
    IF ( ALLOCATED   (this%pvy_) ) DEALLOCATE   (this%pvy_)
    IF ( ALLOCATED   (this%pvz_) ) DEALLOCATE   (this%pvz_)
    IF ( ALLOCATED   (this%lbx_) ) DEALLOCATE   (this%lbx_)
    IF ( ALLOCATED   (this%lby_) ) DEALLOCATE   (this%lby_)
    IF ( ALLOCATED   (this%lbz_) ) DEALLOCATE   (this%lbz_)
    IF ( ALLOCATED   (this%lfx_) ) DEALLOCATE   (this%lfx_)
    IF ( ALLOCATED   (this%lfy_) ) DEALLOCATE   (this%lfy_)
    IF ( ALLOCATED   (this%lfz_) ) DEALLOCATE   (this%lfz_)
    IF ( ALLOCATED (this%ttmp0_) ) DEALLOCATE (this%ttmp0_)

  END SUBROUTINE ChargPIC_dtor
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE ChargPIC_SetStepRKK(this)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : SetStepRKK
!  DESCRIPTION: Initializes an explicit integration timestep for
!               the particle velocity. Must be called at the start
!               of the RK stage execution, together with
!               GPart_SetStepRKK.
!
!  ARGUMENTS  :
!    this     : 'this' class instance
!-----------------------------------------------------------------
    USE fprecision
    USE commtypes
    USE mpivars

    IMPLICIT NONE
    CLASS(ChargPIC)   , INTENT(INOUT) :: this

    INTEGER                           :: j

    ! Initialize solution, v (particle velocity):
    ! v* <-- v:

    ! Cycle over JST loop to update state:
!$omp parallel do
    DO j = 1, this%nparts_
       this%ttmp0_(1,j) = this%pvx_(j)  ! ux_0
       this%ttmp0_(2,j) = this%pvy_(j)  ! uy_0
       this%ttmp0_(3,j) = this%pvz_(j)  ! uz_0
    ENDDO

  END SUBROUTINE ChargPIC_SetStepRKK
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE ChargPIC_StepRKK(this, Ex, Ey, Ez, Bx, By, Bz, dt, xk)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : Step_chargedpic
!  DESCRIPTION: Carries out one stage of explicit RK-like time
!               integration step.  Intended for explicit step within
!               an outer stepper method of the form:
!
!               X = X_0 + dt * V[X(t),t] * xk,
!               V = V_0 + dt * F[V(X(t)),E(X(t)),B(X(t))] * xk,
!
!               where F is the electromagnetic force on the particle.
!               Note that the vx, vy, vz, will be overwritten here.
!  ARGUMENTS  :
!    this    : 'this' class instance
!    Ez,Ey,Ez: compoments of electric field, in real space, partially
!              updated, possibly. These will be overwritten!
!    Bz,By,Bz: compoments of magnetic field in real space
!    dt      : integration timestep
!    xk      : multiplicative RK time stage factor
!-----------------------------------------------------------------
    USE grid
    USE fprecision
    USE commtypes
    USE mpivars

    IMPLICIT NONE
    CLASS(ChargPIC)    ,INTENT(INOUT)                      :: this
    INTEGER                                                :: i,j
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(nx,ny,ksta:kend) :: Bx,By,Bz
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(nx,ny,ksta:kend) :: Ex,Ey,Ez
    REAL(KIND=GP),INTENT   (IN)                            :: dt,xk
    REAL(KIND=GP)                                          :: dtfact
    REAL(KIND=GP)                                          :: dtv
    REAL(KIND=GP), ALLOCATABLE, DIMENSION              (:) :: lid,gid

    dtv    = dt*xk
    CALL GTStart(this%htimers_(GPTIME_STEP))

    ! Find F(u*):
    CALL GPIC_EulerToLag(this,this%lbx_,this%nparts_,Bx,.true. )
    CALL GPIC_EulerToLag(this,this%lby_,this%nparts_,By,.false.)
    CALL GPIC_EulerToLag(this,this%lbz_,this%nparts_,Bz,.false.)
    CALL GPIC_EulerToLag(this,this%lfx_,this%nparts_,Ex,.false.)
    CALL GPIC_EulerToLag(this,this%lfy_,this%nparts_,Ey,.false.)
    CALL GPIC_EulerToLag(this,this%lfz_,this%nparts_,Ez,.false.)
!   Lorentz force
!$omp parallel do
    DO j = 1, this%nparts_
       this%lfx_(j) = this%lfx_(j)+this%pvy_(j)*this%lbz_(j)-this%pvz_(j)*this%lby_(j)
       this%lfy_(j) = this%lfy_(j)+this%pvz_(j)*this%lbx_(j)-this%pvx_(j)*this%lbz_(j)
       this%lfz_(j) = this%lfz_(j)+this%pvx_(j)*this%lby_(j)-this%pvy_(j)*this%lbx_(j)
    ENDDO

    ! ... x:
    dtfact = dt*xk*this%invdel_(1)
!$omp parallel do
    DO j = 1, this%nparts_
      this%px_ (j) = this%ptmp0_(1,j) + dtfact*this%pvx_(j)
      this%pvx_(j) = this%ttmp0_(1,j) + dtv*this%lfx_(j)
    ENDDO

    ! ... y:
    dtfact = dt*xk*this%invdel_(2)
!$omp parallel do
    DO j = 1, this%nparts_
      this%py_ (j) = this%ptmp0_(2,j) + dtfact*this%pvy_(j)
      this%pvy_(j) = this%ttmp0_(2,j) + dtv*this%lfy_(j)
    ENDDO

    ! ... z:
    dtfact = dt*xk*this%invdel_(3)
!$omp parallel do
    DO j = 1, this%nparts_
      this%pz_ (j) = this%ptmp0_(3,j) + dtfact*this%pvz_(j)
      this%pvz_(j) = this%ttmp0_(3,j) + dtv*this%lfz_(j)
    ENDDO

    ! Enforce periodicity in x-y only:
    CALL GPart_MakePeriodicP(this,this%px_,this%py_,this%pz_,this%nparts_,3)

    CALL GTAcc(this%htimers_(GPTIME_STEP))

!   ALLOCATE  (lid(this%maxparts_))
!   ALLOCATE  (gid(this%maxparts_))
!   lid = 0
!   gid = 0
!   do j=1,this%maxparts_
!     lid(this%id_(j)+1) = 1
!   enddo
!   call mpi_allreduce(lid,gid,this%maxparts_,MPI_INTEGER,MPI_SUM,this%comm_,this%ierr_)
!   if ( this%myrank_.eq. 0 ) then
!     do j=1,this%maxparts_
!       if ( gid(j) .gt.1 ) then
!         write(*,*)this%myrank_,': StepRKK: multiplicity > 1: id=',j-1,' p=', &
!                   this%vdb_(1,j),this%vdb_(2,j),this%vdb_(3,j)
!       else if ( gid(j).eq.0 ) then
!         write(*,*)this%myrank_,': StepRKK: particle missing: id=',j-1
!       endif
!     enddo
!   endif
!   DEALLOCATE(lid,gid)

    CALL ChargPIC_EndStageRKK(this,Ex,Ey,Ez,xk)

  END SUBROUTINE ChargPIC_StepRKK
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE ChargPIC_EndStageRKK(this,vx,vy,vz,xk)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : EndStageRKK
!  DESCRIPTION: Called at the end of all RK-like stages to
!               complete inertial particle update.

!  ARGUMENTS  :
!    this    : 'this' class instance
!    Ez,Ey,Ez: compoments of velocity field, in real space, partially
!              updated, possibly. These will be overwritten!
!    xk      : multiplicative RK time stage factor
!-----------------------------------------------------------------
    USE fprecision
    USE commtypes
    USE mpivars
    USE grid

    IMPLICIT NONE
    CLASS(ChargPIC)    ,INTENT(INOUT)                      :: this
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(nx,ny,ksta:kend) :: vx,vy,vz
    REAL(KIND=GP),INTENT   (IN)                            :: xk
    INTEGER                                                :: j,ng

    ! u(t+dt) = u*: done already

    ! IF using nearest-neighbor interface, do particle exchange
    ! between nearest-neighbor tasks BEFORE z-PERIODIZING particle coordinates.
    ! Note this interface has not been tested yet for test particles.
    IF ( this%iexchtype_.EQ.GPEXCHTYPE_NN ) THEN
      CALL GTStart(this%htimers_(GPTIME_COMM))
      CALL this%gpcomm_%PartExchangeV(this%id_,this%px_,this%py_,this%pz_, &
           this%nparts_,this%lxbnds_(3,1),this%lxbnds_(3,2))
      CALL GTAcc(this%htimers_(GPTIME_COMM))
    ENDIF

    ! Enforce periodicity in x, y, & z:
    CALL GPart_MakePeriodicP(this,this%px_,this%py_,this%pz_,this%nparts_,7)

    ! If using VDB interface, do synch-up, and get local work:
    IF ( this%iexchtype_.EQ.GPEXCHTYPE_VDB ) THEN

      IF ( .NOT.GPart_PartNumConsistent(this,this%nparts_) ) THEN
        IF ( this%myrank_.eq.0 ) THEN
          WRITE(*,*) 'ChargPIC_EndStepRKK: Inconsistent particle count'
        ENDIF
      ENDIF
      ! Synch up VDB, if necessary:
      CALL GTStart(this%htimers_(GPTIME_COMM))
      CALL this%gpcomm_%VDBSynch(this%vdb_,this%maxparts_,this%id_, &
                     this%px_,this%py_,this%pz_,this%nparts_,this%ptmp1_)
      CALL this%gpcomm_%VDBSynch(this%gptmp0_,this%maxparts_,this%id_, &
                     this%pvx_,this%pvy_,this%pvz_,this%nparts_,this%ptmp1_)
      CALL GPart_CopyLocalWrk(this,this%pvx_,this%pvy_,this%pvz_, &
                     this%vdb_,this%gptmp0_,this%maxparts_)
      CALL this%gpcomm_%VDBSynch(this%gptmp0_,this%maxparts_,this%id_, &
                     this%ttmp0_(1,:),this%ttmp0_(2,:),this%ttmp0_(3,:),&
                     this%nparts_,this%ptmp1_)
      CALL GPart_CopyLocalWrk(this,this%ttmp0_(1,:),this%ttmp0_(2,:), &
                     this%ttmp0_(3,:),this%vdb_,this%gptmp0_,this%maxparts_)
      CALL this%gpcomm_%VDBSynch(this%gptmp0_,this%maxparts_,this%id_, &
                     this%ptmp0_(1,:),this%ptmp0_(2,:),this%ptmp0_(3,:),&
                     this%nparts_,this%ptmp1_)
      CALL GPart_GetLocalWrk_aux(this,this%id_,this%px_,this%py_,this%pz_,&
                       this%ptmp0_(1,:),this%ptmp0_(2,:),this%ptmp0_(3,:),&
                       this%nparts_,this%vdb_,this%gptmp0_,this%maxparts_)
      CALL GTAcc(this%htimers_(GPTIME_COMM))

      CALL MPI_ALLREDUCE(this%nparts_,ng,1,MPI_INTEGER,   &
                         MPI_SUM,this%comm_,this%ierr_)

    ! vx, vy, vz were not used so far. Can be used in the future
    ! to compute acceleration, or as temporary arrays.

      IF ( this%myrank_.EQ.0 .AND. ng.NE.this%maxparts_) THEN
        WRITE(*,*)'ChargPIC_EndStepRKK: inconsistent d.b.: expected: ', &
                 this%maxparts_, '; found: ',ng
        STOP
      ENDIF

    ENDIF

    RETURN

  END SUBROUTINE ChargPIC_EndStageRKK
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE ChargPIC_InitVel(this,vtherm)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : InitVel
!  DESCRIPTION: Initializes particle velocities with a Gaussian
!               distribution function, with thermal speed vtherm.
!               Other parameters are initialized with GPart_Init.
!  ARGUMENTS:
!    this   : 'this' class instance
!    vtherm : thermal speed
!-----------------------------------------------------------------
    USE random

    IMPLICIT NONE
    CLASS(ChargPIC)   , INTENT(INOUT) :: this
    REAL(KIND=GP),INTENT(IN)        :: vtherm
    REAL(KIND=GP)                   :: low,twopi,u1,u2,gauss
    INTEGER                         :: j

    twopi = 8.0_GP*atan(1.0_GP)
    low   = 1.0e-20_GP

!$omp parallel do
    DO j = 1, this%nparts_

     CALL prandom_number(u1)
     CALL prandom_number(u2)
     u1 =  MAX( u1, low)
     u2 =  twopi * u2
     gauss =  sqrt( -2.0_GP * log( u1)) * CMPLX( cos(u2), sin(u2))

       this%pvx_(j) = vtherm*gauss

     CALL prandom_number(u1)
     CALL prandom_number(u2)
     u1 =  MAX( u1, low)
     u2 =  twopi * u2
     gauss =  sqrt( -2.0_GP * log( u1)) * CMPLX( cos(u2), sin(u2))

       this%pvy_(j) = vtherm*gauss

     CALL prandom_number(u1)
     CALL prandom_number(u2)
     u1 =  MAX( u1, low)
     u2 =  twopi * u2
     gauss =  sqrt( -2.0_GP * log( u1)) * CMPLX( cos(u2), sin(u2))

       this%pvz_(j) = vtherm*gauss

    END DO

  END SUBROUTINE ChargPIC_InitVel
!-----------------------------------------------------------------
!-----------------------------------------------------------------


  SUBROUTINE ChargPIC_GetTemperature(this,T)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : GetTemperature
!  DESCRIPTION: Deposits particle temperature into grid.
!
!  ARGUMENTS  :
!    this     : 'this' class instance
!    T        : Eulerian field containing particle temperature (OUT)
!-----------------------------------------------------------------
    USE grid
    USE mpivars

    IMPLICIT NONE
    CLASS(ChargPIC),INTENT(INOUT)                          :: this
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(nx,ny,ksta:kend) :: T
    INTEGER                                                :: lag

    DO lag=1,this%nparts_
      this%prop_(lag) = this%icv_*(this%pvx_(lag)**2+this%pvy_(lag)**2 &
                                  +this%pvz_(lag)**2)
    END DO

    CALL GPIC_LagToEuler(this,this%prop_,this%nparts_,T,.true.)

    RETURN

  END SUBROUTINE ChargPIC_GetTemperature
!-----------------------------------------------------------------
!-----------------------------------------------------------------
