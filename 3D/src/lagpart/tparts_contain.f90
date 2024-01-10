!=================================================================
! TestGPart SUBROUTINES
!=================================================================

  SUBROUTINE TestGPart_ctor(this)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  Explicit constructor for test particles. Should be called
!  after calling GPart_ctor.
!
!  ARGUMENTS:
!    this    : 'this' class instance
!-----------------------------------------------------------------

    IMPLICIT NONE
    CLASS(TestGPart), INTENT(INOUT)     :: this

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

  END SUBROUTINE TestGPart_ctor
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE TestGPart_dtor(this)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  Explicit destructor for test particles. Should be called
!  after GPart_dtor.
!
!  ARGUMENTS:
!    this    : 'this' class instance
!-----------------------------------------------------------------

    IMPLICIT NONE
    TYPE(TestGPart)   ,INTENT(INOUT)    :: this

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

  END SUBROUTINE TestGPart_dtor
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE TestGPart_InitVel(this,vtherm)
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
    CLASS(TestGPart), INTENT(INOUT) :: this
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

  END SUBROUTINE TestGPart_InitVel
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE TestGPart_SetStepRKK(this)
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
    CLASS(TestGPart), INTENT(INOUT) :: this

    INTEGER                         :: j

    ! Initialize solution, v (particle velocity): 
    ! v* <-- v: 
 
    ! Cycle over JST loop to update state:
!$omp parallel do
    DO j = 1, this%nparts_
       this%ttmp0_(1,j) = this%pvx_(j)  ! ux_0
       this%ttmp0_(2,j) = this%pvy_(j)  ! uy_0
       this%ttmp0_(3,j) = this%pvz_(j)  ! uz_0
    ENDDO
    
  END SUBROUTINE TestGPart_SetStepRKK
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE TestGPart_StepRKK(this, vx, vy, vz, bx, by, bz, jx, &
                jy, jz, dt, xk, tmp1, tmp2, tmp3)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : Step_testp
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
!    vz,vy,vz: compoments of velocity field, in real space, partially
!              updated, possibly. These will be overwritten!
!    bz,by,bz: compoments of magnetic field in real space
!    jz,jy,jz: compoments of current density in real space
!    dt      : integration timestep
!    xk      : multiplicative RK time stage factor
!    tmpX    : temp arrays the same size as vx, vy, vz
!-----------------------------------------------------------------
    USE grid
    USE fprecision
    USE commtypes
    USE mpivars

    IMPLICIT NONE
    CLASS(TestGPart) ,INTENT(INOUT)                        :: this
    INTEGER                                                :: i,j
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(nx,ny,ksta:kend) :: vx,vy,vz
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(nx,ny,ksta:kend) :: tmp1,tmp2,tmp3
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(nx,ny,ksta:kend) :: bx,by,bz
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(nx,ny,ksta:kend) :: jx,jy,jz
    REAL(KIND=GP),INTENT   (IN)                            :: dt,xk
    REAL(KIND=GP)                                          :: dtfact
    REAL(KIND=GP)                                          :: dtv
    REAL(KIND=GP), ALLOCATABLE, DIMENSION              (:) :: lid,gid

    dtv    = dt*xk
    CALL GTStart(this%htimers_(GPTIME_STEP))

    ! Find F(u*):
    CALL GPart_EulerToLag(this,this%lvx_,this%nparts_,vx,.true.,tmp1,tmp2)
    CALL GPart_EulerToLag(this,this%lvy_,this%nparts_,vy,.false.,tmp1,tmp2)
    CALL GPart_EulerToLag(this,this%lvz_,this%nparts_,vz,.false.,tmp1,tmp2)
    CALL GPart_EulerToLag(this,this%lbx_,this%nparts_,bx,.false.,tmp1,tmp2)
    CALL GPart_EulerToLag(this,this%lby_,this%nparts_,by,.false.,tmp1,tmp2)
    CALL GPart_EulerToLag(this,this%lbz_,this%nparts_,bz,.false.,tmp1,tmp2)
    CALL GPart_EulerToLag(this,this%lfx_,this%nparts_,jx,.false.,tmp1,tmp2)
    CALL GPart_EulerToLag(this,this%lfy_,this%nparts_,jy,.false.,tmp1,tmp2)
    CALL GPart_EulerToLag(this,this%lfz_,this%nparts_,jz,.false.,tmp1,tmp2)
!   Lorentz force
!$omp parallel do
    DO j = 1, this%nparts_
       this%lfx_(j) = this%lfx_(j) + (this%pvy_(j)-this%lvy_(j))*this%lbz_(j)- &
                      (this%pvz_(j)-this%lvz_(j))*this%lby_(j)
       this%lfy_(j) = this%lfy_(j) + (this%pvz_(j)-this%lvz_(j))*this%lbx_(j)- &
                      (this%pvx_(j)-this%lvx_(j))*this%lbz_(j)
       this%lfz_(j) = this%lfz_(j) + (this%pvx_(j)-this%lvx_(j))*this%lby_(j)- &
                      (this%pvy_(j)-this%lvy_(j))*this%lbx_(j)
! Lorentz force only with magnetic field (conservative)
!      this%lfx_(j) = this%pvy_(j)*this%lbz_(j)-this%pvz_(j)*this%lby_(j)
!      this%lfy_(j) = this%pvz_(j)*this%lbx_(j)-this%pvx_(j)*this%lbz_(j)
!      this%lfz_(j) = this%pvx_(j)*this%lby_(j)-this%pvy_(j)*this%lbx_(j)
    ENDDO

    ! ... x:
    dtfact = dt*xk*this%invdel_(1)
!$omp parallel do
    DO j = 1, this%nparts_
      this%px_(j) = this%ptmp0_(1,j) + dtfact*this%pvx_(j)
      this%pvx_(j) = this%ttmp0_(1,j) + dtv*this%lfx_(j)
    ENDDO

    ! ... y:
    dtfact = dt*xk*this%invdel_(2)
!$omp parallel do
    DO j = 1, this%nparts_
      this%py_(j) = this%ptmp0_(2,j) + dtfact*this%pvy_(j)
      this%pvy_(j) = this%ttmp0_(2,j) + dtv*this%lfy_(j)
    ENDDO

    ! ... z:
    dtfact = dt*xk*this%invdel_(3)
!$omp parallel do
    DO j = 1, this%nparts_
      this%pz_(j) = this%ptmp0_(3,j) + dtfact*this%pvz_(j)
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

    CALL TestGPart_EndStageRKK(this,vx,vy,vz,xk,tmp1,tmp2)

  END SUBROUTINE TestGPart_StepRKK
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE TestGPart_EndStageRKK(this,vx,vy,vz,xk,tmp1,tmp2)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : EndStageRKK
!  DESCRIPTION: Called at the end of all RK-like stages to
!               complete test particle update.

!  ARGUMENTS  :
!    this    : 'this' class instance
!    vz,vy,vz: compoments of velocity field, in real space, partially
!              updated, possibly. These will be overwritten!
!    xk      : multiplicative RK time stage factor
!    tmpX    : temp arrays the same size as vx, vy, vz
!-----------------------------------------------------------------
    USE fprecision
    USE commtypes
    USE mpivars
    USE grid

    IMPLICIT NONE
    CLASS(TestGPart) ,INTENT(INOUT)                        :: this
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(nx,ny,ksta:kend) :: vx,vy,vz
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(nx,ny,ksta:kend) :: tmp1,tmp2
    REAL(KIND=GP),INTENT   (IN)                            :: xk
    INTEGER                                                :: j,ng

    ! u(t+dt) = u*: done already

    ! IF using nearest-neighbor interface, do particle exchange 
    ! between nearest-neighbor tasks BEFORE z-PERIODIZING particle coordinates.
    ! Note this interface has not been tested yet for test particles.
    IF ( this%iexchtype_.EQ.GPEXCHTYPE_NN ) THEN
      CALL GTStart(this%htimers_(GPTIME_COMM))
      CALL this%gpcomm_%PartExchangeV(this%id_,this%px_,this%py_,this%pz_,  &
           this%nparts_,this%lxbnds_(3,1),this%lxbnds_(3,2),GPEXCH_INIT)
      CALL this%gpcomm_%PartExchangeV(this%id_,this%ptmp0_(1,:),            &
           this%ptmp0_(2,:),this%ptmp0_(3,:),this%nparts_,this%lxbnds_(3,1),&
           this%lxbnds_(3,2),GPEXCH_UPDT)
      CALL this%gpcomm_%PartExchangeV(this%id_,this%ttmp0_(1,:),            &
           this%ttmp0_(2,:),this%ttmp0_(3,:),this%nparts_,                  &
           this%lxbnds_(3,1),this%lxbnds_(3,2),GPEXCH_UPDT)
      CALL this%gpcomm_%PartExchangeV(this%id_,this%pvx_,this%pvy_,this%pvz_,&
           this%nparts_,this%lxbnds_(3,1),this%lxbnds_(3,2),GPEXCH_END)
      CALL GTAcc(this%htimers_(GPTIME_COMM))
      ! Enforce periodicity in x and y:
      CALL GPart_MakePeriodicP(this,this%px_,this%py_,this%pz_,this%nparts_,3)
      ! Enforce periodicity in z and ptmp0(3):
      CALL GPart_MakePeriodicZ(this,this%pz_,this%ptmp0_(3,:),this%nparts_)

    ENDIF

    ! If using VDB interface, do synch-up, and get local work:
    IF ( this%iexchtype_.EQ.GPEXCHTYPE_VDB ) THEN

      ! Enforce periodicity in x, y, & z:
      CALL GPart_MakePeriodicP(this,this%px_,this%py_,this%pz_,this%nparts_,7)

      IF ( .NOT.GPart_PartNumConsistent(this,this%nparts_) ) THEN
        IF ( this%myrank_.eq.0 ) THEN
          WRITE(*,*) 'TestGPart_EndStepRKK: Inconsistent particle count'
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
        WRITE(*,*)'TestGPart_EndStepRKK: inconsistent d.b.: expected: ', &
                 this%maxparts_, '; found: ',ng
        STOP
      ENDIF

    ENDIF

    RETURN

  END SUBROUTINE TestGPart_EndStageRKK
!-----------------------------------------------------------------
!-----------------------------------------------------------------
