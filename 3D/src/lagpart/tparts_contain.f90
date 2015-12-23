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

    DO j = 1, this%nparts_

     u1 = rand()
     u2 = rand()
     u1 =  max( u1, low)
     u2 =  twopi * u2
     gauss =  sqrt( -2.0_GP * log( u1)) * CMPLX( cos(u2), sin(u2))

       this%pvx_(j) = vtherm*gauss

     u1 = rand()
     u2 = rand()
     u1 =  MAX( u1, low)
     u2 =  twopi * u2
     gauss =  sqrt( -2.0_GP * log( u1)) * CMPLX( cos(u2), sin(u2))

       this%pvy_(j) = vtherm*gauss

     u1 = rand()
     u2 = rand()
     u1 =  MAX( u1, low)
     u2 =  twopi * u2
     gauss =  sqrt( -2.0_GP * log( u1)) * CMPLX( cos(u2), sin(u2))

       this%pvz_(j) = vtherm*gauss

    END DO

  END SUBROUTINE TestGPart_InitVel

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

    ! Initialize solution, u: 
    ! u* <-- u: 
 
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
    CLASS(TestGPart) ,INTENT(INOUT)                      :: this
    INTEGER                                              :: i,j
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(n,n,ksta:kend) :: vx,vy,vz
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(n,n,ksta:kend) :: tmp1,tmp2,tmp3
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(n,n,ksta:kend) :: bx,by,bz
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(n,n,ksta:kend) :: jx,jy,jz
    REAL(KIND=GP),INTENT   (IN)                          :: dt,xk
    REAL(KIND=GP)                                        :: dtfact
    REAL(KIND=GP)                                        :: dtv
    REAL(KIND=GP), ALLOCATABLE, DIMENSION            (:) :: lid,gid

    dtfact = dt*xk*real(n,kind=GP)/(8.0_GP*atan(1.0_GP))
    dtv    = dt*xk
    CALL GTStart(this%htimers_(GPTIME_STEP))

    ! Find F(u*):
    CALL GPart_R3toR3(this,tmp3,vx) ! Want vx intact to use later
    CALL GPart_EulerToLag(this,this%lvx_,this%nparts_,tmp3,.true.,tmp1,tmp2)
    CALL GPart_R3toR3(this,tmp3,vy) ! Want vy intact to use later
    CALL GPart_EulerToLag(this,this%lvy_,this%nparts_,tmp3,.false.,tmp1,tmp2)
    CALL GPart_R3toR3(this,tmp3,vz) ! Want vz intact to use later
    CALL GPart_EulerToLag(this,this%lvz_,this%nparts_,tmp3,.false.,tmp1,tmp2)
    CALL GPart_EulerToLag(this,this%lbx_,this%nparts_,bx,.false.,tmp1,tmp2)
    CALL GPart_EulerToLag(this,this%lby_,this%nparts_,by,.false.,tmp1,tmp2)
    CALL GPart_EulerToLag(this,this%lbz_,this%nparts_,bz,.false.,tmp1,tmp2)
    CALL GPart_EulerToLag(this,this%lfx_,this%nparts_,jx,.false.,tmp1,tmp2)
    CALL GPart_EulerToLag(this,this%lfy_,this%nparts_,jy,.false.,tmp1,tmp2)
    CALL GPart_EulerToLag(this,this%lfz_,this%nparts_,jz,.false.,tmp1,tmp2)
!   Lorentz force
    DO j = 1, this%nparts_
       this%lfx_(j) = this%lfx_(j) + (this%pvy_(j)-this%lvy_(j))*this%lbz_(j)-(this%pvz_(j)-this%lvz_(j))*this%lby_(j)
       this%lfy_(j) = this%lfy_(j) + (this%pvz_(j)-this%lvz_(j))*this%lbx_(j)-(this%pvx_(j)-this%lvx_(j))*this%lbz_(j)
       this%lfz_(j) = this%lfz_(j) + (this%pvx_(j)-this%lvx_(j))*this%lby_(j)-(this%pvy_(j)-this%lvy_(j))*this%lbx_(j)
! Lorentz force only with magnetic field (conservative)
!      this%lfx_(j) = this%pvy_(j)*this%lbz_(j)-this%pvz_(j)*this%lby_(j)
!      this%lfy_(j) = this%pvz_(j)*this%lbx_(j)-this%pvx_(j)*this%lbz_(j)
!      this%lfz_(j) = this%pvx_(j)*this%lby_(j)-this%pvy_(j)*this%lbx_(j)
    ENDDO

    ! ... x:
    DO j = 1, this%nparts_
      this%px_(j) = this%ptmp0_(1,j) + dtfact*this%pvx_(j)
      this%pvx_(j) = this%ttmp0_(1,j) + dtv*this%lfx_(j)
    ENDDO

    ! ... y:
    DO j = 1, this%nparts_
      this%py_(j) = this%ptmp0_(2,j) + dtfact*this%pvy_(j)
      this%pvy_(j) = this%ttmp0_(2,j) + dtv*this%lfy_(j)
    ENDDO

    ! ... z:
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

    ! At this point, vx, vy, vz should be intact:
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
    CLASS(TestGPart) ,INTENT(INOUT)             :: this
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(n,n,ksta:kend) :: vx,vy,vz,tmp1,tmp2
    REAL(KIND=GP),INTENT   (IN)                 :: xk
    INTEGER                                     :: j,ng

    ! u(t+dt) = u*: done already

    ! IF using nearest-neighbor interface, do particle exchange 
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
          WRITE(*,*) 'TestGPart_EndStepRKK: Inconsistent particle count'
        ENDIF
      ENDIF
      ! Synch up VDB, if necessary:
      CALL GTStart(this%htimers_(GPTIME_COMM))
      CALL this%gpcomm_%VDBSynch(this%vdb_,this%maxparts_,this%id_, &
           this%px_,this%py_,this%pz_,this%nparts_,this%ptmp1_)
      CALL this%gpcomm_%VDBSynch(this%gptmp0_,this%maxparts_,this%id_, &
                     this%ptmp0_(1,:),this%ptmp0_(2,:),this%ptmp0_(3,:),&
                     this%nparts_,this%ptmp1_)
      CALL GPart_CopyLocalWrk(this,this%ptmp0_(1,:),  &
                        this%ptmp0_(2,:),this%ptmp0_(3,:), &
           this%vdb_,this%gptmp0_,this%maxparts_)
      CALL this%gpcomm_%VDBSynch(this%gptmp0_,this%maxparts_,this%id_, &
           this%pvx_,this%pvy_,this%pvz_,this%nparts_,this%ptmp1_)
      CALL GPart_CopyLocalWrk(this,this%pvx_,  &
           this%pvy_,this%pvz_,this%vdb_,  &
                      this%gptmp0_,this%maxparts_)
      CALL this%gpcomm_%VDBSynch(this%gptmp0_,this%maxparts_,this%id_, &
                     this%ttmp0_(1,:),this%ttmp0_(2,:),this%ttmp0_(3,:),&
                     this%nparts_,this%ptmp1_)
      CALL GPart_CopyLocalWrk(this,this%ttmp0_(1,:),   &
           this%ttmp0_(2,:),this%ttmp0_(3,:),this%vdb_, &
                  this%gptmp0_,this%maxparts_)
      CALL GPart_GetLocalWrk(this,this%id_,this%px_,this%py_,this%pz_, &
      this%nparts_, this%vdb_,this%maxparts_)
      CALL GTAcc(this%htimers_(GPTIME_COMM))

      CALL MPI_ALLREDUCE(this%nparts_,ng,1,MPI_INTEGER,   &
                         MPI_SUM,this%comm_,this%ierr_)

      IF ( this%myrank_.EQ.0 .AND. ng.NE.this%maxparts_) THEN
        WRITE(*,*)'TestGPart_EndStepRKK: inconsistent d.b.: expected: ', &
                 this%maxparts_, '; found: ',ng
!        this%vdb_ = this%vdb_*(8.0_GP*atan(1.0_GP))/real(n,kind=GP)
!        CALL GPART_ascii_write_pdb(this,1,'.','xlgerr','000',0.0_GP,this%vdb_)
!        this%vdb_ = this%vdb_/(8.0_GP*atan(1.0_GP))*real(n,kind=GP)
        STOP
      ENDIF

    ENDIF

    IF ( this%intacc_.EQ.0 ) RETURN

    ! If doing internal acceleration, synch up past time levels:
    CALL GPart_synch_acc(this)

    ! Set t^n+1 velocity based on most recent Lag.particle positions:
    ! NOTE: vx, vy, vz are overwirtten on exit:
    CALL GPart_EulerToLag(this,this%lvx_,this%nparts_,vx,.true. ,tmp1,tmp2)
    CALL GPart_EulerToLag(this,this%lvy_,this%nparts_,vy,.false.,tmp1,tmp2)
    CALL GPart_EulerToLag(this,this%lvz_,this%nparts_,vz,.false.,tmp1,tmp2)

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

  END SUBROUTINE TestGPart_EndStageRKK
!-----------------------------------------------------------------
!-----------------------------------------------------------------

    SUBROUTINE TestGPart_io_write_pdbv(this, iunit, dir, spref, nmb, time)
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
    CLASS(TestGPart) ,INTENT(INOUT)   :: this
    REAL(KIND=GP),INTENT   (IN)       :: time
    REAL(KIND=GP)                     :: prec(3)
    INTEGER,INTENT(IN)                :: iunit
    INTEGER                           :: fh,j,nt
    INTEGER(kind=MPI_OFFSET_KIND)     :: offset
    CHARACTER(len=*),INTENT(IN)       :: dir
    CHARACTER(len=*),INTENT(IN)       :: nmb
    CHARACTER(len=*),INTENT(IN)       :: spref

    IF ( this%iexchtype_.EQ.GPEXCHTYPE_NN ) THEN
       CALL this%gpcomm_%VDBSynch(this%ptmp0_,this%maxparts_,this%id_, &
            this%pvx_,this%pvy_,this%pvz_,this%nparts_,this%ptmp1_)
    ELSE
      ! Store global VDB data into temp array:
!$omp parallel do
      DO j = 1, this%maxparts_
        this%ptmp0_(1,j) = this%pvx_(1,j)
        this%ptmp0_(2,j) = this%pvy_(2,j)
        this%ptmp0_(3,j) = this%pvz_(3,j)
      ENDDO
    ENDIF

    IF ( this%iouttype_ .EQ. 0 ) THEN
      IF ( this%bcollective_.EQ. 1 ) THEN
        ! pass in the current linear _local_ particle velocity arrays
        CALL GPart_binary_write_lag_co(this,iunit,dir,spref,nmb,time,this%nparts_, &
             this%pvx_,this%pvy_,this%pvz_)
      ELSE
        ! pass in the synched-up VDB (copied to ptmp0_):
        CALL GPart_binary_write_lag_t0(this,iunit,dir,spref,nmb,time,this%maxparts_, &
             this%ptmp0_(1,:),this%ptmp0_(2,:),this%ptmp0_(3,:))
     ELSE
      ! pass in the synched-up VDB (copied to ptmp0_):
      CALL GPart_ascii_write_lag(this,iunit,dir,spref,nmb,time,this%maxparts_, &
           this%ptmp0_(1,:),this%ptmp0_(2,:),this%ptmp0_(3,:))
    ENDIF

  END SUBROUTINE TestGPart_io_write_pdbv
!-----------------------------------------------------------------
!-----------------------------------------------------------------
