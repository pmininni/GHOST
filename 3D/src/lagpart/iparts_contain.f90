!=================================================================
! VGPart SUBROUTINES
!=================================================================

  SUBROUTINE VGPart_io_write_pdbv(this, iunit, dir, spref, nmb, time)
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
    CLASS(VGPart),INTENT(INOUT)       :: this
    REAL(KIND=GP),INTENT   (IN)       :: time
    REAL(KIND=GP)                     :: prec(3)
    INTEGER,INTENT(IN)                :: iunit
    INTEGER                           :: fh,j,nt
    INTEGER(kind=MPI_OFFSET_KIND)     :: offset
    CHARACTER(len=*),INTENT(IN)       :: dir
    CHARACTER(len=*),INTENT(IN)       :: nmb
    CHARACTER(len=*),INTENT(IN)       :: spref

!    CALL this%gpcomm_%VDBSynch(this%ptmp0_,this%maxparts_,this%id_, &
!         this%pvx_,this%pvy_,this%pvz_,this%nparts_,this%ptmp1_)

    ! If doing non-collective binary or ascii writes, synch up vector:
    IF ((this%iouttype_.EQ.0.AND.this%bcollective_.EQ.0).OR.this%iouttype_.EQ.1 ) THEN
      CALL this%gpcomm_%VDBSynch_t0(this%ptmp0_,this%maxparts_,this%id_, &
                                 this%pvx_,this%pvy_,this%pvz_,this%nparts_)
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
      ENDIF
    ELSE
      ! pass in the synched-up VDB (copied to ptmp0_):
      CALL GPart_ascii_write_lag(this,iunit,dir,spref,nmb,time,this%maxparts_, &
           this%ptmp0_(1,:),this%ptmp0_(2,:),this%ptmp0_(3,:))
    ENDIF

  END SUBROUTINE VGPart_io_write_pdbv
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE VGPart_io_read_pdbv(this, iunit, dir, spref, nmb)
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
    CLASS(VGPart)   ,INTENT(INOUT)            :: this
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

    IF (this%iexchtype_ .EQ. GPEXCHTYPE_VDB) THEN
      ! Store in particle velocity arrays
      CALL GPart_CopyLocalWrk(this,this%pvx_,this%pvy_,this%pvz_, &
                            this%vdb_,this%ptmp0_,this%maxparts_)
    END IF
    CALL MPI_ALLREDUCE(this%nparts_,ng,1,MPI_INTEGER,   &
                       MPI_SUM,this%comm_,this%ierr_)
    IF ( this%myrank_.EQ.0 .AND. ng.NE.this%maxparts_ ) THEN
      WRITE(*,*)'InerGPart_io_read_pdbv: inconsistent d.b.: expected: ', &
                 this%maxparts_, '; found: ',ng
      STOP
    ENDIF

  END SUBROUTINE VGPart_io_read_pdbv
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE VGPart_ResizeArrays(this,new_size,onlyinc,exc)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : Resize_Arrays
!  DESCRIPTION: Resize all arrays in the VGPart class (including 
!               subclases, i.e. communicator, spline)
!  ARGUMENTS  :
!    this    : 'this' class instance
!    new_size: new number of particles
!    onlyinc : if true, will only resize to increase array size
!-----------------------------------------------------------------
!$  USE threads
 
    IMPLICIT NONE
    CLASS(VGPart) ,INTENT(INOUT)                         :: this
    INTEGER       ,INTENT(IN)                            :: new_size
    LOGICAL       ,INTENT(IN)                            :: onlyinc
    LOGICAL       ,INTENT(IN)   ,OPTIONAL                :: exc
    INTEGER                                              :: n

    CALL GPart_ResizeArrays(this,new_size,onlyinc,exc)

    n = SIZE(this%pvx_)
    IF ((n.lt.new_size).OR.((n.gt.new_size).AND..NOT.onlyinc)) THEN
      CALL Resize_ArrayRank1(this%pvx_,new_size,.true.)
    END IF
    n = SIZE(this%pvy_)
    IF ((n.lt.new_size).OR.((n.gt.new_size).AND..NOT.onlyinc)) THEN
      CALL Resize_ArrayRank1(this%pvy_,new_size,.true.)
    END IF
    n = SIZE(this%pvz_)
    IF ((n.lt.new_size).OR.((n.gt.new_size).AND..NOT.onlyinc)) THEN
      CALL Resize_ArrayRank1(this%pvz_,new_size,.true.)
    END IF

    n = SIZE(this%ttmp0_,2)
    IF ((n.lt.new_size).OR.((n.gt.new_size).AND..NOT.onlyinc)) THEN
      CALL Resize_ArrayRank2(this%ttmp0_,new_size,.true.)
    END IF

    RETURN

  END SUBROUTINE VGPart_ResizeArrays
!-----------------------------------------------------------------
!-----------------------------------------------------------------

!=================================================================
! InerGPart SUBROUTINES
!=================================================================

  SUBROUTINE InerGPart_ctor(this,tau,grav,gamma,nu,donldrag,om,x0)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  Explicit constructor for inertial particles. Should be called
!  after calling GPart_ctor.
!
!  ARGUMENTS:
!    this    : 'this' class instance
!    tau     : Stokes time
!    grav    : gravity acceleration
!    gamma   : mass ratio (= m_f/m_p)
!    nu      : fluid viscosity
!-----------------------------------------------------------------
    USE fprecision

    IMPLICIT NONE
    CLASS(InerGPart), INTENT(INOUT)     :: this
    REAL(KIND=GP),INTENT(IN)            :: tau,grav,gamma,nu
    REAL(KIND=GP),INTENT(IN),OPTIONAL   :: om(3),x0(3)
    INTEGER,      INTENT(IN)            :: donldrag

    this%tau_      = tau
    this%invtau_   = 1.0_GP/tau
    this%grav_     = grav
    this%gamma_    = gamma
    this%nu_       = nu
    this%donldrag_ = donldrag
    IF (PRESENT(om)) THEN
       this%dorotatn_ = 1
       this%omegax_ = om(1)
       this%omegay_ = om(2)
       this%omegaz_ = om(3)
       this%px0_ = x0(1)
       this%py0_ = x0(2)
       this%pz0_ = x0(3)
    ELSE
       this%dorotatn_ = 0
    ENDIF       

    ALLOCATE(this%pvx_     (this%partbuff_))
    ALLOCATE(this%pvy_     (this%partbuff_))
    ALLOCATE(this%pvz_     (this%partbuff_))
    ALLOCATE(this%dfx_     (this%partbuff_))
    ALLOCATE(this%dfy_     (this%partbuff_))
    ALLOCATE(this%dfz_     (this%partbuff_))
    ALLOCATE(this%ttmp0_ (3,this%partbuff_))

  END SUBROUTINE InerGPart_ctor
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE InerGPart_dtor(this)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  Explicit destructor for inertial particles. Should be called
!  after GPart_dtor.
!
!  ARGUMENTS:
!    this    : 'this' class instance
!-----------------------------------------------------------------

    IMPLICIT NONE
    TYPE(InerGPart)   ,INTENT(INOUT)    :: this

    IF ( ALLOCATED   (this%pvx_) ) DEALLOCATE   (this%pvx_)
    IF ( ALLOCATED   (this%pvy_) ) DEALLOCATE   (this%pvy_)
    IF ( ALLOCATED   (this%pvz_) ) DEALLOCATE   (this%pvz_)
    IF ( ALLOCATED   (this%dfx_) ) DEALLOCATE   (this%dfx_)
    IF ( ALLOCATED   (this%dfy_) ) DEALLOCATE   (this%dfy_)
    IF ( ALLOCATED   (this%dfz_) ) DEALLOCATE   (this%dfz_)
    IF ( ALLOCATED (this%ttmp0_) ) DEALLOCATE (this%ttmp0_)

  END SUBROUTINE InerGPart_dtor
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE InerGPart_InitVel(this,vx,vy,vz,tmp1,tmp2)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : InitVel
!  DESCRIPTION: Initializes particle velocities with fluid
!               velocities. Other parameters are initialized
!               with GPart_Init.
!  ARGUMENTS:
!    this    : 'this' class instance
!    vz,vy,vz: compoments of velocity field, in real space, partially
!              updated, possibly. These will be overwritten!
!    tmpX    : temp arrays the same size as vx, vy, vz
!-----------------------------------------------------------------
    USE grid
    USE fprecision
    USE commtypes
    USE mpivars
    
    IMPLICIT NONE
    CLASS(InerGPart),INTENT(INOUT)                         :: this
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(nx,ny,ksta:kend) :: vx,vy,vz
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(nx,ny,ksta:kend) :: tmp1,tmp2
    
    CALL GPart_EulerToLag(this,this%pvx_,this%nparts_,vx,.true. ,tmp1,tmp2)
    CALL GPart_EulerToLag(this,this%pvy_,this%nparts_,vy,.false.,tmp1,tmp2)
    CALL GPart_EulerToLag(this,this%pvz_,this%nparts_,vz,.false.,tmp1,tmp2)

  END SUBROUTINE InerGPart_InitVel
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE InerGPart_SetStepRKK(this)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : SetStepRKK
!  DESCRIPTION: Initializes an explicit integration timestep for
!               the inertial velocity. Must be called at the start
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
    CLASS(InerGPart), INTENT(INOUT) :: this
    INTEGER                          :: j

    ! Initialize solution, v (particle velocity): 
    ! v* <-- v: 
 
    ! Cycle over JST loop to update state:
!$omp parallel do
    DO j = 1, this%nparts_
       this%ttmp0_(1,j) = this%pvx_(j)  ! vx_0
       this%ttmp0_(2,j) = this%pvy_(j)  ! vy_0
       this%ttmp0_(3,j) = this%pvz_(j)  ! vz_0
    ENDDO
    
  END SUBROUTINE InerGPart_SetStepRKK
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE InerGPart_StepRKK(this, vx, vy, vz, dt, xk, tmp1, tmp2)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : Step_testp
!  DESCRIPTION: Carries out one stage of explicit RK-like time
!               integration step.  Intended for explicit step within 
!               an outer stepper method of the form:
!
!               X = X_0 + dt * V[X(t),t] * xk,
!               V = V_0 + dt * (F[V(X(t)),U(X(t))] - g_z) * xk,
!       
!               where F is the drag force, V(X(t)) is the particle
!               velocity, U(X(t)) is the Lagrangian velocity, and
!               g_z is the (positive) z-component of the gravity  
!               acceleration. The drag force is:
!
!               F = 1/tau ( U(X(t)) - V(X(t)) )
!
!               Inertial particles in this method are pointwise and
!               heavy. Any other forces, except for the Stokes drag,
!               are neglected.
!
!               Note that the vx, vy, vz, will be overwritten here.
!  ARGUMENTS  :
!    this     : 'this' class instance
!    vz,vy,vz : compoments of velocity field, in real space, partially
!               updated, possibly. These will be overwritten!
!    dt       : integration timestep
!    xk       : multiplicative RK time stage factor
!    tmpX     : temp arrays the same size as vx, vy, vz
!-----------------------------------------------------------------
    USE grid
    USE fprecision
    USE commtypes
    USE mpivars

    IMPLICIT NONE
    CLASS(InerGPart) ,INTENT(INOUT)                        :: this
    INTEGER                                                :: i,j
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(nx,ny,ksta:kend) :: vx,vy,vz
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(nx,ny,ksta:kend) :: tmp1,tmp2
    REAL(KIND=GP),INTENT   (IN)                            :: dt,xk
    REAL(KIND=GP)                                          :: dtfact
    REAL(KIND=GP)                                          :: dtv
    REAL(KIND=GP), ALLOCATABLE, DIMENSION              (:) :: lid,gid

    dtv    = dt*xk
    CALL GTStart(this%htimers_(GPTIME_STEP))

    ! Find the Lagrangian velocity for F(U*,V*):
    CALL GPart_EulerToLag(this,this%lvx_,this%nparts_,vx,.true.,tmp1,tmp2)
    CALL GPart_EulerToLag(this,this%lvy_,this%nparts_,vy,.false.,tmp1,tmp2)
    CALL GPart_EulerToLag(this,this%lvz_,this%nparts_,vz,.false.,tmp1,tmp2)

    ! Drag force
!$omp parallel do
    DO j = 1, this%nparts_
       this%dfx_(j) = (this%lvx_(j)-this%pvx_(j))*this%invtau_
       this%dfy_(j) = (this%lvy_(j)-this%pvy_(j))*this%invtau_
       this%dfz_(j) = (this%lvz_(j)-this%pvz_(j))*this%invtau_
    ENDDO

    ! ... x:
    dtfact = dt*xk*this%invdel_(1)
!$omp parallel do
    DO j = 1, this%nparts_
      this%px_(j) = this%ptmp0_(1,j) + dtfact*this%pvx_(j)
      this%pvx_(j) = this%ttmp0_(1,j) + dtv*this%dfx_(j)
    ENDDO

    ! ... y:
    dtfact = dt*xk*this%invdel_(2)
!$omp parallel do
    DO j = 1, this%nparts_
      this%py_(j) = this%ptmp0_(2,j) + dtfact*this%pvy_(j)
      this%pvy_(j) = this%ttmp0_(2,j) + dtv*this%dfy_(j)
    ENDDO

    ! ... z:
    dtfact = dt*xk*this%invdel_(3)
!$omp parallel do
    DO j = 1, this%nparts_
      this%pz_(j) = this%ptmp0_(3,j) + dtfact*this%pvz_(j)
      this%pvz_(j) = this%ttmp0_(3,j) + dtv*(this%dfz_(j)-this%grav_)
    ENDDO

    ! Enforce periodicity in x-y only:
    CALL GPart_MakePeriodicP(this,this%px_,this%py_,this%pz_,this%nparts_,3)

    CALL GTAcc(this%htimers_(GPTIME_STEP))

    CALL InerGPart_EndStageRKK(this,vx,vy,vz,xk,tmp1,tmp2)

  END SUBROUTINE InerGPart_StepRKK
!-----------------------------------------------------------------
!-----------------------------------------------------------------

 SUBROUTINE InerGPart_pnlt_StepRKK(this, vx, vy, vz, dt, xk, tmp1, tmp2, x0, y0, z0 , R, shape, n)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : Step_testp
!  DESCRIPTION: Carries out one stage of explicit RK-like time
!               integration step.  Intended for explicit step within 
!               an outer stepper method of the form:
!
!               X = X_0 + dt * V[X(t),t] * xk,
!               V = V_0 + dt * (F[V(X(t)),U(X(t))] - g_z) * xk,
!       
!               where F is the drag force, V(X(t)) is the particle
!               velocity, U(X(t)) is the Lagrangian velocity, and
!               g_z is the (positive) z-component of the gravity  
!               acceleration. The drag force is:
!
!               F = 1/tau ( U(X(t)) - V(X(t)) )
!
!               Inertial particles in this method are pointwise and
!               heavy. Any other forces, except for the Stokes drag,
!               are neglected. This method also computes collitions
!               between the particles and a body submerged in the 
!               fluid using the penalty method.
!
!               Note that the vx, vy, vz, will be overwritten here.
!  ARGUMENTS  :
!    this     : 'this' class instance
!    vx,vy,vz : compoments of velocity field, in real space, partially
!               updated, possibly. These will be overwritten!
!    dt       : integration timestep
!    xk       : multiplicative RK time stage factor
!    tmpX     : temp arrays the same size as vx, vy, vz
!    xc,yc,zc : coordinates of the center of the obstacle
!    R        : radius of the obstacle
!    shape    : shape of the body (only spheres are supported now)
!    n        : absolute time index in the simulation
!-----------------------------------------------------------------
    USE grid
    USE fprecision
    USE commtypes
    USE mpivars

    IMPLICIT NONE
    CLASS(InerGPart) ,INTENT(INOUT)                        :: this
    INTEGER                                                :: i,j
    INTEGER,INTENT         (IN)                            :: n,shape
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(nx,ny,ksta:kend) :: vx,vy,vz
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(nx,ny,ksta:kend) :: tmp1,tmp2
    REAL(KIND=GP),INTENT   (IN)                            :: dt,xk
    REAL(KIND=GP),INTENT   (IN)                            :: x0,y0,z0,R
    REAL(KIND=GP)                                          :: norm_x, norm_y, norm_z, pi, collision_count
    REAL(KIND=GP)                                          :: dtv, nfact, distance, a,b,c, xc,yc,zc, R0, dt_frac
    REAL(KIND=GP)                                          :: x_frac, y_frac, z_frac, vz_temp, size_fact
    REAL(KIND=GP)                                          :: vx_temp, vy_temp, dtfact_x, dtfact_y, dtfact_z
    REAL(KIND=GP), ALLOCATABLE, DIMENSION              (:) :: lid, gid
    
    collision_count = 0 

    dtv    = dt*xk
    CALL GTStart(this%htimers_(GPTIME_STEP))

    ! Find the Lagrangian velocity for F(U*,V*):
    CALL GPart_EulerToLag(this,this%lvx_,this%nparts_,vx,.true.,tmp1,tmp2)
    CALL GPart_EulerToLag(this,this%lvy_,this%nparts_,vy,.false.,tmp1,tmp2)
    CALL GPart_EulerToLag(this,this%lvz_,this%nparts_,vz,.false.,tmp1,tmp2)

    ! Drag force
!$omp parallel do
    DO j = 1, this%nparts_
       this%dfx_(j) = (this%lvx_(j)-this%pvx_(j))*this%invtau_
       this%dfy_(j) = (this%lvy_(j)-this%pvy_(j))*this%invtau_
       this%dfz_(j) = (this%lvz_(j)-this%pvz_(j))*this%invtau_
    ENDDO

! Reflection equations for the particles
   ! Particles parameters
    pi = 4.0_GP*atan(1.0_GP)          
    xc = x0*nx                   ! x coordinate of sphere center in grid units
    yc = y0*ny                   ! y coordinate of sphere center in grid units
    zc = z0*nz                   ! z coordinate of sphere center in grid units
    R0 = R*((ny/2)/pi)           ! Radius of the sphere 
    size_fact = 0.05_GP          ! Ratio of particle to sphere radius (r/R0)

    R0 = R0 * (1 + size_fact)

    ! ... x, y and z:
    dtfact_x = dt*xk*this%invdel_(1)
    dtfact_y = dt*xk*this%invdel_(2)
    dtfact_z = dt*xk*this%invdel_(3)
	
    !$omp parallel do
    DO j = 1, this%nparts_
	
    ! Define the temporary velocties and normal vectors for radial distance calculation	
	vx_temp = this%pvx_(j)
	vy_temp = this%pvy_(j)
	vz_temp = this%pvz_(j)

  	norm_x = (this%ptmp0_(1,j) - xc)/R0
	norm_y = (this%ptmp0_(2,j) - yc)/R0
	norm_z = (this%ptmp0_(3,j) - zc)/R0

    ! Find hypothetical radial distance of particle if it continues along current trajectory for 1 dt
	distance = sqrt((this%ptmp0_(1,j) + dt*xk*this%invdel_(1)*vx_temp - xc)**2 &
			+ (this%ptmp0_(2,j) + dt*xk*this%invdel_(2)*vy_temp - yc)**2 &
			+ (this%ptmp0_(3,j) + dt*xk*this%invdel_(3)*vz_temp - zc)**2)
		
    ! Check for collision condition
	IF ((distance.lt.R0) .AND. ((vx_temp*norm_x + vy_temp*norm_y + vz_temp*norm_z).lt.0))  THEN	
		
    ! Find coefficients of the quadratic 		
		a = (vx_temp*xk*this%invdel_(1))**2 & 
		   +(vy_temp*xk*this%invdel_(2))**2 & 
		   +(vz_temp*xk*this%invdel_(3))**2
		
		b = 2*((this%ptmp0_(1,j) - xc)*vx_temp*xk*this%invdel_(1) &
		     + (this%ptmp0_(2,j) - yc)*vy_temp*xk*this%invdel_(2) &
		     + (this%ptmp0_(3,j) - zc)*vz_temp*xk*this%invdel_(3))        	
		
		c = ((xc-this%ptmp0_(1, j))**2 + (yc-this%ptmp0_(2, j))**2 &
                + (zc-this%ptmp0_(3, j))**2 - R0**2)
		
		! Check for positive solution for fractional timestep		
		IF ((((-b+sqrt(b**2 - 4*a*c))/(2*a)).gt.0) .AND. (((-b+sqrt(b**2 - 4*a*c))/(2*a)).lt.dt)) THEN
      			dt_frac = (-b+sqrt(b**2 - 4*a*c))/(2*a)
    		ELSE
      			dt_frac = (-b-sqrt(b**2 - 4*a*c))/(2*a)
		ENDIF

    ! Determine the fractional position increment
    		x_frac = this%ptmp0_(1,j) + this%ttmp0_(1,j)*dt_frac*xk*this%invdel_(1)  
		y_frac = this%ptmp0_(2,j) + this%ttmp0_(2,j)*dt_frac*xk*this%invdel_(2)   
		z_frac = this%ptmp0_(3,j) + this%ttmp0_(3,j)*dt_frac*xk*this%invdel_(3)   		

    ! Update v_temp as the temporary velocities 
    	        vx_temp = this%ttmp0_(1,j)
        	vy_temp = this%ttmp0_(2,j)
        	vz_temp = this%ttmp0_(3,j)

    ! Calculate the normal vectors at the time of contact
		norm_x = (x_frac - xc)/R0
		norm_y = (y_frac - yc)/R0
		norm_z = (z_frac - zc)/R0

		! Define the common factor (v.n^)
		nfact = vx_temp*norm_x + vy_temp*norm_y + vz_temp*norm_z 
		
		! Update the temporary velocities after collision
		this%ttmp0_(1,j) = this%ttmp0_(1,j) - 2.0*nfact*norm_x
		this%ttmp0_(2,j) = this%ttmp0_(2,j) - 2.0*nfact*norm_y
		this%ttmp0_(3,j) = this%ttmp0_(3,j) - 2.0*nfact*norm_z
		
    ! Update the post collision positions
		this%px_(j) = x_frac + (dt - dt_frac)*this%pvx_(j)*xk*this%invdel_(1)
		this%py_(j) = y_frac + (dt - dt_frac)*this%pvy_(j)*xk*this%invdel_(2)
		this%pz_(j) = z_frac + (dt - dt_frac)*this%pvz_(j)*xk*this%invdel_(3)	
		
    ! Update the post collision velocities 
    		this%pvx_(j) = this%ttmp0_(1,j) 
	  	this%pvy_(j) = this%ttmp0_(2,j) 
		this%pvz_(j) = this%ttmp0_(3,j) 

		collision_count = collision_count + 1

    ! We write data from individual collisions
    	        OPEN(1,file='collisions_individual.txt', position='append')
		  WRITE(1, FMT='(E13.6, I13, E26.18, E26.18, E26.18, E26.18, E26.18, E26.18, &
                   E26.18, E26.18, E26.18, E26.18, E26.18, E26.18, &
                   E26.18, E26.18, E26.18, E26.18)') n*dt, this%id_(j), this%ptmp0_(1,j), &
                                           this%ptmp0_(2,j), this%ptmp0_(3,j), vx_temp, vy_temp, vz_temp, &
                                           this%px_(j), this%py_(j), this%pz_(j), this%pvx_(j), this%pvy_(j), this%pvz_(j), &
                                           dt_frac, x_frac, y_frac, z_frac
	        CLOSE(1)

	ELSE

    ! Normal particle evolution without collisions
		this%px_(j) = this%ptmp0_(1,j) + dtfact_x*this%pvx_(j)
		this%py_(j) = this%ptmp0_(2,j) + dtfact_y*this%pvy_(j)
		this%pz_(j) = this%ptmp0_(3,j) + dtfact_z*this%pvz_(j)

    		this%pvx_(j) = this%ttmp0_(1,j) + dtv*this%dfx_(j)
	  	this%pvy_(j) = this%ttmp0_(2,j) + dtv*this%dfy_(j)
		this%pvz_(j) = this%ttmp0_(3,j) + dtv*(this%dfz_(j) - this%grav_)

	ENDIF	            
    ENDDO

    ! We keep track of the collision count in an external file
    IF (collision_count.gt.0) THEN
      OPEN(1,file='collisions_total.txt', position='append')
      WRITE(1, FMT='(E13.6, E26.18)') n*dt, collision_count
      CLOSE(1)
    ENDIF

    ! Enforce periodicity in x-y only:
    CALL GPart_MakePeriodicP(this,this%px_,this%py_,this%pz_,this%nparts_,3)

    CALL GTAcc(this%htimers_(GPTIME_STEP))

    CALL InerGPart_EndStageRKK(this,vx,vy,vz,xk,tmp1,tmp2)

  END SUBROUTINE InerGPart_pnlt_StepRKK
!-----------------------------------------------------------------
!-----------------------------------------------------------------
  
  SUBROUTINE InerGPart_lite_StepRKK(this, vx, vy, vz, ax, ay, az, dt, xk, tmp1, tmp2)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : Step_testp
!  DESCRIPTION: Carries out one stage of explicit RK-like time
!               integration step.  Intended for explicit step within 
!               an outer stepper method of the form:
!
!               X = X_0 + dt * V[X(t),t] * xk,
!               V = V_0 + dt * (F[V(X(t)),U(X(t))] - G_z + 3/2 R DU/Dt) * xk,
!       
!               where F is the drag force, V(X(t)) is the particle
!               velocity, U(X(t)) is the Lagrangian velocity, G_z is
!               the z-component of the corrected gravity acceleration  
!               (= g*(1-gamma)/(1+gamma/2), with g>0 and gamma=m_f/m_p,
!               where m_f is the fluid mass and m_p the particle mass),
!               DU/dt is the fluid Lagrangian acceleration, and
!               R=gamma/(1+gamma/2).
!               The drag force is:
!
!               F = 1/tau ( U(X(t)) - V(X(t)) )
!
!               or
!
!               F = (1 + 0.15 Re_p^0.687) /tau ( U(X(t)) - V(X(t)) )
!
!               if nonlinear drag is used (see Wang & Maxey 1993),
!               where Re_p = (18 tau gamma/nu)^(1/2) |U - V|, and
!               |U - V| = [(Ux-Vx)^2+(Uy-Vy)^2+(Uz-Vz)^2]^(1/2).
!               If rotation is enabled in the fluid solver, the equation
!               for the velocity also includes the Coriolis force:
!
!               - 2 Omega x [V(X(t)) - 3/2 R U(X(t))]
!
!               and centrifugal force:
!
!               - (1 - 3/2 R) Omega x Omega x [X(t) - X0] 
!
!               where X0 is located in the center of the domain.    
!               This method is intended for light inertial particles, or
!               for inertial particles that include all terms in the
!               Maxey-Riley equations to first order in the particle radius.
!               Note that the vx, vy, vz, will be overwritten here.
!  ARGUMENTS  :
!    this     : 'this' class instance
!    vz,vy,vz : compoments of velocity field, in real space, partially
!               updated, possibly. These will be overwritten!
!    az,ay,az : compoments of Lagrangian acceleration, in real space,
!               partially updated, possibly. These may be overwritten!
!    dt       : integration timestep
!    xk       : multiplicative RK time stage factor
!    tmpX     : temp arrays the same size as vx, vy, vz
!-----------------------------------------------------------------
    USE grid
    USE fprecision
    USE commtypes
    USE mpivars

    IMPLICIT NONE
    CLASS(InerGPart) ,INTENT(INOUT)                        :: this
    INTEGER                                                :: i,j
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(nx,ny,ksta:kend) :: vx,vy,vz
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(nx,ny,ksta:kend) :: ax,ay,az
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(nx,ny,ksta:kend) :: tmp1,tmp2
    REAL(KIND=GP),INTENT   (IN)                            :: dt,xk
    REAL(KIND=GP)                                          :: dtfact
    REAL(KIND=GP)                                          :: dtv
    REAL(KIND=GP)                                          :: rep,cdrag
    REAL(KIND=GP)                                          :: tmparg
    REAL(KIND=GP), ALLOCATABLE, DIMENSION              (:) :: lid,gid

    dtv    = dt*xk
    CALL GTStart(this%htimers_(GPTIME_STEP))

    ! Find the Lagrangian velocity for F(U,V):
    CALL GPart_EulerToLag(this,this%lvx_,this%nparts_,vx,.true.,tmp1,tmp2)
    CALL GPart_EulerToLag(this,this%lvy_,this%nparts_,vy,.false.,tmp1,tmp2)
    CALL GPart_EulerToLag(this,this%lvz_,this%nparts_,vz,.false.,tmp1,tmp2)

    ! Find the Lagrangian acceleration of the fluid:
    CALL GPart_EulerToLag(this,this%dfx_,this%nparts_,ax,.false.,tmp1,tmp2)
    CALL GPart_EulerToLag(this,this%dfy_,this%nparts_,ay,.false.,tmp1,tmp2)
    CALL GPart_EulerToLag(this,this%dfz_,this%nparts_,az,.false.,tmp1,tmp2)

    ! Drag force plus mass ratio term

    tmparg = 1.5_GP*this%gamma_/(1.0_GP+0.5_GP*this%gamma_)
    IF ( this%donldrag_.EQ.0 ) THEN ! Linear drag
!$omp parallel do
    DO j = 1, this%nparts_
       this%dfx_(j) = this%dfx_(j)*tmparg+(this%lvx_(j)-this%pvx_(j))*this%invtau_
       this%dfy_(j) = this%dfy_(j)*tmparg+(this%lvy_(j)-this%pvy_(j))*this%invtau_
       this%dfz_(j) = this%dfz_(j)*tmparg+(this%lvz_(j)-this%pvz_(j))*this%invtau_
    ENDDO
    ELSE IF ( this%donldrag_.EQ.1 ) THEN ! Nonlinear drag
    rep = .15_GP*(18.0_GP*this%tau_*this%gamma_/this%nu_)**0.3435 ! 0.3435 = 0.687/2
!$omp parallel do
    DO j = 1, this%nparts_
       cdrag = 1.0_GP + rep*((this%lvx_(j)-this%pvx_(j))**2 + &
               (this%lvy_(j)-this%pvy_(j))**2 + (this%lvz_(j)-this%pvz_(j))**2)**0.3435
       this%dfx_(j) = this%dfx_(j)*tmparg+(this%lvx_(j)-this%pvx_(j))*cdrag*this%invtau_
       this%dfy_(j) = this%dfy_(j)*tmparg+(this%lvy_(j)-this%pvy_(j))*cdrag*this%invtau_
       this%dfz_(j) = this%dfz_(j)*tmparg+(this%lvz_(j)-this%pvz_(j))*cdrag*this%invtau_
    ENDDO
    ENDIF

    ! Rotation
    IF ( this%dorotatn_.EQ.1 ) THEN
    DO j = 1, this%nparts_
       this%dfx_(j) = this%dfx_(j)-2.0_GP*(this%omegay_*(this%pvz_(j)-tmparg*this%lvz_(j))  - &
                                           this%omegaz_*(this%pvy_(j)-tmparg*this%lvy_(j))) - &
        (1.0_GP-tmparg)*(this%omegay_*(this%omegax_*this%delta_(2)*(this%py_(j)-this%py0_)  - &
                                       this%omegay_*this%delta_(1)*(this%px_(j)-this%px0_)) - &
                         this%omegaz_*(this%omegaz_*this%delta_(1)*(this%px_(j)-this%px0_)  - &
                                       this%omegax_*this%delta_(3)*(this%pz_(j)-this%pz0_)))
       this%dfy_(j) = this%dfy_(j)-2.0_GP*(this%omegaz_*(this%pvx_(j)-tmparg*this%lvx_(j))  - &
                                           this%omegax_*(this%pvz_(j)-tmparg*this%lvz_(j))) - &
        (1.0_GP-tmparg)*(this%omegaz_*(this%omegay_*this%delta_(3)*(this%pz_(j)-this%pz0_)  - &
                                       this%omegaz_*this%delta_(2)*(this%py_(j)-this%py0_)) - &
                         this%omegax_*(this%omegax_*this%delta_(2)*(this%py_(j)-this%py0_)  - &
                                       this%omegay_*this%delta_(1)*(this%px_(j)-this%px0_)))
       this%dfz_(j) = this%dfz_(j)-2.0_GP*(this%omegax_*(this%pvy_(j)-tmparg*this%lvy_(j))  - &
                                           this%omegay_*(this%pvx_(j)-tmparg*this%lvx_(j))) - &
        (1.0_GP-tmparg)*(this%omegax_*(this%omegaz_*this%delta_(1)*(this%px_(j)-this%px0_)  - &
                                       this%omegax_*this%delta_(3)*(this%pz_(j)-this%pz0_)) - &
                         this%omegay_*(this%omegay_*this%delta_(3)*(this%pz_(j)-this%pz0_)  - &
                                       this%omegaz_*this%delta_(2)*(this%py_(j)-this%py0_)))
    ENDDO
    ENDIF

    ! ... x:
    dtfact = dt*xk*this%invdel_(1)
!$omp parallel do
    DO j = 1, this%nparts_
      this%px_(j) = this%ptmp0_(1,j) + dtfact*this%pvx_(j)
      this%pvx_(j) = this%ttmp0_(1,j) + dtv*this%dfx_(j)
    ENDDO

    ! ... y:
    dtfact = dt*xk*this%invdel_(2)
!$omp parallel do
    DO j = 1, this%nparts_
      this%py_(j) = this%ptmp0_(2,j) + dtfact*this%pvy_(j)
      this%pvy_(j) = this%ttmp0_(2,j) + dtv*this%dfy_(j)
    ENDDO

    ! ... z:
    dtfact = dt*xk*this%invdel_(3)
    tmparg = this%grav_*(1.0_GP-this%gamma_)/(1.0_GP+0.5_GP*this%gamma_)
!$omp parallel do
    DO j = 1, this%nparts_
      this%pz_(j) = this%ptmp0_(3,j) + dtfact*this%pvz_(j)
      this%pvz_(j) = this%ttmp0_(3,j) + dtv*(this%dfz_(j)-tmparg)
    ENDDO

    ! Enforce periodicity in x-y only:
    CALL GPart_MakePeriodicP(this,this%px_,this%py_,this%pz_,this%nparts_,3)

    CALL GTAcc(this%htimers_(GPTIME_STEP))

    CALL InerGPart_EndStageRKK(this,vx,vy,vz,xk,tmp1,tmp2)

  END SUBROUTINE InerGPart_lite_StepRKK
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE InerGPart_EndStageRKK(this,vx,vy,vz,xk,tmp1,tmp2)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : EndStageRKK
!  DESCRIPTION: Called at the end of all RK-like stages to
!               complete inertial particle update.

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
    CLASS(InerGPart) ,INTENT(INOUT)                        :: this
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
      CALL this%gpcomm_%IdentifyExchV(this%id_,this%pz_,this%nparts_,ng,    &
                                     this%lxbnds_(3,1),this%lxbnds_(3,2))
      IF (ng.GT.this%partbuff_) THEN
        PRINT *, 'Rank', this%myrank_, 'resizing: nparts=', ng, ' | partbuff=',&
                 this%partbuff_, ' --> ', this%partbuff_ + &
                 (1+(ng-this%partbuff_)/this%partchunksize_)*this%partchunksize_
        this%partbuff_ = this%partbuff_ + &
              (1+(ng-this%partbuff_)/this%partchunksize_)*this%partchunksize_
        CALL this%ResizeArrays(this%partbuff_,.true.)
      END IF
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
      IF (this%stepcounter_.GE.GPSWIPERATE) THEN
        IF ((this%bcollective_.EQ.1).OR.(this%myrank_.NE.0)) THEN
          ng = this%partbuff_ - this%nparts_
          IF (ng.GE.this%partchunksize_) THEN   ! Reduce array size
            PRINT *, 'Rank', this%myrank_, 'resizing: nparts=', this%nparts_, '| partbuff=',&
                      this%partbuff_, ' --> ', this%partbuff_ - &
                      (ng/this%partchunksize_-1)*this%partchunksize_
            this%partbuff_ = this%partbuff_ - &
                         (ng/this%partchunksize_-1)*this%partchunksize_
            CALL this%ResizeArrays(this%partbuff_,.false.)
          END IF
        END IF
        this%stepcounter_ = 1
      ELSE
        this%stepcounter_ = this%stepcounter_ + 1
      END IF
    ENDIF

    ! If using VDB interface, do synch-up, and get local work:
    IF ( this%iexchtype_.EQ.GPEXCHTYPE_VDB ) THEN

      ! Enforce periodicity in x, y, & z:
      CALL GPart_MakePeriodicP(this,this%px_,this%py_,this%pz_,this%nparts_,7)

      IF ( .NOT.GPart_PartNumConsistent(this,this%nparts_) ) THEN 
        IF ( this%myrank_.eq.0 ) THEN
          WRITE(*,*) 'InerGPart_EndStepRKK: Inconsistent particle count'
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
        WRITE(*,*)'InerGPart_EndStepRKK: inconsistent d.b.: expected: ', &
                 this%maxparts_, '; found: ',ng
        STOP
      ENDIF

    ENDIF

    RETURN

  END SUBROUTINE InerGPart_EndStageRKK
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE InerGPart_ResizeArrays(this,new_size,onlyinc,exc)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : Resize_Arrays
!  DESCRIPTION: Resize all arrays in the GPart class (including 
!               subclases, i.e. communicator, spline)
!  ARGUMENTS  :
!    this    : 'this' class instance
!    new_size: new number of particles
!    onlyinc : if true, will only resize to increase array size
!-----------------------------------------------------------------
!$  USE threads
 
    IMPLICIT NONE
    CLASS(InerGPart) ,INTENT(INOUT)                      :: this
    INTEGER          ,INTENT(IN)                         :: new_size
    LOGICAL          ,INTENT(IN)                         :: onlyinc
    LOGICAL          ,INTENT(IN)   ,OPTIONAL             :: exc
    INTEGER                                              :: n

    CALL VGPart_ResizeArrays(this,new_size,onlyinc,exc)

    n = SIZE(this%dfx_)
    IF ((n.lt.new_size).OR.((n.gt.new_size).AND..NOT.onlyinc)) THEN
      CALL Resize_ArrayRank1(this%dfx_,new_size,.false.)
    END IF
    n = SIZE(this%dfy_)
    IF ((n.lt.new_size).OR.((n.gt.new_size).AND..NOT.onlyinc)) THEN
      CALL Resize_ArrayRank1(this%dfy_,new_size,.false.)
    END IF
    n = SIZE(this%dfz_)
    IF ((n.lt.new_size).OR.((n.gt.new_size).AND..NOT.onlyinc)) THEN
      CALL Resize_ArrayRank1(this%dfz_,new_size,.false.)
    END IF

    RETURN

  END SUBROUTINE InerGPart_ResizeArrays
!-----------------------------------------------------------------
!-----------------------------------------------------------------
