!=================================================================
! GFieldComm SUBROUTINES
!=================================================================

  SUBROUTINE GFieldComm_ctor(this)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  Explicit constructor for field communicator. Should be called
!  after calling GPartComm_ctor.
!
!  ARGUMENTS:
!    this    : 'this' class instance
!-----------------------------------------------------------------
    IMPLICIT NONE
    CLASS(GFieldComm), INTENT(INOUT)      :: this
    INTEGER                               :: jf,kg,k,n2p,nt
    INTEGER,ALLOCATABLE,DIMENSION(:)      :: jfwd,kfend,kfsta,nzf

    n2p = this%nd_(3)/this%nprocs_
    nt  = (this%nzghost_+n2p-1)/n2p  ! max no. tasks needed for ghost zones

    ALLOCATE(jfwd (nt))
    ALLOCATE(kfend(nt))
    ALLOCATE(kfsta(nt))
    ALLOCATE(nzf  (nt))

    ALLOCATE(this%ibretp_(nt))
    ALLOCATE(this%itretp_(nt))
    ALLOCATE(this%ibretnz_(nt))
    ALLOCATE(this%itretnz_(nt))
    ALLOCATE(this%ibret_(nt,this%nzghost_+1))
    ALLOCATE(this%itret_(nt,this%nzghost_+1))
    ALLOCATE(this%ibretdst_(nt,this%nzghost_+1))
    ALLOCATE(this%itretdst_(nt,this%nzghost_+1))

    ! *** Find top neighbors to return to:
    DO jf = 1, nt
      jfwd(jf) = modulo(this%myrank_+jf,this%nprocs_)
      CALL range(1,this%nd_(3),this%nprocs_,jfwd(jf),kfsta(jf),kfend(jf))
      nzf (jf) = kfend(jf) - kfsta(jf) + 1
    ENDDO

    this%ntret_ = 1
    this%itretp_(this%ntret_) = jfwd(this%ntret_)
    this%itretnz_(this%ntret_) = 0
    k = 1
    DO kg = 1, this%nzghost_
      IF ( k .GT. nzf(this%ntret_) ) THEN
        k = k - nzf(this%ntret_)
        this%ntret_ = this%ntret_ + 1
        this%itretp_(this%ntret_)  = jfwd(this%ntret_)
        this%itretnz_(this%ntret_) = 0
      ENDIF
      ! local z-index to be returned to top (in _extended_ grid)
      this%itret_   (this%ntret_,k) = this%nzghost_+this%kend_-this%ksta_+1+kg
      ! Destination z-indices in _regular_ grid for top return.
      ! These indices should be in local--not global--form:
      this%itretdst_(this%ntret_,k) = k
      ! Set no. z-indices to return to top task
      this%itretnz_ (this%ntret_  ) = this%itretnz_(this%ntret_) + 1
      k = k + 1
    ENDDO

    ! *** Find bottom neighbors to return to:
    DO jf = 1, nt
      jfwd(jf) = modulo(this%myrank_-jf+this%nprocs_,this%nprocs_)
      CALL range(1,this%nd_(3),this%nprocs_,jfwd(jf),kfsta(jf),kfend(jf))
      nzf (jf) = kfend(jf) - kfsta(jf) + 1
    ENDDO

    this%nbret_ = 1
    this%ibretp_(this%nbret_) = jfwd(this%nbret_)
    this%ibretnz_(this%nbret_) = 0
    k = 1
    DO kg = 1, this%nzghost_
      IF ( k .GT. nzf(this%nbret_) ) THEN
        k = k - nzf(this%nbret_)
        this%nbret_ = this%nbret_ + 1
        this%ibretp_(this%nbret_)  = jfwd(this%nbret_)
        this%ibretnz_(this%nbret_) = 0
      ENDIF
      ! local z-index to be returned to top (in _extended_ grid)
      this%ibret_   (this%nbret_,k) = this%nzghost_ - kg + 1
      ! Destination z-indices in _regular_ grid for bottom return.
      ! These indices should be in local--not global--form:
      this%ibretdst_(this%nbret_,k) = nzf(this%nbret_) + 1 - k
      ! Set no. z-indices to return to bottom task
      this%ibretnz_ (this%nbret_  ) = this%ibretnz_(this%nbret_) + 1
      k = k + 1
    ENDDO
!    PRINT *, this%myrank_, this%nbret_, this%ibretp_(1:this%nbret_), &
!    this%ibretnz_(1:this%nbret_), this%ibret_(1:this%nbret_,1), this%ibretdst_(1:this%nbret_,1)
!    PRINT *, this%myrank_, this%ntsnd_, this%itsndp_(1:this%ntsnd_), &
!    this%itsndnz_(1:this%ntsnd_),this%itsnd_(1:this%ntsnd_,1), this%itsnddst_(1:this%ntsnd_,1)
!    PRINT *, this%myrank_, this%nbrcv_, this%ibrcvp_(1:this%nbrcv_), this%ibrcvnz_(1:this%nbrcv_)

    DEALLOCATE(jfwd,kfend,kfsta,nzf)

  END SUBROUTINE GFieldComm_ctor
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GFieldComm_dtor(this)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  Explicit destructor for field communicator. Should be called
!  after calling GPartComm_dtor.
!
!  ARGUMENTS:
!    this    : 'this' class instance
!-----------------------------------------------------------------
    IMPLICIT NONE
    TYPE(GFieldComm), INTENT(INOUT)           :: this

    IF ( ALLOCATED (this%ibretp_  ) ) DEALLOCATE (this%ibretp_  )
    IF ( ALLOCATED (this%itretp_  ) ) DEALLOCATE (this%itretp_  )
    IF ( ALLOCATED (this%ibretnz_ ) ) DEALLOCATE (this%ibretnz_ )
    IF ( ALLOCATED (this%itretnz_ ) ) DEALLOCATE (this%itretnz_ )
    IF ( ALLOCATED (this%ibret_   ) ) DEALLOCATE (this%ibret_   )
    IF ( ALLOCATED (this%itret_   ) ) DEALLOCATE (this%itret_   )
    IF ( ALLOCATED (this%ibretdst_) ) DEALLOCATE (this%ibretdst_)
    IF ( ALLOCATED (this%itretdst_) ) DEALLOCATE (this%itretdst_)

  END SUBROUTINE GFieldComm_dtor
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GFieldComm_SlabDataExchangeSF(this,v,vext,method)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : GFieldComm_SlabDataExchangeSF
!  DESCRIPTION: Does bdy exchange of field component, v. Output
!               is to data on regular grids, v. 'SF' means
!               that this is the 'single-field' interface.
!  ARGUMENTS  :
!    this      : 'this' class instance (IN)
!    v         : Eulerian velocity component returned on regular
!                grid.
!    vext      : Eulerian velocity components on extended grid.
!    method    : Ghost zone aggregation method
!                   = UNPACK_REP for replacement
!                   = UNPACK_SUM for sum
!
!-----------------------------------------------------------------
    IMPLICIT NONE
    CLASS(GFieldComm),INTENT(INOUT)                     :: this
    INTEGER      ,INTENT   (IN)                         :: method
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(*)            :: v
    REAL(KIND=GP),INTENT   (IN),DIMENSION(*)            :: vext

    INTEGER                                             :: itask,j

    IF ( this%intrfc_ .GE. 1 ) THEN
      WRITE(*,*) 'GFieldComm_SlabDataExchangeSF: SF interface expected'
      STOP
    ENDIF
    IF ( this%nprocs_ .EQ. 1 ) THEN
      CALL GFieldComm_LocalDataRetSF(this,vext,v,method)
      RETURN
    ENDIF

    CALL GFieldComm_Copy2Reg(this,v,vext)

    ! post receives:
    CALL GTStart(this%hcomm_)
    DO j=1,this%nbsnd_  ! from bottom task:
      itask = this%ibsndp_(j)
      CALL MPI_IRECV(this%rbbuff_(:,j),this%nbuff_,GC_REAL,itask, &
                     1,this%comm_,this%ibrh_(j),this%ierr_)
    ENDDO
    CALL GTAcc(this%hcomm_)

    ! return data:
    DO j=1,this%nbret_  ! to bottom task:
      itask = this%ibretp_(j)
!      PRINT *, 'Sending', j, this%myrank_, itask
      CALL GFieldComm_PackSF(this,this%sbbuff_(:,j),vext,j,'b')
      CALL GTStart(this%hcomm_)
      CALL MPI_ISEND(this%sbbuff_,this%nbuff_,GC_REAL,itask, &
                     1,this%comm_,this%ibsh_(j),this%ierr_)
      CALL GTAcc(this%hcomm_)
    ENDDO
!

    CALL GTStart(this%hcomm_)
    DO j=1,this%ntsnd_  ! from top task:
      itask = this%itsndp_(j)
      CALL MPI_IRECV(this%rtbuff_(:,j),this%nbuff_,GC_REAL,itask, &
                     1,this%comm_,this%itrh_(j),this%ierr_)
    ENDDO
    CALL GTAcc(this%hcomm_)


    DO j=1,this%ntret_  ! to top task:
      itask = this%itretp_(j)
      CALL GFieldComm_PackSF(this,this%stbuff_(:,j),vext,j,'t')
      CALL GTStart(this%hcomm_)
      CALL MPI_ISEND(this%stbuff_,this%nbuff_,GC_REAL,itask, &
                     1,this%comm_,this%itsh_(j),this%ierr_)
      CALL GTAcc(this%hcomm_)
    ENDDO

    CALL GTStart(this%hcomm_)
    DO j=1,this%nbret_
      CALL MPI_WAIT(this%ibsh_(j),this%istatus_,this%ierr_)
    ENDDO
    DO j=1,this%nbsnd_
      CALL MPI_WAIT(this%ibrh_(j),this%istatus_,this%ierr_)
    ENDDO


    DO j=1,this%ntret_
      CALL MPI_WAIT(this%itsh_(j),this%istatus_,this%ierr_)
    ENDDO
    DO j=1,this%ntsnd_
      CALL MPI_WAIT(this%itrh_(j),this%istatus_,this%ierr_)
    ENDDO
    CALL GTAcc(this%hcomm_)


    ! Unpack received data:
    DO j=1,this%nbsnd_
       CALL GFieldComm_UnpackSF(this,v,this%rbbuff_(:,j),&
                          this%nbuff_,'b',method,this%ierr_)
    ENDDO
    DO j=1,this%ntsnd_
       CALL GFieldComm_UnpackSF(this,v,this%rtbuff_(:,j),&
                          this%nbuff_,'t',method,this%ierr_)
    ENDDO

  END SUBROUTINE GFieldComm_SlabDataExchangeSF
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GFieldComm_PackSF(this,buff,vext,isnd,sdir)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : PPackSF
!  DESCRIPTION: packs ret buffer with (single) field
!  ARGUMENTS  :
!    this    : 'this' class instance (IN)
!    buff    : packed buffer (returned)
!    vext    : Eulerian velocity component on extended grid
!              in phys. space (IN)
!    isnd    : which send this is
!    sdir    : 't' for top, 'b' for bottom
!
!-----------------------------------------------------------------
    IMPLICIT NONE

    CLASS(GFieldComm),INTENT(INOUT)         :: this
    INTEGER      ,INTENT   (IN)             :: isnd
    INTEGER                                 :: i,j,k,m,nt,nx,ny,nxy
    INTEGER                                 :: jm,km
    REAL(KIND=GP),INTENT  (OUT),DIMENSION(*):: buff
    REAL(KIND=GP),INTENT   (IN),DIMENSION(*):: vext
    CHARACTER*(*),INTENT   (IN)             :: sdir


    IF ( sdir(1:1).NE.'b' .AND. sdir(1:1).NE.'B' &
    .AND.sdir(1:1).NE.'t' .AND. sdir(1:1).NE.'T' ) THEN
      WRITE(*,*) 'GFieldComm_PackSF: Bad direction descriptor'
      STOP
    ENDIF

    nx  = this%nd_(1)
    ny  = this%nd_(2)
    nxy = nx*ny
    IF      ( sdir(1:1) .EQ. 'b' .OR. sdir(1:1) .EQ. 'B' ) THEN
    ! Pack for send to rank at bottom:
    !  ...header
      nt = 1
      buff(1)  = this%ibretnz_(isnd)      ! no. z-indices included
      DO j = 1, this%ibretnz_(isnd)
        nt       = nt + 1
        buff(nt) = this%ibretdst_(isnd,j) ! z-index in regular grid
      ENDDO

    !  ...data
      DO m = 1,this%ibretnz_(isnd)
        k = this%ibret_(isnd,m)
        km = k-1
        DO j = 1, ny
          jm = j-1
          DO i = 1, nx
            nt = nt + 1
            buff(nt) = vext(i+jm*nx+km*nxy)
          ENDDO
        ENDDO
      ENDDO

    ELSE !  Pack for send to rank at top:

      ! ...header
      nt = 1
      buff(1)  = this%itretnz_(isnd)      ! no. z-indices included
      DO j = 1, this%itretnz_(isnd)
        nt       = nt + 1
        buff(nt) = this%itretdst_(isnd,j) ! z-index in regular grid
      ENDDO

      ! ...data
      DO m = 1,this%itretnz_(isnd)
        k = this%itret_(isnd,m)
        km = k-1
        DO j = 1, ny
          jm = j-1
          DO i = 1, nx
            nt = nt + 1
            buff(nt) = vext(i+jm*nx+km*nxy)
          ENDDO
        ENDDO
      ENDDO

    ENDIF

  END SUBROUTINE GFieldComm_PackSF
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GFieldComm_UnpackSF(this,v,buff,nbuff,sb,method,ierr)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : RetUnpackSF
!  DESCRIPTION: Unpacks recv buffer with into regular (single) field.
!               Messages received are self-referential, so contain info
!               on where to 'return' recvd data. So, there is no 't'
!               or 'b' designation required for unpacking.
!  ARGUMENTS  :
!    this        : 'this' class instance (IN)
!    v           : Eulerian velocity component on regular grid
!                  in phys. space (IN)
!    buff        : packed buffer (input) from which to store into
!                  extended grid quantities.
!    nbuff       : maximum buff size
!    sb          : optional buffer name ,'t' or 'b'.
!    method      : data aggregation method
!                   = UNPACK_REP for replacement
!                   = UNPACK_SUM for sum
!    ierr        : err flag: 0 if success; else 1
!
!-----------------------------------------------------------------
    IMPLICIT NONE

    CLASS(GFieldComm),INTENT(INOUT)         :: this
    INTEGER                                 :: i,j,k,m,ngp,nx,nex,nxy,nexy,ny,nz
    INTEGER                                 :: im,ip,ir,ixy,iz,jm,km,nh
    INTEGER      ,INTENT   (IN)             :: nbuff,method
    INTEGER      ,INTENT(INOUT)             :: ierr ! not used now
    REAL(KIND=GP),INTENT   (IN),DIMENSION(*):: buff
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(*):: v
    CHARACTER(len=1),INTENT(IN),OPTIONAL    :: sb

    ierr = 0
    nx   = this%nd_(1)
    ny   = this%nd_(2)
    nxy  = nx*ny
    ngp  = this%nzghost_*this%iextperp_
    nex  = nx+2*ngp
    nexy = (nx+2*ngp)*(ny+2*ngp)

    ! Unpack from either buffer:
    ! For each task, message is of form:
    !     #z-indices:z-index_0:z-index_1:...:nx*ny_0:nx*ny_1: ...
    nz = int(buff(1))
    nh = nz + 1 ! no. items in header
    IF ( method .EQ. UNPACK_REP ) THEN
      DO m = 1, nz
        k   = int(buff(m+1))
        km  = k-1
        ixy = 1
        DO j = 1, ny
          jm = j-1
          DO i = 1, nx
            im = i+ngp+(jm+ngp)*nex+km*nexy
            ir = nh+(m-1)*nxy+ixy
            v(im) = buff(ir)
            ixy = ixy + 1
          ENDDO
        ENDDO
      ENDDO
    ELSE IF ( method .EQ. UNPACK_SUM ) THEN
      DO m = 1, nz
        k   = int(buff(m+1))
        km  = k-1
        ixy = 1
        DO j = 1, ny
          jm = j-1
          DO i = 1, nx
            im = i+ngp+(jm+ngp)*nex+km*nexy
            ir = nh+(m-1)*nxy+ixy
            v(im) = v(im) + buff(ir)
            ixy = ixy + 1
          ENDDO
        ENDDO
      ENDDO
    ENDIF

  END SUBROUTINE GFieldComm_UnpackSF
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GFieldComm_LocalDataRetSF(this,vext,v,method)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : LocalDataExchSF
!  DESCRIPTION: Does 'bdy exchange' of (single) velocity component, when there's
!               only a single MPI task. This is a single-field interface.
!  ARGUMENTS  :
!    this              : 'this' class instance (IN)
!    v                 : Eulerian velocity components returned on regular grid.
!                        Must be of size nd_ set in constructor.
!    vext              : Eulerian velocity components on extended grid (that
!                        used to hold ghost zones).
!    method            : data aggregation method
!                         = UNPACK_REP for replacement
!                         = UNPACK_SUM for sum
!
!-----------------------------------------------------------------
    USE mpivars

    IMPLICIT NONE
    CLASS(GFieldComm),INTENT(INOUT)                     :: this
    INTEGER                                             :: i,j,k,ngp,ngz,nex,nexy,nez
    INTEGER                                             :: nx,nxy,ny,nz
    INTEGER      ,INTENT   (IN)                         :: method
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(*)            :: v
    REAL(KIND=GP),INTENT   (IN),DIMENSION(*)            :: vext

    ngz  = this%nzghost_
    ngp  = ngz * this%iextperp_
    nexy = (this%nd_(1)+2*ngp) * (this%nd_(2)+2*ngp)
    nex  = this%nd_(1)+2*ngp
    nez  = this%nd_(3)+  ngz
    nx   = this%nd_(1)
    ny   = this%nd_(2)
    nz   = this%nd_(3)
    nxy  = nx*ny

    CALL GFieldComm_Copy2Reg(this,v,vext)
    IF ( method .EQ. UNPACK_REP ) THEN
      DO k = 1, ngz  ! extended zones
        DO j=1,ny
          DO i=1,nx
            ! set top bcs:
            v(i+(j-1)*nx+       (k-1)*nxy) = vext(i+ngp+(j+ngp-1)*nex+(nez+k-1)*nexy)

            ! set bottom bcs:
            v(i+(j-1)*nx+(nz-ngz+k-1)*nxy) = vext(i+ngp+(j+ngp-1)*nex+    (k-1)*nexy)

          ENDDO
        ENDDO
      ENDDO
    ELSE IF ( method .EQ. UNPACK_SUM ) THEN
      DO k = 1, ngz  ! extended zones
        DO j=1,ny
          DO i=1,nx
            ! set top bcs:
            v(i+(j-1)*nx+(k-1)*nxy)        = &
                   v(i+(j-1)*nx+(k-1)*nxy) + vext(i+ngp+(j+ngp-1)*nex+(nez+k-1)*nexy)

            ! set bottom bcs:
            v(i+(j-1)*nx+(nz-ngz+k-1)*nxy) = &
                v(i+(j-1)*nx+(nz-ngz+k-1)*nxy) + vext(i+ngp+(j+ngp-1)*nex+(k-1)*nexy)

          ENDDO
        ENDDO
      ENDDO
    ENDIF

  END SUBROUTINE GFieldComm_LocalDataRetSF
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GFieldComm_SlabDataExchangeMF(this,vx,vy,vz,vxext,vyext,vzext,method)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : GFieldComm_SlabDataExchangeMF
!  DESCRIPTION: Does bdy exchange of velocity component, vx,vy,vz. Output
!               is to data on extended grids, vxext, vyext, vzexy. 'MF'
!               means that this is the 'multi-field' interface.
!  ARGUMENTS  :
!    this              : 'this' class instance (IN)
!    vx,vy,vz          : Eulerian velocity components returned on regular
!                        grid
!    vxext,vyext,vzext : Eulerian velocity components on extended grid.
!    method            : Ghost zone aggregation method
!                         = UNPACK_REP for replacement
!                         = UNPACK_SUM for sum
!
!-----------------------------------------------------------------
    IMPLICIT NONE
    CLASS(GFieldComm),INTENT(INOUT)                     :: this
    REAL(KIND=GP),INTENT   (IN),DIMENSION(*)            :: vxext,vyext,vzext
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(*)            :: vx,vy,vz
    INTEGER      ,INTENT   (IN)                         :: method
    INTEGER                                             :: itask,j

    IF ( this%intrfc_ .LT. 1 ) THEN
      WRITE(*,*) 'GFieldComm_SlabDataExchangeMF: SF interface expected'
      STOP
    ENDIF

    IF ( this%nprocs_ .EQ. 1 ) THEN
      CALL GFieldComm_LocalDataRetMF(this,vxext,vyext,vzext,vx,vy,vz,method)
      RETURN
    ENDIF

    ! Post receives:
    CALL GTStart(this%hcomm_)
    DO j=1,this%nbsnd_  ! from bottom task:
      itask = this%ibsndp_(j)
      CALL MPI_IRECV(this%rbbuff_(:,j),this%nbuff_,GC_REAL,itask, &
                     1,this%comm_,this%ibrh_(j),this%ierr_)
    ENDDO
    CALL GTAcc(this%hcomm_)

    CALL GTStart(this%hcomm_)
    DO j=1,this%ntsnd_  ! from top task:
      itask = this%itsndp_(j)
      CALL MPI_IRECV(this%rtbuff_(:,j),this%nbuff_,GC_REAL,itask, &
                     1,this%comm_,this%itrh_(j),this%ierr_)
    ENDDO
    CALL GTAcc(this%hcomm_)

    !
    ! return data:
    DO j=1,this%nbret_  ! to bottom task:
      itask = this%ibretp_(j)
      CALL GFieldComm_PackMF(this,this%sbbuff_(:,j),vxext,vyext,vzext,j,'b')
      CALL GTStart(this%hcomm_)
      CALL MPI_ISEND(this%sbbuff_,this%nbuff_,GC_REAL,itask, &
                     1,this%comm_,this%ibsh_(j),this%ierr_)
      CALL GTAcc(this%hcomm_)
    ENDDO
    DO j=1,this%ntret_  ! to top task:
      itask = this%itretp_(j)
      CALL GFieldComm_PackMF(this,this%stbuff_(:,j),vxext,vyext,vzext,j,'t')
      CALL GTStart(this%hcomm_)
      CALL MPI_ISEND(this%stbuff_,this%nbuff_,GC_REAL,itask, &
                     1,this%comm_,this%itsh_(j),this%ierr_)
      CALL GTAcc(this%hcomm_)

    ENDDO

    CALL GTStart(this%hcomm_)
    DO j=1,this%nbret_
      CALL MPI_WAIT(this%ibsh_(j),this%istatus_,this%ierr_)
    ENDDO
    DO j=1,this%ntret_
      CALL MPI_WAIT(this%itsh_(j),this%istatus_,this%ierr_)
    ENDDO
    DO j=1,this%nbsnd_
      CALL MPI_WAIT(this%ibrh_(j),this%istatus_,this%ierr_)
    ENDDO
    DO j=1,this%ntsnd_
      CALL MPI_WAIT(this%itrh_(j),this%istatus_,this%ierr_)
    ENDDO
    CALL GTAcc(this%hcomm_)

    ! Unpack received data:
    DO j=1,this%nbsnd_
      CALL GFieldComm_UnpackMF(this,vx,vy,vz,this%rbbuff_(:,j),method)
    ENDDO
    DO j=1,this%ntsnd_
      CALL GFieldComm_UnpackMF(this,vx,vy,vz,this%rtbuff_(:,j),method)
    ENDDO

  END SUBROUTINE GFieldComm_SlabDataExchangeMF
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GFieldComm_PackMF(this,buff,vxext,vyext,vzext,isnd,sdir)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : PackMF
!  DESCRIPTION: packs ret buffer with fields; multi-field interface
!  ARGUMENTS  :
!    this             : 'this' class instance (IN)
!    buff             : packed buffer (returned)
!    vxext,vyext,vzext: Eulerian velocity component on extended
!                       grid in phys. space (IN)
!    isnd             : which send this is
!    sdir             : 't' for top, 'b' for bottom
!
!-----------------------------------------------------------------
    IMPLICIT NONE

    CLASS(GFieldComm),INTENT(INOUT)         :: this
    INTEGER      ,INTENT   (IN)             :: isnd
    INTEGER                                 :: i,j,k,m,nt,nx,nxy,ny
    INTEGER                                 :: jm,km
    REAL(KIND=GP),INTENT  (OUT),DIMENSION(*):: buff
    REAL(KIND=GP),INTENT   (IN),DIMENSION(*):: vxext,vyext,vzext
    CHARACTER*(*),INTENT   (IN)             :: sdir


    IF ( sdir(1:1).NE.'b' .AND. sdir(1:1).NE.'B' &
    .AND.sdir(1:1).NE.'t' .AND. sdir(1:1).NE.'T' ) THEN
      WRITE(*,*) 'GFieldComm_PackMF: Bad direction descriptor'
      STOP
    ENDIF

    nx  = this%nd_(1)
    ny  = this%nd_(2)
    nxy = nx*ny
    IF      ( sdir(1:1) .EQ. 'b' .OR. sdir(1:1) .EQ. 'B' ) THEN
    ! Pack for send to rank at bottom:
    !  ...header
      nt = 1
      buff(1)  = this%ibretnz_(isnd)       ! no. z-indices included
      DO j = 1, this%ibretnz_(isnd)
        nt       = nt + 1
        buff(nt) = this%ibretdst_(isnd,j) ! z-index in regular grid
      ENDDO

    !  ...data
      DO m = 1,this%ibretnz_(isnd)
        k = this%ibret_(isnd,m)
        km = k-1
        DO j = 1, ny
          jm = j-1
          DO i = 1, nx
            buff(nt) = vxext(i+jm*nx+km*nxy)
            nt = nt + 1
          ENDDO
        ENDDO
      ENDDO

      DO m = 1, this%ibretnz_(isnd)
        k = this%ibret_(isnd,m)
        km = k-1
        DO j = 1, ny
          jm = j-1
          DO i = 1, nx
            buff(nt) = vyext(i+jm*nx+km*nxy)
            nt = nt + 1
          ENDDO
        ENDDO
      ENDDO

      DO m = 1,this%ibretnz_(isnd)
        k = this%ibret_(isnd,m)
        km = k-1
        DO j = 1, ny
          jm = j-1
          DO i = 1, nx
            buff(nt) = vzext(i+jm*nx+km*nxy)
            nt = nt + 1
          ENDDO
        ENDDO
      ENDDO

    ELSE !  Pack for send to rank at top:

      ! ...header
      nt = 1
      buff(1)  = this%itretnz_(isnd)      ! no. z-indices included
      DO j = 1, this%itretnz_(isnd)
        nt       = nt + 1
        buff(nt) = this%itretdst_(isnd,j) ! z-index in regular grid
      ENDDO

      ! ...data
      DO m = 1,this%itretnz_(isnd)
        k = this%itret_(isnd,m)
        km = k-1
        DO j = 1, ny
          jm = j-1
          DO i = 1, nx
            buff(nt) = vxext(i+jm*nx+km*nxy)
            nt = nt + 1
          ENDDO
        ENDDO
      ENDDO

      DO m = 1,this%itretnz_(isnd)
        k = this%itret_(isnd,m)
        km = k-1
        DO j = 1, ny
          jm = j-1
          DO i = 1, nx
            buff(nt) = vyext(i+jm*nx+km*nxy)
            nt = nt + 1
          ENDDO
        ENDDO
      ENDDO

      DO m = 1,this%itretnz_(isnd)
        k = this%itret_(isnd,m)
        km = k-1
        DO j = 1, ny
          jm = j-1
          DO i = 1, nx
            buff(nt) = vzext(i+jm*nx+km*nxy)
            nt = nt + 1
          ENDDO
        ENDDO
      ENDDO

    ENDIF

  END SUBROUTINE GFieldComm_PackMF
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GFieldComm_UnpackMF(this,vx,vy,vz,buff,method)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : UnpackMF
!  DESCRIPTION: Unpacks recv buffer with into extended (single) field
!               Messages received are self-referential, so contain info
!               on where to 'send' recvd data. So, there is no 't' or
!               'b' designation required for unpacking.
!  ARGUMENTS  :
!    this        : 'this' class instance (IN)
!    buff        : packed buffer (input) from which to store into
!                  extended grid quantities.
!    vx,vy,vz    : Eulerian velocity component on regular grid
!                  in phys. space (IN)
!    method      : data aggregation method
!                   = UNPACK_REP for replacement
!                   = UNPACK_SUM for sum
!
!-----------------------------------------------------------------
    USE mpivars
    IMPLICIT NONE

    CLASS(GFieldComm),INTENT(INOUT)         :: this
    INTEGER                                 :: i,j,k,m,ngp,ngz,nx,nxy,ny,nz
    INTEGER                                 :: ixy,jm,km,nh
    INTEGER      ,INTENT   (IN)             :: method
    REAL(KIND=GP),INTENT   (IN),DIMENSION(*):: buff
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(*):: vx,vy,vz

    nx  = this%nd_(1)
    ny  = this%nd_(2)
    nxy = nx*ny
    ngz = this%nzghost_;
    ngp = ngz*this%iextperp_

  ! Unpack from either top or bottom buffer:
    nz = int(buff(1))
    nh = nz + 1 ! no. items in header
    IF ( method .EQ. UNPACK_REP ) THEN
      DO m = 1,nz
        k   = int(buff(m+1))
        km  = k-1

        ixy = 1
        DO j = 1, ny
          jm = j-1
          DO i = 1, nx
            vx(i+ngp+(jm+ngp)*nx+km*nxy) = buff(nh+(m-1)*nxy+ixy)
            ixy = ixy + 1
          ENDDO
        ENDDO

        DO j = 1, ny
          jm = j-1
          DO i = 1, nx
            vy(i+ngp+(jm+ngp)*nx+km*nxy) = buff(nh+(m-1)*nxy+ixy)
            ixy = ixy + 1
          ENDDO
        ENDDO

        DO j = 1, ny
          jm = j-1
          DO i = 1, nx
            vz(i+ngp+(jm+ngp)*nx+km*nxy) = buff(nh+(m-1)*nxy+ixy)
            ixy = ixy + 1
          ENDDO
        ENDDO

      ENDDO
    ELSE IF ( method .EQ. UNPACK_SUM ) THEN
      DO m = 1,nz
        k   = int(buff(m+1))
        km  = k-1

        ixy = 1
        DO j = 1, ny
          jm = j-1
          DO i = 1, nx
            vx(i+ngp+(jm+ngp)*nx+km*nxy) = vx(i+ngp+(jm+ngp)*nx+km*nxy) &
                                          + buff(nh+(m-1)*nxy+ixy)
            ixy = ixy + 1
          ENDDO
        ENDDO

        DO j = 1, ny
          jm = j-1
          DO i = 1, nx
            vy(i+ngp+(jm+ngp)*nx+km*nxy) = vy(i+ngp+(jm+ngp)*nx+km*nxy) &
                                          + buff(nh+(m-1)*nxy+ixy)
            ixy = ixy + 1
          ENDDO
        ENDDO

        DO j = 1, ny
          jm = j-1
          DO i = 1, nx
            vz(i+ngp+(jm+ngp)*nx+km*nxy) = vz(i+ngp+(jm+ngp)*nx+km*nxy) &
                                          + buff(nh+(m-1)*nxy+ixy)
            ixy = ixy + 1
          ENDDO
        ENDDO

      ENDDO
    ENDIF

  END SUBROUTINE GFieldComm_UnpackMF
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GFieldComm_LocalDataRetMF(this,vxext,vyext,vzext,vx,vy,vz,method)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : LocalDataExch
!  DESCRIPTION: Does 'bdy exchange' of velocity component, when there's
!               only a single MPI task.
!  ARGUMENTS  :
!    this              : 'this' class instance (IN)
!    vxext,vyext,vzext : Eulerian velocity components returned on extended
!                        grid (that used to hold ghost zones). Only z-conditions
!                        are imposed; lateral periodicity is not handled here.
!                        Lateral ghost zones can be accounted for by setting
!                        this%iextperp_=1 in contructor.
!    vx,vy,vz          : Eulerian velocity components returned on regular
!                        grid. Must be of size nd_ set in constructor
!    method            : data aggregation method
!                         = UNPACK_REP for replacement
!                         = UNPACK_SUM for sum
!
!-----------------------------------------------------------------
    USE mpivars

    IMPLICIT NONE
    CLASS(GFieldComm),INTENT(INOUT)                     :: this
    INTEGER                                             :: i,j,k,ngp,ngz,nex,nexy,nez
    INTEGER                                             :: nx,nxy,ny,nz
    INTEGER                                             :: jm,km
    INTEGER      , INTENT  (IN)                         :: method
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(*)            :: vx,vy,vz
    REAL(KIND=GP),INTENT   (IN),DIMENSION(*)            :: vxext,vyext,vzext

    ngz  = this%nzghost_
    ngp  = ngz * this%iextperp_
    nexy = (this%nd_(1)+2*ngp) * (this%nd_(2)+2*ngp)
    nex  = this%nd_(1)+2*ngp
    nez  = this%nd_(3)+  ngz
    nx   = this%nd_(1)
    ny   = this%nd_(2)
    nz   = this%nd_(3)
    nxy  = nx*ny

    CALL GFieldComm_Copy2Reg(this,vx,vxext)
    CALL GFieldComm_Copy2Reg(this,vy,vyext)
    CALL GFieldComm_Copy2Reg(this,vz,vzext)
    IF ( method .EQ. UNPACK_REP ) THEN
      DO k = 1, ngz  ! bottom extended zones
        km = k-1
        DO j=1,ny
          jm = j-1
          DO i=1,nx
            ! set top bcs:
            vx(i+(j-1)*nx+(k-1)*nxy) = vxext(i+ngp+(j+ngp-1)*nex+(nez+k-1)*nexy)
            vy(i+(j-1)*nx+(k-1)*nxy) = vyext(i+ngp+(j+ngp-1)*nex+(nez+k-1)*nexy)
            vz(i+(j-1)*nx+(k-1)*nxy) = vzext(i+ngp+(j+ngp-1)*nex+(nez+k-1)*nexy)

            ! set bottom bcs:
            vx(i+(j-1)*nx+(nz-ngz+k-1)*nxy) = vxext(i+ngp+(j+ngp-1)*nex+    (k-1)*nexy)
            vy(i+(j-1)*nx+(nz-ngz+k-1)*nxy) = vyext(i+ngp+(j+ngp-1)*nex+    (k-1)*nexy)
            vz(i+(j-1)*nx+(nz-ngz+k-1)*nxy) = vzext(i+ngp+(j+ngp-1)*nex+    (k-1)*nexy)

          ENDDO
        ENDDO
      ENDDO
    ELSE IF ( method .EQ. UNPACK_SUM ) THEN
      DO k = 1, ngz  ! bottom extended zones
        km = k-1
        DO j=1,ny
          jm = j-1
          DO i=1,nx
            ! set top bcs:
            vx(i+(j-1)*nx+(k-1)*nxy) = vx(i+(j-1)*nx+(k-1)*nxy) &
                              + vxext(i+ngp+(j+ngp-1)*nex+(nez+k-1)*nexy)
            vy(i+(j-1)*nx+(k-1)*nxy) = vy(i+(j-1)*nx+(k-1)*nxy) &
                              + vyext(i+ngp+(j+ngp-1)*nex+(nez+k-1)*nexy)
            vz(i+(j-1)*nx+(k-1)*nxy) = vz(i+(j-1)*nx+(k-1)*nxy) &
                              + vzext(i+ngp+(j+ngp-1)*nex+(nez+k-1)*nexy)

            ! set bottom bcs:
            vx(i+(j-1)*nx+(nz-ngz+k-1)*nxy) = vx(i+(j-1)*nx+(nz-ngz+k-1)*nxy) &
                              + vxext(i+ngp+(j+ngp-1)*nex+    (k-1)*nexy)
            vy(i+(j-1)*nx+(nz-ngz+k-1)*nxy) = vy(i+(j-1)*nx+(nz-ngz+k-1)*nxy) &
                              + vyext(i+ngp+(j+ngp-1)*nex+    (k-1)*nexy)
            vz(i+(j-1)*nx+(nz-ngz+k-1)*nxy) = vz(i+(j-1)*nx+(nz-ngz+k-1)*nxy) &
                              + vzext(i+ngp+(j+ngp-1)*nex+    (k-1)*nexy)

          ENDDO
        ENDDO
      ENDDO
    ENDIF

  END SUBROUTINE GFieldComm_LocalDataRetMF
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GFieldComm_Copy2Reg(this,v,vext)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : Copy2Reg
!  DESCRIPTION: Copy field from extended to regular grid,
!               without ghost zones.
!
!  ARGUMENTS  :
!    this     : 'this' class instance (IN)
!    v        : regular-grid field
!    vext     : extended-grid field
!    ldims    : local dims of v
!-----------------------------------------------------------------
    IMPLICIT NONE
    CLASS(GFieldComm),INTENT(INOUT)                     :: this
    INTEGER                                             :: i,j,jm,k,km,ngp,ngz,nex,nexy
    INTEGER                                             :: nx,nxy,ny
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(*)            :: v
    REAL(KIND=GP),INTENT   (IN),DIMENSION(*)            :: vext

    ngz  = this%nzghost_
    ngp  = ngz * this%iextperp_
    nexy = (this%nd_(1)+2*ngp) * (this%nd_(2)+2*ngp)
    nex  = this%nd_(1)+2*ngp
    nx   = this%nd_(1)
    ny   = this%nd_(2)
    nxy  = nx*ny

    DO k = 1,this%kend_-this%ksta_+1
      km = k-1
      DO j=1,ny
        jm = j-1
        DO i=1,nx
          v(i+jm*nx+km*nxy) = vext(i+ngp+(jm+ngp)*nex+(km+ngz)*nexy)
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE GFieldComm_Copy2Reg
!-----------------------------------------------------------------
!-----------------------------------------------------------------
