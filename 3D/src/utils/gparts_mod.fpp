!=================================================================
! GHOST GPart Lagrangian particles class
!
! 2013 D. Rosenberg
!      NCCS: ORNL
!
! 15 Jan 2013: Initial version
!=================================================================
MODULE class_GPart
      USE mpivars
      USE fprecision
      USE pdbtypes
      USE GPartComm
      USE class_GSplineInt

      IMPLICIT NONE
      
      ENUM, BIND(C) :: GPINIT
        ENUMERATOR :: GPINIT_RANDLOC=0
        ENUMERATOR :: GPINIT_USERLOC
      ENDENUM GPINIT

      ENUM, BIND(C) :: GPINTRP
        ENUMERATOR :: GPINTRP_CSPLINE=0
        ENUMERATOR :: GPINTRP_LAGINT 
        ENUMERATOR :: GPINTRP_SPECTRAL
      ENDENUM GPINTRP

      TYPE, PUBLIC :: GPart
        PRIVATE
        ! Member data:
        INTEGER, DIMENSION(MPI_STATUS_SIZE)          :: istatus_
        TYPE(GPINIT)                                 :: inittype_
        TYPE(GPINTRP)                                :: iinterp_
        TYPE(GPartComm)                              :: exchop_
        TYPE(GPSplineInt)                            :: intop_
        INTEGER                                      :: intorder_,itorder_,nd_(3),Next_(3),libnds_(3,2)
        INTEGER                                      :: myrank,nprocs_
        INTEGER                                      :: ierr_,nwrk_,seed_
        INTEGER                                      :: maxparts_,nparts_
        INTEGER                                      :: iotype_,sziotype_
!!      TYPE(GPDBrec), ALLOCATABLE, DIMENSION    (:) :: pdb_
        REAL(KIND=GP), ALLOCATABLE, DIMENSION    (:) :: px_,py_,pz_
        REAL(KIND=GP), ALLOCATABLE, DIMENSION  (:,:) :: ptmp0_,ptmp1_
        REAL(KIND=GP), ALLOCATABLE, DIMENSION    (:) :: ltmp0_,ltmp1_
        REAL(KIND=GP), ALLOCATABLE, DIMENSION    (:) :: X_,Y_,Z_,lx_,ly_,lz_
        REAL(KIND=GP)                                :: lxbnds_(3,2),gxbnds_(3,2)
        CHARACTER(len=1024)                          :: seedfile_,serr_,sfile_
      CONTAINS
        ! Public methods:
        PROCEDURE,PUBLIC :: GPart_ctor
        FINAL            :: GPart_dtor
        PROCEDURE,PUBLIC :: Init        => GPart_init
        PROCEDURE,PUBLIC :: Step        => GPart_StepRKK
        PROCEDURE,PUBLIC :: Setstep     => GPart_SetStepRKK
        PROCEDURE,PUBLIC :: Finalize    => GPart_FinalStepRKK
        PROCEDURE,PUBLIC :: io_write    => GPart_io_write_euler,GPart_io_write_pdb
        PROCEDURE,PUBLIC :: io_read     => GPart_io_read
!       PROCEDURE,PUBLIC :: GetPos
        PROCEDURE,PUBLIC :: EulerToLag  => GPart_EulerToLag,GPart_EulerToLocalLag
        PROCEDURE,PUBLIC :: SetInitType => GPart_SetInitType
        PROCEDURE,PUBLIC :: SetSeedFile => GPart_SetSeedFile
        PROCEDURE,PUBLIC :: SetRandSeed => GPart_SetRandSeed
        PROCEDURE,PUBLIC :: SetTimeOrder=> GPart_SetTimeOrder
        PROCEDURE,PUBLIC :: GetSeedFile => GPart_GetSeedFile
        PROCEDURE,PUBLIC :: GetRandSeed => GPart_GetRandSeed
        PROCEDURE,PUBLIC :: GetTimeOrder=> GPart_GetTimeOrder
      END TYPE GPart

      PRIVATE :: GPart_init        , GPart_StepRKK     , GPart_io_write_euler
      PRIVATE :: GPart_SetStepRKK  , GPart_FinalStepRKK
      PRIVATE :: GPart_iowrite_pdb , GPart_io_read     , 
      PRIVATE :: GPart_InitRandSeed, GPart_InitUserSeed, 
      PRIVATE :: GPart_SetInitType , GPart_SetSeedFile
      PRIVATE :: GPart_SetRandSeed , GPart_Delete
      PRIVATE :: GPart_GetWrkInd   , GPart_GetParts    , GPart_MakePeriodicP
      PRIVATE :: GPart_Copy2Ext    , GPart_MakePeriodicExt

! Methods:
  CONTAINS

  SUBROUTINE GPart_ctor(this, comm, mparts, inittype, iinterp, intorder)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  Main explicit constructor
!  ARGUMENTS:
!    this    : 'this' class instance
!    comm    : MPI communicator
!    mparts  : max no. particles allowed on grid. If inittype = GPINIT_USERLOC
!              then this is computed from the user file, and can
!              be any value.
!    inittype: GPINIT-typed quantity to give the type of particle
!              initialization
!    iinterp : GPINTRP-type quantity providing the interpolation
!              scheme
!    intorder: for variable-order (e.g., Lagrange) interpolation, 
!              the order (1, 2, 3...). Sets the number of 'ghost' zones
!              of data transferred between MPI tasks.
!-----------------------------------------------------------------
    USE grid
    USE mpivars

    IMPLICIT NONE
    CLASS(GPart)                      :: this
    INTEGER       , INTENT(IN)        :: mparts,comm
    INTEGER                           :: disp(3),lens(3),types(3),szreal
    TYPE(GPINIT)  , INTENT(IN)        :: inittype
    TYPE(GPINTRP) , INTENT(IN)        :: iinterp

    INTEGER                           :: j,nc

    this%nparts_   = 0 
    this%comm_     = comm
    this%maxparts_ = mparts
    this%nd_(1:3)  = n
    this%nwrk_     = 0
    this%seedfile_ = 'gploc.dat'
    this%iinterp_  = 3          ! fixed for now
    this%inittype_ = inittype
    this%itorder_  = 2
    this%intorder_ = max(intorder,1)
    this%seed_     = 1000

    IF ( this%intorder_ .NE. 3 ) THEN
      WRITE(*,*) 'GPart::ctor: Only 3rd order allowed for now' 
    ENDIF

    CALL MPI_COMM_SIZE(this%comm_,this%nprocs_,this%ierr_)
    CALL MPI_COMM_RANK(this%comm_,this%myrank_,this%ierr_)

    CALL this%exchop_%GPartComm_ctor(1,this%maxparts_,this%nd_,this%intorder_-1,this%comm_)
    CALL this%exchop_%Init()

    this%libnds_(1,1) = 1    ; this%lxbnds_(1,1) = 0.0
    this%libnds_(1,2) = n    ; this%lxbnds_(1,2) = real(n-1,kind=GP)
    this%libnds_(2,1) = 1    ; this%lxbnds_(2,1) = 0.0
    this%libnds_(2,2) = n    ; this%lxbnds_(2,2) = real(n-1,kind=GP)
    this%libnds_(3,1) = kstaa; this%lxbnds_(3,1) = real(ksta-1,kind=GP)
    this%libnds_(3,2) = kend ; this%lxbnds_(3,2) = real(kend-1,kind=GP)

    DO j = 1,3
      this%gxbnds_(j,1) = 0.0_GP 
      this%gxbnds_(j,2) = real(this%nd_(j)-1,kind=GP)
    ENDDO
    CALL intop_%GPSplineInt_ctor(3,this%nd_,ibnds,xbnds,this%maxparts_,this%exchop_)

    ! Create part. d.b. structure type for I/O
    CALL MPI_TYPE_SIZE(GC_REAL,szreal)
    disp (1)=0 
    disp (2)=szreal
    disp (3)=2*szreal
    types(1)=MPI_FLOAT 
    types(2)=MPI_FLOAT 
    types(3)=MPI_FLOAT
    CALL_TYPE_CREATE_STRUCT(3,lens,disp,this%iotype_,this%ierr_)
    CALL_TYPE_COMMIT(this%iotype_,this%ierr_)
    CALL MPI_TYPE_SIZE(this%iotype_,this%sziotype_,this%ierr_)

    DO j=1,2
      this%Next_(j) = this%nd_(j) + 2*(this%intorder_-1)
    ENDDO
    this%Next_(3) = kend-ksta + 1 + 2*(this%intorder_-1)

    ALLOCATE(this%px_      (this%maxparts_))
    ALLOCATE(this%py_      (this%maxparts_))
    ALLOCATE(this%pz_      (this%maxparts_))
    ALLOCATE(this%ptmp0_ (3,this%maxparts_))
    ALLOCATE(this%ptmp1_ (3,this%maxparts_))
    ALLOCATE(this%ltmp1_ (this%maxparts_))
    ALLOCATE(this%ltmp2_ (this%maxparts_))
    ALLOCATE(this%X_ (nd_(1)))
    ALLOCATE(this%Y_ (nd_(2)))
    ALLOCATE(this%Z_ (nd_(3)))
    ALLOCATE(this%lx_(this%iorder_+1))
    ALLOCATE(this%ly_(this%iorder_+1)
    ALLOCATE(this%lz_(this%iorder_+1)

    DO j = 1,nd_(1)
      this%X_(j) = REAL(j-1,KIND=GP)
    ENDDO
    DO j = 1,nd_(2)
      this%Y_(j) = REAL(j-1,KIND=GP)
    ENDDO
    DO j = ksta,ksta+nd_(3)-1
      this%Z_(j) = REAL(j-1,KIND=GP)
    ENDDO

  END SUBROUTINE GPart_ctor
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GPart_dtor(this)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  Main explicit destructor
!  ARGUMENTS:
!    this    : 'this' class instance
!-----------------------------------------------------------------

    IMPLICIT NONE
    TYPE(GPart)                      :: this


    CALL this%exchop_%GPartComm_dtor()
    CALL_TYPE_FREE(this%iotype_,this%ierr_)

    IF ( ALLOCATED    (this%px_) ) DEALLOCATE(this%px_)
    IF ( ALLOCATED    (this%py_) ) DEALLOCATE(this%py_)
    IF ( ALLOCATED    (this%pz_) ) DEALLOCATE(this%pz_)
    IF ( ALLOCATED   (this%pid_) ) DEALLOCATE(this%pid_ )
    IF ( ALLOCATED (this%ptmp0_) ) DEALLOCATE(this%ptmp0_)
    IF ( ALLOCATED  (this%ptmp1) ) DEALLOCATE(this%ptmp1_)
    IF ( ALLOCATED (this%ltmp1_) ) DEALLOCATE(this%ltmp1_)
    IF ( ALLOCATED (this%ltmp2_) ) DEALLOCATE(this%ltmp2_)
    IF ( ALLOCATED     (this%X_) ) DEALLOCATE(this%X_)
    IF ( ALLOCATED     (this%Y_) ) DEALLOCATE(this%Y_)
    IF ( ALLOCATED     (this%Z_) ) DEALLOCATE(this%Z_)
    IF ( ALLOCATED    (this%lx_) ) DEALLOCATE(this%lx_)
    IF ( ALLOCATED    (this%ly_) ) DEALLOCATE(this%ly_)
    IF ( ALLOCATED    (this%lz_) ) DEALLOCATE(this%lz_)
    IF ( ALLOCATED (this%sbbuf_) ) DEALLOCATE(this%sbbuff_)
    IF ( ALLOCATED (this%stbuf_) ) DEALLOCATE(this%stbuff_)
  
  END SUBROUTINE GPart_dtor
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GPart_init(this)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : Init
!  DESCRIPTION: Initializes particle locations before integration.
!               Call after construction.
!  ARGUMENTS  :
!    this    : 'this' class instance
!-----------------------------------------------------------------
    IMPLICIT NONE
    CLASS(GPart)                      :: this

    IF      ( inittype_ .EQ. GPINIT_RANDLOC ) THEN
      CALL GPart_InitRandSeed ()   
    ELSE IF ( inittype_ .EQ. GPINIT_USERLOC ) THEN
      CALL GPart_InitUserSeed ()
    ENDIF

  END SUBROUTINE GPart_init
!-----------------------------------------------------------------
!-----------------------------------------------------------------


  SUBROUTINE GPart_SetTimeOrder(this, iorder)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : SetTimeOrder
!  DESCRIPTION: Sets time step order
!  ARGUMENTS  :
!    this    : 'this' class instance
!    iorder  : value of time step order
!-----------------------------------------------------------------
    IMPLICIT NONE
    CLASS(GPart)                      :: this
    INTEGER, INTENT(IN)               :: iorder

    this%itorder_ = iorder;
   
  END SUBROUTINE GPart_SetTimeOrder
!-----------------------------------------------------------------
!-----------------------------------------------------------------


  INTEGER FUNCTION GPart_GetTimeOrder(this)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : GetTimeOrder
!  DESCRIPTION: Gets time step order
!  ARGUMENTS  :
!    this    : 'this' class instance
!-----------------------------------------------------------------
    IMPLICIT NONE
    CLASS(GPart)                      :: this

    GPart_GetTimeOrder = this%itorder_ 
   
  END SUBROUTINE GPart_GetTimeOrder
!-----------------------------------------------------------------
!-----------------------------------------------------------------


  SUBROUTINE GPart_SetRandSeed(this, iseed)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : SetRandSeed
!  DESCRIPTION: Sets seed for random number generator
!  ARGUMENTS  :
!    this    : 'this' class instance
!    seed    : value of seed
!-----------------------------------------------------------------
    IMPLICIT NONE
    CLASS(GPart)                      :: this
    INTEGER, INTENT(IN)               :: iseed

    this%seed_ = iseed;
   
  END SUBROUTINE GPart_SetRandSeed
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  INTEGER FUNCTION GPart_GetRandSeed(this)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : GetRandSeed
!  DESCRIPTION: Gets random seed
!  ARGUMENTS  :
!    this    : 'this' class instance
!-----------------------------------------------------------------
    IMPLICIT NONE
    CLASS(GPart)                      :: this

    GPart_GetRandSeed = this%seed_
   
  END SUBROUTINE GPart_GetRandSeed
!-----------------------------------------------------------------
!-----------------------------------------------------------------


  SUBROUTINE GPart_SetInitType(this, itype)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : SetInitType
!  DESCRIPTION: Sets particle seed initialization method
!  ARGUMENTS  :
!    this    : 'this' class instance
!    itype   : value of type GPINIT
!-----------------------------------------------------------------
    IMPLICIT NONE
    CLASS(GPart)                      :: this
    TYPE(GPINIT), INTENT(IN)          :: itype

    this%inittype_ = itype;
   
  END SUBROUTINE GPart_SetInitType
!-----------------------------------------------------------------
!-----------------------------------------------------------------


  GPINIT FUNCTION GPart_GetInitType(this, itype)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : GetInitType
!  DESCRIPTION: Gets particle seed initialization method
!  ARGUMENTS  :
!    this    : 'this' class instance
!-----------------------------------------------------------------
    IMPLICIT NONE
    CLASS(GPart)                      :: this

    GPart_GetInitType = this%inittype_ 
   
  END SUBROUTINE GPart_GetInitType
!-----------------------------------------------------------------
!-----------------------------------------------------------------


  SUBROUTINE GPart_SetSeedFile(this, sname)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : SetSeedFile
!  DESCRIPTION: Sets particle seed file name. Note: same file name
!               can be used for user-specified seeds, via InitUserSeed
!               or fixed-location seeds via InitFixedSeed.
!  ARGUMENTS  :
!    this    : 'this' class instance
!    sname   : character file name
!-----------------------------------------------------------------
    USE var
    USE random
    IMPLICIT NONE
    CLASS(GPart)                      :: this
    CHARACTER(len=*) INTENT(IN)       :: sname

    this%seedfile_ = sname;
   
  END SUBROUTINE GPart_SetSeedFile
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  CHARACTER(*) FUNCTION GPart_GetSeedFile(this)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : GetSeedFile
!  DESCRIPTION: Gets particle seed file name. 
!  ARGUMENTS  :
!    this    : 'this' class instance
!-----------------------------------------------------------------
    USE var
    USE random
    IMPLICIT NONE
    CLASS(GPart)                      :: this

    GPart_GetSeedFile = this%seedfile_

  END SUBROUTINE GPart_GetSeedFile
!-----------------------------------------------------------------
!-----------------------------------------------------------------


  SUBROUTINE GPart_InitRandSeed(this)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : InitRandSeed
!  DESCRIPTION: Initializes particle locations by dividing
!               maxparts evenly among MPI slabs, and randomly 
!               selecting positions within a slab.
!  ARGUMENTS  :
!    this    : 'this' class instance
!-----------------------------------------------------------------
    USE var
    USE random
    USE mpivars
    USE fprecision

    IMPLICIT NONE
    CLASS(GPart)                      :: this

    INTEGER                           :: j
    REAL(KIND=GP)                     :: x,r
  
    ! Note: Box is [0,N-1]^3
    this%nparts_ = this%maxparts_/this%nprocs_
    DO j = 1, this%nparts_
       r         = 0.5*(randu(seed_)+1.0)
       px_(j)    = min(int(r*(this%nd_(1)-1)+0.5_GP),real(this%nd_(1)-1,kind=GP) )
       r         = 0.5*(randu(seed_)+1.0)
       py_(j)    = min(int(r*(this%nd_(2)-1)+0.5_GP),real(this%nd_(2)-1,kind=GP) )
       r         = 0.5*(randu(seed_)+1.0)
       x         = real(ksta-1,kind=GP) &
                 + min(int(ksta+r*(kend-ksta+1)+0.5_GP),real(kend-1,kind=GP) )
       pz_(j)    = max(x,ksta-1)
    ENDDO

  END SUBROUTINE GPart_InitRandSeed
!-----------------------------------------------------------------
!-----------------------------------------------------------------


  SUBROUTINE GPart_InitUserSeed(this)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : InitUserSeed
!  DESCRIPTION: Initializes particle locations by letting user
!               specify particle number and locations in
!               file in member data seedfile_.
!  ARGUMENTS  :
!    this    : 'this' class instance
!-----------------------------------------------------------------
    USE mpivars
    USE fprecision
    IMPLICIT NONE
    CLASS(GPart)                      :: this

    INTEGER                           :: navg,nl,nt
    INTEGER,ALLOCATABLE,DIMENSION(:)  :: iproc,ilproc
    REAL(KIND=GP)                     :: x,y,z,zmin,zmax

    zmin  = real(kind=GP,ksta-1)
    zmax  = real(kind=GP,kend-1)

    ! Note: each record (line) consists of x y z real positions
    ! within [0,N-1]^3 box
    OPEN(UNIT=1,FILE=trim(this%seedfile_),STATUS='OLD',ACTION='READ',&
         IOSTAT=this%ierr_, IOMSG=this%serr_);
    IF ( this%ierr_ .NE. 0 ) THEN
      WRITE(*,*)'GPart::InitUserSeed: file:',this%seedfile_,' err: ', trim(this%serr_) 
      STOP
    ENDIF

    ALLOCATE(ilproc(this%nprocs_))
    ALLOCATE(iproc (this%nprocs_))
    ! Find no. parts on each MPI task; only read up to maxparts:
    nt          = 0 
    this%iproc_ = 0
    DO WHILE ( this%ierr_.EQ.0 .AND. nt.LT.this%maxparts_ )
      READ(1,*,IOSTAT=this%ierr_) x, y, z
      IF ( this%ierr_ .NE. 0 ) EXIT
      IF ( z.LE.zmax .AND. z.GE.zmin ) nowned = nowned + 1
      ilproc(this%myrank_+1) = ilproc(this%myrank_+1) + 1 
      nt = nt + 1
    ENDDO
    CLOSE(1)
    CALL MPI_ALLREDUCE(ilproc,iproc,this%nprocs_,GC_INTEGER,      &
                        MPI_SUM,this%comm_,this%ierr_)

    DEALLOCATE(ilproc)
! This check is a dog...
!!   nparts = isumarrp(iprocs,this%nprocs_)
!!   DO WHILE ( nparts.GT.this%maxparts_ )
!!     navg = iavgarrp(iprocs,this%nprocs_) 
!!     j = 1
!!     DO WHILE ( nparts.GT.this%maxparts_ .AND. j.LT.this%nprocs_ )
!!       IF ( iprocs(j).GT.navg .AND. iprocs(j).GE.1 ) THEN
!!         iprocs(j) = iprocs(j) - 1
!!         nparts = isumarrp(iprocs,this%nprocs_)
!!       ENDIF
!!       j = j + 1
!!     ENDDO
!!   ENDDO
    
    OPEN(UNIT=1,FILE=trim(this%seedfile_),STATUS='OLD',ACTION='READ',&
         IOSTAT=this%ierr_, IOMSG=this%serr_);
    IF ( this%ierr_ .NE. 0 ) THEN
      WRITE(*,*)'GPart::InitUserSeed: file:',this%seedfile_,' err: ', trim(this%serr_) 
      STOP
    ENDIF

    nt = 1  ! global particle counter
    nl = 0  ! local particle counter
    DO WHILE ( this%ierr_.EQ.0 .AND. nt.LT.this%maxparts_ )
      READ(1,*,IOSTAT=this%ierr_) x, y, z
      IF ( this%ierr_ .NE. 0 ) EXIT
      IF ( z.LE.zmax .AND. z.GE.zmin ) THEN
        nl = nl + 1
        this%id_(nl) = nt
        this%px_(nl) = x
        this%py_(nl) = y
        this%pz_(nl) = z
      ENDIF
      nt = nt + 1
    ENDDO
    CLOSE(1)
    this%nparts_ = nl;

    DEALLOCATE(iproc )

  END SUBROUTINE GPart_InitUserSeed
!-----------------------------------------------------------------
!-----------------------------------------------------------------


  SUBROUTINE GPart_io_write_euler(this, iunit, dir, fname, nmb, time, evar)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : io_write_euler
!  DESCRIPTION: Converts specified Eulerian real-space variable to
!               a Lagrangian quantity by converting to particle positions;
!               does write of Lagrangian variable to file.
!  ARGUMENTS  :
!    this    : 'this' class instance
!    iunit   : unit number
!    dir     : output directory
!    fname   : filename prefix
!    nmb     : time index
!    time    : real time
!    evar    : Eulerian data from which to compute Lagrangian 
!              quantity: theta(y) = theta(x(y),t). Interpolation
!              of evar is done internally before write.
!    
!-----------------------------------------------------------------
    USE fprecision
    USE commtypes
    USE mpivars

    IMPLICIT NONE
    CLASS(GPart)                      :: this
    REAL(KIND=GP),INTENT(IN)          :: evar(n,n,ksta:kend)
    REAL(KIND=GP),INTENT(IN)          :: time
    INTEGER, INTENT(IN)               :: iunit
    INTEGER                           :: fh,offset,nt,szint,szreal
    CHARACTER(len=100), INTENT(IN)    :: dir
    CHARACTER(len=*), INTENT(IN)      :: nmb
    CHARACTER(len=*), INTENT(IN)      :: fname

    INTEGER                           :: i
    CHARACTER(len=1024)               :: sfile

    CALL GPart_EulerToLag(evar,this%pdb_,this%ltmp1_,this%nparts_)

    CALL MPI_FILE_OPEN(MPI_COMM_WORLD,trim(dir) // '/' // fname // &
         '.' // nmb // '.lag',MPI_MODE_CREATE+MPI_MODE_WRONLY, &
          MPI_INFO_NULL,fh,ioerr)

    ! Must write part. data to correct spot in file:
    CALL MPI_TYPE_SIZE(MPI_INTEGER,szint)
    CALL MPI_TYPE_SIZE(GC_REAL    ,szreal)
    CALL MPI_ALLREDUCE(this%nparts_,nt,this%nprocs_,GC_INTEGER,      &
                        MPI_SUM,this%comm_,this%ierr_)
    IF ( this%myrank_ .EQ. 0 ) THEN
        CALL MPI_FILE_WRITE_AT(fh,0,nt,1,MPI_INTEGER,this%istatus_,this%ierr_)
    ENDIF
    IF ( this%nparts_ .GT. 0 ) THEN
      DO j = 1, this%nparts_
        offset = this%pdb_(j)%id*szreal+szint
        CALL MPI_FILE_WRITE_AT_ALL(fh,offset,this%ltmp1(j),1,GC_REAL,this%istatus_,this%ierr_)
      ENDDO
    ELSE
        CALL MPI_FILE_WRITE_AT_ALL(fh,0,this%ltmp1_(j),0,GC_REAL,this%istatus_,this%ierr_)
    ENDIF
    CALL MPI_FILE_CLOSE(fh,this%ierr_)

    
  END SUBROUTINE GPart_io_write_euler
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GPart_io_write_pdb(this, iunit, dir, spref, nmb, time)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : io_write_pdb
!  DESCRIPTION: Does write of Lagrangian position d.b. to file. 
!               Position of the particle structure in file is the
!               particles id.
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
    CLASS(GPart)                      :: this
    REAL(KIND=GP),INTENT(IN)          :: time
    REAL(KIND=GP)                     :: prec(3)
    INTEGER,INTENT(IN)                :: iunit
    INTEGER                           :: fh,offset,nt,szint
    CHARACTER(len=*), INTENT(IN)      :: dir
    CHARACTER(len=*), INTENT(IN)      :: nmb
    CHARACTER(len=*), INTENT(IN)      :: fname
    TYPE(GPDBrec)                     :: pst

    INTEGER                           :: j

!!  IF ( myrank .EQ. 0 ) THEN
!!    this%sfile_ = trim(dir) // '/' // fname // '.' // &
!!            nmb // '.lag'
!!    OPEN(iunit,file=sfile_,form='unformatted',IOSTAT=ierr_,IOMSG=serr_)

!!    IF ( ierr_ .NE. 0 ) THEN
!!      WRITE(*,*)'GPart::io_write_pdb: file:', sfile_,' err: ', trim(serr_) 
!!    ENDIF

!!    WRITE(iunit) time, real(nparts_,KIND=GP)
!!    DO j = 1, 3
!!      WRITE(iunit) pdb_(1:nparts_,j)
!!    ENDDO

!!    CLOSE(iunit)

!!  ENDIF

    CALL MPI_TYPE_SIZE(MPI_INTEGER,szint)
    CALL MPI_FILE_OPEN(MPI_COMM_WORLD,trim(dir) // '/' // fname // &
         '.' // nmb // '.lag',MPI_MODE_CREATE+MPI_MODE_WRONLY, &
          MPI_INFO_NULL,fh,ioerr)

    ! Must write part. data to correct spot in file:
    CALL MPI_ALLREDUCE(this%nparts_,nt,this%nprocs_,GC_INTEGER,      &
                        MPI_SUM,this%comm_,this%ierr_)
    IF ( this%myrank_ .EQ. 0 ) THEN
        CALL MPI_FILE_WRITE_AT(fh,0,nt,1,MPI_INTEGER,this%istatus_,this%ierr_)
    ENDIF
    IF ( this%nparts_ .GT. 0 ) THEN
      DO j = 1, this%nparts_
        offset = this%pdb_(j)%id*this%this%sziotype_+szint
        prec(1) = this%pdb_(j)%x
        prec(2) = this%pdb_(j)%y
        prec(3) = this%pdb_(j)%z
        CALL MPI_FILE_WRITE_AT_ALL(fh,offset,prec,1,this%iotype_,this%ierr_)
      ENDDO
    ELSE
        CALL MPI_FILE_WRITE_AT_ALL(fh,0,prec,0,this%iotype_,this%ierr_)
    ENDIF
    CALL MPI_FILE_CLOSE(fh,this%ierr_)



  END SUBROUTINE GPart_io_write_pdb
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GPart_io_read(this, iunit, dir, spref, nmb, lvar)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : io_read
!  DESCRIPTION: does read of Lagrangian particle data from file
!  ARGUMENTS  :
!    this    : 'this' class instance
!    iunit   : unit number
!    dir     : input directory
!    spref   : filename prefix
!    nmb     : time index
!    lvar    : Lag. variable, out
!    
!-----------------------------------------------------------------
    USE fprecision
    USE commtypes
    USE mpivars

    IMPLICIT NONE
    CLASS(GPart)                              :: this
    REAL(KIND=GP),INTENT(OUT),DIMENSION(:)    :: lvar
    REAL(KIND=GP)                             :: zmin,zmax
    INTEGER,INTENT(IN)                        :: iunit, bsynch
    INTEGER                                   :: fh,j,nt,szint
    CHARACTER(len=*), INTENT(IN)              :: dir
    CHARACTER(len=*), INTENT(IN)              :: nmb
    CHARACTER(len=*), INTENT(IN)              :: fname


    zmin  = real(kind=GP,ksta-1)
    zmax  = real(kind=GP,kend-1)
    CALL MPI_TYPE_SIZE(MPI_INTEGER,szint)
    CALL MPI_FILE_OPEN(MPI_COMM_WORLD,trim(dir) // '/' // fname // &
         '.' // nmb // '.lag',MPI_MODE_CREATE+MPI_MODE_WRONLY, &
          MPI_INFO_NULL,fh,ioerr)

    ! Must write part. data to correct spot in file:
    CALL MPI_FILE_READ_AT_ALL(fh,0,nt,1,MPI_INTEGER,this%istatus_,this%ierr_)
    IF ( nt.GT.this%maxparts_ ) THEN
      WRITE(*,*) 'GPart_io_read: Attempt to read too many particles'
      STOP
    ENDIF
    CALL MPI_FILE_READ_AT_ALL(fh,szint,this%ptmp0_,nt*3,1,GC_REAL,this%istatus_,this%ierr_)
    CALL MPI_FILE_CLOSE(fh,this%ierr_)
    
    ! Remove from local pdb the parts that don't belong:
    DO j = 1, nt
      this%pdb_(j).id = j-1
      this%pdb_(j).x  = this%ptmp0_(1,j)
      this%pdb_(j).y  = this%ptmp0_(2,j)
      this%pdb_(j).z  = this%ptmp0_(3,j)
      IF ( pdb(j)%z.LT.zmin .OR. pdb(j)%z.GT.zmax) THEN 
        this%pdb_(j).id = GPNULL
      ENDIF
    ENDDO
    CALL GPart_Delete(this,this%pdb_,nt,this%nparts_)

  END SUBROUTINE GPart_io_read
!-----------------------------------------------------------------
!-----------------------------------------------------------------


  SUBROUTINE GPart_SetStepRKK(this)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     :SetStepRKK
!  DESCRIPTION: Initializes an explicit integration timstep.
!               Must be called at the start of the RK stage execution,.

!  ARGUMENTS  :
!    this    : 'this' class instance
!-----------------------------------------------------------------
    USE fprecision
    USE commtypes
    USE mpivars

    IMPLICIT NONE
    CLASS(GPart)                                      :: this

    INTEGER                                           :: i,j

    ! Find array of indices that point into master d.b.
    ! that this MPI task is responsible for (in iwrk_, nwrk_):
!!  CALL GPart_GetWrkInd(pdb_,nparts_)
    
    ! Initialize solution, u: 
    ! u* <-- u: 
 
    ! Cycle over JST loop to update state:
    DO j = 1,3
       DO i = 1, this%nparts_
        ptmp0_(1,i) = pdb_(i)%x  ! u_0
        ptmp0_(2,i) = pdb_(i)%y  ! u_0
        ptmp0_(3,i) = pdb_(i)%z  ! u_0
        ptmp1_(1,i) = pdb_(i)%x  ! u*
        ptmp1_(2,i) = pdb_(i)%y  ! u*
        ptmp1_(3,i) = pdb_(i)%z  ! u*
      ENDDO 
    ENDDO 
    

  END SUBROUTINE GPart_SetStepRKK
!-----------------------------------------------------------------
!-----------------------------------------------------------------


  SUBROUTINE GPart_FinalStepRKK(this)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : FinalStepRKK
!  DESCRIPTION: Called at the end of all RK-like stages to
!               complete Lagrangian particle update.

!  ARGUMENTS  :
!    this    : 'this' class instance
!-----------------------------------------------------------------
    USE fprecision
    USE commtypes
    USE mpivars

    IMPLICIT NONE
    CLASS(GPart)                                      :: this

    ! u(t+dt) = u*:

  END SUBROUTINE GPart_FinalStep
!-----------------------------------------------------------------
!-----------------------------------------------------------------


  SUBROUTINE GPart_stepRKK(this, vx, vy, vz, dt, xk)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : stepRKK
!  DESCRIPTION: Integrates from current state over interval dt
!               for explicit step within an outer stepper
!               method of the form:
!
!               X = X_0 + dt * V(X(t),t) * xk,

!  ARGUMENTS  :
!    this    : 'this' class instance
!    vz,vy,vz: compoments of velocity field, in real space, partially
!              updated, possibly
!    dt      : integration timestep
!    xk      : multiplicative RK factor
!-----------------------------------------------------------------
    USE fprecision
    USE commtypes
    USE mpivars

    IMPLICIT NONE
    CLASS(GPart)                                      :: this

    REAL(KIND=GP),INTENT(IN),DIMENSION(n,n,ksta:kend) :: vx,vy,vz
    REAL(KIND=GP),INTENT(IN)                          :: dt,xk

    INTEGER                                           :: i,j

!   Currently, only allows JST-style explicit RK of 
!   order torder:

    ! Exchange bdy data for velocities, so that we
    ! can perform local interpolations:
    CALL this%exchop_%SlabDataExchange(this,vx)

    ! Find F(u*):
    ! ... x:
    CALL GPart_EulerToLag(vx,ltmp1_)
    ! ux* <-- ux + dt * F(U*)*xk:
    DO j = 1, nwrk_
      ptmp1_(1,j) = ptmp0_(1,j) + dt*ltmp1_(j)*xk
    ENDDO

    ! ... y:
    ! Exchange bdy data for velocities, so that we
    ! can perform local interpolations:
    CALL this%exchop_%SlabDataExchange(this,vz)
    CALL GPart_EulerToLag(vy,ltmp1_)
    ! uy* <-- uy + dt * F(U*)*xk:
    DO j = 1, nwrk_
      ptmp1_(2,j) = ptmp0_(2,j) + dt*ltmp1_(j)*xk
    ENDDO

    ! ... z:
    ! Exchange bdy data for velocities, so that we
    ! can perform local interpolations:
    CALL this%exchop_%SlabDataExchange(this,vz)
    CALL GPart_EulerToLag(vz,ltmp1__)
    ! uz* <-- uz + dt * F(U*)*xk:
    DO j = 1, nwrk_
      ptmp1_(3,j) = ptmp0_(3,j) + dt*ltmp1_(j)*xk
    ENDDO
!
!   in middle of RK stages, should not need to enforce periodicity,
!   or redistribute work, if we're Courant-limted, right? 

    ! Enforce periodicity:
    CALL GPart_MakePeriodicP(this,ptmp1_,nwrk_)

  END SUBROUTINE GPart_stepRKK
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GPart_GetWrkInd(this,pdb,npdb)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : GetWrkInd
!  DESCRIPTION: Finds this MPI task's work, by filling iwrk_ and 
!               nwrk_ data members based on db, pdb
!  ARGUMENTS  :
!    this    : 'this' class instance
!    pdb     : particle d.b. to consider
!    npdb    : no. particles in pdb
!    
!-----------------------------------------------------------------
    USE fprecision
    USE commtypes
    USE mpivars

    IMPLICIT NONE
    CLASS(GPart)                              :: this
    REAL(KIND=GP),INTENT(IN)                  :: pdb(this%maxparts_,3)
    INTEGER,INTENT(IN)                        :: npdb

    INTEGER                                   :: j
   
    ! Cycle through master partle d.b., and locate all particles
    ! that lie in this task's slab:
    nwrk_ = 0
    iwrk_ = 0
    DO j = 1, npdb
      IF ( pdb(j).z.GE.real(ksta,KIND=GP)      &
     .AND. pdb(j).z.LT.real(kend,KIND=GP) ) THEN
        nwrk_ = nwrk_ + 1
        iwrk_(nwrk_) = j
      ENDIF
    ENDDO

  END SUBROUTINE GPart_GetWrkInd
!-----------------------------------------------------------------
!-----------------------------------------------------------------


  SUBROUTINE GPart_EulerToLag(this,evar,lag)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : EulerToLag
!  DESCRIPTION: Computes the Eulerian to Lagrangian
!               transformation by interpolating Eulerian field
!               evar to position of Lagrangian paricles in pdb. 
!               Global set of Lagrangian particles may be considered
!               in this method, as the locations in lag correspond
!               one-to-one with the locations in the pdb d.b..
!  ARGUMENTS  :
!    this    : 'this' class instance (IN)
!    evar    : Eulerian variable (IN)
!    lag     : real array with the containing Lagrangian value
!              of evar field on output (OUT)
!    
!-----------------------------------------------------------------
    USE fprecision
    USE commtypes
    USE mpivars

    IMPLICIT NONE
    CLASS(GPart)                              :: this

    REAL(KIND=GP),INTENT(IN)                  :: evar(n,n,ksta:kend)
    REAL(KIND=GP),INTENT(OUT)                 :: lag(*)

    INTEGER                                   :: j


  END SUBROUTINE GPart_EulerToLag
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GPart_EulerToLocalLag(this,evar,lag)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : EulerToLocalLag
!  DESCRIPTION: Computes the Eulerian to (local) Lagrangian
!               transformation by interpolating Eulerian field
!               evar to position of Lagrangian paricles in pdb. 
!               Only Lagrangian particles owned by this MPI task
!               may be considered here. These are specified by
!               the index array, ilocal.
!  ARGUMENTS  :
!    this    : 'this' class instance (IN)
!    evar    : Eulerian variable (IN)
!    lag     : real array with the containing Lagrangian value
!              of evar field on output (OUT)
!    
!-----------------------------------------------------------------
    USE fprecision
    USE commtypes
    USE mpivars

    IMPLICIT NONE
    CLASS(GPart)                              :: this

    REAL(KIND=GP),INTENT(IN)                  :: evar(n,n,ksta:kend)
    REAL(KIND=GP),INTENT(IN)                  :: gpdb(this%maxparts_,3)
    REAL(KIND=GP),INTENT(OUT)                 :: lag(*)
    INTEGER,INTENT(IN)                        :: ilocal(*)
    INTEGER,INTENT(IN)                        :: np

    INTEGER                                   :: j

  END SUBROUTINE GPart_EulerToLocalLag
!-----------------------------------------------------------------
!-----------------------------------------------------------------


  SUBROUTINE GPart_MakePeriodicP(this,pdb,npdb)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : MakePeriodicP
!  DESCRIPTION: Enforces periodic b.c.'s on particles in pdb
!  ARGUMENTS  :
!    this    : 'this' class instance (IN)
!    pdb     : particle loc. d.b.
!    npdb    : no. particles in pdb
!    
!-----------------------------------------------------------------
    USE fprecision
    USE commtypes
    USE mpivars

    IMPLICIT NONE
    CLASS(GPart)                                     :: this
    REAL(KIND=GP),INTENT(INOUT)                      :: pdb(3,this%maxparts_)
    INTEGER,INTENT(IN)                               :: npdb

    REAL(KIND=GP)                                    :: del
    INTEGER                                          :: j

    del = gxbnds_(i,2)-gxbnds_(i,1)
    DO i = 1, 3
      DO j = 1, npdb
        IF ( pdb(i,j).LT.gxbnds_(i,1) ) pdb(i,j) = del-pdb(i,j)        
        IF ( pdb(i,j).GT.gxbnds_(i,2) ) pdb(i,j) = pdb(i,j)-del
      ENDDO
    ENDDO

  END SUBROUTINE GPart_MakePeriodicP
!-----------------------------------------------------------------
!-----------------------------------------------------------------


  SUBROUTINE GPart_MakePeriodicExt(this,v,nx,ny,kb,ke,nc)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : MakePeriodicExt
!  DESCRIPTION: Enforces periodic b.c.'s on (extended) field
!  ARGUMENTS  :
!    this    : 'this' class instance (IN)
!    v       : real field on extended grid 
!    nx,ny   : global size of field in x-y (including ghost zones)
!    kb,ke   : starting, ending z-indices of slab (including ghost zones)
!    nc      : index in x and y (and z) s.t. f(nc+1,:,:) = f(nx-nc,:,:), etc.
!    
!-----------------------------------------------------------------
    USE fprecision
    USE commtypes

    IMPLICIT NONE
    CLASS(GPart)                                     :: this
    INTEGER,INTENT(IN)                               :: nx,ny,kb,ke
    REAL(KIND=GP),INTENT(INOUT)                      :: v(nx,ny,kb:ke)

    INTEGER                                          :: j

    ! Recall: nx, ny are the dimensions _including_ ghost zones:
    !
    ! Periodicity s.t.:
    !   | | [ | | | | ] | |
    !   a b         a b 
    DO k = kb,ke 
      DO j = 1,ny
        DO i = 1,nc
          v(i,j,k) = v(nx-nc+i,j,k)
        ENDDO
        DO i = nx-nc+1,nx
          v(i,j,k) = v(2*nc+i-nx,j,k)
        ENDDO
      ENDDO
      DO i = 1,nx
        DO j = 1,nc
          v(i,j,k) = v(i,nx-nc+j,k)
        ENDDO
        DO j = ny-nc+1,ny
          v(i,j,k) = v(i,2*nc+j-nx,k)
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE GPart_MakePeriodicExt
!-----------------------------------------------------------------
!-----------------------------------------------------------------


  SUBROUTINE GPart_Copy2Ext(this,ve,v)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : GPart_Copy2Ext
!  DESCRIPTION: Copies field from normal to extended grids.
!  ARGUMENTS  :
!    this    : 'this' class instance (IN)
!    ve      : real field on extended grid 
!    v       : real field on normal grid 
!    nc      : index in x and y (and z) s.t. f(nc+1,:,:) = f(nx-nc,:,:), etc.
!    
!-----------------------------------------------------------------
    USE fprecision
    USE commtypes

    IMPLICIT NONE
    CLASS(GPart)                                     :: this
    INTEGER                                          :: ng,nxy,nexy
    REAL(KIND=GP),INTENT(INOUT)                      :: ve(*),v(*)

    INTEGER                                          :: i,j,k

    ! Recall: nx, ny are the dimensions _including_ ghost zones:
    !
    ! Periodicity s.t.:
    !   | | [ | | | | ] | |
    !   a b         a b 
    ng  = this%intorder_
    nxy = this%nd_(1)*this%nd_(2)
    nexy= (this%nd_(1)+ng)*(this%nd_(2)+ng)

    DO k = 1,this%nd_(3)
      DO j = 1, this%nd_(2)
        DO i = 1, this%nd_(1)
          ve(i+ng+(j+ng)*this%nd_(1)+(k+ng)*nexy) = v(i+j*this%nd_(1)+k*nxy)
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE GPart_Copy2Ext
!-----------------------------------------------------------------
!-----------------------------------------------------------------


  SUBROUTINE GPart_Delete(this,pdb,npdb,nnew)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : GPart_Delete
!  DESCRIPTION: Removes from PDB NULL particles, concatenates list,
!               and sets new number of particles
!  ARGUMENTS  :
!    this    : 'this' class instance (IN)
!    pdb     : part. d.b.
!    npdb    : no. parts. in pdb
!    nnew    : no. non-NULL particles (set in GPartComm class)
!    
!-----------------------------------------------------------------
    USE fprecision
    USE commtypes

    IMPLICIT NONE
    CLASS(GPart)                                     :: this
    TYPE(GPDBrec),INTENT(INOUT),DIMENSION(*)         :: pdb
    INTEGER                                          :: ng,nxy,nexy
    INTEGER,INTENT (IN)                              :: npdb
    INTEGER,INTENT(OUT)                              :: nnew
    INTEGER                                          :: i,j

    nnew = 0
    DO i = 1, npdb
       IF ( pdb(i)%id .NE. GPNULL ) nnew = nnew + 1
    ENDDO

    j     = 1
    DO i = 1, nnew
      DO WHILE ( j.LE.npdb .AND. id(j).EQ.GPNULL )
        j = j + 1
      ENDDO
      IF ( j.LE.npdb .AND. j.NE.i ) THEN
        pdb(i)%id = pdb(j)%id; pdb(j)%id = GPNULL
        pdb(i)%x = pdb(j)%x
        pdb(i)%y = pdb(j)%y
        pdb(i)%z = pdb(j)%z
      ENDIF

    ENDDO

  END SUBROUTINE GPart_Delete
!-----------------------------------------------------------------
!-----------------------------------------------------------------


END MODULE class_GPart
