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
      USE class_GPartComm
      USE class_GPSplineInt

      IMPLICIT NONE
      
!!    ENUM, BIND(C) :: GPINIT
!!      ENUMERATOR :: GPINIT_RANDLOC=0
!!      ENUMERATOR :: GPINIT_USERLOC
!!    ENDENUM GPINIT

!!    ENUM, BIND(C) :: GPINTRP
!!      ENUMERATOR :: GPINTRP_CSPLINE=0
!!      ENUMERATOR :: GPINTRP_LAGINT 
!!      ENUMERATOR :: GPINTRP_SPECTRAL
!!    ENDENUM GPINTRP

!!    ENUM, BIND(C) :: GPEXCHTYPE
!!      ENUMERATOR :: GPEXCHTYPE_NN =0
!!      ENUMERATOR :: GPEXCHTYPE_VDB
!!    ENDENUM GPINTRP

      INTEGER,PARAMETER                              :: GPINIT_RANDLOC =0
      INTEGER,PARAMETER                              :: GPINIT_USERLOC =1

      INTEGER,PARAMETER                              :: GPINTRP_CSPLINE=0
      INTEGER,PARAMETER                              :: GPINTRP_LAGINT =1

      INTEGER,PARAMETER                              :: GPEXCHTYPE_NN  =0
      INTEGER,PARAMETER                              :: GPEXCHTYPE_VDB =1
       

      TYPE, PUBLIC :: GPart
        PRIVATE
        ! Member data:
        INTEGER, DIMENSION(MPI_STATUS_SIZE)          :: istatus_
        INTEGER                                      :: inittype_
        INTEGER                                      :: iinterp_
        INTEGER                                      :: iexchtype_
        TYPE(GPartComm)                              :: exchop_
        TYPE(GPSplineInt)                            :: intop_
        INTEGER                                      :: intorder_,itorder_,nd_(3),libnds_(3,2)
        INTEGER                                      :: myrank_,nprocs_
        INTEGER                                      :: ierr_,iseed_
        INTEGER                                      :: maxparts_,nparts_,nvdb_
        INTEGER                                      :: iotype_,sziotype_
        INTEGER                                      :: comm_
        INTEGER      , ALLOCATABLE, DIMENSION    (:) :: id_
!!      TYPE(GPDBrec), ALLOCATABLE, DIMENSION    (:) :: pdb_
        REAL(KIND=GP), ALLOCATABLE, DIMENSION    (:) :: px_,py_,pz_
        REAL(KIND=GP), ALLOCATABLE, DIMENSION  (:,:) :: ptmp0_,ptmp1_,vdb_
        REAL(KIND=GP), ALLOCATABLE, DIMENSION    (:) :: ltmp0_,ltmp1_
        REAL(KIND=GP)                                :: lxbnds_(3,2),gxbnds_(3,2),gext_(3)
        CHARACTER(len=1024)                          :: seedfile_,serr_,sfile_
      CONTAINS
        ! Public methods:
        PROCEDURE,PUBLIC :: GPart_ctor
        FINAL            :: GPart_dtor
        PROCEDURE,PUBLIC :: Init            => GPart_Init
        PROCEDURE,PUBLIC :: Step            => GPart_StepRKK
        PROCEDURE,PUBLIC :: SetStep         => GPart_SetStepRKK
        PROCEDURE,PUBLIC :: FinalStep       => GPart_FinalizeRKK
        PROCEDURE,PUBLIC :: io_write_euler  => GPart_io_write_euler
        PROCEDURE,PUBLIC :: io_write_pdb    => GPart_io_write_pdb
        PROCEDURE,PUBLIC :: io_read         => GPart_io_read
!       PROCEDURE,PUBLIC :: GetPos
        PROCEDURE,PUBLIC :: EulerToLag      => GPart_EulerToLag
        PROCEDURE,PUBLIC :: SetInitType     => GPart_SetInitType
        PROCEDURE,PUBLIC :: SetSeedFile     => GPart_SetSeedFile
        PROCEDURE,PUBLIC :: SetRandSeed     => GPart_SetRandSeed
        PROCEDURE,PUBLIC :: SetTimeOrder    => GPart_SetTimeOrder
        PROCEDURE,PUBLIC :: GetSeedFile     => GPart_GetSeedFile
        PROCEDURE,PUBLIC :: GetRandSeed     => GPart_GetRandSeed
        PROCEDURE,PUBLIC :: GetTimeOrder    => GPart_GetTimeOrder

        GENERIC  ,PUBLIC :: io_write        => io_write_euler,io_write_pdb
      END TYPE GPart

      PRIVATE :: GPart_Init        , GPart_StepRKK     , GPart_io_write_euler
      PRIVATE :: GPart_SetStepRKK  , GPart_FinalizeRKK
      PRIVATE :: GPart_io_write_pdb, GPart_io_read     
      PRIVATE :: GPart_InitRandSeed, GPart_InitUserSeed 
      PRIVATE :: GPart_SetInitType , GPart_SetSeedFile
      PRIVATE :: GPart_SetRandSeed , GPart_Delete
!     PRIVATE :: GPart_GetParts    , GPart_MakePeriodicP
      PRIVATE :: GPart_MakePeriodicP
      PRIVATE :: GPart_GetLocalWrk , GPart_MakePeriodicExt

! Methods:
  CONTAINS

  SUBROUTINE GPart_ctor(this, comm, mparts, inittype, iinterp, intorder, iexchtyp, csize, nstrip)
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
!    iexchtyp: format for exchanging particle data: GPEXCHTYPE_NN=nearest neighbor, 
!              suggesting that data is passed between neighboring
!              tasks only; GPEXCHTYPE_VDB form suggests that all tasks retain a copy
!              of the 'master' particle d.b., so that passing between
!              neighbors is not necessary. But the VDB form does require
!              a global reduction, and may be more expensive.
!    csize    : cache size param for local transposes
!    nstrip   : 'strip-mining' size for local transposes
!-----------------------------------------------------------------
    USE grid
    USE mpivars
    USE commtypes

    IMPLICIT NONE
    CLASS(GPart)                         :: this
    INTEGER          , INTENT(IN)        :: comm,mparts,csize,nstrip
    INTEGER                              :: disp(3),lens(3),types(3),szreal
    INTEGER          , INTENT(IN)        :: iexchtyp,iinterp,inittype,intorder
    INTEGER                              :: j,nc

    this%nparts_   = 0 
    this%nvdb_     = 0
    this%comm_     = comm
    this%maxparts_ = mparts
    this%nd_(1:3)  = n
    this%seedfile_ = 'gploc.dat'
    this%iinterp_  = 3          ! fixed for now
    this%inittype_ = inittype
    this%itorder_  = 2
    this%intorder_ = max(intorder,1)
    this%iseed_     = 1000
    this%iexchtype_=  iexchtyp

    IF ( this%intorder_ .NE. 3 ) THEN
      WRITE(*,*) 'GPart::ctor: Only 3rd order allowed for now' 
    ENDIF

    CALL MPI_COMM_SIZE(this%comm_,this%nprocs_,this%ierr_)
    CALL MPI_COMM_RANK(this%comm_,this%myrank_,this%ierr_)

    CALL this%exchop_%GPartComm_ctor(1,this%maxparts_,this%nd_,this%intorder_-1,this%comm_)
    CALL this%exchop_%SetCacheParam(csize,nstrip)
    CALL this%exchop_%Init()

    this%libnds_(1,1) = 1    ; this%lxbnds_(1,1) = 0.0
    this%libnds_(1,2) = n    ; this%lxbnds_(1,2) = real(n-1,kind=GP)
    this%libnds_(2,1) = 1    ; this%lxbnds_(2,1) = 0.0
    this%libnds_(2,2) = n    ; this%lxbnds_(2,2) = real(n-1,kind=GP)
    this%libnds_(3,1) = ksta ; this%lxbnds_(3,1) = real(ksta-1,kind=GP)
    this%libnds_(3,2) = kend ; this%lxbnds_(3,2) = real(kend-1,kind=GP)

    DO j = 1,3
      this%gxbnds_(j,1) = 0.0_GP 
      this%gxbnds_(j,2) = real(this%nd_(j)-1,kind=GP)
      this%gext_    (j) = this%gxbnds_(j,2) - this%gxbnds_(j,1)
    ENDDO
    CALL this%intop_%GPSplineInt_ctor(3,this%nd_,this%libnds_,this%lxbnds_,this%maxparts_,this%exchop_)

    ! Create part. d.b. structure type for I/O
    CALL MPI_TYPE_SIZE(GC_REAL,szreal)
    disp (1)=0 
    disp (2)=szreal
    disp (3)=2*szreal
    types(1)=GC_REAL
    types(2)=GC_REAL
    types(3)=GC_REAL
    CALL_TYPE_CREATE_STRUCT(3,lens,disp,this%iotype_,this%ierr_)
    CALL_TYPE_COMMIT(this%iotype_,this%ierr_)
    CALL MPI_TYPE_SIZE(this%iotype_,this%sziotype_,this%ierr_)

    ALLOCATE(this%id_      (this%maxparts_))
    ALLOCATE(this%px_      (this%maxparts_))
    ALLOCATE(this%py_      (this%maxparts_))
    ALLOCATE(this%pz_      (this%maxparts_))
    ALLOCATE(this%ptmp0_ (3,this%maxparts_))
    ALLOCATE(this%ptmp1_ (3,this%maxparts_))
    IF ( this%iexchtype_.EQ.GPEXCHTYPE_VDB ) THEN
      ALLOCATE(this%vdb_(3,this%maxparts_))
    ENDIF
    ALLOCATE(this%ltmp0_ (this%maxparts_))
    ALLOCATE(this%ltmp1_ (this%maxparts_))

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


!!  CALL this%exchop_%GPartComm_dtor()
    CALL MPI_TYPE_FREE(this%iotype_,this%ierr_)

    IF ( ALLOCATED    (this%px_) ) DEALLOCATE(this%px_)
    IF ( ALLOCATED    (this%py_) ) DEALLOCATE(this%py_)
    IF ( ALLOCATED    (this%pz_) ) DEALLOCATE(this%pz_)
    IF ( ALLOCATED    (this%id_) ) DEALLOCATE(this%id_ )
    IF ( ALLOCATED (this%ptmp0_) ) DEALLOCATE(this%ptmp0_)
    IF ( ALLOCATED (this%ptmp1_) ) DEALLOCATE(this%ptmp1_)
    IF ( ALLOCATED   (this%vdb_) ) DEALLOCATE(this%vdb_)
    IF ( ALLOCATED (this%ltmp0_) ) DEALLOCATE(this%ltmp0_)
    IF ( ALLOCATED (this%ltmp1_) ) DEALLOCATE(this%ltmp1_)
  
  END SUBROUTINE GPart_dtor
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GPart_Init(this)
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

    IF      ( this%inittype_ .EQ. GPINIT_RANDLOC ) THEN
      CALL GPart_InitRandSeed (this)   
    ELSE IF ( this%inittype_ .EQ. GPINIT_USERLOC ) THEN
      CALL GPart_InitUserSeed (this)
    ENDIF

  END SUBROUTINE GPart_Init
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


  FUNCTION GPart_GetTimeOrder(this) result(get_res)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : GetTimeOrder
!  DESCRIPTION: Gets time step order
!  ARGUMENTS  :
!    this    : 'this' class instance
!-----------------------------------------------------------------
    IMPLICIT NONE
    CLASS(GPart)                      :: this
    INTEGER                           :: iseed
    INTEGER                           :: get_res

    get_res = this%itorder_ 
   
  END FUNCTION GPart_GetTimeOrder
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

    this%iseed_ = iseed;
   
  END SUBROUTINE GPart_SetRandSeed
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  FUNCTION GPart_GetRandSeed(this) result(get_res)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : GetRandSeed
!  DESCRIPTION: Gets random seed
!  ARGUMENTS  :
!    this    : 'this' class instance
!-----------------------------------------------------------------
    IMPLICIT NONE
    CLASS(GPart)                      :: this
    INTEGER                           :: get_res

    get_res = this%iseed_
   
  END FUNCTION GPart_GetRandSeed
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
    INTEGER     , INTENT(IN)          :: itype

    this%inittype_ = itype;
   
  END SUBROUTINE GPart_SetInitType
!-----------------------------------------------------------------
!-----------------------------------------------------------------


  FUNCTION GPart_GetInitType(this) result(get_res)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : GetInitType
!  DESCRIPTION: Gets particle seed initialization method
!  ARGUMENTS  :
!    this    : 'this' class instance
!-----------------------------------------------------------------
    IMPLICIT NONE
    CLASS(GPart)                      :: this
    INTEGER                           :: get_res

    get_res = this%inittype_ 
   
  END FUNCTION GPart_GetInitType
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
    IMPLICIT NONE
    CLASS(GPart)                      :: this
    CHARACTER(len=*),INTENT(IN)       :: sname

    this%seedfile_ = sname;
   
  END SUBROUTINE GPart_SetSeedFile
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  FUNCTION GPart_GetSeedFile(this) result(get_res)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : GetSeedFile
!  DESCRIPTION: Gets particle seed file name. 
!  ARGUMENTS  :
!    this    : 'this' class instance
!-----------------------------------------------------------------
    IMPLICIT NONE
    CLASS(GPart)                      :: this
    CHARACTER(1024)                   :: get_res

    get_res = this%seedfile_

  END FUNCTION GPart_GetSeedFile
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
    USE random

    IMPLICIT NONE
    CLASS(GPart)                      :: this

    INTEGER                           :: j
    REAL(KIND=GP)                     :: x,r
  
    ! Note: Box is [0,N-1]^3
    this%nparts_ = this%maxparts_/this%nprocs_
    DO j = 1, this%nparts_
       this%id_(j)    = (this%myrank_-1)*this%nparts_ + j - 1
       r         = 0.5*(randu(this%iseed_)+1.0)
       this%px_(j)    = min(r*(this%nd_(1)-1)+0.5_GP,real(this%nd_(1)-1,kind=GP) )
       r         = 0.5*(randu(this%iseed_)+1.0)
       this%py_(j)    = min(r*(this%nd_(2)-1)+0.5_GP,real(this%nd_(2)-1,kind=GP) )
       r         = 0.5*(randu(this%iseed_)+1.0)
       x         = real(ksta-1,kind=GP) &
                 + min(ksta+r*(kend-ksta+1)+0.5_GP,real(kend-1,kind=GP) )
       this%pz_(j)   = max(x,real(ksta-1,kind=GP))
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

    INTEGER                           :: navg,nl,nowned,nt
    INTEGER,ALLOCATABLE,DIMENSION(:)  :: iproc,ilproc
    REAL(KIND=GP)                     :: x,y,z


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
    nt    = 0 
    iproc = 0
    DO WHILE ( this%ierr_.EQ.0 .AND. nt.LT.this%maxparts_ )
      READ(1,*,IOSTAT=this%ierr_) x, y, z
      IF ( this%ierr_ .NE. 0 ) EXIT
      IF ( z.LT.this%lxbnds_(3,2) .AND. z.GT.this%lxbnds_(3,1) ) nowned = nowned + 1
      ilproc(this%myrank_+1) = ilproc(this%myrank_+1) + 1 
      nt = nt + 1
    ENDDO
    CLOSE(1)
    CALL MPI_ALLREDUCE(ilproc,iproc,this%nprocs_,MPI_INTEGER,      &
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

    nt = 0  ! global particle counter
    nl = 0  ! local particle counter
    DO WHILE ( this%ierr_.EQ.0 .AND. nt.LT.this%maxparts_ )
      READ(1,*,IOSTAT=this%ierr_) x, y, z
      IF ( this%ierr_ .NE. 0 ) EXIT
      IF ( z.LT.this%lxbnds_(3,2) .AND. z.GT.this%lxbnds_(3,1) ) THEN
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


  SUBROUTINE GPart_io_write_euler(this, iunit, dir, fname, nmb, &
             time, evar, doupdate, tmp1, tmp2)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : io_write_euler
!  DESCRIPTION: Converts specified Eulerian real-space variable to
!               a Lagrangian quantity by interpolating to particle positions;
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
!              of evar is done internally before write. Note that
!              data in evar is lost on exit.
!    doupdate: if true, do interp point update in interpolator; else don't
!    tmp1/2  : temp arrays of size of evar. Required for interpolation
!-----------------------------------------------------------------
    USE grid
    USE fprecision
    USE commtypes
    USE mpivars

    IMPLICIT NONE
    CLASS(GPart)                                        :: this
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(n,n,ksta:kend):: evar(n,n,ksta:kend)
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(n,n,ksta:kend):: tmp1,tmp2
    REAL(KIND=GP),INTENT   (IN)                         :: time
    INTEGER      ,INTENT   (IN)                         :: iunit
    INTEGER                                             :: fh,offset,nt,szint,szreal
    INTEGER                                             :: j
    LOGICAL      ,INTENT   (IN)                         :: doupdate
    CHARACTER(len=100), INTENT(IN)                      :: dir
    CHARACTER(len=*)  , INTENT(IN)                      :: nmb
    CHARACTER(len=*)  , INTENT(IN)                      :: fname
    CHARACTER(len=1024)                                 :: sfile


    CALL GPart_EulerToLag(this,this%ltmp1_,this%nparts_,evar,doupdate,tmp1,tmp2)

    CALL MPI_FILE_OPEN(MPI_COMM_WORLD,trim(dir) // '/' // fname // &
         '.' // nmb // '.lag',MPI_MODE_CREATE+MPI_MODE_WRONLY, &
          MPI_INFO_NULL,fh,this%ierr_)

    ! Must write part. data to correct spot in file:
    CALL MPI_TYPE_SIZE(MPI_INTEGER,szint)
    CALL MPI_TYPE_SIZE(GC_REAL    ,szreal)
    CALL MPI_ALLREDUCE(this%nparts_,nt,this%nprocs_,MPI_INTEGER,      &
                        MPI_SUM,this%comm_,this%ierr_)
    IF ( this%myrank_ .EQ. 0 ) THEN
        CALL MPI_FILE_WRITE_AT(fh,0,nt,1,MPI_INTEGER,this%istatus_,this%ierr_)
    ENDIF
    IF ( this%nparts_ .GT. 0 ) THEN
      DO j = 1, this%nparts_
        offset = this%id_(j)*szreal+szint
        CALL MPI_FILE_WRITE_AT_ALL(fh,offset,this%ltmp1_(j),1,GC_REAL,this%istatus_,this%ierr_)
      ENDDO
    ELSE
        CALL MPI_FILE_WRITE_AT_ALL(fh,0     ,this%ltmp1_(j),0,GC_REAL,this%istatus_,this%ierr_)
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
!               particle's id.
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
    INTEGER                           :: fh,offset,nt,szreal
    CHARACTER(len=*), INTENT(IN)      :: dir
    CHARACTER(len=*), INTENT(IN)      :: nmb
    CHARACTER(len=*), INTENT(IN)      :: spref
    TYPE(GPDBrec)                     :: pst

    INTEGER                           :: j

!!  IF ( myrank .EQ. 0 ) THEN
!!    this%sfile_ = trim(dir) // '/' // trim(spref) // '.' // &
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

    CALL MPI_TYPE_SIZE(GC_REAL    ,szreal)
    CALL MPI_FILE_OPEN(MPI_COMM_WORLD,trim(dir) // '/' // trim(spref) // &
         '.' // nmb // '.lag',MPI_MODE_CREATE+MPI_MODE_WRONLY, &
          MPI_INFO_NULL,fh,this%ierr_)

    ! Must write part. data to correct position in file:
    CALL MPI_ALLREDUCE(this%nparts_,nt,this%nprocs_,MPI_INTEGER,      &
                        MPI_SUM,this%comm_,this%ierr_)
    IF ( this%myrank_ .EQ. 0 ) THEN
        CALL MPI_FILE_WRITE_AT(fh,     0,real(nt,kind=GP),1,GC_REAL,this%istatus_,this%ierr_)
        CALL MPI_FILE_WRITE_AT(fh,szreal,real(nt,kind=GP),1,GC_REAL,this%istatus_,this%ierr_)
    ENDIF
    IF ( this%nparts_ .GT. 0 ) THEN
      DO j = 1, this%nparts_
        offset  = this%id_(j)*this%sziotype_+2*szreal
        prec(1) = this%px_(j)
        prec(2) = this%py_(j)
        prec(3) = this%pz_(j)
        CALL MPI_FILE_WRITE_AT_ALL(fh,offset,prec,1,this%iotype_,this%ierr_)
      ENDDO
    ELSE
        CALL MPI_FILE_WRITE_AT_ALL(fh,0,prec,0,this%iotype_,this%ierr_)
    ENDIF
    CALL MPI_FILE_CLOSE(fh,this%ierr_)



  END SUBROUTINE GPart_io_write_pdb
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GPart_io_read(this, iunit, dir, spref, nmb)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : io_read
!  DESCRIPTION: Does read of Lagrangian particle data from file
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
    CLASS(GPart)                              :: this
!!  REAL(KIND=GP),INTENT(OUT),DIMENSION(:)    :: lvar
    REAL(KIND=GP)                             :: rvar
    INTEGER,INTENT(IN)                        :: iunit
    INTEGER                                   :: fh,j,nt,szreal
    CHARACTER(len=*), INTENT(IN)              :: dir
    CHARACTER(len=*), INTENT(IN)              :: nmb
    CHARACTER(len=*), INTENT(IN)              :: spref


    CALL MPI_TYPE_SIZE(GC_REAL,szreal)
    CALL MPI_FILE_OPEN(MPI_COMM_WORLD,trim(dir) // '/' // trim(spref) // &
         '.' // nmb // '.lag',MPI_MODE_CREATE+MPI_MODE_WRONLY, &
          MPI_INFO_NULL,fh,this%ierr_)

    ! Must read part. data from correct spot in file:
    CALL MPI_FILE_READ_AT_ALL(fh,0,rvar,1,GC_REAL,this%istatus_,this%ierr_)    !  no.parts
    nt = int(rvar)
    IF ( nt.GT.this%maxparts_ ) THEN
      WRITE(*,*) 'GPart_io_read: Attempt to read too many particles'
      STOP
    ENDIF
    CALL MPI_FILE_READ_AT_ALL(fh,  szreal,rvar,1,GC_REAL,this%istatus_,this%ierr_) ! time
    CALL MPI_FILE_READ_AT_ALL(fh,2*szreal,this%ptmp0_,nt,1,this%iotype_,this%istatus_,this%ierr_)
    CALL MPI_FILE_CLOSE(fh,this%ierr_)
    
    ! From global temp 'd.b.' of particle positions, get local particles to work on:
    CALL GPart_GetLocalWrk(this,this%id_,this%px_,this%py_,this%pz_,this%nparts_,this%ptmp0_,nt)

    ! If there is a global VDB for data 'exchanges', create it here:
    IF ( this%iexchtype_.EQ.GPEXCHTYPE_VDB ) THEN
      this%nvdb_ = nt
      CALL this%exchop_%VDBSynch(this%vdb_,this%maxparts_,this%id_, &
           this%px_,this%py_,this%pz_,this%nparts_,this%ptmp1_)
    ENDIF

!!  CALL GPart_Delete(this,this%pdb_,nt,this%nparts_)

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

    ! Initialize solution, u: 
    ! u* <-- u: 
 
    ! Cycle over JST loop to update state:
    DO i = 1, this%nparts_
       this%ptmp0_(1,i) = this%px_(i)  ! u_0
       this%ptmp0_(2,i) = this%py_(i)  ! u_0
       this%ptmp0_(3,i) = this%pz_(i)  ! u_0
    ENDDO 
    

  END SUBROUTINE GPart_SetStepRKK
!-----------------------------------------------------------------
!-----------------------------------------------------------------


  SUBROUTINE GPart_FinalizeRKK(this)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : FinalizeRKK
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

    ! u(t+dt) = u*: done already
    
    RETURN

  END SUBROUTINE GPart_FinalizeRKK
!-----------------------------------------------------------------
!-----------------------------------------------------------------


  SUBROUTINE GPart_StepRKK(this, vx, vy, vz, dt, xk, tmp1, tmp2)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : StepRKK
!  DESCRIPTION: Carries out one stage of explicit RK-like time
!               inegration step.  Intended for explicit step within 
!               an outer stepper method of the form:
!
!               X = X_0 + dt * V(X(t),t) * xk,

!  ARGUMENTS  :
!    this    : 'this' class instance
!    vz,vy,vz: compoments of velocity field, in real space, partially
!              updated, possibly. These will be overwritten!
!    dt      : integration timestep
!    xk      : multiplicative RK time stage factor
!-----------------------------------------------------------------
    USE grid
    USE fprecision
    USE commtypes
    USE mpivars

    IMPLICIT NONE
    CLASS(GPart)                                         :: this
    INTEGER                                              :: i,j
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(n,n,ksta:kend) :: vx,vy,vz,tmp1,tmp2
    REAL(KIND=GP),INTENT   (IN)                          :: dt,xk
    REAL(KIND=GP)                                        :: dtfact


    dtfact = dt*xk

    ! Find F(u*):
    ! ... x:
    CALL GPart_EulerToLag(this,this%ltmp1_,this%nparts_,vx,.true.,tmp1,tmp2)
    ! ux* <-- ux + dt * F(U*)*xk:
    DO j = 1, this%nparts_
      this%px_(j) = this%ptmp0_(1,j) + dtfact*this%ltmp1_(j)
    ENDDO

    ! ... y:
    ! Exchange bdy data for velocities, so that we
    ! can perform local interpolations:
    CALL GPart_EulerToLag(this,this%ltmp1_,this%nparts_,vy,.false.,tmp1,tmp2)
    ! uy* <-- uy + dt * F(U*)*xk:
    DO j = 1, this%nparts_
      this%py_(j) = this%ptmp0_(2,j) + dtfact*this%ltmp1_(j)
    ENDDO

    ! ... z:
    ! Exchange bdy data for velocities, so that we
    ! can perform local interpolations:
    CALL GPart_EulerToLag(this,this%ltmp1_,this%nparts_,vz,.false.,tmp1,tmp2)
    ! uz* <-- uz + dt * F(U*)*xk:
    DO j = 1, this%nparts_
      this%pz_(j) = this%ptmp0_(3,j) + dtfact*this%ltmp1_(j)
    ENDDO
!

    ! IF using nearest-neighbor interfcae, do particle exchange 
    ! between nearest-neighbor tasks BEFORE PERIODIZING particle coordinates:
    IF ( this%iexchtype_.EQ.GPEXCHTYPE_NN ) THEN
      CALL this%exchop_%PartExchange(this%id_,this%px_,this%py_,this%pz_, &
           this%nparts_,this%lxbnds_(3,1),this%lxbnds_(3,2))
    ENDIF

    ! Enforce periodicity:
    CALL GPart_MakePeriodicP(this,this%px_,this%py_,this%pz_,this%nparts_)

    ! If using VDB interface, do synch-up, and get of local work:
    IF ( this%iexchtype_.EQ.GPEXCHTYPE_VDB ) THEN
      ! Synch up VDB, if necessary:
      CALL this%exchop_%VDBSynch(this%vdb_,this%maxparts_,this%id_, &
           this%px_,this%py_,this%pz_,this%nparts_,this%ptmp1_)
      ! If using VDB, get local particles to work on:
      CALL GPart_GetLocalWrk(this,this%id_,this%px_,this%py_,this%pz_,&
           this%nparts_,this%vdb_,this%maxparts_)
    ENDIF

  END SUBROUTINE GPart_stepRKK
!-----------------------------------------------------------------
!-----------------------------------------------------------------


  SUBROUTINE GPart_EulerToLag(this,lag,nl,evar,doupdate,tmp1,tmp2)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : EulerToLag
!  DESCRIPTION: Computes the Eulerian to Lagrangian
!               transformation by interpolating Eulerian field
!               evar to current position of Lagrangian paricles in 
!               d.b. Array lag must be large enouch to accommodate 
!               max. no. particles; no checking is done. Note
!               that 'evar' array must have local dimensions 
!               for a real array in GHOST (n X n X (kend-ksta+1)).
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
!    tmp1/2  : temp. arrays for interpolation
!    
!-----------------------------------------------------------------
    USE grid
    USE fprecision
    USE commtypes
    USE mpivars

    IMPLICIT NONE
    CLASS(GPart)                                          :: this
    INTEGER      ,INTENT   (IN)                           :: nl
    LOGICAL      ,INTENT   (IN)                           :: doupdate
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(n,n,ksta:kend)  :: evar,tmp1,tmp2
    REAL(KIND=GP),INTENT(INOUT),DIMENSION            (nl)  :: lag
    INTEGER                                               :: j

    IF ( doupdate ) THEN
      CALL this%intop_%PartUpdate3D(this%px_,this%py_,this%pz_,this%nparts_)
    ENDIF
    CALL this%exchop_%CompSpline3D(evar,tmp1,tmp2)
    CALL this%intop_%DoInterp3D(lag,nl)

  END SUBROUTINE GPart_EulerToLag
!-----------------------------------------------------------------
!-----------------------------------------------------------------


  SUBROUTINE GPart_MakePeriodicP(this,px,py,pz,npdb)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : MakePeriodicP
!  DESCRIPTION: Enforces periodic b.c.'s on particles in pdb
!  ARGUMENTS  :
!    this    : 'this' class instance (IN)
!    px      : particle x loc. d.b.
!    py      : particle x loc. d.b.
!    pz      : particle x loc. d.b.
!    npdb    : no. particles in pdb
!    
!-----------------------------------------------------------------
    USE fprecision
    USE commtypes
    USE mpivars

    IMPLICIT NONE
    CLASS(GPart)                                     :: this
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(:)         :: px,py,pz
    INTEGER,INTENT(IN)                               :: npdb

    INTEGER                                          :: j

    DO j = 1, npdb
      px(j) = amod(px(j)+10.0_GP*this%gext_(1),this%gext_(1))
      py(j) = amod(py(j)+10.0_GP*this%gext_(2),this%gext_(2))
      pz(j) = amod(pz(j)+10.0_GP*this%gext_(3),this%gext_(3))
    ENDDO

  END SUBROUTINE GPart_MakePeriodicP
!-----------------------------------------------------------------
!-----------------------------------------------------------------


  SUBROUTINE GPart_MakePeriodicExt(this,v,nx,ny,kb,ke,nc)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : MakePeriodicExt
!  DESCRIPTION: Enforces periodic b.c.'s on extended field
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
    INTEGER,INTENT(IN)                               :: nc,nx,ny,kb,ke
    REAL(KIND=GP),INTENT(INOUT)                      :: v(nx,ny,kb:ke)

    INTEGER                                          :: i,j,k

    ! Recall: nx, ny are the dimensions _including_ ghost zones:
    !
    ! Periodicity s.t.:
    !   | | [ | | | | ] | |
    !   a b       a b 
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


  SUBROUTINE GPart_Delete(this,id,px,py,pz,npdb,nnew)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : GPart_Delete
!  DESCRIPTION: Removes from PDB NULL particles, concatenates list,
!               and sets new number of particles
!  ARGUMENTS  :
!    this    : 'this' class instance (IN)
!    id      : part ids
!    px,py pz: part. d.b.
!    npdb    : no. parts. in pdb
!    nnew    : no. non-NULL particles (set in GPartComm class)
!    
!-----------------------------------------------------------------
    USE fprecision
    USE commtypes

    IMPLICIT NONE
    CLASS(GPart)                               :: this
    INTEGER      ,INTENT (IN)                  :: npdb
    INTEGER      ,INTENT(INOUT),DIMENSION(npdb):: id
    INTEGER      ,INTENT(OUT)                  :: nnew
    INTEGER                                    :: i,j
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(npdb):: px,py,pz

    nnew = 0
    DO i = 1, npdb
       IF ( this%id_(i) .NE. GPNULL ) nnew = nnew + 1
    ENDDO

    j     = 1
    DO i = 1, nnew
      DO WHILE ( j.LE.npdb .AND. id(j).EQ.GPNULL )
        j = j + 1
      ENDDO
      IF ( j.LE.npdb .AND. j.NE.i ) THEN
        id(i) = id(j); id(j) = GPNULL
        px(i) = px(j)
        py(i) = py(j)
        pz(i) = pz(j)
      ENDIF

    ENDDO

  END SUBROUTINE GPart_Delete
!-----------------------------------------------------------------
!-----------------------------------------------------------------


  SUBROUTINE GPart_GetLocalWrk(this,id,lx,ly,lz,nl,gvdb,ngvdb)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : GPart_GetLocalWrk
!  DESCRIPTION: Removes from PDB NULL particles, concatenates list,
!               and sets new number of particles
!  ARGUMENTS  :
!    this    : 'this' class instance (IN)
!    id      : part ids
!    lx,ly,lz: local part. d.b. vectors
!    nl      : no. parts. in local pdb
!    gvdb    : global VDB containing part. position records. Location
!              gives particle id.
!    ngvdb   : no. records in global VDB
!-----------------------------------------------------------------
    USE fprecision
    USE commtypes

    IMPLICIT NONE
    CLASS(GPart)                                  :: this
    INTEGER      ,INTENT(OUT)                     :: nl
    INTEGER      ,INTENT(OUT),DIMENSION(nl)       :: id
    INTEGER      ,INTENT (IN)                     :: ngvdb
    INTEGER                                       :: i,j
    REAL(KIND=GP),INTENT(OUT),DIMENSION(nl)       :: lx,ly,lz
    REAL(KIND=GP),INTENT (IN),DIMENSION(3,ngvdb)  :: gvdb

    nl = 0
    DO j = 1, ngvdb
      id = GPNULL
      IF ( lz(j).GT.this%lxbnds_(3,1) .OR. lz(j).LT.this%lxbnds_(3,2) ) THEN 
        id  (j) = j-1
        lx  (j) = gvdb(1,j)
        ly  (j) = gvdb(2,j)
        lz  (j) = gvdb(3,j)
        nl = nl + 1
      ENDIF
    ENDDO

  END SUBROUTINE GPart_GetLocalWrk
!-----------------------------------------------------------------
!-----------------------------------------------------------------


! SUBROUTINE GPart_GetParts(this,)
!!-----------------------------------------------------------------
!!-----------------------------------------------------------------
!!  METHOD     : GPart_GetParts
!!  DESCRIPTION: Gets particle d.b.
!!  ARGUMENTS  :
!!    this    : 'this' class instance (IN)
!!    id      : part ids
!!    lx,ly,lz: local part. d.b. vectors
!!    nl      : no. parts. in local pdb
!!    gvdb    : global VDB containing part. position records. Location
!!              gives particle id.
!!    ngvdb   : no. records in global VDB
!!-----------------------------------------------------------------
!    USE fprecision
!    USE commtypes
!
!    IMPLICIT NONE 
!    CLASS(GPart)                                  :: this 
!    INTEGER      ,INTENT(OUT)                     :: nl
!    INTEGER      ,INTENT(OUT),DIMENSION(nl)       :: id
!    INTEGER      ,INTENT (IN)                     :: ngvdb
!    INTEGER                                       :: i,j
!    REAL(KIND=GP),INTENT(OUT),DIMENSION(nl)       :: lx,ly,lz
!    REAL(KIND=GP),INTENT (IN),DIMENSION(3,ngvdb)  :: gvdb 
!
!
!  END SUBROUTINE GPart_GetParts
!-----------------------------------------------------------------
!-----------------------------------------------------------------



END MODULE class_GPart
