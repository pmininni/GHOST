! =====================================================================
! NAME       : icp_position.f90
! DESCRIPTION: Initial conditions for the position of all ParticleBase
!              solver classes. These initial conditions can be used to
!              initialize positions of particles in other solver
!              classes, using composition of initial condition
!              operators.
!
! Initial conditions avaliable:
!   read_x   : Reads positions from an input file numbered by stat
!   user_x   : Reads positions from a user defined input file
!   random_x : Random initial positions
!
! DATE       : 03/21/26 (PDM)
! =====================================================================

module icp_position
  use icpbase_mod
  use particlebase_mod

  IMPLICIT NONE

  ! ================= Initial conditions supported ====================
  type, extends(icpBase) :: read_pos
    contains
      procedure :: init_GPState => init_readpos
  end type read_pos 
  type, extends(icpBase) :: user_pos
    contains
      procedure :: init_GPState => init_userpos
  end type user_pos 
  type, extends(icpBase) :: rand_pos
    contains
      procedure :: init_GPState => init_randompos
   end type rand_pos

CONTAINS

  ! ===================================================================
  ! Initial conditions
  ! ===================================================================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Read particle positions from restart files
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_readpos(this,psolver,pstate)
    use gpstate_mod
    use status
    use pstatus
    use filefmt
    use commtypes
    implicit none
    class    (read_pos), intent   (in)   :: this
    class(ParticleBase), intent(inout)   :: psolver
    type  (GPStateComp), intent(inout)   :: pstate(:)

    if ((stat .eq. 0).and.(psolver%myrank_ .eq. 0)) then
      stop 'Cannot read restart files if starting a new run with stat=0'
    endif
    pind = int((stat-1)*lgmult+1)
    WRITE(lgext, lgfmtext) pind
    CALL io_read(psolver,pstate,1,psolver%idir_,'xlg',lgext)
  end subroutine init_readpos


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Read positions from a user defined input file. Pstates
  !! outside the one used for reading must be resized to the
  !! size of arrays in pstate after calling this routine.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_userpos(this,psolver,pstate)
    use gpstate_mod
    use status
    use pstatus
    use filefmt
    use commtypes
    implicit none
    class    (user_pos), intent   (in)   :: this
    class(ParticleBase), intent(inout)   :: psolver
    type  (GPStateComp), intent(inout)   :: pstate(:)
    integer                              :: nl,nt
    real(kind=GP)                        :: x,y,z

    ! Note: each record (line) consists of x y z real positions
    ! within [0,NX-1]x[0,NY-1]x[0,NZ-1] box or the equivalent
    ! in box units depending on wrtunit_ class options.
    OPEN(UNIT=5,FILE=trim(psolver%seedfile_),STATUS='OLD',ACTION='READ',   &
         IOSTAT=psolver%ierr_, IOMSG=psolver%serr_);
    IF ( psolver%ierr_ .NE. 0 ) THEN
      WRITE(*,*)'Init_userpos: file:',trim(psolver%seedfile_),' err: ',    &
         trim(psolver%serr_) 
      STOP
    ENDIF
    READ(5,*,IOSTAT=psolver%ierr_) nt
    IF ( psolver%myrank_.eq.0 .AND. nt.NE.psolver%maxparts_ ) THEN
      WRITE(*,*) 'Init_userpos: Inconsistent seed file: required no. part.=', &
      psolver%maxparts_,' total listed: ',nt,' file:',psolver%seedfile_
      STOP
    ENDIF
    READ(5,*,IOSTAT=psolver%ierr_) x

    nt = 0  ! global part. record counter
    nl = 0  ! local particle counter
    DO WHILE ( psolver%ierr_.EQ.0 )
      READ(5,*,IOSTAT=psolver%ierr_) x, y, z
      IF ( psolver%ierr_ .NE. 0 ) THEN
        WRITE(*,*) 'Init_userpos: terminating read; nt=', nt
        EXIT
      ENDIF
      IF ( psolver%wrtunit_ .EQ. 1 ) THEN ! rescale coordinates to grid units
        x = x*psolver%invdel_(1)
        y = y*psolver%invdel_(2)
        z = z*psolver%invdel_(3)
      ENDIF
      IF ( z.GE.psolver%lxbnds_(3,1) .AND. z.LT.psolver%lxbnds_(3,2) .AND. &
           y.GE.psolver%lxbnds_(2,1) .AND. y.LT.psolver%lxbnds_(2,2) .AND. &
           x.GE.psolver%lxbnds_(1,1) .AND. x.LT.psolver%lxbnds_(1,2) ) THEN
        nl = nl + 1
        IF (nl.GT.psolver%partbuff_) THEN
          psolver%partbuff_ = psolver%partbuff_ + psolver%partchunksize_
          CALL ResizeArrays  (psolver,psolver%partbuff_,.true.)
          call GPState_resize(pstate ,psolver%partbuff_)
        END IF
        psolver%id_(nl)                      = nt
        pstate(psolver%POSITION  )%rcomp(nl) = x
        pstate(psolver%POSITION+1)%rcomp(nl) = y
        pstate(psolver%POSITION+2)%rcomp(nl) = z
      ENDIF
      nt = nt + 1
    ENDDO
    CLOSE(5)

    psolver%nparts_ = nl;
    IF ((psolver%myrank_.NE.0).OR.(psolver%bcollective_.EQ.1)) THEN
      IF (psolver%nparts_.LT.(psolver%partbuff_-psolver%partchunksize_)) THEN
        psolver%partbuff_ = psolver%partbuff_ - ((psolver%partbuff_-psolver%nparts_) &
                        /psolver%partchunksize_-1)*psolver%partchunksize_
        CALL ResizeArrays  (psolver,psolver%partbuff_,.false.) 
        call GPState_resize(pstate ,psolver%partbuff_)
      END IF
    END IF
    CALL MPI_ALLREDUCE(nl,nt,1,MPI_INTEGER,MPI_SUM,psolver%comm_,psolver%ierr_)
    IF ( psolver%myrank_.eq.0 .AND. nt.NE.psolver%maxparts_ ) THEN
      WRITE(*,*) 'Init_userpos: Inconsistent particle count: required no.=',    &
      psolver%maxparts_,' total read: ',nt,' file:',psolver%seedfile_
      STOP
    ENDIF
    IF (psolver%iexchtype_.EQ.GPEXCHTYPE_VDB) THEN
      CALL psolver%gpcomm_%VDBSynch(psolver%vdb_,psolver%maxparts_,psolver%id_, &
           pstate(psolver%POSITION  )%rcomp,                                    &
           pstate(psolver%POSITION+1)%rcomp,                                    &
           pstate(psolver%POSITION+2)%rcomp,                                    &
           psolver%nparts_,psolver%ptmp0_)
    END IF
  end subroutine init_userpos


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Initialize particle positions randomly
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_randompos(this,psolver,pstate)
    use gpstate_mod
    use random
    implicit none
    class    (rand_pos), intent   (in)   :: this
    class(ParticleBase), intent(inout)   :: psolver
    type  (GPStateComp), intent(inout)   :: pstate(:)
    integer                              :: ib,ie,j,iwrk1,iwrk2,nt
    real(kind=GP)                        :: c,r,x1,x2

    ! Note: Box is [0,N-1]^3.
    ! All tasks write to _global_ grid, expecting a VDBSynch afterwards
    ! and a GetLocalWrk call to get the particles a task owns when using VDB
    iwrk1 = psolver%maxparts_/psolver%nprocs_
    iwrk2 = modulo(psolver%maxparts_,psolver%nprocs_)
    IF ( iwrk1.GT.0 ) THEN
      ib = psolver%myrank_*iwrk1+1+min(psolver%myrank_,iwrk2)
      ie = ib+iwrk1-1
      IF ( iwrk2.GT.psolver%myrank_ ) ie = ie + 1
    ELSE
      ib = psolver%myrank_+1
      ie = ib
      IF ( psolver%myrank_+1.GT.psolver%maxparts_ ) THEN
        ie = 0
        ib = 1
      ENDIF
    ENDIF
    ib = ib - 1
    ie = ie - 1
    psolver%nparts_ = ie - ib + 1
    DO j = 1, psolver%nparts_
      psolver%id_(j)    = ib + j - 1
      CALL prandom_number(r)
      x1 = real(psolver%libnds_(1,1)-1,kind=GP)
      x2 = real(psolver%libnds_(1,2),kind=GP)
      c = r*(psolver%nd_(1))
      pstate(psolver%POSITION  )%rcomp(j) = min(max(c,x1),x2)
      CALL prandom_number(r)
      x1 = real(psolver%libnds_(2,1)-1,kind=GP)
      x2 = real(psolver%libnds_(2,2),kind=GP)
      c = r*(psolver%nd_(2))
      pstate(psolver%POSITION+1)%rcomp(j) = min(max(c,x1),x2)
      CALL prandom_number(r)
      x1 = real(psolver%libnds_(3,1)-1,kind=GP)
      x2 = real(psolver%libnds_(3,2),kind=GP);
      pstate(psolver%POSITION+2)%rcomp(j) = min(max(x1+r*(x2-x1),x1),x2)
    ENDDO
    CALL MPI_ALLREDUCE(psolver%nparts_,nt,1,MPI_INTEGER,MPI_SUM,psolver%comm_,  &
                       psolver%ierr_)
    IF ( psolver%myrank_.eq.0 .AND. nt.NE.psolver%maxparts_ ) THEN
      WRITE(*,*) 'init_randompos: Inconsistent particle count: maxparts=',      &
      psolver%maxparts_,' total created: ',nt
      STOP
    ENDIF
    IF ( psolver%iexchtype_.EQ.GPEXCHTYPE_VDB ) THEN
      CALL psolver%gpcomm_%VDBSynch(psolver%vdb_,psolver%maxparts_,psolver%id_, &
                                    pstate(psolver%POSITION  )%rcomp,           &
                                    pstate(psolver%POSITION+1)%rcomp,           &
                                    pstate(psolver%POSITION+2)%rcomp,           &
                                    psolver%nparts_,psolver%ptmp0_)
      CALL GetLocalWrk(psolver,psolver%id_,                                     &
                       pstate(psolver%POSITION  )%rcomp,                        &
                       pstate(psolver%POSITION+1)%rcomp,                        &
                       pstate(psolver%POSITION+2)%rcomp,                        &
                       psolver%nparts_,psolver%vdb_,psolver%maxparts_)  
      IF ( psolver%wrtunit_ .EQ. 1 ) THEN ! rescale coordinates to box units
        psolver%ptmp0_(1,:) = psolver%vdb_(1,:)*psolver%delta_(1)
        psolver%ptmp0_(2,:) = psolver%vdb_(2,:)*psolver%delta_(2)
        psolver%ptmp0_(3,:) = psolver%vdb_(3,:)*psolver%delta_(3)
        CALL ascii_write_lag(psolver,1,'.','xlgInitRndSeed','000',0.0_GP,       &
                             psolver%maxparts_,psolver%ptmp0_(1,:),             &
                             psolver%ptmp0_(2,:),psolver%ptmp0_(3,:))
      ELSE
        CALL ascii_write_lag(psolver,1,'.','xlgInitRndSeed','000',0.0_GP,       &
                             psolver%maxparts_,psolver%vdb_(1,:),               &
                             psolver%vdb_(2,:),psolver%vdb_(3,:))
      ENDIF
    ENDIF
    IF ( .NOT. PartNumConsistent(psolver,psolver%nparts_) ) THEN
      IF ( psolver%myrank_.eq.0 ) THEN
        WRITE(*,*) 'init_randompos: Invalid particle after GetLocalWrk call'
        STOP
      ENDIF
    ENDIF
  end subroutine init_randompos

end module icp_position
