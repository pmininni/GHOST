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
  !! Read the velocity field from restart files
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_readpos(this,psolver,pstate)
    use gpstate_mod
    
    use status
    use pstatus
    use filefmt
    use commtypes
    implicit none
    class    (read_pos), intent   (in)          :: this
    class(ParticleBase), intent(inout)          :: psolver
    type  (GPStateComp), intent(inout)          :: pstate(:)
  end subroutine init_readpos


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Read positions from a user defined input file
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_userpos(this,psolver,pstate)
    use gpstate_mod
    use status
    use pstatus
    use filefmt
    use commtypes
    implicit none
    class    (user_pos), intent   (in)          :: this
    class(ParticleBase), intent(inout)          :: psolver
    type  (GPStateComp), intent(inout)          :: pstate(:)
  end subroutine init_userpos


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Taylor-Green flow
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_randompos(this,psolver,pstate)
    use gpstate_mod
    use random
    implicit none
    class    (rand_pos), intent   (in)          :: this
    class(ParticleBase), intent(inout)          :: psolver
    type  (GPStateComp), intent(inout)          :: pstate(:)
    integer                           :: ib,ie,j,iwrk1,iwrk2,nt
    real(kind=GP)                     :: c,r,x1,x2

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
      WRITE(*,*) 'init_randompos: Inconsistent particle count: maxparts=',  &
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
