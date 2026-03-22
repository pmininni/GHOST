! =====================================================================
! NAME       : icp_position.f90
! DESCRIPTION: Initial conditions for the position of all ParticleBase
!              solver classes. These initial conditions can be used to
!              initialize positions of particles in other solver
!              classes, using composition of initial condition
!              operators.
!
! Initial conditions avaliable:
!   read_position  : Reads positions from an input file numbered by stat
!   user_position  : Reads positions from a user defined input file
!   random_position: Random initial positions
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
    class(ParticleBase), intent   (in)          :: psolver
    type  (GPStateComp), intent(inout)          :: pstate(:)
!!$    real   (kind=GP), pointer, dimension(:,:,:) :: R1
!!$    integer                          :: i
!!$    logical                          :: bret
!!$
!!$    if ((stat .eq. 0).and.(solver%myrank_ .eq. 0)) then
!!$      stop 'Cannot read files if starting a new run with stat=0'
!!$    endif
!!$    call solver%workspace_%get_real_tmp(R1,bret)
!!$    select type (solver)
!!$    class is (VelocityBase)
!!$      tind = int(stat)
!!$      WRITE(ext, fmtext) tind
!!$      do i = solver%VELOCITY,solver%VELOCITY+solver%nc_-1
!!$        call io_read(1,idir,trim(solver%sstate_(i)),ext,solver%planio_,R1)
!!$        call fftp3d_real_to_complex(planrc,R1,state(i)%ccomp,MPI_COMM_WORLD)
!!$      end do
!!$    class default
!!$      error stop "This solver does not support velocity field ICs"
!!$    end select
!!$    call solver%workspace_%free_real_tmp(R1)
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
    class(ParticleBase), intent   (in)          :: psolver
    type  (GPStateComp), intent(inout)          :: pstate(:)
!!$    select type (solver)
!!$    class is (VelocityBase)
!!$!$omp parallel do if (iend-ista.ge.nth) private (j,k)
!!$      DO i = ista,iend
!!$!$omp parallel do if (iend-ista.lt.nth) private (k)
!!$        DO j = 1,ny
!!$          DO k = 1,nz
!!$            state(solver%VELOCITY  )%ccomp(k,j,i) = 0.0_GP
!!$            state(solver%VELOCITY+1)%ccomp(k,j,i) = 0.0_GP
!!$            state(solver%VELOCITY+2)%ccomp(k,j,i) = 0.0_GP
!!$          END DO
!!$        END DO
!!$      END DO
!!$    class default
!!$      error stop "This solver does not support velocity field ICs"
!!$    end select
  end subroutine init_userpos


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Taylor-Green flow
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_randompos(this,psolver,pstate)
    use gpstate_mod
    use status
    use pstatus
    use filefmt
    use commtypes
    implicit none
    class    (rand_pos), intent   (in)          :: this
    class(ParticleBase), intent   (in)          :: psolver
    type  (GPStateComp), intent(inout)          :: pstate(:)
!!$    namelist/ tg_v / u0,kdn,kup
!!$    CALL solver%workspace_%get_real_tmp(R1,bret)
!!$    CALL solver%workspace_%get_real_tmp(R2,bret)
!!$    select type (solver)
!!$    class is (VelocityBase)
!!$    ! Read parameters from a namelist in the input file
!!$    if ( myrank .eq. 0 ) then
!!$      open(1,file=solver%infile_,status='unknown',form="formatted")
!!$      read(1,NML=tg_v)
!!$      close(1)
!!$    endif
!!$    call mpi_bcast(u0 ,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
!!$    call mpi_bcast(kup,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
!!$    call mpi_bcast(kdn,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
!!$    ! Generate a superposition of TG flows
!!$    IF ( abs(Lx-Ly).gt.tiny ) THEN
!!$      IF (myrank.eq.0) error stop "TG initial conditions require Lx=Ly"
!!$    ENDIF
!!$!$omp parallel do if (iend-ista.ge.nth) private (j,k)
!!$    DO i = ista,iend
!!$!$omp parallel do if (iend-ista.lt.nth) private (k)
!!$       DO j = 1,ny
!!$          DO k = 1,nz
!!$            state(solver%VELOCITY+2)%ccomp(k,j,i) = 0.0_GP !vz
!!$          END DO
!!$       END DO
!!$    END DO
!!$!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
!!$    DO k = ksta,kend
!!$!$omp parallel do if (kend-ksta.lt.nth) private (i)
!!$       DO j = 1,ny
!!$          DO i = 1,nx
!!$          R1(i,j,k) = 0.0_GP
!!$          R2(i,j,k) = 0.0_GP
!!$          DO ki = INT(kdn),INT(kup)
!!$             R1(i,j,k) = R1(i,j,k)+SIN(2*pi*ki*(real(i,kind=GP)-1)/ &
!!$                      real(nx,kind=GP))*COS(2*pi*ki*(real(j,kind=GP)-1)/ &
!!$                      real(ny,kind=GP))*COS(2*pi*ki*(real(k,kind=GP)-1)/ &
!!$                      real(nz,kind=GP))
!!$             R2(i,j,k) = R2(i,j,k)-COS(2*pi*ki*(real(i,kind=GP)-1)/ &
!!$                      real(nx,kind=GP))*SIN(2*pi*ki*(real(j,kind=GP)-1)/ &
!!$                      real(ny,kind=GP))*COS(2*pi*ki*(real(k,kind=GP)-1)/ &
!!$                      real(nz,kind=GP))
!!$          END DO
!!$          END DO
!!$       END DO
!!$    END DO
!!$    CALL fftp3d_real_to_complex(planrc,R1, & !vx
!!$                   state(solver%VELOCITY  )%ccomp,MPI_COMM_WORLD)
!!$    CALL fftp3d_real_to_complex(planrc,R2, & !vy
!!$                   state(solver%VELOCITY+1)%ccomp,MPI_COMM_WORLD)
!!$    CALL normalize(state(solver%VELOCITY  )%ccomp, &
!!$                   state(solver%VELOCITY+1)%ccomp, &
!!$                   state(solver%VELOCITY+2)%ccomp, &
!!$                   u0,1,MPI_COMM_WORLD)
!!$    class default
!!$      error stop "This solver does not support velocity field ICs"
!!$    end select
!!$    CALL solver%workspace_%free_real_tmp(R1)
!!$    CALL solver%workspace_%free_real_tmp(R2)
  end subroutine init_randompos

end module icp_position
