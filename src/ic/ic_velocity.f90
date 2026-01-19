! =====================================================================
! NAME       : ic_velocity.f90
! DESCRIPTION: Initial conditions for all VelocityBase solver classes.
!              These initial conditions can be used to initialize
!              velocities on other solver classes, using composition of
!              initial condition operators, as long as no correlation
!              between the fields is needed.
!
! Initial conditions avaliable:
!   null_v   : Null velocity field
!   tg_v     : Taylor-Green flow
!
! DATE       : 01/16/26 (PDM)
! =====================================================================

module ic_velocity
  USE icbase_mod

  IMPLICIT NONE

  ! ================= Initial conditions supported ====================
  type, extends(icBase) :: null_v
    contains
      procedure :: init_GState => init_nullv
  end type null_v
  type, extends(icBase) :: tg_v
    contains
      procedure :: init_GState => init_tgv
  end type tg_v

CONTAINS

  ! ===================================================================
  ! Initial conditions
  ! ===================================================================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Null velocity
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_nullv(this,solver,state)
    use gstate_mod
    use hd_mod
    use fprecision
    use grid
    use mpivars
!$  use threads
    implicit none
    
    class      (null_v), intent   (in) :: this
    class(EquationBase), intent   (in) :: solver
    type       (Gstate), intent(inout) :: state(:)
    integer                            :: i,j,k

    select type (solver)
    class is (VelocityBase)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
        DO j = 1,ny
          DO k = 1,nz
            state(solver%VELOCITY  )%ccomp(k,j,i) = 0.0_GP
            state(solver%VELOCITY+1)%ccomp(k,j,i) = 0.0_GP
            state(solver%VELOCITY+2)%ccomp(k,j,i) = 0.0_GP
          END DO
        END DO
      END DO
    class default
      error stop "This solver does not support velocity field ICs"
    end select
  end subroutine init_nullv


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Taylor-Green flow
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   u0  : amplitude of each Taylor-Green mode
  !   kdn : minimum wave number (rounded to next integer)
  !   kup : maximum wave number (rounded to next integer)
  subroutine init_tgv(this,solver,state)
    use gstate_mod
    use hd_mod
    use fprecision
    use grid
    use boxsize
    use commtypes
    use fftplans
!$  use threads
    implicit none
    
    class        (tg_v), intent   (in)       :: this
    class(EquationBase), intent   (in)       :: solver
    type       (Gstate), intent(inout)       :: state(:)
    real(kind=GP), pointer, dimension(:,:,:) :: R1,R2
    real(kind=GP)                            :: u0,kdn,kup
    integer                                  :: i,j,k,ki
    logical                                  :: bret

    namelist/ tg_v / u0,kdn,kup
    select type (solver)
    class is (VelocityBase)
    ! Read parameters from a namelist in the input file
    if ( myrank .eq. 0 ) then
      open(1,file=solver%infile_,status='unknown',form="formatted")
      read(1,NML=tg_v)
      close(1)
    endif
    call mpi_bcast(u0 ,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(kup,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(kdn,1,GC_REAL,0,MPI_COMM_WORLD,ierr)

    ! Generate a superposition of TG flows
!!$    IF ( abs(Lx-Ly).gt.tiny ) THEN
!!$      IF (myrank.eq.0) error "TG initial conditions require Lx=Ly"
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
!!$    CALL solver%workspace_%get_real_tmp(R1,bret)
!!$    CALL solver%workspace_%get_real_tmp(R2,bret)
!!$!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
!!$    DO k = ksta,kend
!!$!$omp parallel do if (kend-ksta.lt.nth) private (i)
!!$       DO j = 1,ny
!!$          DO i = 1,nx
!!$
!!$          R1(i,j,k) = 0.0_GP
!!$          R2(i,j,k) = 0.0_GP
!!$
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
!!$
!!$          END DO
!!$       END DO
!!$    END DO
!!$    CALL fftp3d_real_to_complex(planrc,R1, & !vx
!!$                    state(solver%VELOCITY  )%ccomp,MPI_COMM_WORLD)
!!$    CALL fftp3d_real_to_complex(planrc,R2, & !vy
!!$                    state(solver%VELOCITY+1)%ccomp,MPI_COMM_WORLD)
!!$    CALL solver%workspace_%free_real_tmp(R1)
!!$    CALL solver%workspace_%free_real_tmp(R2)
    class default
      error stop "This solver does not support velocity field ICs"
    end select
  end subroutine init_tgv

end module ic_velocity
