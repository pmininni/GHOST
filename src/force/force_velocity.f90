! =====================================================================
! NAME       : force_velocity.f90
! DESCRIPTION: Forcing methods for all VelocityBase solver classes.
!              Each force includes an initial set up, and an update.
!              These forcing methods can be used to initialize and
!              update forcings on other solver classes, using
!              composition of operators, as long as no correlation
!              between the forces is needed.
!
! Forces avaliable:
!   null_fv   : Null forcing of the velocity field
!   tg_fv     : Taylor-Green forcing
!
! DATE       : 01/18/26 (PDM)
! =====================================================================

module force_velocity
  USE forcebase_mod

  IMPLICIT NONE

  ! ================= Initial conditions supported ====================
  type, extends(forceBase) :: null_fv
    contains
      procedure ::   init_GForce => init_nullfv
  end type null_fv
  type, extends(forceBase) :: tg_fv
    contains
      procedure ::   init_GForce =>   init_tgfv
  end type tg_fv

CONTAINS

  ! ===================================================================
  ! Forcing functions
  ! ===================================================================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Null forcing
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_nullfv(this,solver,state)
    use gstate_mod
    use hd_mod
    use fprecision
    use grid
    use mpivars
!$  use threads
    implicit none
    
    class     (null_fv), intent   (in) :: this
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
      error stop "This solver does not support velocity forcing"
    end select
  end subroutine init_nullfv


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Taylor-Green forcing
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   f0  : amplitude of the forcing
  !   kdn : minimum wave number (rounded to next integer)
  !   kup : maximum wave number (rounded to next integer)
  subroutine init_tgfv(this,solver,state)
    use gstate_mod
    use hd_mod
    use fprecision
    use grid
    use boxsize
    use commtypes
    use fft
    use ali
    use var
    use pseudospec_norm
!$  use threads
    implicit none
    
    class       (tg_fv), intent   (in)       :: this
    class(EquationBase), intent   (in)       :: solver
    type       (Gstate), intent(inout)       :: state(:)
    real(kind=GP), pointer, dimension(:,:,:) :: R1,R2
    real(kind=GP)                            :: f0,kdn,kup
    integer                                  :: i,j,k,ki
    logical                                  :: bret

    namelist/ tg_fv / f0,kdn,kup
    CALL solver%workspace_%get_real_tmp(R1,bret)
    CALL solver%workspace_%get_real_tmp(R2,bret)
    select type (solver)
    class is (VelocityBase)
    ! Read parameters from a namelist in the input file
    if ( myrank .eq. 0 ) then
      open(1,file=solver%infile_,status='unknown',form="formatted")
      read(1,NML=tg_fv)
      close(1)
    endif
    call mpi_bcast(f0 ,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(kup,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(kdn,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
    ! Generate a superposition of TG flows
    IF ( abs(Lx-Ly).gt.tiny ) THEN
      IF (myrank.eq.0) error stop "TG initial conditions require Lx=Ly"
    ENDIF
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
    DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
       DO j = 1,ny
          DO k = 1,nz
            state(solver%VELOCITY+2)%ccomp(k,j,i) = 0.0_GP !fz
          END DO
       END DO
    END DO
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
    DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
       DO j = 1,ny
          DO i = 1,nx
          R1(i,j,k) = 0.0_GP
          R2(i,j,k) = 0.0_GP
          DO ki = INT(kdn),INT(kup)
             R1(i,j,k) = R1(i,j,k)+SIN(2*pi*ki*(real(i,kind=GP)-1)/ &
                      real(nx,kind=GP))*COS(2*pi*ki*(real(j,kind=GP)-1)/ &
                      real(ny,kind=GP))*COS(2*pi*ki*(real(k,kind=GP)-1)/ &
                      real(nz,kind=GP))
             R2(i,j,k) = R2(i,j,k)-COS(2*pi*ki*(real(i,kind=GP)-1)/ &
                      real(nx,kind=GP))*SIN(2*pi*ki*(real(j,kind=GP)-1)/ &
                      real(ny,kind=GP))*COS(2*pi*ki*(real(k,kind=GP)-1)/ &
                      real(nz,kind=GP))
          END DO
          END DO
       END DO
    END DO
    CALL fftp3d_real_to_complex(planrc,R1, & !fx
                   state(solver%VELOCITY  )%ccomp,MPI_COMM_WORLD)
    CALL fftp3d_real_to_complex(planrc,R2, & !fy
                   state(solver%VELOCITY+1)%ccomp,MPI_COMM_WORLD)
    CALL normalize(state(solver%VELOCITY  )%ccomp, &
                   state(solver%VELOCITY+1)%ccomp, &
                   state(solver%VELOCITY+2)%ccomp, &
                   f0,1,MPI_COMM_WORLD)
    class default
      error stop "This solver does not support velocity forcing"
    end select
    CALL solver%workspace_%free_real_tmp(R1)
    CALL solver%workspace_%free_real_tmp(R2)
  end subroutine init_tgfv


  ! ===================================================================
  ! Forcing update schemes
  ! ===================================================================
  
end module force_velocity
