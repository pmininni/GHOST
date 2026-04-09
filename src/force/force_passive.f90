! =====================================================================
! NAME       : force_passive.f90
! DESCRIPTION: Forcing methods for passive scalars in all VelocityBase
!              solver classes. Each force needs an initial set up, and
!              an update.
!
! Forces avaliable:
!   null_fps    : Null forcing of the passive scalar
!   puff_fps    : Localized forcing of the passive scalar
!   random_fps  : Random forcing
!
! Update methods available:
!   constant_fps: Constant forcing (no method)
!   shift_fps   : Instantaneous phase shift
!
! DATE       : 01/27/26 (PDM)
! =====================================================================

module force_passive
  USE forcebase_mod

  IMPLICIT NONE

  ! ================= Forcing functions supported =====================
  type, extends(forceBase) :: null_fps
    contains
      procedure :: init_GForce => init_nullfps
  end type null_fps
  type, extends(forceBase) :: puff_fps
    contains
      procedure :: init_GForce => init_pufffps
  end type puff_fps
  type, extends(forceBase) :: random_fps
    contains
      procedure :: init_GForce => init_randomfps
  end type random_fps
! type, extends(forceBase) :: userdef_fps
!   contains
!     procedure :: init_GForce => init_userdeffps
! end type userdef_fps

  ! ================= Update methods supported =======================
  type, extends(forceUpdt) :: shiftupdt_fps
    contains
      procedure :: update_GForce => update_shiftfps
  end type shiftupdt_fps
! type, extends(forceUpdt) :: userupdt_fps
!   contains
!     procedure :: init_GForce => init_userupdtfps
! end type userupdt_fps

CONTAINS

  ! ===================================================================
  ! Forcing functions
  ! ===================================================================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Null forcing
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_nullfps(this,solver,state)
    use gstate_mod
    use hd_mod
    use grid
    use mpivars
!$  use threads
    implicit none
    
    class     (null_fps), intent   (in) :: this
    class(EquationBase), intent   (in) :: solver
    type   (GStateComp), intent(inout) :: state(:)
    integer                            :: i,j,k,n

    select type (solver)
    class is (VelocityBase)
      if ( solver%numpassive_ .eq. 0) then
        stop 'Force: Asking for passive scalar forcing with npassive = 0'
      endif
      do n = solver%PASSIVE, solver%PASSIVE+solver%numpassive_-1    
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
        DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
          DO j = 1,ny
            DO k = 1,nz
              state(n)%ccomp(k,j,i) = 0.0_GP
            END DO
          END DO
        END DO
      end do
    class default
      stop "Force: This solver does not support passive scalars"
    end select
  end subroutine init_nullfps


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Forcing in a Gaussian ball ('puff') centered at
  !! (x0, y0, z0) with a FWHM of r0
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   f0        : vector with the amplitudes of the forcing
  !   x0, y0, z0: vectors with the centers of the puffps
  !   r0        : vector with the radii of the puffps
  subroutine init_pufffps(this,solver,state)
    use gstate_mod
    use hd_mod
    use grid
    use mpivars
    use commtypes
    use fft
    use pseudospec_scalar
!$  use threads
    implicit none
    
    class     (puff_fps), intent   (in) :: this
    class(EquationBase), intent   (in) :: solver
    type   (GStateComp), intent(inout) :: state(:)
    real      (kind=GP), pointer       :: R1(:,:,:)
    real      (kind=GP), allocatable, dimension(:)  :: f0,x0,y0,z0,r0
    double precision                   :: tmp
    integer                            :: i,j,k,n
    logical                            :: bret

    namelist/ puff_fps / f0,x0,y0,z0,r0
    CALL solver%workspace_%get_real_tmp(R1,bret)    
    select type (solver)
    class is (VelocityBase)
    if ( solver%numpassive_ .eq. 0) then
      stop 'Force: Asking for passive scalar forcing with npassive = 0'
    endif
    allocate ( f0(solver%numpassive_) )
    allocate ( x0(solver%numpassive_) )
    allocate ( y0(solver%numpassive_) )
    allocate ( z0(solver%numpassive_) )
    allocate ( r0(solver%numpassive_) )
    if ( myrank .eq. 0 ) then
      open(1,file=solver%infile_,status='unknown',form="formatted")
      read(1,NML=puff_fps)
      close(1)
    endif
    call mpi_bcast(f0,solver%numpassive_,GC_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(x0,solver%numpassive_,GC_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(y0,solver%numpassive_,GC_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(z0,solver%numpassive_,GC_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(r0,solver%numpassive_,GC_REAL,0,MPI_COMM_WORLD,ierr)
    do n = solver%PASSIVE, solver%PASSIVE+solver%numpassive_-1    
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      do k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
        do j = 1,ny
          do i = 1,nx
            tmp = (real(i-1,kind=GP)/real(nx-1,kind=GP)-x0(n-solver%PASSIVE+1))**2 &
                + (real(j-1,kind=GP)/real(ny-1,kind=GP)-y0(n-solver%PASSIVE+1))**2 &
                + (real(k-1,kind=GP)/real(nz-1,kind=GP)-z0(n-solver%PASSIVE+1))**2 
            R1(i,j,k) = exp(-tmp**2/r0(n-solver%PASSIVE+1)**2)
          end do
        end do
      end do
      call fftp3d_real_to_complex(planrc,R1,state(n)%ccomp,MPI_COMM_WORLD)
      call variance(state(n)%ccomp,tmp,1)
      call mpi_bcast(tmp,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      do i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
        do j = 1,ny
          do k = 1,nz
            state(n)%ccomp(k,j,i) = state(n)%ccomp(k,j,i)* &
                         f0(n-solver%PASSIVE+1)/sqrt(tmp)
          end do
        end do
      end do
    end do
    class default
      stop "Force: This solver does not support passive scalars"
    end select
  end subroutine init_pufffps

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Random forcing
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   f0 : vector with the amplitudes of the forcing
  !   kdn: vector with minimum forcing wave number
  !   kup: vector with maximum forcing wave number
  subroutine init_randomfps(this,solver,state)
    use gstate_mod
    use hd_mod
    use random
    use grid
    use mpivars
    use commtypes
    use status
    use fft
    use var
    use kes
    use pseudospec_scalar
!$  use threads
    implicit none
    
    class   (random_fps), intent   (in) :: this
    class(EquationBase), intent   (in) :: solver
    type   (GStateComp), intent(inout) :: state(:)
    real      (kind=GP), allocatable, dimension(:)  :: f0,kdn,kup
    real      (kind=GP)                             :: skup,skdn
    real      (kind=GP)                             :: dump,phase
    double precision                   :: tmp
    integer                            :: i,j,k,n
    logical                            :: bret

    namelist/ random_fps / f0,kup,kdn
    select type (solver)
    class is (VelocityBase)
    if ( solver%numpassive_ .eq. 0) then
      stop 'Force: Asking for passive scalar forcing with npassive = 0'
    endif
    allocate ( f0 (solver%numpassive_) )
    allocate ( kdn(solver%numpassive_) )
    allocate ( kup(solver%numpassive_) )
    if ( myrank .eq. 0 ) then
      open(1,file=solver%infile_,status='unknown',form="formatted")
      read(1,NML=random_fps)
      close(1)
    endif
    call mpi_bcast(f0 ,solver%numpassive_,GC_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(kdn,solver%numpassive_,GC_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(kup,solver%numpassive_,GC_REAL,0,MPI_COMM_WORLD,ierr)
    do n = solver%PASSIVE, solver%PASSIVE+solver%numpassive_-1
      skdn = kdn(n-solver%PASSIVE+1)
      skup = kup(n-solver%PASSIVE+1)
      IF (ista.eq.1) THEN
        state(n)%ccomp(1,1,1) = 0.
        DO j = 2,ny/2+1
          IF ((kk2(1,j,1).le.skup**2).and.(kk2(1,j,1).ge.skdn**2)) THEN
            dump = 1./sqrt(kk2(1,j,1))
            phase = 2*pi*randu(seed)
            state(n)%ccomp(1,j,1) = (COS(phase)+im*SIN(phase))*dump
            state(n)%ccomp(1,ny-j+2,1) = conjg(state(n)%ccomp(1,j,1))
          ELSE
            state(n)%ccomp(1,j,1) = 0.
            state(n)%ccomp(1,ny-j+2,1) = 0.
          ENDIF
        END DO
        DO k = 2,nz/2+1
          IF ((kk2(k,1,1).le.skup**2).and.(kk2(k,1,1).ge.skdn**2)) THEN
            dump = 1./sqrt(kk2(k,1,1))
            phase = 2*pi*randu(seed)
            state(n)%ccomp(k,1,1) = (COS(phase)+im*SIN(phase))*dump
            state(n)%ccomp(nz-k+2,1,1) = conjg(state(n)%ccomp(k,1,1))
          ELSE
            state(n)%ccomp(k,1,1) = 0.
            state(n)%ccomp(nz-k+2,1,1) = 0.
          ENDIF
        END DO
        DO j = 2,ny
          DO k = 2,nz/2+1
            IF ((kk2(k,j,1).le.skup**2).and.(kk2(k,j,1).ge.skdn**2)) THEN
              dump = 1./sqrt(kk2(k,j,1))
              phase = 2*pi*randu(seed)
              state(n)%ccomp(k,j,1) = (COS(phase)+im*SIN(phase))*dump
              state(n)%ccomp(nz-k+2,ny-j+2,1) = conjg(state(n)%ccomp(k,j,1))
            ELSE
              state(n)%ccomp(k,j,1) = 0.
              state(n)%ccomp(nz-k+2,ny-j+2,1) = 0.
            ENDIF
          END DO
        END DO
        DO i = 2,iend
          DO j = 1,ny
            DO k = 1,nz
              IF ((kk2(k,j,i).le.skup**2).and.(kk2(k,j,i).ge.skdn**2)) THEN
                dump = 1./sqrt(kk2(k,j,i))
                phase = 2*pi*randu(seed)
                state(n)%ccomp(k,j,i) = 2*(COS(phase)+im*SIN(phase))*dump
              ELSE
                state(n)%ccomp(k,j,i) = 0.
              ENDIF
            END DO
          END DO
        END DO
      ELSE
        DO i = ista,iend
          DO j = 1,ny
            DO k = 1,nz
              IF ((kk2(k,j,i).le.skup**2).and.(kk2(k,j,i).ge.skdn**2)) THEN
                dump = 1./sqrt(kk2(k,j,i))
                phase = 2*pi*randu(seed)
                state(n)%ccomp(k,j,i) = 2*(COS(phase)+im*SIN(phase))*dump
              ELSE
                state(n)%ccomp(k,j,i) = 0.
              ENDIF
            END DO
          END DO
        END DO
      ENDIF
      CALL variance(state(n)%ccomp,tmp,1)
      CALL MPI_BCAST(tmp,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
        DO j = 1,ny
          DO k = 1,nz
            state(n)%ccomp(k,j,i) = state(n)%ccomp(k,j,i)* &
                         f0(n-solver%PASSIVE+1)/sqrt(tmp)
          END DO
        END DO
      END DO
    end do
    class default
      stop "Force: This solver does not support passive scalars"
    end select
  end subroutine init_randomfps


  ! ===================================================================
  ! Update methods
  ! ===================================================================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Instantaneous phase shift
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine update_shiftfps(this, force, solver, state)
    use pseudospec_fluid
    use equationbase_mod
    use commtypes
    use random
    use status
    implicit none
    
    class(shiftupdt_fps),   intent(inout) :: this
    class   (forceBase),   intent   (in) :: force
    class(EquationBase),   intent   (in) :: solver
    type   (GStateComp),   intent(inout) :: state(:)
    complex(kind=GP)                     :: cdump
    real(kind=GP)                        :: phase
    integer                              :: n

    select type (solver)
    class is (VelocityBase)
      if (timef.eq.fstep) then
        do n = solver%PASSIVE, solver%PASSIVE+solver%numpassive_-1
          if (myrank.eq.0) phase = 2*pi*randu(seed)
          call MPI_BCAST(phase,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
          cdump = COS(phase)+im*SIN(phase)
          call phaseshift(state(n)%ccomp,cdump)
        end do
      endif
    class default
      error stop "This solver does not support velocity forcing"
    end select
  end subroutine update_shiftfps
end module force_passive
