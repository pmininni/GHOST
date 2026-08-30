! =====================================================================
! NAME       : force_active.f90
! DESCRIPTION: Forcing methods for active scalars in all ActiveScalarBase
!              solver classes. Each force needs an initial set up, and
!              an update.
!
! Forces avaliable:
!   null_fas    : Null forcing of the active scalar
!   puff_fas    : Localized forcing of the active scalar
!   random_fas  : Random forcing
!
! Update methods available:
!   constant_fas: Constant forcing (no method)
!   shift_fas   : Instantaneous phase shift
!
! DATE       : 04/08/26 (JBG)
! =====================================================================

module force_active
  USE forcebase_mod

  IMPLICIT NONE

  ! ================= Forcing functions supported =====================
  type, extends(forceBase) :: forceNull_fas
    contains
      procedure :: init_GForce => init_nullfas
  end type forceNull_fas 
  type, extends(forceBase) :: forcePuff_fas
    contains
      procedure :: init_GForce => init_pufffas
  end type forcePuff_fas 
  type, extends(forceBase) :: forceRandom_fas
    contains
      procedure :: init_GForce => init_randomfas
  end type forceRandom_fas

  ! ================= Update methods supported =======================
  type, extends(forceUpdt) :: shiftupdt_fas
    contains
      procedure :: update_GForce => update_shiftfas
  end type shiftupdt_fas

CONTAINS

  ! ===================================================================
  ! Forcing functions
  ! ===================================================================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Null forcing
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_nullfas(this,solver,state)
    use gstate_mod
    use hd_mod
    use grid
    use mpivars
!$  use threads
    implicit none
    
    class(forceNull_fas), intent   (in) :: this
    class (EquationBase), intent   (in) :: solver
    type    (GStateComp), intent(inout) :: state(:)
    integer                             :: i,j,k,n

    select type (solver)
    class is (ActiveScalarBase)
      do n = solver%ACTIVESC, solver%ACTIVESC+solver%numactivesc_-1    
!$omp parallel do collapse(2) private (k)
        DO i = ista,iend
          DO j = 1,ny
            DO CONCURRENT (k=1:nz)
              state(n)%ccomp(k,j,i) = 0.0_GP
            END DO
          END DO
        END DO
      end do
    class default
      stop "Force: This solver does not support active scalars"
    end select
  end subroutine init_nullfas


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Forcing in a Gaussian ball ('puff') centered at
  !! (x0, y0, z0) with a FWHM of r0
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   f0        : vector with the amplitudes of the forcing
  !   x0, y0, z0: vectors with the centers of the puffas
  !   r0        : vector with the radii of the puffas
  subroutine init_pufffas(this,solver,state)
    use gstate_mod
    use hd_mod
    use grid
    use mpivars
    use commtypes
    use fft
    use pseudospec_scalar
!$  use threads
    implicit none
    
    class(forcePuff_fas), intent   (in) :: this
    class (EquationBase), intent   (in) :: solver
    type    (GStateComp), intent(inout) :: state(:)
    real       (kind=GP), pointer       :: R1(:,:,:)
    real       (kind=GP), allocatable, dimension(:)  :: f0,x0,y0,z0,r0
    double precision                    :: tmp
    integer                             :: i,j,k,n
    logical                             :: bret

    namelist/ puff_fas / f0,x0,y0,z0,r0
    CALL solver%workspace_%get_real_tmp(R1,bret)    
    select type (solver)
    class is (ActiveScalarBase)
    allocate ( f0(solver%numactivesc_) )
    allocate ( x0(solver%numactivesc_) )
    allocate ( y0(solver%numactivesc_) )
    allocate ( z0(solver%numactivesc_) )
    allocate ( r0(solver%numactivesc_) )
    if ( myrank .eq. 0 ) then
      open(1,file=solver%infile_,status='unknown',form="formatted")
      read(1,NML=puff_fas)
      close(1)
    endif
    call mpi_bcast(f0,solver%numactivesc_,GC_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(x0,solver%numactivesc_,GC_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(y0,solver%numactivesc_,GC_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(z0,solver%numactivesc_,GC_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(r0,solver%numactivesc_,GC_REAL,0,MPI_COMM_WORLD,ierr)
    do n = solver%ACTIVESC, solver%ACTIVESC+solver%numactivesc_-1    
!$omp parallel do collapse(2) private (i)
      do k = ksta,kend
        do j = 1,ny
          do concurrent (i=1:nx)
            tmp = (real(i-1,kind=GP)/real(nx-1,kind=GP)-x0(n-solver%ACTIVESC+1))**2 &
                + (real(j-1,kind=GP)/real(ny-1,kind=GP)-y0(n-solver%ACTIVESC+1))**2 &
                + (real(k-1,kind=GP)/real(nz-1,kind=GP)-z0(n-solver%ACTIVESC+1))**2 
            R1(i,j,k) = exp(-tmp**2/r0(n-solver%ACTIVESC+1)**2)
          end do
        end do
      end do
      call fftp3d_real_to_complex(planrc,R1,state(n)%ccomp,MPI_COMM_WORLD)
      call variance(state(n)%ccomp,tmp,1)
      call mpi_bcast(tmp,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!$omp parallel do collapse(2) private (k)
      do i = ista,iend
        do j = 1,ny
          do concurrent (k=1:nz)
            state(n)%ccomp(k,j,i) = state(n)%ccomp(k,j,i)* &
                         f0(n-solver%ACTIVESC+1)/sqrt(tmp)
          end do
        end do
      end do
    end do
    class default
      stop "Force: This solver does not support active scalars"
    end select
    call solver%workspace_%free_real_tmp(R1)
  end subroutine init_pufffas

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Random forcing
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   f0 : vector with the amplitudes of the forcing
  !   kdn: vector with minimum forcing wave number
  !   kup: vector with maximum forcing wave number
  subroutine init_randomfas(this,solver,state)
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
    
    class(forceRandom_fas), intent   (in) :: this
    class   (EquationBase), intent   (in) :: solver
    type      (GStateComp), intent(inout) :: state(:)
    real         (kind=GP), allocatable, dimension(:) :: f0,kdn,kup
    real         (kind=GP)                            :: skup,skdn
    real         (kind=GP)                            :: dump,phase
    double precision                      :: tmp
    integer                               :: i,j,k,n
    logical                               :: bret

    namelist/ random_fas / f0,kup,kdn
    select type (solver)
    class is (ActiveScalarBase)
    allocate ( f0 (solver%numactivesc_) )
    allocate ( kdn(solver%numactivesc_) )
    allocate ( kup(solver%numactivesc_) )
    if ( myrank .eq. 0 ) then
      open(1,file=solver%infile_,status='unknown',form="formatted")
      read(1,NML=random_fas)
      close(1)
    endif
    call mpi_bcast(f0 ,solver%numactivesc_,GC_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(kdn,solver%numactivesc_,GC_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(kup,solver%numactivesc_,GC_REAL,0,MPI_COMM_WORLD,ierr)
    do n = solver%ACTIVESC, solver%ACTIVESC+solver%numactivesc_-1
      skdn = kdn(n-solver%ACTIVESC+1)
      skup = kup(n-solver%ACTIVESC+1)
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
!$omp parallel do collapse(2) private (k)
      DO i = ista,iend
        DO j = 1,ny
          DO CONCURRENT (k=1:nz)
            state(n)%ccomp(k,j,i) = state(n)%ccomp(k,j,i)* &
                         f0(n-solver%ACTIVESC+1)/sqrt(tmp)
          END DO
        END DO
      END DO
    end do
    class default
      stop "Force: This solver does not support active scalars"
    end select
  end subroutine init_randomfas


  ! ===================================================================
  ! Update methods
  ! ===================================================================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Instantaneous phase shift
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine update_shiftfas(this, force, solver, state)
    use pseudospec_fluid
    use equationbase_mod
    use commtypes
    use random
    use status
    implicit none
    
    class(shiftupdt_fas),   intent(inout) :: this
    class   (forceBase),   intent   (in) :: force
    class(EquationBase),   intent   (in) :: solver
    type   (GStateComp),   intent(inout) :: state(:)
    complex(kind=GP)                     :: cdump
    real(kind=GP)                        :: phase
    integer                              :: n

    select type (solver)
    class is (ActiveScalarBase)
      if (timef.eq.fstep) then
        do n = solver%ACTIVESC, solver%ACTIVESC+solver%numactivesc_-1
          if (myrank.eq.0) phase = 2*pi*randu(seed)
          call MPI_BCAST(phase,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
          cdump = COS(phase)+im*SIN(phase)
          call phaseshift(state(n)%ccomp,cdump)
        end do
      endif
    class default
      error stop "Force_updt: This solver does not support active scalar forcing"
    end select
  end subroutine update_shiftfas
end module force_active
