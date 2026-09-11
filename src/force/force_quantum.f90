! =====================================================================
! NAME       : force_quantum.f90
! DESCRIPTION: Forcing methods for the quantum solvers (all solver
!              classes extending QuantumBase): the forcing of the
!              order parameter (a thermal bath, an additive term in
!              the forcing state), and the external potentials, which
!              are not forces but auxiliary arrays of the solver
!              (solver%vpot_, vlinx_, vliny_) set and updated here.
!
! Forces avaliable:
!   null_fq    : Null forcing of the order parameter
!   thermal_fq : Thermal bath, random forcing with amplitude set by
!                k.T (kttherm) of the solver
!   cyltrap_vq : Harmonic cylindrical (cigar) trapping potential
!                V = V0 (x^2+y^2), with the linear ramps for the
!                angular momentum operator in rotating solvers
!
! Update methods available:
!   constant_fq: Constant forcing (no method)
!   renew_fq   : New random draw of the thermal forcing every fstep
!   constant_vq: Constant potential (no method)
!
! DATE       : 09/11/26 (PDM)
! =====================================================================
module force_quantum
  USE forcebase_mod
  IMPLICIT NONE

  ! ================= Forcing functions supported =====================
  type, extends(forceBase) :: forceNull_fq
    contains
      procedure :: init_GForce => init_nullfq
  end type forceNull_fq
  type, extends(forceBase) :: forceThermal_fq
    contains
      procedure :: init_GForce => init_thermalfq
  end type forceThermal_fq
  type, extends(forceBase) :: forceCyltrap_vq
    contains
      procedure :: init_GForce => init_cyltrapvq
  end type forceCyltrap_vq
  ! ================= Update methods supported =======================
  type, extends(forceUpdt) :: renewupdt_fq
    contains
      procedure :: update_GForce => update_renewfq
  end type renewupdt_fq

CONTAINS

  ! ===================================================================
  ! Forcing functions
  ! ===================================================================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Null forcing
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_nullfq(this,solver,state)
    use gstate_mod
    use equationbase_mod
    use grid
    use mpivars
    implicit none

    class(forceNull_fq), intent   (in) :: this
    class(EquationBase), intent(inout) :: solver
    type   (GStateComp), intent(inout) :: state(:)
    integer                            :: i,j,k

    select type (solver)
    class is (QuantumBase)
!$omp parallel do collapse(2) private (k)
      DO i = ista,iend
         DO j = 1,ny
            DO CONCURRENT (k=1:nz)
              state(solver%ZFUNC  )%ccomp(k,j,i) = 0.0_GP
              state(solver%ZFUNC+1)%ccomp(k,j,i) = 0.0_GP
            END DO
         END DO
      END DO
    class default
      error stop "This solver does not support order parameter forcing"
    end select
  end subroutine init_nullfq


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Thermal forcing: Gaussian random noise in all the modes
  !! (with Hermitian symmetry), with total amplitude
  !! sqrt(kttherm/dt). kttherm is the k.T of the solver (&GL
  !! namelist); a new draw is made when called again.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_thermalfq(this,solver,state)
    use gstate_mod
    use equationbase_mod
    use gl_mod
    use pseudospec_scalar, only: variance
    use random
    use grid
    use kes
    use ali
    use var
    use status
    use commtypes
    use mpivars
    implicit none

    class(forceThermal_fq), intent   (in) :: this
    class   (EquationBase), intent(inout) :: solver
    type      (GStateComp), intent(inout) :: state(:)
    double precision                      :: tmp,tmq,tmr
    real(kind=GP)                         :: kttherm,rmp,rmq,fac
    integer                               :: i,j,k

    select type (solver)
    class is (QuantumBase)
      kttherm = 0.0_GP
      select type (solver)
      type is (GLSolver)
        kttherm = solver%traits_%kttherm
      end select
      associate (fre => state(solver%ZFUNC)%ccomp, fim => state(solver%ZFUNC+1)%ccomp)
      ! Random draws in the order of the modes (deterministic for a
      ! given seed and decomposition)
      IF (ista.eq.1) THEN
         CALL randn_cmplx(rmp,rmq,seed)
         fre(1,1,1) = rmp
         fim(1,1,1) = rmq
         DO j = 2,ny/2+1
            IF ((kn2(1,j,1).le.kmax).and.(kn2(1,j,1).ge.tiny)) THEN
               CALL randn_cmplx(rmp,rmq,seed)
               fre(1,j,1) = (rmp+im*rmq)
               fre(1,ny-j+2,1) = conjg(fre(1,j,1))
               CALL randn_cmplx(rmp,rmq,seed)
               fim(1,j,1) = (rmp+im*rmq)
               fim(1,ny-j+2,1) = conjg(fim(1,j,1))
            ELSE
               fre(1,j,1) = 0.
               fre(1,ny-j+2,1) = 0.
               fim(1,j,1) = 0.
               fim(1,ny-j+2,1) = 0.
            ENDIF
         END DO
         DO k = 2,nz/2+1
            IF ((kn2(k,1,1).le.kmax).and.(kn2(k,1,1).ge.tiny)) THEN
               CALL randn_cmplx(rmp,rmq,seed)
               fre(k,1,1) = (rmp+im*rmq)
               fre(nz-k+2,1,1) = conjg(fre(k,1,1))
               CALL randn_cmplx(rmp,rmq,seed)
               fim(k,1,1) = (rmp+im*rmq)
               fim(nz-k+2,1,1) = conjg(fim(k,1,1))
            ELSE
               fre(k,1,1) = 0.
               fre(nz-k+2,1,1) = 0.
               fim(k,1,1) = 0.
               fim(nz-k+2,1,1) = 0.
            ENDIF
         END DO
         DO j = 2,ny
            DO k = 2,nz/2+1
            IF ((kn2(k,j,1).le.kmax).and.(kn2(k,j,1).ge.tiny)) THEN
               CALL randn_cmplx(rmp,rmq,seed)
               fre(k,j,1) = (rmp+im*rmq)
               fre(nz-k+2,ny-j+2,1) = conjg(fre(k,j,1))
               CALL randn_cmplx(rmp,rmq,seed)
               fim(k,j,1) = (rmp+im*rmq)
               fim(nz-k+2,ny-j+2,1) = conjg(fim(k,j,1))
            ELSE
               fre(k,j,1) = 0.
               fre(nz-k+2,ny-j+2,1) = 0.
               fim(k,j,1) = 0.
               fim(nz-k+2,ny-j+2,1) = 0.
            ENDIF
            END DO
         END DO
         DO i = 2,iend
            DO j = 1,ny
               DO k = 1,nz
               IF ((kn2(k,j,i).le.kmax).and.(kn2(k,j,i).ge.tiny)) THEN
                  CALL randn_cmplx(rmp,rmq,seed)
                  fre(k,j,i) = 2*(rmp+im*rmq)
                  CALL randn_cmplx(rmp,rmq,seed)
                  fim(k,j,i) = 2*(rmp+im*rmq)
               ELSE
                  fre(k,j,i) = 0.
                  fim(k,j,i) = 0.
               ENDIF
               END DO
            END DO
         END DO
      ELSE
         DO i = ista,iend
            DO j = 1,ny
               DO k = 1,nz
               IF ((kn2(k,j,i).le.kmax).and.(kn2(k,j,i).ge.tiny)) THEN
                  CALL randn_cmplx(rmp,rmq,seed)
                  fre(k,j,i) = 2*(rmp+im*rmq)
                  CALL randn_cmplx(rmp,rmq,seed)
                  fim(k,j,i) = 2*(rmp+im*rmq)
               ELSE
                  fre(k,j,i) = 0.
                  fim(k,j,i) = 0.
               ENDIF
               END DO
            END DO
        END DO
      ENDIF
      ! Renormalize to the amplitude sqrt(kttherm/dt)
      CALL variance(fre,tmp,1)
      CALL variance(fim,tmq,1)
      IF (myrank.eq.0) tmr = tmp+tmq
      CALL MPI_BCAST(tmr,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      fac = 0.0_GP
      IF (tmr.gt.0.0D0) fac = sqrt(kttherm/dt)/sqrt(real(tmr,kind=GP))
!$omp parallel do collapse(2) private (k)
      DO i = ista,iend
         DO j = 1,ny
            DO k = 1,nz
               fre(k,j,i) = fre(k,j,i)*fac
               fim(k,j,i) = fim(k,j,i)*fac
            END DO
         END DO
      END DO
      end associate
    class default
      error stop "This solver does not support order parameter forcing"
    end select
  end subroutine init_thermalfq


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Harmonic cylindrical (cigar) trap in x and y, V = V0
  !! (x^2 + y^2), and linear functions of the x and y
  !! coordinates multiplied by the rotation rate omegaz (for
  !! the angular momentum operator Lz = -i hbar (x d/dy - y
  !! d/dx)). The functions are filtered to satisfy the periodic
  !! boundary conditions, and their amplitudes corrected to
  !! have the right derivatives in the center of the domain.
  !! The potential divided by beta and the ramps (not
  !! normalized) are stored in the solver, and the normalized
  !! functions written to files 'Vtrap', 'Vlinx' and 'Vliny'
  !! with extension 'init'. V0 = m.w0^2/(2.hbar) with w0 the
  !! trapping frequency (from the &GL namelist), resulting in a
  !! characteristic lengthscale for the trap a0 = sqrt(hbar/
  !! (m.w0)) = (cspeed.lambda/(V0.sqrt(2)))^(1/4).
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_cyltrapvq(this,solver,state)
    use gstate_mod
    use equationbase_mod
    use iovar
    use grid
    use boxsize
    use kes
    use ali
    use var
    use fft
    use commtypes
    use mpivars
    implicit none

    class(forceCyltrap_vq), intent   (in) :: this
    class   (EquationBase), intent(inout) :: solver
    type      (GStateComp), intent(inout) :: state(:)
    complex(kind=GP), pointer, dimension(:,:,:) :: C1,C2,C3
    real   (kind=GP), pointer, dimension(:,:,:) :: R1,R2,R3
    real(kind=GP)                         :: rmp,rmq,rms,rmt,dump
    real(kind=GP)                         :: V0,beta,omegaz
    integer                               :: i,j,k
    logical                               :: bret

    select type (solver)
    class is (QuantumBase)
      V0     = solver%V0_
      beta   = solver%beta_
      omegaz = solver%omegaz_
      if ( V0 .le. 0.0_GP ) then
        if (myrank.eq.0) error stop 'cyltrap_vq: the trap requires V0 > 0 in the solver namelist'
      endif
      call solver%workspace_%get_complex_tmp(C1,bret)
      call solver%workspace_%get_complex_tmp(C2,bret)
      call solver%workspace_%get_complex_tmp(C3,bret)
      call solver%workspace_%get_real_tmp(R1,bret)
      call solver%workspace_%get_real_tmp(R2,bret)
      call solver%workspace_%get_real_tmp(R3,bret)
      ! We compute the linear functions and the parabola
!$omp parallel do collapse(2) private (i,rms,rmt)
      DO k = ksta,kend
         DO j = 1,ny
            rmt = pi*(2*real(j-1,kind=GP)/real(ny,kind=GP)-1.0_GP)
            DO i = 1,nx
               rms = pi*(2*real(i-1,kind=GP)/real(nx,kind=GP)-1.0_GP)
               solver%vlinx_(i,j,k) = rms              ! omegaz.x
               solver%vliny_(i,j,k) = rmt              ! omegaz.y
               solver%vpot_ (i,j,k) = rms**2 + rmt**2  ! V0.(x^2 + y^2)
            END DO
         END DO
      END DO
      ! We filter these functions to satisfy the periodic boundary conditions
      CALL fftp3d_real_to_complex(planrc,solver%vlinx_,C1,MPI_COMM_WORLD)
      CALL fftp3d_real_to_complex(planrc,solver%vliny_,C2,MPI_COMM_WORLD)
      CALL fftp3d_real_to_complex(planrc,solver%vpot_ ,C3,MPI_COMM_WORLD)
      rmp  = nx*Dkx/17.0_GP ! sigma_x: width of the filter in x
      rmq  = ny*Dky/17.0_GP ! sigma_y: width of the filter in y
      dump = 1.0_GP/(real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))
!$omp parallel do collapse(2) private (k,rms,rmt)
      DO i = ista,iend
         DO j = 1,ny
            rms = exp(-(kx(i)/rmp)**2/2)
            rmt = exp(-(ky(j)/rmq)**2/2)
            DO k = 1,nz
               IF (kn2(k,j,i).le.kmax) THEN
                  C1(k,j,i) = C1(k,j,i) * rms       * dump
                  C2(k,j,i) = C2(k,j,i)       * rmt * dump
                  C3(k,j,i) = C3(k,j,i) * rms * rmt * dump
               ELSE
                  C1(k,j,i) = 0.0_GP
                  C2(k,j,i) = 0.0_GP
                  C3(k,j,i) = 0.0_GP
               END IF
            END DO
         END DO
      END DO
      CALL fftp3d_complex_to_real(plancr,C1,solver%vlinx_,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,C2,solver%vliny_,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,C3,solver%vpot_ ,MPI_COMM_WORLD)
      ! After filtering, the amplitudes may be off. We correct them to
      ! have the right values of the derivatives in the center of the
      ! domain for the trapping potential and for the linear ramps.
      rmp = solver%vlinx_(nx/2+nx/4,ny/2,ksta) - solver%vlinx_(nx/2-nx/4,ny/2,ksta)
      rmq = solver%vliny_(nx/2,ny/2+ny/4,ksta) - solver%vliny_(nx/2,ny/2-ny/4,ksta)
      rmp = real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP)* &
            omegaz*pi*Lx/rmp ! unnormalized omegaz/slope_x
      rmq = real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP)* &
            omegaz*pi*Ly/rmq ! unnormalized omegaz/slope_y
      rms = solver%vpot_(nx/2-nx/4-1,ny/2,ksta) - solver%vpot_(nx/2-nx/4,ny/2,ksta) &
          - solver%vpot_(nx/2+nx/4-1,ny/2,ksta) + solver%vpot_(nx/2+nx/4,ny/2,ksta)
      rmt = solver%vpot_(nx/2,ny/2-ny/4-1,ksta) - solver%vpot_(nx/2,ny/2-ny/4,ksta) &
          - solver%vpot_(nx/2,ny/2+ny/4-1,ksta) + solver%vpot_(nx/2,ny/2+ny/4,ksta)
      rms = real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP)* &
            8*V0*pi**2/((rms*real(nx,kind=GP)/Lx**2 +           &
            rmt*real(ny,kind=GP)/Ly**2)*beta) ! 2V0/(second_deriv*beta)
!$omp parallel do collapse(2) private (i)
      DO k = ksta,kend
         DO j = 1,ny
            DO i = 1,nx
               solver%vlinx_(i,j,k) = solver%vlinx_(i,j,k)*rmp
               solver%vliny_(i,j,k) = solver%vliny_(i,j,k)*rmq
               solver%vpot_ (i,j,k) = solver%vpot_ (i,j,k)*rms
            END DO
         END DO
      END DO
      solver%haspot_ = .true.
      ! We write the normalized potentials and linear ramps to a file
      rmp = 1.0_GP/(real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))
      rmq = beta  /(real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))
!$omp parallel do collapse(2) private (i)
      DO k = ksta,kend
         DO j = 1,ny
            DO i = 1,nx
               R1(i,j,k) = solver%vlinx_(i,j,k)*rmp
               R2(i,j,k) = solver%vliny_(i,j,k)*rmp
               R3(i,j,k) = solver%vpot_ (i,j,k)*rmq
            END DO
         END DO
      END DO
      CALL io_write(1,solver%odir_,'Vlinx','init',solver%planio_,R1)
      CALL io_write(1,solver%odir_,'Vliny','init',solver%planio_,R2)
      CALL io_write(1,solver%odir_,'Vtrap','init',solver%planio_,R3)
      call solver%workspace_%free_real_tmp(R3)
      call solver%workspace_%free_real_tmp(R2)
      call solver%workspace_%free_real_tmp(R1)
      call solver%workspace_%free_complex_tmp(C3)
      call solver%workspace_%free_complex_tmp(C2)
      call solver%workspace_%free_complex_tmp(C1)
    class default
      error stop "This solver does not support external potentials"
    end select
  end subroutine init_cyltrapvq


  ! ===================================================================
  ! Update methods
  ! ===================================================================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! New random draw of the forcing every fstep steps
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine update_renewfq(this, force, solver, state)
    use equationbase_mod
    use gstate_mod
    use status
    implicit none

    class(renewupdt_fq),   intent(inout) :: this
    class   (forceBase),   intent   (in) :: force
    class(EquationBase),   intent(inout) :: solver
    type   (GStateComp),   intent(inout) :: state(:)

    if (timef.eq.fstep) then
      call force%init_GForce(solver,state)
    endif
  end subroutine update_renewfq

end module force_quantum
