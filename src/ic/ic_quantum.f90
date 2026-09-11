! =====================================================================
! NAME       : ic_quantum.f90
! DESCRIPTION: Initial conditions for the order parameter (wave
!              function) z = zre + i zim of the quantum solvers (all
!              solver classes extending QuantumBase). Some ICs also
!              set the advective velocity companion of the solver
!              (solver%vx_, vy_, vz_, constant in time), used by the
!              GL solver; the GPE solver ignores it.
!
! On restarts (stat > 0), list read_z together with the IC that
! created the flow (e.g., "read_z;abc_z"): the ICs with a velocity
! companion then only regenerate the advective velocity, as the
! old code re-created the velocity on restart, and leave z as read.
!
! Initial conditions available:
!   read_z    : Reads zre, zim from input files numbered by stat
!   uniform_z : Uniform condensate, z = sqrt(rho0)
!   gaussian_z: Gaussian perturbation of the uniform condensate
!   abc_z     : Array of vortices following an ABC flow (+ ABC velocity)
!   tg_z      : Array of vortices following a Taylor-Green flow (+ TG
!               velocity)
!   ring_z    : Vortex ring (+ velocity of the ring)
!   trefoil_z : Trefoil vortex knot (+ velocity of the knot)
!   tworings_z: Two linked vortex rings (+ velocity of the rings)
!   trap_z    : Gaussian wave function for a cigar trapping potential
!
! DATE       : 09/11/26 (PDM)
! =====================================================================

module ic_quantum
  use icbase_mod
  use equationbase_mod

  implicit none

  ! ================= Initial conditions supported ====================
  type, extends(icBase) :: icRead_z
    contains
      procedure :: init_GState => init_readz
  end type icRead_z
  type, extends(icBase) :: icUniform_z
    contains
      procedure :: init_GState => init_uniformz
  end type icUniform_z
  type, extends(icBase) :: icGaussian_z
    contains
      procedure :: init_GState => init_gaussianz
  end type icGaussian_z
  type, extends(icBase) :: icAbc_z
    contains
      procedure :: init_GState => init_abcz
  end type icAbc_z
  type, extends(icBase) :: icTg_z
    contains
      procedure :: init_GState => init_tgz
  end type icTg_z
  type, extends(icBase) :: icRing_z
    contains
      procedure :: init_GState => init_ringz
  end type icRing_z
  type, extends(icBase) :: icTrefoil_z
    contains
      procedure :: init_GState => init_trefoilz
  end type icTrefoil_z
  type, extends(icBase) :: icTworings_z
    contains
      procedure :: init_GState => init_tworingsz
  end type icTworings_z
  type, extends(icBase) :: icTrap_z
    contains
      procedure :: init_GState => init_trapz
  end type icTrap_z

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Dealiases the order parameter (zeroes the modes with
  !! kn2 > kmax). The ICs built in real space have modes up to
  !! the Nyquist wavenumber; the solvers evolve only the modes
  !! with kn2 <= kmax, and the others must be removed (the
  !! old code did this at every time step, the new steppers
  !! do not modify the modes with zero time derivative)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine dealias_z(solver,state)
    use grid
    use mpivars
    use kes
    use ali
    class(QuantumBase), intent(in)    :: solver
    type  (GStateComp), intent(inout) :: state(:)
    integer                           :: i,j,k
!$omp parallel do collapse(2) private (k)
    do i = ista,iend
      do j = 1,ny
        do k = 1,nz
          if (kn2(k,j,i).gt.kmax) then
            state(solver%ZFUNC  )%ccomp(k,j,i) = 0.0_GP
            state(solver%ZFUNC+1)%ccomp(k,j,i) = 0.0_GP
          endif
        end do
      end do
    end do
  end subroutine dealias_z


  ! ===================================================================
  ! Initial conditions
  ! ===================================================================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Read the order parameter from restart files
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_readz(this,solver,state)
    use gstate_mod
    use iovar
    use status
    use filefmt
    use fft
    use commtypes
    implicit none

    class    (icRead_z), intent(in)            :: this
    class(EquationBase), intent(inout)         :: solver
    type   (GStateComp), intent(inout)         :: state(:)
    real(kind=GP), pointer, dimension(:,:,:)   :: R1
    logical                                    :: bret

    if ((stat .eq. 0) .and. (solver%myrank_ .eq. 0)) then
      error stop 'Cannot read files if starting a new run with stat=0'
    endif
    call solver%workspace_%get_real_tmp(R1,bret)
    select type (solver)
    class is (QuantumBase)
      tind = int(stat)
      write(ext, fmtext) tind
      call io_read(1,solver%idir_,trim(solver%sstate_(solver%ZFUNC  )),ext, &
                   solver%planio_,R1)
      call fftp3d_real_to_complex(planrc,R1,state(solver%ZFUNC  )%ccomp, &
                   MPI_COMM_WORLD)
      call io_read(1,solver%idir_,trim(solver%sstate_(solver%ZFUNC+1)),ext, &
                   solver%planio_,R1)
      call fftp3d_real_to_complex(planrc,R1,state(solver%ZFUNC+1)%ccomp, &
                   MPI_COMM_WORLD)
    class default
      error stop 'IC: This solver does not support order parameter ICs'
    end select
    call solver%workspace_%free_real_tmp(R1)
  end subroutine init_readz


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Uniform condensate at equilibrium, z = sqrt(rho0)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_uniformz(this,solver,state)
    use gstate_mod
    use grid
    use mpivars
    implicit none

    class (icUniform_z), intent(in)            :: this
    class(EquationBase), intent(inout)         :: solver
    type   (GStateComp), intent(inout)         :: state(:)
    integer                                    :: i,j,k

    select type (solver)
    class is (QuantumBase)
!$omp parallel do collapse(2) private (k)
      do i = ista,iend
        do j = 1,ny
          do concurrent (k=1:nz)
            state(solver%ZFUNC  )%ccomp(k,j,i) = 0.0_GP
            state(solver%ZFUNC+1)%ccomp(k,j,i) = 0.0_GP
          end do
        end do
      end do
      if (myrank .eq. 0) then
        state(solver%ZFUNC)%ccomp(1,1,1) = sqrt(solver%rho0_)* &
              real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP)
      endif
    class default
      error stop 'IC: This solver does not support order parameter ICs'
    end select
  end subroutine init_uniformz


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Gaussian perturbation to the uniform condensate along x
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   z0 : amplitude of the perturbation
  !   l0 : length of the perturbation (in 2.pi units)
  subroutine init_gaussianz(this,solver,state)
    use gstate_mod
    use grid
    use var
    use fft
    use commtypes
    use mpivars
    implicit none

    class(icGaussian_z), intent(in)            :: this
    class(EquationBase), intent(inout)         :: solver
    type   (GStateComp), intent(inout)         :: state(:)
    real(kind=GP), pointer, dimension(:,:,:)   :: R1
    real(kind=GP)                              :: z0,l0,rmp
    integer                                    :: i,j,k
    logical                                    :: bret

    namelist/ gaussian_z / z0,l0
    select type (solver)
    class is (QuantumBase)
      z0 = 0.01_GP; l0 = 1.0_GP
      if ( myrank .eq. 0 ) then
        open(1,file=solver%infile_,status='unknown',form="formatted")
        read(1,NML=gaussian_z)
        close(1)
      endif
      call mpi_bcast(z0,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(l0,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      call solver%workspace_%get_real_tmp(R1,bret)
      rmp = sqrt(solver%rho0_)
!$omp parallel do collapse(2) private (i)
      do k = ksta,kend
        do j = 1,ny
          do i = 1,nx
            R1(i,j,k) = rmp+z0*exp(-(2*pi*(real(i,kind=GP)-1)/ &
                        real(nx,kind=GP)-pi)**2/l0**2)
          end do
        end do
      end do
      call fftp3d_real_to_complex(planrc,R1,state(solver%ZFUNC)%ccomp, &
                                  MPI_COMM_WORLD)
!$omp parallel do collapse(2) private (k)
      do i = ista,iend
        do j = 1,ny
          do concurrent (k=1:nz)
            state(solver%ZFUNC+1)%ccomp(k,j,i) = 0.0_GP
          end do
        end do
      end do
      call dealias_z(solver,state)
      call solver%workspace_%free_real_tmp(R1)
    class default
      error stop 'IC: This solver does not support order parameter ICs'
    end select
  end subroutine init_gaussianz


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Array of vortices following an ABC flow. Also sets the
  !! ABC flow (not normalized) as the advective velocity.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   kdn  : minimum wave number (rounded to next integer)
  !   kup  : maximum wave number (rounded to next integer)
  !   A,B,C: amplitudes of the ABC flow
  subroutine init_abcz(this,solver,state)
    use gstate_mod
    use status, only: stat
    use grid
    use boxsize
    use ali
    use var
    use fft
    use commtypes
    use mpivars
    implicit none

    class     (icAbc_z), intent(in)            :: this
    class(EquationBase), intent(inout)         :: solver
    type   (GStateComp), intent(inout)         :: state(:)
    real(kind=GP), pointer, dimension(:,:,:)   :: R1,R2,R3
    real(kind=GP)                              :: kdn,kup,A,B,C
    real(kind=GP)                              :: alpha,rmp
    complex(kind=GP)                           :: cdump
    integer                                    :: i,j,k,ki
    logical                                    :: bret

    namelist/ abc_z / kdn,kup,A,B,C
    select type (solver)
    class is (QuantumBase)
      if ( (abs(Lx-Ly).gt.tiny).or.(abs(Lx-Lz).gt.tiny) ) then
        if (myrank.eq.0) error stop 'ABC initial conditions require Lx=Ly=Lz'
      endif
      kdn = 1.0_GP; kup = 1.0_GP; A = 1.0_GP; B = 1.0_GP; C = 1.0_GP
      if ( myrank .eq. 0 ) then
        open(1,file=solver%infile_,status='unknown',form="formatted")
        read(1,NML=abc_z)
        close(1)
      endif
      call mpi_bcast(kdn,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(kup,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(A  ,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(B  ,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(C  ,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      alpha = solver%alpha_
      call solver%workspace_%get_real_tmp(R1,bret)
      call solver%workspace_%get_real_tmp(R2,bret)
      call solver%workspace_%get_real_tmp(R3,bret)
      ! Advective velocity: ABC flow (not normalized)
!$omp parallel do collapse(2) private (i,ki)
      do k = ksta,kend
        do j = 1,ny
          do i = 1,nx
            R1(i,j,k) = 0.0_GP
            R2(i,j,k) = 0.0_GP
            R3(i,j,k) = 0.0_GP
            do ki = INT(kdn),INT(kup)
              R1(i,j,k) = R1(i,j,k)+(B*COS(2*pi*ki*(real(j,kind=GP)-1)/  &
                   real(ny,kind=GP))+C*SIN(2*pi*ki*(real(k,kind=GP)-1)/  &
                   real(nz,kind=GP)))
              R2(i,j,k) = R2(i,j,k)+(A*SIN(2*pi*ki*(real(i,kind=GP)-1)/  &
                   real(nx,kind=GP))+C*COS(2*pi*ki*(real(k,kind=GP)-1)/  &
                   real(nz,kind=GP)))
              R3(i,j,k) = R3(i,j,k)+(A*COS(2*pi*ki*(real(i,kind=GP)-1)/  &
                   real(nx,kind=GP))+B*SIN(2*pi*ki*(real(j,kind=GP)-1)/  &
                   real(ny,kind=GP)))
            end do
          end do
        end do
      end do
      call fftp3d_real_to_complex(planrc,R1,solver%vx_,MPI_COMM_WORLD)
      call fftp3d_real_to_complex(planrc,R2,solver%vy_,MPI_COMM_WORLD)
      call fftp3d_real_to_complex(planrc,R3,solver%vz_,MPI_COMM_WORLD)
      solver%hasadv_ = .true.
      if ( stat .gt. 0 ) then ! restart: z is read by read_z
        call solver%workspace_%free_real_tmp(R3)
        call solver%workspace_%free_real_tmp(R2)
        call solver%workspace_%free_real_tmp(R1)
        return
      endif
      ! Order parameter: phase following the flow, quantized circulation
      rmp = sqrt(solver%omegag_/solver%beta_)
!$omp parallel do collapse(2) private (i,ki,cdump)
      do k = ksta,kend
        do j = 1,ny
          do i = 1,nx
            cdump = cmplx(rmp,0.0_GP,kind=GP)
            do ki = INT(kdn),INT(kup)
              cdump = cdump*exp(im*2*pi*(real(i,kind=GP)-1)/real(nx,kind=GP)* &
                (int((B*COS(2*pi*ki*(real(j,kind=GP)-1)/real(ny,kind=GP)))/   &
                (2.0_GP*alpha)+0.5_GP)                                         &
                +int((C*SIN(2*pi*ki*(real(k,kind=GP)-1)/real(nz,kind=GP)))/   &
                (2.0_GP*alpha)+0.5_GP)) )
              cdump = cdump*exp(im*2*pi*(real(j,kind=GP)-1)/real(ny,kind=GP)* &
                (int((A*SIN(2*pi*ki*(real(i,kind=GP)-1)/real(nx,kind=GP)))/   &
                (2.0_GP*alpha)+0.5_GP)                                         &
                +int((C*COS(2*pi*ki*(real(k,kind=GP)-1)/real(nz,kind=GP)))/   &
                (2.0_GP*alpha)+0.5_GP)) )
              cdump = cdump*exp(im*2*pi*(real(k,kind=GP)-1)/real(nz,kind=GP)* &
                (int((A*COS(2*pi*ki*(real(i,kind=GP)-1)/real(nx,kind=GP)))/   &
                (2.0_GP*alpha)+0.5_GP)                                         &
                +int((B*SIN(2*pi*ki*(real(j,kind=GP)-1)/real(ny,kind=GP)))/   &
                (2.0_GP*alpha)+0.5_GP)) )
            end do
            R1(i,j,k) =  real(cdump)
            R2(i,j,k) = aimag(cdump)
          end do
        end do
      end do
      call fftp3d_real_to_complex(planrc,R1,state(solver%ZFUNC  )%ccomp,MPI_COMM_WORLD)
      call fftp3d_real_to_complex(planrc,R2,state(solver%ZFUNC+1)%ccomp,MPI_COMM_WORLD)
      call dealias_z(solver,state)
      call solver%workspace_%free_real_tmp(R3)
      call solver%workspace_%free_real_tmp(R2)
      call solver%workspace_%free_real_tmp(R1)
    class default
      error stop 'IC: This solver does not support order parameter ICs'
    end select
  end subroutine init_abcz


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Array of vortices following a Taylor-Green flow. Also sets
  !! the TG flow (not normalized) as the advective velocity.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   kdn  : minimum wave number (rounded to next integer)
  !   kup  : maximum wave number (rounded to next integer)
  subroutine init_tgz(this,solver,state)
    use gstate_mod
    use status, only: stat
    use grid
    use boxsize
    use ali
    use var
    use fft
    use commtypes
    use mpivars
    implicit none

    class      (icTg_z), intent(in)            :: this
    class(EquationBase), intent(inout)         :: solver
    type   (GStateComp), intent(inout)         :: state(:)
    real(kind=GP), pointer, dimension(:,:,:)   :: R1,R2
    real(kind=GP)                              :: kdn,kup
    real(kind=GP)                              :: alpha,lambda,dump
    real(kind=GP)                              :: rmp,rmq,rms,rmt,rm1,rm2
    complex(kind=GP)                           :: cdump,cdumq
    integer                                    :: i,j,k,ki
    logical                                    :: bret

    namelist/ tg_z / kdn,kup
    select type (solver)
    class is (QuantumBase)
      if ( abs(Lx-Ly).gt.tiny ) then
        if (myrank.eq.0) error stop 'TG initial conditions require at least Lx=Ly'
      endif
      kdn = 1.0_GP; kup = 1.0_GP
      if ( myrank .eq. 0 ) then
        open(1,file=solver%infile_,status='unknown',form="formatted")
        read(1,NML=tg_z)
        close(1)
      endif
      call mpi_bcast(kdn,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(kup,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      alpha  = solver%alpha_
      lambda = sqrt(alpha/solver%omegag_)      ! coherence length xi
      call solver%workspace_%get_real_tmp(R1,bret)
      call solver%workspace_%get_real_tmp(R2,bret)
      ! Advective velocity: TG flow (not normalized)
!$omp parallel do collapse(2) private (i,ki)
      do k = ksta,kend
        do j = 1,ny
          do i = 1,nx
            R1(i,j,k) = 0.0_GP
            R2(i,j,k) = 0.0_GP
            do ki = INT(kdn),INT(kup)
              R1(i,j,k) = R1(i,j,k)+SIN(2*pi*ki*(real(i,kind=GP)-1)/ &
                         real(nx,kind=GP))*COS(2*pi*ki*(real(j,kind=GP)-1)/ &
                         real(ny,kind=GP))*COS(2*pi*ki*(real(k,kind=GP)-1)/ &
                         real(nz,kind=GP))
              R2(i,j,k) = R2(i,j,k)-COS(2*pi*ki*(real(i,kind=GP)-1)/ &
                         real(nx,kind=GP))*SIN(2*pi*ki*(real(j,kind=GP)-1)/ &
                         real(ny,kind=GP))*COS(2*pi*ki*(real(k,kind=GP)-1)/ &
                         real(nz,kind=GP))
            end do
          end do
        end do
      end do
      call fftp3d_real_to_complex(planrc,R1,solver%vx_,MPI_COMM_WORLD)
      call fftp3d_real_to_complex(planrc,R2,solver%vy_,MPI_COMM_WORLD)
!$omp parallel do collapse(2) private (k)
      do i = ista,iend
        do j = 1,ny
          do concurrent (k=1:nz)
            solver%vz_(k,j,i) = 0.0_GP
          end do
        end do
      end do
      solver%hasadv_ = .true.
      if ( stat .gt. 0 ) then ! restart: z is read by read_z
        call solver%workspace_%free_real_tmp(R2)
        call solver%workspace_%free_real_tmp(R1)
        return
      endif
      ! Order parameter: superposition of vortex filaments
      dump = 1.0_GP/sqrt(2.0_GP)
      rmp = sqrt(solver%omegag_/solver%beta_)
!$omp parallel do collapse(2) private (i)
      do k = ksta,kend
        do j = 1,ny
          do i = 1,nx
            R1(i,j,k) = rmp
            R2(i,j,k) = rmp
          end do
        end do
      end do
      do ki = INT(kdn),INT(kup)
!$omp parallel do private (i,j,rmp,rmq,rm1,rm2,rms,rmt,cdump,cdumq)
        do k = ksta,kend
          rmp = sqrt(2*abs(cos(2*pi*ki*(real(k,kind=GP)-1)/real(nz,kind=GP))))
          rmq = rmp*sign(1.0_GP,cos(2*pi*ki*(real(k,kind=GP)-1)/real(nz,kind=GP)))
          do j = 1,ny
            do i = 1,nx
              rm1 = cos(2*pi*ki*(real(i,kind=GP)-1)/real(nx,kind=GP))*rmp
              rm2 = cos(2*pi*ki*(real(j,kind=GP)-1)/real(ny,kind=GP))*rmq
              rms = 1.0_GP/sqrt(((rm1-dump)**2+rm2**2)  &
                               *((rm2-dump)**2+rm1**2)  &
                               *((rm1+dump)**2+rm2**2)  &
                               *((rm2+dump)**2+rm1**2))
              rmt = tanh(dump*sqrt((rm1-dump)**2+rm2**2)/lambda) &
                   *tanh(dump*sqrt((rm2-dump)**2+rm1**2)/lambda) &
                   *tanh(dump*sqrt((rm1+dump)**2+rm2**2)/lambda) &
                   *tanh(dump*sqrt((rm2+dump)**2+rm1**2)/lambda)
              cdump = ((rm1-dump)+im*rm2)  &
                      *(rm1+im*(rm2-dump)) &
                     *((rm1+dump)+im*rm2)  &
                      *(rm1+im*(rm2+dump))
              cdumq = (cdump*rms*rmt)**int(1.0_GP/(2*pi*alpha*ki))
              R1(i,j,k) = R1(i,j,k)*real(cdumq)
              R2(i,j,k) = R2(i,j,k)*aimag(cdumq)
            end do
          end do
        end do
      end do
      call fftp3d_real_to_complex(planrc,R1,state(solver%ZFUNC  )%ccomp,MPI_COMM_WORLD)
      call fftp3d_real_to_complex(planrc,R2,state(solver%ZFUNC+1)%ccomp,MPI_COMM_WORLD)
      call dealias_z(solver,state)
      call solver%workspace_%free_real_tmp(R2)
      call solver%workspace_%free_real_tmp(R1)
    class default
      error stop 'IC: This solver does not support order parameter ICs'
    end select
  end subroutine init_tgz


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Vortex ring, trefoil knot and two linked rings: the three
  !! ICs share the same construction (init_knotz)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_ringz(this,solver,state)
    use gstate_mod
    implicit none
    class    (icRing_z), intent(in)    :: this
    class(EquationBase), intent(inout) :: solver
    type   (GStateComp), intent(inout) :: state(:)
    call init_knotz(1,solver,state)
  end subroutine init_ringz

  subroutine init_trefoilz(this,solver,state)
    use gstate_mod
    implicit none
    class (icTrefoil_z), intent(in)    :: this
    class(EquationBase), intent(inout) :: solver
    type   (GStateComp), intent(inout) :: state(:)
    call init_knotz(2,solver,state)
  end subroutine init_trefoilz

  subroutine init_tworingsz(this,solver,state)
    use gstate_mod
    implicit none
    class(icTworings_z), intent(in)    :: this
    class(EquationBase), intent(inout) :: solver
    type   (GStateComp), intent(inout) :: state(:)
    call init_knotz(3,solver,state)
  end subroutine init_tworingsz


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Position and tangent of the vortex filaments at the
  !! parameter s in [0,2.pi): icurve=1 ring, 2 trefoil, 3 two
  !! linked rings (nc curves). r0 is the size of the knot.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine knot_curve(icurve,r0,s,nc,pos,tan)
    use var
    implicit none
    integer      , intent (in) :: icurve
    real(kind=GP), intent (in) :: r0,s
    integer      , intent(out) :: nc
    real(kind=GP), intent(out) :: pos(3,2),tan(3,2)
    select case (icurve)
    case (1) ! ring in the plane z = pi
      nc = 1
      pos(:,1) = (/ pi+r0*cos(s), pi+r0*sin(s), pi /)
      tan(:,1) = (/ -r0*sin(s), r0*cos(s), 0.0_GP /)
    case (2) ! trefoil
      nc = 1
      pos(:,1) = (/ pi+(r0*(sin(s)+2.0_GP*sin(2.0_GP*s)))/3.0_GP,   &
                    pi+(r0*(cos(s)-2.0_GP*cos(2.0_GP*s)))/3.0_GP,   &
                    pi-(r0*sin(3.0_GP*s))/3.0_GP /)
      tan(:,1) = (/ (r0*(cos(s)+4.0_GP*cos(2.0_GP*s)))/3.0_GP,      &
                    (r0*(-sin(s)+4.0_GP*sin(2.0_GP*s)))/3.0_GP,     &
                    -(r0*cos(3.0_GP*s)) /)
    case default ! two linked rings
      nc = 2
      pos(:,1) = (/ pi-r0/2+r0*cos(s), pi+r0*sin(s), pi /)
      tan(:,1) = (/ -r0*sin(s), r0*cos(s), 0.0_GP /)
      pos(:,2) = (/ pi+r0/2+r0*cos(s), pi, pi+r0*sin(s) /)
      tan(:,2) = (/ -r0*sin(s), 0.0_GP, r0*cos(s) /)
    end select
  end subroutine knot_curve


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Vortex knots: the velocity of the filaments (with a
  !! Gaussian filter) is set as the advective velocity, the
  !! phase of the order parameter is integrated along the
  !! velocity, and the density is set from the distance to
  !! the filaments (Pade approximant of the vortex profile).
  !! Requires Lx=Ly=Lz=2.pi and isotropic grids.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   r0 : radius (size) of the vortex knot
  !   sig: width of the Gaussian filter of the velocity, in
  !        units of the coherence length
  subroutine init_knotz(icurve,solver,state)
    use gstate_mod
    use status, only: stat
    use pseudospec_fluid, only: rotor3, copy3
    use grid
    use boxsize
    use kes
    use ali
    use var
    use fft
    use commtypes
    use mpivars
    implicit none

    integer            , intent(in)            :: icurve
    class(EquationBase), intent(inout)         :: solver
    type   (GStateComp), intent(inout)         :: state(:)
    real   (kind=GP), pointer, dimension(:,:,:) :: R1,R2,R3
    complex(kind=GP), pointer, dimension(:,:,:) :: C1,C2,C3,C4,C5,C6
    real(kind=GP)                              :: r0,sig
    real(kind=GP)                              :: alpha,lambda,rho
    real(kind=GP)                              :: rmp,rmq,rms,rmt,rm1,tmq
    real(kind=GP)                              :: pos(3,2),tan(3,2),dst(3)
    complex(kind=GP)                           :: cdump,cphs
    integer                                    :: i,j,k,ki,ic,nc,nsc
    logical                                    :: bret

    namelist/ ring_z     / r0,sig
    namelist/ trefoil_z  / r0,sig
    namelist/ tworings_z / r0,sig
    select type (solver)
    class is (QuantumBase)
      if ( (abs(Lx-1.0_GP).gt.tiny).or.(abs(Ly-1.0_GP).gt.tiny).or.  &
           (abs(Lz-1.0_GP).gt.tiny).or.(nx.ne.ny).or.(nx.ne.nz).or.(ny.ne.nz) ) then
        if (myrank.eq.0) &
          error stop 'Quantum knots require Lx=Ly=Lz=2.pi and isotropic grids'
      endif
      r0 = 1.0_GP; sig = 1.0_GP
      if ( myrank .eq. 0 ) then
        open(1,file=solver%infile_,status='unknown',form="formatted")
        if (icurve.eq.1) read(1,NML=ring_z)
        if (icurve.eq.2) read(1,NML=trefoil_z)
        if (icurve.eq.3) read(1,NML=tworings_z)
        close(1)
      endif
      call mpi_bcast(r0 ,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(sig,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      alpha  = solver%alpha_
      lambda = sqrt(alpha/solver%omegag_)      ! coherence length xi
      rho    = sqrt(solver%omegag_/solver%beta_)
      call solver%workspace_%get_complex_tmp(C1,bret)
      call solver%workspace_%get_complex_tmp(C2,bret)
      call solver%workspace_%get_complex_tmp(C3,bret)
      call solver%workspace_%get_complex_tmp(C4,bret)
      call solver%workspace_%get_complex_tmp(C5,bret)
      call solver%workspace_%get_complex_tmp(C6,bret)
      call solver%workspace_%get_real_tmp(R1,bret)
      call solver%workspace_%get_real_tmp(R2,bret)
      call solver%workspace_%get_real_tmp(R3,bret)
      nsc = 4*nx                        ! points along the filaments
      ! Vector potential of the filaments, -Lap^(-1) of the vorticity
      rmp = 4.0_GP*pi*alpha/(2.0_GP*pi)**3
      rms = 2.0_GP*pi/real(nsc,kind=GP)
!$omp parallel do collapse(2) private (k,ki,ic,nc,rmq,rmt,pos,tan,cphs)
      do i = ista,iend
        do j = 1,ny
          do k = 1,nz
            C1(k,j,i) = 0.0_GP
            C2(k,j,i) = 0.0_GP
            C3(k,j,i) = 0.0_GP
            if ((i.eq.1).and.(j.eq.1).and.(k.eq.1)) then
              rmt = 0.0_GP              ! 0 for k=0
            else if ( kk2(k,j,i).lt.kmax ) then
              rmt = 1.0_GP/kk2(k,j,i)   ! -Laplac^(-1) in other case
            else
              rmt = 0.0_GP              ! 0 for k>=kmax
            endif
            do ki = 1,nsc
              rmq = rms*real(ki-1,kind=GP)
              call knot_curve(icurve,r0,rmq,nc,pos,tan)
              do ic = 1,nc
                cphs = rmp*rms*rmt*exp(-im*(kx(i)*pos(1,ic)+ky(j)*pos(2,ic)+ &
                                            kz(k)*pos(3,ic)))
                C1(k,j,i) = C1(k,j,i) + cphs*tan(1,ic)
                C2(k,j,i) = C2(k,j,i) + cphs*tan(2,ic)
                C3(k,j,i) = C3(k,j,i) + cphs*tan(3,ic)
              end do
            end do
          end do
        end do
      end do
      call rotor3(C2,C3,C4,1)
      call rotor3(C1,C3,C5,2)
      call rotor3(C1,C2,C6,3)
      ! Filter the velocity by exp(-k^2/(2*k0^2)) and normalize by N^3
      ! (the advective velocity has the normalization of the states,
      ! C4,C5,C6 keep the true Fourier amplitudes)
      rmq = .5_GP*lambda**2/sig**2
      rms = real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP)
!$omp parallel do collapse(2) private (k,rmp)
      do i = ista,iend
        do j = 1,ny
          do k = 1,nz
            rmp = exp(-kk2(k,j,i)*rmq)*rms
            solver%vx_(k,j,i) = C4(k,j,i)*rmp
            solver%vy_(k,j,i) = C5(k,j,i)*rmp
            solver%vz_(k,j,i) = C6(k,j,i)*rmp
          end do
        end do
      end do
      solver%hasadv_ = .true.
      ! Velocity in physical space to correct the mean flow: the line
      ! integrals of the phase along each direction are made multiples
      ! of 2.pi by removing a uniform velocity
      call copy3(C6,C3)
      call fftp3d_complex_to_real(plancr,C6,R3,MPI_COMM_WORLD) ! vz
      do ki = 0,nprocs-1
        if ( ki.eq.myrank ) then
          if ( myrank.eq.0 ) then
            rmp = 0.0_GP
          else
            call MPI_RECV(rmp,1,GC_REAL,myrank-1,1,MPI_COMM_WORLD, &
                          MPI_STATUS_IGNORE,ierr)
          endif
          do k = ksta,kend-1
            rmp = rmp+R3(1,1,k)*pi/(alpha*real(nx,kind=GP))
          end do
          if ( myrank.lt.nprocs-1 ) then
            rmp = rmp+R3(1,1,kend)*pi/(alpha*real(nx,kind=GP))
            call MPI_SEND(rmp,1,GC_REAL,myrank+1,1,MPI_COMM_WORLD,ierr)
          else
            rms = rmp+R3(1,1,nz)*pi/(alpha*real(nx,kind=GP)) ! total phase
          endif
        endif
      end do
      call MPI_BCAST(rms,1,GC_REAL,nprocs-1,MPI_COMM_WORLD,ierr)
      if ( myrank.eq.0 ) then
        solver%vz_(1,1,1) = solver%vz_(1,1,1) - rms/pi*alpha*real(nx,kind=GP)**3
        C3(1,1,1) = C3(1,1,1) - rms/pi*alpha
      endif
      call copy3(C4,C1)
      call fftp3d_complex_to_real(plancr,C4,R3,MPI_COMM_WORLD) ! vx
      rmp = 0.0_GP
      do i = 1,nx-1
        rmp = rmp+R3(i,1,ksta)*pi/(alpha*real(nx,kind=GP))
      end do
      rms = rmp+R3(nx,1,ksta)*pi/(alpha*real(nx,kind=GP))
      if ( myrank.eq.0 ) then
        solver%vx_(1,1,1) = solver%vx_(1,1,1) - rms/pi*alpha*real(nx,kind=GP)**3
        C1(1,1,1) = C1(1,1,1) - rms/pi*alpha
      endif
      call copy3(C5,C2)
      call fftp3d_complex_to_real(plancr,C5,R3,MPI_COMM_WORLD) ! vy
      rmp = 0.0_GP
      do j = 1,ny-1
        rmp = rmp+R3(1,j,ksta)*pi/(alpha*real(nx,kind=GP))
      end do
      rms = rmp+R3(1,ny,ksta)*pi/(alpha*real(nx,kind=GP))
      if ( myrank.eq.0 ) then
        solver%vy_(1,1,1) = solver%vy_(1,1,1) - rms/pi*alpha*real(nx,kind=GP)**3
        C2(1,1,1) = C2(1,1,1) - rms/pi*alpha
      endif
      if ( stat .gt. 0 ) then ! restart: z is read by read_z
        call solver%workspace_%free_real_tmp(R3)
        call solver%workspace_%free_real_tmp(R2)
        call solver%workspace_%free_real_tmp(R1)
        call solver%workspace_%free_complex_tmp(C6)
        call solver%workspace_%free_complex_tmp(C5)
        call solver%workspace_%free_complex_tmp(C4)
        call solver%workspace_%free_complex_tmp(C3)
        call solver%workspace_%free_complex_tmp(C2)
        call solver%workspace_%free_complex_tmp(C1)
        return
      endif
      ! Phase of the order parameter: integrates the velocity from the
      ! point (1,1,1) along z, then along x, and then along y
      call fftp3d_complex_to_real(plancr,C3,R3,MPI_COMM_WORLD) ! vz
      do ki = 0,nprocs-1
        if ( ki.eq.myrank ) then
          if ( myrank.eq.0 ) then
            R1(1,1,1) = 1.0_GP
            R2(1,1,1) = 0.0_GP
            cdump = R1(1,1,1)+im*R2(1,1,1)
          else
            call MPI_RECV(cdump,1,GC_COMPLEX,myrank-1,1,MPI_COMM_WORLD, &
                          MPI_STATUS_IGNORE,ierr)
            R1(1,1,ksta) = real(cdump)
            R2(1,1,ksta) = aimag(cdump)
          endif
          do k = ksta,kend-1
            cdump = cdump*exp(im*R3(1,1,k)*pi/(alpha*real(nx,kind=GP)))
            R1(1,1,k+1) = real(cdump)
            R2(1,1,k+1) = aimag(cdump)
          end do
          if ( myrank.lt.nprocs-1 ) then
            cdump = cdump*exp(im*R3(1,1,kend)*pi/(alpha*real(nx,kind=GP)))
            call MPI_SEND(cdump,1,GC_COMPLEX,myrank+1,1,MPI_COMM_WORLD,ierr)
          endif
        endif
      end do
      call fftp3d_complex_to_real(plancr,C1,R3,MPI_COMM_WORLD) ! vx
      do k = ksta,kend
        do i = 1,nx-1
          cdump = (R1(i,1,k)+im*R2(i,1,k))*exp(im*R3(i,1,k)*pi/(alpha*real(nx,kind=GP)))
          R1(i+1,1,k) = real(cdump)
          R2(i+1,1,k) = aimag(cdump)
        end do
      end do
      call fftp3d_complex_to_real(plancr,C2,R3,MPI_COMM_WORLD) ! vy
      do k = ksta,kend
        do i = 1,nx
          do j = 1,ny-1
            cdump = (R1(i,j,k)+im*R2(i,j,k))*exp(im*R3(i,j,k)*pi/(alpha*real(nx,kind=GP)))
            R1(i,j+1,k) = real(cdump)
            R2(i,j+1,k) = aimag(cdump)
          end do
        end do
      end do
      ! Amplitude of the density from the distance to the filaments
!$omp parallel do collapse(2) private (i,ki,ic,nc,rmt,rm1,pos,tan,dst,tmq)
      do k = ksta,kend
        do j = 1,ny
          do i = 1,nx
            rmt = 2*pi
            do ki = 1,nsc
              rm1 = 2.0_GP*pi*real(ki-1,kind=GP)/real(nsc,kind=GP)
              call knot_curve(icurve,r0,rm1,nc,pos,tan)
              do ic = 1,nc
                dst(1) = 2*pi*real(i-1,kind=GP)/real(nx,kind=GP)-pos(1,ic)
                dst(2) = 2*pi*real(j-1,kind=GP)/real(ny,kind=GP)-pos(2,ic)
                dst(3) = 2*pi*real(k-1,kind=GP)/real(nz,kind=GP)-pos(3,ic)
                rmt = min(rmt,sqrt(dst(1)**2+dst(2)**2+dst(3)**2))
              end do
            end do
            tmq = sqrt((11.0_GP*(rmt/lambda)**2/32.0_GP+11.0_GP*(rmt/lambda)**4/384.0_GP) &
                      /(1.0_GP +(rmt/lambda)**2/3.0_GP +11.0_GP*(rmt/lambda)**4/384.0_GP))
            R1(i,j,k) = R1(i,j,k)*rho*tmq
            R2(i,j,k) = R2(i,j,k)*rho*tmq
          end do
        end do
      end do
      call fftp3d_real_to_complex(planrc,R1,state(solver%ZFUNC  )%ccomp,MPI_COMM_WORLD)
      call fftp3d_real_to_complex(planrc,R2,state(solver%ZFUNC+1)%ccomp,MPI_COMM_WORLD)
      call dealias_z(solver,state)
      call solver%workspace_%free_real_tmp(R3)
      call solver%workspace_%free_real_tmp(R2)
      call solver%workspace_%free_real_tmp(R1)
      call solver%workspace_%free_complex_tmp(C6)
      call solver%workspace_%free_complex_tmp(C5)
      call solver%workspace_%free_complex_tmp(C4)
      call solver%workspace_%free_complex_tmp(C3)
      call solver%workspace_%free_complex_tmp(C2)
      call solver%workspace_%free_complex_tmp(C1)
    class default
      error stop 'IC: This solver does not support order parameter ICs'
    end select
  end subroutine init_knotz


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Real Gaussian wave function for solvers with a cigar
  !! trapping potential V = V0 (x^2+y^2) (use with the forcing
  !! method 'cyltrap_vq'). Uses rho0 and V0 of the solver.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_trapz(this,solver,state)
    use gstate_mod
    use grid
    use kes
    use ali
    use var
    use fft
    use commtypes
    use mpivars
    implicit none

    class    (icTrap_z), intent(in)            :: this
    class(EquationBase), intent(inout)         :: solver
    type   (GStateComp), intent(inout)         :: state(:)
    real(kind=GP), pointer, dimension(:,:,:)   :: R1
    real(kind=GP)                              :: rmp,rmq,rms,rmt
    integer                                    :: i,j,k
    logical                                    :: bret

    select type (solver)
    class is (QuantumBase)
      if ( solver%V0_ .le. 0.0_GP ) then
        if (myrank.eq.0) error stop 'trap_z requires a trapping potential (V0 > 0)'
      endif
      call solver%workspace_%get_real_tmp(R1,bret)
      rmp = .5_GP*sqrt(solver%V0_/solver%alpha_)
      rmq = sqrt(solver%rho0_)*(sqrt(solver%V0_/solver%alpha_)/pi)**(3./4)
!$omp parallel do collapse(2) private (i,rms,rmt)
      do k = ksta,kend
        do j = 1,ny
          do i = 1,nx
            rms = (pi*(2*real(j-1,kind=GP)/real(ny,kind=GP)-1.0_GP))**2 ! y^2
            rmt = (pi*(2*real(i-1,kind=GP)/real(nx,kind=GP)-1.0_GP))**2 ! x^2
            R1(i,j,k) = rmq*exp(-rmp*(rms+rmt))
          end do
        end do
      end do
      call fftp3d_real_to_complex(planrc,R1,state(solver%ZFUNC)%ccomp,MPI_COMM_WORLD)
!$omp parallel do collapse(2) private (k)
      do i = ista,iend
        do j = 1,ny
          do k = 1,nz
            if (kn2(k,j,i).gt.kmax) state(solver%ZFUNC)%ccomp(k,j,i) = 0.0_GP
            state(solver%ZFUNC+1)%ccomp(k,j,i) = 0.0_GP
          end do
        end do
      end do
      call solver%workspace_%free_real_tmp(R1)
    class default
      error stop 'IC: This solver does not support order parameter ICs'
    end select
  end subroutine init_trapz

end module ic_quantum
