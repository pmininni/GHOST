! ===================================================================
! NAME       : equation_base.f90
! DESCRIPTION: Forms base class for all PDEs
!
! DATE       : 11/29/25 (DLR)
! ===================================================================

module equationbase_mod
  use class_GWorkspace3D
  use gstate_mod
  use iovar
  implicit none

  private :: Solver_ctor_interface, init_interface
  private :: dudt_interface       , global_interface
  private :: spectra_interface    , state_size_interface

  ! ================= Base class for all PDEs =======================
  ! Define an abstract base class
  type, abstract :: EquationBase
      type(GWorkspace), pointer     :: workspace_ => null()
      type    (ioplan), pointer     :: planio_    => null()
      integer                       :: myrank_     ! MPI rank
      integer                       :: nprocs_     ! MPI rank 
      character(len=8), allocatable :: sstate_(:)  ! state member nanes
      character(len=128)            :: infile_     ! config file name
      character(len=128)            :: odir_,idir_ ! internal class I/O directories
      character(len=128)            :: todir_      ! internal class TXT output dir
    contains
      procedure(Solver_ctor_interface), deferred :: Solver_ctor ! Constructor
      procedure(init_interface),        deferred :: init        ! init method
      procedure(dudt_interface),        deferred :: dudt        ! RHS method
      procedure(global_interface),      deferred :: global      ! Global qtys
      procedure(spectra_interface),     deferred :: spectra     ! Spectra
      procedure(state_size_interface),  deferred :: state_size  ! Number of states
      procedure, public                          :: write_states
      procedure, public                          :: sync_device ! Aux. arrays to the device
  end type EquationBase

  type, abstract, extends(EquationBase) :: VelocityBase
      integer :: VELOCITY    ! start of velocity sector
      integer :: PASSIVE     ! start of scalar sector
      integer :: numpassive_ ! # passive scalars
      integer :: nd_         ! problem dimension
      integer :: nc_         ! # vector field components
    contains
      procedure, public                 :: rhs_passive
  end type VelocityBase

  type, abstract, extends(VelocityBase) :: ActiveScalarBase
      integer :: ACTIVESC     ! start of active scalar sector
      integer :: numactivesc_ ! # active scalars
  end type ActiveScalarBase

  type, abstract, extends(VelocityBase) :: CompVelocityBase
      integer :: DENSITY     ! start of mass density sector
  end type CompVelocityBase

  type, abstract, extends(VelocityBase) :: MagneticBase
      integer :: MAGNETIC    ! start of magnetic sector
  end type MagneticBase

  type, abstract, extends(MagneticBase) :: CompMagneticBase
      integer :: DENSITY     ! start of mass density sector
  end type CompMagneticBase

  type, abstract, extends(EquationBase) :: QuantumBase
      integer :: ZFUNC       ! start of wavefunction sector
      ! Constants of the GPE-type equations (set by the solvers, read by
      ! the initial conditions and the forcing): alpha = c.xi/sqrt(2),
      ! omegag = c/(xi.sqrt(2)), beta = omegag/rho0, with c the speed of
      ! sound, xi the coherence length and rho0 the equilibrium density
      real(kind=GP) :: alpha_  = 0.0_GP
      real(kind=GP) :: beta_   = 0.0_GP
      real(kind=GP) :: omegag_ = 0.0_GP
      real(kind=GP) :: rho0_   = 1.0_GP
      real(kind=GP) :: V0_     = 0.0_GP ! amplitude of the trapping potential
      real(kind=GP) :: omegaz_ = 0.0_GP ! rotation rate (rotating frame)
      logical       :: hasadv_ = .false. ! advective velocity in use
      logical       :: haspot_ = .false. ! external potential in use
      logical       :: dorot_  = .false. ! rotating frame
      ! Auxiliary fields, constant in time, with device copies: the
      ! advective velocity (Fourier space, unnormalized), |v|^2/(4 alpha
      ! beta) (real space, unnormalized), the external potential divided
      ! by beta, and the linear ramps omegaz.x and omegaz.y for the
      ! angular momentum operator (real space, unnormalized). They are
      ! set by the initial conditions (velocity) and the forcing
      ! (potential), and copied to the device by sync_device.
      complex(kind=GP), allocatable, dimension(:,:,:) :: vx_,vy_,vz_
      real   (kind=GP), allocatable, dimension(:,:,:) :: vsq_,vpot_
      real   (kind=GP), allocatable, dimension(:,:,:) :: vlinx_,vliny_
    contains
      procedure, public :: alloc_aux     => quantum_alloc_aux     ! Allocates aux. arrays
      procedure, public :: free_aux      => quantum_free_aux      ! Frees aux. arrays
      procedure, public :: aux_to_device => quantum_aux_to_device ! Copies them to the device
  end type QuantumBase
  
  abstract interface
     subroutine Solver_ctor_interface(this, infile, workspace, plan)
       use class_GWorkspace3D
       use iovar
       import :: EquationBase
       class(EquationBase), intent(inout)         :: this
       type(GWorkspace)   , intent(inout), target :: workspace
       type(ioplan)       , intent(inout), target :: plan
       character(len=*)   , intent   (in)         :: infile  
     end subroutine Solver_ctor_interface

     subroutine init_interface(this) 
       import :: EquationBase
       class (EquationBase), intent (inout) :: this
     end subroutine init_interface

     subroutine dudt_interface(this, time, uin, uf, dt, dudt) 
       use gstate_mod
       import :: EquationBase
       class(EquationBase), intent(inout)         :: this
       real      (kind=GP), intent   (in)         :: time, dt
       type   (GStateComp), intent(inout), target :: uin(:),uf(:)
       type   (GStateComp), intent(inout), target :: dudt(:) 
     end subroutine dudt_interface

     subroutine global_interface(this, uin, uf, t) 
       use gstate_mod
       use fprecision
       use status
       import :: EquationBase
       class(EquationBase), intent(in)            :: this
       type   (GStateComp), intent(in), target    :: uin(:),uf(:)
       integer            , intent(in)            :: t
     end subroutine global_interface
       
     subroutine spectra_interface(this, uin) 
       use gstate_mod
       use fprecision
       use filefmt
       use status
       import :: EquationBase
       class(EquationBase), intent(in)            :: this
       type   (GStateComp), intent(in), target    :: uin(:)
     end subroutine spectra_interface

     function state_size_interface(this) result(num)
       import :: EquationBase
       class(EquationBase), intent   (in)         :: this
       integer                                    :: num
     end function state_size_interface
  end interface

CONTAINS

  ! ===================================================================
  ! Concrete methods inherited by all solvers
  ! ===================================================================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Concrete method to compute RHS for all passive scalars
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine rhs_passive(this, uin, uf, kappa, dudt)
  !----------------------------------------------------------
  ! Parameters
  !   uin  : current full state
  !   uf   : forces for each state comp
  !   kappa: diffusivities (must be npassive of these)
  !   dudt : computed RHS for full state
  !**********************************************************
    use pseudospec_scalar
    use ali
    use kes
    use var
    use grid
    use mpivars
    use gstate_mod
    implicit none

    class(VelocityBase), intent(in)    :: this
    real      (kind=GP), intent   (in) :: kappa(:)
    type   (GStateComp), intent(inout) :: uin(:)
    type   (GStateComp), intent   (in) :: uf(:)
    type   (GStateComp), intent(inout) :: dudt(:) 
    logical                            :: bret
    integer                            :: n
    complex  (kind=GP), pointer        :: adve(:,:,:),lapl(:,:,:)

    if ( this%numpassive_ .eq. 0 ) return
    call this%workspace_%get_complex_tmp(adve,bret)
    call this%workspace_%get_complex_tmp(lapl,bret)
    if ( .not.bret ) then
      stop 'EquationBase::rhs_passive: workspace get failure'
    endif
    if ( this%nc_ .eq. 3 ) then ! 3d advection
      do n = this%PASSIVE, this%PASSIVE+this%numpassive_-1
        call advect3(uin(this%VELOCITY  )%ccomp, &
                     uin(this%VELOCITY+1)%ccomp, &
                     uin(this%VELOCITY+2)%ccomp, &
                     uin(n)%ccomp,adve)
        call laplak3(uin(n)%ccomp,lapl)
        call rhs_scalar3(lapl,adve,uf(n)%ccomp,kappa(n-this%PASSIVE+1),dudt(n)%ccomp)
      enddo ! end, loop over all scalars
    endif
    call this%workspace_%free_complex_tmp(adve)
    call this%workspace_%free_complex_tmp(lapl)
  end subroutine rhs_passive

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Copies the auxiliary arrays of the solver (if any) to the
  !! device, after they are set by the initial conditions and the
  !! forcing. Default: the solver has no auxiliary arrays.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine sync_device(this)
    class (EquationBase), intent(inout) :: this
    return
  end subroutine sync_device


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Concrete methods of the quantum solvers: allocation (with
  !! device copies, zero by default), deallocation, and copy to
  !! the device of the auxiliary arrays shared by the quantum
  !! solvers (advective velocity, |v|^2, potential, and linear
  !! ramps for the rotation)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine quantum_alloc_aux(this)
    use gmem
    use grid
    use mpivars
    class (QuantumBase), intent(inout) :: this
    call galloc(this%vx_   ,nz,ny,ista,iend)
    call galloc(this%vy_   ,nz,ny,ista,iend)
    call galloc(this%vz_   ,nz,ny,ista,iend)
    call galloc(this%vsq_  ,nx,ny,ksta,kend)
    call galloc(this%vpot_ ,nx,ny,ksta,kend)
    call galloc(this%vlinx_,nx,ny,ksta,kend)
    call galloc(this%vliny_,nx,ny,ksta,kend)
    this%vx_    = 0.0_GP
    this%vy_    = 0.0_GP
    this%vz_    = 0.0_GP
    this%vsq_   = 0.0_GP
    this%vpot_  = 0.0_GP
    this%vlinx_ = 0.0_GP
    this%vliny_ = 0.0_GP
  end subroutine quantum_alloc_aux

  subroutine quantum_free_aux(this)
    use gmem
    class (QuantumBase), intent(inout) :: this
    call gfree(this%vx_)
    call gfree(this%vy_)
    call gfree(this%vz_)
    call gfree(this%vsq_)
    call gfree(this%vpot_)
    call gfree(this%vlinx_)
    call gfree(this%vliny_)
  end subroutine quantum_free_aux

  subroutine quantum_aux_to_device(this)
    use gmem
    class (QuantumBase), intent(inout) :: this
    call gupdate_to(this%vx_)
    call gupdate_to(this%vy_)
    call gupdate_to(this%vz_)
    call gupdate_to(this%vsq_)
    call gupdate_to(this%vpot_)
    call gupdate_to(this%vlinx_)
    call gupdate_to(this%vliny_)
  end subroutine quantum_aux_to_device


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Concrete method to write field states
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine write_states(this, uin, planio)
    use grid
    use iovar
    use status
    use filefmt
    use fft
    use commtypes
    use pseudospec_fluid
    implicit none

    class (EquationBase), intent   (in)             :: this
    type    (GStateComp), intent(inout)             :: uin(:)
    type        (ioplan), intent   (in)             :: planio
    complex    (kind=GP), pointer, dimension(:,:,:) :: C1,C2,C3
    complex    (kind=GP), pointer, dimension(:,:,:) :: C4,C5,C6
    real       (kind=GP), pointer, dimension(:,:,:) :: R1,R2,R3
    real       (kind=GP)                            :: rmp
    integer                          :: i,j,k,o,state_size,nc
    logical                          :: bret

    WRITE(ext, fmtext) tind
    call this%workspace_%get_complex_tmp(C1,bret)
    call this%workspace_%get_real_tmp   (R1,bret)
    state_size = this%state_size()
    do nc = 1,state_size
      rmp = 1.0_GP/(real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))
!$omp parallel do collapse(2) private (k)
      do i = ista,iend
         do j = 1,ny
            do concurrent (k=1:nz)
              C1(k,j,i) = uin(nc)%ccomp(k,j,i)*rmp
            end do
         end do
      end do
      call fftp3d_complex_to_real(plancr,C1,R1,MPI_COMM_WORLD)
      call io_write(1,this%odir_,trim(this%sstate_(nc)),ext,planio,R1)
    end do
    if ( outs .ge. 1) then
      call this%workspace_%get_complex_tmp(C2,bret)
      call this%workspace_%get_complex_tmp(C3,bret)
      call this%workspace_%get_complex_tmp(C4,bret)
      call this%workspace_%get_complex_tmp(C5,bret)
      call this%workspace_%get_complex_tmp(C6,bret)
      call this%workspace_%get_real_tmp   (R2,bret)
      call this%workspace_%get_real_tmp   (R3,bret)
      select type (this)
      class is (VelocityBase)
!$omp parallel do collapse(2) private (k)
        do i = ista,iend
           do j = 1,ny
              do concurrent (k=1:nz)
                C1(k,j,i) = uin(this%VELOCITY  )%ccomp(k,j,i)*rmp
                C2(k,j,i) = uin(this%VELOCITY+1)%ccomp(k,j,i)*rmp
                C3(k,j,i) = uin(this%VELOCITY+2)%ccomp(k,j,i)*rmp
              end do
           end do
        end do
        call rotor3(C2,C3,C4,1) 
        call rotor3(C1,C3,C5,2)
        call rotor3(C1,C2,C6,3)
        call fftp3d_complex_to_real(plancr,C4,R1,MPI_COMM_WORLD)
        call fftp3d_complex_to_real(plancr,C5,R2,MPI_COMM_WORLD)
        call fftp3d_complex_to_real(plancr,C6,R3,MPI_COMM_WORLD)
        call io_write(1,this%odir_,'wx',ext,planio,R1)
        call io_write(1,this%odir_,'wy',ext,planio,R2)
        call io_write(1,this%odir_,'wz',ext,planio,R3)
        select type (this)
        class is (MagneticBase)
!$omp parallel do collapse(2) private (k)
          do i = ista,iend
             do j = 1,ny
                do concurrent (k=1:nz)
                  C1(k,j,i) = uin(this%MAGNETIC  )%ccomp(k,j,i)*rmp
                  C2(k,j,i) = uin(this%MAGNETIC+1)%ccomp(k,j,i)*rmp
                  C3(k,j,i) = uin(this%MAGNETIC+2)%ccomp(k,j,i)*rmp
                end do
             end do
          end do
          call rotor3(C2,C3,C4,1) 
          call rotor3(C1,C3,C5,2)
          call rotor3(C1,C2,C6,3)
          call fftp3d_complex_to_real(plancr,C4,R1,MPI_COMM_WORLD)
          call fftp3d_complex_to_real(plancr,C5,R2,MPI_COMM_WORLD)
          call fftp3d_complex_to_real(plancr,C6,R3,MPI_COMM_WORLD)
          call io_write(1,this%odir_,'bx',ext,planio,R1)
          call io_write(1,this%odir_,'by',ext,planio,R2)
          call io_write(1,this%odir_,'bz',ext,planio,R3)
          if ( outs .eq. 2 ) then
            call laplak3(C1,C4)
            call laplak3(C2,C5)
            call laplak3(C3,C6)
            call fftp3d_complex_to_real(plancr,C4,R1,MPI_COMM_WORLD)
            call fftp3d_complex_to_real(plancr,C5,R2,MPI_COMM_WORLD)
            call fftp3d_complex_to_real(plancr,C6,R3,MPI_COMM_WORLD)
            call io_write(1,this%odir_,'jx',ext,planio,-R1)
            call io_write(1,this%odir_,'jy',ext,planio,-R2)
            call io_write(1,this%odir_,'jz',ext,planio,-R3)
          endif
        end select
      class is (QuantumBase) ! density of the order parameter, |z|^2
!$omp parallel do collapse(2) private (k)
        do i = ista,iend
           do j = 1,ny
              do concurrent (k=1:nz)
                C1(k,j,i) = uin(this%ZFUNC  )%ccomp(k,j,i)*rmp
                C2(k,j,i) = uin(this%ZFUNC+1)%ccomp(k,j,i)*rmp
              end do
           end do
        end do
        call fftp3d_complex_to_real(plancr,C1,R1,MPI_COMM_WORLD)
        call fftp3d_complex_to_real(plancr,C2,R2,MPI_COMM_WORLD)
!$omp parallel do collapse(2) private (i)
        do k = ksta,kend
           do j = 1,ny
              do concurrent (i=1:nx)
                R1(i,j,k) = R1(i,j,k)**2+R2(i,j,k)**2
              end do
           end do
        end do
        call io_write(1,this%odir_,'rho',ext,planio,R1)
      end select 
      call this%workspace_%free_complex_tmp(C2)
      call this%workspace_%free_complex_tmp(C3)
      call this%workspace_%free_complex_tmp(C4)
      call this%workspace_%free_complex_tmp(C5)
      call this%workspace_%free_complex_tmp(C6)
      call this%workspace_%free_real_tmp   (R2)
      call this%workspace_%free_real_tmp   (R3)
    endif
    call this%workspace_%free_complex_tmp(C1)
    call this%workspace_%free_real_tmp   (R1)
  end subroutine write_states

end module equationbase_mod
