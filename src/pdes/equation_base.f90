! ===================================================================
! NAME       : equation_base.f90
! DESCRIPTION: Forms base class for all PDEs
!
! INPUT FILE : All PDEs look for a "&stepper" block that
!              governs configuration of the time stepping
!              scheme. Folliwing are the parameters allowed 
!              in the configuration block:
!                sname : Name of scheme. check stepper_factory_mod 
!                        for valid names
!                norder: Order of scheme
!                nstage: Number of stages (if doing GEXRK)
!                itype : Sub-type of RK scheme is using GEXRK

! DATE       : 11/29/25 (DLR)
! ===================================================================

module equationbase_mod
  use class_GWorkspace3D
  use gstate_mod
  use iovar

  implicit none

  ! ================= Base class for all PDEs =======================
  ! Define an abstract base class
  type, abstract :: EquationBase
      type(GWorkspace), pointer     :: workspace_
      type    (ioplan), pointer     :: planio_
      integer                       :: myrank_   ! MPI rank
      integer                       :: nprocs_   ! MPI rank 
      integer                       :: order_    ! default integration order
      character(len=8), allocatable :: sstate_(:)! state member nanes
      character(len=128)            :: infile_   ! config file name
!     type(StepperBase),allocatable :: stepper_  ! time stepping object
    contains
      procedure(Solver_ctor_interface), deferred :: Solver_ctor ! Constructor
      procedure(init_interface),        deferred :: init        ! init method
      procedure(dudt_interface),        deferred :: dudt        ! RHS method
      procedure(global_interface),      deferred :: global      ! Global qtys
      procedure(spectra_interface),     deferred :: spectra     ! Spectra
      procedure(state_size_interface),  deferred :: state_size  ! Number of states
!     procedure, public                          :: timestep_f, timestep_fp
      procedure, public                          :: write_states
  end type EquationBase

  type, abstract, extends(EquationBase) :: VelocityBase
      integer :: VELOCITY    ! start of velocity sector
      integer :: PASSIVE     ! start of scalar sector
      integer :: numpassive_ ! # passive scalars
      integer :: nd_         ! problem dimension
      integer :: nc_         ! # vector field components
    contains
      procedure, public                          :: rhs_passive
  end type VelocityBase

  type, abstract, extends(VelocityBase) :: ActiveScalarBase
      integer :: ACTIVESC    ! start of active scalar sector
  end type ActiveScalarBase

  type, abstract, extends(VelocityBase) :: MagneticBase
      integer :: MAGNETIC    ! start of magnetic sector
  end type MagneticBase

  type, abstract, extends(EquationBase) :: QuantumBase
      integer :: ZFUNC       ! start of wavefunction sector
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
       class(EquationBase), intent   (in)         :: this
       real      (kind=GP), intent   (in)         :: time, dt
       type   (GStateComp), intent(inout), target :: uin(:),uf(:)
       type   (GStateComp), intent(inout)         :: dudt(:) 
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
  !! Base constructor: 
  !! must be called from each child constructor
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine Solver_base_ctor(this, infile, workspace, plan)
    use class_GWorkspace3D
!   use stepper_factory_mod
    use iovar
    class(EquationBase), intent(inout)         :: this
    type(GWorkspace)   , intent(inout), target :: workspace
    type(ioplan)       , intent(inout), target :: plan
    character(len=*)   , intent   (in)         :: infile

    this%workspace_ => workspace;
!    this%stepper_ = build_stepper_from_file(infile, workspace)
 
  end subroutine Solver_base_ctor

!!$  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$  !! Concrete method to take one time step 
!!$  !! using configured time stepping object.
!!$  !! This method applies only to field evolution.
!!$  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$  subroutine timestep_f(this, time, uin, uf, dt, uout)
!!$!$  use threads
!!$    implicit none
!!$
!!$    class (EquationBase), intent   (in) :: this
!!$    type    (GStateComp), intent(inout) :: uin(:), uf(:), uout(:)
!!$    real       (kind=GP), intent   (in) :: time, dt
!!$
!!$    if ( .not. allocated(this%stepper_) ) then
!!$      stop 'EquationBase::timestep: time stepper object not allocated'
!!$    endif
!!$
!!$    this%stepper(time, uin, uf, dt, uout)
!!$
!!$  end subroutine timestep_f

!!$  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$  !! Concrete method to take one time step 
!!$  !! using configured time stepping object.
!!$  !! This method applies to field and
!!$  !! particle evolution.
!!$  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$  subroutine timestep_fp(this, time, uin, upin, uf, dt, uout, upout)
!!$!$  use threads
!!$    implicit none
!!$
!!$    class (EquationBase), intent   (in) :: this
!!$    type    (GStateComp), intent(inout) :: uin (:), uf(:), uout (:)
!!$    type   (GPStateComp), intent(inout) :: upin(:)       , upout(:)
!!$    real       (kind=GP), intent   (in) :: time, dt
!!$
!!$    if ( .not. allocated(this%stepper_) ) then
!!$      stop 'EquationBase::timestep: time stepper object not allocated'
!!$    endif
!!$
!!$    this%stepper(time, uin, upin, uf, dt, uout, upout)
!!$
!!$  end subroutine timestep_fp

  
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
!$  use threads
    implicit none

    class(VelocityBase), intent(in)   :: this
    real     (kind=GP), intent   (in) :: kappa(:)
    type      (GStateComp), intent(inout) :: uin(:)
    type      (GStateComp), intent   (in) :: uf(:)
    type      (GStateComp), intent(inout) :: dudt(:) 
    logical                           :: bret
    integer                           :: i,j,k,n
    complex  (kind=GP), pointer       :: adve(:,:,:),lapl(:,:,:)

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
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
        do i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
          do j = 1,ny
            do k = 1,nz
              if ((kn2(k,j,i).le.kmax).and.(kn2(k,j,i).ge.tiny)) then
                dudt(n)%ccomp(k,j,i) =                  &
                  kappa(n-this%PASSIVE+1)*lapl(k,j,i) + &
                  adve(k,j,i) + uf(n)%ccomp(k,j,i)
              else if ((kn2(k,j,i).gt.kmax).or.(kn2(k,j,i).lt.tiny)) then
                  dudt(n)%ccomp(k,j,i) = 0.0_GP
              endif
            enddo
          enddo
        enddo
      enddo ! end, loop over all scalars
    endif
    call this%workspace_%free_complex_tmp(adve)
    call this%workspace_%free_complex_tmp(lapl)
  end subroutine rhs_passive

  
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
!$  use threads
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
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      do i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
        do j = 1,ny
          do k = 1,nz
            C1(k,j,i) = uin(nc)%ccomp(k,j,i)*rmp
          end do
        end do
      end do
      call fftp3d_complex_to_real(plancr,C1,R1,MPI_COMM_WORLD)
      call io_write(1,odir,trim(this%sstate_(nc)),ext,planio,R1)
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
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
        do i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
          do j = 1,ny
            do k = 1,nz
              C1(k,j,i) = uin(this%VELOCITY  )%ccomp(k,j,i)*rmp
              C2(k,j,i) = uin(this%VELOCITY+1)%ccomp(k,j,i)*rmp
              C3(k,j,i) = uin(this%VELOCITY+1)%ccomp(k,j,i)*rmp
            end do
          end do
        end do
        call rotor3(C2,C3,C4,1) 
        call rotor3(C1,C3,C5,2)
        call rotor3(C1,C2,C6,3)
        call fftp3d_complex_to_real(plancr,C4,R1,MPI_COMM_WORLD)
        call fftp3d_complex_to_real(plancr,C5,R2,MPI_COMM_WORLD)
        call fftp3d_complex_to_real(plancr,C6,R3,MPI_COMM_WORLD)
        call io_write(1,odir,'wx',ext,planio,R1)
        call io_write(1,odir,'wy',ext,planio,R2)
        call io_write(1,odir,'wz',ext,planio,R3)
        select type (this)
        class is (MagneticBase)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
          do i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
            do j = 1,ny
              do k = 1,nz
                C1(k,j,i) = uin(this%MAGNETIC  )%ccomp(k,j,i)*rmp
                C2(k,j,i) = uin(this%MAGNETIC+1)%ccomp(k,j,i)*rmp
                C3(k,j,i) = uin(this%MAGNETIC+1)%ccomp(k,j,i)*rmp
              end do
            end do
          end do
          call rotor3(C2,C3,C4,1) 
          call rotor3(C1,C3,C5,2)
          call rotor3(C1,C2,C6,3)
          call fftp3d_complex_to_real(plancr,C4,R1,MPI_COMM_WORLD)
          call fftp3d_complex_to_real(plancr,C5,R2,MPI_COMM_WORLD)
          call fftp3d_complex_to_real(plancr,C6,R3,MPI_COMM_WORLD)
          call io_write(1,odir,'bx',ext,planio,R1)
          call io_write(1,odir,'by',ext,planio,R2)
          call io_write(1,odir,'bz',ext,planio,R3)
          if ( outs .eq. 2 ) then
            call laplak3(C1,C4)
            call laplak3(C2,C5)
            call laplak3(C3,C6)
            call fftp3d_complex_to_real(plancr,C4,R1,MPI_COMM_WORLD)
            call fftp3d_complex_to_real(plancr,C5,R2,MPI_COMM_WORLD)
            call fftp3d_complex_to_real(plancr,C6,R3,MPI_COMM_WORLD)
            call io_write(1,odir,'jx',ext,planio,-R1)
            call io_write(1,odir,'jy',ext,planio,-R2)
            call io_write(1,odir,'jz',ext,planio,-R3)
          endif
        end select
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


