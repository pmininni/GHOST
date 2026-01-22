! ===================================================================
! NAME       : equation_base.fpp
! DESCRIPTION: Forms base class for all PDEs
! DATE       : 11/29/25 (DLR)
! ===================================================================

module equationbase_mod
  USE class_GWorkspace3D
  USE gstate_mod
  USE iovar

  IMPLICIT NONE

  ! ================= Base class for all PDEs =======================
  ! Define an abstract base class
  type, abstract :: EquationBase
      type(GWorkspace), pointer     :: workspace_
      type    (ioplan), pointer     :: planio_
      integer                       :: myrank_   ! MPI rank
      integer                       :: nprocs_   ! MPI rank 
      character(len=8), allocatable :: sstate_(:)
      character(len=128)            :: infile_
    contains
      procedure(Solver_ctor_interface), deferred :: Solver_ctor ! Constructor
      procedure(init_interface),        deferred :: init        ! init method
      procedure(dudt_interface),        deferred :: dudt        ! RHS method
      procedure(global_interface),      deferred :: global      ! Global qtys
      procedure(spectra_interface),     deferred :: spectra     ! Spectra
      procedure(state_size_interface),  deferred :: state_size  ! Number of states
      procedure, public                          :: timestep
      procedure, public                          :: read_states
      procedure, public                          :: write_states
  end type EquationBase

  type, abstract, extends(EquationBase) :: VelocityBase
      integer :: VELOCITY    ! start of velocity sector
      integer :: PASSIVE     ! start of scalar sector
      integer :: numpassive_ ! # passive scalars
      integer :: nd_         ! problem dimension
      integer :: nc_         ! # vector field components
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
       USE class_GWorkspace3D
       USE iovar
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
       USE gstate_mod
       import :: EquationBase
       class(EquationBase), intent   (in)         :: this
       real      (kind=GP), intent   (in)         :: time, dt
       type       (GState), intent(inout), target :: uin(:),uf(:)
       type       (GState), intent(inout)         :: dudt(:) 
     end subroutine dudt_interface

     subroutine global_interface(this, uin, uf, t) 
       USE gstate_mod
       USE fprecision
       USE status
       import :: EquationBase
       class(EquationBase), intent(in)            :: this
       type       (GState), intent(in), target    :: uin(:),uf(:)
       integer            , intent(in)            :: t
     end subroutine global_interface
       
     subroutine spectra_interface(this, uin) 
       USE gstate_mod
       USE fprecision
       USE filefmt
       USE status
       import :: EquationBase
       class(EquationBase), intent(in)            :: this
       type       (GState), intent(in), target    :: uin(:)
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
  !! Concrete method to take one time step using Runge-Kutta
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine timestep(this, time, uin, uf, dt, uout)
    use order
!$  use threads
    implicit none

    class (EquationBase), intent   (in) :: this
    type        (GState), intent(inout) :: uin(:), uf(:), uout(:)
    real       (kind=GP), intent   (in) :: time, dt
    real       (kind=GP)                :: eff_dt
    integer                             :: i,j,k,o,state_size,nc

    state_size = this%state_size()
    do o = ord,1,-1
      eff_dt = dt/real(o,kind=GP)
      CALL this%dudt(time, uout, uf, eff_dt, uout)
      do nc = 1,state_size
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
        do i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
          do j = 1,ny
            do k = 1,nz
              uout(nc)%ccomp(k,j,i) = uin(nc)%ccomp(k,j,i) + &
                              eff_dt*uout(nc)%ccomp(k,j,i)
            end do
          end do
        end do
      end do
    end do

  end subroutine timestep

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Concrete method to compute RHS for all passive scalars
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine rhs_passive(this, uin, f, &
               kappa, ivstart, numv, isstart, npassive, dudt)
  !------------------------------------------------------------------
  ! Parameters
  !     uin     : current full state
  !     f       : forces for each state comp
  !     kappa   : diffusivities (must be npassive of these)
  !     ivstart : start index of velocity sector in uin
  !     numv    : number of velocity comps to use for advection
  !     sstart  : start index of scalar sector
  !     npassive: number of passive scalars to advect
  !     dudt    : computed RHS for each scalar, 1-npassive
  !******************************************************************
    use pseudospec_scalar
    use fprecision
    use ali
    use kes
    use var
    use grid
    use mpivars
    use gstate_mod
!$  use threads
    implicit none

    class (EquationBase), intent(in) :: this
    integer          , intent   (in) :: isstart, npassive
    integer          , intent   (in) :: ivstart, numv
    real    (kind=GP), intent   (in) :: kappa(:)
    type     (GState), intent(inout) :: uin(:)
    type     (GState), intent(inout) :: f(:)
    type (GWorkspace), pointer       :: workspace
    type     (GState), intent  (out) :: dudt(:) 
    logical                          :: bret
    integer                          :: i,j,k,n
    complex, pointer                 :: ctmp(:,:,:)
 
    CALL workspace%get_complex_tmp(ctmp,bret)
    if ( .not.bret ) then
      stop 'EquationBase::rhs_passive: workspace get failure'
    endif 

    if ( numv .eq. 3 ) then ! 3d advection
      do n = 1, npassive
         call advect3(uin(ivstart  )%ccomp,uin(ivstart+1  )%ccomp, &
                      uin(ivstart+2)%ccomp,uin(isstart+n-1)%ccomp, &
                      ctmp)
        call laplak3(uin(isstart+n-1)%ccomp,uin(isstart+n-1)%ccomp)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
        do i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
        do j = 1,ny
        do k = 1,nz
          if ((kn2(k,j,i).le.kmax).and.(kn2(k,j,i).ge.tiny)) then
            dudt(n)%ccomp(k,j,i) = &
                kappa(n)*uin(isstart+n-1)%ccomp(k,j,i)+ctmp(k,j,i) &
                +f(isstart)%ccomp(k,j,i)
          else if (kn2(k,j,i).gt.kmax) then
            dudt(n)%ccomp(k,j,i) = 0.0_GP
          else if (kn2(k,j,i).lt.tiny) then
            dudt(n)%ccomp(k,j,i) = 0.0_GP
          endif
        enddo
        enddo
        enddo
      enddo ! end, loop over scalars
!!$    else ! 2d advection:
!!$      do n = 1, npassive
!!$        call advect3(uin(ivstart),uin(ivstart+1),uin(ivstart+2),&
!!$                                 uin(isstart+n-1),ctmp)
!!$        call laplak3(uin(isstart+n-1),uin(isstart+n-1))
!!$!$omp parallel do if (iend-ista.ge.nth) private (j,k)
!!$        do i = ista,iend
!!$!$omp parallel do if (iend-ista.lt.nth) private (k)
!!$        do j = 1,ny
!!$        do k = 1,nz
!!$        if ((kn2(k,j,i).le.kmax).and.(kn2(k,j,i).ge.tiny)) then
!!$          dudt(n+isstart)%ccomp(k,j,i) = &
!!$              kappa(n)*uin(isstart+n-1)%ccomp(k,j,i)+ctmp(k,j,i) &
!!$              +f(isstart)%ccomp(k,j,i)
!!$        else if (kn2(k,j,i).gt.kmax) then
!!$          dudt(n+isstart)%ccomp(k,j,i) = 0.0_GP
!!$        else if (kn2(k,j,i).lt.tiny) then
!!$          dudt(n+isstart)%ccomp(k,j,i) = 0.0_GP
!!$        endif
!!$        enddo
!!$        enddo
!!$        enddo
!!$      enddo ! end, loop over scalars
    endif

    CALL workspace%free_complex_tmp(ctmp)

  end subroutine rhs_passive

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Concrete method to write field states
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine write_states(this, uin, planio)
    use fprecision
    use grid
    use iovar
    use status
    use filefmt
    use fft
    use commtypes
!$  use threads
    implicit none

    class (EquationBase), intent   (in)             :: this
    type        (GState), intent(inout)             :: uin(:)
    type        (ioplan), intent   (in)             :: planio
    complex    (kind=GP), pointer, dimension(:,:,:) :: C1
    real       (kind=GP), pointer, dimension(:,:,:) :: R1
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
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
        DO j = 1,ny
          DO k = 1,nz
            C1(k,j,i) = uin(nc)%ccomp(k,j,i)*rmp
          END DO
        END DO
      END DO
      CALL fftp3d_complex_to_real(plancr,C1,R1,MPI_COMM_WORLD)
      CALL io_write(1,odir,trim(this%sstate_(nc)),ext,planio,R1)
    end do
    call this%workspace_%free_complex_tmp(C1)
    call this%workspace_%free_real_tmp   (R1)
  end subroutine write_states


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Concrete method to read field states
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine read_states(this, uin, planio)
    use fprecision
    use grid
    use iovar
    use status
    use filefmt
    use fft
    use commtypes
!$  use threads
    implicit none

    class (EquationBase), intent   (in)             :: this
    type        (GState), intent(inout)             :: uin(:)
    type        (ioplan), intent   (in)             :: planio
    complex    (kind=GP), pointer, dimension(:,:,:) :: C1
    real       (kind=GP), pointer, dimension(:,:,:) :: R1
    integer                          :: i,j,k,o,state_size,nc
    logical                          :: bret

    WRITE(ext, fmtext) tind
    call this%workspace_%get_complex_tmp(C1,bret)
    call this%workspace_%get_real_tmp(   R1,bret)
    state_size = this%state_size()
    do nc = 1,state_size
      call io_read(1,idir,trim(this%sstate_(nc)),ext,planio,R1)
      call fftp3d_real_to_complex(planrc,R1,uin(nc)%ccomp(k,j,i), &
                                  MPI_COMM_WORLD)
   end do
   call this%workspace_%free_complex_tmp(C1)
   call this%workspace_%free_real_tmp(   R1)
   
  end subroutine read_states
  
end module equationbase_mod


