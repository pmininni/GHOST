!=================================================================
! GPState: particle state components (one rank-1 array per state
! component, e.g. the x, y, z positions). In offload builds every
! component has a device copy created with the state, and the
! device copy is the one the time step works on; the host copy is
! refreshed with GPState_update_from before the particle I/O.
!
! The stepper kernels (GPState_axpy, GPState_upd) update all the
! components of a state on the device while gdev_active is set.
!=================================================================
module gpstate_mod
  use fprecision
  use mpivars
  use grid
  use gmem
  use gdevice, only: gdev_active
  implicit none

  ! Derived type for real space particle field components
  type, public :: GPStateComp
    real(kind=GP)    , allocatable :: rcomp(:) ! real component
  end type GPStateComp

  ! Derived type GPStateArr, whose data is 1D array of GPStateComp's
  type, public :: GPStateArr
    type(GPStateComp), allocatable :: rpstate(:)
  end type GPStateArr

contains

  ! ============ StateComp Allocation routines ==============

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Method to allocate real GPStateComp data types
  !! nc = no.state components
  !! np = no. particles
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine GPState_alloc(pstate, nc, np)
    implicit none
    type(GPStateComp), allocatable, intent(inout) :: pstate(:)
    integer                       , intent   (in) :: nc, np
    integer                                       :: i
    
    call GPState_dealloc(pstate)
    allocate( pstate(nc) )
    do i = 1,nc
      call galloc( pstate(i)%rcomp, np )
    end do
  end subroutine GPState_alloc
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Method to deallocate real GPStateComp data types
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine GPState_dealloc(pstate)
    implicit none
    type(GPStateComp), allocatable, intent(inout) :: pstate(:)
    integer                                       :: i
    
    if ( allocated(pstate) ) then
      do i = 1,size(pstate)
        call gfree( pstate(i)%rcomp )
      end do
      deallocate( pstate )
    endif 
  end subroutine GPState_dealloc

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Method to resize real GPStateComp data types. The
  !! leading min(old,new) values are kept (device and host)
  !! unless keep_data is false.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine GPState_resize(pstate,new_size,keep_data)
    implicit none
    type(GPStateComp), intent(inout), target :: pstate(:)
    integer          , intent(in)            :: new_size
    logical          , intent(in), optional  :: keep_data
    logical                                  :: do_keep
    integer                                  :: i

    do_keep = .TRUE.
    if (present(keep_data)) do_keep = keep_data    
    if (new_size <= 0) stop 'GPState_resize: new_size must be positive.'
    do i = 1,size(pstate)
      call gresize( pstate(i)%rcomp, new_size, do_keep )
    end do
    write(*,*) 'GPstate: GPstate array resized to ', new_size
  end subroutine GPState_resize

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Copies host -> device / device -> host all the components
  !! of a state (no-ops in host builds)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine GPState_update_to(pstate)
    implicit none
    type(GPStateComp), intent(in), target :: pstate(:)
    integer                               :: i
    do i = 1,size(pstate)
      call gupdate_to( pstate(i)%rcomp )
    end do
  end subroutine GPState_update_to

  subroutine GPState_update_from(pstate)
    implicit none
    type(GPStateComp), intent(inout), target :: pstate(:)
    integer                                  :: i
    do i = 1,size(pstate)
      call gupdate_from( pstate(i)%rcomp )
    end do
  end subroutine GPState_update_from

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! dst = src, component by component, in the copy being
  !! worked on (device copy in offload builds)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine GPState_copy(dst,src)
    implicit none
    type(GPStateComp), intent(inout), target :: dst(:)
    type(GPStateComp), intent(in)   , target :: src(:)
    integer                                  :: i
    do i = 1,size(src)
      call gcopy( dst(i)%rcomp, src(i)%rcomp )
    end do
  end subroutine GPState_copy

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! y = y + a*x for the first n entries of every component
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine GPState_axpy(y,a,x,n)
    implicit none
    type(GPStateComp), intent(inout), target :: y(:)
    type(GPStateComp), intent(in)   , target :: x(:)
    real(kind=GP)    , intent(in)            :: a
    integer          , intent(in)            :: n
    integer                                  :: i
    do i = 1,size(y)
      call gp_axpy(n, y(i)%rcomp, a, x(i)%rcomp)
    end do
  end subroutine GPState_axpy

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! y = x + a*y for the first n entries of every component
  !! (the update of the Canuto stepper, y holds the derivative)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine GPState_upd(y,a,x,n)
    implicit none
    type(GPStateComp), intent(inout), target :: y(:)
    type(GPStateComp), intent(in)   , target :: x(:)
    real(kind=GP)    , intent(in)            :: a
    integer          , intent(in)            :: n
    integer                                  :: i
    do i = 1,size(y)
      call gp_upd(n, y(i)%rcomp, a, x(i)%rcomp)
    end do
  end subroutine GPState_upd

  ! ============ Kernels on one component ===================

  subroutine gp_axpy(n,y,a,x)
    implicit none
    integer      , intent(in)    :: n
    real(kind=GP), intent(inout) :: y(n)
    real(kind=GP), intent(in)    :: x(n)
    real(kind=GP), intent(in)    :: a
    integer                      :: k
#if defined(GHOST_GPU)
!$omp target teams distribute parallel do if(target: gdev_active)
    do k = 1,n
#else
!$omp parallel do
    do k = 1,n
#endif
      y(k) = y(k) + a*x(k)
    end do
  end subroutine gp_axpy

  subroutine gp_upd(n,y,a,x)
    implicit none
    integer      , intent(in)    :: n
    real(kind=GP), intent(inout) :: y(n)
    real(kind=GP), intent(in)    :: x(n)
    real(kind=GP), intent(in)    :: a
    integer                      :: k
#if defined(GHOST_GPU)
!$omp target teams distribute parallel do if(target: gdev_active)
    do k = 1,n
#else
!$omp parallel do
    do k = 1,n
#endif
      y(k) = x(k) + a*y(k)
    end do
  end subroutine gp_upd

  ! ============ StateArr Allocation routines ===============

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Method to allocate real GPStateArr data types
  !! narr = no. of GPStateComp's
  !! nc   = no. state components in each GPStateComp's
  !! np   = no. particles
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine GPStateArr_alloc(pstate, narr, nc, np)
    implicit none
    type(GPStateArr), allocatable, intent(inout) :: pstate(:)
    integer                      , intent   (in) :: narr,nc,np
    integer                                      :: i
    
    call GPStateArr_dealloc(pstate)
    allocate( pstate(narr) )
    do i = 1,narr
      call GPState_alloc(pstate(i)%rpstate, nc, np)
    end do
  end subroutine GPStateArr_alloc

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Method to deallocate real GPStateArr data types
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine GPStateArr_dealloc(pstate)
    implicit none
    type(GPStateArr), allocatable, intent(inout) :: pstate(:)
    integer                                      :: i
    
    if ( allocated(pstate) ) then
      do i = 1,size(pstate)
        call GPState_dealloc(pstate(i)%rpstate)
      end do
      deallocate( pstate )
    endif 
  end subroutine GPStateArr_dealloc

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Method to resize real GPStateArr data types
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine GPStateArr_resize(pstate,new_size,keep_data)
    implicit none
    type(GPStateArr), intent(inout), target :: pstate(:)
    integer         , intent(in)            :: new_size
    logical         , intent(in), optional  :: keep_data
    logical                                 :: do_keep
    integer                                 :: i,j

    do_keep = .TRUE.
    if (present(keep_data)) do_keep = keep_data    
    if (new_size <= 0) stop 'GPStateArr_resize: new_size must be positive.'
    do i = 1,size(pstate)
      do j = 1,size(pstate(i)%rpstate)
        call gresize( pstate(i)%rpstate(j)%rcomp, new_size, do_keep )
      end do
    end do
    write(*,*) 'GPstateArr: GPstateArr array resized to ', new_size
  end subroutine GPStateArr_resize

end module gpstate_mod
