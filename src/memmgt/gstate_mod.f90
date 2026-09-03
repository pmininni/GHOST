module gstate_mod
  use fprecision
  use mpivars
  use grid
  use gmem
  implicit none

  ! ================= Base field data types =========================

  ! Derived types for complex (Fourier transformed) field components
  type, public :: GStateComp
    complex(kind=GP), allocatable :: ccomp(:,:,:) ! complex component      
  end type GStateComp
  
  ! Derived type for real space field components
  type, public :: GStateRealComp
    real   (kind=GP), allocatable :: rcomp(:,:,:) ! real component
  end type GStateRealComp

  ! Derived type GCStateArr, whose data is a 1D array of GStateComp's
  type, public :: GCStateArr
    type    (GStateComp), allocatable :: cstate(:)
  end type GCStateArr

  ! Derived type GRStateArr, a 1D array of GStateRealComp's
  type, public :: GRStateArr
    type(GStateRealComp), allocatable :: rstate(:)
 end type GRStateArr

contains

  ! ============ StateComp Allocation routines ==============

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Method to allocate complex GStateComp data types
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine GState_alloc(state, nc)
    use grid
    use mpivars
    implicit none
    type (GStateComp), allocatable, intent(inout) :: state(:)
    integer                       , intent   (in) :: nc
    integer                                       :: i

    if ( allocated(state) ) then
      do i = 1, size(state) 
        call gfree(state(i)%ccomp)
      enddo
    endif
    allocate( state(nc) )
    do i = 1,nc
      call galloc( state(i)%ccomp, nz, ny, ista, iend )
    end do
  end subroutine GState_alloc

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Method to deallocate complex GStateComp data types
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine GState_dealloc(state)
    use grid
    use mpivars
    implicit none
    type (GStateComp), allocatable, intent(inout) :: state(:)
    integer                                       :: i

    if ( allocated(state) ) then
      do i = 1, size(state) 
        call gfree(state(i)%ccomp)
      enddo
    endif
  end subroutine GState_dealloc

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Method to allocate real GStateComp data types
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine GStateReal_alloc(state, nc)
    use grid
    use mpivars
    implicit none
    type (GStateRealComp), allocatable, intent(inout) :: state(:)
    integer                           , intent   (in) :: nc
    integer                                           :: i

    if ( allocated(state) ) then
      do i = 1, size(state) 
        call gfree(state(i)%rcomp)
      enddo
    endif
    allocate( state(nc) )
    do i = 1,nc
      call galloc( state(i)%rcomp, nx, ny, ksta, kend )
    end do
  end subroutine GStateReal_alloc


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Method to deallocate real GStateComp data types
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine GStateReal_dealloc(state)
    use grid
    use mpivars
    implicit none
    type (GStateRealComp), allocatable, intent(inout) :: state(:)
    integer                                           :: i

    if ( allocated(state) ) then
      do i = 1, size(state) 
        call gfree(state(i)%rcomp)
      enddo
    endif
  end subroutine GStateReal_dealloc


  ! ============ Transfers between host and device copies ===========
  ! In host builds these are no-ops (GState_copy is an assignment).

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Copies the host arrays of a state to their device copies
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine GState_update_to(state)
    implicit none
    type (GStateComp), intent(in) :: state(:)
    integer                       :: i
    do i = 1, size(state)
      call gupdate_to(state(i)%ccomp)
    end do
  end subroutine GState_update_to


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Copies the device copies of a state back to the host arrays
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine GState_update_from(state)
    implicit none
    type (GStateComp), intent(inout) :: state(:)
    integer                          :: i
    do i = 1, size(state)
      call gupdate_from(state(i)%ccomp)
    end do
  end subroutine GState_update_from


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Copies one state into another: dst = src. In offload
  !! builds only the device copies are copied (the host copies
  !! are refreshed with GState_update_from when needed), so this
  !! must be used instead of the array assignment dst = src.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine GState_copy(dst, src)
    implicit none
    type (GStateComp), intent(inout) :: dst(:)
    type (GStateComp), intent(in)    :: src(:)
    integer                          :: i
    if ( size(dst) .ne. size(src) ) then
      stop 'GState_copy: inconsistent state sizes'
    endif
    do i = 1, size(src)
      call gcopy(dst(i)%ccomp, src(i)%ccomp)
    end do
  end subroutine GState_copy


  ! ============ StateArr Allocation routines ===============
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Method to allocate complex GCStateArr data types
  !! narr = no. of GStateComp's
  !! nc   = no. state components in each GStateComp's
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine GCStateArr_alloc(state, narr, nc)
    use grid
    use mpivars
    implicit none
    type (GCStateArr), allocatable, intent(inout) :: state(:)
    integer                       , intent   (in) :: nc, narr
    integer                                       :: i

    if ( allocated(state) ) then
      do i = 1, size(state) 
        if ( allocated(state(i)%cstate) ) then
          call GState_dealloc(state(i)%cstate)
        endif
      enddo
      deallocate(state)
    endif
    allocate( state(narr) )
    do i = 1,narr
      call GState_alloc(state(i)%cstate, nc)
    end do
  end subroutine GCStateArr_alloc


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Method to deallocate complex GCStateArr data types
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine GCStateArr_dealloc(state)
    use grid
    use mpivars
    implicit none
    type (GCStateArr), allocatable, intent(inout) :: state(:)
    integer                                       :: i

    if ( allocated(state) ) then
      do i = 1, size(state) 
        if ( allocated(state(i)%cstate) ) then
          call GState_dealloc(state(i)%cstate)
        endif
      enddo
      deallocate(state)
    endif
  end subroutine GCStateArr_dealloc


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Method to allocate real GRStateArr data types
  !! narr = no. of GStateRealComp's
  !! nc   = no. state components in each GStateRealComp's
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine GRStateArr_alloc(state, narr, nc)
    use grid
    use mpivars
    implicit none
    type (GRStateArr), allocatable, intent(inout) :: state(:)
    integer                       , intent   (in) :: nc, narr
    integer                                       :: i

    if ( allocated(state) ) then
      do i = 1, size(state) 
        if ( allocated(state(i)%rstate) ) then
          call GStateReal_dealloc(state(i)%rstate)
        endif
      enddo
      deallocate(state)
    endif
    allocate( state(narr) )
    do i = 1,narr
      call GStateReal_alloc(state(i)%rstate, nc)
    end do
  end subroutine GRStateArr_alloc


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Method to deallocate real GRStateArr data types
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine GRtateArr_dealloc(state)
    use grid
    use mpivars
    implicit none
    type (GRStateArr), allocatable, intent(inout) :: state(:)
    integer                                       :: i

    if ( allocated(state) ) then
      do i = 1, size(state) 
        if ( allocated(state(i)%rstate) ) then
          call GStateReal_dealloc(state(i)%rstate)
        endif
      enddo
      deallocate(state)
    endif
  end subroutine GRtateArr_dealloc

end module gstate_mod
