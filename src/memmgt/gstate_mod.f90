module gstate_mod
  use fprecision
  use mpivars
  use grid
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

  ! Derived type GCState, whose data is a vector of complex components
  type, public :: GCState
    type    (GStateComp), allocatable :: cstate(:)
  end type GCState

  ! Derived type GRState, whose data is a vector of real components
  type, public :: GRState
    type(GStateRealComp), allocatable :: rstate(:)
  end type GRState

contains

  ! ================= Allocation routines ===================

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
        if ( allocated(state(i)%ccomp) ) then
          deallocate(state(i)%ccomp)
        endif
      enddo
    endif
    allocate( state(nc) )
    do i = 1,nc
      allocate( state(i)%ccomp(nz,ny,ista:iend) )
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
        if ( allocated(state(i)%ccomp) ) then
          deallocate(state(i)%ccomp)
        endif
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
        if ( allocated(state(i)%rcomp) ) then
          deallocate(state(i)%rcomp)
        endif
      enddo
    endif
    allocate( state(nc) )
    do i = 1,nc
      allocate( state(i)%rcomp(nx,ny,ksta:kend) )
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
        if ( allocated(state(i)%rcomp) ) then
          deallocate(state(i)%rcomp)
        endif
      enddo
    endif
  end subroutine GStateReal_dealloc

end module gstate_mod
