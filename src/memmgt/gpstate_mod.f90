module gpstate_mod
  use fprecision
  use mpivars
  use grid
  implicit none

  ! ================= Base field data types =========================
  
  ! Derived type for real space particle 
  ! field components
  type, public :: GPStateComp
    real(kind=GP)   , allocatable :: rcomp(:) ! real component
  end type GPStateComp

  ! Derived type GPState, whose data 
  ! is a 1d pointer array of GPStateComp's:
  type, abstract :: GPState
    type(GPStateComp), allocatable, dimension(:) :: rpstate
  end type GPState

contains

  ! ================= Allocation routines ===================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Method to allocate real GPState data types
  !! nc = no.state components
  !! np = no. particles
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine GPState_alloc(pstate, nc, np)
    use grid
    use mpivars
    implicit none
    type(GPStateComp), allocatable, intent(inout) :: pstate(:)
    integer                       , intent   (in) :: nc, np
    integer                                       :: i
    
    if ( allocated(pstate) ) then
      do i = 1,size(pstate)
        if ( allocated(pstate(i)%rcomp) ) then
          deallocate( pstate(i)%rcomp )
        endif
      end do
    endif 
    allocate( pstate(nc) )
    do i = 1,nc
      allocate( pstate(i)%rcomp(np) )
    end do
  end subroutine GPState_alloc
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Method to deallocate real GPState data types
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine GPState_dealloc(pstate)
    use grid
    use mpivars
    implicit none
    type(GPStateComp), allocatable, intent(inout) :: pstate(:)
    integer                                       :: i
    
    if ( allocated(pstate) ) then
      do i = 1,size(pstate)
        if ( allocated(pstate(i)%rcomp) ) then
          deallocate( pstate(i)%rcomp )
        endif
      end do
    endif 
  end subroutine GPState_dealloc

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Method to resize real GPState data types
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine GPState_resize(pstate,new_size,keep_data)
    use grid
    use mpivars
    implicit none
    type(GPStateComp), intent(inout)        :: pstate(:)
    real(kind=GP)    , allocatable          :: tmp(:)
    integer          , intent(in)           :: new_size
    logical          , intent(in), optional :: keep_data
    logical                                 :: do_keep
    integer                                 :: i,copy_n

    do_keep = .TRUE.
    if (present(keep_data)) do_keep = keep_data    
    if (new_size <= 0) stop 'GPState_resize: new_size must be positive.'
    do i = 1,size(pstate)
      if ( allocated(pstate(i)%rcomp) ) then
        copy_n = min(size(pstate(i)%rcomp), new_size)
        allocate(tmp(new_size))
        if (do_keep) tmp(1:copy_n) = pstate(i)%rcomp(1:copy_n)
        call MOVE_ALLOC(tmp, pstate(i)%rcomp)
      endif
    end do
  end subroutine GPState_resize

end module gpstate_mod
