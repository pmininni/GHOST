module gpstate_mod
  use fprecision
  use mpivars
  use grid
  implicit none

  ! ================= Base field data types =========================
  
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
  !! Method to deallocate real GPStateComp data types
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
  !! Method to resize real GPStateComp data types
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine GPState_resize(pstate,new_size,keep_data)
    use grid
    use mpivars
    implicit none
    type(GPStateComp), intent(inout), target :: pstate(:)
    real(kind=GP)    , allocatable           :: tmp(:)
    integer          , intent(in)            :: new_size
    logical          , intent(in), optional  :: keep_data
    logical                                  :: do_keep
    integer                                  :: i,copy_n

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
    write(*,*) 'GPstate: GPstate array resized to ', new_size
  end subroutine GPState_resize


  ! ============ StateArr Allocation routines ===============

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Method to allocate real GPStateArr data types
  !! narr = no. of GPStateComp's
  !! nc   = no. state components in each GPStateComp's
  !! np   = no. particles
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine GPStateArr_alloc(pstate, narr, nc, np)
    use grid
    use mpivars
    implicit none
    type(GPStateArr), allocatable, intent(inout) :: pstate(:)
    integer                      , intent   (in) :: narr,nc,np
    integer                                      :: i
    
    if ( allocated(pstate) ) then
      do i = 1,size(pstate)
        if ( allocated(pstate(i)%rpstate) ) then
          call GPState_dealloc(pstate(i)%rpstate)
        endif
      end do
    endif 
    allocate( pstate(narr) )
    do i = 1,narr
      call GPState_alloc(pstate(i)%rpstate, nc, np)
    end do
  end subroutine GPStateArr_alloc

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Method to deallocate real GPStateArr data types
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine GPStateArr_dealloc(pstate)
    use grid
    use mpivars
    implicit none
    type(GPStateArr), allocatable, intent(inout) :: pstate(:)
    integer                                      :: i
    
    if ( allocated(pstate) ) then
      do i = 1,size(pstate)
        if ( allocated(pstate(i)%rpstate) ) then
          call GPState_dealloc(pstate(i)%rpstate)
        endif
      end do
    endif 
  end subroutine GPStateArr_dealloc

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Method to resize real GPStateArr data types
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine GPStateArr_resize(pstate,new_size,keep_data)
    use grid
    use mpivars
    implicit none
    type(GPStateArr), intent(inout), target :: pstate(:)
    integer         , intent(in)            :: new_size
    logical         , intent(in), optional  :: keep_data
    logical                                 :: do_keep
    real(kind=GP)   , allocatable           :: tmp(:)
    integer                                 :: i,j,copy_n

    do_keep = .TRUE.
    if (present(keep_data)) do_keep = keep_data    
    if (new_size <= 0) then
      stop 'GPStateArr_resize: new_size must be positive.'
    end if
    do i = 1,size(pstate)
      do j = 1,size(pstate(i)%rpstate)
      if ( allocated(pstate(i)%rpstate(j)%rcomp) ) then
        copy_n = min(size(pstate(i)%rpstate(j)%rcomp), new_size)
        allocate(tmp(new_size))
        if (do_keep) tmp(1:copy_n) = pstate(i)%rpstate(j)%rcomp(1:copy_n)
        call MOVE_ALLOC(tmp, pstate(i)%rpstate(j)%rcomp)
      endif
      end do
    end do
    write(*,*) 'GPstateArr: GPstateArr array resized to ', new_size
  end subroutine GPStateArr_resize

end module gpstate_mod
