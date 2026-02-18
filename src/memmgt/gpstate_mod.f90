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
    contains
      procedure, public  :: data => GPState_data
!     procedure, private :: GPState_get_comp
!     generic            :: operator(.get.) => GPState_get_comp
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
      do i = 1,size(state)
        if ( allocated(pstate(i)%rcomp ) then
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
      do i = 1,size(state)
        if ( allocated(pstate(i)%rcomp ) then
          deallocate( pstate(i)%rcomp )
        endif
      end do
    endif 

  end subroutine GPState_dealloc
  
  ! ================= Data access routines  =================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Method to get GPState real data
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function GPState_data(this) result(ret)
    implicit none
    class    (GPState), target,    intent(in) :: this
    type (GPStateComp), pointer, dimension(:) :: ret 
    ret => this%rpstate
  end function GPState_data


! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !! Method to get GPState component data
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! function GPState_get_comp(this, i) result(ret)
!   implicit none
!   class   (GPState), target, intent(in) :: this
!   type(GPStateComp), pointer            :: ret 
!   integer          ,         intent(in) :: i
!
!   if ( (i .lt. 1).or.(i .gt. size(this%rpstate)) ) then
!     stop 'GPState_get_comp: Invalid index'
!   endif
!   ret => this%rpstate(i);
! end function GPState_get_comp

end module gpstate_mod
