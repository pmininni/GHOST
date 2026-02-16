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
    real(kind=GP)   , allocatable :: rcomp(:,:,:) ! real component
  end type GStateRealComp

  ! Derived type GCState, whose data 
  ! is a 1d pointer array of complex components:
  type, abstract :: GCState
    type    (GStateComp), allocatable, dimension(:) :: cstate_
    contains
      procedure, public  :: data => GCState_data
      procedure, private :: GCState_get_comp
      generic            :: operator(.get.) => GCState_get_comp
  end type GCState

  ! Derived type GRState, whose data 
  ! is a 1d pointer array of real components:
  type, abstract :: GRState
    type(GStateRealComp), allocatable, dimension(:) :: rstate_
    contains
      procedure, public  :: data => GRState_data
      procedure, private :: GRState_get_comp
      generic            :: operator(.get.) => GRState_get_comp
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
    allocate( state(nc) )
    do i = 1,nc
      allocate( state(i)%ccomp(nz,ny,ista:iend) )
    end do
  end subroutine GState_alloc

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Method to allocate real GState data types
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine GStateReal_alloc(state, nc)
    use grid
    use mpivars
    implicit none
    type (GStateRealComp), allocatable, intent(inout) :: state(:)
    integer                           , intent   (in) :: nc
    integer                                           :: i
    allocate( state(nc) )
    do i = 1,nc
      allocate( state(i)%rcomp(nx,ny,ksta:kend) )
    end do
  end subroutine GStateReal_alloc
  
  ! ================= Data access routines  =================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Method to get GCState complex data
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function GCState_data(this) result(ret)
    implicit none
    class  (GCState), target, intent(inout) :: this
    type(GStateComp), pointer, dimension(:) :: ret 
    ret => this%cstate_
  end function GCState_data

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Method to get GRState real data
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function GRState_data(this) result(ret)
    implicit none
    class      (GRState), target, intent(inout) :: this
    type(GStateRealComp), pointer, dimension(:) :: ret 
    ret => this%rstate_
  end function GRState_data

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Method to get GCState complex component data
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function GCState_get_comp(this, i) result(ret)
    implicit none
    class  (GCState), target, intent(in) :: this
    type(GStateComp), pointer :: ret 
    integer,       intent(in) :: i

    if ( (i .lt. 1).or.(i .gt. size(this%cstate_)) ) then
      stop 'GCState_comp: Invalid index'
    endif
    ret => this%cstate_(i)
  end function GCState_get_comp

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Method to get GRState real component data
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function GRState_get_comp(this, i) result(ret)
    implicit none
    class      (GRState), target, intent(in) :: this
    type(GStateRealComp), pointer :: ret 
    integer,           intent(in) :: i

    if ( (i .lt. 1).or.(i .gt. size(this%rstate_)) ) then
      stop 'GRState_comp: Invalid index'
    endif
    ret => this%rstate_(i)
  end function GRState_get_comp

end module gstate_mod
