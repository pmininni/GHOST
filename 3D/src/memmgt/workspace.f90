! ===================================================================
! Workspace object, GWorkspace, to checkout, and free
! real and complex tmp arrays. 
! ===================================================================
module class_GWorkspace
  USE fprecision
  USE mpivars
  USE grid

  IMPLICIT NONE
  PRIVATE

  ! Abstract base type for array pool entries
  type, abstract :: ArrayEntry_Base
    integer      :: is_free = 1      ! 1: free, 0: in use
  end type ArrayEntry_Base

  ! Derived type for Real (GP) 3D array entry
  type, EXTendS(ArrayEntry_Base)  :: RealEntry
    real(kind=GP), ALLOCATABLE    :: array(:, :, :)
  end type RealEntry

  ! Derived type for Complex (GP) 3D array entry
  type, EXTendS(ArrayEntry_Base)  :: ComplexEntry
    complex(kind=GP), ALLOCATABLE :: array(:, :, :)
  end type ComplexEntry

  ! Type to hold the collection of arrays (the pool)
  type, public :: GWorkspace 
    ! Use a list of array entries. Note: For a very large pool, a more
    ! sophisticated data structure (like a linked list) might be better,
    ! but an allocatable array of pointers is sufficient for moderate size.
    CLASS   (RealEntry), pointer  :: real_entries_   (:) => NULL()
    CLASS(ComplexEntry), pointer  :: complex_entries_(:) => NULL()
    integer :: real_size_           = 0
    integer :: complex_size_        = 0
    integer :: nreserve_            = 10
    integer :: ncurr_realreserve_   = 10
    integer :: ncurr_complexreserve_= 10

  CONTAINS
    procedure, public :: initialize_pool !, get_real_tmp, get_complex_tmp
!    procedure, public :: free_real_tmp, free_complex_tmp 
    final             :: cleanup_pool

  end type GWorkspace

!  PRIVATE :: get_real_tmp_size, get_complex_tmp_size
  PRIVATE :: add_real_entries !, add_complex_entries
  
CONTAinS

  ! ===================================================================
  ! Pool management methods
  ! ===================================================================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Subroutine to initialize the array pool with a specified 
  ! number of arrays.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine initialize_pool(this, num_real, num_complex)
    class(GWorkspace), intent(inout) :: this
    integer             , intent(in) :: num_real
    integer             , intent(in) :: num_complex
    integer                          :: i
    type   (RealEntry), pointer      :: real_ptr
    type(ComplexEntry), pointer      :: complex_ptr

    if (num_real <= 0) THEN
      stop 'initialize_pool: real pool size must be positive.'
    end if

    if (num_complex <= 0) THEN
      stop 'initialize_pool: complex pool size must be positive.'
    end if

    ! Allocate the pool of pointers
    this%real_size_    = num_real
    this%complex_size_ = num_complex
    ALLOCATE(this%real_entries_   (this%real_size_    + this%ncurr_realreserve_   ))
    ALLOCATE(this%complex_entries_(this%complex_size_ + this%ncurr_complexreserve_))

    ! Initialize Real entries
    do i = 1, num_real
      ! Allocate the derived type contained array
      ALLOCATE(this%real_entries_(i)%array(nx, ny, ksta:kend))
      this%real_entries_(i)%is_free = 1
    end do

    ! Initialize Complex entries
    do i = 1, num_complex
      ! Allocate the derived type contained array
      ALLOCATE(this%complex_entries_(i)%array(nz, ny, ista:iend))
      this%complex_entries_(i)%is_free = 1
    end do

    write(*,*) 'Pool initialized: ', num_real, ' Real and ', num_complex, ' Complex arrays.'
  end subroutine initialize_pool


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Subroutine to clean up and deallocate the entire array 
  !  pool.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine cleanup_pool(this)
    type(GWorkspace), intent(inout) :: this
    integer                          :: i

    if (associated(this%real_entries_)) THEN
      do i = 1, this%real_size_
         DEALLOCATE(this%real_entries_(i)%array)
      end do
      DEALLOCATE(this%real_entries_)
      this%real_size_ = 0
    end if
    if (associated(this%complex_entries_)) THEN
      do i = 1, this%complex_size_
         DEALLOCATE(this%complex_entries_(i)%array)
      end do
      DEALLOCATE(this%complex_entries_)
      this%complex_size_ = 0
    end if
  end subroutine cleanup_pool


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Subroutine to add real arrays to the real array pool
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine add_real_entries(this, num_new)
    class(GWorkspace), intent(inout)          :: this
    integer,           intent(in)             :: num_new
    
    integer                                   :: i
    type(RealEntry)  , ALLOCATABLE, TARGET    :: type_copy(:)

    if (num_new<= 0) THEN
      stop 'add_real_entries: real pool size must be positive.'
    end if

    ! Resize this%real_entries_ array, if necessary:    
    if ( num_new .GT. this%ncurr_realreserve_ ) THEN

      ! Need to extend arrays and copy old data
      ALLOCATE(type_copy(1:this%real_size_))
      type_copy(1:this%real_size_) = this%real_entries_(1:this%real_size_)
      DEALLOCATE(this%real_entries_)
      ALLOCATE(this%real_entries_(1:this%real_size_+num_new+this%nreserve_))
      ! This won't work as target gets deallocated afterwards
      this%real_entries_(1:this%real_size_) => type_copy(1:this%real_size_)
      DEALLOCATE(type_copy)
      
      ! Allocate the remaining arrays
      do i = this%real_size_+1,this%real_size_+num_new
        ALLOCATE(this%real_entries_(i)%array(nx, ny, ksta:kend))
        this%real_entries_(i)%is_free = 1
      end do
      this%ncurr_realreserve_ = this%nreserve_
      
   else
      
      ! Have enough ptr reserves left to fill:
      do i = this%real_size_+1,this%real_size_+num_new
        ALLOCATE(this%real_entries_(i)%array(nx, ny, ksta:kend))
        this%real_entries_(i)%is_free = 1
        this%ncurr_realreserve_ = this%ncurr_realreserve_ - 1
      end do
    endif

    this%real_size_    = this%real_size_ + num_new
  end subroutine add_real_entries


!!$  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$  !  Subroutine to initialize complex array pool with a 
!!$  !  specified number of arrays.
!!$  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$  subroutine add_complex_entries(this, num_new)
!!$    type   (GWorkspace), intent(inout)          :: this
!!$    type(ComplexEntry),  intent(inout), pointer :: complex_ptr
!!$    integer,             intent(in)             :: num_new
!!$
!!$    integer                        :: i
!!$    type(ComplexEntry), pointer    :: ptr_copy(this%complex_size_)
!!$
!!$    if (num_new<= 0) THEN
!!$      stop 'add_complex_entries: real pool size must be positive.'
!!$    end if
!!$
!!$    ! Resize this%complex_entries_ array, if necessary:
!!$    
!!$    if ( num_new .GT. this%ncurr_complexreserve_ ) THEN
!!$      ! Need to extend arrays:
!!$      do i = 1, this%complex_size_
!!$        ptr_copy(i) = this%complex_entries_(i) 
!!$      enddo
!!$      DEALLOCATE(this%complex_entries_)
!!$  
!!$      ALLOCATE(this%complex_entries_(this%complex_size_+num_new+this%nreserve_) )
!!$  
!!$      ! Copy old entries:
!!$      do i = 1, this%complex_size_
!!$        this%complex_entries_(i) = ptr_copy(i)
!!$      end do
!!$  
!!$      do i = 1, num_new
!!$        ! Allocate the derived type and the contained array
!!$        ALLOCATE(complex_ptr)
!!$        ALLOCATE(complex_ptr%array(nz, ny, ista:iend))
!!$        complex_ptr%is_free = 1
!!$        this%complex_entries_(i) => complex_ptr
!!$      end do
!!$      this%ncurr_complexreserve_ = this%nreserve_
!!$    else   
!!$      ! Have enough ptr reserves left to fill:
!!$      do i = 1, num_new
!!$        ! Allocate the derived type and the contained array
!!$        ALLOCATE(complex_ptr)
!!$        ALLOCATE(complex_ptr%array(nx, ny, ksta:kend))
!!$        complex_ptr%is_free = 1
!!$        this%complex_entries_(i) => complex_ptr
!!$        this%ncurr_complexreserve_ = this%ncurr_complexreserve_ - 1
!!$      end do
!!$    endif
!!$
!!$    this%complex_size_    = this%complex_size_ + num_new
!!$  end subroutine add_complex_entries
!!$
!!$
!!$  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$  ! Function to return current number of real entries
!!$  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$  function get_real_tmp_size(this) result(num)
!!$    type   (GWorkspace), intent(inout) :: this
!!$    integer                            :: num
!!$    num = this%real_size_
!!$  end function get_real_tmp_size
!!$
!!$  
!!$  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$  ! Function to return current number of complex entries
!!$  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$  function get_complex_tmp_size(this) result(num)
!!$    type   (GWorkspace), intent(inout) :: this
!!$    integer                            :: num
!!$    num = this%complex_size_
!!$  end function get_complex_tmp_size
!!$
!!$  
!!$  ! ===================================================================
!!$  ! Checkout/get and free methods
!!$  ! ===================================================================
!!$
!!$  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$  ! Function to check out a free real array from the pool.
!!$  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$  function get_real_tmp(this, success) result(ret_ptr)
!!$    type(GWorkspace), intent(inout) :: this
!!$    logical,          intent(out)   :: success
!!$    real(kind=GP),    pointer       :: ret_ptr(:, :, :)
!!$
!!$    integer                         :: i
!!$    type(RealEntry),  pointer       :: real_ptr
!!$
!!$    ret_ptr => NULL()
!!$    success = .FALSE.
!!$
!!$    do i = 1, this%real_size_
!!$      ! Check if the entry is of the correct dynamic type (RealEntry)
!!$      if ( is_same_type(this%real_entries_(i), real_ptr) ) THEN
!!$        ! Downcast the base pointer to the specific derived type
!!$        real_ptr => this%real_entries_(i)
!!$        if ( real_ptr%is_free == 1 ) THEN
!!$          real_ptr%is_free = 0  ! Mark as in use
!!$          ret_ptr => real_ptr%array
!!$          success = .TRUE.
!!$          return
!!$        endif
!!$      endif
!!$    enddo
!!$  end function get_real_tmp
!!$
!!$  
!!$  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$  ! Function to check out a free Complex array from the pool.
!!$  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$  function get_complex_tmp(this, success) result(ret_ptr)
!!$    type(GWorkspace), intent(inout) :: this
!!$    logical,          intent(out)   :: success
!!$    complex(kind=GP), pointer       :: ret_ptr(:, :, :)
!!$
!!$    integer                     :: i
!!$    type(ComplexEntry), pointer :: complex_ptr
!!$
!!$    ret_ptr => NULL()
!!$    success = .FALSE.
!!$
!!$    do i = 1, this%complex_size_
!!$      ! Check if the entry is of the correct dynamic type (ComplexEntry)
!!$      if ( is_same_type(this%complex_entries_(i), complex_ptr) ) THEN
!!$        ! Downcast the base pointer to specific derived type
!!$        complex_ptr => this%complex_entries_(i)
!!$        if ( complex_ptr%is_free == 1 ) THEN
!!$          complex_ptr%is_free = 0  ! Mark as in use
!!$          ret_ptr => complex_ptr%array
!!$          success = .TRUE.
!!$          return
!!$        endif
!!$      endif
!!$    enddo
!!$  end function get_complex_tmp
!!$
!!$
!!$  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$  !  Subroutine to check-in/free a real tmp array, making it 
!!$  !  available again.
!!$  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$  subroutine free_real_tmp(this, in_ptr)
!!$    type(GWorkspace),    intent(inout) :: this
!!$    real(kind=GP), pointer, intent(in) :: in_ptr(:, :, :)
!!$
!!$    integer                  :: i
!!$    type(RealEntry), pointer :: real_ptr
!!$
!!$    do i = 1, this%real_size_
!!$      ! Check if the entry is of the correct dynamic type (RealEntry)
!!$      if ( is_same_type(this%real_entries_(i), real_ptr) ) THEN
!!$        real_ptr => this%real_entries_(i)
!!$        ! Check if the array pointer matches
!!$        if ( associated(real_ptr%array, in_ptr) ) THEN
!!$          if ( real_ptr%is_free == 0 ) THEN
!!$            real_ptr%is_free = 1  ! Mark as free
!!$            ! Optional: Zero the array data for safety/initialization
!!$            real_ptr%array = 0.0_GP
!!$            return
!!$          else
!!$            write(*,*) 'free_real_tmp: Real array already marked as free. Check-in ignored'
!!$            return
!!$          endif
!!$        endif
!!$      endif
!!$    enddo
!!$    stop 'free_real_tmp: Real array not found in pool. Check-in failed'
!!$  end subroutine free_real_tmp
!!$
!!$
!!$  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$  ! Subroutine to free a Complex array, making it 
!!$  ! available again.
!!$  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$  subroutine free_complex_tmp(this, in_ptr)
!!$    type(GWorkspace),       intent(inout) :: this
!!$    complex(kind=GP), pointer, intent(in) :: in_ptr(:, :, :)
!!$
!!$    integer                     :: i
!!$    type(ComplexEntry), pointer :: complex_ptr
!!$
!!$    do i = 1, this%complex_size_
!!$      ! Check if the entry is of the correct dynamic type (ComplexEntry)
!!$      if ( is_same_type(this%complex_entries_(i), complex_ptr) ) THEN
!!$        complex_ptr => this%complex_entries_(i)
!!$        ! Check if the array pointer matches
!!$        if ( associated(complex_ptr%array, in_ptr) ) THEN
!!$          if ( complex_ptr%is_free == 0 ) THEN
!!$            complex_ptr%is_free = 1  ! Mark as free
!!$            ! Optional: Zero the array data for safety/initialization
!!$            complex_ptr%array = (0.0_GP, 0.0_GP)
!!$            return
!!$          else
!!$            write(*,*) 'free_complex_tmp: Complex array already marked as free. Check-in ignored.'
!!$            return
!!$          endif
!!$        endif
!!$      endif
!!$    enddo
!!$    stop 'free_complex_tmp: Complex array not found in pool. Check-in failed'
!!$  end subroutine free_complex_tmp

end module class_GWorkspace
