! ===================================================================
! Workspace object, GWorkspace, to checkout and free real and complex
! tmp arrays.
!
! GWorkspace assumes users never call allocate() themselves, always
! return arrays via free_*_tmp, and do not keep aliases to the arrays
! after returning them. Returning the same array twice or returning
! an array not from the pool results in error.
! ===================================================================

module class_GWorkspace3D
  USE fprecision
  USE mpivars
  USE grid

  IMPLICIT NONE
  PRIVATE

  ! ================= Pool entry types ================================
  ! Abstract base type for array pool entries
  type, abstract :: ArrayEntry_Base
    logical      :: is_free = .TRUE.        ! TRUE: free, FALSE: in use
  end type ArrayEntry_Base

  ! Derived type for Real (GP) 3D array entry
  type, EXTENDS(ArrayEntry_Base)  :: RealEntry
    real(kind=GP), ALLOCATABLE    :: array(:, :, :)
  end type RealEntry

  ! Derived type for Complex (GP) 3D array entry
  type, EXTENDS(ArrayEntry_Base)  :: ComplexEntry
    complex(kind=GP), ALLOCATABLE :: array(:, :, :)
  end type ComplexEntry

  ! Derived type for 'particle' component entry
  type, EXTENDS(ArrayEntry_Base)  :: PCompEntry
    real(kind=GP), ALLOCATABLE    :: array(:)
  end type PCompEntry
  
  ! ================= Workspace =======================================
  ! Type to hold the collection of arrays (the pool)
  type, public :: GWorkspace 
    ! Use a list of array entries. Note: For a very large pool, a more
    ! sophisticated data structure (like a linked list) might be better,
    ! but an allocatable array is sufficient for moderate size.
    CLASS   (RealEntry), ALLOCATABLE  :: real_entries_   (:)
    CLASS(ComplexEntry), ALLOCATABLE  :: complex_entries_(:)
    CLASS  (PCompEntry), ALLOCATABLE  :: pcomp_entries_  (:)
    integer :: real_size_           = 0
    integer :: complex_size_        = 0
    integer :: pcomp_size_          = 0
    integer :: nreserve_            = 8
    integer :: ncurr_realreserve_   = 8
    integer :: ncurr_complexreserve_= 8
    integer :: ncurr_pcompreserve_  = 8
    integer :: nparts_              = 0
  CONTAINS
    procedure, public :: initialize_pool, get_real_tmp, get_complex_tmp, get_pcomp_tmp
    procedure, public :: free_real_tmp    , free_complex_tmp    , free_pcomp_tmp 
    procedure, public :: get_real_tmp_size, get_complex_tmp_size, get_pcomp_tmp_size
    procedure, public :: add_real_entries , add_complex_entries
    procedure, public :: set_nparts, get_nparts
    final             :: cleanup_pool
  end type GWorkspace

CONTAINS

  ! ===================================================================
  ! Pool management methods
  ! ===================================================================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Subroutine to initialize the array pool with a specified 
  ! number of arrays.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine initialize_pool(this, num_real, num_complex, num_pcomp)
    CLASS(GWorkspace), intent(inout) :: this
    integer          , intent(in)    :: num_real
    integer          , intent(in)    :: num_complex
    integer, optional, intent(in)    :: num_pcomp
    integer                          :: i

    if (num_real <= 0) THEN
      stop 'initialize_pool: real pool size must be positive.'
    end if

    if (num_complex <= 0) THEN
      stop 'initialize_pool: complex pool size must be positive.'
    end if
    
    ! Allocate the pool
    this%real_size_    = num_real
    this%complex_size_ = num_complex
    if (present(num_pcomp)) then
      this%pcomp_size_ = num_pcomp
    else
      this%pcomp_size_ = 0
    endif       
    ALLOCATE(this%real_entries_   (this%real_size_   + this%ncurr_realreserve_   ))
    ALLOCATE(this%complex_entries_(this%complex_size_+ this%ncurr_complexreserve_))
    if ( this%pcomp_size_ .gt. 0 ) then
      ALLOCATE(this%pcomp_entries_(this%pcomp_size_  + this%ncurr_pcompreserve_  ))
    endif

    ! Initialize Real entries
    do i = 1, this%real_size_
      ! Allocate the derived type contained array
      ALLOCATE(this%real_entries_(i)%array(nx, ny, ksta:kend))
      this%real_entries_(i)%is_free = .TRUE.
    end do

    ! Initialize Complex entries
    do i = 1, this%complex_size_
      ! Allocate the derived type contained array
      ALLOCATE(this%complex_entries_(i)%array(nz, ny, ista:iend))
      this%complex_entries_(i)%is_free = .TRUE.
    end do

    ! Initialize PComp entries
    if ( this%pcomp_size_ .gt. 0 ) then
      do i = 1, this%pcomp_size_ 
        ! Allocate the derived type contained array
        ALLOCATE(this%pcomp_entries_(i)%array(this%nparts_))
        this%pcomp_entries_(i)%is_free = .TRUE.
      end do
    endif
    
    write(*,*) 'Pool initialized: ', num_real, ' Real and ', & 
               num_complex, ' Complex arrays', num_pcomp, ' PComp arrays'  
  end subroutine initialize_pool


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Subroutine to clean up and deallocate the entire array 
  !  pool.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine cleanup_pool(this)
    type(GWorkspace), intent(inout) :: this

    if (ALLOCATED(this%real_entries_)) THEN
      DEALLOCATE(this%real_entries_)
      this%real_size_ = 0
    end if
    if (ALLOCATED(this%complex_entries_)) THEN
      DEALLOCATE(this%complex_entries_)
      this%complex_size_ = 0
    end if
    if (ALLOCATED(this%pcomp_entries_)) THEN
      DEALLOCATE(this%pcomp_entries_)
      this%pcomp_size_ = 0
    end if
  end subroutine cleanup_pool


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Subroutine to add real arrays to the real array pool.
  !  Note that this routine can make all pointers to the pool
  !  become undefined when the pool is resized with 
  !  MOVE_ALLOC. It should not be used after initiating
  !  arrays in the pool intended for permanent storage.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine add_real_entries(this, num_new)
    CLASS(GWorkspace), intent(inout)          :: this
    integer          , intent(in)             :: num_new
    integer                                   :: i
    type(RealEntry)  , ALLOCATABLE            :: tmp_copy(:)

    if (num_new<= 0) THEN
      stop 'add_real_entries: real pool size must be positive.'
    end if

    ! Resize this%real_entries_ array, if necessary:    
    if ( num_new .GT. this%ncurr_realreserve_ ) THEN

      ! Need to extend arrays and copy old data
      ALLOCATE(tmp_copy(1:this%real_size_+num_new+this%nreserve_))
      tmp_copy(1:this%real_size_) = this%real_entries_(1:this%real_size_)
      call MOVE_ALLOC(tmp_copy, this%real_entries_)
      
      ! Allocate the remaining arrays
      do i = this%real_size_+1,this%real_size_+num_new
        ALLOCATE(this%real_entries_(i)%array(nx, ny, ksta:kend))
        this%real_entries_(i)%is_free = .TRUE.
      end do
      this%ncurr_realreserve_ = this%nreserve_
      
    else
      
      ! Have enough reserves left to fill:
      do i = this%real_size_+1,this%real_size_+num_new
        ALLOCATE(this%real_entries_(i)%array(nx, ny, ksta:kend))
        this%real_entries_(i)%is_free = .TRUE.
        this%ncurr_realreserve_ = this%ncurr_realreserve_ - 1
      end do
    endif

    this%real_size_    = this%real_size_ + num_new
  end subroutine add_real_entries


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Subroutine to add complex arrays to the complex array
  !  pool. Note that this routine can make all pointers to
  !  the pool become undefined when the pool is resized with
  !  MOVE_ALLOC. It should not be used after initiating
  !  arrays in the pool intended for permanent storage.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine add_complex_entries(this, num_new)
    CLASS(GWorkspace), intent(inout)          :: this
    integer          , intent(in)             :: num_new
    integer                                   :: i
    type(ComplexEntry) , ALLOCATABLE          :: tmp_copy(:)

    if (num_new<= 0) THEN
      stop 'add_complex_entries: complex pool size must be positive.'
    end if

    ! Resize this%complex_entries_ array, if necessary:    
    if ( num_new .GT. this%ncurr_complexreserve_ ) THEN

      ! Need to extend arrays and copy old data
      ALLOCATE(tmp_copy(1:this%complex_size_+num_new+this%nreserve_))
      tmp_copy(1:this%complex_size_) = this%complex_entries_(1:this%complex_size_)
      call MOVE_ALLOC(tmp_copy, this%complex_entries_) 
      
      ! Allocate the remaining arrays
      do i = this%complex_size_+1,this%complex_size_+num_new
        ALLOCATE(this%complex_entries_(i)%array(nz, ny, ista:iend))
        this%complex_entries_(i)%is_free = .TRUE.
      end do
      this%ncurr_complexreserve_ = this%nreserve_
      
    else
      
      ! Have enough reserves left to fill:
      do i = this%complex_size_+1,this%complex_size_+num_new
        ALLOCATE(this%complex_entries_(i)%array(nz, ny, ista:iend))
        this%complex_entries_(i)%is_free = .TRUE.
        this%ncurr_complexreserve_ = this%ncurr_complexreserve_ - 1
      end do
    endif

    this%complex_size_    = this%complex_size_ + num_new
  end subroutine add_complex_entries


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Function to return current number of real entries
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PURE function get_real_tmp_size(this) result(num)
    CLASS(GWorkspace), intent(in) :: this
    integer                       :: num
    num = this%real_size_
  end function get_real_tmp_size


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Function to return current number of complex entries
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PURE function get_complex_tmp_size(this) result(num)
    CLASS(GWorkspace), intent(in) :: this
    integer                       :: num
    num = this%complex_size_
  end function get_complex_tmp_size


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Function to return current number of pcomp entries
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PURE function get_pcomp_tmp_size(this) result(num)
    CLASS(GWorkspace), intent(in) :: this
    integer                       :: num
    num = this%pcomp_size_
  end function get_pcomp_tmp_size


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Method to set no. particles
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine set_nparts(this, nparts)
    CLASS(GWorkspace), intent(inout)          :: this
    integer          , intent(in)             :: nparts
    this%nparts_ = nparts
  end subroutine set_nparts


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Function to get no. particles
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PURE function get_nparts(this) result(num)
    CLASS(GWorkspace), intent(in)             :: this
    integer                                   :: num
    num = this%nparts_ 
  end function get_nparts


  ! ===================================================================
  ! Checkout/get and free methods
  ! ===================================================================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Subroutine to check out a free real array from the pool.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_real_tmp(this, ret_ptr, success)
    CLASS(GWorkspace), target, intent(inout) :: this
    real(kind=GP)   , pointer, intent(out)   :: ret_ptr(:,:,:)
    logical        , optional, intent(out)   :: success
    integer                                  :: i

    success = .FALSE.
    do i = 1, this%real_size_
      ! Look for a free entry
      if ( this%real_entries_(i)%is_free ) THEN
        this%real_entries_(i)%is_free = .FALSE. ! Mark as in use
        ret_ptr => this%real_entries_(i)%array
        success = .TRUE.
        return
      endif
    enddo
  end subroutine get_real_tmp


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Subroutine to check out a free complex array from the
  ! pool.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_complex_tmp(this, ret_ptr, success)
    CLASS(GWorkspace), target , intent(inout) :: this
    complex(kind=GP) , pointer, intent(out)   :: ret_ptr(:,:,:)
    logical         , optional, intent(out)   :: success
    integer                                   :: i

    success = .FALSE.
    do i = 1, this%complex_size_
      ! Look for a free entry
      if ( this%complex_entries_(i)%is_free ) THEN
        this%complex_entries_(i)%is_free = .FALSE. ! Mark as in use
        ret_ptr => this%complex_entries_(i)%array
        success = .TRUE.
        return
      endif
    enddo
  end subroutine get_complex_tmp


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Subroutine to check out a free pcomp array from the pool.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_pcomp_tmp(this, ret_ptr, success)
    CLASS(GWorkspace), target, intent(inout) :: this
    real(kind=GP)   , pointer, intent(out)   :: ret_ptr(:)
    logical        , optional, intent(out)   :: success
    integer                                  :: i

    success = .FALSE.
    if ( this%nparts_ .le. 0 ) then
      ret_ptr => null()
      stop 'GWorkspace::get_pcomp_tmp: nparts = 0. Must call set_nparts. '
    endif

    do i = 1, this%pcomp_size_
      ! Look for a free entry
      if ( this%pcomp_entries_(i)%is_free ) THEN
        this%pcomp_entries_(i)%is_free = .FALSE. ! Mark as in use
        ret_ptr => this%pcomp_entries_(i)%array
        success = .TRUE.
        return
      endif
    enddo
  end subroutine get_pcomp_tmp


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Subroutine to check-in/free a real tmp array, making it
  !  available again.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine free_real_tmp(this, in_ptr)
    CLASS(GWorkspace), target , intent(inout) :: this
    real(kind=GP)    , pointer, intent(inout) :: in_ptr(:,:,:)
    integer                                   :: i

    do i = 1, this%real_size_
      if (associated(in_ptr, this%real_entries_(i)%array)) then
        NULLIFY(in_ptr)
        this%real_entries_(i)%is_free = .TRUE. ! Mark as available
        this%real_entries_(i)%array = 0.0_GP   ! Optional: Zero the array
        return
      endif
    enddo
    stop 'free_real_tmp: Real array not found in pool. Check-in failed'
  end subroutine free_real_tmp


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Subroutine to check-in/free a complex tmp array, making
  !  it available again.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine free_complex_tmp(this, in_ptr)
    CLASS(GWorkspace), target, intent(inout) :: this
    complex(kind=GP), pointer, intent(inout) :: in_ptr(:,:,:)
    integer                                  :: i

    do i = 1, this%complex_size_
      if (associated(in_ptr, this%complex_entries_(i)%array)) then
        NULLIFY(in_ptr)
        this%complex_entries_(i)%is_free = .TRUE. ! Mark as available
        this%complex_entries_(i)%array = 0.0_GP   ! Optional: Zero the array
        return
      endif
    enddo
    stop 'free_complex_tmp: Complex array not found in pool. Check-in failed'
  end subroutine free_complex_tmp


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Subroutine to check-in/free a pcomp tmp array, making it
  !  available again.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine free_pcomp_tmp(this, in_ptr)
    CLASS(GWorkspace), target , intent(inout) :: this
    real(kind=GP)    , pointer, intent(inout) :: in_ptr(:)
    integer                                   :: i

    do i = 1, this%pcomp_size_
      if (associated(in_ptr, this%pcomp_entries_(i)%array)) then
        NULLIFY(in_ptr)
        this%pcomp_entries_(i)%is_free = .TRUE. ! Mark as available
        this%pcomp_entries_(i)%array = 0.0_GP   ! Optional: Zero the array
        return
      endif
    enddo
    stop 'free_pcomp_tmp: Real array not found in pool. Check-in failed'
  end subroutine free_pcomp_tmp
  
end module class_GWorkspace3D
