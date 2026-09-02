! ===================================================================
! Workspace object, GWorkspace, to checkout and free real and complex
! tmp arrays.
!
! GWorkspace assumes users never call allocate() themselves, always
! return arrays via free_*_tmp, and do not keep aliases to the arrays
! after returning them. Returning the same array twice or returning
! an array not from the pool results in error. Resizing the pool
! invalidates all pointers to arrays in the pool.
!
! Besides the device-resident real and complex pools there is a host-only
! pool (get_*_htmp / free_*_htmp) for the temporaries of routines that
! run on the host copies of the fields, such as the diagnostics; in
! offload builds those arrays have no device copy.
!
! Arrays handed out by get_*_tmp are NOT zeroed: the free_*_tmp
! methods used to zero each array on release, but that is a full
! field sized memset every time a temporary is returned, sixteen of
! them per time step in the HD solver, and removing it makes the
! whole step 4 to 6% faster. The lines are kept, commented out, in
! the three free_*_tmp methods. Callers must therefore initialize
! a temporary before accumulating into it.
!
! For particles (PComp), the pool life cycle is:
!  1. initialize_pool(..., num_pcomp): Allocates num_pcomp PCompEntry
!     slots but leaves their array(:) components unallocated.
!  2. init_pcomp_arrays(partbuff): Called once partbuff is known
!     (after particle init). Allocates array(partbuff) in every slot.
!  3. resize_pcomp_arrays(new_size [, keep_data]): Reallocates every
!     array in the pool to new_size, optionally preserving existing
!     values up to min(old,new). Pointers obtained via get_pcomp_tmp
!     may remain valid after resize (but this is not fully enforced).
! ===================================================================

module class_GWorkspace3D
  USE fprecision
  USE mpivars
  USE grid
  USE gmem

  IMPLICIT NONE
  PRIVATE

  ! ================= Pool entry types ================================
  ! Abstract base type for array pool entries
  type, abstract :: ArrayEntry_Base
    logical      :: is_free = .TRUE.        ! TRUE: free, FALSE: in use
  end type ArrayEntry_Base

  ! Derived type for Real (GP) 3D array entry
  type, EXTENDS(ArrayEntry_Base)  :: RealEntry
    real   (kind=GP), allocatable :: array(:, :, :)
  end type RealEntry

  ! Derived type for Complex (GP) 3D array entry
  type, EXTENDS(ArrayEntry_Base)  :: ComplexEntry
    complex(kind=GP), allocatable :: array(:, :, :)
  end type ComplexEntry

  ! Derived type for 'particle' component entry
  type, EXTENDS(ArrayEntry_Base)  :: PCompEntry
    real   (kind=GP), allocatable :: array(:)
  end type PCompEntry
  
  ! ================= Workspace =======================================
  ! Type to hold the collection of arrays (the pool)
  type, public :: GWorkspace 
    ! Use a list of array entries. Note: For a very large pool, a more
    ! sophisticated data structure (like a linked list) might be better,
    ! but an allocatable array is sufficient for moderate size.
    TYPE   (RealEntry), ALLOCATABLE  :: real_entries_   (:)
    TYPE(ComplexEntry), ALLOCATABLE  :: complex_entries_(:)
    TYPE  (PCompEntry), ALLOCATABLE  :: pcomp_entries_  (:)
    ! Host-only entries: arrays without device copy, for routines that
    ! run on the host copies of the fields (diagnostics, I/O). They cost
    ! no device memory and are handed out by get_*_htmp.
    TYPE   (RealEntry), ALLOCATABLE  :: hreal_entries_   (:)
    TYPE(ComplexEntry), ALLOCATABLE  :: hcomplex_entries_(:)
    integer :: real_size_           = 0
    integer :: complex_size_        = 0
    integer :: pcomp_size_          = 0
    integer :: hreal_size_          = 0
    integer :: hcomplex_size_       = 0
    ! Peak number of arrays checked out at the same time, per pool,
    ! reported by cleanup_pool (used to size the pools per solver)
    integer :: real_peak_           = 0
    integer :: complex_peak_        = 0
    integer :: hreal_peak_          = 0
    integer :: hcomplex_peak_       = 0
    integer :: nreserve_            = 8
    integer :: ncurr_realreserve_   = 8
    integer :: ncurr_complexreserve_= 8
    integer :: ncurr_pcompreserve_  = 8
    integer :: nparts_              = 0       ! current particle buffer size
    logical :: pcomp_initialised_   = .FALSE. ! TRUE after init_pcomp_arrays
  CONTAINS
    procedure, public :: initialize_pool
    procedure, public :: init_pcomp_arrays    ! allocates pcomp data
    procedure, public :: resize_pcomp_arrays  ! resizes pcomp data in-place
    procedure, public :: get_real_tmp     , get_complex_tmp     , get_pcomp_tmp
    procedure, public :: free_real_tmp    , free_complex_tmp    , free_pcomp_tmp 
    procedure, public :: init_host_entries
    procedure, public :: report_peaks
    procedure, public :: get_real_htmp    , get_complex_htmp
    procedure, public :: free_real_htmp   , free_complex_htmp
    procedure, public :: get_real_tmp_size, get_complex_tmp_size, get_pcomp_tmp_size
    procedure, public :: add_real_entries , add_complex_entries , add_pcomp_entries
    procedure, public :: set_nparts, get_nparts
    final             :: cleanup_pool
  end type GWorkspace

  ! The workspace of the run. GHOST has one pool; this pointer lets
  ! the pseudospectral routines, which have no access to the solver
  ! objects, take their field-sized temporaries from it instead of
  ! declaring automatic arrays (which would not exist on the device).
  ! It is set by initialize_pool.
  CLASS(GWorkspace), POINTER, PUBLIC, SAVE :: gws => null()

CONTAINS

  ! ===================================================================
  ! Pool management methods
  ! ===================================================================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Subroutine to initialize the array pool with a specified 
  ! number of arrays.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine initialize_pool(this, num_real, num_complex, num_pcomp)
    CLASS(GWorkspace), intent(inout), target :: this
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
      ! Allocate the derived type contained array (and its device copy)
      call galloc(this%real_entries_(i)%array, nx, ny, ksta, kend)
      this%real_entries_(i)%is_free = .TRUE.
    end do

    ! Initialize Complex entries
    do i = 1, this%complex_size_
      ! Allocate the derived type contained array (and its device copy)
      call galloc(this%complex_entries_(i)%array, nz, ny, ista, iend)
      this%complex_entries_(i)%is_free = .TRUE.
    end do

    gws => this

    ! PComp: create slots only, no array allocation yet.
    if (this%pcomp_size_ > 0) then
      do i = 1, this%pcomp_size_
        this%pcomp_entries_(i)%is_free = .TRUE.
        ! array intentionally not allocated here
      end do
    endif
    
    write(*,*) 'Pool initialized: ', num_real, ' Real,', & 
               num_complex, ' Complex, and', num_pcomp, ' PComp arrays'  
  end subroutine initialize_pool


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Phase-2 init for PComp entries. Call this once partbuff
  ! is known (after particle ctor). Allocates array(partbuff)
  ! in every slot that does not already have an allocated array.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_pcomp_arrays(this, partbuff)
    CLASS(GWorkspace), intent(inout) :: this
    integer          , intent(in)    :: partbuff
    integer                          :: i
 
    if (partbuff <= 0) &
      stop 'init_pcomp_arrays: partbuff must be positive.'
    if (this%pcomp_size_ == 0) &
      stop 'init_pcomp_arrays: no PComp slots; call initialize_pool with num_pcomp > 0.'
 
    do i = 1, this%pcomp_size_
      if (.NOT. ALLOCATED(this%pcomp_entries_(i)%array)) then
        ALLOCATE(this%pcomp_entries_(i)%array(partbuff))
        this%pcomp_entries_(i)%array   = 0.0_GP
        this%pcomp_entries_(i)%is_free = .TRUE.
      end if
    end do
    this%nparts_             = partbuff
    this%pcomp_initialised_  = .TRUE.
 
    write(*,*) 'PComp arrays initialized: ', this%pcomp_size_, &
               ' arrays of size ', partbuff
  end subroutine init_pcomp_arrays


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Resize all PComp arrays in the pool to new_size. If
  ! keep_data is .TRUE. (default), values up to min(old_size,
  ! new_size) are preserved.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine resize_pcomp_arrays(this, new_size, keep_data)
    CLASS(GWorkspace), intent(inout), target :: this
    integer          , intent(in)            :: new_size
    logical          , intent(in), optional  :: keep_data
    logical                                  :: do_keep
    real(kind=GP)    , allocatable           :: tmp(:)
    integer                                  :: i, copy_n
 
    if (.NOT. this%pcomp_initialised_) &
      stop 'resize_pcomp_arrays: call init_pcomp_arrays first.'
    if (new_size <= 0) &
      stop 'resize_pcomp_arrays: new_size must be positive.'
 
    do_keep = .TRUE.
    if (present(keep_data)) do_keep = keep_data
    do i = 1, this%pcomp_size_
      if (ALLOCATED(this%pcomp_entries_(i)%array)) then
        copy_n = min(this%nparts_, new_size)
        ALLOCATE(tmp(new_size))
        if (do_keep) tmp(1:copy_n) = this%pcomp_entries_(i)%array(1:copy_n)
        call MOVE_ALLOC(tmp, this%pcomp_entries_(i)%array)          
      end if
    end do
    this%nparts_ = new_size
    write(*,*) 'Workspace: PComp arrays resized to ', new_size
  end subroutine resize_pcomp_arrays
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Subroutine to clean up and deallocate the entire array 
  ! pool.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Reports (from task 0) the peak number of arrays checked out
  ! at the same time from each pool, to size the pools per solver
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine report_peaks(this)
    CLASS(GWorkspace), intent(in) :: this
    if ( myrank .eq. 0 ) then
      write(*,'(A,4(I0,A))') 'Workspace peak usage: ', this%real_peak_,     &
        ' of ', this%real_size_, ' real, ', this%complex_peak_, ' of ',     &
        this%complex_size_, ' complex'
      write(*,'(A,4(I0,A))') 'Host pool peak usage: ', this%hreal_peak_,    &
        ' of ', this%hreal_size_, ' real, ', this%hcomplex_peak_, ' of ',   &
        this%hcomplex_size_, ' complex'
    endif
  end subroutine report_peaks


  subroutine cleanup_pool(this)
    TYPE(GWorkspace), intent(inout) :: this
    integer :: i

    if (ALLOCATED(this%real_entries_)) THEN
      do i = 1, this%real_size_
        call gfree(this%real_entries_(i)%array)
      end do
      DEALLOCATE(this%real_entries_)
      this%real_size_ = 0
    end if
    if (ALLOCATED(this%complex_entries_)) THEN
      do i = 1, this%complex_size_
        call gfree(this%complex_entries_(i)%array)
      end do
      DEALLOCATE(this%complex_entries_)
      this%complex_size_ = 0
    end if
    if (ALLOCATED(this%hreal_entries_)) THEN
      do i = 1, this%hreal_size_
        call gfree(this%hreal_entries_(i)%array)
      end do
      DEALLOCATE(this%hreal_entries_)
      this%hreal_size_ = 0
    end if
    if (ALLOCATED(this%hcomplex_entries_)) THEN
      do i = 1, this%hcomplex_size_
        call gfree(this%hcomplex_entries_(i)%array)
      end do
      DEALLOCATE(this%hcomplex_entries_)
      this%hcomplex_size_ = 0
    end if
    if (ALLOCATED(this%pcomp_entries_)) THEN
      DEALLOCATE(this%pcomp_entries_)
      this%pcomp_size_ = 0
    end if
  end subroutine cleanup_pool


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Subroutine to add real arrays to the real array pool.
  ! Note that this routine makes all pointers to the arrays
  ! become undefined. It should not be used after initiating
  ! arrays in the pool intended for permanent storage.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine add_real_entries(this, num_new)
    CLASS(GWorkspace), intent(inout), target  :: this
    integer          , intent   (in)          :: num_new
    integer                                   :: i
    type(RealEntry)  , allocatable            :: tmp_copy(:)

    if (num_new<= 0) THEN
      stop 'add_real_entries: real pool size must be positive.'
    end if

    ! Resize this%real_entries_ array, if necessary:    
    if ( num_new .GT. this%ncurr_realreserve_ ) THEN

      ! Need to extend the list of entries. The arrays are moved, not
      ! copied: a copy would place them at new addresses without device
      ! copies and leave stale device associations behind.
      ALLOCATE(tmp_copy(1:this%real_size_+num_new+this%nreserve_))
      do i = 1, this%real_size_
        call MOVE_ALLOC(this%real_entries_(i)%array, tmp_copy(i)%array)
        tmp_copy(i)%is_free = this%real_entries_(i)%is_free
      end do
      call MOVE_ALLOC(tmp_copy, this%real_entries_)
      
      ! Allocate the remaining arrays
      do i = this%real_size_+1,this%real_size_+num_new
        call galloc(this%real_entries_(i)%array, nx, ny, ksta, kend)
        this%real_entries_(i)%is_free = .TRUE.
      end do
      this%ncurr_realreserve_ = this%nreserve_
      
    else
      
      ! Have enough reserves left to fill:
      do i = this%real_size_+1,this%real_size_+num_new
        call galloc(this%real_entries_(i)%array, nx, ny, ksta, kend)
        this%real_entries_(i)%is_free = .TRUE.
        this%ncurr_realreserve_ = this%ncurr_realreserve_ - 1
      end do
    endif

    this%real_size_    = this%real_size_ + num_new
  end subroutine add_real_entries


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Subroutine to add complex arrays to the complex array pool.
  ! Note that this routine makes all pointers to the arrays
  ! become undefined. It should not be used after initiating
  ! arrays in the pool intended for permanent storage.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine add_complex_entries(this, num_new)
    CLASS (GWorkspace), intent(inout), target  :: this
    integer           , intent   (in)          :: num_new
    integer                                    :: i
    type(ComplexEntry), allocatable            :: tmp_copy(:)

    if (num_new<= 0) THEN
      stop 'add_complex_entries: complex pool size must be positive.'
    end if

    ! Resize this%complex_entries_ array, if necessary:    
    if ( num_new .GT. this%ncurr_complexreserve_ ) THEN

      ! Need to extend the list of entries. The arrays are moved, not
      ! copied (see add_real_entries).
      ALLOCATE(tmp_copy(1:this%complex_size_+num_new+this%nreserve_))
      do i = 1, this%complex_size_
        call MOVE_ALLOC(this%complex_entries_(i)%array, tmp_copy(i)%array)
        tmp_copy(i)%is_free = this%complex_entries_(i)%is_free
      end do
      call MOVE_ALLOC(tmp_copy, this%complex_entries_)
      
      ! Allocate the remaining arrays
      do i = this%complex_size_+1,this%complex_size_+num_new
        call galloc(this%complex_entries_(i)%array, nz, ny, ista, iend)
        this%complex_entries_(i)%is_free = .TRUE.
      end do
      this%ncurr_complexreserve_ = this%nreserve_
      
    else
      
      ! Have enough reserves left to fill:
      do i = this%complex_size_+1,this%complex_size_+num_new
        call galloc(this%complex_entries_(i)%array, nz, ny, ista, iend)
        this%complex_entries_(i)%is_free = .TRUE.
        this%ncurr_complexreserve_ = this%ncurr_complexreserve_ - 1
      end do
    endif

    this%complex_size_    = this%complex_size_ + num_new
  end subroutine add_complex_entries

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Subroutine to add PComp slots to the pool. This routine
  ! makes all pointers to the arrays become undefined. Slots
  ! are added without allocating arrays; call init_pcomp_arrays
  ! afterwards (or the new slots will be allocated by the next
  ! init_pcomp_arrays / resize call).
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine add_pcomp_entries(this, num_new)
    CLASS(GWorkspace), intent(inout)          :: this
    integer          , intent(in)             :: num_new
    integer                                   :: i
    type(PCompEntry) , allocatable            :: tmp_copy(:)
 
    if (num_new <= 0) &
      stop 'add_pcomp_entries: num_new must be positive.'
 
    if (.NOT. ALLOCATED(this%pcomp_entries_)) then
      ! Pool was created without pcomp slots; create from scratch
      ALLOCATE(this%pcomp_entries_(num_new + this%nreserve_))
      do i = 1, num_new
        this%pcomp_entries_(i)%is_free = .TRUE.
      end do
      this%ncurr_pcompreserve_ = this%nreserve_
      this%pcomp_size_         = num_new
      return
    end if
 
    if (num_new > this%ncurr_pcompreserve_) then
      ALLOCATE(tmp_copy(1:this%pcomp_size_+num_new+this%nreserve_))
      ! Preserve existing slots
      tmp_copy(1:this%pcomp_size_) = this%pcomp_entries_(1:this%pcomp_size_)
      call MOVE_ALLOC(tmp_copy, this%pcomp_entries_)
      do i = this%pcomp_size_+1, this%pcomp_size_+num_new
        this%pcomp_entries_(i)%is_free = .TRUE.
        ! If pcomp already initialised, allocate data in new slots too
        if (this%pcomp_initialised_) then
          ALLOCATE(this%pcomp_entries_(i)%array(this%nparts_))
          this%pcomp_entries_(i)%array = 0.0_GP
        end if
      end do
      this%ncurr_pcompreserve_ = this%nreserve_
    else
      do i = this%pcomp_size_+1, this%pcomp_size_+num_new
        this%pcomp_entries_(i)%is_free = .TRUE.
        if (this%pcomp_initialised_) then
          ALLOCATE(this%pcomp_entries_(i)%array(this%nparts_))
          this%pcomp_entries_(i)%array = 0.0_GP
        end if
        this%ncurr_pcompreserve_ = this%ncurr_pcompreserve_ - 1
      end do
    endif
 
    this%pcomp_size_ = this%pcomp_size_ + num_new
  end subroutine add_pcomp_entries
 

  ! ===================================================================
  ! Host-only pool
  ! ===================================================================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Allocates the host-only entries (no device copies). May be
  ! called once, after initialize_pool.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_host_entries(this, num_hreal, num_hcomplex)
    CLASS(GWorkspace), intent(inout), target :: this
    integer          , intent(in)            :: num_hreal, num_hcomplex
    integer                                  :: i

    if ( ALLOCATED(this%hreal_entries_) .or. ALLOCATED(this%hcomplex_entries_) ) then
      stop 'init_host_entries: host entries already initialized.'
    end if
    this%hreal_size_    = max(num_hreal,0)
    this%hcomplex_size_ = max(num_hcomplex,0)
    ALLOCATE(this%hreal_entries_   (this%hreal_size_   ))
    ALLOCATE(this%hcomplex_entries_(this%hcomplex_size_))
    do i = 1, this%hreal_size_
      call galloc_host(this%hreal_entries_(i)%array, nx, ny, ksta, kend)
      this%hreal_entries_(i)%is_free = .TRUE.
    end do
    do i = 1, this%hcomplex_size_
      call galloc_host(this%hcomplex_entries_(i)%array, nz, ny, ista, iend)
      this%hcomplex_entries_(i)%is_free = .TRUE.
    end do
    write(*,*) 'Host pool initialized: ', num_hreal, ' Real and', &
               num_hcomplex, ' Complex arrays'
  end subroutine init_host_entries

  subroutine get_real_htmp(this, ret_ptr, success)
    CLASS(GWorkspace), target , intent(inout) :: this
    real   (kind=GP) , pointer, intent(out)   :: ret_ptr(:,:,:)
    logical         , optional, intent(out)   :: success
    integer                                   :: i

    if (present(success)) success = .FALSE.
    do i = 1, this%hreal_size_
      if ( this%hreal_entries_(i)%is_free ) THEN
        this%hreal_entries_(i)%is_free = .FALSE.
        this%hreal_peak_ = max(this%hreal_peak_, this%hreal_size_ - count(this%hreal_entries_(1:this%hreal_size_)%is_free))
        ret_ptr => this%hreal_entries_(i)%array
        if (present(success)) success = .TRUE.
        return
      endif
    enddo
    ret_ptr => null()
    write(*,*) 'GWorkspace::get_real_htmp: all ', this%hreal_size_, &
               ' host real arrays in the pool are checked out'
    error stop 'GWorkspace::get_real_htmp: host real workspace pool exhausted'
  end subroutine get_real_htmp

  subroutine get_complex_htmp(this, ret_ptr, success)
    CLASS(GWorkspace), target , intent(inout) :: this
    complex(kind=GP) , pointer, intent(out)   :: ret_ptr(:,:,:)
    logical         , optional, intent(out)   :: success
    integer                                   :: i

    if (present(success)) success = .FALSE.
    do i = 1, this%hcomplex_size_
      if ( this%hcomplex_entries_(i)%is_free ) THEN
        this%hcomplex_entries_(i)%is_free = .FALSE.
        this%hcomplex_peak_ = max(this%hcomplex_peak_, this%hcomplex_size_ - count(this%hcomplex_entries_(1:this%hcomplex_size_)%is_free))
        ret_ptr => this%hcomplex_entries_(i)%array
        if (present(success)) success = .TRUE.
        return
      endif
    enddo
    ret_ptr => null()
    write(*,*) 'GWorkspace::get_complex_htmp: all ', this%hcomplex_size_, &
               ' host complex arrays in the pool are checked out'
    error stop 'GWorkspace::get_complex_htmp: host complex workspace pool exhausted'
  end subroutine get_complex_htmp

  subroutine free_real_htmp(this, in_ptr)
    CLASS(GWorkspace), target, intent(inout) :: this
    real   (kind=GP), pointer, intent(inout) :: in_ptr(:,:,:)
    integer                                  :: i

    do i = 1, this%hreal_size_
      if (associated(in_ptr, this%hreal_entries_(i)%array)) then
        NULLIFY(in_ptr)
        this%hreal_entries_(i)%is_free = .TRUE.
        return
      endif
    enddo
    stop 'free_real_htmp: array not found in the host pool. Check-in failed'
  end subroutine free_real_htmp

  subroutine free_complex_htmp(this, in_ptr)
    CLASS(GWorkspace), target, intent(inout) :: this
    complex(kind=GP), pointer, intent(inout) :: in_ptr(:,:,:)
    integer                                  :: i

    do i = 1, this%hcomplex_size_
      if (associated(in_ptr, this%hcomplex_entries_(i)%array)) then
        NULLIFY(in_ptr)
        this%hcomplex_entries_(i)%is_free = .TRUE.
        return
      endif
    enddo
    stop 'free_complex_htmp: array not found in the host pool. Check-in failed'
  end subroutine free_complex_htmp


  ! ===================================================================
  ! Size query and nparts helpers
  ! ===================================================================
  
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

    if (present(success)) success = .FALSE.
    do i = 1, this%real_size_
      ! Look for a free entry
      if ( this%real_entries_(i)%is_free ) THEN
        this%real_entries_(i)%is_free = .FALSE. ! Mark as in use
        this%real_peak_ = max(this%real_peak_, this%real_size_ - count(this%real_entries_(1:this%real_size_)%is_free))
        ret_ptr => this%real_entries_(i)%array
        if (present(success)) success = .TRUE.
        return
      endif
    enddo
    ! The pool is exhausted. Every caller dereferences the array it asks
    ! for, so returning an undefined pointer here only turns a sizing
    ! error into a segfault much later. Fail now, and say how big the
    ! pool is: it is sized at start up from NUMTMPREAL plus whatever the
    ! solver, IC and forcing factories reserve, so this is a request to
    ! raise one of those.
    ret_ptr => null()
    write(*,*) 'GWorkspace::get_real_tmp: all ', this%real_size_, &
               ' real arrays in the pool are checked out'
    error stop 'GWorkspace::get_real_tmp: real workspace pool exhausted'
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

    if (present(success)) success = .FALSE.
    do i = 1, this%complex_size_
      ! Look for a free entry
      if ( this%complex_entries_(i)%is_free ) THEN
        this%complex_entries_(i)%is_free = .FALSE. ! Mark as in use
        this%complex_peak_ = max(this%complex_peak_, this%complex_size_ - count(this%complex_entries_(1:this%complex_size_)%is_free))
        ret_ptr => this%complex_entries_(i)%array
        if (present(success)) success = .TRUE.
        return
      endif
    enddo
    ! See the note in get_real_tmp: an exhausted pool is a sizing error,
    ! and returning an undefined pointer only defers it to a segfault.
    ret_ptr => null()
    write(*,*) 'GWorkspace::get_complex_tmp: all ', this%complex_size_, &
               ' complex arrays in the pool are checked out'
    error stop 'GWorkspace::get_complex_tmp: complex workspace pool exhausted'
  end subroutine get_complex_tmp


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Subroutine to check out a free pcomp array from the pool.
  ! Requires init_pcomp_arrays to have been called first.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_pcomp_tmp(this, ret_ptr, success)
    CLASS(GWorkspace), target, intent(inout) :: this
    real(kind=GP)   , pointer, intent(out)   :: ret_ptr(:)
    logical        , optional, intent(out)   :: success
    integer                                  :: i

    if (present(success)) success = .FALSE.
    if ( this%nparts_ .le. 0 ) then
      ret_ptr => null()
      stop 'GWorkspace::get_pcomp_tmp: nparts = 0. Must call set_nparts. '
    endif
    if (.NOT. this%pcomp_initialised_) then
      ret_ptr => null()
      stop 'GWorkspace::get_pcomp_tmp: call init_pcomp_arrays before first checkout.'
    endif

    do i = 1, this%pcomp_size_
      ! Look for a free entry
      if ( this%pcomp_entries_(i)%is_free ) THEN
        this%pcomp_entries_(i)%is_free = .FALSE. ! Mark as in use
        ret_ptr => this%pcomp_entries_(i)%array
        if (present(success)) success = .TRUE.
        return
      endif
    enddo
    ! See the note in get_real_tmp: an exhausted pool is a sizing error,
    ! and returning an undefined pointer only defers it to a segfault.
    ret_ptr => null()
    write(*,*) 'GWorkspace::get_pcomp_tmp: all ', this%pcomp_size_, &
               ' particle arrays in the pool are checked out'
    error stop 'GWorkspace::get_pcomp_tmp: particle workspace pool exhausted'
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
!       this%real_entries_(i)%array = 0.0_GP   ! Optional: Zero the array
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
!       this%complex_entries_(i)%array = 0.0_GP   ! Optional: Zero the array
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
!       this%pcomp_entries_(i)%array = 0.0_GP   ! Optional: Zero the array
        return
      endif
    enddo
    stop 'free_pcomp_tmp: Real array not found in pool. Check-in failed'
  end subroutine free_pcomp_tmp
  
end module class_GWorkspace3D
