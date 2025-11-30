! ===================================================================
! Workspace object, TmpPool3D, to checkout, and free
! real and complex tmp arrays. 
! ===================================================================
MODULE TmpPool3D
  USE fprecision
  USE mpivars
  USE grid

  IMPLICIT NONE
  PRIVATE


  ! Abstract base type for array pool entries
  TYPE, ABSTRACT :: ArrayEntry_Base
    INTEGER :: is_free = 1  ! 1: free, 0: in use
  END TYPE ArrayEntry_Base

  ! Derived type for Real (GP) 3D array entry
  TYPE, EXTENDS(ArrayEntry_Base) :: RealEntry
    REAL(KIND=GP), ALLOCATABLE :: array(:, :, :)
  END TYPE RealEntry

  ! Derived type for Complex (GP) 3D array entry
  TYPE, EXTENDS(ArrayEntry_Base) :: ComplexEntry
    COMPLEX(KIND=GP), ALLOCATABLE :: array(:, :, :)
  END TYPE ComplexEntry

  ! Type to hold the collection of arrays (the pool)
  TYPE 
    ! Use a list of array entries. Note: For a very large pool, a more
    ! sophisticated data structure (like a linked list) might be better,
    ! but an allocatable array of pointers is sufficient for moderate size.
    CLASS   (RealEntry), POINTER :: real_entries_   (:) => NULL()
    CLASS(ComplexEntry), POINTER :: complex_entries_(:) => NULL()
    INTEGER :: real_size_    = 0
    INTEGER :: complex_size_ = 0
    INTEGER :: nreserve_     = 10
    INTEGER :: ncurr_realreserve_   = 10
    INTEGER :: ncurr_complexreserve_= 10
  END TYPE 

  ! ===================================================================
  ! INTERFACE BLOCK AND PUBLIC ROUTINES
  ! ===================================================================

  PUBLIC :: initialize_pool, get_real_tmp_tmp, get_complex_tmp_tmp, free_real_tmp_tmp, free_complex_tmp_tmp cleanup_pool

  PRIVATE:: get_real_tmp_size, get_complex_tmp_size, add_real_entries, add_complex_entries

! INTERFACE get_real_tmp_tmp
!   MODULE PROCEDURE get_real_tmp, get_complex_tmp
! END INTERFACE get_real_tmp_tmp

! INTERFACE check_in_array
!   MODULE PROCEDURE free_real_tmp, free_complex_tmp
! END INTERFACE check_in_array

CONTAINS

  ! ===================================================================
  ! Pool management methods
  ! ===================================================================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Subroutine to initialize the array pool with a specified 
  ! number of arrays.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE initialize_pool(this, num_real, num_complex)
    TYPE(TmpPool3D), INTENT(INOUT) :: this
    INTEGER        , INTENT(IN) :: num_real
    INTEGER        , INTENT(IN) :: num_complex

    INTEGER                     :: i, total_size
    TYPE   (RealEntry), POINTER :: real_ptr
    TYPE(ComplexEntry), POINTER :: complex_ptr

    IF (num_real <= 0) THEN
      STOP 'initialize_pool: real pool size must be positive.'
    END IF

    IF (num_complex <= 0) THEN
      STOP 'initialize_pool: complex pool size must be positive.'
    END IF

    ! Allocate the pool of pointers
    this%real_size_    = num_real
    this%complex_size_ = num_complex
    ALLOCATE(this%real_entries_   (this%real_size_) )
    ALLOCATE(this%complex_entries_(this%complex_size_) )

    ! Initialize Real entries
    DO i = 1, num_real
      ! Allocate the derived type and the contained array
      ALLOCATE(real_ptr)
      ALLOCATE(real_ptr%array(nx, ny, ksta:kend))
      real_ptr%is_free = 1
      this%real_entries_(i) => real_ptr
    END DO

    ! Initialize Complex entries
    DO i = 1, num_complex
      ! Allocate the derived type and the contained array
      ALLOCATE(complex_ptr)
      ALLOCATE(complex_ptr%array(nz, ny, ista:iend))
      complex_ptr%is_free = 1
      this%complex_entries_(i) => complex_ptr
    END DO

    WRITE(*,*) 'Pool initialized: ', num_real, ' Real and ', num_complex, ' Complex arrays.'
  END SUBROUTINE initialize_pool


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Subroutine to clean up and deallocate the entire array 
  !  pool.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE cleanup_pool(this)
    TYPE(TmpPool3D), INTENT(INOUT) :: this
    INTEGER :: i

    IF (ASSOCIATED(this%real_entries_)) THEN
      DO i = 1, this%real_size_
        IF (ASSOCIATED(this%real_entries_(i))) THEN
          DEALLOCATE(this%real_entries_(i))
        END IF
      END DO
      DO i = 1, this%complex_size_
        IF (ASSOCIATED(this%complex_entries_(i))) THEN
          DEALLOCATE(this%complex_entries_(i))
        END IF
      END DO
      ! Deallocate the array of pointers
      DEALLOCATE(this%real_entries_)
      DEALLOCATE(this%complex_entries_)
      this%real_size_ = 0
      this%complex_size_ = 0
    END IF
  END SUBROUTINE cleanup_pool


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Subroutine to initialize the array pool with a 
  !  specified number of arrays.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE add_real_entries(this, num_new)
    TYPE(TmpPool3D), INTENT(INOUT) :: this
    TYPE(RealEntry), INTENT(INOUT), POINTER :: real_ptr

    INTEGER                     :: i
    TYPE(RealEntry), POINTER    :: ptr_copy(this%real_size_)

    IF (num_new<= 0) THEN
      STOP 'add_real_entries: real pool size must be positive.'
    END IF

    ! Resize this%real_entries array, necessary:
   
    
    IF ( num_new .GT. this%ncurr_realreserve_ ) THEN
      ! Need to extend arrays:
      DO i = 1, this%real_size_
        ptr_copy(i) = this%real_entries_(i) 
      ENDDO
      DEALLOCATE(this%real_entries_)
  
      ALLOCATE(this%real_entries_(this%real_size_+num_new+this%nreserve_) )
  
      ! Copy old entries:
      DO i = 1, this%real_size_
        this%real_entries_(i) = ptr_copy(i)
      END DO
  
      DO i = 1, num_new
        ! Allocate the derived type and the contained array
        ALLOCATE(real_ptr)
        ALLOCATE(real_ptr%array(nx, ny, ksta:kend))
        real_ptr%is_free = 1
        this%real_entries_(i) => real_ptr
      END DO
      this%ncurr_realreserve_ = this%nresrve_
    ELSE   
      ! Have enough ptr reserves left to fill:
      DO i = 1, num_new
        ! Allocate the derived type and the contained array
        ALLOCATE(real_ptr)
        ALLOCATE(real_ptr%array(nx, ny, ksta:kend))
        real_ptr%is_free = 1
        this%real_entries_(i) => real_ptr
        this%ncurr_realreserve_ = this%ncurr_realreserve_ - 1
      END DO
    ENDIF

    this%real_size_    = this%real_size_ + num_new

  END SUBROUTINE add_real_entries


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Subroutine to initialize complex array pool with a 
  !  specified number of arrays.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE add_complex_entries(this, num_new)
    TYPE   (TmpPool3D), INTENT(INOUT) :: this
    TYPE(ComplexEntry), INTENT(INOUT), POINTER :: complex_ptr

    INTEGER                        :: i
    TYPE(ComplexEntry), POINTER    :: ptr_copy(this%complex_size_)

    IF (num_new<= 0) THEN
      STOP 'add_complex_entries: real pool size must be positive.'
    END IF

    ! Resize this%complex_entries array, if necessary:
    
    IF ( num_new .GT. this%ncurr_complexreserve_ ) THEN
      ! Need to extend arrays:
      DO i = 1, this%complex_size_
        ptr_copy(i) = this%complex_entries_(i) 
      ENDDO
      DEALLOCATE(this%complex_entries_)
  
      ALLOCATE(this%complex_entries_(this%complex_size_+num_new+this%nreserve_) )
  
      ! Copy old entries:
      DO i = 1, this%complex_size_
        this%complex_entries_(i) = ptr_copy(i)
      END DO
  
      DO i = 1, num_new
        ! Allocate the derived type and the contained array
        ALLOCATE(complex_ptr)
        ALLOCATE(complex_ptr%array(nz, ny, ista:iend))
        complex_ptr%is_free = 1
        this%complex_entries_(i) => complex_ptr
      END DO
      this%ncurr_complexreserve_ = this%nresrve_
    ELSE   
      ! Have enough ptr reserves left to fill:
      DO i = 1, num_new
        ! Allocate the derived type and the contained array
        ALLOCATE(complex_ptr)
        ALLOCATE(complex_ptr%array(nx, ny, ksta:kend))
        real_ptr%is_free = 1
        this%complex_entries_(i) => complex_ptr
        this%ncurr_complexreserve_ = this%ncurr_complexreserve_ - 1
      END DO
    ENDIF

    this%complex_size_    = this%complex_size_ + num_new

  END SUBROUTINE add_complex_entries

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Function to return current number of real entries
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION get_real_tmp_size(this, num_new) return(num)
    TYPE   (TmpPool3D), INTENT(INOUT) :: this
    INTEGER                           :: num
    num = this%real_size_
  END FUNCTION get_real_tmp_size

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Function to return current number of complex entries
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION get_complex_tmp_size(this, num_new) return(num)
    TYPE   (TmpPool3D), INTENT(INOUT) :: this
    INTEGER                           :: num
    num = this%complex_size_
  END FUNCTION get_complex_tmp_size

  ! ===================================================================
  ! Checkout/get and free methods
  ! ===================================================================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Function to check out a free real array from the pool.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION get_real_tmp(this, success) RESULT(ret_ptr)
    TYPE(TmpPool3D), INTENT(INOUT) :: this
    LOGICAL, INTENT(OUT) :: success
    REAL(KIND=GP), POINTER :: ret_ptr(:, :, :)

    INTEGER :: i
    TYPE(RealEntry), POINTER :: real_ptr

    ret_ptr => NULL()
    success = .FALSE.

    DO i = 1, this%size
      ! Check if the entry is of the correct dynamic type (RealEntry)
      IF ( IS_SAME_TYPE(this%entries(i), real_ptr) ) THEN
        ! Downcast the base pointer to the specific derived type
        real_ptr => this%entries(i)
        IF ( real_ptr%is_free == 1 ) THEN
          real_ptr%is_free = 0  ! Mark as in use
          ret_ptr => real_ptr%array
          success = .TRUE.
          RETURN
        ENDIF
      ENDIF
    ENDDO

  END FUNCTION get_real_tmp

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Function to check out a free Complex array from the pool.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION get_complex_tmp(this, success) RESULT(ret_ptr)
    TYPE(TmpPool3D), INTENT(INOUT) :: this
    LOGICAL, INTENT(OUT) :: success
    COMPLEX(KIND=GP), POINTER :: ret_ptr(:, :, :)

    INTEGER :: i
    TYPE(ComplexEntry), POINTER :: complex_ptr

    ret_ptr => NULL()
    success = .FALSE.

    DO i = 1, this%complex_size_
      ! Check if the entry is of the correct dynamic type (ComplexEntry)
      IF ( IS_SAME_TYPE(this%entries(i), complex_ptr) ) THEN
        ! Downcast the base pointer to specific derived type
        complex_ptr => this%entries(i)
        IF ( complex_ptr%is_free == 1 ) THEN
          complex_ptr%is_free = 0  ! Mark as in use
          ret_ptr => complex_ptr%array
          success = .TRUE.
          RETURN
        ENDIF
      ENDIF
    ENDDO

  END FUNCTION get_complex_tmp


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Subroutine to check-in/free a real tmp array, making it 
  !  available again.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE free_real_tmp(this, in_ptr)
    TYPE(TmpPool3D), INTENT(INOUT) :: this
    REAL(KIND=GP), POINTER, INTENT(IN) :: in_ptr(:, :, :)

    INTEGER :: i
    TYPE(RealEntry), POINTER :: real_ptr

    DO i = 1, this%real_size_
      ! Check if the entry is of the correct dynamic type (RealEntry)
      IF ( IS_SAME_TYPE(this%entries(i), real_ptr) ) THEN
        real_ptr => this%entries(i)
        ! Check if the array pointer matches
        IF ( ASSOCIATED(real_ptr%array, in_ptr) ) THEN
          IF ( real_ptr%is_free == 0 ) THEN
            real_ptr%is_free = 1  ! Mark as free
            ! Optional: Zero the array data for safety/initialization
            real_ptr%array = 0.0_GP
            RETURN
          ELSE
            WRITE(*,*) 'free_real_tmp: Real array already marked as free. Check-in ignored'
            RETURN
          ENDIF
        ENDIF
      ENDIF
    ENDDO

    STOP 'free_real_tmp: Real array not found in pool. Check-in failed'

  END SUBROUTINE free_real_tmp


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Subroutine to free a Complex array, making it 
  ! available again.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE free_complex_tmp(this, in_ptr)
    TYPE(TmpPool3D), INTENT(INOUT) :: this
    COMPLEX(KIND=GP), POINTER, INTENT(IN) :: in_ptr(:, :, :)

    INTEGER :: i
    TYPE(ComplexEntry), POINTER :: complex_ptr

    DO i = 1, this%complex_size_
      ! Check if the entry is of the correct dynamic type (ComplexEntry)
      IF ( IS_SAME_TYPE(this%entries(i), complex_ptr) ) THEN
        complex_ptr => this%entries(i)
        ! Check if the array pointer matches
        IF ( ASSOCIATED(complex_ptr%array, in_ptr) ) THEN
          IF ( complex_ptr%is_free == 0 ) THEN
            complex_ptr%is_free = 1  ! Mark as free
            ! Optional: Zero the array data for safety/initialization
            complex_ptr%array = (0.0_GP, 0.0_GP)
            RETURN
          ELSE
            WRITE(*,*) 'free_complex_tmp: Complex array already marked as free. Check-in ignored.'
            RETURN
          ENDIF
        ENDIF
      ENDIF
    ENDDO

    STOP 'free_complex_tmp: Complex array not found in pool. Check-in failed'
  END SUBROUTINE free_complex_tmp

END MODULE TmpPool3D
