MODULE TmpPool3D
  USE fprecision
  USE mpivars
  USE grid

  IMPLICIT NONE
  PRIVATE


  ! ===================================================================
  ! TYPE DEFINITIONS
  ! ===================================================================

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
  PRIVAT :: get_real_tmp_size, get_complex_tmp_size, add_real_entries, add_complex_entries

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
  ! Checkout/get methods
  ! ===================================================================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Function to check out a free Real array from the pool.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION get_real_tmp(this, success) RESULT(array_ptr)
    TYPE(TmpPool3D), INTENT(INOUT) :: this
    LOGICAL, INTENT(OUT) :: success
    REAL(KIND=GP), POINTER :: array_ptr(:, :, :)

    INTEGER :: i
    TYPE(RealEntry), POINTER :: real_ptr

    array_ptr => NULL()
    success = .FALSE.

    DO i = 1, this%size
      ! Check if the entry is of the correct dynamic type (RealEntry)
      IF (IS_SAME_TYPE(this%entries(i), real_ptr)) THEN
        ! Downcast the base pointer to the specific derived type
        real_ptr => this%entries(i)
        IF (real_ptr%is_free == 1) THEN
          real_ptr%is_free = 0  ! Mark as in use
          array_ptr => real_ptr%array
          success = .TRUE.
          RETURN
        END IF
      END IF
    END DO

  END FUNCTION get_real_tmp

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Function to check out a free Complex array from the pool.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION get_complex_tmp(this, success) RESULT(array_ptr)
    TYPE(TmpPool3D), INTENT(INOUT) :: this
    LOGICAL, INTENT(OUT) :: success
    COMPLEX(KIND=GP), POINTER :: array_ptr(:, :, :)

    INTEGER :: i
    TYPE(ComplexEntry), POINTER :: complex_ptr

    array_ptr => NULL()
    success = .FALSE.

    DO i = 1, this%size
      ! Check if the entry is of the correct dynamic type (ComplexEntry)
      IF (IS_SAME_TYPE(this%entries(i), complex_ptr)) THEN
        ! Downcast the base pointer to the specific derived type
        complex_ptr => this%entries(i)
        IF (complex_ptr%is_free == 1) THEN
          complex_ptr%is_free = 0  ! Mark as in use
          array_ptr => complex_ptr%array
          success = .TRUE.
          RETURN
        END IF
      END IF
    END DO

  END FUNCTION get_complex_tmp

  ! ===================================================================
  ! CHECK-IN ROUTINES
  ! ===================================================================

  !> Subroutine to check in a Real array, making it available again.
  SUBROUTINE free_real_tmp(this, array_ptr)
    TYPE(TmpPool3D), INTENT(INOUT) :: this
    REAL(KIND=GP), POINTER, INTENT(IN) :: array_ptr(:, :, :)

    INTEGER :: i
    TYPE(RealEntry), POINTER :: real_ptr

    DO i = 1, this%size
      ! Check if the entry is of the correct dynamic type (RealEntry)
      IF (IS_SAME_TYPE(this%entries(i), real_ptr)) THEN
        real_ptr => this%entries(i)
        ! Check if the array pointer matches
        IF (ASSOCIATED(real_ptr%array, array_ptr)) THEN
          IF (real_ptr%is_free == 0) THEN
            real_ptr%is_free = 1  ! Mark as free
            ! Optional: Zero the array data for safety/initialization
            real_ptr%array = 0.0_DP
            RETURN
          ELSE
            WRITE(*,*) 'Error: Real array already marked as free. Check-in ignored.'
            RETURN
          END IF
        END IF
      END IF
    END DO

    WRITE(*,*) 'Error: Real array not found in pool. Check-in failed.'
  END SUBROUTINE free_real_tmp

  !> Subroutine to check in a Complex array, making it available again.
  SUBROUTINE free_complex_tmp(this, array_ptr)
    TYPE(TmpPool3D), INTENT(INOUT) :: this
    COMPLEX(KIND=GP), POINTER, INTENT(IN) :: array_ptr(:, :, :)

    INTEGER :: i
    TYPE(ComplexEntry), POINTER :: complex_ptr

    DO i = 1, this%size
      ! Check if the entry is of the correct dynamic type (ComplexEntry)
      IF (IS_SAME_TYPE(this%entries(i), complex_ptr)) THEN
        complex_ptr => this%entries(i)
        ! Check if the array pointer matches
        IF (ASSOCIATED(complex_ptr%array, array_ptr)) THEN
          IF (complex_ptr%is_free == 0) THEN
            complex_ptr%is_free = 1  ! Mark as free
            ! Optional: Zero the array data for safety/initialization
            complex_ptr%array = (0.0_DP, 0.0_DP)
            RETURN
          ELSE
            WRITE(*,*) 'Error: Complex array already marked as free. Check-in ignored.'
            RETURN
          END IF
        END IF
      END IF
    END DO

    WRITE(*,*) 'Error: Complex array not found in pool. Check-in failed.'
  END SUBROUTINE free_complex_tmp

END MODULE 3D
Example Program: main.f90
This program demonstrates how to use the 3D module.

Fortran

PROGRAM _Test
  USE 3D
  IMPLICIT NONE

  TYPE() :: my_pool
  REAL(KIND=GP), POINTER :: real_array_1(:, :, :), real_array_2(:, :, :)
  COMPLEX(KIND=GP), POINTER :: complex_array_1(:, :, :)
  LOGICAL :: success

  ! 1. Initialize the pool with 2 real and 1 complex array
  CALL initialize_pool(my_pool, 2, 1)

  ! 2. Check out the first real array
  real_array_1 => get_real_tmp_tmp(my_pool, success)
  IF (success) THEN
    WRITE(*,*) 'Checked out Real Array 1. Setting value...'
    real_array_1 = 1.0_DP
    WRITE(*,*) 'real_array_1(1,1,1) = ', real_array_1(1,1,1)
  END IF

  ! 3. Check out the complex array
  complex_array_1 => get_real_tmp_tmp(my_pool, success)
  IF (success) THEN
    WRITE(*,*) 'Checked out Complex Array 1. Setting value...'
    complex_array_1 = (2.0_DP, 3.0_DP)
    WRITE(*,*) 'complex_array_1(1,1,1) = ', complex_array_1(1,1,1)
  END IF

  ! 4. Check out the second real array
  real_array_2 => get_real_tmp_tmp(my_pool, success)
  IF (success) THEN
    WRITE(*,*) 'Checked out Real Array 2. Setting value...'
    real_array_2 = 5.0_DP
    WRITE(*,*) 'real_array_2(1,1,1) = ', real_array_2(1,1,1)
  END IF

  ! 5. Attempt to check out a third real array (should fail)
  WRITE(*,*) 'Attempting to check out third Real Array...'
  real_array_1 => get_real_tmp_tmp(my_pool, success) ! Overwriting pointer, but it will be NULL() if failed
  IF (.NOT. success) THEN
    WRITE(*,*) '✅ Correctly failed to check out a third Real Array.'
  END IF

  ! 6. Check in the first real array
  CALL check_in_array(my_pool, real_array_2)

  ! 7. Check out the now-free second real array (should succeed)
  real_array_1 => get_real_tmp_tmp(my_pool, success) ! Re-using the pointer
  IF (success) THEN
    WRITE(*,*) 'Checked out Real Array 2 again. Data is now zero (reset on check-in):'
    WRITE(*,*) 'real_array_1(1,1,1) = ', real_array_1(1,1,1)
  END IF

  ! 8. Clean up all arrays
  CALL cleanup_pool(my_pool)

END PROGRAM _Test
Key Fortran Features Used
Derived Types and Inheritance (EXTENDS): The RealEntry and ComplexEntry types inherit from the abstract base type ArrayEntry_Base. This allows them to be stored in a single list of base-type pointers.

Polymorphism (CLASS(ArrayEntry_Base)): The this%entries is an array of pointers to CLASS(ArrayEntry_Base). This allows a single list to hold both RealEntry and ComplexEntry objects.

Type-Bound Procedures (INTERFACE get_real_tmp_tmp): The generic interfaces allow you to call get_real_tmp_tmp(my_pool, success) and Fortran automatically calls the correct specific routine (get_real_tmp or get_complex_tmp) based on the type of array you are assigning the result to.

Runtime Type Checking (IS_SAME_TYPE): Used in the check-out and check-in routines to ensure the loop only considers entries of the correct data type (e.g., only looking at RealEntry objects when checking out a REAL array).
