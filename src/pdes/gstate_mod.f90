module gstate_mod
  USE fprecision
  USE mpivars
  USE grid

  IMPLICIT NONE

  ! ================= Base field data types =========================
  ! Derived types for complex (Fourier transformed) field components
  type, public :: GStateComp
    complex(kind=GP), allocatable :: ccomp(:,:,:) ! complex component      
  end type GStateComp
  
  ! Derived type for real space field components
  type, public :: GStateCompReal
    real(kind=GP)   , allocatable :: rcomp(:,:,:) ! real component
  end type GStateCompReal

! type, public :: GState
!   type(GStateComp), allocatable :: comp(:) ! real component
! end type GState
  
! type, public :: GStateReal
!   type(GStateCompReal), allocatable :: comp(:) ! real component
! end type GStateReal
  

CONTAINS

  ! ================= Allocation routines ===========================

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
    ALLOCATE( state(nc) )
    DO i = 1,nc
      ALLOCATE( state(i)%ccomp(nz,ny,ista:iend) )
    END DO
  end subroutine GState_alloc

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Method to allocate real GState data types
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine GStateReal_alloc(state, nc)
    use grid
    use mpivars
    implicit none
    type (GStateCompReal), allocatable, intent(inout) :: state(:)
    integer                           , intent   (in) :: nc
    integer                                           :: i
    ALLOCATE( state(nc) )
    DO i = 1,nc
      ALLOCATE( state(i)%rcomp(nx,ny,ksta:kend) )
    END DO
  end subroutine GStateReal_alloc
  
end module gstate_mod
