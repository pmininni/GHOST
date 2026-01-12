module gstate_mod
  USE fprecision
  USE grid
  USE mpivars

  IMPLICIT NONE

  ! ================= Base field data types =========================
  ! Derived types for complex (Fourier transformed) field components
  type, public :: GState
    complex(kind=GP), allocatable :: ccomp(:,:,:) ! complex component      
  end type GState
  
  ! Derived type for real space field components
  type, public :: GStateReal
    complex(kind=GP), allocatable :: rcomp(:,:,:) ! real component
  end type GStateReal

end module gstate_mod
