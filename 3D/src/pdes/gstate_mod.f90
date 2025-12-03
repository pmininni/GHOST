module gstate_types
    use fprecision
    use mpivars

    implicit none

    type, public :: GStateComp 
      complex(kind=GP), pointer, dimension(:,:,:) :: ccomp ! complex component
      real   (kind=GP), pointer, dimension(:,:,:) :: rcomp ! real component
    end type GStateComp


end module gstate_types


