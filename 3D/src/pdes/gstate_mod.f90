module gstate_mod
    use fprecision
    use ali
    use kes
    use var
    use grid
    use mpivars
!$  use threads

    implicit none

    type, public :: GStateComp 
      complex(kind=GP), pointer, dimension(nz,ny,ista:iend) :: ccomp ! complex component
      real   (kind=GP), pointer, dimension(nx,ny,ksta_kend) :: rcomp ! real component
    end type GStateComp

  contains

end module gstate_mod


