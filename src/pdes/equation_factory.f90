! ===================================================================
! Factory for all PDEs.
!
! To add a new solver, "USE" the corresponding module, and add a new
! case clause in the init_pdes_from_file function. The clause should
! allocate the solver class, and declare the number of field
! components the solver needs.
!
! DATE : 03/29/26 (JBG)
! ===================================================================

module equation_factory
  USE equationbase_mod
  USE hd_mod
  USE mhd_mod
! USE userdefinedpde_mod
  
  IMPLICIT NONE

  ! ================= Global parameters ===============================
  integer, public :: NUMTMPCOMP = 0 ! Number of cmplx tmp arrays
  integer, public :: NUMTMPREAL = 0 ! Number of real tmp arrays
  
CONTAINS
  
  ! ================= Factory function ==============================
  function init_pdes_from_file(infile) result(new_object)
    USE commtypes
    class(EquationBase), allocatable  :: new_object
    character(len=*)   , intent(in)   :: infile

    ! Temporary data to read from namelist:
    integer           :: nprocs,myrank,ierr
    character(len=64) :: solver
    ! Required namelist:
    namelist/ pdes / solver

    call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
    if ( myrank .eq. 0 ) then
      open(1,file=infile,status='unknown',form="formatted")
      read(1,NML=pdes)
      close(1)
    endif
    call MPI_BCAST(solver,64,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)

    ! Clauses for each solver class
    select case (trim(adjustl(solver)))
      case ('HD')
        allocate(HDsolver :: new_object)
        NUMTMPCOMP =  8; NUMTMPREAL = 3
      case ('BOUSS')
!       allocate(BOUSSsolver :: new_object)
!       NUMTMPCOMP =  8; NUMTMPREAL = 3
      case ('MOIST')
        allocate(MOISTsolver :: new_object)
        NUMTMPCOMP =  8; NUMTMPREAL = 3 ! TODO check
      case ('MHD')
        allocate(MHDsolver :: new_object)
        NUMTMPCOMP = 12; NUMTMPREAL = 3
!     case ('UserDefined')
!       allocate(UserDefinedsolver :: new_object)
!       NUMFIELDS = 3; NUMTMPCOMP =  8; NUMTMPREAL = 3
      case default
        stop 'Equation factory :: init_pdes_from_file : Unknown solver name'
    end select
  end function init_pdes_from_file

end module equation_factory
