! ===================================================================
! Factory for all PDEs.
!
! To add a new solver, "USE" the corresponding module, and add a new
! case clause in the init_pdes_from_file function.
!
! DATE : 01/13/26 (PDM)
! ===================================================================

module equation_factory
  USE equationbase_mod
  USE hd_mod
! USE userdefinedpdes_mod
  
  IMPLICIT NONE

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
!     case ('UserDefined')
!       allocate(UserDefinedsolver :: new_object)
     case default
        stop 'Unknown solver name'
    end select
  end function init_pdes_from_file

end module equation_factory
