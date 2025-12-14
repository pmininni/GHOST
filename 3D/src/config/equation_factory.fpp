! ===================================================================
! NAME       : equation_factory.fpp
! DESCRIPTION: Forms factory code for all PDEs
! DATE       : 12/12/25 (DLR)
! ===================================================================
module equation_factory_mod
    use hd_solver_mod

    implicit none

    type :: EquationFactory
    contains
        procedure, public build => EquationFactor_build
    end type EquationBase

contains

    subroutine EquationFactory_build(solver_name, nc, workspace, gsolver)
      use gutils
      implicit none

      character(len=*)   , intent(in)  :: solver_name
      integer            , intent(in)  :: nc
      class(EquationBase), allocatable,
                           intent(out) :: gsolver
      type   (GWorkspace), intent(inout), target
                                       :: workspace
 
      character(len=len(solver_name))  :: lc  ! lower case name
      if ( allocated(gsolver) ( then
        deallocate(gsolver)
      endif

      lc = trim(adjust(solver_name))

      select case ( lc ) ! select solver object

        case ('hd')
          allocate(HDSolver :: gsolver)
          gsolver%HDSolver_ctor(traits, workspace, nc)
        case ('roth')
          allocate(HDSolver :: gsolver)
          gsolver%HDSolver_ctor(traits, workspace, nc)
!       case ('bouss')
!         allocate(BoussSolver :: gsolver)
!         gsolver%BoussSolver_ctor(traits, workspace, nc)
!       case ('rotbouss')
!         allocate(BoussSolver :: gsolver)
!         gsolver%BoussSolver_ctor(traits, workspace, nc)
        case default
          resturn

      end select ! end, all solvers

    end subroutine EquationFactory_build


end module equation_factory_mod


