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

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Function to initialize solver
    !! Parameters:
    !!    infile   : input config file (e.g., 'parameter.inp')
    !!    workspace: workspace object
    !!    gsolver  : (output) solver/pde object
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine EquationFactory_build(infile, workspace, gsolver)
      use gutils
      use fprecision
      use ali
      use kes
      use var
      use grid 
      use mpivars
      !  *** NOTE: let's use only those modules needed!!! ***

      implicit none

      character(len=*)   , intent(in)  :: infile
      class(EquationBase), allocatable,
                           intent(out) :: gsolver
      type   (GWorkspace), intent(inout), target
                                       :: workspace
 
      integer                          :: ierr, myrank, nc
      character(len=128)               :: solver_name
      character(len=len(solver_name))  :: lc  ! lower case name

      namelist/ pde    / nc, solver_name

      call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)

      ! Get solver name, other data:
      if ( ny .eq. 1 ) then      ! is 2d
        nc = 2
      else if ( ny .gt. 1 ) then ! is 3d
        nc = 3
      endif
      if ( myrank .eq. 0 ) then
         open(1,file=infile,status='unknown',form="formatted")
         read(1,NML=pde)
         close(1)
         call mpi_bcast(solver_name ,128 ,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
         ! If is 2d, solver may have 2 or 3 vector components:
         if ( ny .eq. 1 ) then
           if ( nc .lt. 2 .or. nc .gt. 3 ) then
             stop 'equation_factory::build: invalid # components, nc'
           endif
           call mpi_bcast(nc          ,1   ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
         endif
      endif
      
      ! Deallocate gsolver if necessary:
      if ( allocated(gsolver) ( then
        deallocate(gsolver)
      endif
      lc = trim(adjust(solver_name))

      ! Enter solver selection list:
      select case ( lc ) ! select solver object

        case ('hd')
          allocate(HDSolver :: gsolver)
          gsolver%HDSolver_ctor(infile, workspace, nc)
        case ('roth')
          allocate(HDSolver :: gsolver)
          gsolver%HDSolver_ctor(infile, workspace, nc)
        case ('phd')
          allocate(HDSolver :: gsolver)
          gsolver%HDSolver_ctor(infile, workspace, nc)
        case ('mphd')
          allocate(HDSolver :: gsolver)
          gsolver%HDSolver_ctor(infile, workspace, nc)
        case ('proth')
          allocate(HDSolver :: gsolver)
          gsolver%HDSolver_ctor(infile, workspace, nc)
        case ('mproth')
          allocate(HDSolver :: gsolver)
          gsolver%HDSolver_ctor(infile, workspace, nc)
!       case ('bouss')
!         allocate(BoussSolver :: gsolver)
!         gsolver%BoussSolver_ctor(infile, workspace, nc)
!       case ('rotbouss')
!         allocate(BoussSolver :: gsolver)
!         gsolver%BoussSolver_ctor(infile, workspace, nc)
!       case ('mprotbouss')
!         allocate(BoussSolver :: gsolver)
!         gsolver%BoussSolver_ctor(infile, workspace, nc)
        case default
          resturn

      end select ! end, all solvers

    end subroutine EquationFactory_build


end module equation_factory_mod


