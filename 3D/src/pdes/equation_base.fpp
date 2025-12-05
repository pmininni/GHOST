! ===================================================================
! NAME       : equation_base.fpp
! DESCRIPTION: Forms base class for all PDEs
! DATE       : 11/29/25 (DLR)
! ===================================================================
module equationbase_mod
    implicit none

    ! Define an abstract base class
    type, abstract :: EquationBase
    contains
        procedure         (step), deferred ::          step ! step method
          generic :: step           => step_impl
        procedure     (tmp_size), deferred ::      tmp_size ! tmp size
          generic :: tmp_size       => tmp_size_impl
        procedure   (state_size), deferred ::    state_size ! state size
          generic :: state_size     => state_size_impl
        procedure(sstate2istate), deferred :: sstate2istate ! state names
          generic :: sstate2istate  => sstate2istate_impl
        procedure   (get_sstate), deferred ::    get_sstate ! get list of state names
          generic :: get_sstate     => get_sstate_impl
        procedure         (init), deferred ::          init ! init method
          generic :: init           => get_sstate_impl
        procedure         (dudt), deferred ::          dudt ! RHS method
          generic :: dudt           => dudt_impl
    end type EquationBase

    ! Abstract interface for deferred (virtual) methods:
    abstract interface
    ! NOTE: DO WE NEED THESE?
        subroutine step(this, time, uin, uf, dt, workspace, uout) 
            import :: EquationBase
            class (EquationBase), intent   (in) :: this
            real    (kind=GP), intent   (in) :: time, dt
            type (GStateComp), intent(inout) :: uin(:)
            type (GWorkspace), intent(inout) :: workspace
            type (GStateComp), intent   (in) :: uout(:)
        end subroutine step

        subroutine sstate2istate(this, sstate, istate) 
            import :: EquationBase
            class (EquationBase), intent   (in) :: this
            character (len=:), intent   (in) :: sstate(:)
            integer          , allocatable &
                             , intent(inout) :: istate(:)
        end subroutine sstate2istate

        subroutine get_sstate(this, sstate) 
            import :: EquationBase
            class (EquationBase), intent   (in) :: this
            class       (pde), intent   (in) :: this
            character (len=:), allocatable, &
                             , intent(inout) :: sstate(:)
        end subroutine get_sstate

        function state_size(this) result(num)
            import :: EquationBase
            class (EquationBase), intent   (in) :: this
            integer :: num
        end function state_size

        function tmp_size(this) result(num)
            import :: EquationBase
            class (EquationBase), intent   (in) :: this
            integer :: num
        end function tmp_size

        subroutine init(this) 
            import :: EquationBase
            class (EquationBase), intent   (in) :: this
        end subroutine init

        subroutine dudt(this, time, uin, uf, dt, dudt) 
            import :: pde
            class       (pde), intent   (in) :: this
            real    (kind=GP), intent   (in) :: time, dt
            type (GStateComp), intent(inout) :: uin(:)
            type (GWorkspace), intent(inout) :: workspace
            type (GStateComp), intent   (in) :: dudt(:) 
        end subroutine dudt
    end interface

contains

    ! Concrete method available to all derived classes here:

end module equationbase_mod


