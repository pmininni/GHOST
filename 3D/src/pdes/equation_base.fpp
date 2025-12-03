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
        procedure         (step), deferred ::          step_impl ! step method
        procedure     (tmp_size), deferred ::      tmp_size_impl ! tmp size
        procedure   (state_size), deferred ::    state_size_impl ! state size
        procedure(sstate2istate), deferred :: sstate2istate_impl ! state names
        procedure   (get_sstate), deferred ::    get_sstate_impl ! get list of state names
        procedure         (init), deferred ::          init_impl ! init method
        procedure         (dudt), deferred ::          dudt_impl ! RHS method
    end type EquationBase

    ! Abstract interface for deferred (virtual) methods:
    abstract interface
        subroutine step(self, time, uin, uf, dt, workspace, uout) 
            import :: pde
            class       (pde), intent   (in) :: self
            real    (kind=GP), intent   (in) :: time, dt
            type (GStateComp), intent(inout) :: uin(:)
            type (GWorkspace), intent(inout) :: workspace
            type (GStateComp), intent   (in) :: uout(:)
        end subroutine step

        subroutine sstate2istate(self, sstate, istate) 
            import :: pde
            class       (pde), intent   (in) :: self
            character (len=:), intent   (in) :: sstate(:)
            integer          , allocatable &
                             , intent(inout) :: istate(:)
        end subroutine sstate2istate

        subroutine get_sstate(self, sstate) 
            import :: pde
            class       (pde), intent   (in) :: self
            character (len=:), allocatable, &
                             , intent(inout) :: sstate(:)
        end subroutine get_sstate

        function state_size(self) result(num)
            import :: pde
            class    (pde), intent   (in) :: self
            integer :: num
        end function state_size

        function tmp_size(self) result(num)
            import :: pde
            class    (pde), intent   (in) :: self
            integer :: num
        end function tmp_size

        subroutine init(self) 
            import :: pde
            class    (pde), intent   (in) :: self
        end subroutine init

        subroutine dudt(self, time, uin, uf, dt, dudt) 
            import :: pde
            class       (pde), intent   (in) :: self
            real    (kind=GP), intent   (in) :: time, dt
            type (GStateComp), intent(inout) :: uin(:)
            type (GWorkspace), intent(inout) :: workspace
            type (GStateComp), intent   (in) :: dudt(:) 
        end subroutine dudt
    end interface

contains

    ! Concrete method available to all derived classes here:

end module equationbase_mod


