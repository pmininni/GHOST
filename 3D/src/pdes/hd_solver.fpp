! ===================================================================
! NAME       : hd_solver.fpp
! DESCRIPTION: Forms class for incompressible HD solver, computing
!
!                dv/dt + v.Grad v 2*Omega x v = -Grad p + nu Del^2 v
!
!                d s_i + v. Grad s_i = kappa_i Del^2 s_i, 
!                                  i = 1, ..., numpassive
!              State ordering is:
!                v1, v2, v3, s1, s2, ..., s_numpassive
!
!              State sector ids are:
!                VELOCITY (VELOCITY+1, VELOCITY+2)
!                PASSIVE  (PASSSIVE+1, PASSIVE+2, ...)
!
! DATE       : 11/30/25 (DLR)
! ===================================================================
module neutralhydro_mod
    use equationbase_mod

    implicit none

    ! Define an abstract base class
    type, extends(EquationBase) :: HDSolver 
        type, public  :: NHTraits
          integer       :: dorot        = 0     ; ! rotation flag
          integer       :: numpassive   = 0     ; ! num passive scalars
          real(kind=GP) :: nu           = 0.0_GP; ! dissipation
          real(kind=GP), allocatable :: passive_diff(:); ! diffusivities
          real(kind=GP), allocatable :: omega( :);! rotation vector
        end type

        ! Member data:
        type (GWorkspace), pointer   :: workspace_
        type   (NHTraits)            :: traits_
        

    contains
        procedure,public  ::      HDSolver_ctor ! constructor
        final             ::      HDSolver_dtor ! desutructor
        procedure,public  ::          init_impl ! init method
        procedure,public  ::          step_impl ! step method
        procedure,public  ::          dudt_impl ! RHS method
        procedure,public  ::      tmp_size_impl ! tmp size
        procedure,public  ::    state_size_impl ! state size
        procedure,public  :: sstate2istate_impl ! state names
        procedure,public  ::    get_sstate_impl ! get list of state names

        procedure,private ::    rhs_norot       ! RHS with no rot
        procedure,private ::    rhs_rot         ! RHS with rot
        procedure,private ::    rhs_passive     ! RHS for passive scalars
    end type HDSolver

    contains

    !! Concrete method available to all derived classes here:
    subroutine HDSolver_ctor(this, traits, workspace) 
      class  (HDSolver), intent   (in) :: this
      real    (kind=GP), intent   (in) :: time, dt
      type   (NHTraits), intent   (in) :: traits
      type         (GWorkspace), intent(inout), target
                                       :: workspace

      this%workspace_ => workspace
      this%traits_.
    end subroutine HDSolver_ctor


    !! Concrete method to take one time step
    subroutine step(this, time, uin, uf, dt, uout) 
      class  (HDSolver), intent   (in) :: this
      real    (kind=GP), intent   (in) :: time, dt
      type (GStateComp), intent(inout) :: uin(:)
      type (GWorkspace), intent(inout) :: workspace
      type (GStateComp), intent   (in) :: uout(:)
    end subroutine step

    subroutine sstate2istate(this, sstate, istate) 
      class  (HDSolver), intent   (in) :: this
      character (len=:), intent   (in) :: sstate(:)
      integer          , allocatable &
                             , intent(inout) :: istate(:)
    end subroutine sstate2istate

    subroutine get_sstate(this, sstate) 
      class  (HDSolver), intent   (in) :: this
      character (len=:), allocatable, &
                             , intent(inout) :: sstate(:)
    end subroutine get_sstate

    !! Function to compute number of state members (equations)
    function state_size(this) result(num)
      class  (HDSolver), intent   (in) :: this
      integer :: num
      !   integer       :: dorot        = 0     ; ! rotation flag
      !   integer       :: dobouss      = 0     ; ! Boussinesq flag
      !   integer       :: domoistbouss = 0     ; ! moist Boussinesq flag
      !   integer       :: numpassive   = 0     ; ! num passive scalars
      num = 0;
      if ( dbouss .eq. 1 ) num = num + 1
    end function state_size

    !! Function to compute number of tmp arrays
    function tmp_size(this) result(num)
      class  (HDSolver), intent   (in) :: this
      integer :: num
    end function tmp_size

    !! Function to initialize solver
    subroutine init(this) 
      class  (HDSolver), intent   (in) :: this
    end subroutine init

    !! Function to compute RHS
    subroutine dudt(this, time, uin, uf, dt, dudt) 
      class  (HDSolver), intent   (in) :: this
      real    (kind=GP), intent   (in) :: time, dt
      type (GStateComp), intent(inout) :: uin(:)
      type (GWorkspace), intent(inout) :: workspace
      type (GStateComp), intent   (in) :: dudt(:) 
    end subroutine dudt


end module neutralhydro_mod


