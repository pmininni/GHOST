! ===================================================================
! NAME       : neutral_hydro.fpp
! DESCRIPTION: Forms class for all neutral hyrdo PDEs
! DATE       : 11/30/25 (DLR)
! ===================================================================
module neutralhydro_mod
    use equationbase_mod

    implicit none

    ! Define an abstract base class
    type, extends(EquationBase) :: NeutralHydro 
    contains
        type, public  :: NHTraits
          integer       :: dorot        = 0     ; ! rotation flag
          integer       :: dobouss      = 0     ; ! Boussinesq flag
          integer       :: domoistbouss = 0     ; ! moist Boussinesq flag
          integer       :: numpassive   = 0     ; ! num passive scalars
          real(kind=GP) :: nu           = 0.0_GP; ! dissipation
          real(kind=GP) :: kappa        = 0.0_GP; ! temp/den diffisuvity 
          real(kind=GP), allocatable :: passive_diff(:); ! diffusivities
        end type

        procedure,public  ::  NeutralHydro_ctor ! constructor
        final             ::  NeutralHydro_dtor ! desutructor
        procedure,public  ::          init_impl ! init method
        procedure,public  ::          step_impl ! step method
        procedure,public  ::          dudt_impl ! RHS method
        procedure,public  ::      tmp_size_impl ! tmp size
        procedure,public  ::    state_size_impl ! state size
        procedure,public  :: sstate2istate_impl ! state names
        procedure,public  ::    get_sstate_impl ! get list of state names
    end type NeutralHydro

    contains

    !! Concrete method available to all derived classes here:
    subroutine NeutralHydro_ctor(this, traits, workspace) 
      class        (NeutralHydro), intent   (in) :: this
      real    (kind=GP), intent   (in) :: time, dt
      type   (NHTraits), intent   (in) :: traits
      type         (GWorkspace), intent(inout) :: workspace
      type (GStateComp), intent   (in) :: uout(:)
    end subroutine NeutralHydro_ctor
    subroutine step(this, time, uin, uf, dt, workspace, uout) 
      class       (NeutralHydro), intent   (in) :: this
      real    (kind=GP), intent   (in) :: time, dt
      type (GStateComp), intent(inout) :: uin(:)
      type (GWorkspace), intent(inout) :: workspace
      type (GStateComp), intent   (in) :: uout(:)
    end subroutine step

    subroutine sstate2istate(this, sstate, istate) 
      class       (NeutralHydro), intent   (in) :: this
      character (len=:), intent   (in) :: sstate(:)
      integer          , allocatable &
                             , intent(inout) :: istate(:)
    end subroutine sstate2istate

    subroutine get_sstate(this, sstate) 
      class       (NeutralHydro), intent   (in) :: this
      character (len=:), allocatable, &
                             , intent(inout) :: sstate(:)
    end subroutine get_sstate

    !! Function to compute number of state members (equations)
    function state_size(this) result(num)
      class    (NeutralHydro), intent   (in) :: this
      integer :: num
      !   integer       :: dorot        = 0     ; ! rotation flag
      !   integer       :: dobouss      = 0     ; ! Boussinesq flag
      !   integer       :: domoistbouss = 0     ; ! moist Boussinesq flag
      !   integer       :: numpassive   = 0     ; ! num passive scalars
      num = 0;
      if ( dbouss .eq. 1 ) num = num + 1
    end function state_size

    function tmp_size(this) result(num)
      class    (NeutralHydro), intent   (in) :: this
      integer :: num
    end function tmp_size

    subroutine init(this) 
      class    (NeutralHydro), intent   (in) :: this
    end subroutine init

    subroutine dudt(this, time, uin, uf, dt, dudt) 
      class       (NeutralHydro), intent   (in) :: this
      real    (kind=GP), intent   (in) :: time, dt
      type (GStateComp), intent(inout) :: uin(:)
      type (GWorkspace), intent(inout) :: workspace
      type (GStateComp), intent   (in) :: dudt(:) 
    end subroutine dudt


end module neutralhydro_mod


