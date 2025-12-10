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
module hd_mod
    use equationbase_mod
    use gstate_mod

    implicit none

    ! Define an abstract base class
    type, extends(EquationBase) :: HDSolver 
        type, public  :: NHTraits
          integer       :: dorot        = 0     ; ! rotation flag
          integer       :: numpassive   = 0     ; ! num passive scalars
          real(kind=GP) :: nu           = 0.0_GP; ! dissipation
          real(kind=GP), allocatable :: passive_diff(:); ! diffusivities
          real(kind=GP), allocatable :: omega(3);! rotation vector
        end type

        ! Member data:
        logical                      :: binit_=.false. 
                                                 ! is initialized?
        integer                      :: VELOCITY ! start of velocity sector
        integer                      :: PASSIVE  ! start of scalar sector
        integer                      :: numpassive_ ! # passive scalars
        integer                      :: nd_      ! problem dimension
        integer                      :: nc_      ! # velocity components
        type (GWorkspace), pointer   :: workspace_
        type   (NHTraits)            :: traits_
        character(len=8), allocatable:: sstate_(:) 
        

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

        procedure,private ::    dudt_norot       ! RHS with no rot
        procedure,private ::    dudt_rot         ! RHS with rot
    end type HDSolver

    contains

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Constructor
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine HDSolver_ctor(this, traits, workspace, nc) 
      class  (HDSolver), intent   (in) :: this
      integer          , intent   (in) :: nc
      real    (kind=GP), intent   (in) :: time, dt
      type   (NHTraits), intent   (in) :: traits
      type         (GWorkspace), intent(inout), target
                                       :: workspace
      this%workspace_ => workspace
      this%nc_        = nc ! # vel. components (usefiul if ny=1)

      ! Make deep copy of traits:
      this%traits_%     dorot = traits%dorot
      this%traits_%numpassive = traits%numpassive
      this%traits_%        nu = traits%nu
      this%traits_%     omega = traits%omega
      if ( traits%numpassive .ne. size(traits%kappa) ) then
        stop 'HDSolver_ctor: inconsistent kappa array size'
      endif
      allocate(this%traits_%kappa(traits%numpassive))
      this%traits_%kappa = traits%kappa
      
    end subroutine HDSolver_ctor

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Destructor
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine HDSolver_dtor(this) 
      class  (HDSolver), intent   (in) :: this

      deallocate(this%traits_%kappa)
      deallocate(this%traits_%sstate_)

      ! If we use persistent workspace, free it here:
      ! this%workspace_%free_complex_tmp(ctmp1)
      ! this%workspace_%free_complex_tmp(ctmp2) ...
      
    end subroutine HDSolver_dtor


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Concrete method to take one time step
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine step_impl(this, time, uin, uf, dt, uout) 
      class  (HDSolver), intent   (in) :: this
      real    (kind=GP), intent   (in) :: time, dt
      type (GState), intent(inout) :: uin(:)
      type (GWorkspace), intent(inout) :: workspace
      type (GState), intent   (in) :: uout(:)

      if ( .not. binit_ ) then
        stop 'HDSolver::step_impl: Solver not initialized'
      endif

    end subroutine step_impl

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Convert input state name to index in state vector
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine sstate2istate_impl(this, sstate, istate) 
      class  (HDSolver), intent   (in) :: this
      character (len=8), intent   (in) :: sstate(:)
      integer          , allocatable &
                       , intent(inout) :: istate(:)
      integer                          :: i,j

      if ( size(sstate) .ne. size(istate) ) then
        stop 'HDSolver::sstate2istate_impl: Incompatible sstate and istate'
      endif
  
      do i = 1, size(sstate)
        istate(i) = -1 ! return unusable index
        do j = 1, size(this%sstate_)
          if ( sstate(i) .eq. this%sstate_(j) ) then
            istate(i) = j
          endif
        enddo
      enddo

    end subroutine sstate2istate_impl

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Get state variable names
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine get_sstate_impl(this, sstate) 
      class  (HDSolver), intent   (in) :: this
      character (len=8), allocatable, &
                       , intent(inout) :: sstate(:)
      character(len=100)               :: snum
      integer                          :: j
      if ( allocated(sstate) ) then
        deallocate(sstate);
        allocate(sstate(this%state_size()))
      endif
        do j = 1:this%nc_
           write(snum,'(I0)') j
           sstate(j) = 'v' // trim(snum)
        enddo
        do j = 1, this%traits_%numpassive
           write(snum,'(I0)') j
           sstate(j) = 's' // trim(snum)
        enddo
      endif
    end subroutine get_sstate_impl

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Function to compute number of state members (equations)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function state_size_impl(this) result(num)
      class  (HDSolver), intent   (in) :: this
      integer :: num

      num = this%nc_ ! # vel. components
      num = num + traits_%numpassive ! # scalars

    end function state_size_impl

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Function to compute number of tmp arrays
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function tmp_size_impl(this) result(num)
      class  (HDSolver), intent   (in) :: this
      integer :: num
    end function tmp_size_impl

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Function to initialize solver
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine init_impl(this) 
      class  (HDSolver), intent   (in) :: this

      this%VELOCITY = 1 ! start of vel sector
      if ( ny .eq. 1 ) then
        this%nd_ =  2 ! 2d
      else 
        this%nd_ =  3 ! 3d
      endif

      this%PASSIVE = this%VELOCITY + 1 ! start of scalar sector

      if ( allocated(this%sstate_) ) then
        deallocate(this%sstate_);
        allocate(this%sstate_(this%get_size()))
      endif
      this%get_sstate(this, this%sstate_)

      binit_ = .true.

    end subroutine init

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Function to compute RHS
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine dudt_impl(this, time, uin, uf, dt, dudt) 
      use fprecision
      use ali
      use kes
      use var
      use grid
      use mpivars
!$    use threads
      implicit none

      class  (HDSolver), intent   (in) :: this
      real    (kind=GP), intent   (in) :: time, dt
      type     (GState), intent(inout) :: uin(:)
      type (GWorkspace), intent(inout) :: workspace
      type     (GState), intent   (in) :: dudt(:) 

      if ( this%traits_%dorot ) then
        call dudt_rot  (this, time, uin, uf, dt, dudt) 
      else
        call dudt_norot(this, time, uin, uf, dt, dudt) 
      endif

    end subroutine dudt_impl


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Function to compute RHS when there's rotation
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine dudt_rot(this, time, uin, uf, dt, dudt) 
      use fprecision
      use ali
      use kes
      use var
      use grid
      use mpivars
!$    use threads
      implicit none

      class  (HDSolver), intent   (in) :: this
      real    (kind=GP), intent   (in) :: time, dt
      type     (GState), intent(inout) :: uin(:)
      type (GWorkspace), intent(inout) :: workspace
      type     (GState), intent   (in) :: dudt(:) 
      complex, pointer, dimension(nz,ny,ista:iend)
                                       :: fx,fy,fz,vx,vy,vz
      integer                          :: i,j,k
      real (kind=GP)                   :: nu,omegaz
      complex, pointer, dimension(nz,ny,ista:iend)
                                       :: c1,c2,c3,c4,c5,c6,c7,c8
       
      nu     = this%traits_.nu
      omegaz = this%traits_.omega(3)

      c1= this%workspace_=>get_complex_tmp(bret)
      c2= this%workspace_=>get_complex_tmp(bret)
      c3= this%workspace_=>get_complex_tmp(bret)
      c4= this%workspace_=>get_complex_tmp(bret)
      c5= this%workspace_=>get_complex_tmp(bret)
      c6= this%workspace_=>get_complex_tmp(bret)
      c7= this%workspace_=>get_complex_tmp(bret)
      c8= this%workspace_=>get_complex_tmp(bret)

      vx => uin  (this%VELOCITY)%ccomp
      vy => uin(this%VELOCITY+1)%ccomp
      vz => uin(this%VELOCITY+2)%ccomp

      fx => uf  (this%VELOCITY)%ccomp
      fy => uf(this%VELOCITY+1)%ccomp
      fz => uf(this%VELOCITY+2)%ccomp
      
      call prodre3(vx,vy,vz,C4,C5,C6)
      call saxpby_c(C4, 1.0_GP, vy, -2.0*omegaz) 
      call saxpby_c(C5, 1.0_GP, vx,  2.0*omegaz) 

      call nonlhd3(C4,C5,C6,C7,1)
      call nonlhd3(C4,C5,C6,C8,2)
      call nonlhd3(C4,C5,C6,C4,3)
      call laplak3(vx,vx)
      call laplak3(vy,vy)
      call laplak3(vz,vz)

!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      do i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
      do j = 1,ny
      do k = 1,nz
         if ((kn2(k,j,i).le.kmax).and.(kn2(k,j,i).ge.tiny)) then
            dudt  (this%VELOCITY)=>ccomp(k,j,i) = &
               nu*vx(k,j,i)+C7(k,j,i) + fx(k,j,i))
            dudt(this%VELOCITY+1)=>ccomp(k,j,i) = &
               nu*vy(k,j,i)+C8(k,j,i) + fy(k,j,i))
            dudt(this%VELOCITY+2)=>ccomp(k,j,i) = &
               nu*vz(k,j,i)+C4(k,j,i) + fz(k,j,i))
         else
            dudt  (this%VELOCITY)=>ccomp(k,j,i) = 0.0_GP
            dudt(this%VELOCITY+1)=>ccomp(k,j,i) = 0.0_GP
            dudt(this%VELOCITY+2)=>ccomp(k,j,i) = 0.0_GP
         endif
      enddo
      enddo
      enddo

      ! Compute passive scalars:
      call this%rhs_passibe(this, uin, uf, this%traits_%kappa, &
              this%VELOCITY, this%nc_, this%PASSIVE, &
              this%numpassive, dudt)

      this%workspace_=>free_complex_tmp(c1)
      this%workspace_=>free_complex_tmp(c2)
      this%workspace_=>free_complex_tmp(c3)
      this%workspace_=>free_complex_tmp(c4)
      this%workspace_=>free_complex_tmp(c5)
      this%workspace_=>free_complex_tmp(c6)
      this%workspace_=>free_complex_tmp(c7)
      this%workspace_=>free_complex_tmp(c8)

    end subroutine dudt_rot


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Function to compute RHS when there's no rotation
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine dudt_norot(this, time, uin, uf, dt, dudt) 
      use fprecision
      use ali
      use kes
      use var
      use grid
      use mpivars
!$    use threads
      implicit none

      class  (HDSolver), intent   (in) :: this
      real    (kind=GP), intent   (in) :: time, dt
      type     (GState), intent(inout) :: uin(:)
      type (GWorkspace), intent(inout) :: workspace
      type     (GState), intent   (in) :: dudt(:) 
      complex, pointer, dimension(nz,ny,ista:iend)
                                       :: fx,fy,fz,vx,vy,vz
      integer                          :: i,j,k
      real (kind=GP)                   :: nu
      complex, pointer, dimension(nz,ny,ista:iend)
                                       :: c1,c2,c3,c4,c5,c6,c7,c8
       
      nu = this%traits_.nu

      c1= this%workspace_=>get_complex_tmp(bret)
      c2= this%workspace_=>get_complex_tmp(bret)
      c3= this%workspace_=>get_complex_tmp(bret)
      c4= this%workspace_=>get_complex_tmp(bret)
      c5= this%workspace_=>get_complex_tmp(bret)
      c6= this%workspace_=>get_complex_tmp(bret)
      c7= this%workspace_=>get_complex_tmp(bret)
      c8= this%workspace_=>get_complex_tmp(bret)

      vx => uin  (this%VELOCITY)%ccomp
      vy => uin(this%VELOCITY+1)%ccomp
      vz => uin(this%VELOCITY+2)%ccomp

      fx => uf  (this%VELOCITY)%ccomp
      fy => uf(this%VELOCITY+1)%ccomp
      fz => uf(this%VELOCITY+2)%ccomp
      
      call prodre3(vx,vy,vz,C4,C5,C6)
      call nonlhd3(C4,C5,C6,C7,1)
      call nonlhd3(C4,C5,C6,C8,2)
      call nonlhd3(C4,C5,C6,C4,3)
      call laplak3(vx,vx)
      call laplak3(vy,vy)
      call laplak3(vz,vz)

!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      do i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
      do j = 1,ny
      do k = 1,nz
         if ((kn2(k,j,i).le.kmax).and.(kn2(k,j,i).ge.tiny)) then
            dudt  (this%VELOCITY)=>ccomp(k,j,i) = &
               nu*vx(k,j,i)+C7(k,j,i) + fx(k,j,i))
            dudt(this%VELOCITY+1)=>ccomp(k,j,i) = &
               nu*vy(k,j,i)+C8(k,j,i) + fy(k,j,i))
            dudt(this%VELOCITY+2)=>ccomp(k,j,i) = &
               nu*vz(k,j,i)+C4(k,j,i) + fz(k,j,i))
         else
            dudt  (this%VELOCITY)=>ccomp(k,j,i) = 0.0_GP
            dudt(this%VELOCITY+1)=>ccomp(k,j,i) = 0.0_GP
            dudt(this%VELOCITY+2)=>ccomp(k,j,i) = 0.0_GP
         endif
      enddo
      enddo
      enddo

      ! Compute passive scalars:
      call this%rhs_passibe(this, uin, uf, this%traits_%kappa, &
              this%VELOCITY, this%nc_, this%PASSIVE, &
              this%numpassive, dudt)


      this%workspace_=>free_complex_tmp(c1)
      this%workspace_=>free_complex_tmp(c2)
      this%workspace_=>free_complex_tmp(c3)
      this%workspace_=>free_complex_tmp(c4)
      this%workspace_=>free_complex_tmp(c5)
      this%workspace_=>free_complex_tmp(c6)
      this%workspace_=>free_complex_tmp(c7)
      this%workspace_=>free_complex_tmp(c8)

    end subroutine dudt_norot


end module hd_mod


