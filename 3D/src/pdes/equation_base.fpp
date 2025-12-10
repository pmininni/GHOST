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

contains

    ! Concrete method available to all derived classes here:
    ! Abstract interface for deferred (virtual) methods:
    ! NOTE: do WE NEED THESE?:
    abstract interface
        subroutine step(this, time, uin, uf, dt, workspace, uout) 
            import :: EquationBase
            class (EquationBase), intent   (in) :: this
            real       (kind=GP), intent   (in) :: time, dt
            type        (GState), intent(inout) :: uin(:)
            type    (GWorkspace), intent(inout) :: workspace
            type        (GState), intent(inout) :: uout(:)
        end subroutine step

        subroutine sstate2istate(this, sstate, istate) 
            import :: EquationBase
            class (EquationBase), intent   (in) :: this
            character (len=:)   , intent   (in) :: sstate(:)
            integer             , allocatable &
                                , intent(inout) :: istate(:)
        end subroutine sstate2istate

        subroutine get_sstate(this, sstate) 
            import :: EquationBase
            class (EquationBase), intent   (in) :: this
            class (EquationBase), intent   (in) :: this
            character (len=:)   , allocatable, &
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
            import :: EquationBase
            class (EquationBase), intent   (in) :: this
            real       (kind=GP), intent   (in) :: time, dt
            type        (GState), intent(inout) :: uin(:)
            type    (GWorkspace), intent(inout) :: workspace
            type        (GState), intent   (in) :: dudt(:) 
        end subroutine dudt
    end interface

  contains

        !  Concrete method to compute RHS for all passive scalars
        subroutine rhs_passive(this, uin, f, &
                   kappa, ivstart, numv, isstart, npassive, dudt)
!-----------------------------------------------------------------
!
! Compute npassive RHS's for passive scalar advection
!
! Parameters
!     uin     : current full state
!     f       : forces for each state comp
!     kappa   : diffusivities (must be npassive of these)
!     ivstart : start index of velocity sector in uin
!     numv    : number of velocity comps to use for advection
!     sstart  : start index of scalar sector
!     npassive: number of passive scalars to advect
!     dudt    : computed RHS for each scalar, 1-npassive

!*****************************************************************
            use fprecision
            use ali
            use kes
            use var
            use grid
            use mpivars
!$          use threads
            implicit none

            class (EquationBase), intent   (in) :: this
            integer          , intent   (in) :: isstart, npassive
            integer          , intent   (in) :: ivstart, numv
            real    (kind=GP), intent   (in) :: time, dt
            real    (kind=GP), intent   (in) :: kappa(:)
            type     (GState), intent(inout) :: uin(:)
            type     (GState), intent(inout) :: f(:)
            type (GWorkspace), intent(inout), pointer &
                                             :: workspace
            type     (GState), intent  (out) :: dudt(:) 

            logical                          :: bret
            complex, pointer, dimension(nz,ny,ista:iend)::ctmp
 
            ctmp = workspace=>get_complex_tmp(bret)
            if ( .not.bret ) then
              stop 'EquationBase::rhs_passive: workspace get failure'
            endif 

            if ( numv .eq. 3 ) then ! 3d advection
              do n = 1, npassive
                call advect3(uin(ivstart),uin(kvstart+1),uin(ivstart+2),&
                             uin(isstart+n-1),ctmp)
                call laplak3(uin(isstart+n-1),uin(isstart+n-1))
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
                do i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
                do j = 1,ny
                do k = 1,nz
                   if ((kn2(k,j,i).le.kmax).and.(kn2(k,j,i).ge.tiny)) then
                     dudt(n)=>ccomp(k,j,i) = &
                      kappa(n)*uin(isstart+n-1)(k,j,i)+ctmp(k,j,i) &
                     +f(isstart)=>ccomp(k,j,i)
                   else if (kn2(k,j,i).gt.kmax) then
                      dudt(n)=>ccomp(k,j,i) = 0.0_GP
                   else if (kn2(k,j,i).lt.tiny) then
                      dudt(n)=>ccomp(k,j,i) = 0.0_GP
                   endif
                enddo
                enddo
                enddo
              enddo ! end, loop over scalars
            else ! 2d advection:
              do n = 1, npassive
                call advect3(uin(ivstart),uin(kvstart+1),uin(ivstart+2),&
                             uin(isstart+n-1),ctmp)
                call laplak3(uin(isstart+n-1),uin(isstart+n-1))
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
                do i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
                do j = 1,ny
                do k = 1,nz
                   if ((kn2(k,j,i).le.kmax).and.(kn2(k,j,i).ge.tiny)) then
                     dudt(n)=>ccomp(k,j,i) = &
                      kappa(n)*uin(isstart+n-1)(k,j,i)+ctmp(k,j,i) &
                     +f(isstart)=>ccomp(k,j,i)
                   else if (kn2(k,j,i).gt.kmax) then
                      dudt(n)=>ccomp(k,j,i) = 0.0_GP
                   else if (kn2(k,j,i).lt.tiny) then
                      dudt(n)=>ccomp(k,j,i) = 0.0_GP
                   endif
                enddo
                enddo
                enddo
              enddo ! end, loop over scalars
            endif

            workspace=>free_complex_tmp(ctmp)

        end subroutine rhs_passive

end module equationbase_mod


