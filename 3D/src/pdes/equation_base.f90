! ===================================================================
! NAME       : equation_base.fpp
! DESCRIPTION: Forms base class for all PDEs
! DATE       : 11/29/25 (DLR)
! ===================================================================

module equationbase_mod
  USE class_GWorkspace3D
  USE gstate_mod

  IMPLICIT NONE

  ! ================= Base class for all PDEs =======================
  ! Define an abstract base class
  type, abstract :: EquationBase
    contains
      procedure(Solver_ctor_interface), deferred ::  Solver_ctor ! Constructor
      procedure(init_interface),        deferred ::  init        ! init method
      procedure(dudt_interface),        deferred ::  dudt        ! RHS method
  end type EquationBase

  abstract interface
     subroutine Solver_ctor_interface(this, infile, workspace, nc)
       USE class_GWorkspace3D
       import :: EquationBase
       class(EquationBase), intent(inout)         :: this
       integer            , intent   (in)         :: nc
       type(GWorkspace)   , intent(inout), target :: workspace
       character(len=*)   , intent   (in)         :: infile  
     end subroutine Solver_ctor_interface

     subroutine init_interface(this) 
       import :: EquationBase
       class (EquationBase), intent (inout) :: this
     end subroutine init_interface

     subroutine dudt_interface(this, time, uin, uf, dt, dudt) 
       USE gstate_mod
       import :: EquationBase
       class(EquationBase), intent   (in)         :: this
       real      (kind=GP), intent   (in)         :: time, dt
       type       (GState), intent(inout), target :: uin(:),uf(:)
       type       (GState), intent(inout)         :: dudt(:) 
     end subroutine dudt_interface
  end interface
     
!!$        subroutine step(this, time, uin, uf, dt, workspace, uout) 
!!$            import :: EquationBase
!!$            class (EquationBase), intent   (in) :: this
!!$            real       (kind=GP), intent   (in) :: time, dt
!!$            type        (GState), intent(inout) :: uin(:)
!!$            type    (GWorkspace), intent(inout) :: workspace
!!$            type        (GState), intent(inout) :: uout(:)
!!$        end subroutine step
!!$
!!$        subroutine sstate2istate(this, sstate, istate) 
!!$            import :: EquationBase
!!$            class (EquationBase), intent   (in) :: this
!!$            character (len=:)   , intent   (in) :: sstate(:)
!!$            integer             , allocatable &
!!$                                , intent(inout) :: istate(:)
!!$        end subroutine sstate2istate
!!$
!!$        subroutine get_sstate(this, sstate) 
!!$            import :: EquationBase
!!$            class (EquationBase), intent   (in) :: this
!!$            class (EquationBase), intent   (in) :: this
!!$            character (len=:)   , allocatable, &
!!$                                , intent(inout) :: sstate(:)
!!$        end subroutine get_sstate
!!$
!!$        function state_size(this) result(num)
!!$            import :: EquationBase
!!$            class (EquationBase), intent   (in) :: this
!!$            integer :: num
!!$        end function state_size
!!$
!!$        function tmp_size(this) result(num)
!!$            import :: EquationBase
!!$            class (EquationBase), intent   (in) :: this
!!$            integer :: num
!!$        end function tmp_size
!!$

CONTAINS
  
  ! ==== Concrete method to compute RHS for all passive scalars =====
  subroutine rhs_passive(this, uin, f, &
                    kappa, ivstart, numv, isstart, npassive, dudt)
  !------------------------------------------------------------------
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
  !******************************************************************
            use fprecision
            use ali
            use kes
            use var
            use grid
            use mpivars
            USE gstate_mod
!$          use threads
            implicit none

            class (EquationBase), intent(in) :: this
            integer          , intent   (in) :: isstart, npassive
            integer          , intent   (in) :: ivstart, numv
!           real    (kind=GP), intent   (in) :: time, dt
            real    (kind=GP), intent   (in) :: kappa(:)
            type     (GState), intent(inout) :: uin(:)
            type     (GState), intent(inout) :: f(:)
            type (GWorkspace), pointer       :: workspace
            type     (GState), intent  (out) :: dudt(:) 

            logical                          :: bret
            integer                          :: i,j,k,n
            complex, pointer                 :: ctmp(:,:,:)
 
            CALL workspace%get_complex_tmp(ctmp,bret)
            if ( .not.bret ) then
              stop 'EquationBase::rhs_passive: workspace get failure'
            endif 

            if ( numv .eq. 3 ) then ! 3d advection
              do n = 1, npassive
                call advect3(uin(ivstart),uin(ivstart+1),uin(ivstart+2),&
                             uin(isstart+n-1),ctmp)
                call laplak3(uin(isstart+n-1),uin(isstart+n-1))
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
                do i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
                do j = 1,ny
                do k = 1,nz
                   if ((kn2(k,j,i).le.kmax).and.(kn2(k,j,i).ge.tiny)) then
                      dudt(n)%ccomp(k,j,i) = &
                           kappa(n)*uin(isstart+n-1)%ccomp(k,j,i)+ctmp(k,j,i) &
                     +f(isstart)%ccomp(k,j,i)
                   else if (kn2(k,j,i).gt.kmax) then
                      dudt(n)%ccomp(k,j,i) = 0.0_GP
                   else if (kn2(k,j,i).lt.tiny) then
                      dudt(n)%ccomp(k,j,i) = 0.0_GP
                   endif
                enddo
                enddo
                enddo
              enddo ! end, loop over scalars
            else ! 2d advection:
              do n = 1, npassive
                call advect3(uin(ivstart),uin(ivstart+1),uin(ivstart+2),&
                             uin(isstart+n-1),ctmp)
                call laplak3(uin(isstart+n-1),uin(isstart+n-1))
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
                do i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
                do j = 1,ny
                do k = 1,nz
                   if ((kn2(k,j,i).le.kmax).and.(kn2(k,j,i).ge.tiny)) then
                     dudt(n+isstart)%ccomp(k,j,i) = &
                      kappa(n)*uin(isstart+n-1)%ccomp(k,j,i)+ctmp(k,j,i) &
                     +f(isstart)%ccomp(k,j,i)
                   else if (kn2(k,j,i).gt.kmax) then
                      dudt(n+isstart)%ccomp(k,j,i) = 0.0_GP
                   else if (kn2(k,j,i).lt.tiny) then
                      dudt(n+isstart)%ccomp(k,j,i) = 0.0_GP
                   endif
                enddo
                enddo
                enddo
              enddo ! end, loop over scalars
            endif

            CALL workspace%free_complex_tmp(ctmp)

        end subroutine rhs_passive

end module equationbase_mod


