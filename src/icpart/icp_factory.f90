! ===================================================================
! Factory for all initial conditions of particles (ICpart)
!
! This factory applies chain operation to create multiple initial
! conditions for all state components in a particle solver.
!
! DATE : 03/21/26 (PDM)
! ===================================================================

module icp_factory
  USE icpbase_mod
  USE icp_position
  USE icp_velocity
  
  IMPLICIT NONE
  
CONTAINS
  
  ! ================= Factory function ==============================
  function init_icp_from_file(infile) result(new_object)
    USE gpstate_mod
    USE commtypes
    USE gutils
    type  (icpChain), allocatable :: new_object(:)
    character(len=*),  intent(in) :: infile

    ! Temporary data to read from namelist:
    integer            :: i,ib,ie,ind,icnmb
    character(len=256) :: iclist
    character(len= 25) :: icname(10)
    ! Required namelist:
    namelist/ part_initial_conditions / iclist

    if ( myrank .eq. 0 ) then
      open(1,file=infile,status='unknown',form="formatted")
      read(1,NML=part_initial_conditions)
      close(1)
    endif
    call MPI_BCAST(iclist,256,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)

    ! Loops through all requested ICs and allocates classes
    call parsestr(iclist,';',icname,25,10,icnmb)
    if ( myrank .eq. 0 ) print *,'Found',icnmb,'particle ICs'
    allocate( new_object(icnmb) )
    do i = 1,icnmb
      if ( myrank .eq. 0 ) print *,'Particle IC',i,':',trim(adjustl(icname(i)))
      select case (trim(adjustl(icname(i))))
      ! Position ICs ---------------------
      case ('read_x')
        allocate( read_pos    :: new_object(i)%icp )
      case ('user_x')
        allocate( user_pos    :: new_object(i)%icp )
      case ('random_x')
        allocate( rand_pos    :: new_object(i)%icp )
      ! Velocity ICs ---------------------    
      case ('read_v')
        allocate( read_vel    :: new_object(i)%icp )
      case ('null_v')
        allocate( null_vel    :: new_object(i)%icp )
      case ('fluid_v')
        allocate( fluid_vel   :: new_object(i)%icp )
      case ('thermal_v')
        allocate( thermal_vel :: new_object(i)%icp )
      case default
        stop 'Unknown initial conditions'
      end select
    end do
  end function init_icp_from_file

end module icp_factory

