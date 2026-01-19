! ===================================================================
! Factory for all forcing functions
!
! To add a new forcing scheme...
! NOTE: Forcing needs two methods. One to generate the initial
! forcing, one to update the forcing. Some could be combined,
! some not. But the actual implementation does not allow this.
! We may have to use two different classes for that.
!
! DATE : 01/18/26 (PDM)
! ===================================================================

module force_factory
  USE gstate_mod
  USE forcebase_mod
  USE force_velocity
  
  IMPLICIT NONE
  
CONTAINS
  
  ! ================= Factory function ==============================
  function init_forcing_from_file(infile) result(new_object)
    USE forcebase_mod
    USE commtypes
    class(forceBase), allocatable :: new_object
    character(len=*),  intent(in) :: infile

    ! Temporary data to read from namelist:
    integer            :: i,ib,ie,ind,force_nmb
    character(len=256) :: forces,updates
    character( len=25) :: forcename(10),updatename(10)
    ! Required namelist:
    namelist/ forcing / forces,updates

    if ( myrank .eq. 0 ) then
      open(1,file=infile,status='unknown',form="formatted")
      read(1,NML=forcing)
      close(1)
    endif
    call MPI_BCAST(forces ,256,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(updates,256,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)

    ! Loops through all requested forces and allocates classes
    ib = 1
    ie = len(forces)
    i  = 1
    do while ( len(trim(forces(ib:ie))) .GT. 0 )  
      ind = index(forces(ib:ie),";")
      if ( ind .eq. 0 ) then
        forcename(i) = trim(adjustl(forces(ib:ie)))
        ib = ie + 1
        i  = i+1
      else
        forcename(i) = trim(adjustl(forces(ib:(ib+ind-2))))
        ib = ib + ind
      endif
    end do
    force_nmb  = i-1
!   do i = 1,force_nmb
    i = 1
    select case (trim(adjustl(forcename(i))))
      case ('null_fv')
        allocate( null_fv :: new_object )
      case ('tg_fv')
        allocate( tg_fv   :: new_object )
      case default
        stop 'Unknown forcing'
    end select
  end function init_forcing_from_file

end module force_factory
