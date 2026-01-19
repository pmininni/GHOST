! ===================================================================
! Factory for all initial conditions (ICs)
!
! To add a new initial condition...
! NOTE: Maybe GState should be part of the IC class. After all, IC
! does an "init" of the states. This would allow us to mark the
! Gstates depending on whether they have velocity, magnetic fields,
! quantum (this is now done with subclasses in the solvers). It would
! also allow for double checks.
! Following the same idea, the forcing vector (which should be more
! flexible and allow for complex and real elements) could be part of
! a "forcing" class.
!
! DATE : 01/16/26 (PDM)
! ===================================================================

module ic_factory
  USE gstate_mod
  USE icbase_mod
  USE ic_velocity
  
  IMPLICIT NONE
  
CONTAINS
  
  ! ================= Factory function ==============================
  function init_ic_from_file(infile) result(new_object)
    USE icbase_mod
    USE commtypes
    class   (icBase), allocatable :: new_object
    character(len=*),  intent(in) :: infile

    ! Temporary data to read from namelist:
    integer            :: i,ib,ie,ind,ic_nmb
    character(len=256) :: iclist
    character( len=25) :: icname(10)
    ! Required namelist:
    namelist/ initial_conditions / iclist

    if ( myrank .eq. 0 ) then
      open(1,file=infile,status='unknown',form="formatted")
      read(1,NML=initial_conditions)
      close(1)
    endif
    call MPI_BCAST(iclist,256,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)

    ! Loops through all requested ICs and allocates classes
    ib = 1
    ie = len(iclist)
    i  = 1
    do while ( len(trim(iclist(ib:ie))) .GT. 0 )  
      ind = index(iclist(ib:ie),";")
      if ( ind .eq. 0 ) then
        icname(i) = trim(adjustl(iclist(ib:ie)))
        ib = ie + 1
        i  = i+1
      else
        icname(i) = trim(adjustl(iclist(ib:(ib+ind-2))))
        ib = ib + ind
      endif
    end do
    ic_nmb  = i-1
!   do i = 1,ic_nmb
    i = 1
    select case (trim(adjustl(icname(i))))
      case ('null_v')
        allocate( null_v :: new_object )
      case ('tg_v')
        allocate( tg_v   :: new_object )
      case default
        stop 'Unknown initial conditions'
    end select
  end function init_ic_from_file

end module ic_factory
