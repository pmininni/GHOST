! ===================================================================
! Factory for all forcing functions
!
! This factory applies chain operation to create multiple forcings
! and forcing update methods for all fields in a pde solver.
!
! DATE : 01/18/26 (PDM)
! ===================================================================

module force_factory
  USE forcebase_mod
  USE force_velocity
  
  IMPLICIT NONE
  
CONTAINS
  
  ! ================= Factory function ==============================
  function init_forcing_from_file(infile) result(new_object)
    USE gstate_mod
    USE commtypes
    USE gutils
    type(forceChain), allocatable :: new_object(:)
    character(len=*),  intent(in) :: infile

    ! Temporary data:
    real(kind=GP)      :: cort
    integer            :: i,ib,ie,ind,fnmb,unmb
    character(len=512) :: forces
    character(len= 50) :: forcelist(10)
    character( len=25) :: forcename(2)
    ! Required namelist:
    namelist/ forcing / forces,cort

    if ( myrank .eq. 0 ) then
      open(1,file=infile,status='unknown',form="formatted")
      read(1,NML=forcing)
      close(1)
    endif
    call MPI_BCAST(forces,512,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(cort  ,  1,GC_REAL      ,0,MPI_COMM_WORLD,ierr)

    ! Loops through all requested forces and allocates classes
    call parsestr(forces,';',forcelist ,50,10,fnmb)
    if ( myrank .eq. 0 ) then
       print *,'Found',fnmb,'force methods'
    endif
    allocate( new_object(fnmb) )
    do i = 1,fnmb
      call parsestr(forcelist(i),',',forcename,25,2,unmb)
      if ( myrank .eq. 0 ) print *,'Forcing',i,':', &
           trim(adjustl(forcename(1))),',',trim(adjustl(forcename(2)))
      select case (trim(adjustl(forcename(1))))
      ! Velocity forcing functions  ------------
      case ('null_fv')
        allocate( null_fv        :: new_object(i)%force )
      case ('tg_fv')
        allocate( tg_fv          :: new_object(i)%force )
      case ('abc_fv')
        allocate( abc_fv         :: new_object(i)%force )
      case ('random_fv')
        allocate( random_fv      :: new_object(i)%force )
!!$      ! Magnetic forcing functions -------------
!!$      case ('read_b')
!!$        allocate( read_b   :: new_object(i)%ic )
!!$      case ('null_b')
!!$        allocate( null_b   :: new_object(i)%ic )
!!$      case ('random_b')
!!$        allocate( random_b :: new_object(i)%ic )
      case default
        stop 'Unknown forcing function'
      end select
      select case (trim(adjustl(forcename(2))))
      ! Velocity forcing update methods --------
      case ('constant_fv')
        if ( allocated(new_object(i)%update) ) deallocate(new_object(i)%update)
      case ('shift_fv')
        allocate( shiftupdt_fv   :: new_object(i)%update )
      case ('shuffle_fv')
        allocate( shuffleupdt_fv :: new_object(i)%update )
        new_object(i)%update%binit_ = .FALSE.
!!$      ! Magnetic field ICs ---------------------    
!!$      case ('read_b')
!!$        allocate( read_b   :: new_object(i)%ic )
!!$      case ('null_b')
!!$        allocate( null_b   :: new_object(i)%ic )
!!$      case ('random_b')
!!$        allocate( random_b :: new_object(i)%ic )
      case default
        stop 'Unknown or undefined forcing update method'
      end select
    end do
  end function init_forcing_from_file

end module force_factory
