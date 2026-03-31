! ===================================================================
! Factory for all initial conditions (ICs)
!
! This factory applies chain operation to create multiple initial
! conditions for all fields in a pde solver.
!
! DATE : 01/16/26 (PDM)
! ===================================================================

module ic_factory
  USE icbase_mod
  USE ic_velocity
  USE ic_magnetic
  USE ic_passive
  USE ic_active
  
  IMPLICIT NONE
  
CONTAINS
  
  ! ================= Factory function ==============================
  function init_ic_from_file(infile) result(new_object)
    USE gstate_mod
    USE commtypes
    USE gutils
    type   (icChain), allocatable :: new_object(:)
    character(len=*),  intent(in) :: infile

    ! Temporary data to read from namelist:
    integer            :: i,ib,ie,ind,icnmb
    character(len=256) :: iclist
    character(len= 25) :: icname(10)
    ! Required namelist:
    namelist/ initial_conditions / iclist

    if ( myrank .eq. 0 ) then
      open(1,file=infile,status='unknown',form="formatted")
      read(1,NML=initial_conditions)
      close(1)
    endif
    call MPI_BCAST(iclist,256,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)

    ! Loops through all requested ICs and allocates classes
    call parsestr(iclist,';',icname,25,10,icnmb)
    if ( myrank .eq. 0 ) print *,'Found',icnmb,'ICs'
    allocate( new_object(icnmb) )
    do i = 1,icnmb
      if ( myrank .eq. 0 ) print *,'IC',i,':',trim(adjustl(icname(i)))
      select case (trim(adjustl(icname(i))))
      ! Velocity field ICs ---------------------
      case ('read_v')
        allocate( read_v   :: new_object(i)%ic )
      case ('null_v')
        allocate( null_v   :: new_object(i)%ic )
      case ('tg_v')
        allocate( tg_v     :: new_object(i)%ic )
      case ('abc_v')
        allocate( abc_v    :: new_object(i)%ic )
      case ('random_v')
        allocate( random_v :: new_object(i)%ic )
      ! Magnetic field ICs ---------------------    
      case ('read_b')
        allocate( read_b   :: new_object(i)%ic )
      case ('null_b')
        allocate( null_b   :: new_object(i)%ic )
      case ('random_b')
        allocate( random_b :: new_object(i)%ic )
      ! Passive scalar ICs ---------------------    
      case ('read_s')
        allocate( read_s     :: new_object(i)%ic )
      case ('constant_s')
        allocate( constant_s :: new_object(i)%ic )
      case ('puff_s')
        allocate( puff_s     :: new_object(i)%ic )
      case ('random_s')
        allocate( random_s   :: new_object(i)%ic )
      case default
        stop 'Unknown initial conditions'
      end select
    end do
  end function init_ic_from_file

end module ic_factory
