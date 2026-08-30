! ===================================================================
! Factory for all initial conditions (ICs)
!
! This factory applies chain operation to create multiple initial
! conditions for all fields in a pde solver.
!
! DATE : 04/08/26 (JBG)
! ===================================================================

module ic_factory
  USE icbase_mod
  USE ic_velocity
  USE ic_magnetic
  USE ic_active
  USE ic_passive
  
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
        allocate( icRead_v      :: new_object(i)%ic )
      case ('null_v')
        allocate( icNull_v      :: new_object(i)%ic )
      case ('tg_v')
        allocate( icTg_v        :: new_object(i)%ic )
      case ('abc_v')
        allocate( icAbc_v       :: new_object(i)%ic )
      case ('random_v')
        allocate( icRandom_v    :: new_object(i)%ic )
      ! Magnetic field ICs ---------------------    
      case ('read_b')
        allocate( icRead_b      :: new_object(i)%ic )
      case ('null_b')
        allocate( icNull_b      :: new_object(i)%ic )
      case ('random_b')
        allocate( icRandom_b    :: new_object(i)%ic )
      ! Active scalar ICs ---------------------    
      case ('read_as')
        allocate( icRead_as     :: new_object(i)%ic )
      case ('uniform_as')
        allocate( icUniform_as :: new_object(i)%ic )
      case ('puff_as')
        allocate( icPuff_as     :: new_object(i)%ic )
      case ('random_as')
        allocate( icRandom_as   :: new_object(i)%ic )
      ! Passive scalar ICs ---------------------    
      case ('read_s')
        allocate( icRead_s      :: new_object(i)%ic )
      case ('uniform_s')
        allocate( icUniform_s  :: new_object(i)%ic )
      case ('puff_s')
        allocate( icPuff_s      :: new_object(i)%ic )
      case ('random_s')
        allocate( icRandom_s    :: new_object(i)%ic )
      case default
        stop 'Unknown initial conditions'
      end select
    end do
  end function init_ic_from_file

end module ic_factory
