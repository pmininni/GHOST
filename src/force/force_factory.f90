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
  USE force_magnetic
  USE force_active
  USE force_passive
  
  IMPLICIT NONE
  
CONTAINS
  
  ! ================= Factory function ==============================
  function init_forcing_from_file(infile,workspace) result(new_object)
    USE class_GWorkspace3D
    USE gstate_mod
    USE commtypes
    USE gutils
    type(forceChain),   allocatable :: new_object(:)
    type(GWorkspace), intent(inout) :: workspace 
    character(len=*), intent   (in) :: infile
    ! Temporary data:
    real(kind=GP)      :: cort
    integer            :: i,ib,ie,ind,fnmb,unmb,poolsz = 0
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

    ! Loops through all requested forces and allocates classes.
    call parsestr(forces,';',forcelist ,50,10,fnmb)
    if ( myrank .eq. 0 ) then
       print *,'Found',fnmb,'force methods'
    endif
    allocate( new_object(fnmb) )
    do i = 1,fnmb
      call parsestr(forcelist(i),',',forcename,25,2,unmb)
      if ( myrank .eq. 0 ) print *,'Forcing',i,':', &
           trim(adjustl(forcename(1))),',',trim(adjustl(forcename(2)))
      ! Forcing functions  ------------------------------------------
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
      ! Electromotive forcing functions --------
      case ('null_fb')
        allocate( null_fb        :: new_object(i)%force )
      case ('random_fb')
        allocate( random_fb      :: new_object(i)%force )
      ! Active scalar forcing functions --------
      case ('null_fas')
        allocate( null_fas       :: new_object(i)%force )
      case ('puff_fas')
        allocate( puff_fas       :: new_object(i)%force )
      case ('random_fas')
        allocate( random_fas     :: new_object(i)%force )
      ! Passive scalar forcing functions ------- 
      case ('null_fs')
        allocate( null_fs        :: new_object(i)%force )
      case ('puff_fs')
        allocate( puff_fs        :: new_object(i)%force )
      case ('random_fs')
        allocate( random_fs      :: new_object(i)%force )
      case default
        stop 'Unknown forcing function'
      end select
      ! Update methods ----------------------------------------------
      select case (trim(adjustl(forcename(2))))
      ! Velocity forcing update methods --------
      case ('constant_fv')
        if ( allocated(new_object(i)%update) ) deallocate(new_object(i)%update)
      case ('shift_fv')
        allocate( shiftupdt_fv   :: new_object(i)%update )
      case ('shuffle_fv')
        allocate( shuffleupdt_fv :: new_object(i)%update )
        new_object(i)%update%binit_ = .FALSE.
        poolsz = poolsz + 6
      ! Electromotive force update methods -----
      case ('constant_fb')
        if ( allocated(new_object(i)%update) ) deallocate(new_object(i)%update)
      case ('shift_fb')
        allocate( shiftupdt_fb   :: new_object(i)%update )
      case ('shuffle_fb')
        allocate( shuffleupdt_fb :: new_object(i)%update )
        new_object(i)%update%binit_ = .FALSE.
        poolsz = poolsz + 6
      ! Active scalar forcing functions --------
      case ('constant_fas')
        if ( allocated(new_object(i)%update) ) deallocate(new_object(i)%update)
      case ('shift_fas')
        allocate( shiftupdt_fas  :: new_object(i)%update )
      ! Passive scalar forcing functions -------
      case ('constant_fs')
        if ( allocated(new_object(i)%update) ) deallocate(new_object(i)%update)
      case ('shift_fs')
        allocate( shiftupdt_fs  :: new_object(i)%update )
      case default
        stop 'Unknown or undefined forcing update method'
      end select
    end do
    ! If needed, we increase the pool size here as it is still safe for 
    ! for permanent storage
    if ( poolsz .gt. 0 ) call workspace%add_complex_entries(poolsz)
  end function init_forcing_from_file

end module force_factory
