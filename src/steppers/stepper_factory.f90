! ===================================================================
! Factory for all allowed time stepper objects
!
! DATE : 02/8/26 (DLR)
! ===================================================================

module stepper_factory
  use gstepperbase_mod
  use canuto_stepper_mod
! use gexrk_stepper_mod
  implicit none

  ! ================= Global parameters ===============================
  
contains
  
  ! ================= Factory function ==============================
  function build_stepper_from_file(infile) result(new_object)
    use gstepperbase_mod
    use commtypes
    class(GStepperBase), allocatable  :: new_object
    character   (len=*), intent(in)   :: infile

    ! Temporary data to read from namelists:
    integer                                :: itype,norder,nstage
    character(len=128)                     :: sname ! stepper name
    type(GStepperTraits)                   :: straits

    ! Required namelists:
    namelist/ stepper    / sname, itype, norder, nstage

    itype          = 1 ! Butcher type if using GEXRK object
    norder         = 2
    nstage         = 2
    sname          = 'TRADITIONAL'

    call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
    if ( myrank .eq. 0 ) then
      open(1,file=infile,status='unknown',form="formatted")
      read(1,NML=stepper)
      close(1)
    endif
    call mpi_bcast(sname    ,128 ,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(itype    ,1   ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(norder   ,1   ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(nstage   ,1   ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)

!!$    ! Build stepper traits structure:
!!$    straits%nstate = this%state_size()
!!$    straits%norder = norder
!!$    straits%nstage = nstage
!!$    straits%sname  = sname

   ! Allocate child object:
   select case (trim(adjustl(sname)))
!     case ('gexrk')
!       allocate( GExRKStepper  :: new_object )
     case ('traditional', 'TRADITIONAL')
       allocate( CanutoStepper :: new_object )
     case default
       stop 'stepper_factory::build_stepper_from_file: Invalid stepper type'
   end select

   ! Call constructor:
   !new_object%GStepper_ctor(straits,workspace)

  end function build_stepper_from_file

end module stepper_factory
