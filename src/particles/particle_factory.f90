! ===================================================================
! Factory for all particles.
!
! To add new particles, "USE" the corresponding module, and add a
! new case clause in the init_particles_from_file function. The
! clause should allocate the solver class and declare the number of
! tmp storage components the particle solver needs.
! DATE : 19/02/26 (PDM)
! ===================================================================

module particle_factory
  USE particlebase_mod
  USE lagpart_mod
  USE inerpart_mod
! USE userdefinedparticle_mod

  IMPLICIT NONE

  ! ================= Global parameters ===============================
  integer, public :: NUMTMPPART = 0   ! Number of tmp arrays for particles

CONTAINS

  ! ================= Factory function ==============================
  function init_particles_from_file(infile) result(new_object)
    USE commtypes
    class(ParticleBase), allocatable  :: new_object
    character(len=*)   , intent(in)   :: infile

    ! Temporary data to read from namelist:
    integer            :: nprocs,myrank,ierr,ios
    character(len=256) :: iomsg
    character(len=64)  :: psolver
    ! Required namelist:
    namelist/ particles / psolver

    call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
    if ( myrank .eq. 0 ) then
      open(1,file=infile,status='unknown',form="formatted")
      read(1,NML=particles,iostat=ios,iomsg=iomsg)
      close(1)
    endif
    call MPI_BCAST(ios,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    if (ios .eq. 0) then
      call MPI_BCAST(psolver,64,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      ! Clauses for each solver class
      select case (trim(adjustl(psolver)))
        case ('lagpart')
          allocate(Gpart  :: new_object)
          NUMTMPPART = 3
        case ('inerpart')
          allocate(Ipart  :: new_object)
          NUMTMPPART = 3
!       case ('maxey')
!         allocate(MRpart :: new_object)
!         NUMTMPPART = 3
!       case ('testpart')
!         allocate(Tpart  :: new_object)
!         NUMTMPPART = 5
!       case ('UserDefined')
!         allocate(UserDefinedpart :: new_object)
!         NUMTMPPART = 2
        case default
          stop 'Particle factory :: init_particles_from_file: Unknown solver name'
      end select
    else ! Not running with particles
      NUMTMPPART = 0
    endif
  end function init_particles_from_file

end module particle_factory
