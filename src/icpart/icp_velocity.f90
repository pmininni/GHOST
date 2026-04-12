! =====================================================================
! NAME       : icp_velocity.f90
! DESCRIPTION: Initial conditions for the velocity of all
!              VelocParticleBase solver classes. These initial
!              conditions can be used to initialize velocities of
!              particles in other solver classes using composition of
!              initial condition operators.
!
! Initial conditions avaliable:
!   read_v   : Reads velocities from an input file numbered by stat
!   null_v   : Null initial velocity
!   fluid_v  : Each particle gets the velocity of its fluid element
!   thermal_v: Random initial velocities with Maxwell-Boltzmann dist
!
! DATE       : 04/11/26 (PDM)
! =====================================================================

module icp_velocity
  use icpbase_mod
  use equationbase_mod
  use particlebase_mod

  IMPLICIT NONE

  ! ================= Initial conditions supported ====================
  type, extends(icpBase) :: read_vel
    contains
      procedure :: init_GPState => init_readvel
  end type read_vel
  type, extends(icpBase) :: null_vel
    contains
      procedure :: init_GPState => init_nullvel
  end type null_vel 
  type, extends(icpBase) :: fluid_vel
    contains
      procedure :: init_GPState => init_fluidvel
  end type fluid_vel
  type, extends(icpBase) :: thermal_vel
    contains
      procedure :: init_GPState => init_thermvel
  end type thermal_vel

CONTAINS

  ! ===================================================================
  ! Initial conditions
  ! ===================================================================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Read particle velocity from restart files.
  !! This method assumes that the chain has read first the
  !! particles positions.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_readvel(this,pde,fluidstate,psolver,pstate)
    use gpstate_mod
    use gstate_mod
    use status
    use pstatus
    use filefmt
    use commtypes
    implicit none
    class    (read_vel),         intent   (in) :: this
    class(EquationBase),         intent   (in) :: pde
    type   (GStateComp), target, intent   (in) :: fluidstate(:) 
    class(ParticleBase),         intent(inout) :: psolver
    type  (GPStateComp),         intent(inout) :: pstate(:)

    if ((stat .eq. 0).and.(psolver%myrank_ .eq. 0)) then
      stop 'Cannot read restart files if starting a new run with stat=0'
    endif
    select type (psolver)
    class is (VelocParticleBase)
      pind = int((stat-1)*lgmult+1)
      write(lgext, lgfmtext) pind
      call io_readvec(psolver,pstate,1,psolver%idir_,trim(psolver%sstate_vel_), &
                      lgext,psolver%VELOCITY)
    class default
      stop "Init_readvel: These ICs do not support particles without velocity"
    end select
  end subroutine init_readvel


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Initialize velocity to zero
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_nullvel(this,pde,fluidstate,psolver,pstate)
    use gpstate_mod
    use gstate_mod
    implicit none
    class    (null_vel),         intent   (in) :: this
    class(EquationBase),         intent   (in) :: pde
    type   (GStateComp), target, intent   (in) :: fluidstate(:) 
    class(ParticleBase),         intent(inout) :: psolver
    type  (GPStateComp),         intent(inout) :: pstate(:)
    integer                                    :: m

    select type (psolver)
    class is (VelocParticleBase)
      do m = 1,psolver%nc_
        pstate(psolver%VELOCITY+m-1)%rcomp = 0.0_GP
      end do
    class default
      stop "Init_fluidvel: These ICs do not support particles without velocity"
    end select
  end subroutine init_nullvel


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Initialize particle velocities using the fluid velocity
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_fluidvel(this,pde,fluidstate,psolver,pstate)
    use gpstate_mod
    use random
    use fft
!$  use threads
    implicit none
    class   (fluid_vel),             intent   (in) :: this
    class(EquationBase),             intent   (in) :: pde
    type   (GStateComp), target,     intent   (in) :: fluidstate(:) 
    class(ParticleBase),             intent(inout) :: psolver
    type  (GPStateComp),             intent(inout) :: pstate(:)
    complex   (kind=GP), pointer, dimension(:,:,:) :: velc
    real      (kind=GP), pointer, dimension(:,:,:) :: velr,tmp1,tmp2
    real      (kind=GP)  :: rmp
    integer              :: i,j,k,m
    logical              :: doupdate(psolver%nc_)
    logical              :: bret

    call psolver%workspace_%get_complex_tmp(velc,bret)
    call psolver%workspace_%get_real_tmp   (velr,bret)
    call psolver%workspace_%get_real_tmp   (tmp1,bret)
    call psolver%workspace_%get_real_tmp   (tmp2,bret)
    call AssignLagPos(psolver, pstate) ! We assign px_,py_,pz_ to the pstate

    select type (pde)
    class is (VelocityBase)
      select type (psolver)
      class is (VelocParticleBase)
      if (psolver%nc_ .ne. pde%nc_) then
        stop "Init_fluidvel: # of components of the particles and pdes must be equal"
      endif
      doupdate    = .false.
      doupdate(1) = .true.
      rmp = 1.0_GP/(real(psolver%nd_(1),kind=GP)*real(psolver%nd_(2),kind=GP)* &
                    real(psolver%nd_(3),kind=GP))
      do m = 1,pde%nc_
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
        do i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
          do j = 1,ny
            do k = 1,nz
              velc(k,j,i) = fluidstate(pde%VELOCITY+m-1)%ccomp(k,j,i)*rmp
            end do
          end do
        end do
        call fftp3d_complex_to_real(plancr,velc,velr,MPI_COMM_WORLD)
        call EulerToLag(psolver,pstate(psolver%VELOCITY+m-1)%rcomp,            &
                        psolver%nparts_,velr,doupdate(m),tmp1,tmp2)
      end do
      class default
        stop "Init_fluidvel: These ICs do not support particles without velocity"
      end select
    class default
      stop "Init_fluidvel: These ICs do not support pdes without a velocity field"
    end select
 
    call psolver%workspace_%free_complex_tmp(velc)
    call psolver%workspace_%free_real_tmp   (velr)
    call psolver%workspace_%free_real_tmp   (tmp1)
    call psolver%workspace_%free_real_tmp   (tmp2)
  end subroutine init_fluidvel


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Initializes particle velocities with a Gaussian 
  !! distribution function with thermal speed vtherm.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_thermvel(this,pde,fluidstate,psolver,pstate)
    use gpstate_mod
    use gstate_mod
    use random
    use status
    use pstatus
    use filefmt
    use commtypes
    implicit none
    class (thermal_vel),         intent   (in) :: this
    class(EquationBase),         intent   (in) :: pde
    type   (GStateComp), target, intent   (in) :: fluidstate(:) 
    class(ParticleBase),         intent(inout) :: psolver
    type  (GPStateComp),         intent(inout) :: pstate(:)
    real      (kind=GP)                        :: vtherm,gauss
    real      (kind=GP)                        :: low,twopi,u1,u2
    integer                                    :: j
    
    namelist/ thermal_v / vtherm
    select type (psolver)
    class is (VelocParticleBase)
      ! Read parameters from a namelist in the input file
      if ( myrank .eq. 0 ) then
        open(1,file=psolver%infile_,status='unknown',form="formatted")
        read(1,NML=thermal_v)
        close(1)
    endif
    call mpi_bcast(vtherm,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
    ! Generate a Gaussian distribution with thermal speed vtherm
    twopi = 8.0_GP*atan(1.0_GP)
    low   = 1.0e-20_GP
!$omp parallel do
    DO j = 1, psolver%nparts_
      CALL prandom_number(u1)
      CALL prandom_number(u2)
      u1 =  MAX( u1, low)
      u2 =  twopi * u2
      gauss =  sqrt( -2.0_GP * log( u1)) * CMPLX( cos(u2), sin(u2))
      pstate(psolver%VELOCITY  )%rcomp = vtherm*gauss
      CALL prandom_number(u1)
      CALL prandom_number(u2)
      u1 =  MAX( u1, low)
      u2 =  twopi * u2
      gauss =  sqrt( -2.0_GP * log( u1)) * CMPLX( cos(u2), sin(u2))
      pstate(psolver%VELOCITY+1)%rcomp = vtherm*gauss
      CALL prandom_number(u1)
      CALL prandom_number(u2)
      u1 =  MAX( u1, low)
      u2 =  twopi * u2
      gauss =  sqrt( -2.0_GP * log( u1)) * CMPLX( cos(u2), sin(u2))
      pstate(psolver%VELOCITY+2)%rcomp = vtherm*gauss
    END DO     
    class default
      stop "Init_thermvel: These ICs do not support particles without velocity"
    end select
  end subroutine init_thermvel
    
end module icp_velocity
