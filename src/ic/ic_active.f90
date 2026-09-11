! =====================================================================
! NAME       : ic_active.f90
! DESCRIPTION: Initial conditions for active scalars in all
!              solver classes supporting ACTIVESC, and for the mass
!              density of the compressible solvers (an active scalar,
!              sector DENSITY): see active_sector below.
!
! Initial conditions available:
!   read_as   : Reads active scalars from input files numbered by stat
!   uniform_as: Uniform active scalar
!   puff_as   : Localized concentration
!   random_as : Random concentration
!
! DATE       : 04/08/26 (JBG)
! =====================================================================

module ic_active
  use icbase_mod
  use equationbase_mod

  implicit none

  ! ================= Initial conditions supported ====================
  type, extends(icBase) :: icRead_as
    contains
      procedure :: init_GState => init_readas
  end type icRead_as  
  type, extends(icBase) :: icUniform_as
    contains
      procedure :: init_GState => init_uniformas
  end type icUniform_as 
  type, extends(icBase) :: icPuff_as
    contains
      procedure :: init_GState => init_puffas
  end type icPuff_as 
  type, extends(icBase) :: icRandom_as
    contains
      procedure :: init_GState => init_randomas
  end type icRandom_as 

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Sector of the active scalars of the solver: the ACTIVESC
  !! sector of the solvers with active scalars, or the mass
  !! density (one scalar) of the compressible solvers
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine active_sector(solver,isc,nsc)
    class(EquationBase), intent(in)  :: solver
    integer            , intent(out) :: isc,nsc
    select type (solver)
    class is (ActiveScalarBase)
      isc = solver%ACTIVESC
      nsc = solver%numactivesc_
    class is (CompMagneticBase)
      isc = solver%DENSITY
      nsc = 1
    class is (CompVelocityBase)
      isc = solver%DENSITY
      nsc = 1
    class default
      error stop 'IC: This solver does not support active scalars'
    end select
  end subroutine active_sector


  ! ===================================================================
  ! Initial conditions
  ! ===================================================================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Read the active scalars from restart files
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_readas(this,solver,state)
    use gstate_mod
    use iovar
    use status
    use filefmt
    use fft
    use commtypes
    implicit none

    class   (icRead_as), intent(in)            :: this
    class(EquationBase), intent(in)            :: solver
    type   (GStateComp), intent(inout)         :: state(:)
    real(kind=GP), pointer, dimension(:,:,:)   :: R1
    integer                                    :: i
    integer                                    :: isc,nsc
    logical                                    :: bret

    if ((stat .eq. 0) .and. (solver%myrank_ .eq. 0)) then
      error stop 'Cannot read files if starting a new run with stat=0'
    endif
    call solver%workspace_%get_real_tmp(R1,bret)
    call active_sector(solver,isc,nsc)
      tind = int(stat)
      write(ext, fmtext) tind
      do i = isc, isc + nsc - 1
        call io_read(1,solver%idir_,trim(solver%sstate_(i)),ext,solver%planio_,R1)
        call fftp3d_real_to_complex(planrc,R1,state(i)%ccomp,MPI_COMM_WORLD)
      end do
    call solver%workspace_%free_real_tmp(R1)
  end subroutine init_readas


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Uniform active scalar
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   c0 : vector with the amplitudes of the active scalars
  subroutine init_uniformas(this,solver,state)
    use gstate_mod
    use grid
    use commtypes
    use mpivars
    implicit none

    class(icUniform_as), intent(in)           :: this
    class(EquationBase), intent(in)           :: solver
    type   (GStateComp), intent(inout)        :: state(:)
    real      (kind=GP), allocatable          :: c0(:)
    integer                                   :: i,j,k,n
    integer                                    :: isc,nsc

    namelist /uniform_as/ c0
    call active_sector(solver,isc,nsc)
      allocate(c0(nsc))
      if (myrank .eq. 0) then
        open(1,file=solver%infile_,status='unknown',form='formatted')
        read(1,nml=uniform_as)
        close(1)
      endif
      call mpi_bcast(c0,nsc,GC_REAL,0,MPI_COMM_WORLD,ierr)
      do n = isc, isc + nsc - 1
!$omp parallel do collapse(2) private (k)
        do i = ista,iend
          do j = 1,ny
            do concurrent (k=1:nz)
              state(n)%ccomp(k,j,i) = 0.0_GP
            end do
          end do
        end do
        if (myrank .eq. 0) then
          state(n)%ccomp(1,1,1) = c0(n-isc+1) * &
                real(nx,kind=GP) * real(ny,kind=GP) * real(nz,kind=GP)
        endif
      end do
  end subroutine init_uniformas


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Active scalar in a Gaussian ball ('puff') centered at
  !! (x0, y0, z0) with a FWHM of r0
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   c0        : vector with the amplitudes
  !   x0, y0, z0: vectors with the centers
  !   r0        : vector with the radii
  subroutine init_puffas(this,solver,state)
    use gstate_mod
    use grid
    use mpivars
    use commtypes
    use fft
    use pseudospec_scalar
    implicit none

    class   (icPuff_as), intent(in)            :: this
    class(EquationBase), intent(in)            :: solver
    type   (GStateComp), intent(inout)         :: state(:)
    real(kind=GP), pointer                     :: R1(:,:,:)
    real(kind=GP), allocatable, dimension(:)   :: c0,x0,y0,z0,r0
    double precision                           :: tmp
    integer                                    :: i,j,k,n
    integer                                    :: isc,nsc
    logical                                    :: bret

    namelist /puff_as/ c0,x0,y0,z0,r0
    call solver%workspace_%get_real_tmp(R1,bret)
    call active_sector(solver,isc,nsc)
      allocate(c0(nsc))
      allocate(x0(nsc))
      allocate(y0(nsc))
      allocate(z0(nsc))
      allocate(r0(nsc))
      if (myrank .eq. 0) then
        open(1,file=solver%infile_,status='unknown',form='formatted')
        read(1,nml=puff_as)
        close(1)
      endif
      call mpi_bcast(c0,nsc,GC_REAL,0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(x0,nsc,GC_REAL,0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(y0,nsc,GC_REAL,0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(z0,nsc,GC_REAL,0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(r0,nsc,GC_REAL,0,MPI_COMM_WORLD,ierr)
      do n = isc, isc + nsc - 1
!$omp parallel do collapse(2) private (i)
        do k = ksta,kend
          do j = 1,ny
            do concurrent (i=1:nx)
              tmp = (real(i-1,kind=GP)/real(nx-1,kind=GP)-x0(n-isc+1))**2 &
                  + (real(j-1,kind=GP)/real(ny-1,kind=GP)-y0(n-isc+1))**2 &
                  + (real(k-1,kind=GP)/real(nz-1,kind=GP)-z0(n-isc+1))**2
              R1(i,j,k) = exp(-tmp**2/r0(n-isc+1)**2)
            end do
          end do
        end do
        call fftp3d_real_to_complex(planrc,R1,state(n)%ccomp,MPI_COMM_WORLD)
        call variance(state(n)%ccomp,tmp,1)
        call mpi_bcast(tmp,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!$omp parallel do collapse(2) private (k)
        do i = ista,iend
          do j = 1,ny
            do concurrent (k=1:nz)
              state(n)%ccomp(k,j,i) = state(n)%ccomp(k,j,i) * &
                                      c0(n-isc+1)/sqrt(tmp)
            end do
          end do
        end do
      end do
    call solver%workspace_%free_real_tmp(R1)
  end subroutine init_puffas


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Random active scalar concentration
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   c0 : vector with the amplitudes
  !   kdn: vector with minimum wave number
  !   kup: vector with maximum wave number
  subroutine init_randomas(this,solver,state)
    use gstate_mod
    use random
    use grid
    use mpivars
    use commtypes
    use status
    use fft
    use var
    use kes
    use pseudospec_scalar
    implicit none

    class (icRandom_as), intent(in)            :: this
    class(EquationBase), intent(in)            :: solver
    type   (GStateComp), intent(inout)         :: state(:)
    real(kind=GP), allocatable, dimension(:)   :: c0,kdn,kup
    real(kind=GP)                              :: skup,skdn
    real(kind=GP)                              :: dump,phase
    double precision                           :: tmp
    integer                                    :: i,j,k,n
    integer                                    :: isc,nsc

    namelist /random_as/ c0,kup,kdn
    call active_sector(solver,isc,nsc)
      allocate(c0 (nsc))
      allocate(kdn(nsc))
      allocate(kup(nsc))
      if (myrank .eq. 0) then
        open(1,file=solver%infile_,status='unknown',form='formatted')
        read(1,nml=random_as)
        close(1)
      endif
      call mpi_bcast(c0 ,nsc,GC_REAL,0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(kdn,nsc,GC_REAL,0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(kup,nsc,GC_REAL,0,MPI_COMM_WORLD,ierr)
      do n = isc, isc + nsc - 1
        skdn = kdn(n-isc+1)
        skup = kup(n-isc+1)
        if (ista .eq. 1) then
          state(n)%ccomp(1,1,1) = 0.0_GP
          do j = 2,ny/2+1
            if ((kk2(1,j,1) .le. skup**2) .and. (kk2(1,j,1) .ge. skdn**2)) then
              dump = 1.0_GP/sqrt(kk2(1,j,1))
              phase = 2*pi*randu(seed)
              state(n)%ccomp(1,j,1) = (cos(phase)+im*sin(phase))*dump
              state(n)%ccomp(1,ny-j+2,1) = conjg(state(n)%ccomp(1,j,1))
            else
              state(n)%ccomp(1,j,1) = 0.0_GP
              state(n)%ccomp(1,ny-j+2,1) = 0.0_GP
            endif
          end do
          do k = 2,nz/2+1
            if ((kk2(k,1,1) .le. skup**2) .and. (kk2(k,1,1) .ge. skdn**2)) then
              dump = 1.0_GP/sqrt(kk2(k,1,1))
              phase = 2*pi*randu(seed)
              state(n)%ccomp(k,1,1) = (cos(phase)+im*sin(phase))*dump
              state(n)%ccomp(nz-k+2,1,1) = conjg(state(n)%ccomp(k,1,1))
            else
              state(n)%ccomp(k,1,1) = 0.0_GP
              state(n)%ccomp(nz-k+2,1,1) = 0.0_GP
            endif
          end do
          do j = 2,ny
            do k = 2,nz/2+1
              if ((kk2(k,j,1) .le. skup**2) .and. (kk2(k,j,1) .ge. skdn**2)) then
                dump = 1.0_GP/sqrt(kk2(k,j,1))
                phase = 2*pi*randu(seed)
                state(n)%ccomp(k,j,1) = (cos(phase)+im*sin(phase))*dump
                state(n)%ccomp(nz-k+2,ny-j+2,1) = conjg(state(n)%ccomp(k,j,1))
              else
                state(n)%ccomp(k,j,1) = 0.0_GP
                state(n)%ccomp(nz-k+2,ny-j+2,1) = 0.0_GP
              endif
            end do
          end do
          do i = 2,iend
            do j = 1,ny
              do k = 1,nz
                if ((kk2(k,j,i) .le. skup**2) .and. (kk2(k,j,i) .ge. skdn**2)) then
                  dump = 1.0_GP/sqrt(kk2(k,j,i))
                  phase = 2*pi*randu(seed)
                  state(n)%ccomp(k,j,i) = 2.0_GP*(cos(phase)+im*sin(phase))*dump
                else
                  state(n)%ccomp(k,j,i) = 0.0_GP
                endif
              end do
            end do
          end do
        else
          do i = ista,iend
            do j = 1,ny
              do k = 1,nz
                if ((kk2(k,j,i) .le. skup**2) .and. (kk2(k,j,i) .ge. skdn**2)) then
                  dump = 1.0_GP/sqrt(kk2(k,j,i))
                  phase = 2*pi*randu(seed)
                  state(n)%ccomp(k,j,i) = 2.0_GP*(cos(phase)+im*sin(phase))*dump
                else
                  state(n)%ccomp(k,j,i) = 0.0_GP
                endif
              end do
            end do
          end do
        endif
        call variance(state(n)%ccomp,tmp,1)
        call mpi_bcast(tmp,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!$omp parallel do collapse(2) private (k)
        do i = ista,iend
          do j = 1,ny
            do concurrent (k=1:nz)
              state(n)%ccomp(k,j,i) = state(n)%ccomp(k,j,i) * &
                                      c0(n-isc+1)/sqrt(tmp)
            end do
          end do
        end do
      end do
  end subroutine init_randomas

end module ic_active
