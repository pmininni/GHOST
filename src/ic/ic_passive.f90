! =====================================================================
! NAME       : ic_passive.f90
! DESCRIPTION: Initial conditions for passive scalars in all
!              VelocityBase solver classes.
!
! Initial conditions avaliable:
!   read_s    : Reads passive scalars from input files numbered by stat
!   constant_s: Uniform passive scalar
!   puff_s    : Localized concentration
!   random_s  : Random concentration
!
! DATE       : 01/26/26 (PDM)
! =====================================================================

module ic_passive
  USE icbase_mod

  IMPLICIT NONE

  ! ================= Initial conditions supported ====================
  type, extends(icBase) :: icRead_s
    contains
      procedure :: init_GState => init_reads
  end type icRead_s 
  type, extends(icBase) :: icConstant_s
    contains
      procedure :: init_GState => init_constants
  end type icConstant_s 
  type, extends(icBase) :: icPuff_s
    contains
      procedure :: init_GState => init_puffs
  end type icPuff_s 
  type, extends(icBase) :: icRandom_s
    contains
      procedure :: init_GState => init_randoms
  end type icRandom_s 
! type, extends(icBase) :: icUserdef_s
!   contains
!     procedure :: init_GState => init_userdefs
! end type icUserdef_s

CONTAINS

  ! ===================================================================
  ! Initial conditions
  ! ===================================================================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Read the scalar from restart files
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_reads(this,solver,state)
    use gstate_mod
    use hd_mod
    use iovar
    use status
    use filefmt
    use fft
    use commtypes
    implicit none

    class    (icRead_s), intent   (in)          :: this
    class(EquationBase), intent   (in)          :: solver
    type(GStateComp), intent(inout)             :: state(:)
    real   (kind=GP), pointer, dimension(:,:,:) :: R1
    integer                                     :: i
    logical                                     :: bret

    if ((stat .eq. 0).and.(solver%myrank_ .eq. 0)) then
      error stop 'Cannot read files if starting a new run with stat=0'
    endif
    call solver%workspace_%get_real_tmp(R1,bret)
    select type (solver)
    class is (VelocityBase)
      if ( solver%numpassive_ .eq. 0) then
        stop 'IC: Asking to read passive scalar ICs with npassive = 0'
      endif
      tind = int(stat)
      WRITE(ext, fmtext) tind
      do i = solver%PASSIVE, solver%PASSIVE+solver%numpassive_-1
        call io_read(1,solver%idir_,trim(solver%sstate_(i)),ext,solver%planio_,R1)
        call fftp3d_real_to_complex(planrc,R1,state(i)%ccomp,MPI_COMM_WORLD)
      end do
    class default
      error stop "This solver does not support passive scalars"
    end select
    call solver%workspace_%free_real_tmp(R1)
  end subroutine init_reads
 

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Constant concentration
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   c0  : vector with the amplitudes of the concentration
  subroutine init_constants(this,solver,state)
    use gstate_mod
    use hd_mod
    use grid
    use commtypes
    use mpivars
!$  use threads
    implicit none
    
    class(icConstant_s), intent   (in) :: this
    class(EquationBase), intent   (in) :: solver
    type   (GstateComp), intent(inout) :: state(:)
    real      (kind=GP), allocatable   :: c0(:)
    integer                            :: i,j,k,n
    
    namelist/ constant_s / c0
    select type (solver)
    class is (VelocityBase)
      if ( solver%numpassive_ .eq. 0) then
        stop 'IC: Asking to read passive scalar ICs with npassive = 0'
      endif
      allocate ( c0(solver%numpassive_) )
      if ( myrank .eq. 0 ) then
        open(1,file=solver%infile_,status='unknown',form="formatted")
        read(1,NML=constant_s)
        close(1)
      endif
      call mpi_bcast(c0,solver%numpassive_,GC_REAL,0,MPI_COMM_WORLD,ierr)
      do n = solver%PASSIVE, solver%PASSIVE+solver%numpassive_-1
        DO CONCURRENT (i=ista:iend, j=1:ny)
           DO CONCURRENT (k=1:nz)
             state(n)%ccomp(k,j,i) = 0.0_GP
           END DO
        END DO
        if ( myrank .eq. 0) then
          state(n)%ccomp(1,1,1) = c0(n-solver%PASSIVE+1)*           &
                real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP)
        endif
      end do
    class default
      error stop "IC: This solver does not support passive scalars"
    end select
  end subroutine init_constants


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Concentration in a Gaussian ball ('puff') centered at
  !! (x0, y0, z0) with a FWHM of r0
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   c0        : vector with the amplitudes of the concentration
  !   x0, y0, z0: vectors with the centers of the puffs
  !   r0        : vector with the radii of the puffs
  subroutine init_puffs(this,solver,state)
    use gstate_mod
    use hd_mod
    use grid
    use mpivars
    use commtypes
    use fft
    use pseudospec_scalar
!$  use threads
    implicit none
    
    class    (icPuff_s), intent   (in) :: this
    class(EquationBase), intent   (in) :: solver
    type   (GStateComp), intent(inout) :: state(:)
    real      (kind=GP), pointer       :: R1(:,:,:)
    real      (kind=GP), allocatable, dimension(:)  :: c0,x0,y0,z0,r0
    double precision                   :: tmp
    integer                            :: i,j,k,n
    logical                            :: bret

    namelist/ puff_s / c0,x0,y0,z0,r0
    CALL solver%workspace_%get_real_tmp(R1,bret)    
    select type (solver)
    class is (VelocityBase)
    if ( solver%numpassive_ .eq. 0) then
      stop 'IC: Asking for passive scalar ICs with npassive = 0'
    endif
    allocate ( c0(solver%numpassive_) )
    allocate ( x0(solver%numpassive_) )
    allocate ( y0(solver%numpassive_) )
    allocate ( z0(solver%numpassive_) )
    allocate ( r0(solver%numpassive_) )
    if ( myrank .eq. 0 ) then
      open(1,file=solver%infile_,status='unknown',form="formatted")
      read(1,NML=puff_s)
      close(1)
    endif
    call mpi_bcast(c0,solver%numpassive_,GC_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(x0,solver%numpassive_,GC_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(y0,solver%numpassive_,GC_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(z0,solver%numpassive_,GC_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(r0,solver%numpassive_,GC_REAL,0,MPI_COMM_WORLD,ierr)
    do n = solver%PASSIVE, solver%PASSIVE+solver%numpassive_-1    
      DO CONCURRENT (k=ksta:kend, j=1:ny) LOCAL(tmp)
         DO CONCURRENT (i=1:nx) LOCAL(tmp)
           tmp = (real(i-1,kind=GP)/real(nx-1,kind=GP)-x0(n-solver%PASSIVE+1))**2 &
               + (real(j-1,kind=GP)/real(ny-1,kind=GP)-y0(n-solver%PASSIVE+1))**2 &
               + (real(k-1,kind=GP)/real(nz-1,kind=GP)-z0(n-solver%PASSIVE+1))**2 
           R1(i,j,k) = exp(-tmp**2/r0(n-solver%PASSIVE+1)**2)
         END DO
      END DO
      call fftp3d_real_to_complex(planrc,R1,state(n)%ccomp,MPI_COMM_WORLD)
      call variance(state(n)%ccomp,tmp,1)
      call mpi_bcast(tmp,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      DO CONCURRENT (i=ista:iend, j=1:ny)
         DO CONCURRENT (k=1:nz)
           state(n)%ccomp(k,j,i) = state(n)%ccomp(k,j,i)* &
                        c0(n-solver%PASSIVE+1)/sqrt(tmp)
         END DO
      END DO
    end do
    class default
      error stop "IC: This solver does not support passive scalars"
    end select
    call solver%workspace_%free_real_tmp(R1)
  end subroutine init_puffs


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Random concentration
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   c0 : vector with the amplitudes of the concentration
  !   kdn: vector with minimum wave number for concentration
  !   kup: vector with maximum wave number for concentration
  subroutine init_randoms(this,solver,state)
    use gstate_mod
    use hd_mod
    use random
    use grid
    use mpivars
    use commtypes
    use status
    use fft
    use var
    use kes
    use pseudospec_scalar
!$  use threads
    implicit none
    
    class  (icRandom_s), intent   (in) :: this
    class(EquationBase), intent   (in) :: solver
    type   (GstateComp), intent(inout) :: state(:)
    real      (kind=GP), allocatable, dimension(:)  :: c0,kdn,kup
    real      (kind=GP)                             :: skup,skdn
    real      (kind=GP)                             :: dump,phase
    double precision                   :: tmp
    integer                            :: i,j,k,n
    logical                            :: bret

    namelist/ random_s / c0,kup,kdn
    select type (solver)
    class is (VelocityBase)
    if ( solver%numpassive_ .eq. 0) then
      stop 'IC: Asking for passive scalar ICs with npassive = 0'
    endif
    allocate ( c0 (solver%numpassive_) )
    allocate ( kdn(solver%numpassive_) )
    allocate ( kup(solver%numpassive_) )
    if ( myrank .eq. 0 ) then
      open(1,file=solver%infile_,status='unknown',form="formatted")
      read(1,NML=random_s)
      close(1)
    endif
    call mpi_bcast(c0 ,solver%numpassive_,GC_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(kdn,solver%numpassive_,GC_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(kup,solver%numpassive_,GC_REAL,0,MPI_COMM_WORLD,ierr)
    do n = solver%PASSIVE, solver%PASSIVE+solver%numpassive_-1
      skdn = kdn(n-solver%PASSIVE+1)
      skup = kup(n-solver%PASSIVE+1)
      IF (ista.eq.1) THEN
        state(n)%ccomp(1,1,1) = 0.
        DO j = 2,ny/2+1
          IF ((kk2(1,j,1).le.skup**2).and.(kk2(1,j,1).ge.skdn**2)) THEN
            dump = 1./sqrt(kk2(1,j,1))
            phase = 2*pi*randu(seed)
            state(n)%ccomp(1,j,1) = (COS(phase)+im*SIN(phase))*dump
            state(n)%ccomp(1,ny-j+2,1) = conjg(state(n)%ccomp(1,j,1))
          ELSE
            state(n)%ccomp(1,j,1) = 0.
            state(n)%ccomp(1,ny-j+2,1) = 0.
          ENDIF
        END DO
        DO k = 2,nz/2+1
          IF ((kk2(k,1,1).le.skup**2).and.(kk2(k,1,1).ge.skdn**2)) THEN
            dump = 1./sqrt(kk2(k,1,1))
            phase = 2*pi*randu(seed)
            state(n)%ccomp(k,1,1) = (COS(phase)+im*SIN(phase))*dump
            state(n)%ccomp(nz-k+2,1,1) = conjg(state(n)%ccomp(k,1,1))
          ELSE
            state(n)%ccomp(k,1,1) = 0.
            state(n)%ccomp(nz-k+2,1,1) = 0.
          ENDIF
        END DO
        DO j = 2,ny
          DO k = 2,nz/2+1
            IF ((kk2(k,j,1).le.skup**2).and.(kk2(k,j,1).ge.skdn**2)) THEN
              dump = 1./sqrt(kk2(k,j,1))
              phase = 2*pi*randu(seed)
              state(n)%ccomp(k,j,1) = (COS(phase)+im*SIN(phase))*dump
              state(n)%ccomp(nz-k+2,ny-j+2,1) = conjg(state(n)%ccomp(k,j,1))
            ELSE
              state(n)%ccomp(k,j,1) = 0.
              state(n)%ccomp(nz-k+2,ny-j+2,1) = 0.
            ENDIF
          END DO
        END DO
        DO i = 2,iend
          DO j = 1,ny
            DO k = 1,nz
              IF ((kk2(k,j,i).le.skup**2).and.(kk2(k,j,i).ge.skdn**2)) THEN
                dump = 1./sqrt(kk2(k,j,i))
                phase = 2*pi*randu(seed)
                state(n)%ccomp(k,j,i) = 2*(COS(phase)+im*SIN(phase))*dump
              ELSE
                state(n)%ccomp(k,j,i) = 0.
              ENDIF
            END DO
          END DO
        END DO
      ELSE
        DO i = ista,iend
          DO j = 1,ny
            DO k = 1,nz
              IF ((kk2(k,j,i).le.skup**2).and.(kk2(k,j,i).ge.skdn**2)) THEN
                dump = 1./sqrt(kk2(k,j,i))
                phase = 2*pi*randu(seed)
                state(n)%ccomp(k,j,i) = 2*(COS(phase)+im*SIN(phase))*dump
              ELSE
                state(n)%ccomp(k,j,i) = 0.
              ENDIF
            END DO
          END DO
        END DO
      ENDIF
      CALL variance(state(n)%ccomp,tmp,1)
      CALL MPI_BCAST(tmp,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      DO CONCURRENT (i=ista:iend, j=1:ny)
         DO CONCURRENT (k=1:nz)
           state(n)%ccomp(k,j,i) = state(n)%ccomp(k,j,i)* &
                        c0(n-solver%PASSIVE+1)/sqrt(tmp)
         END DO
      END DO
    end do
    class default
      error stop "IC: This solver does not support passive scalars"
    end select
  end subroutine init_randoms

end module ic_passive
