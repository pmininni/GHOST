! =====================================================================
! NAME       : ic_velocity.f90
! DESCRIPTION: Initial conditions for all VelocityBase solver classes.
!              These initial conditions can be used to initialize
!              velocities in other solver classes, using composition of
!              initial condition operators, as long as no correlation
!              between the fields is needed.
!
! Initial conditions avaliable:
!   read_v   : Reads velocity from an input file numbered by stat
!   null_v   : Null velocity field
!   tg_v     : Taylor-Green flow
!   abc_v    : ABC flow
!   random_v : Random (Pouquet-Patterson) velocity field
!
! DATE       : 01/16/26 (PDM)
! =====================================================================

module ic_velocity
  USE icbase_mod

  IMPLICIT NONE

  ! ================= Initial conditions supported ====================
  type, extends(icBase) :: read_v
    contains
      procedure :: init_GState => init_readv
  end type read_v
  type, extends(icBase) :: null_v
    contains
      procedure :: init_GState => init_nullv
  end type null_v
  type, extends(icBase) :: tg_v
    contains
      procedure :: init_GState => init_tgv
  end type tg_v
  type, extends(icBase) :: abc_v
    contains
      procedure :: init_GState => init_abcv
  end type abc_v
  type, extends(icBase) :: random_v
    contains
      procedure :: init_GState => init_randomv
  end type random_v
! type, extends(icBase) :: userdef_v
!   contains
!     procedure :: init_GState => init_userdefv
! end type userdef_v

CONTAINS

  ! ===================================================================
  ! Initial conditions
  ! ===================================================================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Read the velocity field from restart files
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_readv(this,solver,state)
    use gstate_mod
    use hd_mod
    use iovar
    use status
    use filefmt
    use fft
    use commtypes
    implicit none

    class      (read_v), intent   (in)          :: this
    class(EquationBase), intent   (in)          :: solver
    type   (GStateComp), intent(inout)          :: state(:)
    real   (kind=GP), pointer, dimension(:,:,:) :: R1
    integer                          :: i
    logical                          :: bret

    if ((stat .eq. 0).and.(solver%myrank_ .eq. 0)) then
      stop 'Cannot read files if starting a new run with stat=0'
    endif
    call solver%workspace_%get_real_tmp(R1,bret)
    select type (solver)
    class is (VelocityBase)
      tind = int(stat)
      WRITE(ext, fmtext) tind
      do i = solver%VELOCITY,solver%VELOCITY+solver%nc_-1
        call io_read(1,solver%idir_,trim(solver%sstate_(i)),ext,solver%planio_,R1)
        call fftp3d_real_to_complex(planrc,R1,state(i)%ccomp,MPI_COMM_WORLD)
      end do
    class default
      error stop "This solver does not support velocity field ICs"
    end select
    call solver%workspace_%free_real_tmp(R1)
  end subroutine init_readv
 

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Null velocity
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_nullv(this,solver,state)
    use gstate_mod
    use hd_mod
    use grid
    use mpivars
!$  use threads
    implicit none
    
    class      (null_v), intent   (in) :: this
    class(EquationBase), intent   (in) :: solver
    type   (GStateComp), intent(inout) :: state(:)
    integer                            :: i,j,k

    select type (solver)
    class is (VelocityBase)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
        DO j = 1,ny
          DO k = 1,nz
            state(solver%VELOCITY  )%ccomp(k,j,i) = 0.0_GP
            state(solver%VELOCITY+1)%ccomp(k,j,i) = 0.0_GP
            state(solver%VELOCITY+2)%ccomp(k,j,i) = 0.0_GP
          END DO
        END DO
      END DO
    class default
      error stop "This solver does not support velocity field ICs"
    end select
  end subroutine init_nullv


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Taylor-Green flow
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   u0  : amplitude of the velocity
  !   kdn : minimum wave number (rounded to next integer)
  !   kup : maximum wave number (rounded to next integer)
  subroutine init_tgv(this,solver,state)
    use gstate_mod
    use hd_mod
    use grid
    use boxsize
    use commtypes
    use fft
    use ali
    use var
    use pseudospec_norm
!$  use threads
    implicit none
    
    class        (tg_v), intent   (in)       :: this
    class(EquationBase), intent   (in)       :: solver
    type   (GStateComp), intent(inout)       :: state(:)
    real(kind=GP), pointer, dimension(:,:,:) :: R1,R2
    real(kind=GP)                            :: u0,kdn,kup
    integer                                  :: i,j,k,ki
    logical                                  :: bret

    namelist/ tg_v / u0,kdn,kup
    CALL solver%workspace_%get_real_tmp(R1,bret)
    CALL solver%workspace_%get_real_tmp(R2,bret)
    select type (solver)
    class is (VelocityBase)
    ! Read parameters from a namelist in the input file
    if ( myrank .eq. 0 ) then
      open(1,file=solver%infile_,status='unknown',form="formatted")
      read(1,NML=tg_v)
      close(1)
    endif
    call mpi_bcast(u0 ,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(kup,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(kdn,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
    ! Generate a superposition of TG flows
    IF ( abs(Lx-Ly).gt.tiny ) THEN
      IF (myrank.eq.0) error stop "TG initial conditions require Lx=Ly"
    ENDIF
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
    DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
       DO j = 1,ny
          DO k = 1,nz
            state(solver%VELOCITY+2)%ccomp(k,j,i) = 0.0_GP !vz
          END DO
       END DO
    END DO
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
    DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
       DO j = 1,ny
          DO i = 1,nx
          R1(i,j,k) = 0.0_GP
          R2(i,j,k) = 0.0_GP
          DO ki = INT(kdn),INT(kup)
             R1(i,j,k) = R1(i,j,k)+SIN(2*pi*ki*(real(i,kind=GP)-1)/ &
                      real(nx,kind=GP))*COS(2*pi*ki*(real(j,kind=GP)-1)/ &
                      real(ny,kind=GP))*COS(2*pi*ki*(real(k,kind=GP)-1)/ &
                      real(nz,kind=GP))
             R2(i,j,k) = R2(i,j,k)-COS(2*pi*ki*(real(i,kind=GP)-1)/ &
                      real(nx,kind=GP))*SIN(2*pi*ki*(real(j,kind=GP)-1)/ &
                      real(ny,kind=GP))*COS(2*pi*ki*(real(k,kind=GP)-1)/ &
                      real(nz,kind=GP))
          END DO
          END DO
       END DO
    END DO
    CALL fftp3d_real_to_complex(planrc,R1, & !vx
                   state(solver%VELOCITY  )%ccomp,MPI_COMM_WORLD)
    CALL fftp3d_real_to_complex(planrc,R2, & !vy
                   state(solver%VELOCITY+1)%ccomp,MPI_COMM_WORLD)
    CALL normalize(state(solver%VELOCITY  )%ccomp, &
                   state(solver%VELOCITY+1)%ccomp, &
                   state(solver%VELOCITY+2)%ccomp, &
                   u0,1,MPI_COMM_WORLD)
    class default
      error stop "This solver does not support velocity field ICs"
    end select
    CALL solver%workspace_%free_real_tmp(R1)
    CALL solver%workspace_%free_real_tmp(R2)
  end subroutine init_tgv


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! ABC flow
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   u0  : normalized amplitude of the velocity
  !   kdn : minimum wave number (rounded to next integer)
  !   kup : maximum wave number (rounded to next integer)
  !   A   : A amplitude (default =1)
  !   B   : B amplitude (default =1)
  !   C   : C amplitude (default =1)
  subroutine init_abcv(this,solver,state)
    use gstate_mod
    use hd_mod
    use grid
    use boxsize
    use commtypes
    use fft
    use ali
    use var
    use pseudospec_norm
!$  use threads
    implicit none

    class       (abc_v), intent   (in)       :: this
    class(EquationBase), intent   (in)       :: solver
    type   (GstateComp), intent(inout)       :: state(:)
    real(kind=GP), pointer, dimension(:,:,:) :: R1,R2,R3
    real(kind=GP)                            :: u0,kdn,kup
    real(kind=GP)                            :: A,B,C
    integer                                  :: i,j,k,ki
    logical                                  :: bret

    namelist/ abc_v / u0,kdn,kup,A,B,C
    CALL solver%workspace_%get_real_tmp(R1,bret)
    CALL solver%workspace_%get_real_tmp(R2,bret)
    CALL solver%workspace_%get_real_tmp(R3,bret)
    select type (solver)
    class is (VelocityBase)
    ! Read parameters from a namelist in the input file
    A = 1.0_GP ; B = 1.0_GP ; C = 1.0_GP
    if ( myrank .eq. 0 ) then
      open(1,file=solver%infile_,status='unknown',form="formatted")
      read(1,NML=abc_v)
      close(1)
    endif
    call mpi_bcast(u0 ,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(kup,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(kdn,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(A  ,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(B  ,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(C  ,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
    ! Generate a superposition of ABC flows
    if ( (abs(Lx-Ly).gt.tiny).or.(abs(Lx-Lz).gt.tiny) ) then
      if (myrank.eq.0) error stop 'ABC initial conditions require Lx=Ly=Lz'
    endif
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
    DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
      DO j = 1,ny
        DO i = 1,nx
          R1(i,j,k) = 0.
          R2(i,j,k) = 0.
          R3(i,j,k) = 0.
          DO ki = INT(kdn),INT(kup)
          R1(i,j,k) = R1(i,j,k)+(B*COS(2*pi*ki*(real(j,kind=GP)-1)/ &
               real(ny,kind=GP))+C*SIN(2*pi*ki*(real(k,kind=GP)-1)/ &
               real(nz,kind=GP)))/ki**2
          R2(i,j,k) = R2(i,j,k)+(A*SIN(2*pi*ki*(real(i,kind=GP)-1)/ &
               real(nx,kind=GP))+C*COS(2*pi*ki*(real(k,kind=GP)-1)/ &
               real(nz,kind=GP)))/ki**2 
          R3(i,j,k) = R3(i,j,k)+(A*COS(2*pi*ki*(real(i,kind=GP)-1)/ &
               real(nx,kind=GP))+B*SIN(2*pi*ki*(real(j,kind=GP)-1)/ &
               real(ny,kind=GP)))/ki**2
          END DO
        END DO
      END DO
    END DO
    call fftp3d_real_to_complex(planrc,R1, & !vx
                   state(solver%VELOCITY  )%ccomp,MPI_COMM_WORLD)
    call fftp3d_real_to_complex(planrc,R2, & !vy
                   state(solver%VELOCITY+1)%ccomp,MPI_COMM_WORLD)
    call fftp3d_real_to_complex(planrc,R3, & !vz
                   state(solver%VELOCITY+2)%ccomp,MPI_COMM_WORLD)
    call normalize(state(solver%VELOCITY  )%ccomp, &
                   state(solver%VELOCITY+1)%ccomp, &
                   state(solver%VELOCITY+2)%ccomp, &
                   u0,1,MPI_COMM_WORLD)
    class default
      error stop "This solver does not support velocity field ICs"
    end select
    call solver%workspace_%free_real_tmp(R1)
    call solver%workspace_%free_real_tmp(R2)
    call solver%workspace_%free_real_tmp(R3)
  end subroutine init_abcv

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Random flow: Superposition of harmonic modes with random
  !! phases correlated to give a tunable total relative
  !! helicity [following Pouquet & Patterson, JFM 1978] with
  !! a ~k^(-4) slope
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   u0   : velocity amplitude
  !   kdn  : minimum wave number (rounded to next integer)
  !   kup  : maximum wave number (rounded to next integer)
  !   alpha: sin(2*alpha) is the total relative helicity
  subroutine init_randomv(this,solver,state)
    use gstate_mod
    use hd_mod
    use random
    use grid
    use boxsize
    use commtypes
    use status
    use fft
    use ali
    use kes
    use pseudospec_norm
!$  use threads
    implicit none

    class    (random_v), intent   (in)          :: this
    class(EquationBase), intent   (in)          :: solver
    type   (GstateComp), intent(inout)          :: state(:)
    complex(kind=GP), pointer, dimension(:,:,:) :: C1,C2,C3,C4
    complex(kind=GP), pointer, dimension(:,:,:) :: C5,C6,C7,C8
    real(kind=GP)                               :: u0,kdn,kup
    real(kind=GP)                               :: alpha,a1,a2
    real(kind=GP)                               :: dump,phase
    integer                                     :: i,j,k
    logical                                     :: bret

    namelist/ random_v / u0,kdn,kup,alpha
    CALL solver%workspace_%get_complex_tmp(C1,bret)
    CALL solver%workspace_%get_complex_tmp(C2,bret)
    CALL solver%workspace_%get_complex_tmp(C3,bret)
    CALL solver%workspace_%get_complex_tmp(C4,bret)
    CALL solver%workspace_%get_complex_tmp(C5,bret)
    CALL solver%workspace_%get_complex_tmp(C6,bret)
    CALL solver%workspace_%get_complex_tmp(C7,bret)
    CALL solver%workspace_%get_complex_tmp(C8,bret)
    select type (solver)
    class is (VelocityBase)
    ! Read parameters from a namelist in the input file
    alpha = 0.0_GP
    if ( myrank .eq. 0 ) then
      open(1,file=solver%infile_,status='unknown',form="formatted")
      read(1,NML=random_v)
      close(1)
    endif
    call mpi_bcast(u0   ,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(kup  ,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(kdn  ,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(alpha,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
    ! Generate the first random velocity field
    IF (ista.eq.1) THEN
      C1(1,1,1) = 0.
      C2(1,1,1) = 0.
      C3(1,1,1) = 0.
      DO j = 2,ny/2+1
        IF ((kk2(1,j,1).le.kup**2).and.(kk2(1,j,1).ge.kdn**2)) THEN
          dump = 1./sqrt(kk2(1,j,1))**4
          phase = 2*pi*randu(seed)
          C1(1,j,1) = (COS(phase)+im*SIN(phase))*dump
          C1(1,ny-j+2,1) = conjg(C1(1,j,1))
          phase = 2*pi*randu(seed)
          C2(1,j,1) = (COS(phase)+im*SIN(phase))*dump
          C2(1,ny-j+2,1) = conjg(C2(1,j,1))
          phase = 2*pi*randu(seed)
          C3(1,j,1) = (COS(phase)+im*SIN(phase))*dump
          C3(1,ny-j+2,1) = conjg(C3(1,j,1))
        ELSE
          C1(1,j,1) = 0.
          C1(1,ny-j+2,1) = 0.
          C2(1,j,1) = 0.
          C2(1,ny-j+2,1) = 0.
          C3(1,j,1) = 0.
          C3(1,ny-j+2,1) = 0.
        ENDIF
      END DO
      DO k = 2,nz/2+1
        IF ((kk2(k,1,1).le.kup**2).and.(kk2(k,1,1).ge.kdn**2)) THEN
          dump = 1./sqrt(kk2(k,1,1))**4
          phase = 2*pi*randu(seed)
          C1(k,1,1) = (COS(phase)+im*SIN(phase))*dump
          C1(nz-k+2,1,1) = conjg(C1(k,1,1))
          phase = 2*pi*randu(seed)
          C2(k,1,1) = (COS(phase)+im*SIN(phase))*dump
          C2(nz-k+2,1,1) = conjg(C2(k,1,1))
          phase = 2*pi*randu(seed)
          C3(k,1,1) = (COS(phase)+im*SIN(phase))*dump
          C3(nz-k+2,1,1) = conjg(C3(k,1,1))
        ELSE
          C1(k,1,1) = 0.
          C1(nz-k+2,1,1) = 0.
          C2(k,1,1) = 0.
          C2(nz-k+2,1,1) = 0.
          C3(k,1,1) = 0.
          C3(nz-k+2,1,1) = 0.
        ENDIF
      END DO
      DO j = 2,ny
        DO k = 2,nz/2+1
          IF ((kk2(k,j,1).le.kup**2).and.(kk2(k,j,1).ge.kdn**2)) THEN
            dump = 1./sqrt(kk2(k,j,1))**4
            phase = 2*pi*randu(seed)
            C1(k,j,1) = (COS(phase)+im*SIN(phase))*dump
            C1(nz-k+2,ny-j+2,1) = conjg(C1(k,j,1))
            phase = 2*pi*randu(seed)
            C2(k,j,1) = (COS(phase)+im*SIN(phase))*dump
            C2(nz-k+2,ny-j+2,1) = conjg(C2(k,j,1))
            phase = 2*pi*randu(seed)
            C3(k,j,1) = (COS(phase)+im*SIN(phase))*dump
            C3(nz-k+2,ny-j+2,1) = conjg(C3(k,j,1))
          ELSE
            C1(k,j,1) = 0.
            C1(nz-k+2,ny-j+2,1) = 0.
            C2(k,j,1) = 0.
            C2(nz-k+2,ny-j+2,1) = 0.
            C3(k,j,1) = 0.
            C3(nz-k+2,ny-j+2,1) = 0.
          ENDIF
        END DO
      END DO
      DO i = 2,iend
        DO j = 1,ny
          DO k = 1,nz
            IF ((kk2(k,j,i).le.kup**2).and.(kk2(k,j,i).ge.kdn**2)) THEN
              dump = 1./sqrt(kk2(k,j,i))**4
              phase = 2*pi*randu(seed)
              C1(k,j,i) = 2*(COS(phase)+im*SIN(phase))*dump
              phase = 2*pi*randu(seed)
              C2(k,j,i) = 2*(COS(phase)+im*SIN(phase))*dump
              phase = 2*pi*randu(seed)
              C3(k,j,i) = 2*(COS(phase)+im*SIN(phase))*dump
            ELSE
              C1(k,j,i) = 0.
              C2(k,j,i) = 0.
              C3(k,j,i) = 0.
            ENDIF
          END DO
        END DO
      END DO
    ELSE
      DO i = ista,iend
        DO j = 1,ny
          DO k = 1,nz
            IF ((kk2(k,j,i).le.kup**2).and.(kk2(k,j,i).ge.kdn**2)) THEN
              dump = 1./sqrt(kk2(k,j,i))**4
              phase = 2*pi*randu(seed)
              C1(k,j,i) = 2*(COS(phase)+im*SIN(phase))*dump
              phase = 2*pi*randu(seed)
              C2(k,j,i) = 2*(COS(phase)+im*SIN(phase))*dump
              phase = 2*pi*randu(seed)
              C3(k,j,i) = 2*(COS(phase)+im*SIN(phase))*dump
            ELSE
              C1(k,j,i) = 0.
              C2(k,j,i) = 0.
              C3(k,j,i) = 0.
            ENDIF
          END DO
        END DO
      END DO
    ENDIF
    CALL rotor3(C2,C3,C4,1)
    CALL rotor3(C1,C3,C5,2)
    CALL rotor3(C1,C2,C6,3)
    CALL normalize(C4,C5,C6,u0,1,MPI_COMM_WORLD)
    ! C4,C5,C6 are the components of the first random vector (say v1)
    ! Repeating this procedure we create the second random vector (v2)
    IF (ista.eq.1) THEN
      C1(1,1,1) = 0.
      C2(1,1,1) = 0.
      C3(1,1,1) = 0.
      DO j = 2,ny/2+1
        IF ((kk2(1,j,1).le.kup**2).and.(kk2(1,j,1).ge.kdn**2)) THEN
          dump = 1./sqrt(kk2(1,j,1))**4
          phase = 2*pi*randu(seed)
          C1(1,j,1) = (COS(phase)+im*SIN(phase))*dump
          C1(1,ny-j+2,1) = conjg(C1(1,j,1))
          phase = 2*pi*randu(seed)
          C2(1,j,1) = (COS(phase)+im*SIN(phase))*dump
          C2(1,ny-j+2,1) = conjg(C2(1,j,1))
          phase = 2*pi*randu(seed)
          C3(1,j,1) = (COS(phase)+im*SIN(phase))*dump
          C3(1,ny-j+2,1) = conjg(C3(1,j,1))
        ELSE
          C1(1,j,1) = 0.
          C1(1,ny-j+2,1) = 0.
          C2(1,j,1) = 0.
          C2(1,ny-j+2,1) = 0.
          C3(1,j,1) = 0.
          C3(1,ny-j+2,1) = 0.
        ENDIF
      END DO
      DO k = 2,nz/2+1
        IF ((kk2(k,1,1).le.kup**2).and.(kk2(k,1,1).ge.kdn**2)) THEN
          dump = 1./sqrt(kk2(k,1,1))**4
          phase = 2*pi*randu(seed)
          C1(k,1,1) = (COS(phase)+im*SIN(phase))*dump
          C1(nz-k+2,1,1) = conjg(C1(k,1,1))
          phase = 2*pi*randu(seed)
          C2(k,1,1) = (COS(phase)+im*SIN(phase))*dump
          C2(nz-k+2,1,1) = conjg(C2(k,1,1))
          phase = 2*pi*randu(seed)
          C3(k,1,1) = (COS(phase)+im*SIN(phase))*dump
          C3(nz-k+2,1,1) = conjg(C3(k,1,1))
        ELSE
          C1(k,1,1) = 0.
          C1(nz-k+2,1,1) = 0.
          C2(k,1,1) = 0.
          C2(nz-k+2,1,1) = 0.
          C3(k,1,1) = 0.
          C3(nz-k+2,1,1) = 0.
        ENDIF
      END DO
      DO j = 2,ny
        DO k = 2,nz/2+1
          IF ((kk2(k,j,1).le.kup**2).and.(kk2(k,j,1).ge.kdn**2)) THEN
            dump = 1./sqrt(kk2(k,j,1))**4
            phase = 2*pi*randu(seed)
            C1(k,j,1) = (COS(phase)+im*SIN(phase))*dump
            C1(nz-k+2,ny-j+2,1) = conjg(C1(k,j,1))
            phase = 2*pi*randu(seed)
            C2(k,j,1) = (COS(phase)+im*SIN(phase))*dump
            C2(nz-k+2,ny-j+2,1) = conjg(C2(k,j,1))
            phase = 2*pi*randu(seed)
            C3(k,j,1) = (COS(phase)+im*SIN(phase))*dump
            C3(nz-k+2,ny-j+2,1) = conjg(C3(k,j,1))
          ELSE
            C1(k,j,1) = 0.
            C1(nz-k+2,ny-j+2,1) = 0.
            C2(k,j,1) = 0.
            C2(nz-k+2,ny-j+2,1) = 0.
            C3(k,j,1) = 0.
            C3(nz-k+2,ny-j+2,1) = 0.
          ENDIF
        END DO
      END DO
      DO i = 2,iend
        DO j = 1,ny
          DO k = 1,nz
            IF ((kk2(k,j,i).le.kup**2).and.(kk2(k,j,i).ge.kdn**2)) THEN
              dump = 1./sqrt(kk2(k,j,i))**4
              phase = 2*pi*randu(seed)
              C1(k,j,i) = 2*(COS(phase)+im*SIN(phase))*dump
              phase = 2*pi*randu(seed)
              C2(k,j,i) = 2*(COS(phase)+im*SIN(phase))*dump
              phase = 2*pi*randu(seed)
              C3(k,j,i) = 2*(COS(phase)+im*SIN(phase))*dump
            ELSE
              C1(k,j,i) = 0.
              C2(k,j,i) = 0.
              C3(k,j,i) = 0.
            ENDIF
          END DO
        END DO
      END DO
    ELSE
      DO i = ista,iend
        DO j = 1,ny
          DO k = 1,nz
            IF ((kk2(k,j,i).le.kup**2).and.(kk2(k,j,i).ge.kdn**2)) THEN
              dump = 1./sqrt(kk2(k,j,i))**4
              phase = 2*pi*randu(seed)
              C1(k,j,i) = 2*(COS(phase)+im*SIN(phase))*dump
              phase = 2*pi*randu(seed)
              C2(k,j,i) = 2*(COS(phase)+im*SIN(phase))*dump
              phase = 2*pi*randu(seed)
              C3(k,j,i) = 2*(COS(phase)+im*SIN(phase))*dump
            ELSE
              C1(k,j,i) = 0.
              C2(k,j,i) = 0.
              C3(k,j,i) = 0.
            ENDIF
          END DO
        END DO
      END DO
    ENDIF
    CALL rotor3(C2,C3,C7,1)
    CALL rotor3(C1,C3,C8,2)
    CALL rotor3(C1,C2,C3,3)
    CALL normalize(C7,C8,C3,u0,1,MPI_COMM_WORLD)
    ! So far, v1 = (C4,C5,C6) and v2 = (C7,C8,C3) are two random 
    ! normalized vectors in Fourier space. Correlating vectors 
    ! v1 and v2 we create the velocities vectors vx, vy, and vz.
    a1=SIN(alpha)
    a2=COS(alpha)
    ! vx: components y and z of u2 are calculated in order to get 
    ! the x component of curl(u2)
    DO i = ista,iend
      DO j = 1,ny
        DO k = 1,nz
          IF ((kk2(k,j,i).le.kup**2).and.(kk2(k,j,i).ge.kdn**2)) THEN
            C1(k,j,i) = a1*C5(k,j,i)+a2*C8(k,j,i)
            C2(k,j,i) = a1*C6(k,j,i)+a2*C3(k,j,i)
          ENDIF
        END DO
      END DO
    END DO
    CALL rotor3(C1,C2,C1,1)
    DO i = ista,iend
      DO j = 1,ny
        DO k = 1,nz
          IF ((kk2(k,j,i).le.kup**2).and.(kk2(k,j,i).ge.kdn**2)) THEN
            dump = 1./sqrt(kk2(k,j,i))
            state(solver%VELOCITY  )%ccomp(k,j,i) = a2*C4(k,j,i) +  &
                        a1*C7(k,j,i) + C1(k,j,i)*dump
          ENDIF
        END DO
      END DO
    END DO
    ! vy: components x and z of u2 are calculated in order to get 
    ! the y component of curl(u2)
    DO i = ista,iend
      DO j = 1,ny
        DO k = 1,nz
          IF ((kk2(k,j,i).le.kup**2).and.(kk2(k,j,i).ge.kdn**2)) THEN
            C1(k,j,i) = a1*C4(k,j,i)+a2*C7(k,j,i)
            C2(k,j,i) = a1*C6(k,j,i)+a2*C3(k,j,i)
          ENDIF
        END DO
      END DO
    END DO
    CALL rotor3(C1,C2,C1,2)
    DO i = ista,iend
      DO j = 1,ny
        DO k = 1,nz
          IF ((kk2(k,j,i).le.kup**2).and.(kk2(k,j,i).ge.kdn**2)) THEN
            dump = 1./sqrt(kk2(k,j,i))
            state(solver%VELOCITY+1)%ccomp(k,j,i) = a2*C5(k,j,i) +  &
                        a1*C8(k,j,i) + C1(k,j,i)*dump
          ENDIF
        END DO
      END DO
    END DO
    ! vz: components x and y of u2 are calculated in order to get
    ! the z component of curl(u2)
    DO i = ista,iend
      DO j = 1,ny
        DO k = 1,nz
          IF ((kk2(k,j,i).le.kup**2).and.(kk2(k,j,i).ge.kdn**2)) THEN
            C1(k,j,i) = a1*C4(k,j,i)+a2*C7(k,j,i)
            C2(k,j,i) = a1*C5(k,j,i)+a2*C8(k,j,i)
          ENDIF
        END DO
      END DO
    END DO
    CALL rotor3(C1,C2,C1,3)
    DO i = ista,iend
      DO j = 1,ny
        DO k = 1,nz
          IF ((kk2(k,j,i).le.kup**2).and.(kk2(k,j,i).ge.kdn**2)) THEN
            dump = 1./sqrt(kk2(k,j,i))
            state(solver%VELOCITY+2)%ccomp(k,j,i) = a2*C6(k,j,i) +  &
                        a1*C3(k,j,i) + C1(k,j,i)*dump
          ENDIF
        END DO
      END DO
    END DO
    call normalize(state(solver%VELOCITY  )%ccomp, &
                   state(solver%VELOCITY+1)%ccomp, &
                   state(solver%VELOCITY+2)%ccomp, &
                   u0,1,MPI_COMM_WORLD)
    class default
      error stop "This solver does not support velocity field ICs"
    end select
    call solver%workspace_%free_complex_tmp(C1)
    call solver%workspace_%free_complex_tmp(C2)
    call solver%workspace_%free_complex_tmp(C3)
    call solver%workspace_%free_complex_tmp(C4)
    call solver%workspace_%free_complex_tmp(C5)
    call solver%workspace_%free_complex_tmp(C6)
    call solver%workspace_%free_complex_tmp(C7)
    call solver%workspace_%free_complex_tmp(C8)
  end subroutine init_randomv

end module ic_velocity
