! =====================================================================
! NAME       : ic_magnetic.f90
! DESCRIPTION: Initial conditions for all MagneticBase solver classes.
!              These initial conditions can be used to initialize
!              magnetic fields in other solver classes, using
!              composition of initial condition operators, as long as
!              no correlation between the fields is needed. For v-b
!              correlated cases, this file can provide methods to
!              initialize both fields simultaneously.
!
! Initial conditions avaliable:
!   read_b   : Reads vector potential from input file numbered by stat
!   null_b   : Null vector potential
!   random_b : Random (Pouquet-Patterson) vector potential
!
! DATE       : 01/22/26 (PDM)
! =====================================================================

module ic_magnetic
  USE icbase_mod

  IMPLICIT NONE

  ! ================= Initial conditions supported ====================
  type, extends(icBase) :: read_b
    contains
      procedure :: init_GState => init_readb
  end type read_b 
  type, extends(icBase) :: null_b
    contains
      procedure :: init_GState => init_nullb
  end type null_b
  type, extends(icBase) :: random_b
    contains
      procedure :: init_GState => init_randomb
  end type random_b
! type, extends(icBase) :: userdef_b
!   contains
!     procedure :: init_GState => init_userdefb
! end type random_b

CONTAINS

  ! ===================================================================
  ! Initial conditions
  ! ===================================================================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Read the vector potential from restart files
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_readb(this,solver,state)
    use gstate_mod
    use hd_mod
    use iovar
    use status
    use filefmt
    use fft
    use commtypes
!$  use threads
    implicit none

    class      (read_b), intent   (in)          :: this
    class(EquationBase), intent   (in)          :: solver
    type    (GState), intent(inout)             :: state(:)
    complex(kind=GP), pointer, dimension(:,:,:) :: C1
    real   (kind=GP), pointer, dimension(:,:,:) :: R1
    integer                          :: i
    logical                          :: bret

    if ((stat .eq. 0).and.(solver%myrank_ .eq. 0)) then
       error stop 'Cannot read files if starting a new run with stat=0'
    endif
    call solver%workspace_%get_complex_tmp(C1,bret)
    call solver%workspace_%get_real_tmp(   R1,bret)
    select type (solver)
    class is (MagneticBase)
      tind = int(stat)
      WRITE(ext, fmtext) tind
      do i = solver%MAGNETIC,solver%MAGNETIC+solver%nc_-1
        call io_read(1,idir,trim(solver%sstate_(i)),ext,solver%planio_,R1)
        call fftp3d_real_to_complex(planrc,R1,state(i)%ccomp,MPI_COMM_WORLD)
      end do
    class default
      error stop "This solver does not support magnetic field ICs"
    end select
    call solver%workspace_%free_complex_tmp(C1)
    call solver%workspace_%free_real_tmp(   R1)
  end subroutine init_readb
 

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Null vector potential
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_nullb(this,solver,state)
    use gstate_mod
    use hd_mod
    use fprecision
    use grid
    use mpivars
!$  use threads
    implicit none
    
    class      (null_b), intent   (in) :: this
    class(EquationBase), intent   (in) :: solver
    type       (Gstate), intent(inout) :: state(:)
    integer                            :: i,j,k

    select type (solver)
    class is (MagneticBase)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
        DO j = 1,ny
          DO k = 1,nz
            state(solver%MAGNETIC  )%ccomp(k,j,i) = 0.0_GP
            state(solver%MAGNETIC+1)%ccomp(k,j,i) = 0.0_GP
            state(solver%MAGNETIC+2)%ccomp(k,j,i) = 0.0_GP
          END DO
        END DO
      END DO
    class default
      error stop "This solver does not support magnetic field ICs"
    end select
  end subroutine init_nullb

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Random vector potential: Superposition of harmonic modes
  !! with random phases correlated to give a tunable total
  !! relative helicity [following Pouquet & Patterson, JFM
  !! 1978] with a ~k^(-4) slope
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   b0   : magnetic field amplitude
  !   kdn  : minimum wave number (rounded to next integer)
  !   kup  : maximum wave number (rounded to next integer)
  !   alpha: sin(2*alpha) is the total relative helicity
  subroutine init_randomb(this,solver,state)
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

    class    (random_b), intent   (in)          :: this
    class(EquationBase), intent   (in)          :: solver
    type       (Gstate), intent(inout)          :: state(:)
    complex(kind=GP), pointer, dimension(:,:,:) :: C1,C2,C3,C4
    complex(kind=GP), pointer, dimension(:,:,:) :: C5,C6,C7,C8
    real(kind=GP)                               :: b0,kdn,kup
    real(kind=GP)                               :: alpha,a1,a2
    real(kind=GP)                               :: dump,phase
    integer                                     :: i,j,k
    logical                                     :: bret

    namelist/ random_b / b0,kdn,kup,alpha
    CALL solver%workspace_%get_complex_tmp(C1,bret)
    CALL solver%workspace_%get_complex_tmp(C2,bret)
    CALL solver%workspace_%get_complex_tmp(C3,bret)
    CALL solver%workspace_%get_complex_tmp(C4,bret)
    CALL solver%workspace_%get_complex_tmp(C5,bret)
    CALL solver%workspace_%get_complex_tmp(C6,bret)
    CALL solver%workspace_%get_complex_tmp(C7,bret)
    CALL solver%workspace_%get_complex_tmp(C8,bret)
    select type (solver)
    class is (MagneticBase)
    ! Read parameters from a namelist in the input file
    alpha = 0.0_GP
    if ( myrank .eq. 0 ) then
      open(1,file=solver%infile_,status='unknown',form="formatted")
      read(1,NML=random_b)
      close(1)
    endif
    call mpi_bcast(b0   ,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(kup  ,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(kdn  ,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(alpha,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
    ! Generate the first random field
    IF (ista.eq.1) THEN
      C1(1,1,1) = 0.
      C2(1,1,1) = 0.
      C3(1,1,1) = 0.
      DO j = 2,ny/2+1
        IF ((kk2(1,j,1).le.kup**2).and.(kk2(1,j,1).ge.kdn**2)) THEN
          dump = 1./sqrt(kk2(1,j,1))**5
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
          dump = 1./sqrt(kk2(k,1,1))**5
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
            dump = 1./sqrt(kk2(k,j,1))**5
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
              dump = 1./sqrt(kk2(k,j,i))**5
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
              dump = 1./sqrt(kk2(k,j,i))**5
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
    CALL normalize(C4,C5,C6,b0,0,MPI_COMM_WORLD)
    ! C4,C5,C6 are the components of the first random vector (say v1)
    ! Repeating this procedure we create the second random vector (v2)
    IF (ista.eq.1) THEN
      C1(1,1,1) = 0.
      C2(1,1,1) = 0.
      C3(1,1,1) = 0.
      DO j = 2,ny/2+1
        IF ((kk2(1,j,1).le.kup**2).and.(kk2(1,j,1).ge.kdn**2)) THEN
          dump = 1./sqrt(kk2(1,j,1))**5
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
          dump = 1./sqrt(kk2(k,1,1))**5
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
            dump = 1./sqrt(kk2(k,j,1))**5
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
              dump = 1./sqrt(kk2(k,j,i))**5
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
              dump = 1./sqrt(kk2(k,j,i))**5
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
    CALL normalize(C7,C8,C3,b0,0,MPI_COMM_WORLD)
    ! So far, v1 = (C4,C5,C6) and v2 = (C7,C8,C3) are two random 
    ! normalized vectors in Fourier space. Correlating vectors 
    ! v1 and v2 we create the vector fields vx, vy, and vz.
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
            state(solver%MAGNETIC  )%ccomp(k,j,i) = a2*C4(k,j,i) +  &
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
            state(solver%MAGNETIC+1)%ccomp(k,j,i) = a2*C5(k,j,i) +  &
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
            state(solver%MAGNETIC+2)%ccomp(k,j,i) = a2*C6(k,j,i) +  &
                        a1*C3(k,j,i) + C1(k,j,i)*dump
          ENDIF
        END DO
      END DO
    END DO
    call normalize(state(solver%MAGNETIC  )%ccomp, &
                   state(solver%MAGNETIC+1)%ccomp, &
                   state(solver%MAGNETIC+2)%ccomp, &
                   b0,0,MPI_COMM_WORLD)
    class default
      error stop "This solver does not support magnetic field ICs"
    end select
    call solver%workspace_%free_complex_tmp(C1)
    call solver%workspace_%free_complex_tmp(C2)
    call solver%workspace_%free_complex_tmp(C3)
    call solver%workspace_%free_complex_tmp(C4)
    call solver%workspace_%free_complex_tmp(C5)
    call solver%workspace_%free_complex_tmp(C6)
    call solver%workspace_%free_complex_tmp(C7)
    call solver%workspace_%free_complex_tmp(C8)
  end subroutine init_randomb

end module ic_magnetic

