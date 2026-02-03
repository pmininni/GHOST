! =====================================================================
! NAME       : gexrk.f90
! DESCRIPTION: Forms class for all methods required to take a step
!              using explicit RK of specified order and number of
!              stages
!
! DATE       : 2/2/26 (DLR)
! =====================================================================

module gexrk_mod
  USE gstate_mod

  IMPLICIT NONE

  ! ================= Solver traits ===================================
  type, public  :: GexRKTraits
!   logical       :: bSSP       = .FALSE. ! strong stability flag?
    integer       :: norder     = 2  
    integer       :: nstage     = 2 
  end type


  ! ================= Stepper ==========================================
  ! Define class:
  type :: GExRKStepper
    ! Member data:
    logical              :: binit_=.false. ! is initialized?
    type  (GexRKTraits)  :: traits_
  CONTAINS
    procedure, public :: set_order     ! set order, nstages
    procedure, public :: step          ! take one timestep
    procedure, public :: init          ! initialize
    final             :: GExRKStepper_dtor
  end type GExRKStepper

CONTAINS

  ! ===================================================================
  ! Solver initialization, this is where parameter files are read
  ! ===================================================================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Subroutine to initialize the solver
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init(this)
    USE commtypes
    class  (GExRKStepper), intent (inout) :: this

    ! Temporary data to read from namelists:
    logical                    :: dorot
    integer                    :: npassive
    integer                    :: ierr
    real(kind=GP)              :: nu, bkappa, bvfreq, omegax, omegay, omegaz
    real(kind=GP), allocatable :: kappa(:)

    call MPI_COMM_SIZE(MPI_COMM_WORLD,this%nprocs_,ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,this%myrank_,ierr)

  end subroutine init

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Function to one step
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine step(this, time, uin, uf, dt, dudt) 
    use pseudospec_fluid
    use ali
    use kes
    use var
    use grid
    use mpivars
!$  use threads
    implicit none

    class (GExRKStepper), intent   (in)             :: this
    real   (kind=GP), intent   (in)             :: time, dt
    type    (GState), intent(inout), target     :: uin(:),uf(:)
    type    (GState), intent(inout)             :: dudt(:) 
    complex(kind=GP), pointer, dimension(:,:,:) :: fx,fy,fz,vx,vy,vz
    complex(kind=GP), pointer, dimension(:,:,:) :: fth,th
    complex(kind=GP), pointer, dimension(:,:,:) :: C1,C2,C3,C4,C5,C6
    complex(kind=GP), pointer, dimension(:,:,:) :: C7,C8
    real   (kind=GP)                            :: bkappa,bvfreq,nu
    real   (kind=GP)                            :: xmom,xtemp
    real   (kind=GP)                            :: omegax,omegay,omegaz
    integer                                     :: i,j,k
    logical                                     :: bret
       
    nu     = this%traits_%nu
    bkappa = this%traits_%bkappa
    bvfreq = this%traits_%bvfreq
    xmom   = this%traits_%xmom  * this%bvfreq
    xtemp  = this%traits_%xtemp * this%traits_%bvfreq

    CALL this%workspace_%get_complex_tmp(C1,bret)
    CALL this%workspace_%get_complex_tmp(C2,bret)
    CALL this%workspace_%get_complex_tmp(C3,bret)
    CALL this%workspace_%get_complex_tmp(C4,bret)
    CALL this%workspace_%get_complex_tmp(C5,bret)
    CALL this%workspace_%get_complex_tmp(C6,bret)
    CALL this%workspace_%get_complex_tmp(C7,bret)
    CALL this%workspace_%get_complex_tmp(C8,bret)

    vx  => uin(this%VELOCITY  )%ccomp
    vy  => uin(this%VELOCITY+1)%ccomp
    vz  => uin(this%VELOCITY+2)%ccomp
    th  => uin(this%TEMP)      %ccomp
    fx  => uf (this%VELOCITY  )%ccomp
    fy  => uf (this%VELOCITY+1)%ccomp
    fz  => uf (this%VELOCITY+2)%ccomp
    fth => uf (this%TEMP)      %ccomp
      
    call prodre3(vx,vy,vz,C4,C5,C6)                    ! w x v
    if ( this%traits_%dorot ) then
      omegax = this%traits_%omega(1)
      omegay = this%traits_%omega(2)
      omegaz = this%traits_%omega(3)
      call saxpby_c(C1, vz, 2*omegay, vy, -2.0*omegaz) ! 2 Omega x v
      call saxpby_c(C2, vx, 2*omegaz, vz, -2.0*omegax)
      call saxpby_c(C3, vy, 2*omegax, vx, -2.0*omegay)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      do i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
      do j = 1,ny
      do k = 1,nz
         C4(k,j,i) = C4(k,j,i) + C1(k,j,i) ! (w x v + 2 Omega x v)_x
         C5(k,j,i) = C5(k,j,i) + C2(k,j,i) ! (w x v + 2 Omega x v)_y
         C6(k,j,i) = C6(k,j,i) + C3(k,j,i) ! (w x v + 2 Omega x v)_z
      end do
      end do
      end do
    endif

!$omp parallel do if (iend-ista.ge.nth) private (j,k)
    do i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
    do j = 1,ny
    do k = 1,nz
       C6(k,j,i) = C6(k,j,i) - xmom*th(k,j,i) ! buoyancy term
    end do
    end do
    end do

    call nonlhd3(C4,C5,C6,C1,1)  ! -[(w + 2 Omega) x v + Grad p]_x
    call nonlhd3(C4,C5,C6,C2,2)  ! -[(w + 2 Omega) x v + Grad p]_y
    call nonlhd3(C4,C5,C6,C3,3)  ! -[(w + 2 Omega) x v + Grad p]_z
    call laplak3(vx,C4)          ! Del^2 vx
    call laplak3(vy,C5)          ! Del^2 vy
    call laplak3(vz,C6)          ! Del^2 vz

    CALL advect3(vx,vy,vz,th,C7) ! -(v.Grad) th
    call laplak3(th,C8)          ! Del^2 th

!$omp parallel do if (iend-ista.ge.nth) private (j,k)
    DO i = ista,iend               ! heat 'currrent':
!$omp parallel do if (iend-ista.lt.nth) private (k)
      DO j = 1,ny
        DO k = 1,nz
          C7(k,j,i) = C7(k,j,i) + xtemp*vz(k,j,i) ! add N vz term
        END DO
      END DO
    END DO


!$omp parallel do if (iend-ista.ge.nth) private (j,k)
    do i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
    do j = 1,ny
    do k = 1,nz
      if ((kn2(k,j,i).le.kmax).and.(kn2(k,j,i).ge.tiny)) then
        dudt(this%VELOCITY  )%ccomp(k,j,i) = nu*C4(k,j,i) + C1(k,j,i) + fx(k,j,i)
        dudt(this%VELOCITY+1)%ccomp(k,j,i) = nu*C5(k,j,i) + C2(k,j,i) + fy(k,j,i)
        dudt(this%VELOCITY+2)%ccomp(k,j,i) = nu*C6(k,j,i) + C3(k,j,i) + fz(k,j,i)
        dudt(this%TEMP)%ccomp(k,j,i) = bkappa*C8(k,j,i) + C7(k,j,i) + fth(k,j,i)
      else
        dudt(this%VELOCITY  )%ccomp(k,j,i) = 0.0_GP
        dudt(this%VELOCITY+1)%ccomp(k,j,i) = 0.0_GP
        dudt(this%VELOCITY+2)%ccomp(k,j,i) = 0.0_GP
        dudt(this%TEMP)      %ccomp(k,j,i) = 0.0_GP
      endif
    enddo
    enddo
    enddo

    CALL this%workspace_%free_complex_tmp(C1)
    CALL this%workspace_%free_complex_tmp(C2)
    CALL this%workspace_%free_complex_tmp(C3)
    CALL this%workspace_%free_complex_tmp(C4)
    CALL this%workspace_%free_complex_tmp(C5)
    CALL this%workspace_%free_complex_tmp(C6)
    CALL this%workspace_%free_complex_tmp(C7)
    CALL this%workspace_%free_complex_tmp(C8)

    ! Compute passive scalars:
    call this%rhs_passive(uin, uf, this%traits_%kappa, dudt)

  end subroutine dudt_impl


  ! ===================================================================
  ! Computation of global quantities and spectra
  ! ===================================================================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Function to compute and write global quantities
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine global_impl(this, uin, uf, t) 
    use pseudospec_hd
    use pseudospec_phd
    use status
    implicit none

    class (GExRKStepper), intent(in)             :: this
    type    (GState), intent(in), target        :: uin(:), uf(:)
    integer         , intent(in)                :: t
    complex(kind=GP), pointer, dimension(:,:,:) :: fs,fx,fy,fz,vx,vy,vz,th
    complex(kind=GP), pointer, dimension(:,:,:) :: c1,c2,c3
    real   (kind=GP)                            :: rmp,rmq
    integer                                     :: i

    CALL this%workspace_%get_complex_tmp(C1,bret)

    vx => uin(this%VELOCITY  )%ccomp
    vy => uin(this%VELOCITY+1)%ccomp
    vz => uin(this%VELOCITY+2)%ccomp
    th => uin(this%TEMP)      %ccomp
    fx => uf (this%VELOCITY  )%ccomp
    fy => uf (this%VELOCITY+1)%ccomp
    fz => uf (this%VELOCITY+2)%ccomp
    fs => uf (this%TEMP)      %ccomp


    CALL hdcheck(vx,vy,vz,fx,fy,fz,t,dt,1,0)
    CALL maxabs(vx,vy,vz,rmp,0)

    CALL pscheck(th,fs,t,dt)
    CALL derivk3(th,c1,1)
    CALL derivk3(th,c2,2)
    CALL derivk3(th,c3,3)
    CALL maxabs(c1,c2,c3,rmq,2)

    IF (myrank.eq.0) THEN
      OPEN(1,file='maximum.txt',position='append')
      WRITE(1,FMT='(E13.6,E13.6,E13.6)') (t-1)*dt,rmp,rmq
      CLOSE(1)
    ENDIF
    do i = this%PASSIVE, this%PASSIVE+this%numpassive_-1
      call pscheck(uin(i)%ccomp,uf(i)%ccomp,t,dt,trim(this%sstate_(i)))
    end do   

    CALL this%workspace_%free_complex_tmp(c1)
    CALL this%workspace_%free_complex_tmp(c2)
    CALL this%workspace_%free_complex_tmp(c3)

  end subroutine global_impl

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Function to compute and write spectra
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine spectra_impl(this, uin) 
    use pseudospec_aniso
    use pseudospec_scalar
    use pseudospec_anisca
    use filefmt
    use status
    use iovar
    implicit none

    class (GExRKStepper), intent(in)             :: this
    type    (GState), intent(in), target        :: uin(:)
    complex(kind=GP), pointer, dimension(:,:,:) :: vx,vy,vz,th
    complex(kind=GP), pointer, dimension(:,:,:) :: c1,c2,c3
    integer                                     :: i
    logical                                     :: bret

    WRITE(ext, fmtext) sind
    vx => uin(this%VELOCITY  )%ccomp
    vy => uin(this%VELOCITY+1)%ccomp
    vz => uin(this%VELOCITY+2)%ccomp
    th => uin(this%TEMP)      %ccomp
    call spectrum(vx,vy,vz,ext,1,1)
    call this%workspace_%get_complex_tmp(c1,bret)
    call this%workspace_%get_complex_tmp(c2,bret)
    call this%workspace_%get_complex_tmp(c3,bret)
    call gradre3(vx,vy,vz,c1,c2,c3)            ! Computes v.Grad(v)
    call entrans(vx,vy,vz,-c1,-c2,-c3,ext,1)   ! Writes the energy flux
    ! Uncomment the following line to compute the helicity flux
    ! call heltrans(vx,vy,vz,-c1,-c2,-c3,ext,1)
    if ( this%traits_%dorot ) then
      call specpara(vx,vy,vz,ext,1,1)
      call specperp(vx,vy,vz,ext,1,1)
      call spectrsc(th,ext,0)
      call specscpa(th,ext,0)
      call specscpe(th,ext,0)

      call entpara(vx,vy,vz,-c1,-c2,-c3,ext,1) ! Writes the energy flux
      call entperp(vx,vy,vz,-c1,-c2,-c3,ext,1) ! Writes the energy fluxq
      ! The following two lines compute anisotropic helicity fluxes
      ! call helpara(vx,vy,vz,-c1,-c2,-c3,ext,1)
      ! call helperp(vx,vy,vz,-c1,-c2,-c3,ext,1)
      ! Uncomment the following line to compute 2D spectra
      ! CALL spec2D(vx,vy,vz,ext,odir,1,1)
      ! Uncomment the following lines to compute spatio-temporal spectra
      ! CALL write_fourier(vx,'vx',ext,odir)
      ! CALL write_fourier(vy,'vy',ext,odir)
      ! CALL write_fourier(vz,'vz',ext,odir)

      ! Uncomment the following line to compute vert. spectrum of pot'l vorticity
      !  CALL spectpv(vx,vy,vz,th,ext)

    endif
    call this%workspace_%free_complex_tmp(c1)
    call this%workspace_%free_complex_tmp(c2)
    call this%workspace_%free_complex_tmp(c3)
    if ( this%numpassive_ .gt. 0) then
      do i = this%PASSIVE, this%PASSIVE+this%numpassive_-1
        call spectrsc(uin(i)%ccomp,ext,0,trim(this%sstate_(i)))
        if ( this%traits_%dorot ) then
          call specscpa(uin(i)%ccomp,ext,0,trim(this%sstate_(i)))
          call specscpe(uin(i)%ccomp,ext,0,trim(this%sstate_(i)))
        endif
      end do
    endif
  end subroutine spectra_impl
  

  ! ===================================================================
  ! Solver specific methods
  ! ===================================================================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Constructor
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine GExRKStepper_ctor(this, infile, workspace, plan)
    use iovar
    class  (GExRKStepper), intent(inout)      :: this
    type(GWorkspace) , intent(inout), target :: workspace
    type(ioplan)     , intent(inout), target :: plan
    character(len=*) , intent   (in)         :: infile
    this%infile_    =  infile    ! input file
    this%workspace_ => workspace
    this%planio_    => plan
    call this%init()
  end subroutine GExRKStepper_ctor


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Destructor
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine GExRKStepper_dtor(this) 
    type  (GExRKStepper), intent(inout) :: this
    if (associated(this%workspace_))   nullify(this%workspace_)
    if (associated(this%planio_))      nullify(this%planio_)
    if (allocated(this%sstate_))       deallocate(this%sstate_)
    if (allocated(this%traits_%kappa)) deallocate(this%traits_%kappa)
  end subroutine GExRKStepper_dtor


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Convert input state name to index in state vector
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine sstate2istate_impl(this, sstate, istate) 
    class (GExRKStepper)              , intent   (in) :: this
    character(len=8)              , intent   (in) :: sstate(:)
    integer         , allocatable , intent(inout) :: istate(:)
    integer                                       :: i,j
    if ( size(sstate) .ne. size(istate) ) then
      stop 'GExRKStepper::sstate2istate_impl: Incompatible sstate and istate'
    endif  
    do i = 1, size(sstate)
      istate(i) = -1 ! return unusable index
      do j = 1, size(this%sstate_)
        if ( sstate(i) .eq. this%sstate_(j) ) then
          istate(i) = j
        endif
      enddo
    enddo
  end subroutine sstate2istate_impl
  

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Get state variable names
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_sstate_impl(this, sstate) 
    class  (GExRKStepper)             , intent   (in) :: this
    character (len=8), allocatable, intent(inout) :: sstate(:)
    character(len=100)                            :: snum
    character(len=1)                              :: comp(3)
    integer                                       :: j

    if ( size(sstate) .lt. this%state_size() ) then
      deallocate(sstate)
      allocate(sstate,this%state_size())
    endif

    comp = ['x', 'y', 'z']
    do j = this%VELOCITY,this%VELOCITY+this%nc_-1
       sstate(j) = 'v' // comp(j-this%VELOCITY+1)
    enddo
    sstate(this%TEMP) = 'th' 
    do j = this%PASSIVE,this%PASSIVE+this%numpassive_-1
       write(snum,'(I0)') j-this%PASSIVE+1
       sstate(j) = 's' // trim(snum)
    enddo
  end subroutine get_sstate_impl
  

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Function to compute number of state members (equations)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PURE function state_size_impl(this) result(num)
    class(GExRKStepper), intent(in) :: this
    integer                     :: num
    num = this%nc_ + 1           ! # vel. components + theta
    num = num + this%numpassive_ ! # passive scalars
  end function state_size_impl

end module gexrk_mod
