!=================================================================
! GHOST computation of SGS  ML data model. Assumes
! trained model is in ONNx format, this class 
! uses the ECMWF Infero engine to make a pred-
! iction for all SGS terms 
!   (hat(N)(hat(u),hat(u)) - hat(N(u,u))
! given feature input data from truncated
! numerical simulation
! 
! Requires Infero built with the ONNx Runtime 
! backend and fckit available.
!
! 2025 D. Rosenberg
!      CIRA/ NOAA
!
!=================================================================
MODULE class_GSGSmodel
 !    USE kes
 !    USE mpivars
      USE fprecision
      USE commtypes
      USE fftplans
      USE inferof     , only: infero_model, &
          infero_initialise, infero_finalise, infero_check
      USE fckit_module, only: fckit_map, fckit_tensor_real32
 !    USE, intrinsic :: iso_c_binding
      USE  iso_c_binding

      IMPLICIT NONE
   !  INCLUDE 'mpif.h' 


      INTERFACE
        TYPE(C_PTR) FUNCTION allocate_c_array(size_in_bytes) BIND(C, NAME='malloc')
          USE iso_c_binding
          IMPLICIT NONE

          INTEGER(C_INT), INTENT(IN), VALUE :: size_in_bytes
        END FUNCTION allocate_c_array

        SUBROUTINE free_c_array(c_ptr_to_free) BIND(C, NAME='free')
          USE iso_c_binding
          IMPLICIT NONE

          TYPE(C_PTR), INTENT(IN), VALUE :: c_ptr_to_free
        END SUBROUTINE free_c_array
      END INTERFACE
!-----------------------------------------------------------------
!-----------------------------------------------------------------

      TYPE, PUBLIC :: GSGSmodelTraits

        INTEGER            :: nx, ny, nz
        INTEGER            :: nchannel

        CHARACTER(len=1024) :: model_path, model_type
        CHARACTER(len=1024) :: in_name, out_name
        CHARACTER(len=:), ALLOCATABLE :: yaml_config

      END TYPE GSGSmodelTraits
!-----------------------------------------------------------------
!-----------------------------------------------------------------


      TYPE, PUBLIC :: GSGSmodel
        PRIVATE
        ! Member data:
        INTEGER, DIMENSION(MPI_STATUS_SIZE)          :: istatus_
        INTEGER                                      :: myrank_,nprocs_
        INTEGER                                      :: nx,ny,nz,ntot
        INTEGER                                      :: comm_
        INTEGER                                      :: ista,iend,ksta,kend
        INTEGER                                      :: ierr_, rlen_
        INTEGER, ALLOCATABLE, DIMENSION(:)           :: sndtype_, rcvtype_
        TYPE(GSGSmodelTraits)                        :: modelTraits_

        REAL(KIND=GP), ALLOCATABLE, DIMENSION(:,:,:) :: kk2
        REAL(KIND=GP), ALLOCATABLE, DIMENSION    (:) :: kx,ky,kz
        TYPE(FFTPLAN), POINTER                       :: plancr, planrc

        ! Infero data:
        REAL(KIND=GP), POINTER, DIMENSION(:,:)       :: t_in_
        REAL(KIND=GP), POINTER, DIMENSION(:,:)       :: t_out_
        TYPE(C_PTR)                                  :: c_ptr_t_in_
        TYPE(C_PTR)                                  :: c_ptr_t_out_
        ! fckit wrappers and name->tensor maps
        TYPE(fckit_tensor_real32)                    :: tin_wrapped_, tout_wrapped_
        TYPE(fckit_map)                              :: imap_, omap_

        ! Inference model
        TYPE(infero_model)                           :: infmodel_

        CHARACTER(len=MPI_MAX_ERROR_STRING)          :: serr_

      CONTAINS
        ! Public methods:
        PROCEDURE,PUBLIC :: GSGS_ctor
        FINAL            :: GSGS_dtor
        PROCEDURE,PUBLIC :: sgs_model         => GSGS_compute_model

      END TYPE GSGSmodel

      PRIVATE :: GSGS_compute_model, GSGS_init_infero, GSGS_pack, GSGS_unpack 
      PRIVATE :: GSGS_real_exch_types, GSGS_real_exch


! Methods:
  CONTAINS

  SUBROUTINE GSGS_ctor(this, comm, ngrid, bnds, arbsz, Dk, plancr, planrc, modtraits)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  Main explicit constructor
!  ARGUMENTS:
!    this    : 'this' class instance
!    comm    : MPI communicator
!    ngrid   : array of size 3 giving grid size
!    bnds    : array of [ista, iend, ksta, kend]
!    arbsz   : arbitrary size flag (0, 1)
!    Dk      : array of size 3 giving Fourier shell widths
!    plancr,
!     planrc : FFT plans
!    modraits: GSGSmodelTraits structure
!-----------------------------------------------------------------
 !  USE var
 !  USE grid
 !  USE boxsize
 !  USE mpivars
 !  USE random
    USE commtypes
    USE fftplans

    IMPLICIT NONE
    CLASS(GSGSmodel)      ,INTENT(INOUT)     :: this
    INTEGER          ,INTENT   (IN)     :: comm
    INTEGER          ,INTENT   (IN)     :: arbsz
    INTEGER          ,INTENT   (IN)     :: ngrid(3)
    INTEGER          ,INTENT   (IN)     :: bnds(4)
    INTEGER                             :: aniso,i,j,k,n(3)
    INTEGER                             :: ierr
    INTEGER                             :: nprocs, myrank
    REAL(KIND=GP)    ,INTENT   (IN)     :: Dk(3)
    REAL(KIND=GP)                       :: rmp,rmq,rms
    TYPE(FFTPLAN)    ,INTENT   (IN), TARGET     :: plancr, planrc
    TYPE(GSGSmodelTraits)               :: modtraits

    this%comm_   = MPI_COMM_NULL
    this%nprocs_ = 0
    this%myrank_ = -1
    this%nx = ngrid(1)
    this%ny = ngrid(2)
    this%nz = ngrid(3)
!   CALL range(1,this%nx/2+1,this%nprocs_,this%myrank_,this%ista,this%iend)
!   CALL range(1,this%nz,this%nprocs_,this%myrank_,this%ksta,this%kend)
    this%ista = bnds(1)
    this%iend = bnds(2)
    this%ksta = bnds(3)
    this%kend = bnds(4)

    IF ( comm .NE. MPI_COMM_NULL ) THEN

    CALL MPI_COMM_SIZE(comm,nprocs,this%ierr_)
    CALL MPI_COMM_RANK(comm,myrank,this%ierr_)

    CALL MPI_COMM_DUP(comm, this%comm_, this%ierr_)
    IF (this%ierr_ .NE. MPI_SUCCESS) THEN
      CALL MPI_Error_string(this%ierr_, this%serr_, this%rlen_, this%ierr_)
      WRITE(*,*) myrank, ' GSGS_ctor:', TRIM(this%serr_(:this%rlen_))
      STOP
    END IF

    CALL MPI_COMM_SIZE(this%comm_,this%nprocs_,this%ierr_)
    IF (this%ierr_ .NE. MPI_SUCCESS) THEN
      CALL MPI_Error_string(this%ierr_, this%serr_, this%rlen_, this%ierr_)
      WRITE(*,*) 'GSGS_ctor:', TRIM(this%serr_(:this%rlen_))
      STOP
    END IF

    CALL MPI_COMM_RANK(this%comm_,this%myrank_,this%ierr_)
    IF (this%ierr_ .NE. MPI_SUCCESS) THEN
      CALL MPI_Error_string(this%ierr_, this%serr_, this%rlen_, this%ierr_)
      WRITE(*,*) 'GSGS_ctor:', TRIM(this%serr_(:this%rlen_))
      STOP
    END IF
    ENDIF

    this%nx = ngrid(1)
    this%ny = ngrid(2)
    this%nz = ngrid(3)
    this%ntot = ngrid(1)*ngrid(2)*ngrid(3)
!   CALL range(1,this%nx/2+1,this%nprocs_,this%myrank_,this%ista,this%iend)
!   CALL range(1,this%nz,this%nprocs_,this%myrank_,this%ksta,this%kend)
    this%ista = bnds(1)
    this%iend = bnds(2)
    this%ksta = bnds(3)
    this%kend = bnds(4)


!   n = ngrid
!   CALL fftp3d_create_plan_comm(this%planrc,n,FFTW_REAL_TO_COMPLEX, &
!                                FFTW_ESTIMATE, this%comm_)
!   CALL fftp3d_create_plan_comm(this%plancr,n,FFTW_COMPLEX_TO_REAL, &
!                                FFTW_ESTIMATE, this%comm_)

    this%plancr => plancr
    this%planrc => planrc

    ALLOCATE( this%kx(this%nx), this%ky(this%ny), this%kz(this%nz) )
    ALLOCATE( this%kk2(this%nz,this%ny,this%ista:this%iend) )
    ALLOCATE( this%sndtype_(0:nprocs-1) )
    ALLOCATE( this%rcvtype_(0:nprocs-1) )
    if ( arbsz .EQ. 1 ) THEN
      aniso = 1
    ELSE
      IF ((this%nx.ne.this%ny).or.(this%ny.ne.this%nz)) THEN
         aniso = 1
      ELSE
         aniso = 0
      ENDIF
    ENDIF

    DO i = 1,this%nx/2
       this%kx(i) = real(i-1,kind=GP)
       this%kx(i+this%nx/2) = real(i-this%nx/2-1,kind=GP)
     END DO
     DO j = 1,this%ny/2
        this%ky(j) = real(j-1,kind=GP)
        this%ky(j+this%ny/2) = real(j-this%ny/2-1,kind=GP)
     END DO
     IF (this%ny.eq.1) THEN
        this%ky(1) = 0.0_GP
     ENDIF
     DO k = 1,this%nz/2
        this%kz(k) = real(k-1,kind=GP)
        this%kz(k+this%nz/2) = real(k-this%nz/2-1,kind=GP)
     END DO
     IF (aniso.eq.1) THEN
        rmp = 1.0_GP/real(this%nx,kind=GP)**2
        rmq = 1.0_GP/real(this%ny,kind=GP)**2
        rms = 1.0_GP/real(this%nz,kind=GP)**2
     ELSE
        rmp = 1.0_GP
        rmq = 1.0_GP
        rms = 1.0_GP
     ENDIF

!$omp parallel do if (this%iend-this%ista.ge.nth) private (j,k)
     DO i = this%ista,this%iend
!$omp parallel do if (this%iend-this%ista.lt.nth) private (k)
        DO j = 1,this%ny
           DO k = 1,this%nz
              this%kk2(k,j,i) = rmp*this%kx(i)**2+rmq*this%ky(j)**2+rms*this%kz(k)**2
           END DO
        END DO
     END DO

     IF ( arbsz .eq. 1 ) THEN
       this%kx = this%kx*Dk(1)
       this%ky = this%ky*Dk(2)
       this%kz = this%kz*Dk(3)
     ENDIF

     IF (aniso.eq.1) THEN
!$omp parallel do if (this%iend-this%ista.ge.nth) private (j,k)
        DO i = this%ista,this%iend
!$omp parallel do if (this%iend-this%ista.lt.nth) private (k)
           DO j = 1,this%ny
              DO k = 1,this%nz
                 this%kk2(k,j,i) = this%kx(i)**2+this%ky(j)**2+this%kz(k)**2
              END DO
           END DO
        END DO
     ENDIF 

    CALL GSGS_init_infero(this, modtraits)
  !   write(*,*)this%myrank_, ' GSGS_ctor: ksta=', this%ksta, ' kend=', this%kend, ' ista=', this%ista, ' iend=', this%iend

    RETURN
  END SUBROUTINE GSGS_ctor
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GSGS_dtor(this)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  Main explicit destructor
!  ARGUMENTS:
!    this    : 'this' class instance
!-----------------------------------------------------------------

    IMPLICIT NONE
    TYPE(GSGSmodel)   ,INTENT(INOUT)             :: this
    INTEGER                                 :: j

    IF ( ALLOCATED    (this%kk2) ) DEALLOCATE   (this%kk2)
    IF ( ALLOCATED    (this%kx)  ) DEALLOCATE   (this%kx)
    IF ( ALLOCATED    (this%ky)  ) DEALLOCATE   (this%ky)
    IF ( ALLOCATED    (this%kz)  ) DEALLOCATE   (this%kz)
    IF ( ALLOCATED(this%sndtype_)) DEALLOCATE   (this%sndtype_)
    IF ( ALLOCATED(this%rcvtype_)) DEALLOCATE   (this%rcvtype_)

!   CALL fftp3d_destroy_plan(this%plancr)
!   CALL fftp3d_destroy_plan(this%planrc)

    ! Clean up Infero:
    CALL infero_check( this%infmodel_%free() )
    CALL this%tin_wrapped_%final()
    CALL this%tout_wrapped_%final()
    CALL this%imap_%final()
    CALL this%omap_%final()
    CALL infero_check( infero_finalise() )

    CALL free_c_array(this%c_ptr_t_in_)
    CALL free_c_array(this%c_ptr_t_out_)

    CALL MPI_Comm_free(this%comm_, this%ierr_)

    RETURN
  END SUBROUTINE GSGS_dtor
!-----------------------------------------------------------------
!-----------------------------------------------------------------



  SUBROUTINE GSGS_init_infero(this, modtraits)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  Performs inference of ML model and returns the SGS terms
!  ARGUMENTS:
!    this    : 'this' class instance
!-----------------------------------------------------------------
    USE  iso_c_binding
    IMPLICIT NONE
    TYPE(GSGSmodel)      ,INTENT(INOUT)           :: this
    INTEGER(C_INT)                                :: nb
    INTEGER                                       :: nc, nn
    TYPE(GSGSmodelTraits), INTENT  (IN)           :: modtraits

!   INTEGER, PARAMETER                            :: KIND(0.0D0)

    this%modelTraits_ = modtraits

    IF ( modtraits%nx .NE. this%nx &
    .OR. modtraits%ny .NE. this%ny &
    .OR. modtraits%nz .NE. this%nz ) THEN
      WRITE(*,*) 'GSGS_init_infero: Incompatible grid sizes'
      STOP 
    ENDIF
    
    CALL infero_check( infero_initialise() )
  
    ! Wrap Fortran arrays into fckit tensors and map them by layer names
    this%tin_wrapped_  = fckit_tensor_real32(this%t_in_)
    this%tout_wrapped_ = fckit_tensor_real32(this%t_out_)

    this%imap_ = fckit_map()
    CALL this%imap_%insert(trim(this%modelTraits_%in_name),  this%tin_wrapped_%c_ptr())

    this%omap_ = fckit_map()
    CALL this%omap_%insert(trim(this%modelTraits_%out_name), this%tout_wrapped_%c_ptr())

    ! Allocate C arrays 
    nb = this%modelTraits_%nchannel * this%modelTraits_%nx * this%modelTraits_%ny * this%modelTraits_%nz * SIZEOF(1.0_GP)
    this%c_ptr_t_in_  = allocate_c_array(nb)

    nb = 3 * this%modelTraits_%nx * this%modelTraits_%ny * this%modelTraits_%nz * SIZEOF(1.0_GP)
    this%c_ptr_t_out_ = allocate_c_array(nb)

    ! Associate Fortran pointers with C memory:
    nc = this%modelTraits_%nchannel
    nn = this%ntot
    CALL C_F_POINTER(this%c_ptr_t_in_ , this%t_in_ , SHAPE=[nc, nn])
    CALL C_F_POINTER(this%c_ptr_t_out_, this%t_out_, SHAPE=[3 , nn])

    ! Build a YAML config to point Infero at the 
    ! model and backend type:
    this%modelTraits_%yaml_config = &
        '---'//new_line('a')                            // &
        '  path: ' // trim(this%modelTraits_%model_path) // new_line('a') // &
        '  type: ' // trim(this%modelTraits_%model_type) // c_null_char

    ! Initialize model from YAML string
    CALL infero_check( this%infmodel_%initialise_from_yaml_string(this%modelTraits_%yaml_config) )


    RETURN

  END SUBROUTINE GSGS_init_infero
!-----------------------------------------------------------------
!-----------------------------------------------------------------


  SUBROUTINE GSGS_compute_model(this, vx, vy, vz, th, C1, R1,  SGS1, SGS2, SGS3, SGSth)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  Performs inference of ML model and returns the SGS terms
!  ARGUMENTS:
!    this    : 'this' class instance
!    vx,vy,vz: input velocities
!    th      : input pot'l temp
!    C1      : complex tmp array(s)
!    R1      : real tmp array(s)
!    SGSi    : output SGS components
!-----------------------------------------------------------------

    IMPLICIT NONE
    class(GSGSmodel),INTENT (INOUT)         :: this
    COMPLEX(KIND=GP), INTENT   (IN), DIMENSION(this%nz,this%ny,this%ista:this%iend) :: vx,vy,vz,th
    COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(this%nz,this%ny,this%ista:this%iend) :: C1
    REAL(KIND=GP)   , INTENT(INOUT), DIMENSION(this%nx,this%ny,this%ksta:this%kend) :: R1
    COMPLEX(KIND=GP), INTENT  (OUT), DIMENSION(this%nz,this%ny,this%ista:this%iend) :: SGS1,SGS2,SGS3,SGSth

    ! Pack model input layer:
    ! shape = ("time", "channel", "x0", "x1", "x2")
    CALL GSGS_pack(this, vx, 0, C1, R1, this%t_in_)
    CALL GSGS_pack(this, vy, 1, C1, R1, this%t_in_)
    CALL GSGS_pack(this, vz, 2, C1, R1, this%t_in_)
    CALL GSGS_pack(this, th, 3, C1, R1, this%t_in_)

    ! Run inference
    CALL infero_check( this%infmodel_%infer(this%imap_, this%omap_) )

    ! (Optional) print stats/config
    CALL infero_check( this%infmodel_%print_statistics() )
    CALL infero_check( this%infmodel_%print_config() )

    ! Show output
!   PRINT *, 'Inference output (first batch row):'
!   PRINT '(100(1x,f10.6))', this%t_out_(1, :)

    ! Unpack model output and compute FFTs:
    CALL GSGS_unpack(this, this%t_out_, 0, R1, SGS1)
    CALL GSGS_unpack(this, this%t_out_, 1, R1, SGS2)
    CALL GSGS_unpack(this, this%t_out_, 2, R1, SGS3)


    RETURN

  END SUBROUTINE GSGS_compute_model
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GSGS_pack(this, cvar, ivar, C1, R1, itensor)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  Packs tensor argument for inference
!  ARGUMENTS:
!    cvar    : channel/feature to pack into itensor
!    ivar    : channel/feature  id (0, nchannel-1)
!    C1      : complex tmp array
!    R1      : real tmp array
!    itensor : input tensor to pack into
!-----------------------------------------------------------------
    IMPLICIT NONE
    class(GSGSmodel), INTENT(INOUT)         :: this
    INTEGER         , INTENT   (IN)         :: ivar
    INTEGER                                 :: i, j, k
    COMPLEX(KIND=GP), INTENT   (IN), DIMENSION(this%nz,this%ny,this%ista:this%iend) :: cvar
    COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(this%nz,this%ny,this%ista:this%iend) :: C1
    REAL(KIND=GP)   , INTENT(INOUT), DIMENSION(this%nx,this%ny,this%ksta:this%kend) :: R1
    REAL(KIND=GP)   , INTENT  (OUT), DIMENSION(this%modelTraits_%nchannel,this%nx,this%ny,this%ksta:this%kend) :: itensor

!$omp parallel do if (this%iend-2.ge.nth) private (j,k)
    DO i = this%ista,this%iend
!$omp parallel do if (this%iend-2.lt.nth) private (k)
       DO j = 1,this%ny
          DO k = 1,this%nz
             C1(k,j,i) = cvar(k,j,i)
          END DO
       END DO
    END DO

    CALL fftp3d_complex_to_real(this%plancr,C1,R1)

    CALL GSGS_real_exch(R1, ivar, itensor) 

#if 0
!$omp parallel do if (this%kend-this%ksta.ge.nth) private (j,i)
    DO k = this%ksta,this%kend
!$omp parallel do if (this%kend-this%ksta.lt.nth) private (i)
       DO j = 1,this%ny
          DO i = 1,this%nx
             itensor(ivar+1,i,j,k) = R1(i,j,k)
          END DO
       END DO
    END DO
#endif
       
    RETURN
  END SUBROUTINE GSGS_pack
!-----------------------------------------------------------------
!-----------------------------------------------------------------


  SUBROUTINE GSGS_unpack(this, otensor, ivar, R1, sgs)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  Packs tensor argument for inference
!  ARGUMENTS:
!    otensor : tensor output to read
!    ivar    : channel/feature  id (0, nchannel-1)
!    C1      : complex tmp array
!    R1      : real tmp array
!    sgs     : return complex SGS data
!-----------------------------------------------------------------
    IMPLICIT NONE      
    class(GSGSmodel), INTENT(INOUT)         :: this
    INTEGER         , INTENT   (IN)         :: ivar
    INTEGER                                 :: i, j, k
    COMPLEX(KIND=GP), INTENT   (IN), DIMENSION(this%nz,this%ny,this%ista:this%iend) :: sgs
    REAL(KIND=GP)   , INTENT(INOUT), DIMENSION(this%nx,this%ny,this%ksta:this%kend) :: R1
    REAL(KIND=GP)   , INTENT  (OUT), DIMENSION(3, this%nx,this%ny,this%ksta:this%kend) :: otensor
    REAL(KIND=GP)                            :: tmp

    tmp = 1.0_GP/ &
            (real(this%nx,kind=GP)*real(this%ny,kind=GP)*real(this%nz,kind=GP))

!$omp parallel do if (this%kend-this%ksta.ge.nth) private (j,i)
    DO k = this%ksta,this%kend
!$omp parallel do if (this%kend-this%ksta.lt.nth) private (i)
       DO j = 1,this%ny
          DO i = 1,this%nx
             R1(i,j,k) = otensor(ivar+1,i,j,k) * tmp
          END DO
       END DO
    END DO

    CALL fftp3d_real_to_complex(this%planrc,R1,sgs)
       

    RETURN
  END SUBROUTINE GSGS_unpack
!-----------------------------------------------------------------
!-----------------------------------------------------------------

      SUBROUTINE GSGS_real_exch_types(n,nprocs,myrank,sndtype,rcvtype)
!-----------------------------------------------------------------
!
! Defines derived data types for sending and receiving 
! blocks of the 3D matrix between processors. The data 
! types are used to transpose the matrix during the FFT.
!
! Parameters
!     n      : the size of the dimensions of the input array [IN]
!     nprocs : the number of processors [IN]
!     myrank : the rank of the processor [IN]
!     sndtype: contains a derived data type for sending [OUT]
!     rcvtype: contains a derived data type for receiving [OUT]
!-----------------------------------------------------------------
      USE commtypes
      IMPLICIT NONE

      INTEGER, INTENT(OUT), DIMENSION(0:nprocs-1) :: sndtype,rcvtype
      INTEGER, INTENT(IN) :: n(3),nprocs
      INTEGER, INTENT(IN) :: myrank

      INTEGER :: ista,iend
      INTEGER :: ksta,kend
      INTEGER :: irank,krank
      INTEGER :: itemp1,itemp2

      CALL range(1,n(3),nprocs,myrank,ksta,kend)
      DO irank = 0,nprocs-1
         ista = 1
         iend = n(1)
         CALL block3d(1,n(1),1,n(2),ksta,ista,iend,1,n(2), &
                     ksta,kend,GC_REAL,itemp1)
         sndtype(irank) = itemp1
      END DO
      CALL range(1,n(3),nprocs,myrank,ksta,kend)
      DO krank = 0,nprocs-1
         ista = 1
         iend = n(1)
         CALL block3d(1,n(1),1,n(2),ksta,ista,iend,1,n(2), &
                     ksta,kend,GC_REAL,itemp1)
         rcvype(krank) = itemp2
      END DO

      RETURN
      END SUBROUTINE GSGS_real_exch_types
!-----------------------------------------------------------------
!-----------------------------------------------------------------

      SUBROUTINE GSGS_real_exch(R1, ivar, t_in)
!-----------------------------------------------------------------
!
! Parameters
!     R1     : Real field
!     ivar   : which channel/feature (0, ... nchannel-1)
!     t_in   : input tensor that will contain all other
!              tasks' field data
!-----------------------------------------------------------------

      USE commtypes
      IMPLICIT NONE

      INTEGER         , INTENT   (IN)     :: ivar
      REAL(KIND=GP)   , INTENT   (IN), DIMENSION(this%nx,this%ny,this%ksta:this%kend) :: R1
      REAL(KIND=GP)   , INTENT(INOUT), DIMENSION(this%modelTraits_%nchannel,this%ntot):: t_in

      INTEGER                             :: iprocs, irank, istrip, nstrip
      INTEGER                             :: isendTo, igetFrom
      INTEGER, DIMENSION(0:this%nprocs_-1) :: ireq1,ireq2
      INTEGER, DIMENSION(MPI_STATUS_SIZE) :: istatus

      nstrip = 1
      do iproc = 0, this%nprocs_-1, nstrip
         do istrip=0, nstrip-1
            irank = iproc + istrip

            isendTo = this%myrank_ + irank
            if ( isendTo .ge. this%nprocs_ ) isendTo = isendTo - this%nprocs_

            igetFrom = this%myrank_ - irank
            if ( igetFrom .lt. 0 ) igetFrom = igetFrom + this%nprocs_
            CALL MPI_IRECV(t_in(ivar+1,:),1,this%rcvtype(igetFrom),igetFrom,      &
                          1,this%comm_,ireq2(irank),this%ierr_)

            CALL MPI_ISEND(R1,1,this%sndtype(isendTo),isendTo, &
                          1,this%comm_,ireq1(irank),this%ierr_)
         enddo

         do istrip=0, nstrip-1
            irank = iproc + istrip
            CALL MPI_WAIT(ireq1(irank),istatus,this%ierr_)
            CALL MPI_WAIT(ireq2(irank),istatus,this%ierr_)
         enddo
      enddo

      RETURN
      END SUBROUTINE GSGS_real_exch
!-----------------------------------------------------------------
!-----------------------------------------------------------------



END MODULE class_GSGSmodel
