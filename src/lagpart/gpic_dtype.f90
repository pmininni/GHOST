!=================================================================
! GHOST particle-in-cell subclasses
!
! 2023 F. Pugliese
!      Department of Physics,
!      Facultad de Ciencias Exactas y Naturales.
!      Universidad de Buenos Aires.
!=================================================================

  TYPE, PUBLIC, EXTENDS ( GPart ) :: GPIC

        PRIVATE
        ! Member data:
!        TYPE(GFieldComm)                             :: gfcomm_
        TYPE(GPICSplineInt)                          :: picspl_
        REAL(KIND=GP), ALLOCATABLE, DIMENSION    (:) :: prop_,weight_
        REAL(KIND=GP)                                :: icv_
      CONTAINS
        ! Public methods:
        PROCEDURE,PUBLIC :: GPIC_ctor
        PROCEDURE,PUBLIC :: Init              => GPIC_Init
        PROCEDURE,PUBLIC :: GPICStep          => GPIC_StepRKK
        PROCEDURE,PUBLIC :: GetDensity        => GPIC_GetDensity
        PROCEDURE,PUBLIC :: io_write_wgt      => GPIC_io_write_wgt
        PROCEDURE,PUBLIC :: io_read_wgt       => GPIC_io_read_wgt
        PROCEDURE,PUBLIC :: InitUserSeed      => GPIC_InitUserSeed
        PROCEDURE,PUBLIC :: ResizeArrays      => GPIC_ResizeArrays
      END TYPE GPIC

  TYPE, PUBLIC, EXTENDS ( GPIC ) :: VGPIC
        PRIVATE
        ! Member data:
        REAL(KIND=GP), ALLOCATABLE, DIMENSION    (:) :: pvx_,pvy_,pvz_
        REAL(KIND=GP), ALLOCATABLE, DIMENSION  (:,:) :: ttmp0_
      CONTAINS
        ! Public methods:
        PROCEDURE,PUBLIC :: io_write_pdbv     => VGPIC_io_write_pdbv
        PROCEDURE,PUBLIC :: io_readv          => VGPIC_io_read_pdbv
        PROCEDURE,PUBLIC :: GetFlux           => VGPIC_GetFlux
        PROCEDURE,PUBLIC :: GetMoment         => VGPIC_GetMoment
        PROCEDURE,PUBLIC :: InitUserSeed      => VGPIC_InitUserSeed
        PROCEDURE,PUBLIC :: ResizeArrays      => VGPIC_ResizeArrays
      END TYPE VGPIC
  
  TYPE, PUBLIC, EXTENDS ( VGPIC ) :: ChargPIC
        PRIVATE
        ! Member data:
        REAL(KIND=GP), ALLOCATABLE, DIMENSION    (:) :: lbx_,lby_,lbz_
        REAL(KIND=GP), ALLOCATABLE, DIMENSION    (:) :: lfx_,lfy_,lfz_
      CONTAINS
        PROCEDURE,PUBLIC :: ChargPIC_ctor
        FINAL            :: ChargPIC_dtor
        PROCEDURE,PUBLIC :: SetStepVel        => ChargPIC_SetStepRKK
        PROCEDURE,PUBLIC :: StepChargedPIC    => ChargPIC_StepRKK
        PROCEDURE,PUBLIC :: StepChargedPICBor => ChargPIC_StepBoris
        PROCEDURE,PUBLIC :: EndStageChargedPIC=> ChargPIC_EndStageRKK
        PROCEDURE,PUBLIC :: GetTemperature    => ChargPIC_GetTemperature
        PROCEDURE,PUBLIC :: GetTemperatureAnis=> ChargPIC_GetTemperatureAnis
        PROCEDURE,PUBLIC :: InitFromFields    => ChargPIC_InitFromFields
        PROCEDURE,PUBLIC :: ResizeArrays      => ChargPIC_ResizeArrays
      END TYPE ChargPIC
