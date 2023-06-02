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
        TYPE(GFieldComm)                             :: gfcomm_
        TYPE(GPICSplineInt)                          :: picspl_
        REAL(KIND=GP), ALLOCATABLE, DIMENSION    (:) :: prop_
        REAL(KIND=GP)                                :: icv_
      CONTAINS
        ! Public methods:
        PROCEDURE,PUBLIC :: GPIC_ctor
        PROCEDURE,PUBLIC :: Step              => GPIC_StepRKK
        PROCEDURE,PUBLIC :: GetDensity        => GPIC_GetDensity
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
      END TYPE VGPIC
  
  TYPE, PUBLIC, EXTENDS ( VGPIC ) :: ChargPIC
        PRIVATE
        ! Member data:
        REAL(KIND=GP), ALLOCATABLE, DIMENSION    (:) :: lbx_,lby_,lbz_
        REAL(KIND=GP), ALLOCATABLE, DIMENSION    (:) :: lfx_,lfy_,lfz_
      CONTAINS
        PROCEDURE,PUBLIC :: ChargPIC_ctor
        FINAL            :: ChargPIC_dtor
        PROCEDURE,PUBLIC :: InitVel           => ChargPIC_InitVel
        PROCEDURE,PUBLIC :: SetStepVel        => ChargPIC_SetStepRKK
        PROCEDURE,PUBLIC :: StepChargedPIC    => ChargPIC_StepRKK
        PROCEDURE,PUBLIC :: EndStage          => ChargPIC_EndStageRKK
        PROCEDURE,PUBLIC :: GetTemperature    => ChargPIC_GetTemperature
      END TYPE ChargPIC
