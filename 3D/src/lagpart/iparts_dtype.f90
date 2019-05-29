!=================================================================
! GHOST InerGPart test particles subclass
! Also creates the VGPart subclass, for all particles that
! need a velocity different from the Lagrangian velocity. This
! subclass must be declared before any other subclass that uses
! particle velocities pvx_,pvy_,pvz_.
!
! 2018 F. Falkinhoff and P.D. Mininni
!      Department of Physics, 
!      Facultad de Ciencias Exactas y Naturales.
!      Universidad de Buenos Aires.
!=================================================================

  TYPE, PUBLIC, EXTENDS ( GPart ) :: VGPart
        PRIVATE
        ! Member data:
        REAL(KIND=GP), ALLOCATABLE, DIMENSION    (:) :: pvx_,pvy_,pvz_
        REAL(KIND=GP), ALLOCATABLE, DIMENSION  (:,:) :: ttmp0_
      CONTAINS
        ! Public methods:
        PROCEDURE,PUBLIC :: io_write_pdbv     => VGPart_io_write_pdbv
        PROCEDURE,PUBLIC :: io_readv          => VGPart_io_read_pdbv
      END TYPE VGPart

  TYPE, PUBLIC, EXTENDS ( VGPart ) :: InerGPart
        PRIVATE
        ! Member data:
        REAL(KIND=GP), ALLOCATABLE, DIMENSION    (:) :: dfx_,dfy_,dfz_
        REAL(KIND=GP)    :: invtau_,grav_,gamma_
      CONTAINS
        ! Public methods:
        PROCEDURE,PUBLIC :: InerGPart_ctor
        FINAL            :: InerGPart_dtor
        PROCEDURE,PUBLIC :: InitVel           => InerGPart_InitVel
        PROCEDURE,PUBLIC :: SetStepVel        => InerGPart_SetStepRKK
        PROCEDURE,PUBLIC :: StepInerp         => InerGPart_StepRKK
        PROCEDURE,PUBLIC :: StepLitep         => InerGPart_lite_StepRKK
        PROCEDURE,PUBLIC :: EndStage          => InerGPart_EndStageRKK
      END TYPE InerGPart
