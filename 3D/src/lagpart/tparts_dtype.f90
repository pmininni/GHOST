!=================================================================
! GHOST TestGPart test particles subclass
!
! 2015 P. Dmitruk
!      Department of Physics, 
!      Facultad de Ciencias Exactas y Naturales.
!      Universidad de Buenos Aires.
!=================================================================

  TYPE, PUBLIC, EXTENDS ( GPart ) :: TestGPart
        PRIVATE
        ! Member data:
        REAL(KIND=GP), ALLOCATABLE, DIMENSION    (:) :: pvx_,pvy_,pvz_
        REAL(KIND=GP), ALLOCATABLE, DIMENSION    (:) :: lbx_,lby_,lbz_
        REAL(KIND=GP), ALLOCATABLE, DIMENSION    (:) :: lfx_,lfy_,lfz_
        REAL(KIND=GP), ALLOCATABLE, DIMENSION  (:,:) :: ttmp0_
      CONTAINS
        ! Public methods:
        PROCEDURE,PUBLIC :: TestGPart_ctor
        FINAL            :: TestGPart_dtor
        PROCEDURE,PUBLIC :: InitVel           => TestGPart_InitVel
        PROCEDURE,PUBLIC :: SetStepVel        => TestGPart_SetStepRKK
        PROCEDURE,PUBLIC :: StepTestp         => TestGPart_StepRKK
        PROCEDURE,PUBLIC :: EndStage          => TestGPart_EndStageRKK
        PROCEDURE,PUBLIC :: io_write_pdbv     => TestGPart_io_write_pdbv
      END TYPE TestGPart
