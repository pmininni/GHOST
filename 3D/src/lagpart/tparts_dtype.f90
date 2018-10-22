!=================================================================
! GHOST TestGPart test particles subclass
!
! 2015 P. Dmitruk
!      Department of Physics, 
!      Facultad de Ciencias Exactas y Naturales.
!      Universidad de Buenos Aires.
!=================================================================

  TYPE, PUBLIC, EXTENDS ( VGPart ) :: TestGPart
        PRIVATE
        ! Member data:
        REAL(KIND=GP), ALLOCATABLE, DIMENSION    (:) :: lbx_,lby_,lbz_
        REAL(KIND=GP), ALLOCATABLE, DIMENSION    (:) :: lfx_,lfy_,lfz_
      CONTAINS
        ! Public methods:
        PROCEDURE,PUBLIC :: TestGPart_ctor
        FINAL            :: TestGPart_dtor
        PROCEDURE,PUBLIC :: InitVel           => TestGPart_InitVel
        PROCEDURE,PUBLIC :: SetStepVel        => TestGPart_SetStepRKK
        PROCEDURE,PUBLIC :: StepTestp         => TestGPart_StepRKK
        PROCEDURE,PUBLIC :: EndStage          => TestGPart_EndStageRKK
      END TYPE TestGPart
