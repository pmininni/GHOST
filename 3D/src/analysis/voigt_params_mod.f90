!=================================================================

  MODULE voigt_params
!
      USE fprecision
 
      TYPE GVoigtParamType
         REAL(KIND=GP)    :: nu
         REAL(KIND=GP)    :: kappa
         REAL(KIND=GP)    :: bvfreq
         REAL(KIND=GP)    :: rotf
      END TYPE GVoigtParamType
      SAVE

  END MODULE voigt_params
