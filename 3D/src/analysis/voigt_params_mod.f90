!=================================================================

  MODULE voigt_params
!
      USE fprecision
 
      TYPE GVoigtParamType
         INTEGER           :: doSGSinj
         INTEGER           :: prtbin
         INTEGER           :: icycle ! continuous time icycle
         REAL(KIND=GP)     :: nu
         REAL(KIND=GP)     :: kappa
         REAL(KIND=GP)     :: bvfreq
         REAL(KIND=GP)     :: rotf
         REAL(KIND=GP)     :: ttime  ! time stamp
         REAL(KIND=GP)     :: dt
         CHARACTER(len=64) :: ext
      END TYPE GVoigtParamType
      SAVE

  END MODULE voigt_params
