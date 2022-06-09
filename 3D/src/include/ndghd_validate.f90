! This file verifies that all parameters and compilation options
! are compatible with the solver used. If empty, the solver is
! compatible with all possible configurations of GHOST.

      IF ( stat.eq.0 ) THEN
        IF (myrank.eq.0) THEN
           PRINT *,'MAIN: Nudging can only be used with stat >= 0'
           PRINT *,'The solver must read fxold and fxnew when it'
           PRINT *,'starts to initialize the nudging fields'
        ENDIF
        STOP
      ENDIF

      IF ( rand.ne.2 ) THEN
        IF (myrank.eq.0) THEN
           PRINT *,'MAIN: Nudging can only be used with rand=2'
           PRINT *,'The solver reads the nudging fields using'
           PRINT *,'initialfv.f90_nudging'
        ENDIF
        STOP
      ENDIF
