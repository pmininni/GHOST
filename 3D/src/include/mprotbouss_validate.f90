! This file verifies that all parameters and compilation options
! are compatible with the solver used. If empty, the solver is
! compatible with all possible configurations of GHOST.

      IF ( (abs(omegax).gt.tiny).or.(abs(omegay).gt.tiny ) ) THEN
        IF (myrank.eq.0) &
           PRINT *,'This solver has rotation only in the z-direction.'
           PRINT *,'Please set omegax and omegay = 0, or remove them'
           PRINT *,'from the input file.'
        STOP
      ENDIF

