! This file verifies that all parameters and compilation options
! are compatible with the solver used. If empty, the solver is
! compatible with all possible configurations of GHOST.

      PRINT *,'When using particles, penalty will only compute'
      PRINT *,'collisions between the obstacle and the particles'
      PRINT *,'if heavy particles are used.'
      
      IF (shape.ne.1) THEN
        IF (myrank.eq.0) THEN
           PRINT *,'Penalty for the moment only supports spherical'
           PRINT *,'objects. Please set shape = 1 or remove it from'
           PRINT *,'the input file.'
        ENDIF
        STOP
      ENDIF

      IF ((u0.lt.tiny).or.(vparam8.lt.tiny)) THEN
        IF (myrank.eq.0) THEN
           PRINT *,'Penalty requires setting a mean velocity (u0),'
           PRINT *,'and a rate at which the flow mean velocity is'
           PRINT *,'initially increased (vparam0). Please set these'
           PRINT *,'values in the input file.'
        ENDIF
        STOP
      ENDIF

