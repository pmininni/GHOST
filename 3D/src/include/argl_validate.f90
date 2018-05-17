! This file verifies that all parameters and compilation options
! are compatible with the solver used. If empty, the solver is
! compatible with all possible configurations of GHOST.

      IF ( ord.ne.1 ) THEN ! Check if the order is correct for ARGL
         WRITE(*,*)'MAIN: ARGL solver must be compiled with ord=1'
         STOP
      ENDIF
