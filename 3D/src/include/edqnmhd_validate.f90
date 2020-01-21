! This file verifies that all parameters and compilation options
! are compatible with the solver used. If empty, the solver is
! compatible with all possible configurations of GHOST.

      IF ( (anis.eq.1).or.(nx.ne.ny).or.(nx.ne.nz).or.(ny.ne.nz) ) THEN
        IF (myrank.eq.0) THEN
           PRINT *,'EDQNM solvers require (2.pi)^3 cubic boxes.'
           PRINT *,'Please compile with ARBSIZE = no, and with nx = ny = nz.'
        ENDIF
        STOP
      ENDIF
     
