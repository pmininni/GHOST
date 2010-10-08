! Initial condition for the velocity.
! This file contains the expression used for the initial
! streamfunction or velocity field (depending on the solver). 
! You can use the temporary real array R1 of size 
! (n,jsta:jend), and temporary complex arrays C1, C2 of size 
! (n,ista:iend) to do intermediate computations. The variable 
! u0 should control the global amplitude of the field, and 
! variables vparam0-9 can be used to control the amplitudes 
! of individual terms. At the end, the streamfunction should 
! be stored in the array ps, or the velocity field components 
! in the arrays vx and vy (plus vz in 2.5D solvers).

! Superposition of harmonic modes with random phases
! (streamfunction, 2D)
!     kdn : minimum wave number
!     kup : maximum wave number

      IF (ista.eq.1) THEN
         ps(1,1) = 0.0_GP
         DO j = 2,n/2+1
            IF ((ka2(j,1).le.kup**2).and.(ka2(j,1).ge.kdn**2)) THEN
               phase = 2*pi*randu(seed)
               ps(j,1) = (COS(phase)+im*SIN(phase))/sqrt(ka2(j,1))
               ps(n-j+2,1) = conjg(ps(j,1))
            ELSE
               ps(j,1) = 0.0_GP
               ps(n-j+2,1) = 0.0_GP
            ENDIF
         END DO
         DO j = 1,n
            DO i = 2,iend
               IF ((ka2(j,i).le.kup**2).and.(ka2(j,i).ge.kdn**2)) THEN
                  phase = 2*pi*randu(seed)
                  ps(j,i) = 2*(COS(phase)+im*SIN(phase))/sqrt(ka2(j,i))
               ELSE
                  ps(j,i) = 0.0_GP
               ENDIF
            END DO
         END DO
      ELSE
         DO j = 1,n
            DO i = ista,iend
               IF ((ka2(j,i).le.kup**2).and.(ka2(j,i).ge.kdn**2)) THEN
                  phase = 2*pi*randu(seed)
                  ps(j,i) = 2*(COS(phase)+im*SIN(phase))/sqrt(ka2(j,i))
               ELSE
                  ps(j,i) = 0.0_GP
               ENDIF
            END DO
         END DO
      ENDIF
      CALL normalize(ps,u0,1,MPI_COMM_WORLD)
