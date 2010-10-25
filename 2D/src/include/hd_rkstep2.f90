! Step 2 of Runge-Kutta for the HD equations
! Computes the nonlinear terms and evolves the equations in dt/o
         CALL laplak2(ps,C2)
         CALL poisson(ps,C2,ps)

         rmp = 1.0_GP/real(o,kind=GP)
         DO i = ista,iend
         DO j = 1,n
            IF ((ka2(j,i).le.kmax).and.(ka2(j,i).ge.tiny)) THEN
               ps(j,i) = C1(j,i)+dt*(nu*C2(j,i)-ps(j,i)/ka2(j,i)   &
              +fk(j,i))*rmp
            ELSE
               ps(j,i) = 0.0_GP
            ENDIF
         END DO
         END DO
