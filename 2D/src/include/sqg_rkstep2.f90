! Step 2 of Runge-Kutta for the SQG equations
! Computes the nonlinear terms and evolves the equations in dt/o
         CALL scalarder(ps,C2)
         CALL poisson(ps,C2,ps)

         rmp = 1.0_GP/real(o,kind=GP)
         DO i = ista,iend
         DO j = 1,n
            IF ((ka2(j,i).le.kmax).and.(ka2(j,i).ge.tiny)) THEN
               ps(j,i) = C1(j,i)+dt*(sqrt(ka2(j,i))*nu*C2(j,i)+ps(j,i) &
              /sqrt(ka2(j,i))+fk(j,i))*rmp
            ELSE
               ps(j,i) = 0.0_GP
            ENDIF
         END DO
         END DO
