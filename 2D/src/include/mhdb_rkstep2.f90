! Step 2 of Runge-Kutta for the MHD equations
! Computes the nonlinear terms and evolves the equations in dt/o
         CALL laplak2(ps,C3)
         CALL laplak2(az,C4)
         CALL poissonb0(ps,az,C5,by0)
         IF ((trans.eq.1).and.(times.eq.0).and.(bench.eq.0).and.(o.eq.ord)) &
            CALL vectrans(az,C5,ext)
         CALL poissonb0(C4,az,az,by0)
         CALL poisson(ps,C3,ps)

         rmp = 1.0_GP/real(o,kind=GP)
         DO i = ista,iend
            DO j = 1,n
            
            IF ((ka2(j,i).le.kmax).and.(ka2(j,i).ge.tiny)) THEN
               ps(j,i) = C1(j,i)+dt*(nu*C3(j,i)+(az(j,i)-ps(j,i))    &
              /ka2(j,i)+fk(j,i))*rmp
               az(j,i) = C2(j,i)+dt*(mu*C4(j,i)+C5(j,i)+mk(j,i))*rmp
            ELSE
               ps(j,i) = 0.0_GP
               az(j,i) = 0.0_GP
            ENDIF
            
            END DO
         END DO

