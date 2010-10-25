! Step 2 of Runge-Kutta for the HMHD equations
! Computes the nonlinear terms and evolves the equations in dt/o
         CALL laplak2(ps,C5)
         CALL laplak2(az,C6)
         CALL poisson(ps,az,C7)
         CALL poisson(ps,vz,C8)
         CALL poisson(ps,bz,C9)
         CALL poisson(az,vz,C10)
         CALL poisson(az,bz,C11)
         CALL poisson(ps,C5,ps)
         CALL poisson(az,C6,az)
         CALL laplak2(vz,vz)
         CALL laplak2(bz,bz)

         DO i = ista,iend
            DO j = 1,n

            IF ((ka2(j,i).le.kmax).and.(ka2(j,i).ge.tiny)) THEN
               ps(j,i) = C1(j,i)+dt*(nu*C5(j,i)+(az(j,i)-ps(j,i))    &
              /ka2(j,i)+fk(j,i))/float(o)
               az(j,i) = C2(j,i)+dt*(mu*C6(j,i)+C7(j,i)+ep*C11(j,i)  &
              +mk(j,i))/float(o)
               vz(j,i) = C3(j,i)+dt*(nu*vz(j,i)+C8(j,i)-C11(j,i)     &
              +fz(j,i))/float(o)
               bz(j,i) = C4(j,i)+dt*(mu*bz(j,i)+C9(j,i)-C10(j,i)     &
              -ep*az(j,i)+mz(j,i))/float(o)
            ELSE
               ps(j,i) = 0.0_GP
               az(j,i) = 0.0_GP
               vz(j,i) = 0.0_GP
               bz(j,i) = 0.0_GP
            ENDIF

            END DO
         END DO
