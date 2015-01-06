! Step 2 of Runge-Kutta for the SWHD equations
! Computes the nonlinear terms and evolves the equations in dt/o
         DO i = ista,iend
            DO j = 1,n
               C13(j,i) = th(j,i)-fs(j,i)
            END DO
         END DO
         CALL buhler(C13,vx,vy,C17,C18)
         CALL elemwise2(vx,C13,C14)
         CALL elemwise2(vy,C13,C13)
         CALL gradre2(vx,vy,C15,C16)
         CALL laplak2(vx,vx)
         CALL laplak2(vy,vy)

        IF ((trans.eq.1).and.(times.eq.0).and.(bench.eq.0).and.(o.eq.ord)) &
           CALL v_entrans(C1,C2,C15,C16,ext)

         rmp = 1.0_GP/real(o,kind=GP)
         DO i = ista,iend
         DO j = 1,n
            IF ((ka2(j,i).le.kmax).and.(ka2(j,i).ge.tiny)) THEN
               tmp     = 1/(1+switch*(c0**2)*ka2(j,i)/3)
               vx(j,i) = C1(j,i) + dt*(nu*vx(j,i) - C15(j,i) &
                        - g*im*ka(i)*th(j,i) + fx(j,i)       &
                        + nu*C17(j,i))*rmp*tmp
               vy(j,i) = C2(j,i) + dt*(nu*vy(j,i) - C16(j,i) &
                        - g*im*ka(j)*th(j,i) + fy(j,i)       &
                        + nu*C18(j,i))*rmp*tmp
               th(j,i) = C12(j,i) - dt*im*(ka(i)*C14(j,i)   &
                        + ka(j)*C13(j,i))*rmp
            ELSE IF (ka2(j,i).gt.kmax) THEN
               vx(j,i) = 0.0_GP
               vy(j,i) = 0.0_GP
               th(j,i) = 0.0_GP
            ELSE IF (ka2(j,i).lt.tiny) THEN
               vx(j,i) = 0.0_GP
               vy(j,i) = 0.0_GP
               th(j,i) = C12(j,i)               
            ENDIF
         END DO
         END DO
         
