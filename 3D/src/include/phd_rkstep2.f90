! Step 2 of Runge-Kutta for the HD equations
! Computes the nonlinear terms and evolves the equations in dt/o

         CALL prodre3(vx,vy,vz,C4,C5,C6)
         CALL nonlhd3(C4,C5,C6,C7,1)
         CALL nonlhd3(C4,C5,C6,C8,2)
         CALL nonlhd3(C4,C5,C6,C4,3)
         CALL advect3(vx,vy,vz,th,C5)
         CALL laplak3(vx,vx)
         CALL laplak3(vy,vy)
         CALL laplak3(vz,vz)
         CALL laplak3(th,th)
         IF ((trans.eq.1).and.(times.eq.0).and.(bench.eq.0).and.(o.eq.ord)) &
            THEN
            CALL entrans(C1,C2,C3,C7,C8,C4,ext,1)
            CALL sctrans(C20,C5,ext)
         ENDIF

         rmp = 1./float(o)
         DO i = ista,iend 
         DO j = 1,n
         DO k = 1,n
            IF ((ka2(k,j,i).le.kmax).and.(ka2(k,j,i).ge.tiny)) THEN
               vx(k,j,i) = C1(k,j,i)+dt*(nu*vx(k,j,i)+C7(k,j,i) &
              +fx(k,j,i))*rmp
               vy(k,j,i) = C2(k,j,i)+dt*(nu*vy(k,j,i)+C8(k,j,i) &
              +fy(k,j,i))*rmp
               vz(k,j,i) = C3(k,j,i)+dt*(nu*vz(k,j,i)+C4(k,j,i) &
              +fz(k,j,i))*rmp
               th(k,j,i) = C20(k,j,i)+dt*(kappa*th(k,j,i)+C5(k,j,i) &
              +fs(k,j,i))*rmp
            ELSE IF (ka2(k,j,i).gt.kmax) THEN
               vx(k,j,i) = 0.
               vy(k,j,i) = 0.
               vz(k,j,i) = 0.
               th(k,j,i) = 0.
            ELSE IF (ka2(k,j,i).lt.tiny) THEN
               vx(k,j,i) = 0.
               vy(k,j,i) = 0.
               vz(k,j,i) = 0.
            ENDIF
         END DO
         END DO
         END DO
