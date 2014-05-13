! Step 2 of Runge-Kutta for the Lagrangian averaged HD equations

         CALL aprodre3(vx,vy,vz,C4,C5,C6,alpk)
         CALL nonlhd3(C4,C5,C6,C7,1)
         CALL nonlhd3(C4,C5,C6,C8,2)
         CALL nonlhd3(C4,C5,C6,C4,3)
         CALL laplak3(vx,vx)
         CALL laplak3(vy,vy)
         CALL laplak3(vz,vz)
         IF ((trans.eq.1).and.(times.eq.sstep).and.(bench.eq.0).and.(o.eq.ord)) &
            CALL entrans(C1,C2,C3,C7,C8,C4,ext,1)

         DO i = ista,iend
            DO j = 1,n
               DO k = 1,n

            IF ((ka2(k,j,i).le.kmax).and.(ka2(k,j,i).ge.tiny)) THEN
               vx(k,j,i) = C1(k,j,i)+dt*(nu*vx(k,j,i)+C7(k,j,i) &
              +fx(k,j,i))/real(o,kind=GP)
               vy(k,j,i) = C2(k,j,i)+dt*(nu*vy(k,j,i)+C8(k,j,i) &
              +fy(k,j,i))/real(o,kind=GP)
               vz(k,j,i) = C3(k,j,i)+dt*(nu*vz(k,j,i)+C4(k,j,i) &
              +fz(k,j,i))/real(o,kind=GP)
            ELSE
               vx(k,j,i) = 0.
               vy(k,j,i) = 0.
               vz(k,j,i) = 0.
            ENDIF

               END DO
            END DO
         END DO
