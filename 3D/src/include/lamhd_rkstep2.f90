! Step 2 of Runge-Kutta for the Lagrangian averaged MHD equations

         CALL rotor3(ay,az,C7,1)
         CALL rotor3(ax,az,C8,2)
         CALL rotor3(ax,ay,C9,3)             ! b is now C7-9
         CALL aprodre3(vx,vy,vz,C10,C11,C12,alpk) ! omega X u in C10-2
         CALL aprodre3(C7,C8,C9,C13,C14,C15,alpm) ! j X b_s in C13-5
         IF ((trans.ge.1).and.(times.eq.0).and.(bench.eq.0).and.(o.eq.ord)) &
            CALL aentrans(vx,vy,vz,C13,C14,C15,alpk,ext,2)
         CALL nonlin3(C10,C11,C12,C13,C14,C15,C16,1)
         CALL nonlin3(C10,C11,C12,C13,C14,C15,C17,2)
         CALL nonlin3(C10,C11,C12,C13,C14,C15,C10,3) ! C16-7,10 include divP
         CALL avector3(vx,vy,vz,C7,C8,C9,C11,C12,C13,alpk,alpm) ! u X b_s C11-3
         CALL gauge3(C11,C12,C13,C7,1)
         CALL gauge3(C11,C12,C13,C8,2)
         CALL gauge3(C11,C12,C13,C9,3) ! u X b_s so div a_s=0 in C7-9
         CALL laplak3(vx,vx)
         CALL laplak3(vy,vy)
         CALL laplak3(vz,vz)   ! vx-z now Laplacian of vx-z
         CALL laplak3(ax,ax)
         CALL laplak3(ay,ay)
         CALL laplak3(az,az)   ! ax-z now Laplacian of ax-z
         IF ((trans.ge.1).and.(times.eq.0).and.(bench.eq.0).and.(o.eq.ord)) &
            THEN
            CALL aentrans(vx,vy,vz,C16,C17,C10,alpk,ext,1)
            CALL aentrans(ax,ay,az,C7,C8,C9,alpm,ext,0)
         ENDIF

         DO i = ista,iend
            DO j = 1,n
               DO k = 1,n

            IF ((ka2(k,j,i).le.kmax).and.(ka2(k,j,i).ge.tiny)) THEN
               vx(k,j,i) = C1(k,j,i)+dt*(nu*vx(k,j,i)+C16(k,j,i) &
              +fx(k,j,i))/real(o,kind=GP)
               vy(k,j,i) = C2(k,j,i)+dt*(nu*vy(k,j,i)+C17(k,j,i) &
              +fy(k,j,i))/real(o,kind=GP)
               vz(k,j,i) = C3(k,j,i)+dt*(nu*vz(k,j,i)+C10(k,j,i) &
              +fz(k,j,i))/real(o,kind=GP)
               ax(k,j,i) = C4(k,j,i)+dt*(mu*ax(k,j,i)+C7(k,j,i)  &
              +mx(k,j,i))*(1+alpm**2*ka2(k,j,i))/real(o,kind=GP)
               ay(k,j,i) = C5(k,j,i)+dt*(mu*ay(k,j,i)+C8(k,j,i)  &
              +my(k,j,i))*(1+alpm**2*ka2(k,j,i))/real(o,kind=GP)
               az(k,j,i) = C6(k,j,i)+dt*(mu*az(k,j,i)+C9(k,j,i)  &
              +mz(k,j,i))*(1+alpm**2*ka2(k,j,i))/real(o,kind=GP)
             ELSE
               vx(k,j,i) = 0.
               vy(k,j,i) = 0.
               vz(k,j,i) = 0.
               ax(k,j,i) = 0.
               ay(k,j,i) = 0.
               az(k,j,i) = 0.
            ENDIF

                END DO
            END DO
         END DO
