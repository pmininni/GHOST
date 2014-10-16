! Step 2 of Runge-Kutta for the MHD equations
! Computes the nonlinear terms and evolves the equations in dt/o

         CALL rotor3(ay,az,C7,1)
         CALL rotor3(ax,az,C8,2)
         CALL rotor3(ax,ay,C9,3)
         CALL prodre3(vx,vy,vz,C10,C11,C12)
         CALL prodre3(C7,C8,C9,C13,C14,C15)
         IF ((trans.eq.1).and.(times.eq.0).and.(bench.eq.0).and.(o.eq.ord)) &
            CALL entrans(C1,C2,C3,C13,C14,C15,ext,2)
         CALL nonlin3(C10,C11,C12,C13,C14,C15,C16,1)
         CALL nonlin3(C10,C11,C12,C13,C14,C15,C17,2)
         CALL nonlin3(C10,C11,C12,C13,C14,C15,C10,3)
         IF ((trans.eq.1).and.(times.eq.0).and.(bench.eq.0).and.(o.eq.ord)) &
            CALL entrans(C7,C8,C9,C16,C17,C10,ext,4)
         CALL vector3(vx,vy,vz,C7,C8,C9,C11,C12,C13)
         CALL gauge3(C11,C12,C13,C7,1)
         CALL gauge3(C11,C12,C13,C8,2)
         CALL gauge3(C11,C12,C13,C9,3)
         CALL laplak3(vx,vx)
         CALL laplak3(vy,vy)
         CALL laplak3(vz,vz)
         CALL laplak3(ax,ax)
         CALL laplak3(ay,ay)
         CALL laplak3(az,az)
         IF ((trans.eq.1).and.(times.eq.0).and.(bench.eq.0).and.(o.eq.ord)) &
            THEN
            CALL entrans(C1,C2,C3,C16,C17,C10,ext,1)
            CALL entrans(C4,C5,C6,C7,C8,C9,ext,0)
            CALL rotor3(C8,C9,C11,1)
            CALL rotor3(C7,C9,C12,2)
            CALL rotor3(C7,C8,C13,3)
            CALL entrans(C1,C2,C3,C11,C12,C13,ext,3)
         ENDIF

         rmp = 1./real(o,kind=GP)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,n
         DO k = 1,n
            IF ((ka2(k,j,i).le.kmax).and.(ka2(k,j,i).ge.tiny)) THEN
               vx(k,j,i) = C1(k,j,i)+dt*(nu*vx(k,j,i)+C16(k,j,i) &
              +fx(k,j,i))*rmp
               vy(k,j,i) = C2(k,j,i)+dt*(nu*vy(k,j,i)+C17(k,j,i) &
              +fy(k,j,i))*rmp
               vz(k,j,i) = C3(k,j,i)+dt*(nu*vz(k,j,i)+C10(k,j,i) &
              +fz(k,j,i))*rmp
               ax(k,j,i) = C4(k,j,i)+dt*(mu*ax(k,j,i)+C7(k,j,i)  &
              +mx(k,j,i))*rmp
               ay(k,j,i) = C5(k,j,i)+dt*(mu*ay(k,j,i)+C8(k,j,i)  &
              +my(k,j,i))*rmp
               az(k,j,i) = C6(k,j,i)+dt*(mu*az(k,j,i)+C9(k,j,i)  &
              +mz(k,j,i))*rmp
            ELSE
               vx(k,j,i) = 0.0_GP
               vy(k,j,i) = 0.0_GP
               vz(k,j,i) = 0.0_GP
               ax(k,j,i) = 0.0_GP
               ay(k,j,i) = 0.0_GP
               az(k,j,i) = 0.0_GP
            ENDIF
         END DO
         END DO
         END DO
