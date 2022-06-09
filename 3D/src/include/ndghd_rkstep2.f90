! Step 2 of Runge-Kutta for the HD equations
! Computes the nonlinear terms and evolves the equations in dt/o
         CALL prodre3(vx,vy,vz,C4,C5,C6)
         CALL nonlhd3(C4,C5,C6,C7,1)
         CALL nonlhd3(C4,C5,C6,C8,2)
         CALL nonlhd3(C4,C5,C6,C4,3)
         CALL laplak3(vx,vx)
         CALL laplak3(vy,vy)
         CALL laplak3(vz,vz)
         IF ((trans.eq.1).and.(times.eq.0).and.(bench.eq.0).and.(o.eq.ord)) &
            CALL entrans(C1,C2,C3,C7,C8,C4,ext,1)
         rmp = 1./real(o,kind=GP)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
         DO k = 1,nz
            IF ((kn2(k,j,i).le.kmax).and.(kn2(k,j,i).ge.tiny)) THEN
               ! Apply nudging only in forcing region
               IF ((kk2(k,j,i).le.kup**2).and.(kk2(k,j,k).ge.kdn**2)) THEN
                  vx(k,j,i) = C1(k,j,i)+dt*(nu*vx(k,j,i)+C7(k,j,i) &
                 +ndgamp*(vx(k,j,i)-fx(k,j,i)))*rmp
                  vy(k,j,i) = C2(k,j,i)+dt*(nu*vy(k,j,i)+C8(k,j,i) &
                 +ndgamp*(vx(k,j,i)-fy(k,j,i)))*rmp
                  vz(k,j,i) = C3(k,j,i)+dt*(nu*vz(k,j,i)+C4(k,j,i) &
                 +ndgamp*(vx(k,j,i)-fz(k,j,i)))*rmp
               ELSE
                  vx(k,j,i) = C1(k,j,i)+dt*(nu*vx(k,j,i)+C7(k,j,i))*rmp
                  vy(k,j,i) = C2(k,j,i)+dt*(nu*vy(k,j,i)+C8(k,j,i))*rmp
                  vz(k,j,i) = C3(k,j,i)+dt*(nu*vz(k,j,i)+C4(k,j,i))*rmp
               ENDIF
            ELSE
               vx(k,j,i) = 0.0_GP
               vy(k,j,i) = 0.0_GP
               vz(k,j,i) = 0.0_GP
            ENDIF
         END DO
         END DO
         END DO
