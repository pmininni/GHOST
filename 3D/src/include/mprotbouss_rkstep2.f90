! Step 2 of Runge-Kutta for the Boussinesq equations 
! Computes the nonlinear terms and evolves the equations in dt/o

         CALL prodre3(vx,vy,vz,C4,C5,C6)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend               ! Coriolis force
!$omp parallel do if (iend-ista.lt.nth) private (k)
            DO j = 1,ny
               DO k = 1,nz
                  C4(k,j,i) = C4(k,j,i)-2*omega*vy(k,j,i)
                  C5(k,j,i) = C5(k,j,i)+2*omega*vx(k,j,i)
                  C6(k,j,i) = C6(k,j,i)+xmom*th(k,j,i)
               END DO
            END DO
         END DO
         CALL nonlhd3(C4,C5,C6,C7,1)
         CALL nonlhd3(C4,C5,C6,C8,2)
         CALL nonlhd3(C4,C5,C6,C4,3)
         CALL advect3(vx,vy,vz,th,C5)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend               ! heat 'currrent':
!$omp parallel do if (iend-ista.lt.nth) private (k)
            DO j = 1,ny
               DO k = 1,nz
                  C5(k,j,i) = C5(k,j,i)+xtemp*vz(k,j,i)
               END DO
            END DO
         END DO
         CALL advect3(vx,vy,vz,th1,C6)
         CALL advect3(vx,vy,vz,th2,C24)
         CALL laplak3(vx,vx)
         CALL laplak3(vy,vy)
         CALL laplak3(vz,vz)
         CALL laplak3(th,th)
         CALL laplak3(th1,th1)
         CALL laplak3(th2,th2)

         IF ((trans.eq.1).and.(times.eq.0).and.(bench.eq.0).and.(o.eq.ord)) &
            THEN
            CALL entrans (C1,C2,C3,C7,C8,C4,ext,1)
            CALL entpara (C1,C2,C3,C7,C8,C4,ext,1)
            CALL entperp (C1,C2,C3,C7,C8,C4,ext,1)
            CALL heltrans(C1,C2,C3,C7,C8,C4,ext,1)
            CALL heltpara(C1,C2,C3,C7,C8,C4,ext,1)
            CALL heltperp(C1,C2,C3,C7,C8,C4,ext,1)
            CALL sctrans (C20,C5,ext,0)
            CALL sctpara (C20,C5,ext,0)
            CALL sctperp (C20,C5,ext,0)
            CALL sctrans (C21,C6,ext,1)
            CALL sctpara (C21,C6,ext,1)
            CALL sctperp (C21,C6,ext,1)
            CALL sctrans (C22,C24,ext,2)
            CALL sctpara (C22,C24,ext,2)
            CALL sctperp (C22,C24,ext,2)
         ENDIF

         rmp = 1.0_GP/(real(o,kind=GP))
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
         DO k = 1,nz
            IF ((kn2(k,j,i).le.kmax).and.(kn2(k,j,i).ge.tiny)) THEN
               vx(k,j,i) = C1 (k,j,i)+dt*(nu   *vx(k,j,i)+C7(k,j,i) &
              +fx(k,j,i))*rmp
               vy(k,j,i) = C2 (k,j,i)+dt*(nu   *vy(k,j,i)+C8(k,j,i) &
              +fy(k,j,i))*rmp
               vz(k,j,i) = C3 (k,j,i)+dt*(nu   *vz(k,j,i)+C4(k,j,i) &
              +fz(k,j,i))*rmp
               th(k,j,i) = C20(k,j,i)+dt*(kappa*th(k,j,i)+C5(k,j,i) &
              +fs(k,j,i))*rmp
               th1(k,j,i) = C21(k,j,i)+dt*(kappa1*th1(k,j,i)+C6(k,j,i) &
              +fs1(k,j,i))*rmp
               th2(k,j,i) = C22(k,j,i)+dt*(kappa2*th2(k,j,i)+C24(k,j,i) &
              +fs2(k,j,i))*rmp
            ELSE IF (kn2(k,j,i).gt.kmax) THEN
               vx (k,j,i) = 0.0_GP
               vy (k,j,i) = 0.0_GP
               vz (k,j,i) = 0.0_GP
               th (k,j,i) = 0.0_GP
               th1(k,j,i) = 0.0_GP
               th2(k,j,i) = 0.0_GP
            ELSE IF (kn2(k,j,i).lt.tiny) THEN
               vx (k,j,i) = 0.0_GP
               vy (k,j,i) = 0.0_GP
               vz (k,j,i) = 0.0_GP
               th (k,j,i) = C20(k,j,i)
               th1(k,j,i) = C21(k,j,i)
               th2(k,j,i) = C22(k,j,i)
            ENDIF
         END DO
         END DO
         END DO
