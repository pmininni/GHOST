! Step 2 of Runge-Kutta for the EDQNMROTH equations (with eddy viscosity)
!
! Compute eddy viscosities, transfer functions from partially-updated v_i:
         CALL spectrumc(vx,vy,vz,1,heli,Eold,Hold)
         CALL evedqnm(Eold,Hold,nu,(t-1)*dt,heli,tepq,thpq,tve,tvh,Eext,Hext)
         rmp = 1./real(n,KIND=GP)**6
         DO i = ista,iend 
         DO j = 1,n
         DO k = 1,n
            Eden(k,j,i) = (abs(vx(k,j,i))**2+abs(vy(k,j,i))**2+           &
                          abs(vz(k,j,i))**2)*rmp
            Hden(k,j,i) = 2*rmp*(ka(k)*(real(vx(k,j,i))*aimag(vy(k,j,i))- &
                          real(vy(k,j,i))*aimag(vx(k,j,i)))               &
                          +ka(j)*(real(vz(k,j,i))*aimag(vx(k,j,i))-       &
                          real(vx(k,j,i))*aimag(vz(k,j,i)))               &
                          +ka(i)*(real(vy(k,j,i))*aimag(vz(k,j,i))-       &
                          real(vz(k,j,i))*aimag(vy(k,j,i))))
         END DO
         END DO
         END DO
!
! Compute the nonlinear terms and evolve the equations in dt/o
         CALL prodre3(vx,vy,vz,C4,C5,C6)
         DO i = ista,iend               ! Coriolis force
            DO j = 1,n
               DO k = 1,n
                  C4(k,j,i) = C4(k,j,i)+2*omega*vy(k,j,i)
                  C5(k,j,i) = C5(k,j,i)-2*omega*vx(k,j,i)
               END DO
            END DO
         END DO
         CALL nonlhd3(C4,C5,C6,C7,1)
         CALL nonlhd3(C4,C5,C6,C8,2)
         CALL nonlhd3(C4,C5,C6,C4,3)
         CALL rotor3 (vy,vz,C5 ,1) ! w_1=curl(v)_1
         CALL rotor3 (vx,vz,C6 ,2) ! w_2=curl(v)_2
         CALL rotor3 (vx,vy,C19,3) ! w_3=curl(v)_3
         CALL laplak3(vx,vx)
         CALL laplak3(vy,vy)
         CALL laplak3(vz,vz)
         IF ((trans.eq.1).and.(times.eq.sstep).and.(bench.eq.0).and.(o.eq.ord)) &
            THEN
            CALL entrans(C1,C2,C3,C7,C8,C4,ext,1)
            CALL entpara(C1,C2,C3,C7,C8,C4,ext,1)
            CALL entperp(C1,C2,C3,C7,C8,C4,ext,1)
            CALL heltrans(C1,C2,C3,C7,C8,C4,ext,1)
            CALL heltpara(C1,C2,C3,C7,C8,C4,ext,1)
            CALL heltperp(C1,C2,C3,C7,C8,C4,ext,1)
         ENDIF
!
! Advance first using only kinematic and energy and helicity eddy viscosities
         rmp = 1./real(o,kind=GP)
         DO i = ista,iend 
         DO j = 1,n
         DO k = 1,n
            IF ((ka2(k,j,i).le.kmax).and.(ka2(k,j,i).ge.tiny)) THEN
               ki = int(sqrt(ka2(k,j,i)) + 0.501)
               vx(k,j,i) = C1(k,j,i)+dt*((nu+tve(ki))*vx(k,j,i)+C7(k,j,i) &
              +tvh(ki)*C5(k,j,i)+fx(k,j,i))*rmp
               vy(k,j,i) = C2(k,j,i)+dt*((nu+tve(ki))*vy(k,j,i)+C8(k,j,i) &
              +tvh(ki)*C6(k,j,i)+fy(k,j,i))*rmp
               vz(k,j,i) = C3(k,j,i)+dt*((nu+tve(ki))*vz(k,j,i)+C4(k,j,i) &
              +tvh(ki)*C19(k,j,i)+fz(k,j,i))*rmp
            ELSE
               vx(k,j,i) = 0.
               vy(k,j,i) = 0.
               vz(k,j,i) = 0.
            ENDIF
         END DO
         END DO
         END DO
!
! Finally, must reconstruct the velocity components based on corrected
! energy and helicity spectra:
         CALL vcorrect(vx,vy,vz,Eden,Hden,tepq,thpq,Eold,Hold,dt*rmp,heli)
