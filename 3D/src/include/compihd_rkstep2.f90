! Step 2 of Runge-Kutta for the compressible HD equations
! Computes the nonlinear terms and evolves the equations in dt/o

         CALL gradre3(vx,vy,vz,C4,C5,C6)          ! v.Grad v
         CALL gradpressi(gam1,th,C31,C32,C33)     ! Grad p term

         CALL divide(rho,C31,C32,C33)             ! divide Grad p by rho
         CALL vdiss(nu,nu2,vx,vy,vz)              ! viscous term
         CALL divide(rho,vx,vy,vz)                ! divide viscous term by rho

         CALL divrhov(rho,vx,vy,vz,C7)            ! div(rho.v)

         CALL divrhov(th,vx,vy,vz,C8)             ! div(e.v)
         CALL pdVwork(gam1,th,vx,vy,vz,C34)       ! p.div(v) = (gamma-1) e.div(v)
         CALL viscHeatRayleigh(vx,vy,vz,C36)      ! phi/mu, visc. heat
         CALL laplak3(th,th)                      ! laplacian(e)
        
         rmp = 1./real(o,kind=GP)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
         DO k = 1,nz

            IF (kn2(k,j,i).le.kmax) THEN
               vx (k,j,i) = C1(k,j,i)+dt*(vx(k,j,i)-C4(k,j,i)-C31(k,j,i) &
              +fx(k,j,i))*rmp
               vy (k,j,i) = C2(k,j,i)+dt*(vy(k,j,i)-C5(k,j,i)-C32(k,j,i) &
              +fy(k,j,i))*rmp
               vz (k,j,i) = C3(k,j,i)+dt*(vz(k,j,i)-C6(k,j,i)-C33(k,j,i) &
              +fz(k,j,i))*rmp
               rho(k,j,i) = C20(k,j,i)-dt*C7(k,j,i)*rmp
               th (k,j,i) = C35(k,j,i)+dt*(mu*C36(k,j,i)-C8(k,j,i)-C34(k,j,i) &
              +fs(k,j,i))*rmp
            ELSE IF (kn2(k,j,i).gt.kmax) THEN
               vx (k,j,i) = 0.0_GP
               vy (k,j,i) = 0.0_GP
               vz (k,j,i) = 0.0_GP
               rho(k,j,i) = 0.0_GP
               th (k,j,i) = 0.0_GP
            ELSE IF (kn2(k,j,i).lt.tiny) THEN
               vx (k,j,i) = 0.0_GP
               vy (k,j,i) = 0.0_GP
               vz (k,j,i) = 0.0_GP
               rho(k,j,i) = C20(k,j,i)
               th (k,j,i) = C35(k,j,i)
            ENDIF

         END DO
         END DO
         END DO
