! Step 2 of Runge-Kutta for the compressible HD equations
! Computes the nonlinear terms and evolves the equations in dt/o

         CALL prodre3(vx,vy,vz,C4,C5,C6)          ! om x v
         CALL gradpressi(gam1,th,vx,vy,vz,C21,C22,C23) ! Grad (p+0.5v^2) term
         CALL vdiss(nu,nu2,vx,vy,vz)              ! viscous term
         CALL divide(rho,vx,vy,vz)                ! viscous term with rho
         CALL divide(rho,C21,C22,C23)             ! divide Grad(p_tot) by rho

         CALL divrhov(rho,vx,vy,vz,C7)            ! div(rho.v)

         CALL divrhov(th,vx,vy,vz,C8)             ! div(e.v)
         CALL pdVwork(gam1,th,vx,vy,vz,C24)       ! p.div(v) = (gamma-1) e.div(v)
         CALL laplak3(th,th)                      ! laplacian(e)
        
         rmp = 1./real(o,kind=GP)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
         DO k = 1,nz

            IF (kn2(k,j,i).le.kmax) THEN
               vx (k,j,i) = C1(k,j,i)+dt*(vx(k,j,i)-C4(k,j,i)-C21(k,j,i) &
              +fx(k,j,i))*rmp
               vy (k,j,i) = C2(k,j,i)+dt*(vy(k,j,i)-C5(k,j,i)-C22(k,j,i) &
              +fy(k,j,i))*rmp
               vz (k,j,i) = C3(k,j,i)+dt*(vz(k,j,i)-C6(k,j,i)-C23(k,j,i) &
              +fz(k,j,i))*rmp
               rho(k,j,i) = C20(k,j,i)-dt*C7(k,j,i)*rmp
               th (k,j,i) = C21(k,j,i)+dt*(kappa*th(k,j,i)-C8(k,j,i)-C24(k,j,i))*rmp
            ELSE
               vx (k,j,i) = 0.
               vy (k,j,i) = 0.
               vz (k,j,i) = 0.
               rho(k,j,i) = 0.
               th (k,j,i) = 0.
            ENDIF

         END DO
         END DO
         END DO
