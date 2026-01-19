! Step 2 of Runge-Kutta for the compressible HD equations
! Computes the nonlinear terms and evolves the equations in dt/o

         CALL prodre3(vx,vy,vz,C4,C5,C6)          ! om x v
         CALL gradpress(cp1,gam1,th,vx,vy,vz,C25,C26,C27) ! grad pressure term
         CALL divrhov(th,vx,vy,vz,0,C7)           ! div(rho.v)
         CALL vdiss(nu,nu2,vx,vy,vz)              ! viscous term
         CALL divide(th,vx,vy,vz)                 ! viscous term with rho
        
         rmp = 1./real(o,kind=GP)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
         DO k = 1,nz

            IF (kn2(k,j,i).le.kmax) THEN
               vx(k,j,i) = C1(k,j,i)+dt*(vx(k,j,i)-C4(k,j,i)-C25(k,j,i) &
              +fx(k,j,i))*rmp
               vy(k,j,i) = C2(k,j,i)+dt*(vy(k,j,i)-C5(k,j,i)-C26(k,j,i) &
              +fy(k,j,i))*rmp
               vz(k,j,i) = C3(k,j,i)+dt*(vz(k,j,i)-C6(k,j,i)-C27(k,j,i) &
              +fz(k,j,i))*rmp
               th(k,j,i) = C20(k,j,i)-dt*C7(k,j,i)*rmp
            ELSE
               vx(k,j,i) = 0.
               vy(k,j,i) = 0.
               vz(k,j,i) = 0.
               th(k,j,i) = 0.
            ENDIF

         END DO
         END DO
         END DO
