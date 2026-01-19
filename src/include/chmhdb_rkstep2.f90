! Step 2 of Runge-Kutta for the CHMHDB equations with uniform B_0
! Computes the nonlinear terms and evolves the equations in dt/o

         CALL rotor3(ay,az,C7,1)
         CALL rotor3(ax,az,C8,2)
         CALL rotor3(ax,ay,C9,3)
         IF (myrank.eq.0) THEN          ! b = b + B_0
            C7(1,1,1) = bx0*real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP)
            C8(1,1,1) = by0*real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP)
            C9(1,1,1) = bz0*real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP)
         ENDIF
         CALL vector3(vx,vy,vz,C7,C8,C9,C10,C11,C12)  ! v x B
         CALL prodre3(C7,C8,C9,C13,C14,C15)           ! j x B
         CALL divide(th,C13,C14,C15)                  ! (j x B)/rho
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend                             ! vxB - ep.(jxB)/rho
!$omp parallel do if (iend-ista.ge.nth) private (k)
         DO j = 1,ny
         DO k = 1,nz
            C10(k,j,i) = C10(k,j,i)-ep*C13(k,j,i)
            C11(k,j,i) = C11(k,j,i)-ep*C14(k,j,i)
            C12(k,j,i) = C12(k,j,i)-ep*C15(k,j,i)
         END DO
         END DO
         END DO
         CALL gauge3(C10,C11,C12,C7,1)               
         CALL gauge3(C10,C11,C12,C8,2)
         CALL gauge3(C10,C11,C12,C9,3)
         CALL prodre3(vx,vy,vz,C10,C11,C12)           ! om x v
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend                             ! c2*jxB/rho + v x om
!$omp parallel do if (iend-ista.ge.nth) private (k)
         DO j = 1,ny
         DO k = 1,nz
            C10(k,j,i) = cp2*C13(k,j,i)-C10(k,j,i)
            C11(k,j,i) = cp2*C14(k,j,i)-C11(k,j,i)
            C12(k,j,i) = cp2*C15(k,j,i)-C12(k,j,i)
         END DO
         END DO
         END DO
         CALL laplak3(ax,ax)
         CALL laplak3(ay,ay)
         CALL laplak3(az,az)
         CALL gradpress(cp1,gam1,th,vx,vy,vz,C13,C14,C15) ! grad pressure term
         PRINT *, myrank, SHAPE(th), SHAPE(C16)
         CALL divrhov(th,vx,vy,vz,0,C16)          ! div(rho.v)
         CALL vdiss(nu,nu2,vx,vy,vz)              ! viscous term
         CALL divide(th,vx,vy,vz)                 ! viscous term with rho
        
         rmp = 1./real(o,kind=GP)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.ge.nth) private (k)
         DO j = 1,ny
         DO k = 1,nz

            IF ((kn2(k,j,i).le.kmax)) THEN
               vx(k,j,i) = C1(k,j,i)+dt*(vx(k,j,i)+C10(k,j,i)-C13(k,j,i) &
              +fx(k,j,i))*rmp
               vy(k,j,i) = C2(k,j,i)+dt*(vy(k,j,i)+C11(k,j,i)-C14(k,j,i) &
              +fy(k,j,i))*rmp
               vz(k,j,i) = C3(k,j,i)+dt*(vz(k,j,i)+C12(k,j,i)-C15(k,j,i) &
              +fz(k,j,i))*rmp
               ax(k,j,i) = C4(k,j,i)+dt*(mu*ax(k,j,i)+C7(k,j,i)  &
              +mx(k,j,i))*rmp
               ay(k,j,i) = C5(k,j,i)+dt*(mu*ay(k,j,i)+C8(k,j,i)  &
              +my(k,j,i))*rmp
               az(k,j,i) = C6(k,j,i)+dt*(mu*az(k,j,i)+C9(k,j,i)  &
              +mz(k,j,i))*rmp
               th(k,j,i) = C20(k,j,i)-dt*C16(k,j,i)*rmp
            ELSE
               vx(k,j,i) = 0.0_GP
               vy(k,j,i) = 0.0_GP
               vz(k,j,i) = 0.0_GP
               ax(k,j,i) = 0.0_GP
               ay(k,j,i) = 0.0_GP
               az(k,j,i) = 0.0_GP
               th(k,j,i) = 0.0_GP
            ENDIF

         END DO
         END DO
         END DO
