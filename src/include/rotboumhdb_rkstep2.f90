! Step 2 of Runge-Kutta for the Boussinesq + MHD equations 
! Computes the nonlinear terms and evolves the equations in dt/o

         CALL rotor3(ay,az,C7,1)        ! b = curl(a)
         CALL rotor3(ax,az,C8,2)
         CALL rotor3(ax,ay,C9,3)
         IF (myrank.eq.0) THEN          ! b = b + B_0
            C7(1,1,1) = bx0*real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP)
            C8(1,1,1) = by0*real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP)
            C9(1,1,1) = bz0*real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP)
         ENDIF  
         CALL prodre3(C7,C8,C9,C13,C14,C15) ! j x b
         CALL prodre3(vx,vy,vz,C10,C11,C12) ! w x v
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend               ! NL = w x v + 2 Om x v + N.theta
!$omp parallel do if (iend-ista.lt.nth) private (k)
          DO j = 1,ny
           DO k = 1,nz
            C10(k,j,i) = C10(k,j,i)+2*(omegay*vz(k,j,i)-omegaz*vy(k,j,i))
            C11(k,j,i) = C11(k,j,i)+2*(omegaz*vx(k,j,i)-omegax*vz(k,j,i))
            C12(k,j,i) = C12(k,j,i)+2*(omegax*vy(k,j,i)-omegay*vx(k,j,i)) &
                        +xmom*th(k,j,i)
           END DO
          END DO
         END DO
         CALL nonlin3(C10,C11,C12,C13,C14,C15,C16,1) ! -NL -j x b -grad(p)
         CALL nonlin3(C10,C11,C12,C13,C14,C15,C17,2)
         CALL nonlin3(C10,C11,C12,C13,C14,C15,C10,3)
         CALL vector3(vx,vy,vz,C7,C8,C9,C10,C11,C12) ! v x b
         CALL gauge3(C10,C11,C12,C7,1)  ! v x b - grad(phi)
         CALL gauge3(C10,C11,C12,C8,2)
         CALL gauge3(C10,C11,C12,C9,3)
         CALL advect3(vx,vy,vz,th,C11)  ! - v.grad(theta)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend               ! heat 'currrent'
!$omp parallel do if (iend-ista.lt.nth) private (k)
            DO j = 1,ny
               DO k = 1,nz
                  C11(k,j,i) = C11(k,j,i)+xtemp*vz(k,j,i)
               END DO
            END DO
         END DO
         CALL laplak3(vx,vx)            ! laplacian(v)
         CALL laplak3(vy,vy)
         CALL laplak3(vz,vz)
         CALL laplak3(ax,ax)            ! laplacian(a)
         CALL laplak3(ay,ay)
         CALL laplak3(az,az)
         CALL laplak3(th,th)            ! laplacial(theta)
         
         IF ((trans.eq.1).and.(times.eq.0).and.(bench.eq.0).and.(o.eq.ord)) &
            THEN
            CALL entrans(C1,C2,C3,C16,C17,C10,ext,1)
            CALL entrans(C4,C5,C6, C7, C8, C9,ext,0)
            CALL sctrans(C20,C11,ext,0)
         ENDIF

         rmp = 1.0_GP/(real(o,kind=GP))
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
         DO k = 1,nz
            IF ((kn2(k,j,i).le.kmax).and.(kn2(k,j,i).ge.tiny)) THEN
               vx(k,j,i) = C1 (k,j,i)+dt*(nu   *vx(k,j,i)+C16(k,j,i) &
              +fx(k,j,i))*rmp
               vy(k,j,i) = C2 (k,j,i)+dt*(nu   *vy(k,j,i)+C17(k,j,i) &
              +fy(k,j,i))*rmp
               vz(k,j,i) = C3 (k,j,i)+dt*(nu   *vz(k,j,i)+C10(k,j,i) &
              +fz(k,j,i))*rmp
               ax(k,j,i) = C4 (k,j,i)+dt*(mu   *ax(k,j,i)+C7 (k,j,i) &
              +mx(k,j,i))*rmp
               ay(k,j,i) = C5 (k,j,i)+dt*(mu   *ay(k,j,i)+C8 (k,j,i) &
              +my(k,j,i))*rmp
               az(k,j,i) = C6 (k,j,i)+dt*(mu   *az(k,j,i)+C9 (k,j,i) &
              +mz(k,j,i))*rmp
               th(k,j,i) = C20(k,j,i)+dt*(kappa*th(k,j,i)+C11(k,j,i) &
              +fs(k,j,i))*rmp
            ELSE IF (kn2(k,j,i).gt.kmax) THEN
               vx(k,j,i) = 0.0_GP
               vy(k,j,i) = 0.0_GP
               vz(k,j,i) = 0.0_GP
               ax(k,j,i) = 0.0_GP
               ay(k,j,i) = 0.0_GP
               az(k,j,i) = 0.0_GP
               th(k,j,i) = 0.0_GP
            ELSE IF (kn2(k,j,i).lt.tiny) THEN
               vx(k,j,i) = 0.0_GP
               vy(k,j,i) = 0.0_GP
               vz(k,j,i) = 0.0_GP
               ax(k,j,i) = 0.0_GP
               ay(k,j,i) = 0.0_GP
               az(k,j,i) = 0.0_GP
               th(k,j,i) = C20(k,j,i)
            ENDIF
         END DO
         END DO
         END DO
