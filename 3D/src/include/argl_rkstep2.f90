! Step 2 of Runge-Kutta for the ARGL equation:
!    dz/dt = omegag.z - beta.|z|^2 z + alpha.Lap(z) - i v.grad(z) - |v^2|.z/(4.alpha) + f
! Computes the nonlinear terms and evolves the equations using first 
! order forward implicit Euler. This solver only works with order ord=1.
! The real and imaginary parts of the equations are solved separately:
!    dzre/dt - alpha.Lap(zre) = omegag.zre - beta.|z|^2 zre + v.grad(zim) - |v^2|.zre/(4.alpha) + fre
!    dzim/dt - alpha.Lap(zre) = omegag.zim - beta.|z|^2 zim - v.grad(zre) - |v^2|.zim/(4.alpha) + fim

         IF (o.eq.1) THEN ! Only iterate once (first order)

         CALL squareabs(zre,zim,R1,1)
         IF (cflow.eq.0) THEN ! If not doing counterflow we have the |v^2| term
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
            DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
               DO j = 1,n
                  DO i = 1,n
                     R1(i,j,k) = -R1(i,j,k)-vsq(i,j,k)
                  END DO
               END DO
            END DO
         ENDIF
         CALL nonlgpe(R1,zre,C3)
         CALL nonlgpe(R1,zim,C4)
         CALL advect3(vx,vy,vz,zre,C5) ! -v.grad(zre)
         CALL advect3(vx,vy,vz,zim,C6) ! -v.grad(zim)

         rmq = 1.0_GP+dt*omegag
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,n
         DO k = 1,n
            IF (ka2(k,j,i).le.kmax) THEN
               rmp = 1.0_GP/(1.0_GP+alpha*ka2(k,j,i)*dt)
               zre(k,j,i) = (rmq*zre(k,j,i)+dt*(beta*C3(k,j,i) &
                           -C6(k,j,i)+fre(k,j,i)))*rmp
               zim(k,j,i) = (rmq*zim(k,j,i)+dt*(beta*C4(k,j,i) &
                           +C5(k,j,i)+fim(k,j,i)))*rmp
            ELSE
               zre(k,j,i) = 0.0_GP
               zim(k,j,i) = 0.0_GP
            ENDIF
         END DO
         END DO
         END DO

         ENDIF
