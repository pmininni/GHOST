! Step 2 of Runge-Kutta for the rotating GPE (RGPE) equation:
!
!    dz/dt = i [omegag.z - beta.|z|^2 z + alpha.Lap(z)
!               - i omegaz.(x.dz/dy - y.dz/dx) - Vtrap.z]
!
! where Vtrap(x,y) = V0.(x^2 + y^2) with V0 = m.w0^2/(2.hbar) and where w0 is
! the trapping frequency, resulting in a characteristic lengthscale for the
! trap a0 = sqrt(hbar/(m.w0)) = (cspeed.lambda/(V0.sqrt(2)))^(1/4). 
!
! Computes the nonlinear terms and evolves the equations at intermediate 
! steps of the R-K method. This solver only works with order ord=2 or 4.
! The real and imaginary parts of the equations are solved separately:
!    dzre/dt = - omegag.zim + beta.|z|^2 zim - alpha.Lap(zim)
!              + omegaz.(x.dzre/dy - y.dzre/dx) + Vtrap.zim
!    dzim/dt =   omegag.zre - beta.|z|^2 zre + alpha.Lap(zre)
!              + omegaz.(x.dzim/dy - y.dzim/dx) - Vtrap.zre

        IF ((trans.eq.1).and.(times.eq.0).and.(bench.eq.0).and.(o.eq.ord)) &
           CALL gperealspecc(zre,zim,iold,qold,kold,cold)
        
        CALL squareabs(zre,zim,R1,1)
        CALL nonlgpe(R1,zre,C3)
        CALL nonlgpe(R1,zim,C4)
        CALL laplak3(zre,C5)
        CALL laplak3(zim,C6)
        CALL derivk3(zre,C28,1)       ! dzre/dx
        CALL nonlgpe(Vliny,C28,C28)   ! omegaz*y*dzre/dx
        CALL derivk3(zre,C29,2)       ! dzre/dy
        CALL nonlgpe(Vlinx,C29,C29)   ! omegaz*x*dzre/dy
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
        DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
           DO j = 1,ny
              DO k = 1,nz             ! omegaz*(x.dzre/dy - y.dzre/dx)
                 C29(k,j,i) = C29(k,j,i)-C28(k,j,i)
              END DO
           END DO
        END DO
        CALL derivk3(zim,C28,1)       ! dzim/dx
        CALL nonlgpe(Vliny,C28,C28)   ! omegaz*y*dzim/dx
        CALL derivk3(zim,C30,2)       ! dzim/dy
        CALL nonlgpe(Vlinx,C30,C30)   ! omegaz*x*dzim/dy
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
        DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
           DO j = 1,ny
              DO k = 1,nz             ! omegaz*(x.dzim/dy - y.dzim/dx)
                 C30(k,j,i) = C30(k,j,i)-C28(k,j,i)
              END DO
           END DO
        END DO
        CALL nonlgpe(Vtrap,zre,C28)   ! Vtrap.zre
        CALL nonlgpe(Vtrap,zim,C31)   ! Vtrap.zim

        IF ( (ord-o).eq.0 ) THEN ! First RK iteration
           rmp = 0.5_GP
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
           DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
           DO j = 1,ny
           DO k = 1,nz
              IF (kn2(k,j,i).le.kmax) THEN
                 C7(k,j,i) = - omegag*zim(k,j,i) + beta*C4(k,j,i)  &
                   - alpha*C6(k,j,i) + C29(k,j,i) + C31(k,j,i) ! store K1
                 C8(k,j,i) =   omegag*zre(k,j,i) - beta*C3(k,j,i)  &
                   + alpha*C5(k,j,i) + C30(k,j,i) - C28(k,j,i) ! store K1
                 zre(k,j,i) = C1(k,j,i) + dt*C7(k,j,i)*rmp
                 zim(k,j,i) = C2(k,j,i) + dt*C8(k,j,i)*rmp
               ELSE
                 zre(k,j,i) = 0.0_GP
                 zim(k,j,i) = 0.0_GP
               ENDIF
           END DO
           END DO
           END DO
        ENDIF

        IF ( (ord-o).eq.1 ) THEN ! Second RK iteration
           IF ( ord.eq.2 ) rmp = 1.0_GP
           IF ( ord.eq.4 ) rmp = 0.5_GP
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
           DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
           DO j = 1,ny
           DO k = 1,nz
              IF (kn2(k,j,i).le.kmax) THEN
                 C4(k,j,i) = - omegag*zim(k,j,i) + beta*C4(k,j,i)  &
                   - alpha*C6(k,j,i) + C29(k,j,i) + C31(k,j,i)  ! K2
                 C5(k,j,i) =   omegag*zre(k,j,i) - beta*C3(k,j,i)  &
                   + alpha*C5(k,j,i) + C30(k,j,i) - C28(k,j,i)  ! K2
                 zre(k,j,i) = C1(k,j,i) + dt*C4(k,j,i)*rmp
                 zim(k,j,i) = C2(k,j,i) + dt*C5(k,j,i)*rmp  
              ELSE
                 zre(k,j,i) = 0.0_GP
                 zim(k,j,i) = 0.0_GP
              ENDIF
           END DO
           END DO
           END DO
           IF ( ord.eq.4 ) THEN
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
              DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
              DO j = 1,ny
              DO k = 1,nz
                 IF (kn2(k,j,i).le.kmax) THEN
                    C7(k,j,i) = C7(k,j,i) + 2*C4(k,j,i) ! K1+2.K2
                    C8(k,j,i) = C8(k,j,i) + 2*C5(k,j,i) ! K1+2.K2
                 ENDIF
              END DO
              END DO
              END DO
           ENDIF
        ENDIF

        IF ( (ord-o).eq.2 ) THEN ! Third RK iteration
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
           DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
           DO j = 1,ny
           DO k = 1,nz
              IF (kn2(k,j,i).le.kmax) THEN
                 C4(k,j,i) = - omegag*zim(k,j,i) + beta*C4(k,j,i)  &
                   - alpha*C6(k,j,i) + C29(k,j,i) + C31(k,j,i)  ! K3
                 C5(k,j,i) =   omegag*zre(k,j,i) - beta*C3(k,j,i)  &
                   + alpha*C5(k,j,i) + C30(k,j,i) - C28(k,j,i)  ! K3
                 C7(k,j,i) = C7(k,j,i) + 2*C4(k,j,i) ! K1+2.K2+2.K3
                 C8(k,j,i) = C8(k,j,i) + 2*C5(k,j,i) ! K1+2.K2+2.K3
                 zre(k,j,i) = C1(k,j,i) + dt*C4(k,j,i)
                 zim(k,j,i) = C2(k,j,i) + dt*C5(k,j,i)
              ELSE
                 zre(k,j,i) = 0.0_GP
                 zim(k,j,i) = 0.0_GP
              ENDIF
           END DO
           END DO
           END DO
        ENDIF

        IF ( (ord-o).eq.3 ) THEN ! Fourth RK iteration
           rmp = 1.0_GP/6
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
           DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
           DO j = 1,ny
           DO k = 1,nz
              IF (kn2(k,j,i).le.kmax) THEN
                 C7(k,j,i) = C7(k,j,i) + (- omegag*zim(k,j,i) + beta*C4(k,j,i) &
                   -alpha*C6(k,j,i) + C29(k,j,i) + C31(k,j,i)) ! K1+2.K2+2.K3+K4
                 C8(k,j,i) = C8(k,j,i) + (  omegag*zre(k,j,i) - beta*C3(k,j,i) &
                   +alpha*C5(k,j,i) + C30(k,j,i) - C28(k,j,i)) ! K1+2.K2+2.K3+K4
                 zre(k,j,i) = C1(k,j,i) + dt*C7(k,j,i)*rmp
                 zim(k,j,i) = C2(k,j,i) + dt*C8(k,j,i)*rmp
              ELSE
                 zre(k,j,i) = 0.0_GP
                 zim(k,j,i) = 0.0_GP
              ENDIF
           END DO
           END DO
           END DO
        ENDIF

        IF ((trans.eq.1).and.(times.eq.0).and.(bench.eq.0).and.(o.eq.1)) &
           THEN
           CALL gperealspecc(zre,zim,inew,qnew,knew,cnew)
           CALL gperealtrans(dt,iold,qold,kold,cold,inew,qnew,knew,cnew,ext)
        ENDIF
