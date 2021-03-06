! Step 2 of Runge-Kutta for the rotating ARGL (RARGL) equation:
!
!    dz/dt = omegag.z - beta.|z|^2 z + alpha.Lap(z) - i v.grad(z) -
!          - |v^2|.z/(4.alpha) + f - i omegaz.(x.dz/dy - y.dz/dx) - Vtrap.z
!
! where Vtrap(x,y) = V0.(x^2 + y^2) with V0 = m.w0^2/(2.hbar) and where w0 is
! the trapping frequency, resulting in a characteristic lengthscale for the
! trap a0 = sqrt(hbar/(m.w0)) = (cspeed.lambda/(V0.sqrt(2)))^(1/4). 
!
! Computes the nonlinear terms and evolves the equations using first 
! order forward implicit Euler. This solver only works with order ord=1.
! The real and imaginary parts of the equations are solved separately:
!    dzre/dt - alpha.Lap(zre) = omegag.zre - beta.|z|^2 zre + v.grad(zim)
!            - |v^2|.zre/(4.alpha) + fre + omegaz.(x.dzim/dy - y.dzim/dx)
!            - Vtrap(x,y).zre
!    dzim/dt - alpha.Lap(zre) = omegag.zim - beta.|z|^2 zim - v.grad(zre)
!            - |v^2|.zim/(4.alpha) + fim - omegaz.(x.dzre/dy - y.dzre/dx)
!            - Vtrap(x,y).zim

        IF (o.eq.1) THEN ! Only iterate once (first order)

        CALL squareabs(zre,zim,R1,1)
        rmq = real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP)*omegag/beta
        IF (cflow.eq.0) THEN ! If not doing counterflow we have the |v^2| term
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
           DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
              DO j = 1,ny
                 DO i = 1,nx
                    R1(i,j,k) = rmq - R1(i,j,k) - Vtrap(i,j,k) - vsq(i,j,k)
                 END DO
              END DO
           END DO
        ELSE
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
           DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
              DO j = 1,ny
                 DO i = 1,nx
                    R1(i,j,k) = rmq - R1(i,j,k) - Vtrap(i,j,k)
                 END DO
               END DO
           END DO
        ENDIF
        CALL nonlgpe(R1,zre,C3)
        CALL nonlgpe(R1,zim,C4)
        CALL advect3(vx,vy,vz,zre,C5) ! -v.grad(zre)
        CALL advect3(vx,vy,vz,zim,C6) ! -v.grad(zim)
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

!$omp parallel do if (iend-ista.ge.nth) private (j,k)
        DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
        DO j = 1,ny
        DO k = 1,nz
           IF (kn2(k,j,i).le.kmax) THEN
              rmp = 1.0_GP/(1.0_GP+alpha*kk2(k,j,i)*dt)
              zre(k,j,i) = (zre(k,j,i) + dt*(beta*C3(k,j,i) - C6(k,j,i) &
                          + fre(k,j,i) + C30(k,j,i)))*rmp
              zim(k,j,i) = (zim(k,j,i) + dt*(beta*C4(k,j,i) + C5(k,j,i) &
                          + fim(k,j,i) - C29(k,j,i)))*rmp
           ELSE
              zre(k,j,i) = 0.0_GP
              zim(k,j,i) = 0.0_GP
           ENDIF
        END DO
        END DO
        END DO

        IF ((t.eq.step).and.(iter_max_newt.gt.0)) THEN ! Do Newton at the end
           n_dim_1d = 4*(iend-ista+1)*ny*nz
           rmp = 1.0_GP/(real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))
           !$omp parallel do if (kend-ksta.ge.nth) private (j,i)
           DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
              DO j = 1,ny
                 DO i = 1,nx
                    R1(i,j,k) = Vlinx(i,j,k)*rmp ! Normalized omegaz.x
                    R2(i,j,k) = Vliny(i,j,k)*rmp ! Normalized omegaz.y
                END DO
              END DO
           END DO
           CALL fftp3d_real_to_complex(planrc,R1,C7,MPI_COMM_WORLD)
           CALL fftp3d_real_to_complex(planrc,R2,C8,MPI_COMM_WORLD)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
           DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
              DO j = 1,ny
                 DO k = 1,nz
                     C1(k,j,i) = vx(k,j,i) - C8(k,j,i) ! vx' = vx - omegaz.y
                     C2(k,j,i) = vy(k,j,i) + C7(k,j,i) ! vy' = vy + omegaz.x
                 END DO
              END DO
           END DO
           IF (cflow_newt.eq.0) THEN ! If not doing counterflow we have |v^2|
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
              DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
                 DO j = 1,ny
                    DO i = 1,nx
                       R2(i,j,k) = Vtrap(i,j,k) + vsq(i,j,k) ! |v^2| + Vtrap
                    END DO
                 END DO
              END DO
           ELSE
     !$omp parallel do if (kend-ksta.ge.nth) private (j,i)
              DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
                 DO j = 1,ny
                    DO i = 1,nx
                       R2(i,j,k) = Vtrap(i,j,k)              ! vsq' = Vtrap
                    END DO
                 END DO
              END DO
           ENDIF
           CALL newton(zre,zim,0,dt_newt,C1,C2,vz,R2,R1,C3,C4,C5,C6)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
           DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
              DO j = 1,ny
                 DO k = 1,nz
                     C3(k,j,i) = zre(k,j,i)*rmp
                     C4(k,j,i) = zim(k,j,i)*rmp
                 END DO
              END DO
           END DO
           CALL fftp3d_complex_to_real(plancr,C3,R1,MPI_COMM_WORLD)
           CALL fftp3d_complex_to_real(plancr,C4,R2,MPI_COMM_WORLD)
           CALL io_write(1,odir,'phi_re','NWT',planio,R1) ! Writes output
           CALL io_write(1,odir,'phi_im','NWT',planio,R2) ! in 'NWT' files
           IF (outs.ge.1) THEN
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
              DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
                 DO j = 1,ny
                    DO i = 1,nx
                       R1(i,j,k) = R1(i,j,k)**2+R2(i,j,k)**2
                    END DO
                 END DO
              END DO
              CALL io_write(1,odir,'rho','NWT',planio,R1)
           ENDIF
        ENDIF
        
        ! Renormalization for finite temperature runs
        IF (kttherm.gt.0) THEN
           ! Computes the mass
           CALL variance(zre,tmp,1)
           CALL variance(zim,tmq,1)
           IF (myrank.eq.0) tmr = tmp+tmq
           CALL MPI_BCAST(tmr,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

           ! Renormalization factor
           tmr = sqrt(omegag/beta)/sqrt(tmr)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
           DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
              DO j = 1,ny
                 DO k = 1,nz
                    zre(k,j,i) = zre(k,j,i)*tmr
                    zim(k,j,i) = zim(k,j,i)*tmr
                 END DO
              END DO
           END DO
        ENDIF

        ENDIF

