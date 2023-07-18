         CALL picpart%GetDensity(R1)
         CALL picpart%GetFlux(Re1,Re2,Re3)        
         CALL fftp3d_real_to_complex(planrc,Re1,C4,MPI_COMM_WORLD)
         CALL fftp3d_real_to_complex(planrc,Re2,C5,MPI_COMM_WORLD)
         CALL fftp3d_real_to_complex(planrc,Re3,C6,MPI_COMM_WORLD)
         CALL fftp3d_real_to_complex(planrc,R1 ,C7,MPI_COMM_WORLD)

         CALL rotor3(ax,ay,C8 ,1) ! bx
         CALL rotor3(ay,az,C9 ,2) ! by
         CALL rotor3(az,ax,C10,3) ! bz
         IF (myrank.eq.0) THEN          ! b = b + B_0
             C7(1,1,1) = bx0*real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP)
             C8(1,1,1) = by0*real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP)
             C9(1,1,1) = bz0*real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP)
         ENDIF
         CALL laplak3(ax,ax) ! -jx
         CALL laplak3(ay,ay) ! -jy
         CALL laplak3(az,az) ! -jz

!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.ge.nth) private (k)
         DO j = 1,ny
         DO k = 1,nz
            IF ((kn2(k,j,i).le.kmax)) THEN
               C4(k,j,i) = ax(k,j,i) + C4(k,j,i)
               C5(k,j,i) = ay(k,j,i) + C5(k,j,i)
               C6(k,j,i) = az(k,j,i) + C6(k,j,i)
               C7(k,j,i) = C7(k,j,i)
            ELSE
               C4(k,j,i) = 0.0_GP
               C5(k,j,i) = 0.0_GP
               C6(k,j,i) = 0.0_GP
               C7(k,j,i) = 0.0_GP
            ENDIF
         END DO
         END DO
         END DO
         
         CALL divide(C7,C4,C5,C6)  ! -u_e
         CALL vector3(C4,C5,C6,C8,C9,C10,C11,C12,C13) ! -u_exB
         CALL gauge3(C11,C12,C13,C4,1)
         CALL gauge3(C11,C12,C13,C5,2)
         CALL gauge3(C11,C12,C13,C6,3) ! -u_exB + grad(Phi)

         rmp = 1./real(o,kind=GP)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.ge.nth) private (k)
         DO j = 1,ny
         DO k = 1,nz
            IF ((kn2(k,j,i).le.kmax)) THEN
               ax(k,j,i) = C1(k,j,i) + dt*(mu*ax(k,j,i) - C4(k,j,i))*rmp
               ay(k,j,i) = C2(k,j,i) + dt*(mu*ay(k,j,i) - C5(k,j,i))*rmp
               az(k,j,i) = C3(k,j,i) + dt*(mu*az(k,j,i) - C6(k,j,i))*rmp
            ELSE
               ax(k,j,i) = 0.0_GP
               ay(k,j,i) = 0.0_GP
               az(k,j,i) = 0.0_GP
            ENDIF
         END DO
         END DO
         END DO


!$omp parallel do if (iend-ista.ge.nth) private (j,i)
         DO k = ksta,kend
!$omp parallel do if (iend-ista.lt.nth) private (i)
         DO j = 1,ny
         DO i = 1,nx
            R1(i,j,k) = cp1*R1(i,j,k)**gam1 ! P_e
         END DO
         END DO
         END DO

         CALL fftp3d_real_to_complex(planrc,R1,C7,MPI_COMM_WORLD)
         CALL derivk3(C7,C8 ,1)
         CALL derivk3(C7,C9 ,2)
         CALL derivk3(C7,C10,3)
   
         tmp = 1.0_GP/(real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.ge.nth) private (k)
         DO j = 1,ny
         DO k = 1,nz
            IF ((kn2(k,j,i).le.kmax)) THEN
               C7(k,j,i)  = C7(k,j,i)*tmp
               C8(k,j,i)  = C8(k,j,i)*tmp
               C9(k,j,i)  = C9(k,j,i)*tmp
               C11(k,j,i) = (C11(k,j,i) - mu*ax(k,j,i) -  C8(k,j,i))*tmp
               C12(k,j,i) = (C12(k,j,i) - mu*ay(k,j,i) -  C9(k,j,i))*tmp
               C13(k,j,i) = (C13(k,j,i) - mu*az(k,j,i) - C10(k,j,i))*tmp
            ELSE
               C7(k,j,i) = 0.0_GP
               C8(k,j,i) = 0.0_GP
               C9(k,j,i) = 0.0_GP
               C11(k,j,i) = 0.0_GP
               C12(k,j,i) = 0.0_GP
               C13(k,j,i) = 0.0_GP 
            ENDIF
         END DO
         END DO
         END DO

         CALL fftp3d_complex_to_real(plancr,C7 ,Rb1,MPI_COMM_WORLD)
         CALL fftp3d_complex_to_real(plancr,C8 ,Rb2,MPI_COMM_WORLD)
         CALL fftp3d_complex_to_real(plancr,C9 ,Rb3,MPI_COMM_WORLD)
         CALL fftp3d_complex_to_real(plancr,C11,Re1,MPI_COMM_WORLD)
         CALL fftp3d_complex_to_real(plancr,C12,Re2,MPI_COMM_WORLD)
         CALL fftp3d_complex_to_real(plancr,C13,Re3,MPI_COMM_WORLD)

         CALL picpart%StepChargedPIC(Re1,Re2,Re3,Rb1,Rb2,Rb3,dt,rmp)
