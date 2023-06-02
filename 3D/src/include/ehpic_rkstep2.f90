         CALL derivk3(phi,C1,1)
         CALL derivk3(phi,C2,2)
         CALL derivk3(phi,C3,3)
         rmp = 1.0_GP/(real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.ge.nth) private (k)
         DO j = 1,ny
         DO k = 1,nz
            C1(k,j,i) = C1(k,j,i)*rmp
            C2(k,j,i) = C2(k,j,i)*rmp
            C3(k,j,i) = C3(k,j,i)*rmp
         END DO
         END DO
         END DO
         CALL fftp3d_complex_to_real(plancr,C1,Re1,MPI_COMM_WORLD)
         CALL fftp3d_complex_to_real(plancr,C2,Re2,MPI_COMM_WORLD)
         CALL fftp3d_complex_to_real(plancr,C3,Re3,MPI_COMM_WORLD)
         
         rmp = 1.0_GP/real(o,kind=GP)         
         CALL picpart%StepChargedPIC(Re1,Re2,Re3,Rb1,Rb2,Rb3,dt,rmp)

         CALL picpart%GetDensity(R1)
         CALL fftp3d_real_to_complex(planrc,R1,rhoc,MPI_COMM_WORLD)
         IF ( myrank.EQ.0 ) THEN
            rhoc(1,1,1) = 0.0_GP
         ENDIF

         CALL poisson_elecstat(rhoc,kde2,phi)
