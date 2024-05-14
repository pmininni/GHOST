            CALL picpart%GetDensity(R1)
            CALL picpart%GetFlux(Rj1,Rj2,Rj3)                
            CALL picpart%GetTemperature(R2)
!$omp parallel do if (iend-ista.ge.nth) private (j,i)
            DO k = ksta,kend
!$omp parallel do if (iend-ista.lt.nth) private (i)
            DO j = 1,ny
            DO i = 1,nx
               Rj1(i,j,k) = Rj1(i,j,k)/R1(i,j,k)
               Rj2(i,j,k) = Rj2(i,j,k)/R1(i,j,k)
               Rj3(i,j,k) = Rj3(i,j,k)/R1(i,j,k)
            END DO
            END DO
            END DO
            CALL hpiccheck(R2,Rj1,Rj2,Rj3,R1,ax,ay,az,gammae,betae,t,dt)
