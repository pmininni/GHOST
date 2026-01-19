! Spectra computed in Hybrid-PIC runs
            CALL spectrum(ax,ay,az,ext,0,0)
            CALL specpara(ax,ay,az,ext,0,0)
            CALL specperp(ax,ay,az,ext,0,0)
            
            CALL picpart%GetDensity(R1)
            CALL picpart%GetFlux(Rj1,Rj2,Rj3)  ! ni*u
!$omp parallel do if (iend-ista.ge.nth) private (j,i)
            DO k = ksta,kend
!$omp parallel do if (iend-ista.lt.nth) private (i)
            DO j = 1,ny
            DO i = 1,nx
               tmp = 1.0_GP/SQRT(R1(i,j,k))
               Rj1(i,j,k) = Rj1(i,j,k)*tmp
               Rj2(i,j,k) = Rj2(i,j,k)*tmp
               Rj3(i,j,k) = Rj3(i,j,k)*tmp ! sqrt(ni)*u
            END DO
            END DO
            END DO

            CALL fftp3d_real_to_complex(planrc,Rj1,ux,MPI_COMM_WORLD)
            CALL fftp3d_real_to_complex(planrc,Rj2,uy,MPI_COMM_WORLD)
            CALL fftp3d_real_to_complex(planrc,Rj3,uz,MPI_COMM_WORLD)
            CALL dealias(ux)
            CALL dealias(uy)
            CALL dealias(uz)
            CALL spectrum(ux,uy,uz,ext,1,0) ! sqrt(ni)*u spectra
            CALL specpara(ux,uy,uz,ext,1,0)
            CALL specperp(ux,uy,uz,ext,1,0)

            CALL fftp3d_real_to_complex(planrc,R1,rhoc,MPI_COMM_WORLD)
            CALL dealias(rhoc)
            CALL spectrsc(rhoc,ext,0)
            CALL specscpa(rhoc,ext,0)
            CALL specscpe(rhoc,ext,0)
