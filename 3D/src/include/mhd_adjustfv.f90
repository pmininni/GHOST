! Adjust the amplitude of the mechanical forcing
! to keep the energy constant in a MHD or Hall-MHD run

               CALL energy(vx,vy,vz,tmp,0)
               CALL energy(ax,ay,az,tmq,2)
               CALL cross(vx,vy,vz,fx,fy,fz,eps,1)
               CALL cross(ax,ay,az,mx,my,mz,epm,0)
               CALL MPI_BCAST(tmp,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
               CALL MPI_BCAST(tmq,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
               CALL MPI_BCAST(eps,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
               CALL MPI_BCAST(epm,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
               Faux1(1:9) = Faux1(2:10)
               Faux2(1:9) = Faux2(2:10)
               ampl = ampl*(nu*tmp+mu*tmq-epm)/eps
               Faux1(10) = ampl
               Faux2(10) = nu*tmp+mu*tmq-epm-eps
               ampl = 0.9_GP*ampl
               DO i = 1,9
                  ampl = ampl+(0.1_GP*Faux1(i)+0.01_GP*Faux2(i))/9.0_GP
               END DO
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
               DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
                  DO j = 1,n
                     DO k = 1,n
                        fx(k,j,i) = fx(k,j,i)*ampl
                        fy(k,j,i) = fy(k,j,i)*ampl
                        fz(k,j,i) = fz(k,j,i)*ampl
                     END DO
                  END DO
               END DO

