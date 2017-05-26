! Write scalar
           tmp = 1.0D0/(real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
           DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
             DO j = 1,ny
               DO k = 1,nz
                 C7(k,j,i) = th(k,j,i)*tmp
               END DO
             END DO
           END DO
           CALL fftp3d_complex_to_real(plancr,C7,R1,MPI_COMM_WORLD)
           CALL lagpart%io_write_euler(1,odir,'thlg',lgext,(t-1)*dt,R1,.true.,R2,R3)
           
           ! nonlinear term:
           CALL advect3(vx,vy,vz,th,C7)
           CALL fftp3d_complex_to_real(plancr,C7,R1,MPI_COMM_WORLD)
           CALL lagpart%io_write_euler(1,odir,'thnllg',lgext,(t-1)*dt,R1,.true.,R2,R3)
