! Write scalar
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
           DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
             DO j = 1,n
               DO k = 1,n
                 C7(k,j,i) = th(k,j,i)/real(n,kind=GP)**3
               END DO
             END DO
           END DO
           CALL fftp3d_complex_to_real(plancr,C7,R1,MPI_COMM_WORLD)
           CALL lagpart%io_write_euler(1,odir,'thlg',ext,(t-1)*dt,R1,.true.,R2,R3)
