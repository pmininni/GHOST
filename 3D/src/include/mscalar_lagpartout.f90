! Write scalar
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
           DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
             DO j = 1,n
               DO k = 1,n
                 C7(k,j,i) = th1(k,j,i)/real(n,kind=GP)**3
               END DO
             END DO
           END DO
           CALL fftp3d_complex_to_real(plancr,C7,R1,MPI_COMM_WORLD)
           CALL lagpart%io_write_euler(1,odir,'th1lg',ext,(t-1)*dt,R1,.true.,R2,R3)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
           DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
             DO j = 1,n
               DO k = 1,n
                 C7(k,j,i) = th2(k,j,i)/real(n,kind=GP)**3
               END DO
             END DO
           END DO
           CALL fftp3d_complex_to_real(plancr,C7,R1,MPI_COMM_WORLD)
           CALL lagpart%io_write_euler(1,odir,'th2lg',ext,(t-1)*dt,R1,.true.,R2,R3)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
           DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
             DO j = 1,n
               DO k = 1,n
                 C7(k,j,i) = th3(k,j,i)/real(n,kind=GP)**3
               END DO
             END DO
           END DO
           CALL fftp3d_complex_to_real(plancr,C7,R1,MPI_COMM_WORLD)
           CALL lagpart%io_write_euler(1,odir,'th3lg',ext,(t-1)*dt,R1,.true.,R2,R3)

           ! nonlinear terms:
           CALL advect3(vx,vy,vz,th1,C7)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
           DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
             DO j = 1,n
               DO k = 1,n
                 C5(k,j,i) = th1(k,j,i)/real(n,kind=GP)**3
                 C7(k,j,i) = C7(k,j,i)/real(n,kind=GP)**3
               END DO
             END DO
           END DO
           CALL fftp3d_complex_to_real(plancr,C5,R1,MPI_COMM_WORLD)
           CALL fftp3d_complex_to_real(plancr,C7,R2,MPI_COMM_WORLD)
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
           DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
             DO j = 1,n
                DO i = 1,n
                   R1(i,j,k) = R1(i,j,k)*R2(i,j,k)
                END DO
             END DO
           END DO
           CALL lagpart%io_write_euler(1,odir,'s1translg',ext,(t-1)*dt,R1,.true.,R2,R3)


!           CALL advect3(vx,vy,vz,th2,C7)
!!$omp parallel do if (iend-ista.ge.nth) private (j,k)
!           DO i = ista,iend
!!$omp parallel do if (iend-ista.lt.nth) private (k)
!             DO j = 1,n
!               DO k = 1,n
!                 C5(k,j,i) = th2(k,j,i)/real(n,kind=GP)**3
!                 C7(k,j,i) = C7(k,j,i)/real(n,kind=GP)**3
!               END DO
!             END DO
!           END DO
!           CALL fftp3d_complex_to_real(plancr,C5,R1,MPI_COMM_WORLD)
!           CALL fftp3d_complex_to_real(plancr,C7,R2,MPI_COMM_WORLD)
!!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
!           DO k = ksta,kend
!!$omp parallel do if (kend-ksta.lt.nth) private (i)
!             DO j = 1,n
!                DO i = 1,n
!                   R1(i,j,k) = R1(i,j,k)*R2(i,j,k)
!                END DO
!             END DO
!           END DO
!           CALL lagpart%io_write_euler(1,odir,'s2translg',ext,(t-1)*dt,R1,.true.,R2,R3)
!
!
!           CALL advect3(vx,vy,vz,th3,C7)
!!$omp parallel do if (iend-ista.ge.nth) private (j,k)
!           DO i = ista,iend
!!$omp parallel do if (iend-ista.lt.nth) private (k)
!             DO j = 1,n
!               DO k = 1,n
!                 C5(k,j,i) = th3(k,j,i)/real(n,kind=GP)**3
!                 C7(k,j,i) = C7(k,j,i)/real(n,kind=GP)**3
!               END DO
!             END DO
!           END DO
!           CALL fftp3d_complex_to_real(plancr,C5,R1,MPI_COMM_WORLD)
!           CALL fftp3d_complex_to_real(plancr,C7,R2,MPI_COMM_WORLD)
!!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
!           DO k = ksta,kend
!!$omp parallel do if (kend-ksta.lt.nth) private (i)
!             DO j = 1,n
!                DO i = 1,n
!                   R1(i,j,k) = R1(i,j,k)*R2(i,j,k)
!                END DO
!             END DO
!           END DO
!           CALL lagpart%io_write_euler(1,odir,'s3translg',ext,(t-1)*dt,R1,.true.,R2,R3)
!
