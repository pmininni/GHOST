!!!!!!! Write positions and velocities: !!!!!!
!
!!!!!!! Write positions and velocities: !!!!!!

! Set the Lagrangian velocities so output doesn't
! give 0: 
!
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
           DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
             DO j = 1,n
               DO k = 1,n
                 C7(k,j,i) = vx(k,j,i)/real(n,kind=GP)**3
               END DO
             END DO
           END DO
           CALL fftp3d_complex_to_real(plancr,C7,R1,MPI_COMM_WORLD)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
           DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
             DO j = 1,n
               DO k = 1,n
                 C7(k,j,i) = vy(k,j,i)/real(n,kind=GP)**3
               END DO
             END DO
           END DO
           CALL fftp3d_complex_to_real(plancr,C7,R2,MPI_COMM_WORLD)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
           DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
             DO j = 1,n
               DO k = 1,n
                 C7(k,j,i) = vz(k,j,i)/real(n,kind=GP)**3
               END DO
             END DO
           END DO
           CALL fftp3d_complex_to_real(plancr,C7,R3,MPI_COMM_WORLD)
           CALL lagpart%SetLagVec(R1,R2,R3,.true.,R4,R5)

           timep = 0
           pind = pind+1
           WRITE(lgext,lgfmtext) pind
           CALL lagpart%io_write_pdb(1,odir,'xlg',lgext,(t-1)*dt)
           CALL lagpart%io_write_vec(1,odir,'vlg',lgext,(t-1)*dt)

           ! Write Lagrangian vorticity components:
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
           DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
              DO j = 1,n
                 DO k = 1,n
                    C1(k,j,i) = vx(k,j,i)/real(n,kind=GP)**3
                    C2(k,j,i) = vy(k,j,i)/real(n,kind=GP)**3
                    C3(k,j,i) = vz(k,j,i)/real(n,kind=GP)**3
                 END DO
              END DO
           END DO
           CALL rotor3(C2,C3,C4,1)
           CALL rotor3(C1,C3,C5,2)
           CALL rotor3(C1,C2,C6,3)
           CALL fftp3d_complex_to_real(plancr,C4,R1,MPI_COMM_WORLD)
           CALL fftp3d_complex_to_real(plancr,C5,R2,MPI_COMM_WORLD)
           CALL fftp3d_complex_to_real(plancr,C6,R3,MPI_COMM_WORLD)


           CALL lagpart%SetLagVec(R1,R2,R3,.true.,R4,R5)
           CALL lagpart%io_write_vec(1,odir,'wlg',lgext,(t-1)*dt)

           nwpart = nwpart + 1

!
!!!!!! Write nonlinear transfer:!!!!!!
!
           CALL prodre3(vx,vy,vz,C4,C5,C6)
           CALL fftp3d_complex_to_real(plancr,C4,R4,MPI_COMM_WORLD)
           R5 = 0.0_GP
           tmp = 1.0D0/real(n,kind=GP)**3
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
           DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
             DO j = 1,n
                DO i = 1,n
                   R5(i,j,k) = R5(i,j,k)+ R1(i,j,k)*R4(i,j,k)*tmp
                END DO
             END DO
           END DO

           CALL fftp3d_complex_to_real(plancr,C5,R4,MPI_COMM_WORLD)
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
           DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
             DO j = 1,n
                DO i = 1,n
                   R5(i,j,k) = R5(i,j,k)+ R2(i,j,k)*R4(i,j,k)*tmp
                END DO
             END DO
           END DO

           CALL fftp3d_complex_to_real(plancr,C6,R4,MPI_COMM_WORLD)
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
           DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
             DO j = 1,n
                DO i = 1,n
                   R5(i,j,k) = R5(i,j,k)+ R3(i,j,k)*R4(i,j,k)*tmp
                END DO
             END DO
           END DO
           CALL lagpart%io_write_euler(1,odir,'etranslg',lgext,(t-1)*dt,R5,.false.,R2,R3)
!
!!!!!!! Write strain-rate tensor components: !!!!!!
!
           CALL derivk3(vx,C1,1)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
           DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
              DO j = 1,n
                 DO k = 1,n
                    C1(k,j,i) = C1(k,j,i)/real(n,kind=GP)**3
                 END DO
              END DO
           END DO
           CALL fftp3d_complex_to_real(plancr,C1,R1,MPI_COMM_WORLD)
           CALL lagpart%io_write_euler(1,odir,'s11',lgext,(t-1)*dt,R1,.false.,R2,R3)
           CALL derivk3(vx,C1,2)
           CALL derivk3(vy,C2,1)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
           DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
              DO j = 1,n
                 DO k = 1,n
                    C1(k,j,i) = 0.5*(C1(k,j,i)+C2(k,j,i))/real(n,kind=GP)**3
                 END DO
              END DO
           END DO
           CALL fftp3d_complex_to_real(plancr,C1,R1,MPI_COMM_WORLD)
           CALL lagpart%io_write_euler(1,odir,'s12',lgext,(t-1)*dt,R1,.false.,R2,R3)
           CALL derivk3(vx,C1,3)
           CALL derivk3(vz,C2,1)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
           DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
              DO j = 1,n
                 DO k = 1,n
                    C1(k,j,i) = 0.5*(C1(k,j,i)+C2(k,j,i))/real(n,kind=GP)**3
                 END DO
              END DO
           END DO
           CALL fftp3d_complex_to_real(plancr,C1,R1,MPI_COMM_WORLD)
           CALL lagpart%io_write_euler(1,odir,'s13',lgext,(t-1)*dt,R1,.false.,R2,R3)
           CALL derivk3(vy,C1,2)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
           DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
              DO j = 1,n
                 DO k = 1,n
                    C1(k,j,i) = C1(k,j,i)/real(n,kind=GP)**3
                 END DO
              END DO
           END DO
           CALL fftp3d_complex_to_real(plancr,C1,R1,MPI_COMM_WORLD)
           CALL lagpart%io_write_euler(1,odir,'s22',lgext,(t-1)*dt,R1,.false.,R2,R3)
           CALL derivk3(vy,C1,3)
           CALL derivk3(vz,C2,2)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
           DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
              DO j = 1,n
                 DO k = 1,n
                    C1(k,j,i) = 0.5*(C1(k,j,i)+C2(k,j,i))/real(n,kind=GP)**3
                 END DO
              END DO
           END DO
           CALL fftp3d_complex_to_real(plancr,C1,R1,MPI_COMM_WORLD)
           CALL lagpart%io_write_euler(1,odir,'s23',lgext,(t-1)*dt,R1,.false.,R2,R3)
