! Write velocity of test particles
           CALL lagpart%io_write_pdbv(1,odir,'vtp',lgext,(t-1)*dt)
! Write magnetic field and current density:
           tmp = 1.0E0 / real(n,kind=GP)**3
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
           DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
              DO j = 1,n
                 DO k = 1,n
                    C1(k,j,i) = ax(k,j,i)*tmp
                    C2(k,j,i) = ay(k,j,i)*tmp
                    C3(k,j,i) = az(k,j,i)*tmp
                 END DO
              END DO
           END DO
           CALL rotor3(C2,C3,C4,1)
           CALL fftp3d_complex_to_real(plancr,C4,R1,MPI_COMM_WORLD)
           CALL rotor3(C1,C3,C5,2)
           CALL fftp3d_complex_to_real(plancr,C5,R2,MPI_COMM_WORLD)
           CALL rotor3(C1,C2,C6,3)
           CALL fftp3d_complex_to_real(plancr,C6,R3,MPI_COMM_WORLD)
           CALL lagpart%SetLagVec(R1,R2,R3,.false.,R4,R5)
           CALL lagpart%io_write_vec(1,odir,'blg',lgext,(t-1)*dt)

           CALL laplak3(C1,C4)
           CALL laplak3(C2,C5)
           CALL laplak3(C3,C6)
           CALL fftp3d_complex_to_real(plancr,C4,R1,MPI_COMM_WORLD)
           CALL fftp3d_complex_to_real(plancr,C5,R2,MPI_COMM_WORLD)
           CALL fftp3d_complex_to_real(plancr,C6,R3,MPI_COMM_WORLD)
           CALL lagpart%SetLagVec(R1,R2,R3,.false.,R4,R5)
           CALL lagpart%io_write_vec(1,odir,'jlg',lgext,(t-1)*dt)
