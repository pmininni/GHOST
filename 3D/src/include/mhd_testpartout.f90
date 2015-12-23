! Write positions and velocities:
!
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
           WRITE(lgext, lgfmtext) pind
           CALL lagpart%io_write_pdb (1,odir,'xlg',lgext,(t-1)*dt)
           CALL lagpart%io_write_pdbv(1,odir,'vtp',lgext,(t-1)*dt)
           CALL lagpart%io_write_vec (1,odir,'vlg',lgext,(t-1)*dt)
           CALL lagpart%io_write_pdbm1(1,odir,'xlgm1',lgext,(t-2)*dt)
           CALL lagpart%io_write_vecm1(1,odir,'vlgm1',lgext,(t-2)*dt)

!!!!!! Write internal Lagrangian acceleration components: !!!!!!
!      NOTE: if dopacc > 0, then time centering of 'vlgm1', xlgm1' are
!            the same as 'algm1', while the other quantities are at the
!            most recent time.] 
           CALL lagpart%io_write_acc(tbeta,1,odir,'algm1',lgext,(t-2)*dt)
