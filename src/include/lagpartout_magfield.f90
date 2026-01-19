! Write positions and velocities:
!
! Set the Lagrangian velocities so output doesn't
! give 0: 
!
           rmp = 1.0_GP/(real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
           DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
             DO j = 1,ny
               DO k = 1,nz
                 C7(k,j,i) = vx(k,j,i)*rmp
               END DO
             END DO
           END DO
           CALL fftp3d_complex_to_real(plancr,C7,R1,MPI_COMM_WORLD)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
           DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
             DO j = 1,ny
               DO k = 1,nz
                 C7(k,j,i) = vy(k,j,i)*rmp
               END DO
             END DO
           END DO
           CALL fftp3d_complex_to_real(plancr,C7,R2,MPI_COMM_WORLD)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
           DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
             DO j = 1,ny
               DO k = 1,nz
                 C7(k,j,i) = vz(k,j,i)*rmp
               END DO
             END DO
           END DO
           CALL fftp3d_complex_to_real(plancr,C7,R3,MPI_COMM_WORLD)
           CALL lagpart%SetLagVec(R1,R2,R3,.true.,R4,R5)

           timep = 0
           pind = pind+1
           WRITE(lgext, lgfmtext) pind
           CALL lagpart%io_write_pdb (1,odir,'xlg',lgext,(t-1)*dt)
           CALL lagpart%io_write_vec (1,odir,'vlg',lgext,(t-1)*dt)

