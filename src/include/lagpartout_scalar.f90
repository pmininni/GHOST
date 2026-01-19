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
         
!!
!!!!!!!! Write temperature gradient components: !!!!!!
!!
!           rmp = 1.0_GP/(real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))
!           CALL derivk3(th,C1,1)
!!$omp parallel do if (iend-ista.ge.nth) private (j,k)
!           DO i = ista,iend
!!$omp parallel do if (iend-ista.lt.nth) private (k)
!              DO j = 1,ny
!                 DO k = 1,nz
!                    C1(k,j,i) = C1(k,j,i)*rmp
!                END DO
!              END DO
!           END DO
!           CALL fftp3d_complex_to_real(plancr,C1,R1,MPI_COMM_WORLD)
!           CALL lagpart%io_write_euler(1,odir,'dth1'  ,lgext,(t-1)*dt,R1,.false.,R2,R3)
!
!           CALL derivk3(th,C2,2)
!!$omp parallel do if (iend-ista.ge.nth) private (j,k)
!           DO i = ista,iend
!!$omp parallel do if (iend-ista.lt.nth) private (k)
!              DO j = 1,ny
!                 DO k = 1,nz
!                    C2(k,j,i) = C2(k,j,i)*rmp
!                 END DO
!              END DO
!           END DO
!           CALL fftp3d_complex_to_real(plancr,C2,R1,MPI_COMM_WORLD)
!           CALL lagpart%io_write_euler(1,odir,'dth2'  ,lgext,(t-1)*dt,R1,.false.,R2,R3)
!
!           CALL derivk3(th,C3,3)
!!$omp parallel do if (iend-ista.ge.nth) private (j,k)
!           DO i = ista,iend
!!$omp parallel do if (iend-ista.lt.nth) private (k)
!              DO j = 1,ny
!                 DO k = 1,nz
!                    C3(k,j,i) = C3(k,j,i)*rmp
!                 END DO
!              END DO
!           END DO
!           CALL fftp3d_complex_to_real(plancr,C3,R1,MPI_COMM_WORLD)
!           CALL lagpart%io_write_euler(1,odir,'dth3'  ,lgext,(t-1)*dt,R1,.false.,R2,R3)
