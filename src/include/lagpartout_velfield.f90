!
!!!!!!! Write positions and velocities: !!!!!!

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
           WRITE(lgext,lgfmtext) pind

!       NOTE: if dopacc > 0 (set in solver), write out position and velocity corresp. 
!             to acceleration time stamp; change name to reflect that these lag by one 
!             dt from the current time:
           CALL lagpart%io_write_pdb  (1,odir,'xlg'  ,lgext,(t-1)*dt)
           CALL lagpart%io_write_vec  (1,odir,'vlg'  ,lgext,(t-1)*dt)
           CALL lagpart%io_write_pdbm1(1,odir,'xlgm1',lgext,(t-2)*dt)
           CALL lagpart%io_write_vecm1(1,odir,'vlgm1',lgext,(t-2)*dt)

! If doing 'fixed point' particles, set the Lagrangian velocities
! and repeat the process:
!
           IF ( blgdofp.GT.0 ) THEN

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
           CALL lagfp%SetLagVec  (R1,R2,R3,.true.,R4,R5,1)
           CALL lagfp%io_write_vec  (1,odir,'fpvlg'  ,lgext,(t-1)*dt)
           CALL lagfp%io_write_vecm1(1,odir,'fpvlgm1',lgext,(t-2)*dt)
           CALL lagfp%io_write_acc(tbeta,1,odir,'fpalgm1',lgext,(t-2)*dt)

           ENDIF

!!!!!! Write internal Lagrangian acceleration components: !!!!!!
!      NOTE: if dopacc > 0, then time centering of 'vlgm1', xlgm1' are
!            the same as 'algm1', while the other quantities are at the
!            most recent time.] 
           CALL lagpart%io_write_acc(tbeta,1,odir,'algm1',lgext,(t-2)*dt)
!
!!!!!! Write Lagrangian vorticity components: !!!!!!
!
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
           DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
              DO j = 1,ny
                 DO k = 1,nz
                    C1(k,j,i) = vx(k,j,i)*rmp
                    C2(k,j,i) = vy(k,j,i)*rmp
                    C3(k,j,i) = vz(k,j,i)*rmp
                 END DO
              END DO
           END DO
           CALL rotor3(C2,C3,C4,1)
           CALL rotor3(C1,C3,C5,2)
           CALL rotor3(C1,C2,C6,3)
           CALL fftp3d_complex_to_real(plancr,C4,R1,MPI_COMM_WORLD)
           CALL fftp3d_complex_to_real(plancr,C5,R2,MPI_COMM_WORLD)
           CALL fftp3d_complex_to_real(plancr,C6,R3,MPI_COMM_WORLD)
           CALL lagpart%SetLagVec(R1,R2,R3,.false.,R4,R5)
           CALL lagpart%io_write_vec(1,odir,'wlg'  ,lgext,(t-1)*dt)

           nwpart = nwpart + 1

!
!!!!!! Write nonlinear terms:!!!!!!
!
           if ( .false. ) then
           CALL prodre3(vx,vy,vz,C4,C5,C6)
           CALL fftp3d_complex_to_real(plancr,C4,R4,MPI_COMM_WORLD)
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
           DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
             DO j = 1,ny
                DO i = 1,nz
                   R4(i,j,k) = R4(i,j,k)*rmp
                END DO
             END DO
           END DO
           CALL lagpart%io_write_euler(1,odir,'v1nllg'  ,lgext,(t-1)*dt,R4,.false.,R2,R3)

           CALL fftp3d_complex_to_real(plancr,C5,R4,MPI_COMM_WORLD)
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
           DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
             DO j = 1,ny
                DO i = 1,nz
                   R4(i,j,k) = R4(i,j,k)*rmp
                END DO
             END DO
           END DO
           CALL lagpart%io_write_euler(1,odir,'v2nllg'  ,lgext,(t-1)*dt,R4,.false.,R2,R3)

           CALL fftp3d_complex_to_real(plancr,C6,R4,MPI_COMM_WORLD)
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
           DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
             DO j = 1,ny
                DO i = 1,nz
                   R4(i,j,k) = R4(i,j,k)*rmp
                END DO
             END DO
           END DO
           CALL lagpart%io_write_euler(1,odir,'v3nllg'  ,lgext,(t-1)*dt,R4,.false.,R2,R3)
           endif

!
!!!!!!! Write strain-rate tensor components: !!!!!!
!
           CALL derivk3(vx,C1,1)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
           DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
              DO j = 1,ny
                 DO k = 1,nz
                    C1(k,j,i) = C1(k,j,i)*rmp
                 END DO
              END DO
           END DO
           CALL fftp3d_complex_to_real(plancr,C1,R1,MPI_COMM_WORLD)
           CALL lagpart%io_write_euler(1,odir,'s11'  ,lgext,(t-1)*dt,R1,.false.,R2,R3)
           CALL derivk3(vx,C1,2)
           CALL derivk3(vy,C2,1)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
           DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
              DO j = 1,ny
                 DO k = 1,nz
                    C1(k,j,i) = 0.5*(C1(k,j,i)+C2(k,j,i))*rmp
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
              DO j = 1,ny
                 DO k = 1,nz
                    C1(k,j,i) = 0.5*(C1(k,j,i)+C2(k,j,i))*rmp
                 END DO
              END DO
           END DO
           CALL fftp3d_complex_to_real(plancr,C1,R1,MPI_COMM_WORLD)
           CALL lagpart%io_write_euler(1,odir,'s13',lgext,(t-1)*dt,R1,.false.,R2,R3)
           CALL derivk3(vy,C1,2)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
           DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
              DO j = 1,ny
                 DO k = 1,nz
                    C1(k,j,i) = C1(k,j,i)*rmp
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
              DO j = 1,ny
                 DO k = 1,nz
                    C1(k,j,i) = 0.5*(C1(k,j,i)+C2(k,j,i))*rmp
                 END DO
              END DO
           END DO
           CALL fftp3d_complex_to_real(plancr,C1,R1,MPI_COMM_WORLD)
           CALL lagpart%io_write_euler(1,odir,'s23',lgext,(t-1)*dt,R1,.false.,R2,R3)
