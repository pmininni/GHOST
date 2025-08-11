! Step 2 of Runge-Kutta for the Boussinesq equations 
! Computes the nonlinear terms and evolves the equations in dt/o

         CALL prodre3(vx,vy,vz,C4,C5,C6)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend               ! Gravity:
!$omp parallel do if (iend-ista.lt.nth) private (k)
            DO j = 1,ny
               DO k = 1,nz
                  C6(k,j,i) = C6(k,j,i)+xmom*th(k,j,i)
               END DO
            END DO
         END DO
         CALL nonlhd3(C4,C5,C6,C7,1)
         CALL nonlhd3(C4,C5,C6,C8,2)
         CALL nonlhd3(C4,C5,C6,C4,3)
         CALL advect3(vx,vy,vz,th,C5)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend               ! heat 'currrent':
!$omp parallel do if (iend-ista.lt.nth) private (k)
            DO j = 1,ny
               DO k = 1,nz
                  C5(k,j,i) = C5(k,j,i)+xtemp*vz(k,j,i)
               END DO
            END DO
         END DO
         CALL laplak3(vx,vx)
         CALL laplak3(vy,vy)
         CALL laplak3(vz,vz)
         CALL laplak3(th,th)

         IF ((trans.eq.1).and.(times.eq.0).and.(bench.eq.0).and.(o.eq.ord)) &
            THEN
            times = 0
            CALL entrans(C1,C2,C3,C7,C8,C4,ext,1)
            CALL entpara(C1,C2,C3,C7,C8,C4,ext,1)
            CALL entperp(C1,C2,C3,C7,C8,C4,ext,1)
            CALL sctrans(C20,C5,ext,0)
            CALL sctpara(C20,C5,ext,0)
            CALL sctperp(C20,C5,ext,0)
         ENDIF

         rmp = 1.0_GP/(real(o,kind=GP))

         if ( .NOT. use_voigt ) THEN

!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
         DO k = 1,nz
            IF ((kn2(k,j,i).le.kmax).and.(kn2(k,j,i).ge.tiny)) THEN
               vx(k,j,i) = C1 (k,j,i)+dt*(nu   *vx(k,j,i)+C7(k,j,i) &
              +fx(k,j,i))*rmp
               vy(k,j,i) = C2 (k,j,i)+dt*(nu   *vy(k,j,i)+C8(k,j,i) &
              +fy(k,j,i))*rmp
               vz(k,j,i) = C3 (k,j,i)+dt*(nu   *vz(k,j,i)+C4(k,j,i) &
              +fz(k,j,i))*rmp
               th(k,j,i) = C20(k,j,i)+dt*(kappa*th(k,j,i)+C5(k,j,i) &
              +fs(k,j,i))*rmp
            ELSE IF (kn2(k,j,i).gt.kmax) THEN
               vx(k,j,i) = 0.0_GP
               vy(k,j,i) = 0.0_GP
               vz(k,j,i) = 0.0_GP
               th(k,j,i) = 0.0_GP
            ELSE IF (kn2(k,j,i).lt.tiny) THEN
               vx(k,j,i) = 0.0_GP
               vy(k,j,i) = 0.0_GP
               vz(k,j,i) = 0.0_GP
               th(k,j,i) = C20(k,j,i)
            ENDIF
         END DO
         END DO
         END DO

         ELSE ! using Voigt

!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
         DO k = 1,nz
            IF ((kn2(k,j,i).le.kmax).and.(kn2(k,j,i).ge.tiny)) THEN
               vx(k,j,i) = C1 (k,j,i)+dt*(nu   *vx(k,j,i)+C7(k,j,i) &
              +fx(k,j,i))*rmp * Hinv(k,j,i)
               vy(k,j,i) = C2 (k,j,i)+dt*(nu   *vy(k,j,i)+C8(k,j,i) &
              +fy(k,j,i))*rmp * Hinv(k,j,i)
               vz(k,j,i) = C3 (k,j,i)+dt*(nu   *vz(k,j,i)+C4(k,j,i) &
              +fz(k,j,i))*rmp * Hinv(k,j,i)
               th(k,j,i) = C20(k,j,i)+dt*(kappa*th(k,j,i)+C5(k,j,i) &
              +fs(k,j,i))*rmp
            ELSE IF (kn2(k,j,i).gt.kmax) THEN
               vx(k,j,i) = 0.0_GP
               vy(k,j,i) = 0.0_GP
               vz(k,j,i) = 0.0_GP
               th(k,j,i) = 0.0_GP
            ELSE IF (kn2(k,j,i).lt.tiny) THEN
               vx(k,j,i) = 0.0_GP
               vy(k,j,i) = 0.0_GP
               vz(k,j,i) = 0.0_GP
               th(k,j,i) = C20(k,j,i)
            ENDIF
         END DO
         END DO
         END DO

         ENDIF ! Voigt check

         IF ( ord.EQ.3 .AND. o.LE.1 ) THEN
 
          IF ( use_voigt ) THEN
            WRITE(*,*) 'Cannot use higher order integration with Voigt'
            STOP
          ENDIF 

! First, compute u <-- L(u) + N(u,u):
          CALL prodre3(C1,C2,C3,C4,C5,C6)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
          DO i = ista,iend               ! Gravity:
!$omp parallel do if (iend-ista.lt.nth) private (k)
             DO j = 1,ny
                DO k = 1,nz
                   C6(k,j,i) = C6(k,j,i)+xmom*C20(k,j,i)
                END DO
            END DO
          END DO

          CALL nonlhd3(C4,C5,C6,C7,1)
          CALL nonlhd3(C4,C5,C6,C8,2)
          CALL nonlhd3(C4,C5,C6,C4,3)
          CALL laplak3(C1 ,C1 )
          CALL laplak3(C2 ,C2 )
          CALL laplak3(C3 ,C3 )

!$omp parallel do if (iend-ista.ge.nth) private (j,k)
          ! Compute L(u) + N(u,u):
          DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
            DO j = 1,ny
             DO k = 1,nz
              IF ((kn2(k,j,i).le.kmax).and.(kn2(k,j,i).ge.tiny)) THEN
                 C1 (k,j,i) = C1 (k,j,i)+C7(k,j,i)
                 C2 (k,j,i) = C2 (k,j,i)+C8(k,j,i)
                 C3 (k,j,i) = C3 (k,j,i)+C4(k,j,i)
              ELSE 
                 C1 (k,j,i) = 0.0
                 C2 (k,j,i) = 0.0
                 C3 (k,j,i) = 0.0
              ENDIF
             END DO
            END DO
          END DO

! Next, compute u <-- 2N(u,u)
          CALL prodre3(C1,C2,C3,C4,C5,C6)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
          DO i = ista,iend               ! Gravity:
!$omp parallel do if (iend-ista.lt.nth) private (k)
             DO j = 1,ny
                DO k = 1,nz
                   C6(k,j,i) = C6(k,j,i)+xmom*C20(k,j,i)
                END DO
            END DO
          END DO

          CALL nonlhd3(C4,C5,C6,C7,1)
          CALL nonlhd3(C4,C5,C6,C8,2)
          CALL nonlhd3(C4,C5,C6,C4,3)

          rmp = dt**3/24.0_GP

! Finally, u = u_* + dt^3/24 * u; where u_* is from the solution
! from the uncorrected JST above:
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
          DO i = ista,iend               ! Gravity:
!$omp parallel do if (iend-ista.lt.nth) private (k)
             DO j = 1,ny
                DO k = 1,nz
                   vx(k,j,i) = vx(k,j,i)+2.0_GP*C7(k,j,i)*rmp
                   vy(k,j,i) = vy(k,j,i)+2.0_GP*C8(k,j,i)*rmp
                   vz(k,j,i) = vz(k,j,i)+2.0_GP*C4(k,j,i)*rmp
                END DO
            END DO
          END DO

         ENDIF
