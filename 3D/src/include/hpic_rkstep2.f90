         CALL picpart%GetDensity(R1)
         CALL picpart%GetFlux(Rj1,Rj2,Rj3)                

         CALL fftp3d_real_to_complex(planrc,Rj1,C4,MPI_COMM_WORLD)
         CALL fftp3d_real_to_complex(planrc,Rj2,C5,MPI_COMM_WORLD)
         CALL fftp3d_real_to_complex(planrc,Rj3,C6,MPI_COMM_WORLD)
         CALL fftp3d_real_to_complex(planrc,R1 ,C7,MPI_COMM_WORLD)
         CALL rotor3(ay,az,C8 ,1) ! bx
         CALL rotor3(ax,az,C9 ,2) ! by
         CALL rotor3(ax,ay,C10,3) ! bz
         tmp = real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP)
         IF (myrank.eq.0) THEN          ! b = b + B_0
            C8(1,1,1) = bx0*tmp
            C9(1,1,1) = by0*tmp
            C10(1,1,1)= bz0*tmp
         ENDIF
         CALL laplak3(ax,ax) ! -jx
         CALL laplak3(ay,ay) ! -jy
         CALL laplak3(az,az) ! -jz

         tmp = 1.0_GP/dii
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.ge.nth) private (k)
         DO j = 1,ny
         DO k = 1,nz
            IF ((kn2(k,j,i).le.kmax)) THEN
               C4(k,j,i) = dii*ax(k,j,i) + C4(k,j,i)!*C17(k,j,i)
               C5(k,j,i) = dii*ay(k,j,i) + C5(k,j,i)!*C17(k,j,i)
               C6(k,j,i) = dii*az(k,j,i) + C6(k,j,i)!*C17(k,j,i)
               ax(k,j,i) = C4(k,j,i)*tmp  ! Store electron current j_e
               ay(k,j,i) = C5(k,j,i)*tmp
               az(k,j,i) = C6(k,j,i)*tmp
               C7(k,j,i) = C7(k,j,i)!*C17(k,j,i)
            ELSE
               C4(k,j,i) = 0.0_GP
               C5(k,j,i) = 0.0_GP
               C6(k,j,i) = 0.0_GP
               C7(k,j,i) = 0.0_GP
            ENDIF
         END DO
         END DO
         END DO

         IF (gammae.EQ.1) THEN
            CALL derivk3(C7,C14,1)
            CALL derivk3(C7,C15,2)
            CALL derivk3(C7,C16,3)     ! grad(P_e)
            CALL divide(C7,C14,C15,C16)
         ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,i)
            DO k = ksta,kend
!$omp parallel do if (iend-ista.lt.nth) private (i)
            DO j = 1,ny
            DO i = 1,nx
               R2(i,j,k) = R1(i,j,k)**gam1
            END DO
            END DO
            END DO
            CALL fftp3d_real_to_complex(planrc,R2,C11,MPI_COMM_WORLD)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
            DO i = ista,iend
!$omp parallel do if (iend-ista.ge.nth) private (k)
            DO j = 1,ny
            DO k = 1,nz
               IF ((kn2(k,j,i).le.kmax)) THEN
                  C11(k,j,i) = C11(k,j,i)!*C17(k,j,i)
               ELSE
                  C11(k,j,i) = 0.0_GP
               ENDIF
            END DO
            END DO
            END DO
            CALL derivk3(C11,C14,1)
            CALL derivk3(C11,C15,2)
            CALL derivk3(C11,C16,3) ! grad(P_e)/n_e
         END IF       
        
         CALL divide(C7,C4,C5,C6)  ! u_e
         CALL dealias(C4)
         CALL dealias(C5)
         CALL dealias(C6)
         CALL vector3(C4,C5,C6,C8,C9,C10,C11,C12,C13) ! u_exB
         CALL gauge3(C11,C12,C13,C4,1)
         CALL gauge3(C11,C12,C13,C5,2)
         CALL gauge3(C11,C12,C13,C6,3) ! u_exB - grad(Phi)

         rmp = 1.0_GP/real(o*Bmult,kind=GP)
         tmp = 1/(dii*real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.ge.nth) private (k)
         DO j = 1,ny
         DO k = 1,nz
            IF ((kn2(k,j,i).le.kmax)) THEN
               C8(k,j,i) = C8(k,j,i) *tmp
               C9(k,j,i) = C9(k,j,i) *tmp
               C10(k,j,i)= C10(k,j,i)*tmp
               C11(k,j,i)= (-C11(k,j,i) - mu*ax(k,j,i) - cp1*C14(k,j,i)&
                            - 0*mx(k,j,i))*tmp
               C12(k,j,i)= (-C12(k,j,i) - mu*ay(k,j,i) - cp1*C15(k,j,i)&
                            - 0*my(k,j,i))*tmp
               C13(k,j,i)= (-C13(k,j,i) - mu*az(k,j,i) - cp1*C16(k,j,i)&
                            - 0*mz(k,j,i))*tmp
               ax(k,j,i) = C1(k,j,i) + dt*(mu*ax(k,j,i) +              &
                                           C4(k,j,i) + mx(k,j,i))*rmp
               ay(k,j,i) = C2(k,j,i) + dt*(mu*ay(k,j,i) +              &
                                           C5(k,j,i) + my(k,j,i))*rmp
               az(k,j,i) = C3(k,j,i) + dt*(mu*az(k,j,i) +              &
                                           C6(k,j,i) + mz(k,j,i))*rmp
            ELSE
               C8(k,j,i) = 0.0_GP
               C9(k,j,i) = 0.0_GP
               C10(k,j,i)= 0.0_GP
               C11(k,j,i)= 0.0_GP
               C12(k,j,i)= 0.0_GP
               C13(k,j,i)= 0.0_GP 
               ax(k,j,i) = 0.0_GP
               ay(k,j,i) = 0.0_GP
               az(k,j,i) = 0.0_GP
            ENDIF
         END DO
         END DO
         END DO

         CALL fftp3d_complex_to_real(plancr,C8 ,Rb1,MPI_COMM_WORLD)
         CALL fftp3d_complex_to_real(plancr,C9 ,Rb2,MPI_COMM_WORLD)
         CALL fftp3d_complex_to_real(plancr,C10,Rb3,MPI_COMM_WORLD)
         CALL fftp3d_complex_to_real(plancr,C11,Re1,MPI_COMM_WORLD)
         CALL fftp3d_complex_to_real(plancr,C12,Re2,MPI_COMM_WORLD)
         CALL fftp3d_complex_to_real(plancr,C13,Re3,MPI_COMM_WORLD)

         CALL picpart%StepChargedPICBor(Re1,Re2,Re3,Rb1,Rb2,Rb3,dt,o)

         DO iB = 2,Bmult
            CALL fftp3d_real_to_complex(planrc,Rj1,C4,MPI_COMM_WORLD)
            CALL fftp3d_real_to_complex(planrc,Rj2,C5,MPI_COMM_WORLD)
            CALL fftp3d_real_to_complex(planrc,Rj3,C6,MPI_COMM_WORLD)
            CALL rotor3(ay,az,C8 ,1) ! bx
            CALL rotor3(ax,az,C9 ,2) ! by
            CALL rotor3(ax,ay,C10,3) ! bz
            tmp = real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP)
            IF (myrank.eq.0) THEN          ! b = b + B_0
               C8(1,1,1) = bx0*tmp
               C9(1,1,1) = by0*tmp
               C10(1,1,1)= bz0*tmp
            ENDIF
            CALL laplak3(ax,C14) ! -jx
            CALL laplak3(ay,C15) ! -jy
            CALL laplak3(az,C16) ! -jz
   
            tmp = 1.0_GP/dii
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
            DO i = ista,iend
!$omp parallel do if (iend-ista.ge.nth) private (k)
            DO j = 1,ny
            DO k = 1,nz
               IF ((kn2(k,j,i).le.kmax)) THEN
                  C4(k,j,i)  = dii*C14(k,j,i) + C4(k,j,i)!*C17(k,j,i)
                  C5(k,j,i)  = dii*C15(k,j,i) + C5(k,j,i)!*C17(k,j,i)
                  C6(k,j,i)  = dii*C16(k,j,i) + C6(k,j,i)!*C17(k,j,i)
                  C14(k,j,i) = C4(k,j,i)*tmp  ! Store electron current -j_e
                  C15(k,j,i) = C5(k,j,i)*tmp
                  C16(k,j,i) = C6(k,j,i)*tmp
               ELSE
                  C4(k,j,i) = 0.0_GP
                  C5(k,j,i) = 0.0_GP
                  C6(k,j,i) = 0.0_GP
               ENDIF
            END DO
            END DO
            END DO
            
            CALL divide(C7,C4,C5,C6)  ! u_e
            CALL dealias(C4)
            CALL dealias(C5)
            CALL dealias(C6)
            CALL vector3(C4,C5,C6,C8,C9,C10,C11,C12,C13) ! u_exB
            CALL gauge3(C11,C12,C13,C4,1)
            CALL gauge3(C11,C12,C13,C5,2)
            CALL gauge3(C11,C12,C13,C6,3) ! u_exB - grad(Phi)
   
            rmp = 1.0_GP/real(o*Bmult,kind=GP)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
            DO i = ista,iend
!$omp parallel do if (iend-ista.ge.nth) private (k)
            DO j = 1,ny
            DO k = 1,nz
               IF ((kn2(k,j,i).le.kmax)) THEN
                  ax(k,j,i) = ax(k,j,i) + dt*(mu*C14(k,j,i) +          &
                                              C4(k,j,i) + mx(k,j,i))*rmp
                  ay(k,j,i) = ay(k,j,i) + dt*(mu*C15(k,j,i) +          &
                                              C5(k,j,i) + my(k,j,i))*rmp
                  az(k,j,i) = az(k,j,i) + dt*(mu*C16(k,j,i) +          &
                                              C6(k,j,i) + mz(k,j,i))*rmp
               ELSE
                  ax(k,j,i) = 0.0_GP
                  ay(k,j,i) = 0.0_GP
                  az(k,j,i) = 0.0_GP
               ENDIF
            END DO
            END DO
            END DO
         END DO 

