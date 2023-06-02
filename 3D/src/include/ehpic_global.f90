            CALL picpart%GetTemperature(R1)

            rmp = 1/(real(nx,KIND=GP)*real(ny,KIND=GP)*real(nz,KIND=GP))
            DO k=ksta,kend
              DO j=1,ny
                DO i=1,nx
                  ekin = ekin + R1(k,j,i)*rmp
                ENDDO
              ENDDO
            ENDDO

            CALL derivk3(phi,C1,1)
            CALL derivk3(phi,C2,2)
            CALL derivk3(phi,C3,3)
            CALL energy(C1,C2,C3,epot,1)

            IF (myrank.eq.0) THEN
               OPEN(1,file='energy.txt',position='append')
               WRITE(1,FMT='(E13.6,E22.14,E22.14)') (t-1)*dt,ekin,epot
               CLOSE(1)
            ENDIF    
