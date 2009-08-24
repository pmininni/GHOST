! Spectra computed in Hall-MHD runs

            CALL spectrum(vx,vy,vz,ext,1,1)
            CALL spectrum(ax,ay,az,ext,0,1)
            IF (gspe.eq.1) THEN
               DO i = ista,iend
                  DO j = 1,n
                     DO k = 1,n
                        C1(k,j,i) = ax(k,j,i)+ep*vx(k,j,i)
                        C2(k,j,i) = ay(k,j,i)+ep*vy(k,j,i)
                        C3(k,j,i) = az(k,j,i)+ep*vz(k,j,i)
                     END DO
                  END DO
               END DO
               CALL spectrum(C1,C2,C3,ext,2,1)
            ENDIF
