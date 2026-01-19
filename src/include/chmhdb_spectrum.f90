! Spectra computed in compressible Hall-MHD runs with B_0

            CALL spectrum(vx,vy,vz,ext,1,1)
            CALL specpara(vx,vy,vz,ext,1,1)
            CALL specperp(vx,vy,vz,ext,1,1)
            CALL spectrum(ax,ay,az,ext,0,1)
            CALL specpara(ax,ay,az,ext,0,1)
            CALL specperp(ax,ay,az,ext,0,1)
            CALL spectrsc(th,ext,0)
            CALL specscpa(th,ext,0)
            CALL specscpe(th,ext,0)
            IF (gspe.eq.1) THEN
               DO i = ista,iend
                  DO j = 1,ny
                     DO k = 1,nz
                        C1(k,j,i) = ax(k,j,i)+ep*vx(k,j,i)
                        C2(k,j,i) = ay(k,j,i)+ep*vy(k,j,i)
                        C3(k,j,i) = az(k,j,i)+ep*vz(k,j,i)
                     END DO
                  END DO
               END DO
               CALL spectrum(C1,C2,C3,ext,2,1)
               CALL specpara(C1,C2,C3,ext,2,1)
               CALL specperp(C1,C2,C3,ext,2,1)
            ENDIF

