! Spectra computed in HD runs in the rotating frame

            CALL spectrum(vx,vy,vz,ext,1,1)
            CALL specpara(vx,vy,vz,ext,1,1)
            CALL specperp(vx,vy,vz,ext,1,1)
!           CALL spec2D(vx,vy,vz,ext,odir,1,1)
!           CALL write_fourier(vx,'vx',ext,odir)
!           CALL write_fourier(vy,'vy',ext,odir)
!           CALL write_fourier(vz,'vz',ext,odir)
