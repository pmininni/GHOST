! Spectra computed in compressible MHD runs with B_0

            CALL spectrum(vx,vy,vz,ext,1,1)
            CALL specpara(vx,vy,vz,ext,1,1)
            CALL specperp(vx,vy,vz,ext,1,1)
            CALL spectrum(ax,ay,az,ext,0,1)
            CALL specpara(ax,ay,az,ext,0,1)
            CALL specperp(ax,ay,az,ext,0,1)
            CALL spectrsc(th,ext,0)
            CALL specscpa(th,ext)
            CALL specscpe(th,ext)
