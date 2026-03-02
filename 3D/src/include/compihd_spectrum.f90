! Spectra computed in compressible HD runs

            CALL spectrum(vx,vy,vz,ext,1,1)
            CALL specpara(vx,vy,vz,ext,1,1)
            CALL specperp(vx,vy,vz,ext,1,1)

            CALL spectrsc(th,ext,0)


            CALL spectrsc(th,ext,0)
            CALL specscpa(th,ext,0)
            CALL specscpe(th,ext,0)

            CALL spectrsc(rho,ext,-1)
            CALL specscpa(rho,ext,-1)
            CALL specscpe(rho,ext,-1)
            CALL spec2D(vx,vy,vz,ext,odir,1,1)
            CALL specsc2D(th,ext,odir,0)

