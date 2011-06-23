! Spectra computed in BOUSS runs 

            CALL spectrum(vx,vy,vz,ext,1,1)
            CALL specpara(vx,vy,vz,ext,1,1)
            CALL specperp(vx,vy,vz,ext,1,1)
            CALL spectrsc(th,ext)
            CALL specscpa(th,ext)
            CALL specscpe(th,ext)
!           CALL spec2D(vx,vy,vz,ext,odir,1,1)
!           CALL specsc2D(th,ext,odir)
            CALL havghvel(vx,vy,ext)
            CALL spectpv(vx,vy,vz,th,ext)
