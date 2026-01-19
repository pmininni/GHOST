! Spectra computed in PROTH runs (HD rotation with passive scalar)

            CALL spectrum(vx,vy,vz,ext,1,1)
            CALL specpara(vx,vy,vz,ext,1,1)
            CALL specperp(vx,vy,vz,ext,1,1)
            CALL spectrsc(th,ext,0)
            CALL specscpa(th,ext,0)
            CALL specscpe(th,ext,0)

! Uncomment the following lines to compute 1D mean profiles
!           CALL difucx(vx,th,ext)
!           CALL difucz(vz,th,ext)

! Uncomment the following lines to compute 2D axisymmetric spectra
!           CALL spec2D(vx,vy,vz,ext,odir,1,1)
!           CALL specsc2D(th,ext,odir,0)
