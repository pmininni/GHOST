! Spectra computed in HD runs

            CALL spectrum(vx,vy,vz,ext,1,1)
            CALL spectrsc(th,ext,0)

! Uncomment the following lines to compute 1D mean profiles
!           CALL difucx(vx,th,ext)
!           CALL difucz(vz,th,ext)
