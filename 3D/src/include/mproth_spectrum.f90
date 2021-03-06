! Spectra computed in MPROTH runs (HD rotation with 2 or 3 passive scalar)

            CALL spectrum(vx,vy,vz,ext,1,1)
            CALL specpara(vx,vy,vz,ext,1,1)
            CALL specperp(vx,vy,vz,ext,1,1)
            CALL spectrsc(th1,ext,1)
            CALL specscpa(th1,ext,1)
            CALL specscpe(th1,ext,1)
            CALL spectrsc(th2,ext,2)
            CALL specscpa(th2,ext,2)
            CALL specscpe(th2,ext,2)
!!          CALL spectrsc(th3,ext,3)
!!          CALL specscpa(th3,ext,3)
!!          CALL specscpe(th3,ext,3)

! Uncomment the following lines to compute 2D axisymmetric spectra
!           CALL spec2D(vx,vy,vz,ext,odir,1,1)
!           CALL specsc2D(th1,ext,odir,1)
!           CALL specsc2D(th2,ext,odir,2)
!           CALL specsc2D(th3,ext,odir,3)
