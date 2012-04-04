! Spectra computed in SWHD runs

            CALL spectrsc(th,ext)
            CALL swspectrum(vx,vy,th,fs,g,ext,1)
            CALL swspectrum(vx,vy,th,fs,g,ext,2)
            CALL vspectrum(vx,vy,ext)
