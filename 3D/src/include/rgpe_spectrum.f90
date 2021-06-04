! Spectra computed in RGPE runs

            CALL gpemassspec(zre,zim,ext)
            CALL gperealspec(zre,zim,ext)
            CALL gpehelspec(zre,zim,ext)
            CALL gpemomtspec(zre,zim,ext)
            CALL gperealspecperp(zre,zim,ext)

! Uncomment the following lines to compute spatio-temporal spectra
!           CALL write_fourier(zre,'zre',ext,odir)
!           CALL write_fourier(zim,'zim',ext,odir)
