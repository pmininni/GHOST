! Spectra computed in EDQNMHD runs

            CALL spectrum(vx,vy,vz,ext,1,1)
            CALL specpara(vx,vy,vz,ext,1,1)
            CALL specperp(vx,vy,vz,ext,1,1)
            call edqnmout(Eext,3*(n/2+1),dt*(t-1),'eext.'//ext//'.txt')
            call edqnmout(tepq,n/2+1    ,dt*(t-1),'tepq.'//ext//'.txt')
            call edqnmout(tve ,n/2+1    ,dt*(t-1),'tve.'//ext//'.txt')
            IF (heli.eq.1) THEN
            call edqnmout(Hext,3*(n/2+1),dt*(t-1),'hext.'//ext//'.txt')
            call edqnmout(thpq,n/2+1    ,dt*(t-1),'thpq.'//ext//'.txt')
            call edqnmout(tvh ,n/2+1    ,dt*(t-1),'tvh.'//ext//'.txt')
            ENDIF
