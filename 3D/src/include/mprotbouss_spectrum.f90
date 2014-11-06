! Spectra computed in BOUSS runs 

            CALL spectrum(vx,vy,vz,ext,1,1)
            CALL specpara(vx,vy,vz,ext,1,1)
            CALL specperp(vx,vy,vz,ext,1,1)

            CALL spectrsc(th,ext)
            CALL specscpa(th,ext)
            CALL specscpe(th,ext)

            CALL spectrsc(th1,ext,1)
            CALL specscpa(th1,ext,1)
            CALL specscpe(th1,ext,1)

            CALL spectrsc(th2,ext,2)
            CALL specscpa(th2,ext,2)
            CALL specscpe(th2,ext,2)

!!          CALL spectrsc(th3,ext,3)
!!          CALL specscpa(th3,ext,3)
!!          CALL specscpe(th3,ext,3)

            CALL spec2D(vx,vy,vz,ext,odir,1,1)
            CALL specsc2D(th,ext,odir)
            CALL specsc2D(th1,ext,odir,1)
            CALL specsc2D(th2,ext,odir,2)
!!          CALL specsc2D(th3,ext,odir,3)
!
            CALL havgwrite(0,'shear'  ,ext,vx,vy,vz,th,omega,bvfreq) ! hor. avg shear
            CALL havgwrite(1,'tgradz' ,ext,vx,vy,vz,th,omega,bvfreq) ! hor. avg dtheta/dz
            CALL havgwrite(2,'hawdtdz',ext,vx,vy,vz,th,omega,bvfreq) ! hor. avg u_z * dtheta/dz
            CALL havgwrite(3,'hahke'  ,ext,vx,vy,vz,th,omega,bvfreq) ! hor. avg hor. k.e.
            CALL havgwrite(4,'havke'  ,ext,vx,vy,vz,th,omega,bvfreq) ! hor. avg vert. k.e.
            CALL havgwrite(5,'haphel' ,ext,vx,vy,vz,th,omega,bvfreq) ! hor. avg perp. helicity
            CALL havgwrite(6,'haomzt' ,ext,vx,vy,vz,th,omega,bvfreq) ! hor. avg ometa_z * theta
            CALL havgwrite(7,'hapv2'  ,ext,vx,vy,vz,th,omega,bvfreq) ! hor. avg pot'l vorticity^2
            CALL havgwrite(8,'hasuph' ,ext,vx,vy,vz,th,omega,bvfreq) ! hor. avg super-helicity
            CALL havgwrite(9,'hari'   ,ext,vx,vy,vz,th,omega,bvfreq) ! hor. avg Richardson no.
!           CALL spectpv(vx,vy,vz,th,ext)
