! Spectra computed in BOUSS runs 

            CALL spectrum(vx,vy,vz,ext,1,1)
            CALL specpara(vx,vy,vz,ext,1,1)
            CALL specperp(vx,vy,vz,ext,1,1)
            CALL spectrsc(th,ext)
            CALL specscpa(th,ext)
            CALL specscpe(th,ext)
            CALL spec2D(vx,vy,vz,ext,odir,1,1)
            CALL specsc2D(th,ext,odir)
!
            CALL havgwrite(0,'hashear',ext,vx,vy,vz,th,omega,bvfreq) ! hor. avg shear
            CALL havgwrite(1,'hadtdz' ,ext,vx,vy,vz,th,omega,bvfreq) ! hor. avg dtheta/dz
            CALL havgwrite(2,'hawdtdz',ext,vx,vy,vz,th,omega,bvfreq) ! hor. avg u_z * dtheta/dz
            CALL havgwrite(3,'hahke'  ,ext,vx,vy,vz,th,omega,bvfreq) ! hor. avg hor. k.e.
            CALL havgwrite(4,'havke'  ,ext,vx,vy,vz,th,omega,bvfreq) ! hor. avg vert. k.e.
            CALL havgwrite(5,'haphel' ,ext,vx,vy,vz,th,omega,bvfreq) ! hor. avg perp. helicity
            CALL havgwrite(6,'haomzt' ,ext,vx,vy,vz,th,omega,bvfreq) ! hor. avg ometa_z * theta
            CALL havgwrite(7,'hapv2'  ,ext,vx,vy,vz,th,omega,bvfreq) ! hor. avg pot'l vorticity^2
            CALL havgwrite(8,'hasuph' ,ext,vx,vy,vz,th,omega,bvfreq) ! hor. avg super-helicity
!           CALL spectpv(vx,vy,vz,th,ext)
