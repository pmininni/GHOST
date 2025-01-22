! Spectra computed in MOIST runs 

            CALL spectrum(vx,vy,vz,ext,1,1)
            CALL specpara(vx,vy,vz,ext,1,1)
            CALL specperp(vx,vy,vz,ext,1,1)

            CALL spectrsc(th,ext,0)
            CALL specscpa(th,ext,0)
            CALL specscpe(th,ext,0)

            CALL spectrsc(th1,ext,1)
            CALL specscpa(th1,ext,1)
            CALL specscpe(th1,ext,1)

!!          CALL spectrsc(th2,ext,2)
!!          CALL specscpa(th2,ext,2)
!!          CALL specscpe(th2,ext,2)

!!          CALL spectrsc(th3,ext,3)
!!          CALL specscpa(th3,ext,3)
!!          CALL specscpe(th3,ext,3)

            CALL spec2D(vx,vy,vz,ext,odir,1,1)
            CALL specsc2D(th ,ext,odir,0)
            CALL specsc2D(th1,ext,odir,1)
!!          CALL specsc2D(th2,ext,odir,2)
!!          CALL specsc2D(th3,ext,odir,3)

! Uncomment the following lines to compute hor. averaged quantities
            CALL havgwrite(0,'shear'  ,ext,vx,vy,vz,th,0.0_GP,bvfreq) ! shear
            CALL havgwrite(1,'tgradz' ,ext,vx,vy,vz,th,0.0_GP,bvfreq) ! dtheta/dz
            CALL havgwrite(2,'hawdtdz',ext,vx,vy,vz,th,0.0_GP,bvfreq) ! u_z*dtheta/dz
            CALL havgwrite(3,'hahke'  ,ext,vx,vy,vz,th,0.0_GP,bvfreq) ! hor. k.e.
            CALL havgwrite(4,'havke'  ,ext,vx,vy,vz,th,0.0_GP,bvfreq) ! vert. k.e.
            CALL havgwrite(5,'haphel' ,ext,vx,vy,vz,th,0.0_GP,bvfreq) ! perp. helicity
            CALL havgwrite(6,'haomzt' ,ext,vx,vy,vz,th,0.0_GP,bvfreq) ! ometa_z*theta
            CALL havgwrite(7,'hapv2'  ,ext,vx,vy,vz,th,0.0_GP,bvfreq) ! pot'l vorticity^2
            CALL havgwrite(8,'hasuph' ,ext,vx,vy,vz,th,0.0_GP,bvfreq) ! super-helicity
            CALL havgwrite(9,'hari'   ,ext,vx,vy,vz,th,0.0_GP,bvfreq) ! Richardson no.

! Uncomment the following line to compute vert. spectrum of pot'l vorticity 
!           CALL spectpv(vx,vy,vz,th,ext)
