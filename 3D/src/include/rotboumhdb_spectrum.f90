! Spectra computed in BOUSS runs 

            CALL spectrum(vx,vy,vz,ext,1,1) ! Kinetic energy spectra
            CALL specpara(vx,vy,vz,ext,1,1)
            CALL specperp(vx,vy,vz,ext,1,1)
            CALL spectrum(ax,ay,az,ext,0,1) ! Magnetic energy spectra
            CALL specpara(ax,ay,az,ext,0,1)
            CALL specperp(ax,ay,az,ext,0,1)
            CALL spectrsc(th,ext,0)         ! Potential energy spectra
            CALL specscpa(th,ext,0)
            CALL specscpe(th,ext,0)
! Examples of how to write 1D longitudinal or transverse spectra
!           CALL spectr1d(vx,ext,'k',1,1)
!           CALL spectr1d(vy,ext,'k',2,1)
!           CALL spectr1d(vz,ext,'k',3,1)
!
! Uncomment to write 2D axisymmetric spectra
!           CALL spec2D(vx,vy,vz,ext,odir,1,1)
!           CALL spec2D(ax,ay,az,ext,odir,0,0)
!           CALL specsc2D(th,ext,odir,0)
!
! Uncomment the following lines to compute hor. averaged quantities
!  CALL havgwrite(0,'shear'  ,ext,vx,vy,vz,th,omega,bvfreq) ! shear
!  CALL havgwrite(1,'tgradz' ,ext,vx,vy,vz,th,omega,bvfreq) ! dtheta/dz
!  CALL havgwrite(2,'hawdtdz',ext,vx,vy,vz,th,omega,bvfreq) ! u_z*dtheta/dz
!  CALL havgwrite(3,'hahke'  ,ext,vx,vy,vz,th,omega,bvfreq) ! hor. k.e.
!  CALL havgwrite(4,'havke'  ,ext,vx,vy,vz,th,omega,bvfreq) ! vert. k.e.
!  CALL havgwrite(5,'haphel' ,ext,vx,vy,vz,th,omega,bvfreq) ! perp. helicity
!  CALL havgwrite(6,'haomzt' ,ext,vx,vy,vz,th,omega,bvfreq) ! ometa_z*theta
!  CALL havgwrite(7,'hapv2'  ,ext,vx,vy,vz,th,omega,bvfreq) ! pot'l vorticity^2
!  CALL havgwrite(8,'hasuph' ,ext,vx,vy,vz,th,omega,bvfreq) ! super-helicity
!  CALL havgwrite(9,'hari'   ,ext,vx,vy,vz,th,omega,bvfreq) ! Richardson no.
!
! Uncomment the following line to compute vert. spectrum of pot'l vorticity
!  CALL spectpv(vx,vy,vz,th,ext)
