!=================================================================
! GHOST code: Geophysical High Order Suite for Turbulence
!
! Header file with definitions for conditional compilation.
! Main code is in main3D.fpp. See that file for more details.
!
! 2015 Pablo D. Mininni.
!      Department of Physics, 
!      Facultad de Ciencias Exactas y Naturales.
!      Universidad de Buenos Aires.
!      e-mail: mininni@df.uba.ar
!=================================================================

! The lines below define preprocessor variables for each solver.
! New solvers can be created by adding new definitions here, 
! adding the necessary solver files in 'include', and listing
! all objects in 'Makefile'.

! Fluid solvers

#ifdef HD_SOL
#define DNS_
#define VELOC_
#define INCLUDEFNAME_ 'hd_
#endif

#ifdef PHD_SOL
#define DNS_
#define VELOC_
#define SCALAR_
#define INCLUDEFNAME_ 'phd_
#endif

#ifdef MPHD_SOL
#define DNS_
#define VELOC_
#define MULTISCALAR_
#define INCLUDEFNAME_ 'mphd_
#endif

#ifdef MHD_SOL
#define DNS_
#define VELOC_
#define MAGFIELD_
#define INCLUDEFNAME_ 'mhd_
#endif

#ifdef MHDB_SOL
#define DNS_
#define VELOC_
#define MAGFIELD_
#define UNIFORMB_
#define INCLUDEFNAME_ 'mhdb_
#endif

#ifdef RMHDB_SOL
#define DNS_
#define VELOC_
#define MAGFIELD_
#define UNIFORMB_
#define ROTATION_
#define INCLUDEFNAME_ 'rmhdb_
#endif

#ifdef HMHD_SOL
#define DNS_
#define VELOC_
#define MAGFIELD_
#define HALLTERM_
#define INCLUDEFNAME_ 'hmhd_
#endif

#ifdef HMHDB_SOL
#define DNS_
#define VELOC_
#define MAGFIELD_
#define HALLTERM_
#define UNIFORMB_
#define INCLUDEFNAME_ 'hmhdb_
#endif

#ifdef COMPRHD_SOL
#define DNS_
#define VELOC_
#define SCALAR_
#define COMPRESSIBLE_
#define COMPR_AUX_ARR_
#define INCLUDEFNAME_ 'comprhd_
#endif

#ifdef COMPIHD_SOL
#define DNS_
#define VELOC_
#define MOM_
#define DENSITY_
#define SCALAR_
#define COMPRESSIBLE_
#define COMPI_AUX_ARR_
#define INCLUDEFNAME_ 'compihd_
#endif

#ifdef CMHD_SOL
#define DNS_
#define CMHD_
#define VELOC_
#define SCALAR_
#define MAGFIELD_
#define COMPRESSIBLE_
#define INCLUDEFNAME_ 'cmhd_
#endif

#ifdef CMHDB_SOL
#define DNS_
#define CMHD_
#define VELOC_
#define SCALAR_
#define MAGFIELD_
#define UNIFORMB_
#define COMPRESSIBLE_
#define INCLUDEFNAME_ 'cmhdb_
#endif

#ifdef CHMHD_SOL
#define DNS_
#define CMHD_
#define VELOC_
#define SCALAR_
#define MAGFIELD_
#define HALLTERM_
#define COMPRESSIBLE_
#define INCLUDEFNAME_ 'chmhd_
#endif

#ifdef CHMHDB_SOL
#define DNS_
#define CMHD_
#define VELOC_
#define SCALAR_
#define MAGFIELD_
#define HALLTERM_
#define UNIFORMB_
#define COMPRESSIBLE_
#define INCLUDEFNAME_ 'chmhdb_
#endif

#ifdef ROTH_SOL
#define DNS_
#define VELOC_
#define ROTATION_
#define INCLUDEFNAME_ 'roth_
#endif 

#ifdef PROTH_SOL
#define DNS_
#define VELOC_
#define SCALAR_
#define ROTATION_
#define INCLUDEFNAME_ 'proth_
#endif

#ifdef MPROTH_SOL
#define DNS_
#define VELOC_
#define ROTATION_
#define MULTISCALAR_
#define INCLUDEFNAME_ 'mproth_
#endif

#ifdef BOUSS_SOL
#define DNS_
#define VELOC_
#define SCALAR_
#define BOUSSINESQ_
#define INCLUDEFNAME_ 'bouss_
#endif

#ifdef MPBOUSS_SOL
#define DNS_
#define VELOC_
#define SCALAR_
#define BOUSSINESQ_
#define MULTISCALAR_
#define INCLUDEFNAME_ 'mpbouss_
#endif

#ifdef ROTBOUSS_SOL
#define DNS_
#define VELOC_
#define SCALAR_
#define ROTATION_
#define BOUSSINESQ_
#define INCLUDEFNAME_ 'rotbouss_
#endif

#ifdef ROTBOUSSSGS_SOL
#define DNS_
#define VELOC_
#define VELOCSGS_
#define SCALAR_
#define SCALARSGS_
#define ROTATION_
#define BOUSSINESQ_
#define INCLUDEFNAME_ 'rotbouss_
#endif

#ifdef ROTBOUMHDB_SOL
#define DNS_
#define VELOC_
#define SCALAR_
#define MAGFIELD_
#define UNIFORMB_
#define ROTATION_
#define BOUSSINESQ_
#define INCLUDEFNAME_ 'rotboumhdb_
#endif

#ifdef MPROTBOUSS_SOL
#define DNS_
#define VELOC_
#define SCALAR_
#define ROTATION_
#define BOUSSINESQ_
#define MULTISCALAR_
#define INCLUDEFNAME_ 'mprotbouss_
#endif

#ifdef MOIST_SOL
#define DNS_
#define VELOC_
#define SCALAR_
#define BOUSSINESQ_
#define MULTISCALAR_
#define INCLUDEFNAME_ 'moist_
#endif 

#ifdef HDPNLT_SOL
#define DNS_
#define VELOC_
#define PENALTY_
#define INCLUDEFNAME_ 'hdpnlt_
#endif

#ifdef GPE_SOL
#define DNS_
#define WAVEFUNCTION_
#define INCLUDEFNAME_ 'gpe_
#endif

#ifdef ARGL_SOL
#define DNS_
#define NEWTON_
#define ADVECT_
#define QFORCE_
#define WAVEFUNCTION_
#define INCLUDEFNAME_ 'argl_
#endif

#ifdef RGPE_SOL
#define DNS_
#define WAVEFUNCTION_
#define TRAP_
#define ROTATION_
#define INCLUDEFNAME_ 'rgpe_
#endif

#ifdef RARGL_SOL
#define DNS_
#define ADVECT_
#define QFORCE_
#define WAVEFUNCTION_
#define TRAP_
#define ROTATION_
#define INCLUDEFNAME_ 'rargl_
#endif

#ifdef LAHD_SOL
#define ALPHAV_
#define VELOC_
#define INCLUDEFNAME_ 'lahd_
#endif

#ifdef CAHD_SOL
#define ALPHAV_
#define VELOC_
#define INCLUDEFNAME_ 'cahd_
#endif

#ifdef LHD_SOL
#define ALPHAV_
#define VELOC_
#define INCLUDEFNAME_ 'lhd_

#endif

#ifdef LAMHD_SOL
#define ALPHAV_
#define ALPHAB_
#define VELOC_
#define MAGFIELD_
#define INCLUDEFNAME_ 'lamhd_
#endif

#ifdef EDQNMHD_SOL
#define DNS_
#define EDQNM_
#define VELOC_
#define INCLUDEFNAME_ 'edqnmhd_
#endif

#ifdef EDQNMROTH_SOL
#define DNS_
#define EDQNM_
#define VELOC_
#define ROTATION_
#define INCLUDEFNAME_ 'edqnmroth_
#endif

#ifdef EHPIC_SOL
#define DNS_
#define ELECSTAT_
#define ELECFIELD_
#define UNIFORMB_
#define CPIC_
#define INCLUDEFNAME_ 'ehpic_
#endif

#ifdef HPIC_SOL
#define DNS_
#define HYBPIC_
#define MAGFIELD_
#define UNIFORMB_
#define CPIC_
#define INCLUDEFNAME_ 'hpic_
#endif

! Particles subclasses

#ifdef DEF_GHOST_LAGP
#define PART_
#define LAGPART_
#endif

#ifdef DEF_GHOST_INERP
#define PART_
#define VPART_
#define INERPART_ 
#endif

#ifdef DEF_GHOST_TESTP
#define PART_ 
#define VPART_
#define TESTPART_ 
#endif

#ifdef DEF_GHOST_PICP
#define PIC_
#endif

! Do not edit below this line!
! Builds the names of all files to include for each solver
#define STRINGIFY(a) a
#define SOLVERCHECK_ STRINGIFY(INCLUDEFNAME_)validate.f90'
#define GLOBALOUTPUT_ STRINGIFY(INCLUDEFNAME_)global.f90'
#define SPECTROUTPUT_ STRINGIFY(INCLUDEFNAME_)spectrum.f90'
#define RKSTEP1_ STRINGIFY(INCLUDEFNAME_)rkstep1.f90'
#define RKSTEP2_ STRINGIFY(INCLUDEFNAME_)rkstep2.f90'
