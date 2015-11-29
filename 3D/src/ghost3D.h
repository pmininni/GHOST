!=================================================================
! GHOST code: Geophysical High Order Suite for Turbulence
!
! Header file with definitions for conditional compilation.
! Main code is in main2D.fpp. See this file for more details.
!
! 2015 Pablo D. Mininni.
!      Department of Physics, 
!      Facultad de Ciencias Exactas y Naturales.
!      Universidad de Buenos Aires.
!      e-mail: mininni@df.uba.ar
!=================================================================

#ifdef HD_SOL
#define DNS_
#define VELOC_
#endif

#ifdef PHD_SOL
#define DNS_
#define VELOC_
#define SCALAR_
#endif

#ifdef MPHD_SOL
#define DNS_
#define VELOC_
#define MULTISCALAR_
#endif

#ifdef MHD_SOL
#define DNS_
#define VELOC_
#define MAGFIELD_
#endif

#ifdef MHDB_SOL
#define DNS_
#define VELOC_
#define MAGFIELD_
#define UNIFORMB_
#endif

#ifdef HMHD_SOL
#define DNS_
#define VELOC_
#define MAGFIELD_
#define HALLTERM_
#endif

#ifdef HMHDB_SOL
#define DNS_
#define VELOC_
#define MAGFIELD_
#define HALLTERM_
#define UNIFORMB_
#endif

#ifdef ROTH_SOL
#define DNS_
#define VELOC_
#define ROTATION_
#endif 

#ifdef PROTH_SOL
#define DNS_
#define VELOC_
#define SCALAR_
#define ROTATION_
#endif

#ifdef MPROTH_SOL
#define DNS_
#define VELOC_
#define ROTATION_
#define MULTISCALAR_
#endif

#ifdef BOUSS_SOL
#define DNS_
#define VELOC_
#define SCALAR_
#define BOUSSINESQ_
#endif

#ifdef MPBOUSS_SOL
#define DNS_
#define VELOC_
#define SCALAR_
#define BOUSSINESQ_
#define MULTISCALAR_
#endif

#ifdef ROTBOUSS_SOL
#define DNS_
#define VELOC_
#define SCALAR_
#define ROTATION_
#define BOUSSINESQ_
#endif

#ifdef MPROTBOUSS_SOL
#define DNS_
#define VELOC_
#define SCALAR_
#define ROTATION_
#define BOUSSINESQ_
#define MULTISCALAR_
#endif

#ifdef GPE_SOL
#define DNS_
#define WAVEFUNCTION_
#endif

#ifdef ARGL_SOL
#define DNS_
#define NEWTON_
#define ADVECT_
#define QFORCE_
#define WAVEFUNCTION_
#endif

#ifdef LAHD_SOL
#define ALPHAV_
#define VELOC_
#endif

#ifdef CAHD_SOL
#define ALPHAV_
#define VELOC_
#endif

#ifdef LHD_SOL
#define ALPHAV_
#define VELOC_
#endif

#ifdef LAMHD_SOL
#define ALPHAV_
#define ALPHAB_
#define VELOC_
#define MAGFIELD_
#endif

#ifdef EDQNMHD_SOL
#define DNS_
#define EDQNM_
#define VELOC_
#endif

#ifdef EDQNMROTH_SOL
#define DNS_
#define EDQNM_
#define VELOC_
#define ROTATION_
#endif

#ifdef DEF_GHOST_LAGP
#define LAGPART_
#endif
