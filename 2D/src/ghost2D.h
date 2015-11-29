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
#define STREAM_
#endif

#ifdef PHD_SOL
#define DNS_
#define STREAM_
#define SCALAR_
#endif

#ifdef MHD_SOL
#define DNS_
#define STREAM_
#define VECPOT_
#endif

#ifdef MHDB_SOL
#define DNS_
#define STREAM_
#define VECPOT_
#define UNIFORMB_
#endif

#ifdef HMHD_SOL
#define DNS_
#define D25_
#define STREAM_
#define VECPOT_
#define HALLTERM_
#endif

#ifdef SQG_SOL
#define DNS_
#define STREAM_
#endif

#ifdef SWHD_SOL
#define DNS_
#define VELOC_
#define SCALAR_
#define SW_
#endif
