!=================================================================
! MODULES for 3D codes
!
! 2003 Pablo D. Mininni.
!      Department of Physics, 
!      Facultad de Ciencias Exactas y Naturales.
!      Universidad de Buenos Aires.
!      e-mail: mininni@df.uba.ar 
!=================================================================

!=================================================================

  MODULE grid
!
! n: number of points in the spatial grid
      INTEGER :: n = N_
      SAVE

  END MODULE grid
!=================================================================

  MODULE order
!
! ord: number of iterations in the Runge-Kutta method
      INTEGER :: ord = ORD_
      SAVE

  END MODULE order
!=================================================================

  MODULE filefmt
!
! Change the length of the string 'ext' to change the number 
! of characters used to number the binary files and files with 
! the spectra. The format fmtext should be consistent with the 
! length of the string, e.g. if len=5 then fmtext = '(i5.5)'.
      CHARACTER(len=3) :: ext
      CHARACTER(len=6),SAVE :: fmtext = '(i3.3)'

  END MODULE filefmt
!=================================================================

  MODULE fft
!
      USE fftplans
      TYPE(FFTPLAN) :: planrc, plancr
      SAVE

  END MODULE fft
!=================================================================

  MODULE ali
      REAL :: kmax
      REAL :: tiny
      REAL :: tinyf
      SAVE

  END MODULE ali
!=================================================================

  MODULE var
      REAL    :: pi = 3.14159265
      COMPLEX :: im = (0.,1.)
      SAVE

  END MODULE var
!=================================================================

  MODULE hall
      REAL    :: ep
      INTEGER :: gspe
      SAVE

  END MODULE hall
!=================================================================

  MODULE kes
      REAL, ALLOCATABLE, DIMENSION (:)     :: ka
      REAL, ALLOCATABLE, DIMENSION (:,:,:) :: ka2
      SAVE

  END MODULE kes
!=================================================================

  MODULE random
      CONTAINS
       REAL FUNCTION randu(idum)
!
! Uniform distributed random numbers between -1 and 
! 1. The seed idum must be between 0 and the value 
! of mask

       INTEGER, PARAMETER :: iq=127773,ir=2836,mask=123459876
       INTEGER, PARAMETER :: ia=16807,im=2147483647
       INTEGER            :: k,idum
       REAL, PARAMETER    :: am=1./im

       idum = ieor(idum,mask)
       k = idum/iq
       idum = ia*(idum-k*iq)-ir*k
       IF (idum.lt.0) idum = idum+im
       randu = am*idum
       randu = (randu-.5)*2
       idum = ieor(idum,mask)
       END FUNCTION randu

       REAL FUNCTION randn(idum)
!
! Normally distributed random numbers with zero mean 
! and unit variance. The seed idum must be between 0 
! and the value of mask in randu.

       REAL          :: v1,v2,ran1
       REAL          :: fac,rsq
       REAL,SAVE     :: gset
       INTEGER       :: idum
       INTEGER, SAVE :: iset

       IF ((iset.ne.0).or.(iset.ne.1)) iset=0
       IF (idum.lt.0) iset=0
       IF (iset.eq.0) THEN
          rsq = 2.
          DO WHILE ((rsq.ge.1.).or.(rsq.eq.0.))
             v1 = randu(idum)
             v2 = randu(idum)
             rsq = v1**2+v2**2
          END DO
          fac = sqrt(-2.*log(rsq)/rsq)
          gset = v1*fac
          randn = v2*fac
          iset = 1
       ELSE
          randn = gset
          iset = 0
       ENDIF
       END FUNCTION randn

  END MODULE random
!=================================================================
