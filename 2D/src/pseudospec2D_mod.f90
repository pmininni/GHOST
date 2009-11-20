!=================================================================
! MODULES for 2D codes
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
      INTEGER :: n = 256
      SAVE

  END MODULE grid
!=================================================================

  MODULE fft
!
      USE fftplans
      TYPE (FFTPLAN) :: planrc, plancr
      SAVE

  END MODULE fft
!=================================================================

  MODULE ali
      REAL :: kmax
      REAL :: tiny
      SAVE

  END MODULE ali
!=================================================================

  MODULE var
      REAL    :: pi = 3.14159265
      COMPLEX :: im = (0.,1.)
      SAVE

  END MODULE var
!=================================================================

  MODULE kes
      REAL, ALLOCATABLE, DIMENSION (:)   :: ka
      REAL, ALLOCATABLE, DIMENSION (:,:) :: ka2
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
!
       REAL FUNCTION randn(idum, sig)
!
! Normally-distributed random number, given variance, sig.
! Uses randu. Taken from Numerical Recipies, and modified to 
! include variance.
       INTEGER:: iset
       REAL   :: x, y, rsq, rkeep, ssq, w
       SAVE   :: iset
       
       ssq = sig*sig 
       IF ( idum.lt.0 ) iset = 0
       IF ( iset .eq. 0 ) THEN
         rsq = 0.0
         w   = 10.0 
         DO WHILE ( rsq .ge. 1.0 .or. rsq .eq. 0.0 )
           x = 2.0*ssq* randu(idum) - 1.0
           y = 2.0*ssq* randu(idum) - 1.0
           rsq = x*x + y*y
         ENDDO
         w = sqrt(-2.0*log(rsq)/rsq)
         rkeep = w * x
         iset = 1
         randn = y * w
       ELSE
         iset = 0
         randn = rkeep
       ENDIF

       END FUNCTION randn

  END MODULE random
!=================================================================
