!=================================================================
! GHOST suite Fortran utilities 
!
! 2010 D. Rosenberg
!      NCAR
!
! 17 Nov 2010: Initial version
!=================================================================
MODULE gutils
! Utils data, if any:
!
! ...
! end, member data
!
!
! Methods:
      CONTAINS

      SUBROUTINE rarray_byte_swap(Rin, nin)
!-----------------------------------------------------------------
!
! Performs endian conversion on array of floats
!
! Parameters
!     Rin : artibrary-rank array whose values will be endian-swapped, returned
!     nin : dimension of Rin
!-----------------------------------------------------------------
      USE mpivars
      USE threads
      USE fprecision

      IMPLICIT NONE

      INTEGER, INTENT(IN)                        :: nin
      REAL(KIND=GP), INTENT(INOUT), DIMENSION(*) :: Rin

      INTEGER(KIND=GP) :: ie0, ie1
      INTEGER          :: i, j, k, m, nb

      nb = 8  ! no. bits per byte

      DO k = 1, nin
          ie0 = TRANSFER(Rin(k), 0_GP)
          DO m = 1, GP
             CALL MVBITS( ie0, (GP-m)*nb, nb, ie1, (m-1)*nb  )
          END DO
          Rin(k) = TRANSFER(ie1, 0.0_GP)
       END DO

      RETURN

      END SUBROUTINE rarray_byte_swap

!
!
      SUBROUTINE parseind(sind, sep, ind, nmax, nind) 
!-----------------------------------------------------------------
!
! Parses string of integers that are ';' separated, and stores, and
! stores them in the integer array, specifying how many integers
! were found in the string. A test is made against the specified
! max number of integer indices, nmax, before attempting to add
! any more indices to this array.
!
! Parameters
!     sind : ';'-separated string of integers (IN)
!     sep  : separator string (IN)
!     ind  : integer array contining integers found in sind (IN)
!     nmax : max size of 'ind' array (OUT)
!     nind : num. integers found in 'sind' (OUT)
!-----------------------------------------------------------------
      INTEGER, INTENT(OUT)         :: ind(*), nind
      INTEGER, INTENT(IN)          :: nmax
      CHARACTER(len=*), INTENT(IN) :: sind, sep

      INTEGER                      :: i
      CHARACTER(len=1024)          :: sint

      ib = 1;
      ie = len(sind)
      nind = 0
      DO WHILE ( len(trim(sind(ib:ie))) .GT. 0 )
        i = index(sind(ib:ie),sep)
        IF ( i .eq. 0 ) THEN
          sint = trim(adjustl(sind(ib:ie)))
          ib = ie + 1
        ELSE
          sint = trim(adjustl(sind(ib:(ib+i-2))))
          ib = ib + i
        ENDIF
        nind = nind + 1
        IF ( nind.GT.nmax ) RETURN
        READ(sint,'(I10)') ind(nind)
      ENDDO
      RETURN

      END SUBROUTINE parseind

END MODULE gutils
