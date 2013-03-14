!=================================================================
! GHOST particle database types (class)
!
! 2011 D. Rosenberg
!      NCAR
!
! 15 Aug 2011: Initial version
!=================================================================
MODULE pdbtypes
      USE fprecision
      IMPLICIT NONE

      TYPE, PUBLIC :: GPDBrec
        SEQUENCE
        INTEGER      : id
        REAL(KIND=GP): x
        REAL(KIND=GP): y
        REAL(KIND=GP): z
      END TYPE GPDBrec

!-----------------------------------------------------------------
!-----------------------------------------------------------------
END MODULE pdbtypes
