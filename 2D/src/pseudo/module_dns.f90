!=================================================================
! MODULES to normalize fields in DNS
!
! 2007 Jonathan Pietarila Graham and Pablo D. Mininni
!      National Center for Atmospheric Research.
!=================================================================

!=================================================================

MODULE dns

CONTAINS
  SUBROUTINE normalize(fs,f0,kin,comm)
    USE fprecision
    USE commtypes
    USE grid
    USE mpivars
!$  USE threads
    IMPLICIT NONE

    DOUBLE PRECISION    :: tmp
    COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(n,ista:iend) :: fs
    REAL(KIND=GP), INTENT(IN) :: f0
    INTEGER, INTENT(IN) :: kin
    INTEGER, INTENT(IN) :: comm
    INTEGER :: i,j

    CALL energy(fs,tmp,kin)
    CALL MPI_BCAST(tmp,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    DO i = ista,iend
       DO j = 1,n
          fs(j,i) = fs(j,i)*f0/sqrt(tmp)
       END DO
    END DO
  END SUBROUTINE normalize

END MODULE dns
!=================================================================
