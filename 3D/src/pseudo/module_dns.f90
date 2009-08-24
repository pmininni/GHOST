!=================================================================
! MODULES to normalize fields in DNS
!
! 2007 Jonathan Pietarila Graham and Pablo D. Mininni
!      National Center for Atmospheric Research.
!=================================================================

!=================================================================

MODULE dns

CONTAINS
  SUBROUTINE normalize(fx,fy,fz,f0,kin,comm)
    USE grid
    USE mpivars
    IMPLICIT NONE

    DOUBLE PRECISION    :: tmp
    COMPLEX, INTENT(INOUT), DIMENSION(n,n,ista:iend) :: fx,fy,fz
    REAL, INTENT(IN)    :: f0
    INTEGER, INTENT(IN) :: kin
    INTEGER, INTENT(IN) :: comm
    
    CALL energy(fx,fy,fz,tmp,kin)
    CALL MPI_BCAST(tmp,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    fx = fx*f0/sqrt(tmp)
    fy = fy*f0/sqrt(tmp)
    fz = fz*f0/sqrt(tmp)
  END SUBROUTINE normalize

END MODULE dns
!=================================================================
