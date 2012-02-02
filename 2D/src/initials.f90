! Initial condition for the passive/active scalar.
! The scalar can be passive (e.g., in the PHD solver) or active (as 
! in the shallow-water solvers, where it represents the height of 
! the column of fluid). 
! This file contains the expression used for the initial scalar. 
! You can use temporary real arrays R1 of size (n,jsta:jend), and 
! temporary complex arrays C1, C2, C12 and C13 of size (n,ista:iend)
! to do intermediate computations. 
! The variable c0 should control the global concentration of the 
! scalar, and variables cparam0-9 can be used to control the 
! amplitudes of individual terms. At the end, the height in 
! spectral space should be stored in the array th.

! Null initial scalar

    DO i=ista,iend
        DO j=1,n
            th(j,i) = 0.0_GP
        END DO
    END DO
