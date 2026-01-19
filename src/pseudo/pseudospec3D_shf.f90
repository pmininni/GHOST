!=================================================================
! PSEUDOSPECTRAL subroutines
!
! Subroutines to compute structure functions in the SO(2) and 
! SO(3) decompositions in HD, MHD, and Hall-MHD simulations 
! with the GHOST code. This file provides functions to perform 
! shifts in the three spatial directions.
!
! NOTATION: index 'i' is 'x' 
!           index 'j' is 'y'
!           index 'k' is 'z'
!
! 2008 Luis Martin and Pablo D. Mininni.
!      Department of Physics, 
!      Facultad de Ciencias Exactas y Naturales.
!      Universidad de Buenos Aires.
!=================================================================

!*****************************************************************
      SUBROUTINE shiftz(a,dir,comm)
!-----------------------------------------------------------------
!
! Shifts a three-dimensional array one slice in the z-direction 
!
! Parameters
!     a    : input matrix (rewriten with the shifted array)
!     dir  : displacement (+1 or -1)
!     comm : the MPI communicator (handle)
!
      USE fprecision
      USE commtypes
      USE grid
      USE mpivars
      USE threads
      IMPLICIT NONE

      REAL(KIND=GP), INTENT(INOUT), DIMENSION(nx,ny,ksta:kend) :: a
      REAL(KIND=GP), DIMENSION(nx,ny) :: buffer1,buffer2
      INTEGER, DIMENSION(MPI_STATUS_SIZE) :: istatus
      INTEGER, INTENT(IN)  :: dir
      INTEGER, INTENT(IN)  :: comm
      INTEGER :: isendTo,igetFrom
      INTEGER :: ireq1,ireq2
      INTEGER :: zsend,zrecv
      INTEGER :: i,j,k

      igetFrom = myrank+dir
      isendTo = myrank-dir
      IF ( igetFrom .ge. nprocs ) igetFrom = igetFrom - nprocs
      IF ( isendTo .ge. nprocs ) isendTo = isendTo - nprocs
      IF ( igetFrom .lt. 0 ) igetFrom = igetFrom + nprocs
      IF ( isendTo .lt. 0 ) isendTo = isendTo + nprocs

      IF ( dir .eq. 1 ) THEN
         zsend = ksta
         zrecv = kend
      ELSE
         zsend = kend
         zrecv = ksta
      ENDIF
!$omp parallel do private (i)
      DO j = 1,ny
         DO i = 1,nx
            buffer2(i,j) = a(i,j,zsend)
         END DO
      END DO

      CALL MPI_IRECV(buffer1,nx*ny,GC_REAL,igetFrom,1,comm,ireq1,ierr)
      CALL MPI_ISEND(buffer2,nx*ny,GC_REAL,isendTo,1,comm,ireq2,ierr)

!$omp parallel do if ((zrecv-dir-zsend)/dir.ge.nth) private (j,i)
      DO k = zsend,zrecv-dir,dir
!$omp parallel do if ((zrecv-dir-zsend)/dir.ge.nth) private (i)
         DO j = 1,ny
            DO i = 1,nx
               a(i,j,k) = a(i,j,k+dir)
            END DO
         END DO
      END DO

      CALL MPI_WAIT(ireq1,istatus,ierr)
      CALL MPI_WAIT(ireq2,istatus,ierr)

!$omp parallel do private (i)
      DO j = 1,ny
         DO i = 1,nx
            a(i,j,zrecv) = buffer1(i,j)
         END DO
      END DO

      RETURN
      END SUBROUTINE shiftz

!*****************************************************************
      SUBROUTINE shiftx(a,dis)
!-----------------------------------------------------------------
!
! Shifts a three-dimensional array in the x-direction
!
! Parameters
!     a   : input matrix (rewriten with the shifted array)
!     dis : displacement in x
!
      USE fprecision
      USE grid
      USE mpivars
      USE threads
      IMPLICIT NONE

      REAL(KIND=GP), INTENT(INOUT), DIMENSION(nx,ny,ksta:kend) :: a
      INTEGER, INTENT(IN)                   :: dis
      REAL(KIND=GP), DIMENSION(abs(dis),ny,ksta:kend) :: buffer
      INTEGER :: d,s,inib,endb,alef,arig
      INTEGER :: i,j,k

      d = ABS(dis)
      s = dis/d
      IF ( s.gt.0 ) THEN
         inib = 0
         endb = nx-d
         alef = 0
         arig = d
      ELSE
         inib = nx-d
         endb = 0
         alef = d
         arig = 0
      ENDIF

!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
         DO j = 1,ny
            DO i = 1,d
               buffer(i,j,k) = a(i+inib,j,k)
            END DO
         END DO
      END DO
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
         DO j = 1,ny
            DO i = 1,nx-d,s
               a(i+alef,j,k) = a(i+arig,j,k)
            END DO
         END DO
      END DO
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
         DO j = 1,ny
            DO i = 1,d
               a(i+endb,j,k) = buffer(i,j,k)
            END DO
         END DO
      END DO

      RETURN
      END SUBROUTINE shiftx

!*****************************************************************
      SUBROUTINE shifty(a,dis)
!-----------------------------------------------------------------
!
! Shifts a three-dimensional array in the y-direction
!
! Parameters
!     a   : input matrix (rewriten with the shifted array)
!     dis : displacement in y
!
      USE fprecision
      USE grid
      USE mpivars
      USE threads
      IMPLICIT NONE

      REAL(KIND=GP), INTENT(INOUT), DIMENSION(nx,ny,ksta:kend) :: a
      INTEGER, INTENT(IN)                   :: dis
      REAL(KIND=GP), DIMENSION(nx,abs(dis),ksta:kend) :: buffer
      INTEGER :: d,s,inib,endb,alef,arig
      INTEGER :: i,j,k

      d = ABS(dis)
      s = dis/d
      IF ( s.gt.0 ) THEN
         inib = 0
         endb = ny-d
         alef = 0
         arig = d
      ELSE
         inib = ny-d
         endb = 0
         alef = d
         arig = 0
      ENDIF

!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
         DO j = 1,d
            DO i = 1,nx
               buffer(i,j,k) = a(i,j+inib,k)
            END DO
         END DO
      END DO
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
         DO j = 1,ny-d,s
            DO i = 1,nx
               a(i,j+alef,k) = a(i,j+arig,k)
            END DO
         END DO
      END DO
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
         DO j = 1,d
            DO i = 1,nx
               a(i,j+endb,k) = buffer(i,j,k)
            END DO
         END DO
      END DO

      RETURN
      END SUBROUTINE shifty

