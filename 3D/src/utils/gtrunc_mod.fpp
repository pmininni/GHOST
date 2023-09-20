!=================================================================
! GHOST suite Fortran utilities 
!
! 2023 D. Rosenberg
!      CIRA
!
! 13 Sep 2023: Initial version
!=================================================================
MODULE gtrunc
!
      USE fprecision
      USE fftplans
      USE iovar
!
! Utils data, if any:
!
     ! Input config data for call to bouss_lescomp method
     TYPE TRUNCDAT
         TYPE(IOPLAN)    :: planiot
         TYPE(FFTPLAN)   :: plancrt
         INTEGER         :: commtrunc
         COMPLEX(KIND=GP), DIMENSION (:,:,:), POINTER :: C1,C2,C3,C4,C5
         COMPLEX(KIND=GP), DIMENSION (:,:,:), POINTER :: CT1
         REAL   (KIND=GP), DIMENSION (:,:,:), POINTER :: R1,R2,R3
         REAL   (KIND=GP), DIMENSION (:,:,:), POINTER :: RT1
         REAL(KIND=GP)   :: ktrunc = 1.0_GP / 9.0_GP
         CHARACTER(LEN=1024) &
                         :: odir
      END TYPE TRUNCDAT

! end, member data
!
!
! Methods:
      CONTAINS

      SUBROUTINE trunc(Cin, n, nt, kmax, inorm, Ctmp, Ctr) 
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!
! Truncates input spectral array, Cin, to maximum wave number, kmax,
! and stores in reduced-size array Ctr. Draws starting/ending indices,
! ista/iend, itsta/itend, ksta/kend, ktsta/ktend from mod files.
!
! Parameters
!     Cin  : complex input array, on original grid.
!     n    : full grid dimensions. NOTE: these must be
!            consistent with (ista,jsta,ksta)-(iend,jend,kend)
!     nt   : reduced grid dimensions. NOTE: these must be
!            consistent with (itsta,jtsta,ktsta)-(itend,jtend,ktend)
!            indices in mpivars module
!     kmax : maximun kn2
!     inorm: normalize (1), or not (0)
!     Ctmp : complex tmp array on original grid
!     Ctr  : truncated complex array on reduced grid
!-----------------------------------------------------------------
      USE grid
      USE mpivars
      USE kes

      IMPLICIT NONE
      INTEGER,INTENT   (IN)                                   :: inorm, n(3), nt(3)
      COMPLEX,INTENT   (IN) ,DIMENSION(n(3),n(2),ista:iend)   :: Cin
      COMPLEX,INTENT(INOUT) ,DIMENSION(n(3),n(2),ista:iend)   :: Ctmp
      REAL,INTENT      (IN)                                   :: kmax
      COMPLEX,INTENT  (OUT),DIMENSION(nt(3),nt(2),itsta:itend):: Ctr

      REAL(kind=GP)                                           :: fact
      INTEGER                                                 :: i, j, k
!
! Truncate in Fourier space:
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
            DO j = 1,n(2)
               DO k = 1,n(3)
                  Ctmp(k,j,i) = Cin(k,j,i)
                  IF (  kn2(k,j,i).GT.kmax ) Ctmp(k,j,i) = 0.0
               END DO
            END DO
         END DO
!
         fact = 1.0_GP / (real(n(1),kind=GP)*real(n(2),kind=GP)*real(n(3),kind=GP))
         if ( inorm .LE. 0 ) fact = 1.0_GP
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
            DO j = 1,n(2)
               DO k = 1,n(3)
                  Ctmp(k,j,i) = Ctmp(k,j,i)*fact
               END DO
            END DO
         END DO
!
! Store in coeffs in reduced storage for inverse transform:
         Ctr = 0.0
         DO i = itsta,itend
            DO j = 1,nt(2)/2
               DO k = 1,nt(3)/2
                  Ctr(k,j,i) = Ctmp(k,j,i)
               END DO
               DO k = nt(3)/2+1,nt(3)
                  Ctr(k,j,i) = Ctmp(k-nt(3)+n(3),j,i)
               END DO
            END DO
            DO j = nt(2)/2+1,nt(2)
               DO k = 1,nt(3)/2
                  Ctr(k,j,i) = Ctmp(k,j-nt(2)+n(2),i)
               END DO
               DO k = nt(3)/2+1,nt(3)
                  Ctr(k,j,i) = Ctmp(k-nt(3)+n(3),j-nt(2)+n(2),i)
               END DO
            END DO
         END DO
!

      END SUBROUTINE trunc
!-----------------------------------------------------------------
!-----------------------------------------------------------------


      SUBROUTINE trunc_dd(Cin, n,ind, nt, indt, kmax, inorm, Ctmp, Ctr) 
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!
! Truncates input spectral array, Cin, to maximum wave number, kmax,
! and stores in reduced-size array Ctr. Uses specified domain
! decomposition (dd) indices, instead of drawing them from the mod files
!
! Parameters
!     Cin  : complex input array, on original grid.
!     n    : full grid dimensions. NOTE: these must be
!            consistent with (ista,jsta,ksta)-(iend,jend,kend)
!     ind  : array [ista,iend] corresponding to n, and
!            that define the domain decomposition of Cin
!     nt   : reduced grid dimensions. NOTE: these must be
!            consistent with (itsta,jtsta,ktsta)-(itend,jtend,ktend)
!            indices in mpivars module
!     indt : array [ista,iend] corresponding to nt, and
!            that define the domain decomposition of Ctr 
!     kmax : maximun kn2
!     inorm: normalize (1), or not (0)
!     Ctmp : complex tmp array on original grid
!     Ctr  : truncated complex array on reduced grid
!-----------------------------------------------------------------
      USE grid
!     USE mpivars
      USE kes

      IMPLICIT NONE
      INTEGER,INTENT   (IN)                                       :: inorm, ind(2), indt(2)
      INTEGER,INTENT   (IN)                                       :: n(3), nt(3)
      COMPLEX,INTENT   (IN) ,DIMENSION(n(3),n(2),ind(1):ind(2))   :: Cin
      COMPLEX,INTENT(INOUT) ,DIMENSION(n(3),n(2),ind(1):ind(2))   :: Ctmp
      REAL,INTENT      (IN)                                       :: kmax
      COMPLEX,INTENT  (OUT),DIMENSION(nt(3),nt(2),indt(1):indt(2)):: Ctr

      REAL(kind=GP)                                               :: fact
      INTEGER                                                     :: i, j, k
!
! Truncate in Fourier space:
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ind(1), ind(2)
!$omp parallel do if (iend-ista.lt.nth) private (k)
            DO j = 1,n(2)
               DO k = 1,n(3)
                  Ctmp(k,j,i) = Cin(k,j,i)
                  IF (  kn2(k,j,i).GT.kmax ) Ctmp(k,j,i) = 0.0
               END DO
            END DO
         END DO
!
         fact = 1.0_GP / (real(n(1),kind=GP)*real(n(2),kind=GP)*real(n(3),kind=GP))
         if ( inorm .LE. 0 ) fact = 1.0_GP
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ind(1), ind(2)
!$omp parallel do if (iend-ista.lt.nth) private (k)
            DO j = 1,n(2)
               DO k = 1,n(3)
                  Ctmp(k,j,i) = Ctmp(k,j,i)*fact
               END DO
            END DO
         END DO
!
! Store in coeffs in reduced storage for inverse transform:
         Ctr = 0.0
         DO i = indt(1), indt(2)
            DO j = 1,nt(2)/2
               DO k = 1,nt(3)/2
                  Ctr(k,j,i) = Ctmp(k,j,i)
               END DO
               DO k = nt(3)/2+1,nt(3)
                  Ctr(k,j,i) = Ctmp(k-nt(3)+n(3),j,i)
               END DO
            END DO
            DO j = nt(2)/2+1,nt(2)
               DO k = 1,nt(3)/2
                  Ctr(k,j,i) = Ctmp(k,j-nt(2)+n(2),i)
               END DO
               DO k = nt(3)/2+1,nt(3)
                  Ctr(k,j,i) = Ctmp(k-nt(3)+n(3),j-nt(2)+n(2),i)
               END DO
            END DO
         END DO
!

      END SUBROUTINE trunc_dd
!-----------------------------------------------------------------
!-----------------------------------------------------------------

END MODULE gtrunc

