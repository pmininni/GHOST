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

      SUBROUTINE trunc(Cin, n, nt, kmax, Ctr) 
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!
! Truncates input spectral array, Cin, to maximum wave number, kmax,
! and stores in reduced-size array Ctr.
!
! Parameters
!     Cin : complex input array
!     n   : full grid dimensions. NOTE: these must be
!           consistent with (ista,jsta,ksta)-(iend,jend,kend)
!     nt  : reduced grid dimensions. NOTE: these must be
!           consistent with (itsta,jtsta,ktsta)-(itend,jtend,ktend)
!           indices in mpivars module
!     Ctr : truncated complex array on reduced grid
!-----------------------------------------------------------------
      USE grid
      USE mpivars
      USE kes

      IMPLICIT NONE
      INTEGER,INTENT   (IN)                                   :: n(3), nt(3)
      COMPLEX,INTENT(INOUT) ,DIMENSION(ista:iend,n(2),n(3))   :: Cin
      REAL,INTENT      (IN)                                   :: kmax
      COMPLEX,INTENT  (OUT),DIMENSION(itsta:itend,nt(2),nt(3)):: Ctr
      REAL(kind=GP)                                           :: fact
      INTEGER                                                 :: i, j, k
!
! Truncate in Fourier space:
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
            DO j = 1,n(2)
               DO k = 1,n(3)

                  IF (  kn2(k,j,i).GT.kmax ) Cin(k,j,i) = 0.0
               END DO
            END DO
         END DO
!
         fact = 1.0_GP / (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
            DO j = 1,n(2)
               DO k = 1,n(3)
                  Cin(k,j,i) = Cin(k,j,i)*fact
               END DO
            END DO
         END DO
!
! Store in coeffs in reduced storage for inverse transform:
         Ctr = 0.0
         DO i = itsta,itend
            DO j = 1,nt(2)
               DO k = 1,nt(3)
                  Ctr(k,j,i) = Cin(k,j,i)
               END DO
               DO k = nt(3)/2+1,nt(3)
                  Ctr(k,j,i) = Cin(k-nt(3)+n(3),j,i)
               END DO
            END DO
            DO j = nt(2)/2+1,nt(2)
               DO k = 1,nt(3)/2
                  Ctr(k,j,i) = Cin(k,j-nt(2)+n(2),i)
               END DO
               DO k = nt(3)/2+1,nt(3)
                  Ctr(k,j,i) = Cin(k-nt(3)+n(3),j-nt(2)+n(2),i)
               END DO
            END DO
         END DO
!
!!Compute inverse FT of truncated variable:
!!
!!       CALL trrange(1,nz    ,nzt    ,nprocs,myrank,ksta,kend)
!!       CALL trrange(1,nx/2+1,nxt/2+1,nprocs,myrank,ista,iend)
!!       IF ( myrank.eq.0 ) write(*,*)'main: complex_to_real...'
!!       CALL fftp3d_complex_to_real(plancrt,Ctr,tr,MPI_COMM_WORLD)
!!       IF ( myrank.eq.0 ) write(*,*)'main: complex_to_real done.'
!!       CALL range(1,nx/2+1,nprocs,myrank,ista,iend)
!!       CALL range(1,nz,    nprocs,myrank,ksta,kend)

      END SUBROUTINE trunc
!-----------------------------------------------------------------
!-----------------------------------------------------------------

END MODULE gtrunc

