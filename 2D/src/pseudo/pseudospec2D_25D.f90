!=================================================================
! PSEUDOSPECTRAL subroutines
!
! Subroutines for computing spatial derivatives and nonlinear 
! terms in incompressible equations in 2.5D using a 
! pseudo-spectral method. You should use the FFTPLANS and 
! MPIVARS modules (see the file 'fftp2D_mod.f90') in each 
! program that call any of the subroutines in this file. 
!
! NOTATION: index 'i' is 'x' 
!           index 'j' is 'y'
!
! 2003 Pablo D. Mininni.
!      Department of Physics, 
!      Facultad de Ciencias Exactas y Naturales.
!      Universidad de Buenos Aires.
!      e-mail: mininni@df.uba.ar 
!=================================================================

!*****************************************************************
      SUBROUTINE energy25(a,b,c,kin)
!-----------------------------------------------------------------
!
! Computes the mean kinetic or magnetic energy in 2.5D,
! and the mean square current density or vorticity. 
! The output is valid only in the first node.
!
! Parameters
!     a  : input matrix with the scalar field
!     b  : input matrix with the z component of the field
!     c  : at the output contains the energy
!     kin: =1 computes the energy
!          =0 computes the current or vorticity
!
      USE fprecision
      USE commtypes
      USE mpivars
      USE grid
      USE kes
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,ista:iend) :: a,b
      DOUBLE PRECISION, INTENT(OUT) :: c
      DOUBLE PRECISION              :: cloc
      REAL(KIND=GP)                 :: tmp
      INTEGER, INTENT(IN) :: kin
      INTEGER             :: i,j

      cloc = 0.0_GP
      tmp = 1.0_GP/real(n,kind=GP)**4
!
! Computes the energy
!
      IF (kin.eq.1) THEN
         IF (ista.eq.1) THEN
            DO j = 1,n
               cloc = cloc+(ka2(j,1)*abs(a(j,1))**2+abs(b(j,1))**2)*tmp
            END DO
            DO i = 2,iend
               DO j = 1,n
                  cloc = cloc+2*(ka2(j,i)*abs(a(j,i))**2+abs(b(j,i))**2)*tmp
               END DO
            END DO
          ELSE
            DO i = ista,iend
               DO j = 1,n
                  cloc = cloc+2*(ka2(j,i)*abs(a(j,i))**2+abs(b(j,i))**2)*tmp
               END DO
            END DO
          ENDIF
!
! Computes the current or vorticity
!
      ELSE
         IF (ista.eq.1) THEN
            DO j = 1,n
               cloc = cloc+ka2(j,1)*(ka2(j,1)*abs(a(j,1))**2+ &
                      abs(b(j,1))**2)*tmp
            END DO
            DO i = 2,iend
               DO j = 1,n
                  cloc = cloc+2*ka2(j,i)*(ka2(j,i)*abs(a(j,i))**2+ &
                         abs(b(j,i))**2)*tmp
               END DO
            END DO
          ELSE
            DO i = ista,iend
               DO j = 1,n
                  cloc = cloc+2*ka2(j,i)*(ka2(j,i)*abs(a(j,i))**2+ &
                         abs(b(j,i))**2)*tmp
               END DO
            END DO
          ENDIF
      ENDIF
!
! Computes the reduction between nodes
!
      CALL MPI_REDUCE(cloc,c,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)

      RETURN
      END SUBROUTINE energy25

!*****************************************************************
      SUBROUTINE helicity25(a,b,c)
!-----------------------------------------------------------------
!
! Computes the mean kinetic or magnetic helicity in 2.5D.
! The output is only valid in the first node.
!
! Parameters
!     a  : input matrix with the scalar field
!     b  : input matrix with the z component of the field
!     c  : at the output contains the helicity
!
      USE fprecision
      USE commtypes
      USE mpivars
      USE grid
      USE kes
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,ista:iend) :: a,b
      DOUBLE PRECISION, INTENT(OUT) :: c
      DOUBLE PRECISION              :: cloc
      REAL(KIND=GP)  :: tmp
      INTEGER        :: i,j

      cloc = 0.0_GP
      tmp = 1.0_GP/real(n,kind=GP)**4

      IF (ista.eq.1) THEN
         DO j = 1,n
            cloc = cloc+2*ka2(j,1)*real(b(j,1)*conjg(a(j,1)))*tmp
         END DO
         DO i = 2,iend
            DO j = 1,n
               cloc = cloc+4*ka2(j,i)*real(b(j,i)*conjg(a(j,i)))*tmp
            END DO
         END DO
      ELSE
         DO i = ista,iend
            DO j = 1,n
               cloc = cloc+4*ka2(j,i)*real(b(j,i)*conjg(a(j,i)))*tmp
            END DO
         END DO
      ENDIF
      CALL MPI_REDUCE(cloc,c,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)

      RETURN
      END SUBROUTINE helicity25

!*****************************************************************
      SUBROUTINE mhdcheck25(a,b,c,d,e,f,g,h,t,hel)
!-----------------------------------------------------------------
!
! Consistency check for the conservation of energy in 2.5D
!
! Parameters
!     a  : streamfunction
!     b  : vector potential
!     c  : velocity field in the z-direction
!     d  : magnetic field in the z-direction
!     e  : external mechanic force
!     f  : external magnetic force
!     g  : external mechanic force in the z-direction
!     h  : external magnetic force in the z-direction
!     t  : time
!     hel: =1 computes the magnetic helicity
!          =0 skips helicity computation
!
      USE fprecision
      USE commtypes
      USE mpivars
      USE grid
      USE kes
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,ista:iend) :: a,b,c,d
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,ista:iend) :: e,f,g,h
      DOUBLE PRECISION :: engk,engm,eng,tmp
      DOUBLE PRECISION :: potk,potm,ens,cur
      DOUBLE PRECISION :: helk,helm
      REAL(KIND=GP)    :: t
      REAL(KIND=GP)    :: tmq
      INTEGER, INTENT(IN) :: hel
      INTEGER             :: i,j

      tmq = 1.0_GP/real(n,kind=GP)**4
!
! Computes the mean energy, enstrophy, square current
! and kinetic and magnetic helicity
!
      CALL energy25(a,c,engk,1)
      CALL energy25(b,d,engm,1)
      CALL energy25(a,c,ens,0)
      CALL energy25(b,d,cur,0)
      eng = engk+engm
      IF (hel.eq.1) THEN
         CALL helicity25(a,c,helk)
         CALL helicity25(b,d,helm)
      ENDIF
!
! Computes the  kinetic energy injection rate
!
      IF (ista.eq.1) THEN
         DO j = 1,n
            tmp = tmp+(ka2(j,1)*(real(e(j,1)*conjg(a(j,1)))      &
                  +real(f(j,1)*conjg(b(j,1)))))*tmq
         END DO
         DO i = 2,iend
            DO j = 1,n
               tmp = tmp+2*(ka2(j,i)*(real(e(j,i)*conjg(a(j,i))) &
                     +real(f(j,i)*conjg(b(j,i)))))*tmq
            END DO
         END DO
      ELSE
         DO i = ista,iend
            DO j = 1,n
               tmp = tmp+2*(ka2(j,i)*(real(e(j,i)*conjg(a(j,i))) &
                     +real(f(j,i)*conjg(b(j,i)))))*tmq
            END DO
         END DO
      ENDIF
      CALL MPI_REDUCE(tmp,potk,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
!
! Computes the magnetic energy injection rate
!
      IF (ista.eq.1) THEN
         DO j = 1,n
            tmp = tmp+(ka2(j,1)*(real(g(j,1)*conjg(c(j,1)))      &
                  +real(h(j,1)*conjg(d(j,1)))))*tmq
         END DO
         DO i = 2,iend
            DO j = 1,n
               tmp = tmp+2*(ka2(j,i)*(real(g(j,i)*conjg(c(j,i))) &
                     +real(h(j,i)*conjg(d(j,i)))))*tmq
            END DO
         END DO
      ELSE
         DO i = ista,iend
            DO j = 1,n
               tmp = tmp+2*(ka2(j,i)*(real(g(j,i)*conjg(c(j,i))) &
                     +real(h(j,i)*conjg(d(j,i)))))*tmq
            END DO
         END DO
      ENDIF
      CALL MPI_REDUCE(tmp,potm,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
!
! Creates external files to store the results
!
      IF (myrank.eq.0) THEN
         OPEN(1,file='balance.txt',position='append')
         WRITE(1,10) t,eng,ens,cur
   10    FORMAT( E13.6,E22.14,E22.14,E22.14 )
         CLOSE(1)
         OPEN(1,file='energy.txt',position='append')
         WRITE(1,20) t,engk,engm
   20    FORMAT( E13.6,E22.14,E22.14 )
         CLOSE(1)
         OPEN(1,file='inject.txt',position='append')
         WRITE(1,20) t,potk,potm
         CLOSE(1)
         IF (hel.eq.1) THEN
            OPEN(1,file='helicity.txt',position='append')
            WRITE(1,20) t,helk,helm
            CLOSE(1)
         ENDIF
      ENDIF

      RETURN
      END SUBROUTINE mhdcheck25

!*****************************************************************
      SUBROUTINE spectrum25(a,b,nmb,kin)
!-----------------------------------------------------------------
!
! Computes the energy power spectrum in 2.5D. 
! The output is written to a file by the first node.
!
! Parameters
!     a  : streamfunction or vector potential
!     c  : field component in the z-direction
!     ext: the extension used when writting the file
!     kin: =1 computes the kinetic spectrum
!          =0 computes the magnetic spectrum
!
      USE fprecision
      USE commtypes
      USE mpivars
      USE grid
      USE kes
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,ista:iend) :: a,b
      DOUBLE PRECISION, DIMENSION(n/2+1)  :: Ekp,Ekptot
      DOUBLE PRECISION, DIMENSION(n/2+1)  :: Ekn,Ekntot
      REAL(KIND=GP)       :: tmp
      INTEGER, INTENT(IN) :: kin
      INTEGER             :: kmn
      INTEGER             :: i,j
      CHARACTER(len=*), INTENT(IN) :: nmb

!
! Sets Ekp and Ekn to zero
!
      DO i = 1,n/2+1
         Ekp(i) = 0.0_GP
         Ekn(i) = 0.0_GP
      END DO
!
! Computes the energy spectrum
!
      tmp = 1./real(n,kind=GP)**4
      IF (ista.eq.1) THEN
         DO j = 1,n
            kmn = int(sqrt(ka2(j,1))+.501)
            IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
               Ekp(kmn) = Ekp(kmn)+ka2(j,1)*abs(a(j,1))**2*tmp
               Ekn(kmn) = Ekn(kmn)+abs(b(j,1))**2*tmp
            ENDIF
         END DO
         DO i = 2,iend
            DO j = 1,n
               kmn = int(sqrt(ka2(j,i))+.501)
               IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                  Ekp(kmn) = Ekp(kmn)+2*ka2(j,i)*abs(a(j,i))**2*tmp
                  Ekn(kmn) = Ekn(kmn)+2*abs(b(j,i))**2*tmp
               ENDIF
            END DO
         END DO
      ELSE
         DO i = ista,iend
            DO j = 1,n
               kmn = int(sqrt(ka2(j,i))+.501)
               IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                  Ekp(kmn) = Ekp(kmn)+2*ka2(j,i)*abs(a(j,i))**2*tmp
                  Ekn(kmn) = Ekn(kmn)+2*abs(b(j,i))**2*tmp
               ENDIF
            END DO
         END DO
      ENDIF
!
! Computes the reduction between nodes
! and exports the result to a file
!
      CALL MPI_REDUCE(Ekp,Ekptot,n/2+1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(Ekn,Ekntot,n/2+1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
      IF (myrank.eq.0) THEN
         IF (kin.eq.1) THEN
            OPEN(1,file='psspectrum.' // nmb // '.txt')
         ELSE
            OPEN(1,file='azspectrum.' // nmb // '.txt')
         ENDIF
         WRITE(1,30) Ekptot
   30    FORMAT( E23.15 ) 
         CLOSE(1)
      ENDIF
      IF (myrank.eq.0) THEN
         IF (kin.eq.1) THEN
            OPEN(1,file='vzspectrum.' // nmb // '.txt')
         ELSE
            OPEN(1,file='bzspectrum.' // nmb // '.txt')
         ENDIF
         WRITE(1,30) Ekntot
         CLOSE(1)
      ENDIF

      RETURN
      END SUBROUTINE spectrum25
