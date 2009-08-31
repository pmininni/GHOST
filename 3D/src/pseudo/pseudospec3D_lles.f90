!=================================================================
! PSEUDOSPECTRAL subroutines
!
! Extra subroutines to compute spatial nonlinear terms common 
! to several Lagrangian averaged subgrid models. You should use 
! the FFTPLANS and MPIVARS modules (see the file 'fftp_mod.f90') 
! in each program that call any of the subroutines in this file. 
!
! NOTATION: index 'i' is 'x' 
!           index 'j' is 'y'
!           index 'k' is 'z'
!
! 2007 Jonathan Pietarila Graham and Pablo D. Mininni
!      National Center for Atmospheric Research.
!=================================================================

!*****************************************************************
      SUBROUTINE avector3(a,b,c,d,e,f,x,y,z,alpk,alpm)
!-----------------------------------------------------------------
!
! Computes the product AsxBs in real space. The components 
! of the vector fields A and B are given by the matrixes 
! a,b,c,d,e and f, following the right hand convention.
!
! Parameters
!     a   : input matrix with A_x
!     b   : input matrix with A_y
!     c   : input matrix with A_z
!     d   : input matrix with B_x
!     e   : input matrix with B_y
!     f   : input matrix with B_z
!     x   : at the output contains (AsxBs)_x in Fourier space
!     y   : at the output contains (AsxBs)_y in Fourier space
!     z   : at the output contains (AsxBs)_z in Fourier space
!     alpk: value of alpha for the field A
!     alpk: value of alpha for the field B
!
      USE fprecision
      USE commtypes
      USE mpivars
      USE grid
      USE fft
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,n,ista:iend)  :: a,b,c
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,n,ista:iend)  :: d,e,f
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(n,n,ista:iend) :: x,y,z
      REAL(KIND=GP), DIMENSION(n,n,ksta:kend) :: r1,r2
      REAL(KIND=GP), DIMENSION(n,n,ksta:kend) :: r3,r4
      REAL(KIND=GP), DIMENSION(n,n,ksta:kend) :: r5,r6
      REAL(KIND=GP), DIMENSION(n,n,ksta:kend) :: r7
      REAL(KIND=GP), INTENT(IN) :: alpk,alpm
      REAL(KIND=GP)             :: tmp
      INTEGER          :: i,j,k

      CALL smooth3(a,b,c,x,y,z,alpk)
      CALL fftp3d_complex_to_real(plancr,x,r1,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,y,r2,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,z,r3,MPI_COMM_WORLD)
      CALL smooth3(d,e,f,x,y,z,alpm)
      CALL fftp3d_complex_to_real(plancr,x,r4,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,y,r5,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,z,r6,MPI_COMM_WORLD)

      tmp = 1./real(n,kind=GP)**6
      DO k = ksta,kend
         DO j = 1,n
            DO i = 1,n
               r7(i,j,k) = (r2(i,j,k)*r6(i,j,k)-r5(i,j,k)*r3(i,j,k))*tmp
               r3(i,j,k) = (r3(i,j,k)*r4(i,j,k)-r6(i,j,k)*r1(i,j,k))*tmp
               r1(i,j,k) = (r1(i,j,k)*r5(i,j,k)-r4(i,j,k)*r2(i,j,k))*tmp
            END DO
         END DO
      END DO

      CALL fftp3d_real_to_complex(planrc,r7,x,MPI_COMM_WORLD)
      CALL fftp3d_real_to_complex(planrc,r3,y,MPI_COMM_WORLD)
      CALL fftp3d_real_to_complex(planrc,r1,z,MPI_COMM_WORLD)

      RETURN
      END SUBROUTINE avector3

!*****************************************************************
      SUBROUTINE aenergy(a,b,c,d,alp,kin)
!-----------------------------------------------------------------
!
! Computes the mean alpha-energy of a vector 
! field. The output is valid only in the first 
! node.
!
! Parameters
!     a  : input matrix in the x-direction
!     b  : input matrix in the y-direction
!     c  : input matrix in the z-direction
!     d  : at the output contains the energy
!     alp: value of alpha
!     kin: =0 computes the magnetic energy
!          =1 computes the kinetic energy

      USE fprecision
      USE commtypes
      USE kes
      USE grid
      USE mpivars
      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(OUT) :: d
      DOUBLE PRECISION              :: dloc
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,n,ista:iend) :: a,b,c
      COMPLEX(KIND=GP), DIMENSION(n,n,ista:iend)             :: c1,c2,c3
      REAL(KIND=GP), INTENT(IN)    :: alp
      INTEGER, INTENT(IN) :: kin
      INTEGER             :: i,j,k

      dloc = 0.

!
! Computes the kinetic energy
!
      IF (kin.eq.1) THEN
         IF (ista.eq.1) THEN
            DO j = 1,n
               DO k = 1,n
                  dloc = dloc+(abs(a(k,j,1))**2+abs(b(k,j,1))**2+        &
                         abs(c(k,j,1))**2)/((1.+alp**2*ka2(k,j,1))*      &
                         real(n,kind=GP)**6)
               END DO
            END DO
            DO i = 2,iend
               DO j = 1,n
                  DO k = 1,n
                     dloc = dloc+2*(abs(a(k,j,i))**2+abs(b(k,j,i))**2+   &
                            abs(c(k,j,i))**2)/((1.+alp**2*ka2(k,j,i))*   &
                            real(n,kind=GP)**6)
                  END DO
               END DO
            END DO
          ELSE
            DO i = ista,iend
               DO j = 1,n
                  DO k = 1,n
                     dloc = dloc+2*(abs(a(k,j,i))**2+abs(b(k,j,i))**2+   &
                            abs(c(k,j,i))**2)/((1.+alp**2*ka2(k,j,i))*   &
                            real(n,kind=GP)**6)
                  END DO
               END DO
            END DO
          ENDIF
!
! Computes the magnetic energy
!
      ELSE
         CALL rotor3(b,c,c1,1)
         CALL rotor3(a,c,c2,2)
         CALL rotor3(a,b,c3,3)
         IF (ista.eq.1) THEN
            DO j = 1,n
               DO k = 1,n
                  dloc = dloc+(abs(c1(k,j,1))**2+abs(c2(k,j,1))**2+      &
                         abs(c3(k,j,1))**2)/((1.+alp**2*ka2(k,j,1))*     &
                         real(n,kind=GP)**6)
               END DO
            END DO
            DO i = 2,iend
               DO j = 1,n
                  DO k = 1,n
                     dloc = dloc+2*(abs(c1(k,j,i))**2+abs(c2(k,j,i))**2+ &
                            abs(c3(k,j,i))**2)/((1.+alp**2*ka2(k,j,i))*  &
                            real(n,kind=GP)**6)
                  END DO
               END DO
            END DO
         ELSE
            DO i = ista,iend
               DO j = 1,n
                  DO k = 1,n
                     dloc = dloc+2*(abs(c1(k,j,i))**2+abs(c2(k,j,i))**2+ &
                            abs(c3(k,j,i))**2)/((1.+alp**2*ka2(k,j,i))*  &
                            real(n,kind=GP)**6)
                  END DO
               END DO
            END DO
         ENDIF
      ENDIF
!
! Computes the reduction between nodes
!
      CALL MPI_REDUCE(dloc,d,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)

      RETURN
      END SUBROUTINE aenergy

!*****************************************************************
      SUBROUTINE across(a,b,c,d,e,f,g,alp,kin)
!-----------------------------------------------------------------
!
! Computes the alpha-inner product of two vector 
! fields. The output is valid only in the first 
! node.
!
! Parameters
!     a  : first field x-component
!     b  : first field y-component
!     c  : first field z-component
!     d  : second field x-component
!     e  : second field y-component
!     f  : second field z-component
!     g  : at the output contains the inner product
!     alp: value of alpha
!     kin: =0 computes the inner product of the curls
!          =1 computes the inner product of the fields
!
      USE fprecision
      USE commtypes
      USE kes
      USE grid
      USE mpivars
      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(OUT) :: g
      DOUBLE PRECISION              :: gloc
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,n,ista:iend) :: a,b,c
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,n,ista:iend) :: d,e,f
      REAL(KIND=GP), INTENT(IN)    :: alp
      INTEGER, INTENT(IN) :: kin
      INTEGER             :: i,j,k

      gloc = 0.
!
! Computes the averaged inner product between the fields
!
      IF (kin.eq.1) THEN
         IF (ista.eq.1) THEN
            DO j = 1,n
               DO k = 1,n
                  gloc = gloc+real(a(k,j,1)*conjg(d(k,j,1))+                 &
                        b(k,j,1)*conjg(e(k,j,1))+c(k,j,1)*                   &
                        conjg(f(k,j,1)))/((1.+alp**2*ka2(k,j,1))*            &
                        real(n)**6)
               END DO
            END DO
            DO i = 2,iend
               DO j = 1,n
                  DO k = 1,n
                     gloc = gloc+2*real(a(k,j,i)*conjg(d(k,j,i))+            &
                           b(k,j,i)*conjg(e(k,j,i))+c(k,j,i)*                &
                           conjg(f(k,j,i)))/((1.+alp**2*ka2(k,j,i))*         &
                           real(n,kind=GP)**6)
                  END DO
               END DO
            END DO
         ELSE
            DO i = ista,iend
               DO j = 1,n
                  DO k = 1,n
                     gloc = gloc+2*real(a(k,j,i)*conjg(d(k,j,i))+            &
                           b(k,j,i)*conjg(e(k,j,i))+c(k,j,i)*                &
                           conjg(f(k,j,i)))/((1.+alp**2*ka2(k,j,i))*         &
                           real(n,kind=GP)**6)
                  END DO
               END DO
            END DO
         ENDIF
!
! Computes the averaged inner product between 
! the curl of the fields
!
      ELSE
         IF (ista.eq.1) THEN
            DO j = 1,n
               DO k = 1,n
                  gloc = gloc+ka2(k,j,1)*real(a(k,j,1)*conjg(d(k,j,1))+      &
                        b(k,j,1)*conjg(e(k,j,1))+c(k,j,1)*                   &
                        conjg(f(k,j,1)))/((1.+alp**2*ka2(k,j,1))*            &
                        real(n,kind=GP)**6)
               END DO
            END DO
            DO i = 2,iend
               DO j = 1,n
                  DO k = 1,n
                     gloc = gloc+2*ka2(k,j,i)*real(a(k,j,i)*conjg(d(k,j,i))+ &
                           b(k,j,i)*conjg(e(k,j,i))+c(k,j,i)*                &
                           conjg(f(k,j,i)))/((1.+alp**2*ka2(k,j,i))*         &
                           real(n,kind=GP)**6)
                  END DO
               END DO
            END DO
         ELSE
            DO i = ista,iend
               DO j = 1,n
                  DO k = 1,n
                     gloc = gloc+2*ka2(k,j,i)*real(a(k,j,i)*conjg(d(k,j,i))+ &
                           b(k,j,i)*conjg(e(k,j,i))+c(k,j,i)*                &
                           conjg(f(k,j,i)))/((1.+alp**2*ka2(k,j,i))*         &
                           real(n)**6)
                  END DO
               END DO
            END DO
         ENDIF
      ENDIF
!
! Computes the reduction between nodes
!
      CALL MPI_REDUCE(gloc,g,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,    &
                        MPI_COMM_WORLD,ierr)

      RETURN
      END SUBROUTINE across

!*****************************************************************
      SUBROUTINE ahelicity(a,b,c,d,alp)
!-----------------------------------------------------------------
!
! Computes the mean alpha-helicity of a vector 
! field. The output is valid only in the first 
! node.
!
! Parameters
!     a  : input matrix in the x-direction
!     b  : input matrix in the y-direction
!     c  : input matrix in the z-direction
!     d  : at the output contains the helicity
!     alp: value of alpha
!
      USE fprecision
      USE commtypes
      USE kes
      USE grid
      USE mpivars
      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(OUT) :: d
      DOUBLE PRECISION              :: dloc
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,n,ista:iend) :: a,b,c
      COMPLEX(KIND=GP), DIMENSION(n,n,ista:iend)             :: c1
      REAL(KIND=GP), INTENT(IN)    :: alp
      INTEGER             :: i,j,k

      dloc = 0.

      CALL rotor3(b,c,c1,1)
      IF (ista.eq.1) THEN
         DO j = 1,n
            DO k = 1,n
                 dloc = dloc+real(a(k,j,1)*conjg(c1(k,j,1)))      &
                        /((1.+alp**2*ka2(k,j,1))**2*real(n,kind=GP)**6)
            END DO
         END DO
         DO i = 2,iend
            DO j = 1,n
               DO k = 1,n
                  dloc = dloc+2*real(a(k,j,i)*conjg(c1(k,j,i)))   &
                         /((1.+alp**2*ka2(k,j,i))**2*real(n,kind=GP)**6)
               END DO
            END DO
         END DO
      ELSE
         DO i = ista,iend
            DO j = 1,n
               DO k = 1,n
                  dloc = dloc+2*real(a(k,j,i)*conjg(c1(k,j,i)))   &
                         /((1.+alp**2*ka2(k,j,i))**2*real(n,kind=GP)**6)
               END DO
            END DO
         END DO
      ENDIF
      CALL rotor3(a,c,c1,2)
      IF (ista.eq.1) THEN
         DO j = 1,n
            DO k = 1,n
               dloc = dloc+real(b(k,j,1)*conjg(c1(k,j,1)))        &
                      /((1.+alp**2*ka2(k,j,1))**2*real(n,kind=GP)**6)
            END DO
         END DO
         DO i = 2,iend
            DO j = 1,n
               DO k = 1,n
                  dloc = dloc+2*real(b(k,j,i)*conjg(c1(k,j,i)))   &
                         /((1.+alp**2*ka2(k,j,i))**2*real(n,kind=GP)**6)
               END DO
            END DO
         END DO
      ELSE
         DO i = ista,iend
            DO j = 1,n
               DO k = 1,n
                  dloc = dloc+2*real(b(k,j,i)*conjg(c1(k,j,i)))   &
                         /((1.+alp**2*ka2(k,j,i))**2*real(n,kind=GP)**6)
               END DO
            END DO
         END DO
      ENDIF
      CALL rotor3(a,b,c1,3)
      IF (ista.eq.1) THEN
         DO j = 1,n
            DO k = 1,n
               dloc = dloc+real(c(k,j,1)*conjg(c1(k,j,1)))        &
                      /((1.+alp**2*ka2(k,j,1))**2*real(n,kind=GP)**6)
            END DO
         END DO
         DO i = 2,iend
            DO j = 1,n
               DO k = 1,n
                  dloc = dloc+2*real(c(k,j,i)*conjg(c1(k,j,i)))   &
                         /((1.+alp**2*ka2(k,j,i))**2*real(n,kind=GP)**6)
               END DO
            END DO
         END DO
      ELSE
         DO i = ista,iend
            DO j = 1,n
               DO k = 1,n
                  dloc = dloc+2*real(c(k,j,i)*conjg(c1(k,j,i)))   &
                         /((1.+alp**2*ka2(k,j,i))**2*real(n,kind=GP)**6)
               END DO
            END DO
         END DO
      ENDIF

      CALL MPI_REDUCE(dloc,d,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,    &
                      MPI_COMM_WORLD,ierr)

      RETURN
      END SUBROUTINE ahelicity

!*****************************************************************
      SUBROUTINE ahdcheck(a,b,c,d,e,f,t,dt,alp,hel)
!-----------------------------------------------------------------
!
! Consistency check for the conservation of energy 
! and helicity in the alpha-model Navier-Stokes 
! equations.
!
! Parameters
!     a  : velocity field in the x-direction
!     b  : velocity field in the y-direction
!     c  : velocity field in the z-direction
!     d  : force in the x-direction
!     e  : force in the y-direction
!     f  : force in the z-direction
!     t  : number of time steps made
!     dt : time step
!     alp: value of alpha
!     hel: =1 computes the kinetic helicity
!          =0 skips kinetic helicity computation
!
      USE fprecision
      USE kes
      USE grid
      USE mpivars
      IMPLICIT NONE

      DOUBLE PRECISION :: eng,ens
      DOUBLE PRECISION :: pot,khe
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,n,ista:iend) :: a,b,c
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,n,ista:iend) :: d,e,f
      COMPLEX(KIND=GP), DIMENSION(n,n,ista:iend)             :: c1,c2,c3
      REAL(KIND=GP), INTENT(IN)    :: alp
      REAL(KIND=GP), INTENT(IN)    :: dt
      INTEGER, INTENT(IN) :: hel
      INTEGER             :: i,j,k
      INTEGER, INTENT(IN) :: t

!
! Computes the mean energy, enstrophy, 
! and kinetic helicity
!
      CALL aenergy(a,b,c,eng,alp,1)
      CALL aenergy(a,b,c,ens,alp,0)
      IF (hel.eq.1) THEN
         CALL ahelicity(a,b,c,khe,alp)
      ENDIF
!
! Computes the energy injection rate
!
      CALL across(a,b,c,d,e,f,pot,alp,1)
!
! Creates external files to store the results
!
      IF (myrank.eq.0) THEN
         OPEN(1,file='abalance.txt',position='append')
         WRITE(1,10) (t-1)*dt,eng,ens,pot
   10    FORMAT( E13.6,E22.14,E22.14,E22.14 )
         CLOSE(1)
         IF (hel.eq.1) THEN
            OPEN(1,file='ahelicity.txt',position='append')
            WRITE(1,*) (t-1)*dt,khe
            CLOSE(1)
         ENDIF
      ENDIF

      RETURN
      END SUBROUTINE ahdcheck

!*****************************************************************
      SUBROUTINE aspectrum(a,b,c,alp,nmb,kin,hel)
!-----------------------------------------------------------------
!
! Computes the alpha energy and alpha 
! helicity power spectrum. The output 
! is written to a file by the first node.
!
! Parameters
!     a  : input matrix in the x-direction
!     b  : input matrix in the y-direction
!     c  : input matrix in the z-direction
!     alp: value of alpha
!     nmb: the extension used when writting the file
!     kin: =1 computes the kinetic spectrum
!          =0 computes the magnetic spectrum
!     hel: =1 computes the helicity spectrum
!          =0 skips helicity spectrum computation
!
      USE fprecision
      USE commtypes
      USE kes
      USE grid
      USE mpivars
      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(n/2+1)            :: Ek,Ektot
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,n,ista:iend) :: a,b,c
      COMPLEX(KIND=GP), DIMENSION(n,n,ista:iend)             :: c1,c2,c3
      REAL(KIND=GP), INTENT(IN)    :: alp
      INTEGER, INTENT(IN) :: kin,hel
      INTEGER             :: i,j,k
      INTEGER             :: kmn
      CHARACTER(len=*), INTENT(IN) :: nmb

!
! Sets Ek to zero
!
      DO i = 1,n/2+1
         Ek(i) = 0.
      END DO
!
! Computes the curl of the field if needed
!
      IF ((kin.eq.0).or.(hel.eq.1)) THEN
         CALL rotor3(b,c,c1,1)
         CALL rotor3(a,c,c2,2)
         CALL rotor3(a,b,c3,3)
      ENDIF
!
! Computes the kinetic energy spectrum
!
      IF (kin.eq.1) THEN
         IF (ista.eq.1) THEN
            DO j = 1,n
               DO k = 1,n
                  kmn = int(sqrt(ka2(k,j,1))+.501)
                  IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                     Ek(kmn) = Ek(kmn)+                               &
                       (abs(a(k,j,1))**2+abs(b(k,j,1))**2+            &
                       abs(c(k,j,1))**2)/((1.+alp**2*ka2(k,j,1))*     &
                       real(n,kind=GP)**6)
                  ENDIF
               END DO
            END DO
            DO i = 2,iend
               DO j = 1,n
                  DO k = 1,n
                     kmn = int(sqrt(ka2(k,j,i))+.501)
                     IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                        Ek(kmn) = Ek(kmn)+                            &
                          2*(abs(a(k,j,i))**2+abs(b(k,j,i))**2+       &
                          abs(c(k,j,i))**2)/((1.+alp**2*ka2(k,j,i))*  &
                          real(n,kind=GP)**6)
                     ENDIF
                  END DO
               END DO
            END DO
          ELSE
            DO i = ista,iend
               DO j = 1,n
                  DO k = 1,n
                     kmn = int(sqrt(ka2(k,j,i))+.501)
                     IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                        Ek(kmn) = Ek(kmn)+                            &
                          2*(abs(a(k,j,i))**2+abs(b(k,j,i))**2+       &
                          abs(c(k,j,i))**2)/((1.+alp**2*ka2(k,j,i))*  &
                          real(n,kind=GP)**6)
                     ENDIF
                  END DO
               END DO
            END DO
          ENDIF
!
! Computes the magnetic energy spectrum
!
      ELSE
         IF (ista.eq.1) THEN
            DO j = 1,n
               DO k = 1,n
                  kmn = int(sqrt(ka2(k,j,1))+.501)
                  IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                     Ek(kmn) = Ek(kmn)+                               &
                       (abs(c1(k,j,1))**2+abs(c2(k,j,1))**2+          &
                       abs(c3(k,j,1))**2)/((1.+alp**2*ka2(k,j,1))*    &
                       real(n,kind=GP)**6)
                  ENDIF
               END DO
            END DO
            DO i = 2,iend
               DO j = 1,n
                  DO k = 1,n
                     kmn = int(sqrt(ka2(k,j,i))+.501)
                     IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                        Ek(kmn) = Ek(kmn)+                            &
                          2*(abs(c1(k,j,i))**2+abs(c2(k,j,i))**2+     &
                          abs(c3(k,j,i))**2)/((1.+alp**2*ka2(k,j,i))* &
                          real(n,kind=GP)**6)
                     ENDIF
                  END DO
               END DO
            END DO
         ELSE
            DO i = ista,iend
               DO j = 1,n
                  DO k = 1,n
                     kmn = int(sqrt(ka2(k,j,i))+.501)
                     IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                        Ek(kmn) = Ek(kmn)+                            &
                          2*(abs(c1(k,j,i))**2+abs(c2(k,j,i))**2+     &
                          abs(c3(k,j,i))**2)/((1.+alp**2*ka2(k,j,i))* &
                          real(n,kind=GP)**6)
                     ENDIF
                  END DO
               END DO
            END DO
         ENDIF
      ENDIF
!
! Computes the reduction between nodes
! and exports the result to a file
!
      CALL MPI_REDUCE(Ek,Ektot,n/2+1,MPI_DOUBLE_PRECISION,MPI_SUM,0,  &
                      MPI_COMM_WORLD,ierr)
      IF (myrank.eq.0) THEN
         IF (kin.eq.1) THEN
            OPEN(1,file='kaspectrum.' // nmb // '.txt')
         ELSE
            OPEN(1,file='maspectrum.' // nmb // '.txt')
         ENDIF
         WRITE(1,20) Ektot
   20    FORMAT( E23.15 ) 
         CLOSE(1)
      ENDIF
!
! Computes the helicity spectrum
!
      IF (hel.eq.1) THEN
         DO i = 1,n/2+1
            Ek(i) = 0.
         END DO
         IF (ista.eq.1) THEN
            DO j = 1,n
               DO k = 1,n
                  kmn = int(sqrt(ka2(k,j,1))+.501)
                  IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                     Ek(kmn) = Ek(kmn)+                               &
                       (real(a(k,j,1)*conjg(c1(k,j,1)))+              &
                       real(b(k,j,1)*conjg(c2(k,j,1)))+               &
                       real(c(k,j,1)*conjg(c3(k,j,1))))/              &
                       ((1.+alp**2*ka2(k,j,1))**2*real(n,kind=GP)**6)
                  ENDIF
               END DO
            END DO
            DO i = 2,iend
               DO j = 1,n
                  DO k = 1,n
                     kmn = int(sqrt(ka2(k,j,i))+.501)
                     IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                        Ek(kmn) = Ek(kmn)+                            &
                          2*(real(a(k,j,i)*conjg(c1(k,j,i)))+         &
                          real(b(k,j,i)*conjg(c2(k,j,i)))+            &
                          real(c(k,j,i)*conjg(c3(k,j,i))))/           &
                          ((1.+alp**2*ka2(k,j,i))**2*real(n,kind=GP)**6)
                     ENDIF
                 END DO
               END DO
            END DO
         ELSE
            DO i = ista,iend
               DO j = 1,n
                  DO k = 1,n
                     kmn = int(sqrt(ka2(k,j,i))+.501)
                     IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                        Ek(kmn) = Ek(kmn)+                            &
                          2*(real(a(k,j,i)*conjg(c1(k,j,i)))+         &
                          real(b(k,j,i)*conjg(c2(k,j,i)))+            &
                          real(c(k,j,i)*conjg(c3(k,j,i))))/           &
                          ((1.+alp**2*ka2(k,j,i))**2*real(n,kind=GP)**6)
                     ENDIF
                  END DO
               END DO
            END DO
         ENDIF
!
! Computes the reduction between nodes
! and exports the result to a file
!
         CALL MPI_REDUCE(Ek,Ektot,n/2+1,MPI_DOUBLE_PRECISION,MPI_SUM, &
                         0,MPI_COMM_WORLD,ierr)
         IF (myrank.eq.0) THEN
            IF (kin.eq.1) THEN
               OPEN(1,file='kahelicity.' // nmb // '.txt')
            ELSE
               OPEN(1,file='mahelicity.' // nmb // '.txt')
            ENDIF
            WRITE(1,30) Ektot
   30       FORMAT( E23.15 )
            CLOSE(1)
         ENDIF
      ENDIF

      RETURN
      END SUBROUTINE aspectrum

!*****************************************************************
      SUBROUTINE amhdcheck(a,b,c,ma,mb,mc,t,dt,alpk,alpm,hel,crs)
!-----------------------------------------------------------------
!
! Consistency check for the conservation of the total
! energy and helicity in alpha-model MHD
!
! Parameters
!     a   : velocity field in the x-direction
!     b   : velocity field in the y-direction
!     c   : velocity field in the z-direction
!     ma  : vector potential in the x-direction
!     mb  : vector potential in the y-direction
!     mc  : vector potential in the z-direction
!     t   : number of time steps made
!     dt  : time step
!     alpk: kinetic alpha
!     alpm: magnetic alpha
!     hel : =0 skips helicity computation
!           =1 computes the helicity
!     crs : =0 skips cross helicity computation
!           =1 computes cross helicity
!
      USE fprecision
      USE kes
      USE grid
      USE mpivars
      IMPLICIT NONE

      DOUBLE PRECISION :: engk,engm,eng,ens
      DOUBLE PRECISION :: helk,helm,cur
      DOUBLE PRECISION :: asq,crh
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,n,ista:iend) :: a,b,c
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,n,ista:iend) :: ma,mb,mc
      COMPLEX(KIND=GP), DIMENSION(n,n,ista:iend)             :: c1,c2,c3
      REAL(KIND=GP), INTENT(IN)    :: dt
      REAL(KIND=GP), INTENT(IN)    :: alpk,alpm
      INTEGER, INTENT(IN) :: hel,crs
      INTEGER             :: i,j,k
      INTEGER, INTENT(IN) :: t

!
! Computes the mean energy, enstrophy 
! and square current, and the kinetic 
! and magnetic helicity
!
      call aenergy(a,b,c,engk,alpk,1)
      call aenergy(a,b,c,ens,alpk,0)
      call rotor3(mb,mc,c1,1)
      call rotor3(ma,mc,c2,2)
      call rotor3(ma,mb,c3,3)
      call aenergy(c1,c2,c3,engm,alpm,1)
      call aenergy(c1,c2,c3,cur,alpm,0)
      eng = engk+engm
      IF (hel.eq.1) THEN
         CALL ahelicity(a,b,c,helk,alpk)
         CALL ahelicity(ma,mb,mc,helm,alpm)
      ENDIF
!
! Computes the square vector potential 
! and the cross helicity
!
      IF (crs.eq.1) THEN
         CALL aenergy(ma,mb,mc,asq,alpm,1)
         CALL across(a,b,c,c1,c2,c3,crh,alpm,1)
      ENDIF
!
! Creates external files to store the results
!
      IF (myrank.eq.0) THEN
         OPEN(1,file='abalance.txt',position='append')
         WRITE(1,10) (t-1)*dt,eng,ens,cur
   10    FORMAT( E13.6,E22.14,E22.14,E22.14,E22.14 )
         CLOSE(1)
         OPEN(1,file='aenergy.txt',position='append')
         WRITE(1,*) (t-1)*dt,engk,engm
         CLOSE(1)
         IF (hel.eq.1) THEN
            OPEN(1,file='ahelicity.txt',position='append')
            WRITE(1,*) (t-1)*dt,helk,helm
            CLOSE(1)
         ENDIF
         IF (crs.eq.1) THEN
            OPEN(1,file='across.txt',position='append')
            WRITE(1,*) (t-1)*dt,crh,asq
            CLOSE(1)
         ENDIF
      ENDIF

      RETURN
      END SUBROUTINE amhdcheck

!*****************************************************************
      SUBROUTINE aentrans(a,b,c,d,e,f,alp,nmb,kin)
!-----------------------------------------------------------------
!
! Computes the alpha energy transfer in Fourier 
! space in 3D. The output is written to a file 
! by the first node.
!
! Parameters
!     a  : field component in the x-direction
!     b  : field component in the y-direction
!     c  : field component in the z-direction
!     d  : nonlinear term in the x-direction
!     e  : nonlinear term in the y-direction
!     f  : nonlinear term in the z-direction
!     alp: value of alpha for the first field
!     nmb: the extension used when writting the file
!     kin: =0 computes the magnetic energy transfer
!          =1 computes the kinetic energy transfer
!          =2 computes the Lorentz force transfer
!
      USE fprecision
      USE commtypes
      USE kes
      USE grid
      USE mpivars
      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(n/2+1)            :: Ek,Ektot
      COMPLEX(KIND=GP), DIMENSION(n,n,ista:iend), INTENT(IN) :: a,b,c
      COMPLEX(KIND=GP), DIMENSION(n,n,ista:iend), INTENT(IN) :: d,e,f
      REAL(KIND=GP), INTENT(IN)    :: alp
      REAL(KIND=GP)                :: tmp
      INTEGER, INTENT(IN) :: kin
      INTEGER             :: i,j,k
      INTEGER             :: kmn
      CHARACTER(len=*), INTENT(IN) :: nmb

!
! Sets Ek to zero
!
      DO i = 1,n/2+1
         Ek(i) = 0.
      END DO
!
! Computes the kinetic energy transfer
!
      IF ((kin.eq.1).or.(kin.eq.2)) THEN
         IF (ista.eq.1) THEN
            DO j = 1,n
               DO k = 1,n
                  kmn = int(sqrt(ka2(k,j,1))+.501)
                  IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                     Ek(kmn) = Ek(kmn)+                             &
                       (real(a(k,j,1)*conjg(d(k,j,1)))+             &
                       real(b(k,j,1)*conjg(e(k,j,1)))+              &
                       real(c(k,j,1)*conjg(f(k,j,1))))/             &
                       ((1.+alp**2*ka2(k,j,1))*real(n,kind=GP)**6)
                  ENDIF
               END DO
            END DO
            DO i = 2,iend
               DO j = 1,n
                  DO k = 1,n
                     kmn = int(sqrt(ka2(k,j,i))+.501)
                     IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                        Ek(kmn) = Ek(kmn)+                          &
                          2*(real(a(k,j,i)*conjg(d(k,j,i)))+        &
                          real(b(k,j,i)*conjg(e(k,j,i)))+           &
                          real(c(k,j,i)*conjg(f(k,j,i))))/          &
                          ((1.+alp**2*ka2(k,j,i))*real(n,kind=GP)**6)
                     ENDIF
                 END DO
               END DO
            END DO
         ELSE
            DO i = ista,iend
               DO j = 1,n
                  DO k = 1,n
                     kmn = int(sqrt(ka2(k,j,i))+.501)
                     IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                        Ek(kmn) = Ek(kmn)+                          &
                          2*(real(a(k,j,i)*conjg(d(k,j,i)))+        &
                          real(b(k,j,i)*conjg(e(k,j,i)))+           &
                          real(c(k,j,i)*conjg(f(k,j,i))))/          &
                          ((1.+alp**2*ka2(k,j,i))*real(n,kind=GP)**6)
                     ENDIF
                  END DO
               END DO
            END DO
         ENDIF
!
! Computes the magnetic energy transfer
!
      ELSE
         tmp = 1./real(n,kind=GP)**6
         IF (ista.eq.1) THEN
            DO j = 1,n
               DO k = 1,n
                  kmn = int(sqrt(ka2(k,j,1))+.501)
                  IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                     Ek(kmn) = Ek(kmn)+ka2(k,j,1)*                  &
                       (real(a(k,j,1)*conjg(d(k,j,1)))+             &
                       real(b(k,j,1)*conjg(e(k,j,1)))+              &
                       real(c(k,j,1)*conjg(f(k,j,1))))*tmp
                  ENDIF
               END DO
            END DO
            DO i = 2,iend
               DO j = 1,n
                  DO k = 1,n
                     kmn = int(sqrt(ka2(k,j,i))+.501)
                     IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                        Ek(kmn) = Ek(kmn)+2*ka2(k,j,i)*             &
                          (real(a(k,j,i)*conjg(d(k,j,i)))+          &
                          real(b(k,j,i)*conjg(e(k,j,i)))+           &
                          real(c(k,j,i)*conjg(f(k,j,i))))*tmp
                     ENDIF
                 END DO
               END DO
            END DO
         ELSE
            DO i = ista,iend
               DO j = 1,n
                  DO k = 1,n
                     kmn = int(sqrt(ka2(k,j,i))+.501)
                     IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                        Ek(kmn) = Ek(kmn)+2*ka2(k,j,i)*             &
                          (real(a(k,j,i)*conjg(d(k,j,i)))+          &
                          real(b(k,j,i)*conjg(e(k,j,i)))+           &
                          real(c(k,j,i)*conjg(f(k,j,i))))*tmp
                     ENDIF
                  END DO
               END DO
            END DO
         ENDIF
      ENDIF
!
! Computes the reduction between nodes
! and exports the result to a file
!
      CALL MPI_REDUCE(Ek,Ektot,n/2+1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
      IF (myrank.eq.0) THEN
         IF (kin.eq.0) THEN
            OPEN(1,file='amtransfer.' // nmb // '.txt')
         ELSEIF (kin.eq.1) THEN
            OPEN(1,file='aktransfer.' // nmb // '.txt')
         ELSE
            OPEN(1,file='ajtransfer.' // nmb // '.txt')
         ENDIF
         WRITE(1,40) Ektot
   40    FORMAT( E23.15 ) 
         CLOSE(1)
      ENDIF

      RETURN
      END SUBROUTINE aentrans
