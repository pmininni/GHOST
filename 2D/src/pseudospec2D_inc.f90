!=================================================================
! PSEUDOSPECTRAL subroutines
!
! Subroutines for computing spatial derivatives and nonlinear 
! terms in incompressible HD, MHD and Hall-MHD equations in 2D 
! using a pseudo-spectral method. You should use the FFTPLANS 
! and MPIVARS modules (see the file 'fftp2D_mod.f90') in each 
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
      SUBROUTINE derivk2(a,b,dir)
!-----------------------------------------------------------------
!
! Two-dimensional derivative of the matrix 'a'
!
! Parameters
!     a  : input matrix
!     b  : at the output contains the derivative da/dk_dir
!     dir: =1 derivative in the x-direction
!          =2 derivative in the y-direction
!
      USE ali
      USE kes
      USE var
      USE grid
      USE mpivars
      IMPLICIT NONE

      COMPLEX, DIMENSION(n,ista:iend) :: a,b
      INTEGER :: dir
      INTEGER :: i,j

!
! Derivative in the x-direction
!
      IF (dir.eq.1) THEN
         DO i = ista,iend
            DO j = 1,n
               IF ((ka2(j,i).le.kmax).and.(ka2(j,i).ge.tiny)) THEN
                  b(j,i) = im*ka(i)*a(j,i)
               ELSE
                  b(j,i) = 0.
               ENDIF
            END DO
         END DO
!
! Derivative in the y-direction
!
      ELSE
         DO i = ista,iend
            DO j = 1,n
               IF ((ka2(j,i).le.kmax).and.(ka2(j,i).ge.tiny)) THEN
                  b(j,i) = im*ka(j)*a(j,i)
               ELSE
                  b(j,i) = 0.
               ENDIF
            END DO
         END DO
      ENDIF

      RETURN
      END SUBROUTINE derivk2

!*****************************************************************
      SUBROUTINE laplak2(a,b)
!-----------------------------------------------------------------
!
! Two-dimensional Laplacian of the matrix 'a'
!
! Parameters
!     a: input matrix
!     b: at the output contains the Laplacian d2a/dka2
!
      USE kes
      USE grid
      USE mpivars
      IMPLICIT NONE

      COMPLEX, DIMENSION(n,ista:iend) :: a,b
      INTEGER :: i,j

      DO i = ista,iend
         DO j = 1,n
            b(j,i) = -ka2(j,i)*a(j,i)
         END DO
      END DO

      RETURN
      END SUBROUTINE laplak2

!*****************************************************************
      SUBROUTINE poisson(a,b,c)
!-----------------------------------------------------------------
!
! Poisson bracket of the scalar fields A and B 
! in real space.
!
! Parameters
!     a: input matrix
!     b: input matrix
!     c: Poisson bracket {A,B} [output]
!
      USE mpivars
      USE grid
      USE fft
      IMPLICIT NONE

      COMPLEX, DIMENSION(n,ista:iend) :: a,b,c
      COMPLEX, DIMENSION(n,ista:iend) :: c1,c2
      REAL, DIMENSION(n,jsta:jend)    :: r1,r2,r3
      INTEGER :: i,j

!
! Computes dA/dx.dB/dy
!
      CALL derivk2(a,c1,1)
      CALL derivk2(b,c2,2)
      CALL fftp2d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      CALL fftp2d_complex_to_real(plancr,c2,r2,MPI_COMM_WORLD)
      DO j = jsta,jend
         DO i = 1,n
            r3(i,j) = r1(i,j)*r2(i,j)
         END DO
      END DO
!
! Computes dA/dy.dB/dx
!
      CALL derivk2(a,c1,2)
      CALL derivk2(b,c2,1)
      CALL fftp2d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      CALL fftp2d_complex_to_real(plancr,c2,r2,MPI_COMM_WORLD)
      DO j = jsta,jend
         DO i = 1,n
            r3(i,j) = (r3(i,j)-r1(i,j)*r2(i,j))/float(n)**4
         END DO
      END DO

      CALL fftp2d_real_to_complex(planrc,r3,c,MPI_COMM_WORLD)

      RETURN
      END SUBROUTINE poisson

!*****************************************************************
      SUBROUTINE energy(a,b,kin)
!-----------------------------------------------------------------
!
! Computes the mean kinetic or magnetic energy in 2D,
! and the mean square current density or vorticity. 
! The output is valid only in the first node.
!
! Parameters
!     a  : input matrix with the scalar field
!     b  : at the output contains the energy
!     kin: =2 computes the square of the scalar field
!          =1 computes the energy
!          =0 computes the current or vorticity
!
      USE kes
      USE grid
      USE mpivars
      IMPLICIT NONE

      COMPLEX, DIMENSION(n,ista:iend) :: a
      REAL*8  :: b
      REAL*8  :: bloc
      REAL    :: tmp
      INTEGER :: kin
      INTEGER :: i,j

      bloc = 0.
      tmp = 1./float(n)**4

!
! Computes the square of the scalar field
!
      IF (kin.eq.2) THEN
         IF (ista.eq.1) THEN
            DO j = 1,n
               bloc = bloc+abs(a(j,1))**2*tmp
            END DO
            DO i = 2,iend
               DO j = 1,n
                  bloc = bloc+2*abs(a(j,i))**2*tmp
               END DO
            END DO
         ELSE
            DO i = ista,iend
               DO j = 1,n
                  bloc = bloc+2*abs(a(j,i))**2*tmp
               END DO
            END DO
         ENDIF
!
! Computes the energy
!
      ELSE IF (kin.eq.1) THEN
         IF (ista.eq.1) THEN
            DO j = 1,n
               bloc = bloc+ka2(j,1)*abs(a(j,1))**2*tmp
            END DO
            DO i = 2,iend
               DO j = 1,n
                  bloc = bloc+2*ka2(j,i)*abs(a(j,i))**2*tmp
               END DO
            END DO
         ELSE
            DO i = ista,iend
               DO j = 1,n
                  bloc = bloc+2*ka2(j,i)*abs(a(j,i))**2*tmp
               END DO
            END DO
         ENDIF
!
! Computes the current or vorticity
!
      ELSE
         IF (ista.eq.1) THEN
            DO j = 1,n
               bloc = bloc+ka2(j,1)**2*abs(a(j,1))**2*tmp
            END DO
            DO i = 2,iend
               DO j = 1,n
                  bloc = bloc+2*ka2(j,i)**2*abs(a(j,i))**2*tmp
               END DO
            END DO
         ELSE
            DO i = ista,iend
               DO j = 1,n
                  bloc = bloc+2*ka2(j,i)**2*abs(a(j,i))**2*tmp
               END DO
            END DO
         ENDIF
      ENDIF
!
! Computes the reduction between nodes
!
      CALL MPI_REDUCE(bloc,b,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)

      RETURN
      END SUBROUTINE energy

!*****************************************************************
      SUBROUTINE maxabs(a,b,c,d)
!-----------------------------------------------------------------
!
! Computes the maximum absolute value of the 
! vorticity or current density. The output is 
! only valid in the first node.
!
! Parameters
!     a: input field
!     b: at the output contains the maximum value
!     c: at the output contains the i-index where the maximum is
!     d: at the output contains the j-index where the maximum is
!
      USE fft
      USE grid
      USE mpivars
      IMPLICIT NONE

      COMPLEX, DIMENSION(n,ista:iend) :: a
      COMPLEX, DIMENSION(n,ista:iend) :: c1
      REAL, DIMENSION(n,jsta:jend)    :: r1
      REAL    :: b,bloc
      INTEGER :: c,d
      INTEGER :: i,j

      CALL laplak2(a,c1)
      CALL fftp2d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      bloc = 0.
      DO j = jsta,jend
         DO i = 1,n
            bloc = max(bloc,abs(r1(i,j)))
         END DO
      END DO
      bloc = bloc/float(n)**2
      CALL MPI_REDUCE(bloc,b,1,MPI_REAL,MPI_MAX,0,MPI_COMM_WORLD,ierr)

      RETURN
      END SUBROUTINE maxabs

!*****************************************************************
      SUBROUTINE hdcheck(a,b,t)
!-----------------------------------------------------------------
!
! Consistency check for the conservation of energy in HD 2D
!
! Parameters
!     a  : streamfunction
!     b  : external force
!     t  : time
!
      USE kes
      USE grid
      USE mpivars

      IMPLICIT NONE

      COMPLEX, DIMENSION(n,ista:iend) :: a,b
      REAL*8  :: eng,ens,pot,tmp
      REAL    :: t
      REAL    :: tmq
      INTEGER :: i,j

      pot = 0.
      tmp = 0.
      tmq = 1./float(n)**4

!
! Computes the mean energy and enstrophy
!
      CALL energy(a,eng,1)
      CALL energy(a,ens,0)
!
! Computes the energy injection rate
!
      IF (ista.eq.1) THEN
         DO j = 1,n
            tmp = tmp+ka2(j,1)*real(b(j,1)*conjg(a(j,1)))*tmq
         END DO
         DO i = 2,iend
            DO j = 1,n
               tmp = tmp+2*ka2(j,i)*real(b(j,i)*conjg(a(j,i)))*tmq
            END DO
         END DO
      ELSE
         DO i = ista,iend
            DO j = 1,n
               tmp = tmp+2*ka2(j,i)*real(b(j,i)*conjg(a(j,i)))*tmq
            END DO
         END DO
      ENDIF
      CALL MPI_REDUCE(tmp,pot,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
!
! Creates external files to store the results
!
      IF (myrank.eq.0) THEN
         OPEN(1,file='balance.txt',position='append')
         WRITE(1,10) t,eng,ens,pot
   10    FORMAT( E13.6,D22.14,D22.14,D22.14 )
         CLOSE(1)
      ENDIF      

      RETURN
      END SUBROUTINE hdcheck

!*****************************************************************
      SUBROUTINE mhdcheck(a,b,c,d,t)
!-----------------------------------------------------------------
!
! Consistency check for the conservation of energy in MHD 2D
!
! Parameters
!     a  : streamfunction
!     b  : vector potential
!     c  : external kinetic force
!     d  : external magnetic force
!     t  : time
!
      USE kes
      USE grid
      USE mpivars

      IMPLICIT NONE

      COMPLEX, DIMENSION(n,ista:iend) :: a,b,c,d
      REAL*8  :: engk,engm,eng,asq,udb
      REAL*8  :: ens,cur,pot,tmp
      REAL    :: t
      REAL    :: tmq
      INTEGER :: i,j

      pot = 0.
      tmp = 0.
      tmq = 1./float(n)**4

!
! Computes the mean energy, enstrophy, square 
! current, and square vector potential.
!
      CALL energy(a,engk,1)
      CALL energy(b,engm,1)
      CALL energy(a,ens,0)
      CALL energy(b,cur,0)
      CALL energy(b,asq,2)
      eng = engk+engm
!
! Computes the total energy injection rate
!
      IF (ista.eq.1) THEN
         DO j = 1,n
            tmp = tmp+ka2(j,1)*(real(c(j,1)*conjg(a(j,1)))+      &
                  real(d(j,1)*conjg(b(j,1))))*tmq
         END DO
         DO i = 2,iend
            DO j = 1,n
               tmp = tmp+2*ka2(j,i)*(real(c(j,i)*conjg(a(j,i)))+ &
                     real(d(j,i)*conjg(b(j,i))))*tmq
            END DO
         END DO
      ELSE
         DO i = ista,iend
            DO j = 1,n
               tmp = tmp+2*ka2(j,i)*(real(c(j,i)*conjg(a(j,i)))+ &
                     real(d(j,i)*conjg(b(j,i))))*tmq
            END DO
         END DO
      ENDIF
      CALL MPI_REDUCE(tmp,pot,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,  &
                      MPI_COMM_WORLD,ierr)

!
! Computes the cross correlation between
! velocity and magnetic fields
!
      IF (ista.eq.1) THEN
         DO j = 1,n
            tmp = tmp+ka2(j,1)*real(b(j,1)*conjg(a(j,1)))*tmq
         END DO
         DO i = 2,iend
            DO j = 1,n
               tmp = tmp+2*ka2(j,i)*real(b(j,i)*conjg(a(j,i)))*tmq
            END DO
         END DO
      ELSE
         DO i = ista,iend
            DO j = 1,n
               tmp = tmp+2*ka2(j,i)*real(b(j,i)*conjg(a(j,i)))*tmq
            END DO
         END DO
      ENDIF
      CALL MPI_REDUCE(tmp,udb,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
!
! Creates external files to store the results
!
      IF (myrank.eq.0) THEN
         OPEN(1,file='balance.txt',position='append')
         WRITE(1,10) t,eng,ens,cur,pot
   10    FORMAT( E13.6,E22.14,E22.14,E22.14,E22.14 )
         CLOSE(1)
         OPEN(1,file='cross.txt',position='append')
         WRITE(1,*) udb,asq
         CLOSE(1)
         OPEN(1,file='energy.txt',position='append')
         WRITE(1,*) engk,engm
         CLOSE(1)
      ENDIF      

      RETURN
      END SUBROUTINE mhdcheck

!*****************************************************************
      SUBROUTINE spectrum(a,ext,kin)
!-----------------------------------------------------------------
!
! Computes the energy power spectrum in 2D. 
! The output is written to a file by the first node.
!
! Parameters
!     a  : streamfunction or vector potential
!     ext: the extension used when writting the file
!     kin: =1 computes the kinetic spectrum
!          =0 computes the magnetic spectrum
!
      USE kes
      USE grid
      USE mpivars
      IMPLICIT NONE

      REAL*8, DIMENSION(n/2+1)        :: Ek,Ektot
      COMPLEX, DIMENSION(n,ista:iend) :: a
      REAL        :: tmp
      INTEGER     :: kin
      INTEGER     :: kmn
      INTEGER     :: i,j
      CHARACTER*3 :: ext

!
! Sets Ek to zero
!
      DO i = 1,n/2+1
         Ek(i) = 0.
      END DO
!
! Computes the energy spectrum
!
      tmp = 1./float(n)**4
      IF (ista.eq.1) THEN
         DO j = 1,n
            kmn = int(sqrt(ka2(j,1))+.5)
            IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
               Ek(kmn) = Ek(kmn)+ka2(j,1)*abs(a(j,1))**2*tmp
            ENDIF
         END DO
         DO i = 2,iend
            DO j = 1,n
               kmn = int(sqrt(ka2(j,i))+.5)
               IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                  Ek(kmn) = Ek(kmn)+2*ka2(j,i)*abs(a(j,i))**2*tmp
               ENDIF
            END DO
         END DO
      ELSE
         DO i = ista,iend
            DO j = 1,n
               kmn = int(sqrt(ka2(j,i))+.5)
               IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                  Ek(kmn) = Ek(kmn)+2*ka2(j,i)*abs(a(j,i))**2*tmp
               ENDIF
            END DO
         END DO
      ENDIF
!
! Computes the reduction between nodes
! and exports the result to a file
!
      CALL MPI_REDUCE(Ek,Ektot,n/2+1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
      IF (myrank.eq.0) THEN
         IF (kin.eq.1) THEN
            OPEN(1,file='kspectrum.' // ext // '.txt')
         ELSE
            OPEN(1,file='mspectrum.' // ext // '.txt')
         ENDIF
         WRITE(1,20) Ektot
   20    FORMAT( E23.15 ) 
         CLOSE(1)
      ENDIF

      RETURN
      END SUBROUTINE spectrum

!*****************************************************************
      SUBROUTINE vectrans(a,b,ext)
!-----------------------------------------------------------------
!
! Computes the square vector potential transfer in 
! Fourier space in 2D MHD. The output is written 
! to a file by the first node.
!
! Parameters
!     a  : streamfunction
!     b  : vector potential
!     ext: the extension used when writting the file
!
      USE kes
      USE grid
      USE mpivars
      IMPLICIT NONE

      REAL*8, DIMENSION(n/2+1)        :: Ek,Ektot
      COMPLEX, DIMENSION(n,ista:iend) :: a,b,c1
      REAL        :: tmp
      INTEGER     :: kmn
      INTEGER     :: i,j
      CHARACTER*3 :: ext

!
! Sets Ek to zero
!
      DO i = 1,n/2+1
         Ek(i) = 0.
      END DO
!
! Computes the square vector potential flux
!
      tmp = 1./float(n)**4
      CALL poisson(a,b,c1)
      IF (ista.eq.1) THEN
         DO j = 1,n
            kmn = int(sqrt(ka2(j,1))+.5)
            Ek(kmn) = Ek(kmn)+real(c1(j,1)*conjg(b(j,1)))*tmp
         END DO
         DO i = 2,iend
            DO j = 1,n
               kmn = int(sqrt(ka2(j,i))+.5)
               IF (kmn.le.n/2+1) THEN
                  Ek(kmn) = Ek(kmn)+2*real(c1(j,i)*conjg(b(j,i)))*tmp
               ENDIF
            END DO
         END DO
      ELSE
         DO i = ista,iend
            DO j = 1,n
               kmn = int(sqrt(ka2(j,i))+.5)
               IF (kmn.le.n/2+1) THEN
                  Ek(kmn) = Ek(kmn)+2*real(c1(j,i)*conjg(b(j,i)))*tmp
               ENDIF
            END DO
         END DO
      ENDIF
!
! Computes the reduction between nodes
! and exports the result to a file
!
      CALL MPI_REDUCE(Ek,Ektot,n/2+1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
      IF (myrank.eq.0) THEN
         OPEN(1,file='vectrans.' // ext // '.txt')
         WRITE(1,20) Ektot
   20    FORMAT( E23.15 ) 
         CLOSE(1)
      ENDIF

      RETURN
      END SUBROUTINE vectrans
