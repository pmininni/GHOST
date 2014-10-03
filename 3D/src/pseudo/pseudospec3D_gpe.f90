!=================================================================
! PSEUDOSPECTRAL subroutines
!
! Subroutines to compute spatial derivatives and nonlinear 
! terms in the GPE and ARLG equations in 3D using a 
! pseudo-spectral method. You should use the FFTPLANS 
! and MPIVARS modules (see the file 'fftp_mod.f90') in each 
! program that calls any of the subroutines in this file. 
!
! NOTATION: index 'i' is 'x' 
!           index 'j' is 'y'
!           index 'k' is 'z'
!
! 2014 Pablo D. Mininni.
!      Department of Physics,
!      Facultad de Ciencias Exactas y Naturales.
!      Universidad de Buenos Aires.
!      e-mail: mininni@df.uba.ar 
!=================================================================

!*****************************************************************
      SUBROUTINE squareabs(a,b,r,dealias)
!-----------------------------------------------------------------
!
! Pointwise squared absolute value of a complex wavefunction Z.
! Note the output is not normalized (i.e., not divided by N^3).
!
! Parameters
!     a : real part of the wavefunction in Fourier space
!     b : imaginary part of the wavefunction in Fourier space
!     r : |Z|^2 in real space [output]
!     dealias: =0 does not dealias the result
!              =1 dealiases the result

      USE fprecision
      USE commtypes
      USE mpivars
      USE grid
      USE kes
      USE ali
      USE fft
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,n,ista:iend)  :: a,b
      COMPLEX(KIND=GP), DIMENSION(n,n,ista:iend) :: c1
      REAL(KIND=GP), INTENT(OUT), DIMENSION(n,n,ksta:kend)    :: r
      REAL(KIND=GP), DIMENSION(n,n,ksta:kend)    :: r1
      REAL(KIND=GP)       :: rmp
      INTEGER, INTENT(IN) :: dealias
      INTEGER :: i,j,k

!
! Computes the square of the real part of the wavefunction
!
      c1 = a
      CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
         DO j = 1,n
            DO i = 1,n
               r(i,j,k) = r1(i,j,k)**2
            END DO
         END DO
      END DO
!
! Computes the square of the imaginary part of the wavefunction
!
      c1 = b
      CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      rmp = 1.0_GP/real(n,kind=GP)**6
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
         DO j = 1,n
            DO i = 1,n
               r(i,j,k) = (r(i,j,k)+r1(i,j,k)**2)*rmp
            END DO
         END DO
      END DO
!
! Dealiases the result and returs to real space
!
      IF (dealias.eq.1) THEN
         CALL fftp3d_real_to_complex(planrc,r,c1,MPI_COMM_WORLD)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
            DO j = 1,n
               DO k = 1,n
                  IF (ka2(k,j,i).gt.kmax) THEN
                     c1(k,j,i) = 0.0_GP
                  ENDIF
               END DO
            END DO
         END DO
         CALL fftp3d_complex_to_real(plancr,c1,r,MPI_COMM_WORLD)
      ENDIF

      RETURN
      END SUBROUTINE squareabs

!*****************************************************************
      SUBROUTINE nonlgpe(r,a,b)
!-----------------------------------------------------------------
!
! Computes Z.|Z|^2 in real space.
!
! Parameters
!     r  : input matrix with |Z|^2 in real space (not normalized)
!     a  : input with real or imaginary part of Z in Fourier space
!     b  : Z.|Z|^2 in Fourier space [output]
!
      USE fprecision
      USE commtypes
      USE mpivars
      USE grid
      USE fft
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN),  DIMENSION(n,n,ista:iend) :: a
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(n,n,ista:iend) :: b
      REAL(KIND=GP), INTENT(IN), DIMENSION(n,n,ksta:kend)     :: r
      REAL(KIND=GP), DIMENSION(n,n,ksta:kend)    :: r1
      REAL(KIND=GP)    :: rmp
      INTEGER :: i,j,k

!
! Computes Z.|Z|^2
!
      b = a
      CALL fftp3d_complex_to_real(plancr,b,r1,MPI_COMM_WORLD)
      rmp = 1.0_GP/real(n,kind=GP)**6
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
         DO j = 1,n
            DO i = 1,n
               r1(i,j,k) = r(i,j,k)*r1(i,j,k)*rmp
            END DO
         END DO
      END DO
      CALL fftp3d_real_to_complex(planrc,r1,b,MPI_COMM_WORLD)

      RETURN
      END SUBROUTINE nonlgpe

!*****************************************************************
      SUBROUTINE gpecheck(a,b,t,dt)
!-----------------------------------------------------------------
!
! Computes the mass, the kinetic plus quantum energy Ekq, and 
! the quartic term in the energy:
!    Ekq    = 2.alpha^2.|z|^2
!    Equart = alpha.beta.|z|^4
! The potential energy then is
!    Epot = Equart-2*alpha.omegag.mass+alpha.omegag^2/beta
! This quantity should be zero in the condensate. The total 
! energy is simply
!    E = Ekq+Epot
! The results are written to a file by the first node.
!
! Parameters
!     a : input matrix with the real part of the wavefunction
!     b : input matrix with the imaginary part of the wavefunction
!     t : number of time steps made
!     dt: time step
!
      USE fprecision
      USE commtypes
      USE grid
      USE hbar
      USE mpivars
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,n,ista:iend) :: a,b
      DOUBLE PRECISION    :: mass,ekq
      DOUBLE PRECISION    :: tmp,tmq
      REAL(KIND=GP), DIMENSION(n,n,ksta:kend) :: r1
      REAL(KIND=GP), INTENT(IN)               :: dt
      INTEGER, INTENT(IN) :: t
      INTEGER             :: i,j,k

!
! Computes the mass
!
      CALL variance(a,tmp,1)
      CALL variance(b,tmq,1)
      IF (myrank.eq.0) mass = tmp+tmq
!
! Computes the kinetic + quantum energy
!
      CALL variance(a,tmp,0)
      CALL variance(b,tmq,0)
      IF (myrank.eq.0) ekq = 2*alpha**2*(tmp+tmq)
!
! Computes the quartic energy and then the potential energy
!
      CALL squareabs(a,b,r1,1)
      tmp = 0.D0
!$omp parallel do if (kend-ksta.ge.nth) private (j,i) reduction(+:tmp)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i) reduction(+:tmp)
         DO j = 1,n
            DO i = 1,n
               tmp = tmp+r1(i,j,k)**2
            END DO
         END DO
      END DO
      CALL MPI_REDUCE(tmp,tmq,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
      IF (myrank.eq.0) THEN
         tmq = alpha*(beta*tmq/real(n,kind=GP)**9-2*omegag*mass &
              +omegag**2/beta)
      ENDIF
!
! Creates a external file to store the results
!
      IF (myrank.eq.0) THEN
         OPEN(1,file='balance.txt',position='append')
         WRITE(1,10) (t-1)*dt,mass,ekq,tmq
   10    FORMAT( E13.6,E22.14,E22.14,E22.14 )
         CLOSE(1)
      ENDIF

      RETURN
      END SUBROUTINE gpecheck

!*****************************************************************
      SUBROUTINE momentum(a,b,t,dt)
!-----------------------------------------------------------------
!
! Computes the three components of the total momentum
!    p = 2.alpha[zbar.grad(z)-z.grad(zbar)]
! The result is written to a file by the first node.
!
! Parameters
!     a : input matrix with the real part of the wavefunction
!     b : input matrix with the imaginary part of the wavefunction
!     t : number of time steps made
!     dt: time step
!
      USE fprecision
      USE grid
      USE hbar
      USE mpivars
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,n,ista:iend) :: a,b
      COMPLEX(KIND=GP), DIMENSION(n,n,ista:iend) :: C1,C2
      DOUBLE PRECISION    :: tmp,tmq
      DOUBLE PRECISION    :: jx,jy,jz
      REAL(KIND=GP), INTENT(IN)    :: dt
      INTEGER, INTENT(IN) :: t

      CALL derivk3(a,C1,1)
      CALL derivk3(b,C2,1)
      CALL product(a,C2,tmp)
      CALL product(b,C1,tmq)
      IF (myrank.eq.0) THEN
         jx = 2*alpha*(tmp-tmq)
      ENDIF
      CALL derivk3(a,C1,2)
      CALL derivk3(b,C2,2)
      CALL product(a,C2,tmp)
      CALL product(b,C1,tmq)
      IF (myrank.eq.0) THEN
         jy = 2*alpha*(tmp-tmq)
      ENDIF
      CALL derivk3(a,C1,3)
      CALL derivk3(b,C2,3)
      CALL product(a,C2,tmp)
      CALL product(b,C1,tmq)
      IF (myrank.eq.0) THEN
         jz = 2*alpha*(tmp-tmq)
         OPEN(1,file='momentum.txt',position='append')
         WRITE(1,FMT='(E13.6,E22.14,E22.14,E22.14)') (t-1)*dt,jx,jy,jz
         CLOSE(1)
      ENDIF

      RETURN
      END SUBROUTINE momentum

!*****************************************************************
      SUBROUTINE gpemassspec(a,b,nmb)
!-----------------------------------------------------------------
!
! Computes the spectrum of mass, which can be computed in 
! spectral space directly as mass is quadratic in the 
! wavefunction. The spectrum starts at k=0 to preserve 
! information of the total mass (k=0,1,...,N/2). The output 
! is written to a file by the first node.
!
! Parameters
!     a : real part of the wavefunction in Fourier space
!     b : imaginary part of the wavefunction in Fourier space
!     nmb: the extension used when writting the file
!
      USE fprecision
      USE commtypes
      USE kes
      USE grid
      USE mpivars
      USE filefmt
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,n,ista:iend) :: a,b
      DOUBLE PRECISION, DIMENSION(n/2+1)                :: Ek,Ektot
      DOUBLE PRECISION :: tmq
      REAL(KIND=GP)    :: rmp
      INTEGER          :: i,j,k
      INTEGER          :: kmn
      CHARACTER(len=*), INTENT(IN) :: nmb

!
! Sets Ek to zero
!
      DO i = 1,n/2+1
         Ek(i) = 0.0D0
      END DO
!
! Computes the power spectrum
!
      rmp = 1./real(n,kind=GP)**6
      IF (ista.eq.1) THEN
!$omp parallel do private (k,kmn,tmq)
         DO j = 1,n
            DO k = 1,n
               kmn = int(sqrt(ka2(k,j,1))+1.501)
               IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                  tmq = rmp*(abs(a(k,j,1))**2+abs(b(k,j,1))**2)
!$omp atomic
                  Ek(kmn) = Ek(kmn)+tmq
               ENDIF
            END DO
         END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k,kmn,tmq)
         DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k,kmn,tmq)
            DO j = 1,n
               DO k = 1,n
                  kmn = int(sqrt(ka2(k,j,i))+1.501)
                  IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                     tmq = 2*rmp*(abs(a(k,j,i))**2+abs(b(k,j,i))**2)
!$omp atomic
                     Ek(kmn) = Ek(kmn)+tmq
                  ENDIF
               END DO
            END DO
         END DO
      ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k,kmn,tmq)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,kmn,tmq)
            DO j = 1,n
               DO k = 1,n
                  kmn = int(sqrt(ka2(k,j,i))+1.501)
                  IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                     tmq = 2*rmp*(abs(a(k,j,i))**2+abs(b(k,j,i))**2)
!$omp atomic
                     Ek(kmn) = Ek(kmn)+tmq
                  ENDIF
               END DO
            END DO
         END DO
      ENDIF
!
! Computes the reduction between nodes
!
      CALL MPI_REDUCE(Ek,Ektot,n/2+1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
!
! Exports the spectrum to a file
!
      IF (myrank.eq.0) THEN
         OPEN(1,file='massspectrum.' // nmb // '.txt')
         WRITE(1,FMT='(E23.15)') Ektot
         CLOSE(1)
      ENDIF

      RETURN
      END SUBROUTINE gpemassspec

!*****************************************************************
      SUBROUTINE gperealspec(a,b,nmb)
!-----------------------------------------------------------------
!
! Computes the spectrum of kinetic, quantum, and potential (or 
! internal) energy. These quantities must be computed in real 
! space first, and then transformed to Fourier space to compute 
! the spectrum. The spectra start at k=0 to preserve information
! of the energy in the condensate, and are not dealiased 
! (k = 0,1,...,N/2). The output is written to files by the first
! node.
!
! Parameters
!     a : real part of the wavefunction in Fourier space
!     b : imaginary part of the wavefunction in Fourier space
!     nmb: the extension used when writting the file
!
      USE fprecision
      USE commtypes
      USE kes
      USE fft
      USE ali
      USE grid
      USE hbar
      USE mpivars
      USE filefmt
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,n,ista:iend) :: a,b
      COMPLEX(KIND=GP), DIMENSION(n,n,ista:iend) :: c1,c2,c3,c4
      DOUBLE PRECISION, DIMENSION(n/2+1)         :: Ek,Ektot,Ec,Ectot
      REAL(KIND=GP), DIMENSION(n,n,ksta:kend)    :: r1,r2,r3
      REAL(KIND=GP), DIMENSION(n,n,ksta:kend)    :: qua,kin
      REAL(KIND=GP)    :: rmp,rmq
      INTEGER          :: i,j,k
      INTEGER          :: kmn
      CHARACTER(len=*), INTENT(IN) :: nmb

!
! Transforms the wavefunction to real space
!
!$omp parallel do if (iend-2.ge.nth) private (j,k,kmn,tmq)
      DO i = ista,iend
!$omp parallel do if (iend-2.lt.nth) private (k,kmn,tmq)
         DO j = 1,n
            DO k = 1,n
               c1(k,j,i) = a(k,j,i)
               c2(k,j,i) = b(k,j,i)
            END DO
         END DO
      END DO
      CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,c2,r2,MPI_COMM_WORLD)
!
! Computes and writes the internal energy spectrum
! Eint = 2.alpha.beta.(|z|^2-rho0)^2
!
      rmp = sqrt(alpha*beta)*omegag/beta
      rmq = sqrt(alpha*beta)/real(n,kind=GP)**6
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
         DO j = 1,n
            DO i = 1,n
               r3(i,j,k) = (r1(i,j,k)**2+r2(i,j,k)**2)*rmq-rmp
            END DO
         END DO
      END DO
      CALL fftp3d_real_to_complex(planrc,r3,c1,MPI_COMM_WORLD)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend   ! This spectrum must be dealiased
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,n
            DO k = 1,n
               IF (ka2(k,j,i).gt.kmax) THEN
                  c1(k,j,i) = 0.0_GP
               ENDIF
            END DO
         END DO
      END DO
      CALL spectrscc(c1,Ek,1.0_GP)
      IF (myrank.eq.0) THEN
         OPEN(1,file='intspectrum.' // nmb // '.txt')
         WRITE(1,FMT='(E23.15)') Ek
         CLOSE(1)
      ENDIF
!
! Computes z/sqrt(|z|^2) inplace
!
      CALL zturn(r1,r2)
!
! Computes and writes the quantum energy spectrum, and 
! prepares to compute the kinetic energy spectra
! Equa ~ (zturnre*grad(zre)+zturnim*grad(zim))^2
! Ekin ~ (zturnre*grad(zim)-zturnim*grad(zre))^2
!
      CALL derivk3(a,c1,1)   ! x component
      CALL fftp3d_complex_to_real(plancr,c1,r3,MPI_COMM_WORLD)
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
         DO j = 1,n
            DO i = 1,n
               qua(i,j,k) =  r1(i,j,k)*r3(i,j,k)
               kin(i,j,k) = -r2(i,j,k)*r3(i,j,k)
            END DO
         END DO
      END DO
      CALL derivk3(b,c1,1)
      CALL fftp3d_complex_to_real(plancr,c1,r3,MPI_COMM_WORLD)
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
         DO j = 1,n
            DO i = 1,n
               qua(i,j,k) = qua(i,j,k)+r2(i,j,k)*r3(i,j,k)
               kin(i,j,k) = kin(i,j,k)+r1(i,j,k)*r3(i,j,k)
            END DO
         END DO
      END DO
      CALL fftp3d_real_to_complex(planrc,qua,c1,MPI_COMM_WORLD)
      CALL fftp3d_real_to_complex(planrc,kin,c2,MPI_COMM_WORLD)
      CALL spectrscc(c1,Ektot,1.0_GP)

      CALL derivk3(a,c1,2)   ! y component
      CALL fftp3d_complex_to_real(plancr,c1,r3,MPI_COMM_WORLD)
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
         DO j = 1,n
            DO i = 1,n
               qua(i,j,k) =  r1(i,j,k)*r3(i,j,k)
               kin(i,j,k) = -r2(i,j,k)*r3(i,j,k)
            END DO
         END DO
      END DO
      CALL derivk3(b,c1,2)
      CALL fftp3d_complex_to_real(plancr,c1,r3,MPI_COMM_WORLD)
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
         DO j = 1,n
            DO i = 1,n
               qua(i,j,k) = qua(i,j,k)+r2(i,j,k)*r3(i,j,k)
               kin(i,j,k) = kin(i,j,k)+r1(i,j,k)*r3(i,j,k)
            END DO
         END DO
      END DO
      CALL fftp3d_real_to_complex(planrc,qua,c1,MPI_COMM_WORLD)
      CALL fftp3d_real_to_complex(planrc,kin,c3,MPI_COMM_WORLD)
      CALL spectrscc(c1,Ek,1.0_GP)
      IF (myrank.eq.0) Ektot = Ektot+Ek

      CALL derivk3(a,c1,3)   ! z component
      CALL fftp3d_complex_to_real(plancr,c1,r3,MPI_COMM_WORLD)
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
         DO j = 1,n
            DO i = 1,n
               qua(i,j,k) =  r1(i,j,k)*r3(i,j,k)
               kin(i,j,k) = -r2(i,j,k)*r3(i,j,k)
            END DO
         END DO
      END DO
      CALL derivk3(b,c1,3)
      CALL fftp3d_complex_to_real(plancr,c1,r3,MPI_COMM_WORLD)
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
         DO j = 1,n
            DO i = 1,n
               qua(i,j,k) = qua(i,j,k)+r2(i,j,k)*r3(i,j,k)
               kin(i,j,k) = kin(i,j,k)+r1(i,j,k)*r3(i,j,k)
            END DO
         END DO
      END DO
      CALL fftp3d_real_to_complex(planrc,qua,c1,MPI_COMM_WORLD)
      CALL fftp3d_real_to_complex(planrc,kin,c4,MPI_COMM_WORLD)
      CALL spectrscc(c1,Ek,1.0_GP)
      IF (myrank.eq.0) THEN
         rmq = 2*alpha**2/real(n,kind=GP)**6
         Ektot = (Ektot+Ek)*rmq
         OPEN(1,file='qspectrum.' // nmb // '.txt')
         WRITE(1,FMT='(E23.15)') Ektot
         CLOSE(1)
      ENDIF
!
! Computes the compressible and incompressible kinetic energy spectra
!
      CALL gauge3(c2,c3,c4,c1,1)      ! x component 
      CALL spectrscc(c1,Ektot,1.0_GP) ! incompressible
      c1 = c2-c1
      CALL spectrscc(c1,Ectot,1.0_GP) ! compressible

      CALL gauge3(c2,c3,c4,c1,2)      ! y component
      CALL spectrscc(c1,Ek,1.0_GP)    ! incompressible
      c1 = c3-c1
      CALL spectrscc(c1,Ec,1.0_GP)    ! compressible
      IF (myrank.eq.0) THEN
         Ektot = Ektot+Ek
         Ectot = Ectot+Ec
      ENDIF

      CALL gauge3(c2,c3,c4,c1,3)      ! z component
      CALL spectrscc(c1,Ek,1.0_GP)    ! incompressible
      c1 = c4-c1
      CALL spectrscc(c1,Ec,1.0_GP)    ! compressible
      IF (myrank.eq.0) THEN
         rmq = 2*alpha**2/real(n,kind=GP)**6
         Ektot = (Ektot+Ek)*rmq
         Ectot = (Ectot+Ec)*rmq
         OPEN(1,file='kincspectrum.' // nmb // '.txt')
         WRITE(1,FMT='(E23.15)') Ektot
         CLOSE(1)
         OPEN(1,file='kcomspectrum.' // nmb // '.txt')
         WRITE(1,FMT='(E23.15)') Ectot
         CLOSE(1)
      ENDIF

      RETURN
      END SUBROUTINE gperealspec

!*****************************************************************
      SUBROUTINE zturn(ra,rb)
!-----------------------------------------------------------------
!
! Computes the real and imaginary parts of z/sqrt(|z|^2) in 
! place (i.e., the input is destroyed, and replaced by the 
! output). Note zturn in GHOST is the complex conjugate of the
! quantity called zturn in Brachet's GPE codes (including TYGRES).
!
! Parameters
!     a : real part of the wavefunction in real space
!     b : imaginary part of the wavefunction in real space
!
      USE fprecision
      USE grid
      USE hbar
      USE mpivars
!$    USE threads
      IMPLICIT NONE

      REAL(KIND=GP), INTENT(INOUT), DIMENSION(n,n,ksta:kend) :: ra,rb
      REAL(KIND=GP)    :: rmp
      INTEGER          :: i,j,k

!
! Computes zbar/sqrt(|z|^2)
!
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
         DO j = 1,n
            DO i = 1,n
               rmp = 1.0_GP/sqrt(ra(i,j,k)**2+rb(i,j,k)**2+ &
                     real(n,kind=GP)**6*regu*omegag/beta)
               ra(i,j,k) = ra(i,j,k)*rmp
               rb(i,j,k) = rb(i,j,k)*rmp
            END DO
         END DO
      END DO

      RETURN
      END SUBROUTINE zturn
