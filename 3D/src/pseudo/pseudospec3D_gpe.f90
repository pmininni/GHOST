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
!
! 17 May 2018: Support for elongated box (N.Muller & P.D.Mininni) 
!=================================================================

!*****************************************************************
      SUBROUTINE squareabs(a,b,r,dealias)
!-----------------------------------------------------------------
!
! Pointwise squared absolute value of a complex wavefunction Z.
! Note that when dealiased, the output is not normalized (i.e., 
! not divided by N^3).
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

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend)  :: a,b
      COMPLEX(KIND=GP), DIMENSION(nz,ny,ista:iend) :: c1
      REAL(KIND=GP), INTENT(OUT), DIMENSION(nx,ny,ksta:kend)    :: r
      REAL(KIND=GP), DIMENSION(nx,ny,ksta:kend)    :: r1
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
         DO j = 1,ny
            DO i = 1,nx
               r(i,j,k) = r1(i,j,k)**2
            END DO
         END DO
      END DO
!
! Computes the square of the imaginary part of the wavefunction
!
      c1 = b
      CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      rmp = 1.0_GP/ &
            (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
         DO j = 1,ny
            DO i = 1,nx
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
            DO j = 1,ny
               DO k = 1,nz
                  IF (kn2(k,j,i).gt.kmax) THEN
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

      COMPLEX(KIND=GP), INTENT(IN),  DIMENSION(nz,ny,ista:iend) :: a
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(nz,ny,ista:iend) :: b
      REAL(KIND=GP), INTENT(IN), DIMENSION(nx,ny,ksta:kend)     :: r
      REAL(KIND=GP), DIMENSION(nx,ny,ksta:kend)    :: r1
      REAL(KIND=GP)    :: rmp
      INTEGER :: i,j,k

!
! Computes Z.|Z|^2
!
      b = a
      CALL fftp3d_complex_to_real(plancr,b,r1,MPI_COMM_WORLD)
      rmp = 1.0_GP/ &
            (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
         DO j = 1,ny
            DO i = 1,nx
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
!    Ekq    = 2.alpha^2.|grad(z)|^2
!    Equart = alpha.beta.|z|^4
! The potential energy then is
!    Epot = Equart-2*alpha.omegag.mass+alpha.omegag^2/beta
! This quantity should be zero in the condensate. The total 
! energy is simply
!    E = Ekq+Epot
! The results are written to a file by the first node.
!
! Output files contain:
! 'balance.txt': time, mass, kinetic+quantum energy, quartic energy
!   [Ekq = 2.alpha^2.|grad(z)|^2, Equart = alpha.beta.|z|^4, and the      ]
!   [pot. energy is Epot = Equart-2*alpha.omegag.mass+alpha.omegag^2/beta.]
!   [Note this output replaces all 'balance.txt' files in quantum solvers.]
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
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a,b
      DOUBLE PRECISION    :: mass,ekq
      DOUBLE PRECISION    :: tmp,tmq
      REAL(KIND=GP), DIMENSION(nx,ny,ksta:kend) :: r1
      REAL(KIND=GP), INTENT(IN)                 :: dt
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
         DO j = 1,ny
            DO i = 1,nx
               tmp = tmp+r1(i,j,k)**2
            END DO
         END DO
      END DO
      CALL MPI_REDUCE(tmp,tmq,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
      IF (myrank.eq.0) THEN
         tmq = alpha*(beta*tmq/ &
          (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**3 - &
          2*omegag*mass + omegag**2/beta)
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
! Output files contain:
! 'momentum.txt': time, momentum_x, momentum_y, momentum_z
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

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a,b
      COMPLEX(KIND=GP), DIMENSION(nz,ny,ista:iend) :: C1,C2
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
! Output files contain:
! 'massspectrum.XXX.txt': k, mass(k)
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
      USE boxsize
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a,b
      DOUBLE PRECISION, DIMENSION(nmax/2+1)               :: Ek,Ektot
      DOUBLE PRECISION :: tmq
      REAL(KIND=GP)    :: rmp
      INTEGER          :: i,j,k
      INTEGER          :: kmn
      CHARACTER(len=*), INTENT(IN) :: nmb

!
! Sets Ek to zero
!
      DO i = 1,nmax/2+1
         Ek(i) = 0.0D0
      END DO
!
! Computes the power spectrum
!
      rmp = 1./ &
            (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2
      IF (ista.eq.1) THEN
!$omp parallel do private (k,kmn,tmq)
         DO j = 1,ny
            DO k = 1,nz
               kmn = int(sqrt(kk2(k,j,1))/Dkk+1.501)
               IF ((kmn.gt.0).and.(kmn.le.nmax/2+1)) THEN
                  tmq = rmp*(abs(a(k,j,1))**2+abs(b(k,j,1))**2)
!$omp atomic
                  Ek(kmn) = Ek(kmn)+tmq
               ENDIF
            END DO
         END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k,kmn,tmq)
         DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k,kmn,tmq)
            DO j = 1,ny
               DO k = 1,nz
                  kmn = int(sqrt(kk2(k,j,i))/Dkk+1.501)
                  IF ((kmn.gt.0).and.(kmn.le.nmax/2+1)) THEN
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
            DO j = 1,ny
               DO k = 1,nz
                  kmn = int(sqrt(kk2(k,j,i))/Dkk+1.501)
                  IF ((kmn.gt.0).and.(kmn.le.nmax/2+1)) THEN
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
      CALL MPI_REDUCE(Ek,Ektot,nmax/2+1,MPI_DOUBLE_PRECISION,MPI_SUM,0,&
                      MPI_COMM_WORLD,ierr)
!
! Exports the spectrum to a file
!
      IF (myrank.eq.0) THEN
         OPEN(1,file='massspectrum.' // nmb // '.txt')
         DO i=1,nmax/2+1
            WRITE(1,FMT='(E13.6,E23.15)') Dkk*(i-1),Ektot(i)/Dkk
         END DO
         CLOSE(1)
      ENDIF
! Igual que en HD

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
! Output files contain:
! 'intspectrum.XXX.txt' : k, Eint(k) [Eint = 2.alpha.beta.(|z|^2-rho0)^2]
! 'qspectrum.XXX.txt'   : k, Equa(k)
!   [Equa ~ (zturnre*grad(zre)+zturnim*grad(zim))^2]
! 'kincspectrum.XXX.txt': k, Einc(k)
! 'kcomspectrum.XXX.txt': k, Ecom(k)
!   [Ekin ~ (zturnre*grad(zim)-zturnim*grad(zre))^2, decomposed into]
!   [incompressible and compressible components.]
!
! Parameters
!     a : real part of the wavefunction in Fourier space
!     b : imaginary part of the wavefunction in Fourier space
!     nmb: the extension used when writting the file
!
      USE kes
      USE grid
      USE mpivars
      USE filefmt
      USE boxsize
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a,b
      DOUBLE PRECISION, DIMENSION(nmax/2+1)        :: Eint,Equa,Einc,Ecom
      INTEGER                      :: i
      CHARACTER(len=*), INTENT(IN) :: nmb

!
! Computes all the energy  spectra
!
      CALL gperealspecc(a,b,Eint,Equa,Einc,Ecom)
!
! Exports the energy spectrum to a file
!
      IF (myrank.eq.0) THEN
         OPEN(1,file='intspectrum.' // nmb // '.txt')
         DO i=1,nmax/2+1
            WRITE(1,FMT='(E13.6,E23.15)') Dkk*(i-1),Eint(i)/Dkk
         END DO
         CLOSE(1)
         OPEN(1,file='qspectrum.' // nmb // '.txt')
         DO i=1,nmax/2+1                                       
            WRITE(1,FMT='(E13.6,E23.15)') Dkk*(i-1),Equa(i)/Dkk   
         END DO  
         CLOSE(1)
         OPEN(1,file='kincspectrum.' // nmb // '.txt')
         DO i=1,nmax/2+1
            WRITE(1,FMT='(E13.6,E23.15)') Dkk*(i-1),Einc(i)/Dkk
         END DO
         CLOSE(1)
         OPEN(1,file='kcomspectrum.' // nmb // '.txt')
         DO i=1,nmax/2+1
            WRITE(1,FMT='(E13.6,E23.15)') Dkk*(i-1),Ecom(i)/Dkk
         END DO
         CLOSE(1)
      ENDIF

      RETURN
      END SUBROUTINE gperealspec
        
!*****************************************************************
      SUBROUTINE gperealspecc(a,b,Eint,Equa,Einc,Ecom)
!-----------------------------------------------------------------
!
! Computes the spectrum of kinetic, quantum, and potential (or 
! internal) energy, returning them.
!
! Parameters
!     a   : real part of the wavefunction in Fourier space
!     b   : imaginary part of the wavefunction in Fourier space
!     Eint: at the output contains the internal energy spectrum
!     Equa: at the output contains the quantum energy spectrum
!     Einc: at the output contains the incompressible energy spec.
!     Ecom: at the output contains the compressible energy spec.
!
      USE fprecision
      USE commtypes
      USE kes
      USE fft
      USE ali
      USE grid
      USE hbar
      USE mpivars
      USE boxsize
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a,b
      COMPLEX(KIND=GP), DIMENSION(nz,ny,ista:iend) :: c1,c2,c3,c4
      DOUBLE PRECISION, INTENT(OUT), DIMENSION(nmax/2+1) :: Eint,Equa,Einc,Ecom
      DOUBLE PRECISION, DIMENSION(nmax/2+1)        :: Ek1,Ek2
      REAL(KIND=GP), DIMENSION(nx,ny,ksta:kend)    :: r1,r2,r3
      REAL(KIND=GP), DIMENSION(nx,ny,ksta:kend)    :: qua,kin
      REAL(KIND=GP)    :: rmp,rmq
      INTEGER          :: i,j,k
      INTEGER          :: kmn

!
! Transforms the wavefunction to real space
!
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
            DO k = 1,nz
               c1(k,j,i) = a(k,j,i)
               c2(k,j,i) = b(k,j,i)
            END DO
         END DO
      END DO
      CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,c2,r2,MPI_COMM_WORLD)
!
! Computes the internal energy spectrum
! Eint = 2.alpha.beta.(|z|^2-rho0)^2
!
      rmp = sqrt(alpha*beta)*omegag/beta
      rmq = sqrt(alpha*beta)/ &
            (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
         DO j = 1,ny
            DO i = 1,nx
               r3(i,j,k) = (r1(i,j,k)**2+r2(i,j,k)**2)*rmq-rmp
            END DO
         END DO
      END DO
      CALL fftp3d_real_to_complex(planrc,r3,c1,MPI_COMM_WORLD)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend   ! This spectrum must be dealiased
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
            DO k = 1,nz
               IF (kn2(k,j,i).gt.kmax) THEN
                  c1(k,j,i) = 0.0_GP
               ENDIF
            END DO
         END DO
      END DO
      CALL spectrscc(c1,Eint,1.0_GP)
!
! Computes z/sqrt(|z|^2) inplace
!
      CALL zturn(r1,r2)
!
! Computes the quantum energy spectrum, and 
! prepares to compute the kinetic energy spectra
! Equa ~ (zturnre*grad(zre)+zturnim*grad(zim))^2
! Ekin ~ (zturnre*grad(zim)-zturnim*grad(zre))^2
!
      CALL derivk3(a,c1,1)   ! x component
      CALL fftp3d_complex_to_real(plancr,c1,r3,MPI_COMM_WORLD)
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
         DO j = 1,ny
            DO i = 1,nx
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
         DO j = 1,ny
            DO i = 1,nx
               qua(i,j,k) = qua(i,j,k)+r2(i,j,k)*r3(i,j,k)
               kin(i,j,k) = kin(i,j,k)+r1(i,j,k)*r3(i,j,k)
            END DO
         END DO
      END DO
      CALL fftp3d_real_to_complex(planrc,qua,c1,MPI_COMM_WORLD)
      CALL fftp3d_real_to_complex(planrc,kin,c2,MPI_COMM_WORLD)
      CALL spectrscc(c1,Ek1,1.0_GP)

      CALL derivk3(a,c1,2)   ! y component
      CALL fftp3d_complex_to_real(plancr,c1,r3,MPI_COMM_WORLD)
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
         DO j = 1,ny
            DO i = 1,nx
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
         DO j = 1,ny
            DO i = 1,nx
               qua(i,j,k) = qua(i,j,k)+r2(i,j,k)*r3(i,j,k)
               kin(i,j,k) = kin(i,j,k)+r1(i,j,k)*r3(i,j,k)
            END DO
         END DO
      END DO
      CALL fftp3d_real_to_complex(planrc,qua,c1,MPI_COMM_WORLD)
      CALL fftp3d_real_to_complex(planrc,kin,c3,MPI_COMM_WORLD)
      CALL spectrscc(c1,Ek2,1.0_GP)
      Ek1 = Ek1+Ek2

      CALL derivk3(a,c1,3)   ! z component
      CALL fftp3d_complex_to_real(plancr,c1,r3,MPI_COMM_WORLD)
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
         DO j = 1,ny
            DO i = 1,nx
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
         DO j = 1,ny
            DO i = 1,nx
               qua(i,j,k) = qua(i,j,k)+r2(i,j,k)*r3(i,j,k)
               kin(i,j,k) = kin(i,j,k)+r1(i,j,k)*r3(i,j,k)
            END DO
         END DO
      END DO
      CALL fftp3d_real_to_complex(planrc,qua,c1,MPI_COMM_WORLD)
      CALL fftp3d_real_to_complex(planrc,kin,c4,MPI_COMM_WORLD)
      CALL spectrscc(c1,Ek2,1.0_GP)
      rmq = 2*alpha**2/ &
            (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2
      Equa = (Ek1+Ek2)*rmq
!
! Computes the compressible and incompressible kinetic energy spectra
!
      CALL gauge3(c2,c3,c4,c1,1)      ! x component 
      CALL spectrscc(c1,Einc,1.0_GP) ! incompressible
      c1 = c2-c1
      CALL spectrscc(c1,Ecom,1.0_GP) ! compressible

      CALL gauge3(c2,c3,c4,c1,2)      ! y component
      CALL spectrscc(c1,Ek1,1.0_GP)    ! incompressible
      c1 = c3-c1
      CALL spectrscc(c1,Ek2,1.0_GP)    ! compressible
      Einc = Einc+Ek1
      Ecom = Ecom+Ek2

      CALL gauge3(c2,c3,c4,c1,3)      ! z component
      CALL spectrscc(c1,Ek1,1.0_GP)    ! incompressible
      c1 = c4-c1
      CALL spectrscc(c1,Ek2,1.0_GP)    ! compressible
      rmq = 2*alpha**2/ &
            (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2
      Einc = (Einc+Ek1)*rmq
      Ecom = (Ecom+Ek2)*rmq

      RETURN
      END SUBROUTINE gperealspecc

!*****************************************************************
      SUBROUTINE gperealtrans(dt,io,qo,ko,co,in,qn,kn,cn,nmb)
!-----------------------------------------------------------------
!
! Computes the energy transfers in Fourier space for the GPE
! equations in 3D. Normalization of the transfer function is such
! that the fluxes are Pi = -sum[T(k).Dkk], where Dkk is the width
! of the Fourier shells. The output is written to files by the
! first node.
!
! Output files contain:
! 'inttransfer.XXX.txt' : k, Ti(k) (internal energy transfer function)
! 'qtransfer.XXX.txt'   : k, Tq(k) (quantum energy transfer)
! 'kinctransfer.XXX.txt': k, Tk(k) (incompressible kin. energy transfer)
! 'kcomtransfer.XXX.txt': k, Tc(k) (compressible kin. energy transfer)
!   [Each transfer function is computed as to T_x(k) = dE_x(k)/dt]
!
! Parameters
!     dt : time step
!     io : spectrum of internal energy at t-dt
!     qo : spectrum of quantum energy at t-dt
!     ko : spectrum of incompressible kin. energy at t-dt
!     co : spectrum of compressible kin. energy at t-dt
!     in : spectrum of internal energy at t
!     qn : spectrum of quantum energy at t
!     kn : spectrum of incompressible kin. energy at t
!     cn : spectrum of compressible kin. energy at t-dt
!     nmb: nmb: the extension used when writting the file
!
      USE kes
      USE grid
      USE mpivars
      USE filefmt
      USE boxsize
      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN), DIMENSION(nmax/2+1)    :: io,qo,ko,co
      DOUBLE PRECISION, INTENT(INOUT), DIMENSION(nmax/2+1) :: in,qn,kn,cn
      REAL(KIND=GP),INTENT(IN)     :: dt
      INTEGER                      :: i
      CHARACTER(len=*), INTENT(IN) :: nmb

!
! Computes time derivatives
!
      IF (myrank.eq.0) THEN
         DO i=1,nmax/2+1
            in(i) = (in(i)-io(i))/dt
            qn(i) = (qn(i)-qo(i))/dt
            kn(i) = (kn(i)-ko(i))/dt
            cn(i) = (cn(i)-co(i))/dt
         END DO
!
! Exports the transfer functions to files
!
         OPEN(1,file='inttransfer.' // nmb // '.txt')
         DO i=1,nmax/2+1
            WRITE(1,FMT='(E13.6,E23.15)') Dkk*(i-1),in(i)/Dkk
         END DO
         CLOSE(1)
         OPEN(1,file='qtransfer.' // nmb // '.txt')
         DO i=1,nmax/2+1                                       
            WRITE(1,FMT='(E13.6,E23.15)') Dkk*(i-1),qn(i)/Dkk   
         END DO  
         CLOSE(1)
         OPEN(1,file='kinctransfer.' // nmb // '.txt')
         DO i=1,nmax/2+1
            WRITE(1,FMT='(E13.6,E23.15)') Dkk*(i-1),kn(i)/Dkk
         END DO
         CLOSE(1)
         OPEN(1,file='kcomtransfer.' // nmb // '.txt')
         DO i=1,nmax/2+1
            WRITE(1,FMT='(E13.6,E23.15)') Dkk*(i-1),cn(i)/Dkk
         END DO
         CLOSE(1)
      ENDIF

      RETURN
      END SUBROUTINE gperealtrans

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

      REAL(KIND=GP), INTENT(INOUT), DIMENSION(nx,ny,ksta:kend) :: ra,rb
      REAL(KIND=GP)    :: rmp
      INTEGER          :: i,j,k

!
! Computes zbar/sqrt(|z|^2)
!
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
         DO j = 1,ny
            DO i = 1,nx
               rmp = 1.0_GP/sqrt(ra(i,j,k)**2+rb(i,j,k)**2+ &
                     (real(nx,kind=GP)*real(ny,kind=GP)*    &
                     real(nz,kind=GP))**2*regu*omegag/beta)
               ra(i,j,k) = ra(i,j,k)*rmp
               rb(i,j,k) = rb(i,j,k)*rmp
            END DO
         END DO
      END DO

      RETURN
      END SUBROUTINE zturn

!*******************************************************************
      SUBROUTINE gpekfield(a,b,c,d,e,f,g,h)
!-------------------------------------------------------------------
!
! Computes all the components of the compressible and incompressible
! parts of sqrt(rho)*v. These quantities must be computed in real
! space and then transformed to Fourier space.
!
! Parameters
!     a : real part of the wavefunction in Fourier space
!     b : imaginary part of the wavefunction in Fourier space
!     c : x-component of the incompressible part [output]
!     d : x-component of the compressible part   [output]
!     e : y-component of the incompressible part [output]
!     f : y-component of the compressible part   [output]
!     g : z-component of the incompressible part [output]
!     h : z-component of the compressible part   [output]
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

      COMPLEX(KIND=GP), INTENT(IN),  DIMENSION(nz,ny,ista:iend) :: a,b
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(nz,ny,ista:iend) :: c,d
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(nz,ny,ista:iend) :: e,f
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(nz,ny,ista:iend) :: g,h
      COMPLEX(KIND=GP), DIMENSION(nz,ny,ista:iend) :: c1,c2,c3,c4
      DOUBLE PRECISION, DIMENSION(nmax/2+1)        :: Ek,Ektot,Ec,Ectot
      REAL(KIND=GP), DIMENSION(nx,ny,ksta:kend)    :: r1,r2,r3
      REAL(KIND=GP), DIMENSION(nx,ny,ksta:kend)    :: qua,kin
      REAL(KIND=GP)    :: rmp,rmq
      INTEGER          :: i,j,k
      INTEGER          :: kmn

!
! Transforms the wavefunction to real space
!
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
            DO k = 1,nz
               c1(k,j,i) = a(k,j,i)
               c2(k,j,i) = b(k,j,i)
            END DO
         END DO
      END DO
      CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,c2,r2,MPI_COMM_WORLD)
!
! Computes z/sqrt(|z|^2) inplace
!
      CALL zturn(r1,r2)
!
! Prepares to compute the kinetic energy spectra
! Ekin ~ (zturnre*grad(zim)-zturnim*grad(zre))^2
!
      CALL derivk3(a,c1,1)   ! x component
      CALL fftp3d_complex_to_real(plancr,c1,r3,MPI_COMM_WORLD)
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
         DO j = 1,ny
            DO i = 1,nx
               kin(i,j,k) = -r2(i,j,k)*r3(i,j,k)
            END DO
         END DO
      END DO
      CALL derivk3(b,c1,1)
      CALL fftp3d_complex_to_real(plancr,c1,r3,MPI_COMM_WORLD)
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
         DO j = 1,ny
            DO i = 1,nx
               kin(i,j,k) = kin(i,j,k)+r1(i,j,k)*r3(i,j,k)
            END DO
         END DO
      END DO
      CALL fftp3d_real_to_complex(planrc,kin,c2,MPI_COMM_WORLD)

      CALL derivk3(a,c1,2)   ! y component
      CALL fftp3d_complex_to_real(plancr,c1,r3,MPI_COMM_WORLD)
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
         DO j = 1,ny
            DO i = 1,nx
               kin(i,j,k) = -r2(i,j,k)*r3(i,j,k)
            END DO
         END DO
      END DO
      CALL derivk3(b,c1,2)
      CALL fftp3d_complex_to_real(plancr,c1,r3,MPI_COMM_WORLD)
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
         DO j = 1,ny
            DO i = 1,nx
               kin(i,j,k) = kin(i,j,k)+r1(i,j,k)*r3(i,j,k)
            END DO
         END DO
      END DO
      CALL fftp3d_real_to_complex(planrc,kin,c3,MPI_COMM_WORLD)

      CALL derivk3(a,c1,3)   ! z component
      CALL fftp3d_complex_to_real(plancr,c1,r3,MPI_COMM_WORLD)
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
         DO j = 1,ny
            DO i = 1,nx
               kin(i,j,k) = -r2(i,j,k)*r3(i,j,k)
            END DO
         END DO
      END DO
      CALL derivk3(b,c1,3)
      CALL fftp3d_complex_to_real(plancr,c1,r3,MPI_COMM_WORLD)
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
         DO j = 1,ny
            DO i = 1,nx
               kin(i,j,k) = kin(i,j,k)+r1(i,j,k)*r3(i,j,k)
            END DO
         END DO
      END DO
      CALL fftp3d_real_to_complex(planrc,kin,c4,MPI_COMM_WORLD)
!
! Computes the compressible and incompressible kinetic energy spectra
!
      CALL gauge3(c2,c3,c4,c1,1)      ! x component 
      c  = c1    ! incompressible
      c1 = c2-c1
      d  = c1    ! compressible

      CALL gauge3(c2,c3,c4,c1,2)      ! y component
      e = c1     ! incompressible
      c1 = c3-c1
      f = c1     ! compressible

      CALL gauge3(c2,c3,c4,c1,3)      ! z component
      g = c1     ! incompressible
      c1 = c4-c1
      h = c1     ! compressible

      RETURN
      END SUBROUTINE gpekfield


!***********************************************************************
      SUBROUTINE gpehelicity(a,b,t,dt)
!-----------------------------------------------------------------------
!
! Computes two measurements of helicity, based on a regularized parallel
! velocity and on v (usual definition of the helicity for a classical
! fluid). The output is written to a file by the first node.
!
! Output files contain:
! 'helicity.txt': time, regularized helicity, classical helicity
!
! Parameters
!     a    : real part of the wavefunction in Fourier space
!     b    : imaginary part of the wavefunction in Fourier space
!     t    : number of time steps made
!     dt   : time step
!
      USE fprecision
      USE commtypes
      USE kes
      USE fft
      USE ali
      USE grid
      USE hbar
      USE mpivars
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a,b
      COMPLEX(KIND=GP), DIMENSION(nz,ny,ista:iend) :: c1,c2,c3
      COMPLEX(KIND=GP), DIMENSION(nz,ny,ista:iend) :: c4,c5,c6
      COMPLEX(KIND=GP), DIMENSION(nz,ny,ista:iend) :: c7,c8,c9
      REAL(KIND=GP), DIMENSION(nx,ny,ksta:kend)    :: r1,r2,r3,r4
      REAL(KIND=GP), DIMENSION(nx,ny,ksta:kend)    :: r5,r6,r7,r8
      REAL(KIND=GP), DIMENSION(nx,ny,ksta:kend)    :: r9,r10,r11,r12
      REAL(KIND=GP), INTENT(IN)          :: dt
      DOUBLE PRECISION    :: Htot1,Htot2
      INTEGER, INTENT(IN) :: t
      REAL(KIND=GP)       :: rmp,rmq
      REAL(KIND=GP)       :: tmp
      INTEGER             :: i,j,k
      INTEGER             :: kmn

!
! Transforms the wavefunction to real space
!
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
            DO k = 1,nz
               c1(k,j,i) = a(k,j,i)
               c2(k,j,i) = b(k,j,i)
            END DO
         END DO
      END DO
      CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,c2,r2,MPI_COMM_WORLD)
!
! Obtains e_parallel by calculating the cross product of
! grad(zre) and grad(zim). Then it normalizes.
!
      CALL derivk3(a,c1,1)
      CALL fftp3d_complex_to_real(plancr,c1,r3,MPI_COMM_WORLD)
      CALL derivk3(b,c1,1)
      CALL fftp3d_complex_to_real(plancr,c1,r4,MPI_COMM_WORLD)
      CALL derivk3(a,c1,2)
      CALL fftp3d_complex_to_real(plancr,c1,r5,MPI_COMM_WORLD)
      CALL derivk3(b,c1,2)
      CALL fftp3d_complex_to_real(plancr,c1,r6,MPI_COMM_WORLD)
      CALL derivk3(a,c1,3)  
      CALL fftp3d_complex_to_real(plancr,c1,r7,MPI_COMM_WORLD)
      CALL derivk3(b,c1,3)
      CALL fftp3d_complex_to_real(plancr,c1,r8,MPI_COMM_WORLD)
      tmp = 1.0_GP/ &
            (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
         DO j = 1,ny
            DO i = 1,nx
               r9(i,j,k)  = (r5(i,j,k)*r8(i,j,k)-r7(i,j,k)*r6(i,j,k))* &
                            tmp
               r10(i,j,k) = (r7(i,j,k)*r4(i,j,k)-r3(i,j,k)*r8(i,j,k))* &
                            tmp
               r11(i,j,k) = (r3(i,j,k)*r6(i,j,k)-r5(i,j,k)*r4(i,j,k))* &
                            tmp
               r12(i,j,k) = (real(nx,kind=GP)*real(ny,kind=GP)  *      &
                            real(nz,kind=GP))**2/(r3(i,j,k)**2  +      &
                            r4(i,j,k)**2 + r5(i,j,k)**2         +      &
                            r6(i,j,k)**2 + r7(i,j,k)**2         +      &
                            r8(i,j,k)**2 + (real(nx,kind=GP)*          &
                            real(ny,kind=GP)*real(nz,kind=GP))**2*     &
                            regu*omegag/beta)
            END DO
         END DO
      END DO

!
! Calculates v_parallel by doing regularizing and the projecting
! The magnitud calculated is:
! 2*alpha e_parallel.((grad(zre).grad)(grad(zim)) -
! (grad(zim).grad)(grad(zre)))/(grad(zre)**2+grad(zim)**2)
! 
      CALL derivk3(b,c1,1)
      CALL derivk3(b,c2,2)
      CALL derivk3(b,c3,3)
!
! Computes d_x(zre) e_parallel.d_x(grad(zim))
!
      CALL derivk3(a,c7,1)
      CALL derivk3(c1,c4,1)
      CALL derivk3(c2,c5,1)
      CALL derivk3(c3,c6,1)
      CALL fftp3d_complex_to_real(plancr,c7,r3,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,c4,r5,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,c5,r6,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,c6,r7,MPI_COMM_WORLD)
      r3 = r3*tmp
      r5 = r5*tmp
      r6 = r6*tmp
      r7 = r7*tmp

!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
         DO j = 1,ny
            DO i = 1,nx
               r4(i,j,k) = r3(i,j,k)*(r9(i,j,k)*r5(i,j,k) +   &  
                           r10(i,j,k)*r6(i,j,k)           +   &  
                           r11(i,j,k)*r7(i,j,k))
            END DO
         END DO
      END DO
!
! Computes d_y(zre) e_parallel.d_y(grad(zim))
!
      CALL derivk3(a,c7,2)
      CALL derivk3(c1,c4,2)
      CALL derivk3(c2,c5,2)
      CALL derivk3(c3,c6,2)
      CALL fftp3d_complex_to_real(plancr,c7,r3,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,c4,r5,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,c5,r6,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,c6,r7,MPI_COMM_WORLD)
      r3 = r3*tmp
      r5 = r5*tmp
      r6 = r6*tmp
      r7 = r7*tmp

!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
         DO j = 1,ny
            DO i = 1,nx
               r4(i,j,k) = r4(i,j,k)                      +   &
                           r3(i,j,k)*(r9(i,j,k)*r5(i,j,k) +   &  
                           r10(i,j,k)*r6(i,j,k)           +   &  
                           r11(i,j,k)*r7(i,j,k))
            END DO
         END DO
      END DO
!
! Computes d_z(zre) e_parallel.d_z(grad(zim))
!
      CALL derivk3(a,c7,3)
      CALL derivk3(c1,c4,3)
      CALL derivk3(c2,c5,3)
      CALL derivk3(c3,c6,3)
      CALL fftp3d_complex_to_real(plancr,c7,r3,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,c4,r5,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,c5,r6,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,c6,r7,MPI_COMM_WORLD)
      r3 = r3*tmp
      r5 = r5*tmp
      r6 = r6*tmp
      r7 = r7*tmp

!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
         DO j = 1,ny
            DO i = 1,nx
               r4(i,j,k) = r4(i,j,k)                      +   &
                           r3(i,j,k)*(r9(i,j,k)*r5(i,j,k) +   &  
                           r10(i,j,k)*r6(i,j,k)           +   &  
                           r11(i,j,k)*r7(i,j,k))
            END DO
         END DO
      END DO

!
! Computes -e_parallel.grad(zim).grad(grad(zre))
! 
      CALL derivk3(a,c1,1)
      CALL derivk3(a,c2,2)
      CALL derivk3(a,c3,3)
!
! Computes - d_x(zim) e_parallel.d_x(grad(zre))
!
      CALL derivk3(b,c7,1)
      CALL derivk3(c1,c4,1)
      CALL derivk3(c2,c5,1)
      CALL derivk3(c3,c6,1)
      CALL fftp3d_complex_to_real(plancr,c7,r3,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,c4,r5,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,c5,r6,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,c6,r7,MPI_COMM_WORLD)
      r3 = r3*tmp
      r5 = r5*tmp
      r6 = r6*tmp
      r7 = r7*tmp

!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
         DO j = 1,ny
            DO i = 1,nx
               r4(i,j,k) = r4(i,j,k)                      -   &
                           r3(i,j,k)*(r9(i,j,k)*r5(i,j,k) +   &  
                           r10(i,j,k)*r6(i,j,k)           +   &  
                           r11(i,j,k)*r7(i,j,k))
            END DO
         END DO
      END DO
!
! Computes - d_y(zim) e_parallel.d_y(grad(zre))
!
      CALL derivk3(b,c7,2)
      CALL derivk3(c1,c4,2)
      CALL derivk3(c2,c5,2)
      CALL derivk3(c3,c6,2)
      CALL fftp3d_complex_to_real(plancr,c7,r3,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,c4,r5,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,c5,r6,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,c6,r7,MPI_COMM_WORLD)
      r3 = r3*tmp
      r5 = r5*tmp
      r6 = r6*tmp
      r7 = r7*tmp

!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
         DO j = 1,ny
            DO i = 1,nx
               r4(i,j,k) = r4(i,j,k)                      -   &
                           r3(i,j,k)*(r9(i,j,k)*r5(i,j,k) +   &  
                           r10(i,j,k)*r6(i,j,k)           +   &  
                           r11(i,j,k)*r7(i,j,k))
            END DO
         END DO
      END DO
!
! Computes - d_z(zim) e_parallel.d_z(grad(zre))
!
      CALL derivk3(b,c7,3)
      CALL derivk3(c1,c4,3)
      CALL derivk3(c2,c5,3)
      CALL derivk3(c3,c6,3)
      CALL fftp3d_complex_to_real(plancr,c7,r3,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,c4,r5,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,c5,r6,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,c6,r7,MPI_COMM_WORLD)
      r3 = r3*tmp
      r5 = r5*tmp
      r6 = r6*tmp
      r7 = r7*tmp

!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
         DO j = 1,ny
            DO i = 1,nx     
               r4(i,j,k) = 2*alpha*(r4(i,j,k)             -   &
                           r3(i,j,k)*(r9(i,j,k)*r5(i,j,k) +   &  
                           r10(i,j,k)*r6(i,j,k)           +   &  
                           r11(i,j,k)*r7(i,j,k)))!*r12(i,j,k)
            END DO
         END DO
      END DO
!
! Does v_parallel*e_parallel
!
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
         DO j = 1,ny
            DO i = 1,nx
               rmp = (real(nx,kind=GP)*real(ny,kind=GP)* &
                     real(nz,kind=GP))**2/(r9(i,j,k)**2 + &
                     r10(i,j,k)**2 + r11(i,j,k)**2 + &
                     (real(nx,kind=GP)*real(ny,kind=GP)* &
                     real(nz,kind=GP))**2*regu*omegag/beta)
               r9(i,j,k)  = r4(i,j,k)*r9(i,j,k)*rmp*r12(i,j,k)
               r10(i,j,k) = r4(i,j,k)*r10(i,j,k)*rmp*r12(i,j,k)
               r11(i,j,k) = r4(i,j,k)*r11(i,j,k)*rmp*r12(i,j,k)
            END DO
         END DO
      END DO
      CALL fftp3d_real_to_complex(planrc,r9, c2,MPI_COMM_WORLD)
      CALL fftp3d_real_to_complex(planrc,r10,c3,MPI_COMM_WORLD)
      CALL fftp3d_real_to_complex(planrc,r11,c4,MPI_COMM_WORLD)

!
! Obtains v by computing 2.alpha[zre*grad(zim)-zim*grad(zre)]
! and dividing by rho.
!
      CALL derivk3(a,c1,1)   ! x component
      CALL fftp3d_complex_to_real(plancr,c1,r3,MPI_COMM_WORLD)
      CALL derivk3(b,c1,1)
      CALL fftp3d_complex_to_real(plancr,c1,r4,MPI_COMM_WORLD)
!$omp parallel do if (kend-ksta.ge.nth) private (j,i,rmp)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i,rmp)
         DO j = 1,ny
            DO i = 1,nx
               rmp = 1.0_GP/(r1(i,j,k)**2+r2(i,j,k)**2+       &
                     (real(nx,kind=GP)*real(ny,kind=GP)*      &
                     real(nz,kind=GP))**2*regu*omegag/beta)
               r3(i,j,k) = 2*alpha*(r1(i,j,k)*r4(i,j,k)-      &
                           r2(i,j,k)*r3(i,j,k))*rmp       ! v_x
            END DO
         END DO
      END DO
      CALL fftp3d_real_to_complex(planrc,r3,c5,MPI_COMM_WORLD)

      CALL derivk3(a,c1,2)   ! y component
      CALL fftp3d_complex_to_real(plancr,c1,r3,MPI_COMM_WORLD)
      CALL derivk3(b,c1,2)
      CALL fftp3d_complex_to_real(plancr,c1,r4,MPI_COMM_WORLD)
!$omp parallel do if (kend-ksta.ge.nth) private (j,i,rmp)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i,rmp)
         DO j = 1,ny
            DO i = 1,nx
               rmp = 1.0_GP/(r1(i,j,k)**2+r2(i,j,k)**2+       &
                     (real(nx,kind=GP)*real(ny,kind=GP)*      &
                     real(nz,kind=GP))**2*regu*omegag/beta)
               r3(i,j,k) = 2*alpha*(r1(i,j,k)*r4(i,j,k)-      &
                           r2(i,j,k)*r3(i,j,k))*rmp       ! v_y
            END DO
         END DO
      END DO
      CALL fftp3d_real_to_complex(planrc,r3,c6,MPI_COMM_WORLD)

      CALL derivk3(a,c1,3)   ! z component
      CALL fftp3d_complex_to_real(plancr,c1,r3,MPI_COMM_WORLD)
      CALL derivk3(b,c1,3)
      CALL fftp3d_complex_to_real(plancr,c1,r4,MPI_COMM_WORLD)
!$omp parallel do if (kend-ksta.ge.nth) private (j,i,rmp)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i,rmp)
         DO j = 1,ny
            DO i = 1,nx
               rmp = 1.0_GP/(r1(i,j,k)**2+r2(i,j,k)**2+       &
                     (real(nx,kind=GP)*real(ny,kind=GP)*      &
                     real(nz,kind=GP))**2*regu*omegag/beta)
               r3(i,j,k) = 2*alpha*(r1(i,j,k)*r4(i,j,k)-      &
                           r2(i,j,k)*r3(i,j,k))*rmp       ! v_z
            END DO
         END DO
      END DO
      CALL fftp3d_real_to_complex(planrc,r3,c1,MPI_COMM_WORLD)
!
! Computes the helicity
!
      CALL helicity(c5,c6,c1,Htot2)
!
! Computes the regularized helicity
!
      CALL rotor3(c6,c1,c7,1)
      CALL rotor3(c5,c1,c8,2)
      CALL rotor3(c5,c6,c9,3)
      CALL cross(c2,c3,c4,c7,c8,c9,Htot1,1)
!
! Writes the result to a file
!
      IF (myrank.eq.0) THEN
         OPEN(1,file='helicity.txt',position='append')
         WRITE(1,20) (t-1)*dt,Htot1,Htot2
20    FORMAT( E13.6,E22.14,E22.14 )
         CLOSE(1)
      ENDIF

      RETURN
      END SUBROUTINE gpehelicity

!***********************************************************************
      SUBROUTINE gpehelspec(a,b,nmb)
!-----------------------------------------------------------------------
!
! Computes two spectra of helicity, one based on a regularized parallel
! velocity and one on v (usual definition of the helicity for a
! classical fluid). The output is written to a file by the first node.
!
! Output files contain:
! 'hspectrum.XXX.txt': k, H_regularized(k), H_classical(k)
!
! Parameters
!     a    : real part of the wavefunction in Fourier space
!     b    : imaginary part of the wavefunction in Fourier space
!     nmb  : the extension used when writting the file
!
      USE fprecision
      USE commtypes
      USE kes
      USE fft
      USE ali
      USE grid
      USE hbar
      USE boxsize
      USE mpivars
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a,b
      COMPLEX(KIND=GP), DIMENSION(nz,ny,ista:iend) :: c1,c2,c3
      COMPLEX(KIND=GP), DIMENSION(nz,ny,ista:iend) :: c4,c5,c6
      COMPLEX(KIND=GP), DIMENSION(nz,ny,ista:iend) :: c7,c8,c9
      REAL(KIND=GP), DIMENSION(nx,ny,ksta:kend)    :: r1,r2,r3,r4
      REAL(KIND=GP), DIMENSION(nx,ny,ksta:kend)    :: r5,r6,r7,r8
      REAL(KIND=GP), DIMENSION(nx,ny,ksta:kend)    :: r9,r10,r11,r12
      DOUBLE PRECISION, DIMENSION(nmax/2+1) :: Htot3,Htot4
      CHARACTER(len=*), INTENT(IN)       :: nmb
      REAL(KIND=GP)       :: rmp,rmq
      REAL(KIND=GP)       :: tmp
      INTEGER             :: i,j,k
      INTEGER             :: kmn

!
! Transforms the wavefunction to real space
!
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
            DO k = 1,nz
               c1(k,j,i) = a(k,j,i)
               c2(k,j,i) = b(k,j,i)
            END DO
         END DO
      END DO
      CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,c2,r2,MPI_COMM_WORLD)
!
! Obtains e_parallel by calculating the cross product of
! grad(zre) and grad(zim). Then it normalizes.
!
      CALL derivk3(a,c1,1)
      CALL fftp3d_complex_to_real(plancr,c1,r3,MPI_COMM_WORLD)
      CALL derivk3(b,c1,1)
      CALL fftp3d_complex_to_real(plancr,c1,r4,MPI_COMM_WORLD)
      CALL derivk3(a,c1,2)
      CALL fftp3d_complex_to_real(plancr,c1,r5,MPI_COMM_WORLD)
      CALL derivk3(b,c1,2)
      CALL fftp3d_complex_to_real(plancr,c1,r6,MPI_COMM_WORLD)
      CALL derivk3(a,c1,3)  
      CALL fftp3d_complex_to_real(plancr,c1,r7,MPI_COMM_WORLD)
      CALL derivk3(b,c1,3)
      CALL fftp3d_complex_to_real(plancr,c1,r8,MPI_COMM_WORLD)
      tmp = 1.0_GP/ &
            (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
         DO j = 1,ny
            DO i = 1,nx
               r9(i,j,k)  = (r5(i,j,k)*r8(i,j,k)-r7(i,j,k)*          &
                            r6(i,j,k))*tmp
               r10(i,j,k) = (r7(i,j,k)*r4(i,j,k)-r3(i,j,k)*          &
                            r8(i,j,k))*tmp
               r11(i,j,k) = (r3(i,j,k)*r6(i,j,k)-r5(i,j,k)*          &
                            r4(i,j,k))*tmp
               r12(i,j,k) = (real(nx,kind=GP)*real(ny,kind=GP)*      &
                            real(nz,kind=GP))**2/(r3(i,j,k)**2 +     &
                            r4(i,j,k)**2 +                           &
                            r5(i,j,k)**2 + r6(i,j,k)**2 +            &
                            r7(i,j,k)**2 + r8(i,j,k)**2 +            &
                            (real(nx,kind=GP)*real(ny,kind=GP)*      &
                            real(nz,kind=GP))**2*regu*omegag/beta)
            END DO
         END DO
      END DO

!
! Calculates v_parallel by doing regularizing and the projecting
! The magnitud calculated is:
! 2*alpha e_parallel.((grad(zre).grad)(grad(zim)) -
! (grad(zim).grad)(grad(zre)))/(grad(zre)**2+grad(zim)**2)
! 
      CALL derivk3(b,c1,1)
      CALL derivk3(b,c2,2)
      CALL derivk3(b,c3,3)
!
! Computes d_x(zre) e_parallel.d_x(grad(zim))
!
      CALL derivk3(a,c7,1)
      CALL derivk3(c1,c4,1)
      CALL derivk3(c2,c5,1)
      CALL derivk3(c3,c6,1)
      CALL fftp3d_complex_to_real(plancr,c7,r3,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,c4,r5,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,c5,r6,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,c6,r7,MPI_COMM_WORLD)
      r3 = r3*tmp
      r5 = r5*tmp
      r6 = r6*tmp
      r7 = r7*tmp

!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
         DO j = 1,ny
            DO i = 1,nx
               r4(i,j,k) = r3(i,j,k)*(r9(i,j,k)*r5(i,j,k) +   &  
                           r10(i,j,k)*r6(i,j,k)           +   &  
                           r11(i,j,k)*r7(i,j,k))
            END DO
         END DO
      END DO
!
! Computes d_y(zre) e_parallel.d_y(grad(zim))
!
      CALL derivk3(a,c7,2)
      CALL derivk3(c1,c4,2)
      CALL derivk3(c2,c5,2)
      CALL derivk3(c3,c6,2)
      CALL fftp3d_complex_to_real(plancr,c7,r3,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,c4,r5,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,c5,r6,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,c6,r7,MPI_COMM_WORLD)
      r3 = r3*tmp
      r5 = r5*tmp
      r6 = r6*tmp
      r7 = r7*tmp

!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
         DO j = 1,ny
            DO i = 1,nx
               r4(i,j,k) = r4(i,j,k)                      +   &
                           r3(i,j,k)*(r9(i,j,k)*r5(i,j,k) +   &  
                           r10(i,j,k)*r6(i,j,k)           +   &  
                           r11(i,j,k)*r7(i,j,k))
            END DO
         END DO
      END DO
!
! Computes d_z(zre) e_parallel.d_z(grad(zim))
!
      CALL derivk3(a,c7,3)
      CALL derivk3(c1,c4,3)
      CALL derivk3(c2,c5,3)
      CALL derivk3(c3,c6,3)
      CALL fftp3d_complex_to_real(plancr,c7,r3,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,c4,r5,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,c5,r6,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,c6,r7,MPI_COMM_WORLD)
      r3 = r3*tmp
      r5 = r5*tmp
      r6 = r6*tmp
      r7 = r7*tmp

!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
         DO j = 1,ny
            DO i = 1,nx
               r4(i,j,k) = r4(i,j,k)                      +   &
                           r3(i,j,k)*(r9(i,j,k)*r5(i,j,k) +   &  
                           r10(i,j,k)*r6(i,j,k)           +   &  
                           r11(i,j,k)*r7(i,j,k))
            END DO
         END DO
      END DO

!
! Computes -e_parallel.grad(zim).grad(grad(zre))
! 
      CALL derivk3(a,c1,1)
      CALL derivk3(a,c2,2)
      CALL derivk3(a,c3,3)
!
! Computes - d_x(zim) e_parallel.d_x(grad(zre))
!
      CALL derivk3(b,c7,1)
      CALL derivk3(c1,c4,1)
      CALL derivk3(c2,c5,1)
      CALL derivk3(c3,c6,1)
      CALL fftp3d_complex_to_real(plancr,c7,r3,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,c4,r5,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,c5,r6,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,c6,r7,MPI_COMM_WORLD)
      r3 = r3*tmp
      r5 = r5*tmp
      r6 = r6*tmp
      r7 = r7*tmp

!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
         DO j = 1,ny
            DO i = 1,nx
               r4(i,j,k) = r4(i,j,k)                      -   &
                           r3(i,j,k)*(r9(i,j,k)*r5(i,j,k) +   &  
                           r10(i,j,k)*r6(i,j,k)           +   &  
                           r11(i,j,k)*r7(i,j,k))
            END DO
         END DO
      END DO
!
! Computes - d_y(zim) e_parallel.d_y(grad(zre))
!
      CALL derivk3(b,c7,2)
      CALL derivk3(c1,c4,2)
      CALL derivk3(c2,c5,2)
      CALL derivk3(c3,c6,2)
      CALL fftp3d_complex_to_real(plancr,c7,r3,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,c4,r5,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,c5,r6,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,c6,r7,MPI_COMM_WORLD)
      r3 = r3*tmp
      r5 = r5*tmp
      r6 = r6*tmp
      r7 = r7*tmp

!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
         DO j = 1,ny
            DO i = 1,nx
               r4(i,j,k) = r4(i,j,k)                      -   &
                           r3(i,j,k)*(r9(i,j,k)*r5(i,j,k) +   &  
                           r10(i,j,k)*r6(i,j,k)           +   &  
                           r11(i,j,k)*r7(i,j,k))
            END DO
         END DO
      END DO
!
! Computes - d_z(zim) e_parallel.d_z(grad(zre))
!
      CALL derivk3(b,c7,3)
      CALL derivk3(c1,c4,3)
      CALL derivk3(c2,c5,3)
      CALL derivk3(c3,c6,3)
      CALL fftp3d_complex_to_real(plancr,c7,r3,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,c4,r5,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,c5,r6,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,c6,r7,MPI_COMM_WORLD)
      r3 = r3*tmp
      r5 = r5*tmp
      r6 = r6*tmp
      r7 = r7*tmp

!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
         DO j = 1,ny
            DO i = 1,nx     
               r4(i,j,k) = 2*alpha*(r4(i,j,k)             -   &
                           r3(i,j,k)*(r9(i,j,k)*r5(i,j,k) +   &  
                           r10(i,j,k)*r6(i,j,k)           +   &  
                           r11(i,j,k)*r7(i,j,k)))!*r12(i,j,k)
            END DO
         END DO
      END DO
!
! Does v_parallel*e_parallel
!
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
         DO j = 1,ny
            DO i = 1,nx
               rmp = (real(nx,kind=GP)*real(ny,kind=GP)*        &
                     real(nz,kind=GP))**2/(r9(i,j,k)**2 +       &
                     r10(i,j,k)**2 + r11(i,j,k)**2 +            &
                     (real(nx,kind=GP)*real(ny,kind=GP)*        &
                     real(nz,kind=GP))**2*regu*omegag/beta)
               r9(i,j,k)  = r4(i,j,k)*r9(i,j,k)*rmp*r12(i,j,k)
               r10(i,j,k) = r4(i,j,k)*r10(i,j,k)*rmp*r12(i,j,k)
               r11(i,j,k) = r4(i,j,k)*r11(i,j,k)*rmp*r12(i,j,k)
            END DO
         END DO
      END DO
      CALL fftp3d_real_to_complex(planrc,r9, c2,MPI_COMM_WORLD)
      CALL fftp3d_real_to_complex(planrc,r10,c3,MPI_COMM_WORLD)
      CALL fftp3d_real_to_complex(planrc,r11,c4,MPI_COMM_WORLD)

!
! Obtains v by computing 2.alpha[zre*grad(zim)-zim*grad(zre)]
! and dividing by rho.
!
      CALL derivk3(a,c1,1)   ! x component
      CALL fftp3d_complex_to_real(plancr,c1,r3,MPI_COMM_WORLD)
      CALL derivk3(b,c1,1)
      CALL fftp3d_complex_to_real(plancr,c1,r4,MPI_COMM_WORLD)
!$omp parallel do if (kend-ksta.ge.nth) private (j,i,rmp)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i,rmp)
         DO j = 1,ny
            DO i = 1,nx
               rmp = 1.0_GP/(r1(i,j,k)**2+r2(i,j,k)**2+       &
                     (real(nx,kind=GP)*real(ny,kind=GP)*      &
                     real(nz,kind=GP))**2*regu*omegag/beta)
               r3(i,j,k) = 2*alpha*(r1(i,j,k)*r4(i,j,k)-      &
                           r2(i,j,k)*r3(i,j,k))*rmp       ! v_x
            END DO
         END DO
      END DO
      CALL fftp3d_real_to_complex(planrc,r3,c5,MPI_COMM_WORLD)

      CALL derivk3(a,c1,2)   ! y component
      CALL fftp3d_complex_to_real(plancr,c1,r3,MPI_COMM_WORLD)
      CALL derivk3(b,c1,2)
      CALL fftp3d_complex_to_real(plancr,c1,r4,MPI_COMM_WORLD)
!$omp parallel do if (kend-ksta.ge.nth) private (j,i,rmp)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i,rmp)
         DO j = 1,ny
            DO i = 1,nx
               rmp = 1.0_GP/(r1(i,j,k)**2+r2(i,j,k)**2+       &
                     (real(nx,kind=GP)*real(ny,kind=GP)*      &
                     real(nz,kind=GP))**2*regu*omegag/beta)
               r3(i,j,k) = 2*alpha*(r1(i,j,k)*r4(i,j,k)-      &
                           r2(i,j,k)*r3(i,j,k))*rmp       ! v_y
            END DO
         END DO
      END DO
      CALL fftp3d_real_to_complex(planrc,r3,c6,MPI_COMM_WORLD)

      CALL derivk3(a,c1,3)   ! z component
      CALL fftp3d_complex_to_real(plancr,c1,r3,MPI_COMM_WORLD)
      CALL derivk3(b,c1,3)
      CALL fftp3d_complex_to_real(plancr,c1,r4,MPI_COMM_WORLD)
!$omp parallel do if (kend-ksta.ge.nth) private (j,i,rmp)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i,rmp)
         DO j = 1,ny
            DO i = 1,nx
               rmp = 1.0_GP/(r1(i,j,k)**2+r2(i,j,k)**2+       &
                     (real(nx,kind=GP)*real(ny,kind=GP)*      &
                     real(nz,kind=GP))**2*regu*omegag/beta)
               r3(i,j,k) = 2*alpha*(r1(i,j,k)*r4(i,j,k)-      &
                           r2(i,j,k)*r3(i,j,k))*rmp       ! v_z
            END DO
         END DO
      END DO
      CALL fftp3d_real_to_complex(planrc,r3,c1,MPI_COMM_WORLD)
!
! Computes the helicity spectra
!
      CALL gpespectrumc(c5,c6,c1,Htot4)
!
! Computes the regularized helicity spectra
!
      CALL crosspecc(c2,c3,c4,c5,c6,c1,Htot3,1.0_GP)
      IF (myrank.eq.0) THEN
         OPEN(1,file='hspectrum.' // nmb // '.txt')
         DO i = 1,nmax/2+1
            WRITE(1,30) Dkk*(i-1),Htot3(i)/Dkk,Htot4(i)/Dkk
         END DO
30       FORMAT( E13.6,E23.15,E23.15 )
         CLOSE(1)
      ENDIF

      RETURN
      END SUBROUTINE

!***********************************************************************
      SUBROUTINE gpespectrumc(a,b,c,Hktot)
!-----------------------------------------------------------------------
!
! Computes the helicity power spectra, returning it. The spectra start
! at k=0 to preserve information of the energy in the condensate, and
! are not dealiased (k = 0,1,...,N/2). 
!
! Parameters
!     a    : input matrix in the x-direction
!     b    : input matrix in the y-direction
!     c    : input matrix in the z-direction
!     Hktot: output helicity spectrum
!
      USE fprecision
      USE commtypes
      USE kes
      USE grid
      USE mpivars
      USE boxsize
!$    USE threads
      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(nmax/2+1) :: Ek
      DOUBLE PRECISION, INTENT(OUT), DIMENSION(nmax/2+1) :: Hktot
      DOUBLE PRECISION    :: tmq
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a,b,c
      COMPLEX(KIND=GP), DIMENSION(nz,ny,ista:iend)             :: c1,c2,c3
      REAL(KIND=GP)       :: tmp
      INTEGER             :: i,j,k
      INTEGER             :: kmn

!
! Computes the curl of the field
!
      CALL rotor3(b,c,c1,1)
      CALL rotor3(a,c,c2,2)
      CALL rotor3(a,b,c3,3)
!
! Computes the kinetic energy spectrum
!
      tmp = 1.0_GP/ &
            (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2
!
! Computes the helicity spectrum
!
      DO i = 1,nmax/2+1
         Ek(i) = 0.0D0
         Hktot(i) = 0.0D0
      END DO
      IF (ista.eq.1) THEN
!$omp parallel do private (k,kmn,tmq)
         DO j = 1,ny
            DO k = 1,nz
               kmn = int(sqrt(kk2(k,j,1))/Dkk+1.501)
               IF ((kmn.gt.0).and.(kmn.le.nmax/2+1)) THEN
                  tmq = (real(a(k,j,1)*conjg(c1(k,j,1)))+          &
                         real(b(k,j,1)*conjg(c2(k,j,1)))+          &
                         real(c(k,j,1)*conjg(c3(k,j,1))))*tmp
!$omp atomic
                  Ek(kmn) = Ek(kmn)+tmq
               ENDIF
            END DO
         END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k,kmn,tmq)
         DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k,kmn,tmq)
            DO j = 1,ny
               DO k = 1,nz
                  kmn = int(sqrt(kk2(k,j,i))/Dkk+1.501)
                  IF ((kmn.gt.0).and.(kmn.le.nmax/2+1)) THEN
                     tmq = 2*(real(a(k,j,i)*conjg(c1(k,j,i)))+     &
                              real(b(k,j,i)*conjg(c2(k,j,i)))+     &
                              real(c(k,j,i)*conjg(c3(k,j,i))))*tmp
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
            DO j = 1,ny
               DO k = 1,nz
                  kmn = int(sqrt(kk2(k,j,i))/Dkk+1.501)
                  IF ((kmn.gt.0).and.(kmn.le.nmax/2+1)) THEN
                     tmq = 2*(real(a(k,j,i)*conjg(c1(k,j,i)))+     &
                              real(b(k,j,i)*conjg(c2(k,j,i)))+     &
                              real(c(k,j,i)*conjg(c3(k,j,i))))*tmp
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
      CALL MPI_ALLREDUCE(Ek,Hktot,nmax/2+1,MPI_DOUBLE_PRECISION,   &
                      MPI_SUM,MPI_COMM_WORLD,ierr)

      RETURN
      END SUBROUTINE gpespectrumc

!***********************************************************************
      SUBROUTINE gpemomtspec(a,b,nmb)
!-----------------------------------------------------------------------
!
! Computes the spectrum of momentum. These quantities must be computed
! in real space first, and then transformed to Fourier space to compute
! the spectrum. The spectra start at k=0 to preserve information of the
! energy in the condensate, and are not dealiased (k = 0,1,...,N/2). The
! output is written to files by the first node.
!
! Output files contain:
! 'momtspectrum.XXX.txt': k,P(k) [spectrum of incompressible momentum]
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
      USE var
      USE grid
      USE hbar
      USE mpivars
      USE filefmt
      USE boxsize
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a,b
      COMPLEX(KIND=GP), DIMENSION(nz,ny,ista:iend) :: c1,c2,c3,c4
      DOUBLE PRECISION, DIMENSION(nmax/2+1)        :: Ek,Ektot
      REAL(KIND=GP), DIMENSION(nx,ny,ksta:kend)    :: r1,r2,r3,r4,r5
      REAL(KIND=GP)    :: rmq
      INTEGER          :: i,j,k
      INTEGER          :: kmn
      CHARACTER(len=*), INTENT(IN) :: nmb

!
! Transforms the wavefunction to real space
!
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
            DO k = 1,nz
               c1(k,j,i) = a(k,j,i)
               c2(k,j,i) = b(k,j,i)
            END DO
         END DO
      END DO
      CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,c2,r2,MPI_COMM_WORLD)
!
! Do derivatives and multiply in real space
!
      rmq = 2*alpha/(real(nx,kind=GP)*real(ny,kind=GP)* &
            real(nz,kind=GP))**2
      CALL derivk3(a,c1,1)   ! x component
      CALL fftp3d_complex_to_real(plancr,c1,r3,MPI_COMM_WORLD)
      CALL derivk3(b,c1,1)
      CALL fftp3d_complex_to_real(plancr,c1,r5,MPI_COMM_WORLD)
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
         DO j = 1,ny
            DO i = 1,nx
               r4(i,j,k) = (r2(i,j,k)*r3(i,j,k) - &
                            r1(i,j,k)*r5(i,j,k))*rmq
            END DO
         END DO
      END DO
      CALL fftp3d_real_to_complex(planrc,r4,c2,MPI_COMM_WORLD)

      CALL derivk3(a,c1,2)   ! y component
      CALL fftp3d_complex_to_real(plancr,c1,r3,MPI_COMM_WORLD)
      CALL derivk3(b,c1,2)
      CALL fftp3d_complex_to_real(plancr,c1,r5,MPI_COMM_WORLD)
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
         DO j = 1,ny
            DO i = 1,nx
               r4(i,j,k) = (r2(i,j,k)*r3(i,j,k) - &
                            r1(i,j,k)*r5(i,j,k))*rmq
            END DO
         END DO
      END DO
      CALL fftp3d_real_to_complex(planrc,r4,c3,MPI_COMM_WORLD)

      CALL derivk3(a,c1,3)   ! z component
      CALL fftp3d_complex_to_real(plancr,c1,r3,MPI_COMM_WORLD)
      CALL derivk3(b,c1,3)
      CALL fftp3d_complex_to_real(plancr,c1,r5,MPI_COMM_WORLD)
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
         DO j = 1,ny
            DO i = 1,nx
               r4(i,j,k) = (r2(i,j,k)*r3(i,j,k) - &
                            r1(i,j,k)*r5(i,j,k))*rmq
            END DO
         END DO
      END DO
      CALL fftp3d_real_to_complex(planrc,r4,c4,MPI_COMM_WORLD)
!
! Computes the incompressible momentum spectrum
!
      CALL gauge3(c2,c3,c4,c1,1)      ! x component 
      CALL spectrscc(c1,Ektot,1.0_GP) ! incompressible

      CALL gauge3(c2,c3,c4,c1,2)      ! y component
      CALL spectrscc(c1,Ek,1.0_GP)    ! incompressible
      IF (myrank.eq.0) THEN
         Ektot = Ektot+Ek
      ENDIF

      CALL gauge3(c2,c3,c4,c1,3)      ! z component
      CALL spectrscc(c1,Ek,1.0_GP)    ! incompressible
      IF (myrank.eq.0) THEN
         Ektot = (Ektot+Ek)
         OPEN(1,file='momtspectrum.' // nmb // '.txt')
         DO i=1,nmax/2+1
            WRITE(1,FMT='(E13.6,E23.15)') Dkk*(i-1),Ektot(i)/Dkk
         END DO
         CLOSE(1)
      ENDIF

      RETURN
      END SUBROUTINE gpemomtspec

!**********************************************************************
      SUBROUTINE combine(r1,r2,dir,nmb,plan)
!----------------------------------------------------------------------
!
! Combine wavefunctions from ARGL and SGLE for finite temperature runs.
! Output is in real space.
!
! Parameters
!     r1 : real part of the wavefunction in real space [OUT]
!     r2 : imaginary part of the wavefunction in real space [OUT]
!     dir   : directory from which the files are read [IN]
!     nmb   : extension with the time label [IN]. Should be zero.
!     plan  : I/O plan [IN]
!

      USE fprecision
      USE commtypes
      USE mpivars
      USE iovar
      USE iompi
      USE grid
      USE kes
      USE ali
      USE fft
!$    USE threads
      IMPLICIT NONE

      REAL(KIND=GP), INTENT(OUT), DIMENSION(nx,ny,ksta:kend)  :: r1,r2
      CHARACTER(len=100), INTENT(IN) :: dir
      CHARACTER(len=*), INTENT(IN)   :: nmb
      TYPE(IOPLAN),INTENT  (IN)      :: plan

      REAL(KIND=GP), DIMENSION(nx,ny,ksta:kend)    :: r3,r4,r5,r6
      INTEGER :: i,j,k


      ! ARGL field
      CALL io_read(1,dir,'argl_re',nmb,plan,r3)
      CALL io_read(1,dir,'argl_im',nmb,plan,r4)

      ! SGLE field
      CALL io_read(1,dir,'sgle_re',nmb,plan,r5)
      CALL io_read(1,dir,'sgle_im',nmb,plan,r6)

      ! Combine and transform
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
         DO j = 1,ny
            DO i = 1,nx
               R1(i,j,k) = R3(i,j,k)*R5(i,j,k) - R4(i,j,k)*R6(i,j,k)
               R2(i,j,k) = R3(i,j,k)*R6(i,j,k) + R4(i,j,k)*R5(i,j,k)
            END DO
         END DO
      END DO

      RETURN
      END SUBROUTINE combine
