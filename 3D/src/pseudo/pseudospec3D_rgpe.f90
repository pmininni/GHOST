!=================================================================
! PSEUDOSPECTRAL subroutines
!
! Subroutines to compute energies and spectra in the GPE and
! ARLG equations with rotation and/or with a trapping potential.
! You should use the FFTPLANS and MPIVARS modules (see the file
! 'fftp_mod.f90') in each program that calls any of the
! subroutines in this file. 
!
! NOTATION: index 'i' is 'x' 
!           index 'j' is 'y'
!           index 'k' is 'z'
!
! 2020 Julián Amette Estrada
!      Department of Physics,
!      Facultad de Ciencias Exactas y Naturales.
!      Universidad de Buenos Aires.
!=================================================================

!*****************************************************************
      SUBROUTINE trapenergy(a,b,v,alpha,t,dt)
!-----------------------------------------------------------------
!
! Computes the potential energy associated to the trapping
! potential, Etrap = < V(x) |z|^2 >, and saves the output to
! a file.
!
! Output file contains:
! 'trenergy.txt': time, harmonic potential energy <V(x)|z|^2>
!
! Parameters
!     a    : real part of the wavefunction in Fourier space
!     b    : imaginary part of the wavefunction in Fourier space
!     v    : potential of the harmonic trap in Fourier space
!     alpha: amplitude of the alpha coefficient in GPE
!     t    : number of time steps made
!     dt   : time step
!
      USE fprecision
      USE commtypes
      USE mpivars
      USE grid
      USE kes
      USE ali
      USE fft
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a,b
      DOUBLE PRECISION          :: tmp,ene,enet
      REAL(KIND=GP), INTENT(IN), DIMENSION(nx,ny,ksta:kend)    :: v
      REAL(KIND=GP), DIMENSION(nx,ny,ksta:kend)                :: r
      REAL(KIND=GP), INTENT(IN) :: dt,alpha
      INTEGER, INTENT(IN)       :: t
      INTEGER                   :: i,j,k

!
! Computes the square of the wavefunction
!
      CALL squareabs(a,b,r,0) ! Non-dealiased and normalized
! 
! Computes product with the potential in real space
!
      ene = 0.0D0
      tmp = 1.0D0/ &
            (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
         DO j = 1,ny
            DO i = 1,nx
               ene  = ene + v(i,j,k)*r(i,j,k)*tmp
            END DO
         END DO
      END DO
      ene = 2.0D0*ene*alpha
      CALL MPI_REDUCE(ene,enet,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
!
! Writes the output to a file
!
      IF (myrank.eq.0) THEN
         OPEN(1,file='trenergy.txt',position='append')
         WRITE(1,FMT='(E13.6,E22.14)') (t-1)*dt,enet
         CLOSE(1)
      ENDIF
      END SUBROUTINE trapenergy

!*****************************************************************
      SUBROUTINE rotenergy(a,b,vlinx,vliny,omegaz,alpha,t,dt)
!-----------------------------------------------------------------
!
! Computes the energy associated to the rotation of the
! condensante (including the effective repulsive centrifugal
! potential), in the rotating frame of reference,
! Erot = < z* Omega . (r x P) z >, where P is the momentum
! operator. This subroutine writes the output to a file.
!
! Output file contains:
! 'rotenergy.txt': time, rot. energy < z* Omega.(r x P) z >
!
! Parameters
!     a     : real part of the wavefunction in Fourier space
!     b     : imaginary part of the wavefunction in Fourier space
!     vlinx : linear x-coordinate (with the fix for the b.c.)
!     vliny : linear y-coordinate (with the fix for the b.c.)
!     omegaz: condensate rotation rate
!     alpha : amplitude of the alpha coefficient in GPE
!     t     : number of time steps made
!     dt    : time step
!
      USE fprecision
      USE commtypes
      USE mpivars
      USE grid
      USE kes
      USE ali
      USE fft
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a,b
      COMPLEX(KIND=GP), DIMENSION(nz,ny,ista:iend)             :: c1
      DOUBLE PRECISION                          :: tmp,tmq,tmr,ene,enet
      REAL(KIND=GP), INTENT(IN), DIMENSION(nx,ny,ksta:kend) :: vlinx,vliny
      REAL(KIND=GP), DIMENSION(nx,ny,ksta:kend) :: r1,r2
      REAL(KIND=GP), DIMENSION(nx,ny,ksta:kend) :: dx1,dx2,dy1,dy2
      REAL(KIND=GP), INTENT(IN)                 :: omegaz,dt,alpha
      INTEGER, INTENT(IN) :: t
      INTEGER             :: i,j,k

!
! Computes < z* Omega.(r x P) z >
!
      c1 = a
      CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      c1 = b
      CALL fftp3d_complex_to_real(plancr,c1,r2,MPI_COMM_WORLD)
      CALL derivk3(a,c1,1) ! dz /dx
      CALL fftp3d_complex_to_real(plancr,c1,dx1,MPI_COMM_WORLD)
      CALL derivk3(b,c1,1) ! dz*/dx
      CALL fftp3d_complex_to_real(plancr,c1,dx2,MPI_COMM_WORLD)
      CALL derivk3(a,c1,2) ! dz /dy
      CALL fftp3d_complex_to_real(plancr,c1,dy1,MPI_COMM_WORLD)
      CALL derivk3(b,c1,2) ! dz*/dy
      CALL fftp3d_complex_to_real(plancr,c1,dy2,MPI_COMM_WORLD)
      ene = 0.0D0
      tmp = 2.0D0*omegaz*alpha/ &
            (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**4
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
         DO j = 1,ny
            DO i = 1,nx
               tmq = r2(i,j,k)*(vliny(i,j,k)*dx1(i,j,k)-vlinx(i,j,k) &
                     *dy1(i,j,k))
               tmr = r1(i,j,k)*(vlinx(i,j,k)*dy2(i,j,k)-vliny(i,j,k) &
                     *dx2(i,j,k))
               ene  = ene + tmp*(tmq+tmr) 
            END DO
         END DO
      END DO
      CALL MPI_REDUCE(ene,enet,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
!
! Writes the result to a file
!
      IF (myrank.eq.0) THEN
         OPEN(1,file='rotenergy.txt',position='append')
         WRITE(1,FMT='(E13.6,E22.14)') (t-1)*dt,enet
         CLOSE(1)
      ENDIF

      END SUBROUTINE rotenergy

!*****************************************************************
      SUBROUTINE gperealspecperp(a,b,nmb)
!-----------------------------------------------------------------
!
! Computes the spectrum of kinetic, quantum, and potential (or 
! internal) energy as a function of the wavenumber perpendicular
! to z. The k-shells are cylindrical surfaces with
! kperp = Dkk*(0,...,max{nx*Dkx/Dkk,nyDky/Dkk}/2). The spectra
! start at k = 0 to preserve information of the energy in the
! condensate, and are not dealiased. The output is written to
! files by the first node.
!
! Output files contain:
! 'intspecperp.XXX.txt'  : k, Eint(k) [Eint = 2.alpha.beta.(|z|^2-rho0)^2]
! 'qspectrumperp.XXX.txt': k, Equa(k)
!   [Equa ~ (zturnre*grad(zre)+zturnim*grad(zim))^2]
! 'kincspecperp.XXX.txt' : k, Einc(k)
! 'kcomspecperp.XXX.txt' : k, Ecom(k)
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
      DOUBLE PRECISION, DIMENSION(nmax/2+1)    :: Eint,Equa,Einc,Ecom
      INTEGER                      :: i
      CHARACTER(len=*), INTENT(IN) :: nmb

!
! Computes all the energy  spectra
!
      CALL gperealspecperpc(a,b,Eint,Equa,Einc,Ecom)
!
! Exports the energy spectrum to a file
!
      IF (myrank.eq.0) THEN
         OPEN(1,file='intspecperp.' // nmb // '.txt')
         DO i=1,nmax/2+1
            WRITE(1,FMT='(E13.6,E23.15)') Dkk*(i-1),Eint(i)/Dkk
         END DO
         CLOSE(1)
         OPEN(1,file='qspecperp.' // nmb // '.txt')
         DO i=1,nmax/2+1                                       
            WRITE(1,FMT='(E13.6,E23.15)') Dkk*(i-1),Equa(i)/Dkk   
         END DO  
         CLOSE(1)
         OPEN(1,file='kincspecperp.' // nmb // '.txt')
         DO i=1,nmax/2+1
            WRITE(1,FMT='(E13.6,E23.15)') Dkk*(i-1),Einc(i)/Dkk
         END DO
         CLOSE(1)
         OPEN(1,file='kcomspecperp.' // nmb // '.txt')
         DO i=1,nmax/2+1
            WRITE(1,FMT='(E13.6,E23.15)') Dkk*(i-1),Ecom(i)/Dkk
         END DO
         CLOSE(1)
      ENDIF

      RETURN
      END SUBROUTINE gperealspecperp

!*****************************************************************
      SUBROUTINE gperealspecperpc(a,b,Eint,Equa,Einc,Ecom)
!-----------------------------------------------------------------
!
! Computes the reduced perpendicular spectrum of kinetic, quantum,
! and potential (or internal) energy, as a function of the
! perpendicular wavenumber, returning them.
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
      DOUBLE PRECISION, INTENT(OUT), DIMENSION(nmax/2+1) ::Eint,Equa,Einc,Ecom
      DOUBLE PRECISION, DIMENSION(nmax/2+1)        ::Etmp
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
      CALL specscpec(c1,Eint,Etmp)
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
      CALL specscpec(c1,Ek1,Etmp)

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
      CALL specscpec(c1,Ek2,Etmp)
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
      CALL specscpec(c1,Ek2,Etmp)
      rmq = 2*alpha**2/ &
            (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2
      Equa = (Ek1+Ek2)*rmq
!
! Computes the compressible and incompressible kinetic energy spectra
!
      CALL gauge3(c2,c3,c4,c1,1)     ! x component 
      CALL specscpec(c1,Einc,Etmp)   ! incompressible
      c1 = c2-c1
      CALL specscpec(c1,Ecom,Etmp)   ! compressible

      CALL gauge3(c2,c3,c4,c1,2)     ! y component
      CALL specscpec(c1,Ek1,Etmp)    ! incompressible
      c1 = c3-c1
      CALL specscpec(c1,Ek2,Etmp)    ! compressible
      Einc = Einc+Ek1
      Ecom = Ecom+Ek2

      CALL gauge3(c2,c3,c4,c1,3)     ! z component
      CALL specscpec(c1,Ek1,Etmp)    ! incompressible
      c1 = c4-c1
      CALL specscpec(c1,Ek2,Etmp)    ! compressible
      rmq = 2*alpha**2/ &
            (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2
      Einc = (Einc+Ek1)*rmq
      Ecom = (Ecom+Ek2)*rmq

      RETURN
      END SUBROUTINE gperealspecperpc
