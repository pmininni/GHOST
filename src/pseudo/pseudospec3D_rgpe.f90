!=================================================================
! PSEUDOSPECTRAL modules
!
! CONTAINS:
!      MODULE pseudospec_rgpe
!
! Subroutines to compute energies and spectra in the GPE and
! ARGL equations with rotation and/or with a trapping potential.
! You should use the FFTPLANS and MPIVARS modules (see the file
! 'fftp_mod.f90') in each program that calls any of the
! subroutines in this file. All the routines are diagnostics
! (computed in the host, with host-only temporaries); the terms
! of the equations with rotation and the potential are built in
! the solvers with the routines of pseudospec_quantum.
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

MODULE pseudospec_rgpe
   USE fprecision
   USE pseudospec_fluid
   USE pseudospec_magnetic
   USE pseudospec_anisca
   USE pseudospec_quantum, ONLY: squareabs
   USE pseudospec_gpe, ONLY: zturn, regu
   USE class_GWorkspace3D, ONLY: gws
   IMPLICIT NONE
   CONTAINS

!*****************************************************************
      SUBROUTINE trapenergy(a,b,v,alpha,beta,t,dt,path)
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
!     v    : potential of the trap divided by beta, in real space
!            (not normalized)
!     alpha: amplitude of the alpha coefficient in GPE
!     t    : number of time steps made
!     dt   : time step
!     path : path for the output
!
      USE fprecision
      USE commtypes
      USE mpivars
      USE grid
      USE kes
      USE ali
      USE fft
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a,b
      REAL(KIND=GP), INTENT(IN), DIMENSION(nx,ny,ksta:kend)    :: v
      REAL(KIND=GP), INTENT(IN)    :: dt,alpha,beta
      INTEGER, INTENT(IN)          :: t
      CHARACTER(len=*), INTENT(IN) :: path
      REAL(KIND=GP), POINTER, DIMENSION(:,:,:) :: r
      DOUBLE PRECISION          :: tmp,ene,enet
      INTEGER                   :: i,j,k
      LOGICAL                   :: bret

      CALL gws%get_real_htmp(r,bret)
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
!$omp parallel do collapse(2) private (i) reduction(+:ene)
      DO k = ksta,kend
         DO j = 1,ny
            DO i = 1,nx
               ene  = ene + v(i,j,k)*r(i,j,k)*tmp
            END DO
         END DO
      END DO
      ene = 2.0D0*ene*alpha*beta ! v is the potential divided by beta
      CALL MPI_REDUCE(ene,enet,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
!
! Writes the output to a file
!
      IF (myrank.eq.0) THEN
         OPEN(1,file=trim(path) // '/trenergy.txt',position='append')
         WRITE(1,FMT='(E13.6,E22.14)') (t-1)*dt,enet
         CLOSE(1)
      ENDIF

      CALL gws%free_real_htmp(r)

      RETURN
      END SUBROUTINE trapenergy

!*****************************************************************
      SUBROUTINE rotenergy(a,b,vlinx,vliny,alpha,t,dt,path)
!-----------------------------------------------------------------
!
! Computes the energy associated to the rotation of the
! condensate (including the effective repulsive centrifugal
! potential), in the rotating frame of reference,
! Erot = - < z* Omega . (r x P) z > = - Omega . L, where P is
! the momentum operator and L the angular momentum. With this
! sign, the energy in the rotating frame H + Etrap + Erot is
! conserved (H is the energy in 'balance.txt', Etrap the energy
! in 'trenergy.txt'). This subroutine writes the output to a file.
!
! Output file contains:
! 'rotenergy.txt': time, rot. energy - < z* Omega.(r x P) z >
!
! Parameters
!     a     : real part of the wavefunction in Fourier space
!     b     : imaginary part of the wavefunction in Fourier space
!     vlinx : omegaz times the linear x-coordinate (with the fix
!             for the b.c.), in real space (not normalized)
!     vliny : omegaz times the linear y-coordinate
!     alpha : amplitude of the alpha coefficient in GPE
!     t     : number of time steps made
!     dt    : time step
!     path  : path for the output
!
      USE fprecision
      USE commtypes
      USE mpivars
      USE grid
      USE kes
      USE ali
      USE fft
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a,b
      REAL(KIND=GP), INTENT(IN), DIMENSION(nx,ny,ksta:kend) :: vlinx,vliny
      REAL(KIND=GP), INTENT(IN)    :: dt,alpha
      INTEGER, INTENT(IN)          :: t
      CHARACTER(len=*), INTENT(IN) :: path
      COMPLEX(KIND=GP), POINTER, DIMENSION(:,:,:) :: c1
      REAL(KIND=GP), POINTER, DIMENSION(:,:,:) :: r1,r2,dx1,dx2,dy1,dy2
      DOUBLE PRECISION    :: tmp,tmq,tmr,ene,enet
      INTEGER             :: i,j,k
      LOGICAL             :: bret

      CALL gws%get_complex_htmp(c1,bret)
      CALL gws%get_real_htmp(r1,bret)
      CALL gws%get_real_htmp(r2,bret)
      CALL gws%get_real_htmp(dx1,bret)
      CALL gws%get_real_htmp(dx2,bret)
      CALL gws%get_real_htmp(dy1,bret)
      CALL gws%get_real_htmp(dy2,bret)
!
! Computes < z* Omega.(r x P) z > (the linear ramps already
! include the rotation rate omegaz)
!
      CALL copy3(a,c1)
      CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      CALL copy3(b,c1)
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
      tmp = -2.0D0*alpha/ &
            (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**4
!$omp parallel do collapse(2) private (i,tmq,tmr) reduction(+:ene)
      DO k = ksta,kend
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
         OPEN(1,file=trim(path) // '/rotenergy.txt',position='append')
         WRITE(1,FMT='(E13.6,E22.14)') (t-1)*dt,enet
         CLOSE(1)
      ENDIF

      CALL gws%free_real_htmp(dy2)
      CALL gws%free_real_htmp(dy1)
      CALL gws%free_real_htmp(dx2)
      CALL gws%free_real_htmp(dx1)
      CALL gws%free_real_htmp(r2)
      CALL gws%free_real_htmp(r1)
      CALL gws%free_complex_htmp(c1)

      END SUBROUTINE rotenergy

!*****************************************************************
      SUBROUTINE gperealspecperp(a,b,alpha,beta,omegag,path,nmb)
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
! 'qspecperp.XXX.txt'    : k, Equa(k)
!   [Equa ~ (zturnre*grad(zre)+zturnim*grad(zim))^2]
! 'kincspecperp.XXX.txt' : k, Einc(k)
! 'kcomspecperp.XXX.txt' : k, Ecom(k)
!   [Ekin ~ (zturnre*grad(zim)-zturnim*grad(zre))^2, decomposed into]
!   [incompressible and compressible components.]
!
! Parameters
!     a : real part of the wavefunction in Fourier space
!     b : imaginary part of the wavefunction in Fourier space
!     alpha,beta,omegag: constants of the GPE
!     path: path for the output
!     nmb: the extension used when writting the file
!
      USE kes
      USE grid
      USE mpivars
      USE filefmt
      USE boxsize
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a,b
      REAL(KIND=GP), INTENT(IN)    :: alpha,beta,omegag
      DOUBLE PRECISION, DIMENSION(nmax/2+1)    :: Eint,Equa,Einc,Ecom
      INTEGER                      :: i
      CHARACTER(len=*), INTENT(IN) :: path,nmb

!
! Computes all the energy  spectra
!
      CALL gperealspecperpc(a,b,alpha,beta,omegag,Eint,Equa,Einc,Ecom)
!
! Exports the energy spectrum to a file
!
      IF (myrank.eq.0) THEN
         OPEN(1,file=trim(path) // '/intspecperp.' // nmb // '.txt')
         DO i=1,nmax/2+1
            WRITE(1,FMT='(E13.6,E23.15)') Dkk*(i-1),Eint(i)/Dkk
         END DO
         CLOSE(1)
         OPEN(1,file=trim(path) // '/qspecperp.' // nmb // '.txt')
         DO i=1,nmax/2+1
            WRITE(1,FMT='(E13.6,E23.15)') Dkk*(i-1),Equa(i)/Dkk
         END DO
         CLOSE(1)
         OPEN(1,file=trim(path) // '/kincspecperp.' // nmb // '.txt')
         DO i=1,nmax/2+1
            WRITE(1,FMT='(E13.6,E23.15)') Dkk*(i-1),Einc(i)/Dkk
         END DO
         CLOSE(1)
         OPEN(1,file=trim(path) // '/kcomspecperp.' // nmb // '.txt')
         DO i=1,nmax/2+1
            WRITE(1,FMT='(E13.6,E23.15)') Dkk*(i-1),Ecom(i)/Dkk
         END DO
         CLOSE(1)
      ENDIF

      RETURN
      END SUBROUTINE gperealspecperp

!*****************************************************************
      SUBROUTINE gperealspecperpc(a,b,alpha,beta,omegag,Eint,Equa,Einc,Ecom)
!-----------------------------------------------------------------
!
! Computes the reduced perpendicular spectrum of kinetic, quantum,
! and potential (or internal) energy, as a function of the
! perpendicular wavenumber, returning them.
!
! Parameters
!     a   : real part of the wavefunction in Fourier space
!     b   : imaginary part of the wavefunction in Fourier space
!     alpha,beta,omegag: constants of the GPE
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
      USE mpivars
      USE boxsize
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a,b
      REAL(KIND=GP), INTENT(IN) :: alpha,beta,omegag
      DOUBLE PRECISION, INTENT(OUT), DIMENSION(nmax/2+1) ::Eint,Equa,Einc,Ecom
      DOUBLE PRECISION, DIMENSION(nmax/2+1)        :: Etmp
      DOUBLE PRECISION, DIMENSION(nmax/2+1)        :: Ek1,Ek2
      COMPLEX(KIND=GP), POINTER, DIMENSION(:,:,:) :: c1,c2,c3,c4
      REAL(KIND=GP), POINTER, DIMENSION(:,:,:) :: r1,r2,r3,qua,kin
      REAL(KIND=GP)    :: rmp,rmq
      INTEGER          :: i,j,k,m
      LOGICAL          :: bret

      CALL gws%get_complex_htmp(c1,bret)
      CALL gws%get_complex_htmp(c2,bret)
      CALL gws%get_complex_htmp(c3,bret)
      CALL gws%get_complex_htmp(c4,bret)
      CALL gws%get_real_htmp(r1,bret)
      CALL gws%get_real_htmp(r2,bret)
      CALL gws%get_real_htmp(r3,bret)
      CALL gws%get_real_htmp(qua,bret)
      CALL gws%get_real_htmp(kin,bret)
!
! Transforms the wavefunction to real space
!
      CALL copy3(a,c1)
      CALL copy3(b,c2)
      CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,c2,r2,MPI_COMM_WORLD)
!
! Computes the internal energy spectrum
! Eint = 2.alpha.beta.(|z|^2-rho0)^2
!
      rmp = sqrt(alpha*beta)*omegag/beta
      rmq = sqrt(alpha*beta)/ &
            (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2
!$omp parallel do collapse(2) private (i)
      DO k = ksta,kend
         DO j = 1,ny
            DO i = 1,nx
               r3(i,j,k) = (r1(i,j,k)**2+r2(i,j,k)**2)*rmq-rmp
            END DO
         END DO
      END DO
      CALL fftp3d_real_to_complex(planrc,r3,c1,MPI_COMM_WORLD)
!$omp parallel do collapse(2) private (k)
      DO i = ista,iend   ! This spectrum must be dealiased
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
      CALL zturn(r1,r2,regu*omegag/beta)
!
! Computes the quantum energy spectrum, and
! prepares to compute the kinetic energy spectra
! Equa ~ (zturnre*grad(zre)+zturnim*grad(zim))^2
! Ekin ~ (zturnre*grad(zim)-zturnim*grad(zre))^2
!
      DO m = 1,3
         CALL derivk3(a,c1,m)
         CALL fftp3d_complex_to_real(plancr,c1,r3,MPI_COMM_WORLD)
!$omp parallel do collapse(2) private (i)
         DO k = ksta,kend
            DO j = 1,ny
               DO i = 1,nx
                  qua(i,j,k) =  r1(i,j,k)*r3(i,j,k)
                  kin(i,j,k) = -r2(i,j,k)*r3(i,j,k)
               END DO
            END DO
         END DO
         CALL derivk3(b,c1,m)
         CALL fftp3d_complex_to_real(plancr,c1,r3,MPI_COMM_WORLD)
!$omp parallel do collapse(2) private (i)
         DO k = ksta,kend
            DO j = 1,ny
               DO i = 1,nx
                  qua(i,j,k) = qua(i,j,k)+r2(i,j,k)*r3(i,j,k)
                  kin(i,j,k) = kin(i,j,k)+r1(i,j,k)*r3(i,j,k)
               END DO
            END DO
         END DO
         CALL fftp3d_real_to_complex(planrc,qua,c1,MPI_COMM_WORLD)
         IF (m.eq.1) THEN
            CALL fftp3d_real_to_complex(planrc,kin,c2,MPI_COMM_WORLD)
            CALL specscpec(c1,Ek1,Etmp)
         ELSE IF (m.eq.2) THEN
            CALL fftp3d_real_to_complex(planrc,kin,c3,MPI_COMM_WORLD)
            CALL specscpec(c1,Ek2,Etmp)
            Ek1 = Ek1+Ek2
         ELSE
            CALL fftp3d_real_to_complex(planrc,kin,c4,MPI_COMM_WORLD)
            CALL specscpec(c1,Ek2,Etmp)
         ENDIF
      END DO
      rmq = 2*alpha**2/ &
            (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2
      Equa = (Ek1+Ek2)*rmq
!
! Computes the compressible and incompressible kinetic energy spectra
!
      CALL gauge3(c2,c3,c4,c1,1)     ! x component
      CALL specscpec(c1,Einc,Etmp)   ! incompressible
      CALL saxpby_c(c1,c2,1.0_GP,c1,-1.0_GP)
      CALL specscpec(c1,Ecom,Etmp)   ! compressible

      CALL gauge3(c2,c3,c4,c1,2)     ! y component
      CALL specscpec(c1,Ek1,Etmp)    ! incompressible
      CALL saxpby_c(c1,c3,1.0_GP,c1,-1.0_GP)
      CALL specscpec(c1,Ek2,Etmp)    ! compressible
      Einc = Einc+Ek1
      Ecom = Ecom+Ek2

      CALL gauge3(c2,c3,c4,c1,3)     ! z component
      CALL specscpec(c1,Ek1,Etmp)    ! incompressible
      CALL saxpby_c(c1,c4,1.0_GP,c1,-1.0_GP)
      CALL specscpec(c1,Ek2,Etmp)    ! compressible
      Einc = (Einc+Ek1)*rmq
      Ecom = (Ecom+Ek2)*rmq

      CALL gws%free_real_htmp(kin)
      CALL gws%free_real_htmp(qua)
      CALL gws%free_real_htmp(r3)
      CALL gws%free_real_htmp(r2)
      CALL gws%free_real_htmp(r1)
      CALL gws%free_complex_htmp(c4)
      CALL gws%free_complex_htmp(c3)
      CALL gws%free_complex_htmp(c2)
      CALL gws%free_complex_htmp(c1)

      RETURN
      END SUBROUTINE gperealspecperpc

END MODULE pseudospec_rgpe
