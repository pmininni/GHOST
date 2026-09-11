!=================================================================
! PSEUDOSPECTRAL modules
!
! CONTAINS:
!      MODULE pseudospec_quantum
!      MODULE pseudospec_gpe
!
! Subroutines to compute spatial derivatives and nonlinear
! terms in the GPE and ARGL equations in 3D using a
! pseudo-spectral method. You should use the FFTPLANS
! and MPIVARS modules (see the file 'fftp_mod.f90') in each
! program that calls any of the subroutines in this file.
!
! The module pseudospec_quantum has the routines used in the
! evolution equations (kernels with the dual host/device
! directives, temporaries from the workspace pool), and the
! module pseudospec_gpe the diagnostics (computed in the host,
! with host-only temporaries). The constants of the equations
! (alpha, beta, omegag) are passed as arguments; the old 'hbar'
! module is not used. See pseudospec3D_hd.f90 for the methodology.
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

MODULE pseudospec_quantum
   USE pseudospec_fluid
   USE class_GWorkspace3D, ONLY: gws
   USE gdevice, ONLY: gdev_active
   CONTAINS

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
!
      USE fprecision
      USE commtypes
      USE mpivars
      USE grid
      USE kes
      USE ali
      USE fft
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend)  :: a,b
      REAL(KIND=GP), INTENT(OUT), DIMENSION(nx,ny,ksta:kend)    :: r
      COMPLEX(KIND=GP), POINTER, DIMENSION(:,:,:) :: c1
      REAL(KIND=GP), POINTER, DIMENSION(:,:,:) :: r1
      REAL(KIND=GP)       :: rmp
      INTEGER, INTENT(IN) :: dealias
      INTEGER :: i,j,k
      LOGICAL :: bret

      CALL gws%get_complex_tmp(c1,bret)
      CALL gws%get_real_tmp(r1,bret)
!
! Computes the square of the real part of the wavefunction
!
      CALL copy3(a,c1)
      CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
#if defined(GHOST_GPU)
!$omp target teams distribute parallel do collapse(3) if(target: gdev_active)
#else
!$omp parallel do collapse(2) private (i)
#endif
      DO k = ksta,kend
         DO j = 1,ny
            DO i = 1,nx
               r(i,j,k) = r1(i,j,k)**2
            END DO
         END DO
      END DO
!
! Computes the square of the imaginary part of the wavefunction
!
      CALL copy3(b,c1)
      CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      rmp = 1.0_GP/ &
            (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2
#if defined(GHOST_GPU)
!$omp target teams distribute parallel do collapse(3) if(target: gdev_active)
#else
!$omp parallel do collapse(2) private (i)
#endif
      DO k = ksta,kend
         DO j = 1,ny
            DO i = 1,nx
               r(i,j,k) = (r(i,j,k)+r1(i,j,k)**2)*rmp
            END DO
         END DO
      END DO
!
! Dealiases the result and returns to real space
!
      IF (dealias.eq.1) THEN
         CALL fftp3d_real_to_complex(planrc,r,c1,MPI_COMM_WORLD)
#if defined(GHOST_GPU)
!$omp target teams distribute parallel do collapse(3) if(target: gdev_active)
#else
!$omp parallel do collapse(2) private (k)
#endif
         DO i = ista,iend
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

      CALL gws%free_real_tmp(r1)
      CALL gws%free_complex_tmp(c1)

      RETURN
      END SUBROUTINE squareabs

!*****************************************************************
      SUBROUTINE nonlgpe(r,a,b)
!-----------------------------------------------------------------
!
! Computes Z.|Z|^2 in real space (or, in general, the product of
! the scalar r in real space by the field a).
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
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN),  DIMENSION(nz,ny,ista:iend) :: a
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(nz,ny,ista:iend) :: b
      REAL(KIND=GP), INTENT(IN), DIMENSION(nx,ny,ksta:kend)     :: r
      REAL(KIND=GP), POINTER, DIMENSION(:,:,:) :: r1
      REAL(KIND=GP)    :: rmp
      INTEGER :: i,j,k
      LOGICAL :: bret

      CALL gws%get_real_tmp(r1,bret)
!
! Computes Z.|Z|^2
!
      CALL copy3(a,b)
      CALL fftp3d_complex_to_real(plancr,b,r1,MPI_COMM_WORLD)
      rmp = 1.0_GP/ &
            (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2
#if defined(GHOST_GPU)
!$omp target teams distribute parallel do collapse(3) if(target: gdev_active)
#else
!$omp parallel do collapse(2) private (i)
#endif
      DO k = ksta,kend
         DO j = 1,ny
            DO i = 1,nx
               r1(i,j,k) = r(i,j,k)*r1(i,j,k)*rmp
            END DO
         END DO
      END DO
      CALL fftp3d_real_to_complex(planrc,r1,b,MPI_COMM_WORLD)

      CALL gws%free_real_tmp(r1)

      RETURN
      END SUBROUTINE nonlgpe

!*****************************************************************
      SUBROUTINE gpe_mass(a,b,mass)
!-----------------------------------------------------------------
!
! Computes the mass (mean density) of the wavefunction,
! <|Z|^2> = sum_k (|a_k|^2+|b_k|^2)/N^6, in Fourier space, and
! returns it in all the tasks. Same normalization as the sum of
! the variances of a and b (routine 'variance'). Used in the
! evolution (renormalization in finite temperature runs), so
! the reduction runs on the device in offload builds.
!
! Parameters
!     a   : real part of the wavefunction in Fourier space
!     b   : imaginary part of the wavefunction in Fourier space
!     mass: mean density [output, in all tasks]
!
      USE fprecision
      USE commtypes
      USE mpivars
      USE grid
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a,b
      DOUBLE PRECISION, INTENT(OUT) :: mass
      DOUBLE PRECISION              :: bloc
      REAL(KIND=GP)                 :: tmp,w
      INTEGER                       :: i,j,k

      bloc = 0.0D0
      tmp = 1.0_GP/ &
            (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2
      ! The plane kx = 0 (i = 1) is counted once, the others twice
#if defined(GHOST_GPU)
!$omp target teams distribute parallel do collapse(3) if(target: gdev_active) &
!$omp   private(w) reduction(+:bloc)
#else
!$omp parallel do collapse(2) private (k,w) reduction(+:bloc)
#endif
      DO i = ista,iend
         DO j = 1,ny
            DO k = 1,nz
               w = 2.0_GP
               IF (i.eq.1) w = 1.0_GP
               bloc = bloc + w*tmp*(real(a(k,j,i))**2+aimag(a(k,j,i))**2 + &
                                    real(b(k,j,i))**2+aimag(b(k,j,i))**2)
            END DO
         END DO
      END DO
      CALL MPI_ALLREDUCE(bloc,mass,1,MPI_DOUBLE_PRECISION,MPI_SUM, &
                         MPI_COMM_WORLD,ierr)

      RETURN
      END SUBROUTINE gpe_mass

END MODULE pseudospec_quantum


!=================================================================
! Diagnostics of the GPE and ARGL equations (computed in the host)
!=================================================================
MODULE pseudospec_gpe
   USE fprecision
   USE pseudospec_fluid
   USE pseudospec_scalar
   USE pseudospec_magnetic
   USE pseudospec_quantum, ONLY: squareabs
   USE class_GWorkspace3D, ONLY: gws
   IMPLICIT NONE
   ! Regularization of the divisions by |z|^2
   REAL(KIND=GP), PARAMETER, PUBLIC :: regu = 1.e-20_GP
   CONTAINS

!*****************************************************************
      SUBROUTINE gpecheck(a,b,alpha,beta,omegag,t,dt,path)
!-----------------------------------------------------------------
!
! Computes the mass, the kinetic plus quantum energy Ekq, the
! potential energy, and the quartic term in the energy:
!    Ekq    = 2.alpha^2.|grad(z)|^2
!    Equart = alpha.beta.|z|^4
! The potential energy then is
!    Epot = Equart-2*alpha.omegag.mass+alpha.omegag^2/beta
! This quantity should be zero in the condensate. The total
! energy is simply:
!    E = Ekq + Epot
! The energy from the Hamiltonian is:
!    H = Ekq + Equart
! and the free energy can be computed as H - mu.N where
!    mu.N/V = 2.alpha.omegag.mass = c^2 mass
! where mass is the mass density.
! The results are written to a file by the first node.
!
! Output files contain:
! 'balance.txt': time, mass, kinetic+quantum en., potential en., quartic en.
!   [Ekq = 2.alpha^2.|grad(z)|^2, Equart = alpha.beta.|z|^4, and the      ]
!   [pot. energy is Epot = Equart-2*alpha.omegag.mass+alpha.omegag^2/beta.]
!   [Note this output replaces all 'balance.txt' files in quantum solvers.]
!
! Parameters
!     a : input matrix with the real part of the wavefunction
!     b : input matrix with the imaginary part of the wavefunction
!     alpha,beta,omegag: constants of the GPE
!     t : number of time steps made
!     dt: time step
!     path: path for the output
!
      USE fprecision
      USE commtypes
      USE grid
      USE mpivars
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a,b
      REAL(KIND=GP), INTENT(IN)     :: alpha,beta,omegag,dt
      INTEGER, INTENT(IN)           :: t
      CHARACTER(len=*), INTENT(IN)  :: path
      REAL(KIND=GP), POINTER, DIMENSION(:,:,:) :: r1
      DOUBLE PRECISION    :: mass,ekq
      DOUBLE PRECISION    :: tmp,tmq
      INTEGER             :: i,j,k
      LOGICAL             :: bret

      CALL gws%get_real_htmp(r1,bret)
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
!$omp parallel do collapse(2) private (i) reduction(+:tmp)
      DO k = ksta,kend
         DO j = 1,ny
            DO i = 1,nx
               tmp = tmp+r1(i,j,k)**2
            END DO
         END DO
      END DO
      CALL MPI_REDUCE(tmp,tmq,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
      IF (myrank.eq.0) THEN
         tmq = alpha*beta*tmq/ &
          (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**3
         tmp = tmq + alpha*(omegag**2/beta - 2*omegag*mass)
      ENDIF
!
! Creates a external file to store the results
!
      IF (myrank.eq.0) THEN
         OPEN(1,file=trim(path) // '/balance.txt',position='append')
         WRITE(1,10) (t-1)*dt,mass,ekq,tmp,tmq
   10    FORMAT( E13.6,E22.14,E22.14,E22.14,E22.14 )
         CLOSE(1)
      ENDIF

      CALL gws%free_real_htmp(r1)

      RETURN
      END SUBROUTINE gpecheck

!*****************************************************************
      SUBROUTINE momentum(a,b,alpha,t,dt,path)
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
!     alpha: constant of the GPE
!     t : number of time steps made
!     dt: time step
!     path: path for the output
!
      USE fprecision
      USE grid
      USE mpivars
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a,b
      REAL(KIND=GP), INTENT(IN)    :: alpha,dt
      INTEGER, INTENT(IN)          :: t
      CHARACTER(len=*), INTENT(IN) :: path
      COMPLEX(KIND=GP), POINTER, DIMENSION(:,:,:) :: C1,C2
      DOUBLE PRECISION    :: tmp,tmq
      DOUBLE PRECISION    :: jx,jy,jz
      LOGICAL             :: bret

      CALL gws%get_complex_htmp(C1,bret)
      CALL gws%get_complex_htmp(C2,bret)

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
         OPEN(1,file=trim(path) // '/momentum.txt',position='append')
         WRITE(1,FMT='(E13.6,E22.14,E22.14,E22.14)') (t-1)*dt,jx,jy,jz
         CLOSE(1)
      ENDIF

      CALL gws%free_complex_htmp(C2)
      CALL gws%free_complex_htmp(C1)

      RETURN
      END SUBROUTINE momentum

!*****************************************************************
      SUBROUTINE gpemassspec(a,b,path,nmb)
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
!     path: path for the output
!     nmb: the extension used when writting the file
!
      USE fprecision
      USE commtypes
      USE kes
      USE grid
      USE mpivars
      USE filefmt
      USE boxsize
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a,b
      DOUBLE PRECISION, DIMENSION(nmax/2+1)               :: Ek,Ektot
      DOUBLE PRECISION :: tmq
      REAL(KIND=GP)    :: rmp
      INTEGER          :: i,j,k
      INTEGER          :: kmn
      CHARACTER(len=*), INTENT(IN) :: path,nmb

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
!$omp parallel private (k,kmn,tmq) reduction(+:Ek)
!$omp do
         DO j = 1,ny
            DO k = 1,nz
               kmn = int(sqrt(kk2(k,j,1))/Dkk+1.501)
               IF ((kmn.gt.0).and.(kmn.le.nmax/2+1)) THEN
                  tmq = rmp*(abs(a(k,j,1))**2+abs(b(k,j,1))**2)
                  Ek(kmn) = Ek(kmn)+tmq
               ENDIF
            END DO
         END DO
!$omp end do
!$omp do collapse(2)
         DO i = 2,iend
            DO j = 1,ny
               DO k = 1,nz
                  kmn = int(sqrt(kk2(k,j,i))/Dkk+1.501)
                  IF ((kmn.gt.0).and.(kmn.le.nmax/2+1)) THEN
                     tmq = 2*rmp*(abs(a(k,j,i))**2+abs(b(k,j,i))**2)
                     Ek(kmn) = Ek(kmn)+tmq
                  ENDIF
               END DO
            END DO
         END DO
!$omp end do
!$omp end parallel
      ELSE
!$omp parallel do collapse(2) private (k,kmn,tmq) reduction(+:Ek)
         DO i = ista,iend
            DO j = 1,ny
               DO k = 1,nz
                  kmn = int(sqrt(kk2(k,j,i))/Dkk+1.501)
                  IF ((kmn.gt.0).and.(kmn.le.nmax/2+1)) THEN
                     tmq = 2*rmp*(abs(a(k,j,i))**2+abs(b(k,j,i))**2)
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
         OPEN(1,file=trim(path) // '/massspectrum.' // nmb // '.txt')
         DO i=1,nmax/2+1
            WRITE(1,FMT='(E13.6,E23.15)') Dkk*(i-1),Ektot(i)/Dkk
         END DO
         CLOSE(1)
      ENDIF

      RETURN
      END SUBROUTINE gpemassspec

!*****************************************************************
      SUBROUTINE gperealspec(a,b,alpha,beta,omegag,path,nmb)
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
      DOUBLE PRECISION, DIMENSION(nmax/2+1)        :: Eint,Equa,Einc,Ecom
      INTEGER                      :: i
      CHARACTER(len=*), INTENT(IN) :: path,nmb

!
! Computes all the energy  spectra
!
      CALL gperealspecc(a,b,alpha,beta,omegag,Eint,Equa,Einc,Ecom)
!
! Exports the energy spectrum to a file
!
      IF (myrank.eq.0) THEN
         OPEN(1,file=trim(path) // '/intspectrum.' // nmb // '.txt')
         DO i=1,nmax/2+1
            WRITE(1,FMT='(E13.6,E23.15)') Dkk*(i-1),Eint(i)/Dkk
         END DO
         CLOSE(1)
         OPEN(1,file=trim(path) // '/qspectrum.' // nmb // '.txt')
         DO i=1,nmax/2+1
            WRITE(1,FMT='(E13.6,E23.15)') Dkk*(i-1),Equa(i)/Dkk
         END DO
         CLOSE(1)
         OPEN(1,file=trim(path) // '/kincspectrum.' // nmb // '.txt')
         DO i=1,nmax/2+1
            WRITE(1,FMT='(E13.6,E23.15)') Dkk*(i-1),Einc(i)/Dkk
         END DO
         CLOSE(1)
         OPEN(1,file=trim(path) // '/kcomspectrum.' // nmb // '.txt')
         DO i=1,nmax/2+1
            WRITE(1,FMT='(E13.6,E23.15)') Dkk*(i-1),Ecom(i)/Dkk
         END DO
         CLOSE(1)
      ENDIF

      RETURN
      END SUBROUTINE gperealspec

!*****************************************************************
      SUBROUTINE gperealspecc(a,b,alpha,beta,omegag,Eint,Equa,Einc,Ecom)
!-----------------------------------------------------------------
!
! Computes the spectrum of kinetic, quantum, and potential (or
! internal) energy, returning them.
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
      DOUBLE PRECISION, INTENT(OUT), DIMENSION(nmax/2+1) :: Eint,Equa,Einc,Ecom
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
      CALL spectrscc(c1,Eint,1.0_GP)
!
! Computes z/sqrt(|z|^2) inplace
!
      CALL zturn(r1,r2,regu*omegag/beta)
!
! Computes the quantum energy spectrum, and
! prepares to compute the kinetic energy spectra
! Equa ~ (zturnre*grad(zre)+zturnim*grad(zim))^2
! Ekin ~ (zturnre*grad(zim)-zturnim*grad(zre))^2
! (the kinetic energy components go to c2, c3 and c4)
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
            CALL spectrscc(c1,Ek1,1.0_GP)
         ELSE IF (m.eq.2) THEN
            CALL fftp3d_real_to_complex(planrc,kin,c3,MPI_COMM_WORLD)
            CALL spectrscc(c1,Ek2,1.0_GP)
            Ek1 = Ek1+Ek2
         ELSE
            CALL fftp3d_real_to_complex(planrc,kin,c4,MPI_COMM_WORLD)
            CALL spectrscc(c1,Ek2,1.0_GP)
         ENDIF
      END DO
      rmq = 2*alpha**2/ &
            (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2
      Equa = (Ek1+Ek2)*rmq
!
! Computes the compressible and incompressible kinetic energy spectra
!
      CALL gauge3(c2,c3,c4,c1,1)      ! x component
      CALL spectrscc(c1,Einc,1.0_GP)  ! incompressible
      CALL saxpby_c(c1,c2,1.0_GP,c1,-1.0_GP)
      CALL spectrscc(c1,Ecom,1.0_GP)  ! compressible

      CALL gauge3(c2,c3,c4,c1,2)      ! y component
      CALL spectrscc(c1,Ek1,1.0_GP)   ! incompressible
      CALL saxpby_c(c1,c3,1.0_GP,c1,-1.0_GP)
      CALL spectrscc(c1,Ek2,1.0_GP)   ! compressible
      Einc = Einc+Ek1
      Ecom = Ecom+Ek2

      CALL gauge3(c2,c3,c4,c1,3)      ! z component
      CALL spectrscc(c1,Ek1,1.0_GP)   ! incompressible
      CALL saxpby_c(c1,c4,1.0_GP,c1,-1.0_GP)
      CALL spectrscc(c1,Ek2,1.0_GP)   ! compressible
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
      END SUBROUTINE gperealspecc

!*****************************************************************
      SUBROUTINE gperealtrans(dt,io,qo,ko,co,in,qn,kn,cn,path,nmb)
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
!     path: path for the output
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
      CHARACTER(len=*), INTENT(IN) :: path,nmb

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
         OPEN(1,file=trim(path) // '/inttransfer.' // nmb // '.txt')
         DO i=1,nmax/2+1
            WRITE(1,FMT='(E13.6,E23.15)') Dkk*(i-1),in(i)/Dkk
         END DO
         CLOSE(1)
         OPEN(1,file=trim(path) // '/qtransfer.' // nmb // '.txt')
         DO i=1,nmax/2+1
            WRITE(1,FMT='(E13.6,E23.15)') Dkk*(i-1),qn(i)/Dkk
         END DO
         CLOSE(1)
         OPEN(1,file=trim(path) // '/kinctransfer.' // nmb // '.txt')
         DO i=1,nmax/2+1
            WRITE(1,FMT='(E13.6,E23.15)') Dkk*(i-1),kn(i)/Dkk
         END DO
         CLOSE(1)
         OPEN(1,file=trim(path) // '/kcomtransfer.' // nmb // '.txt')
         DO i=1,nmax/2+1
            WRITE(1,FMT='(E13.6,E23.15)') Dkk*(i-1),cn(i)/Dkk
         END DO
         CLOSE(1)
      ENDIF

      RETURN
      END SUBROUTINE gperealtrans

!*****************************************************************
      SUBROUTINE zturn(ra,rb,reg)
!-----------------------------------------------------------------
!
! Computes the real and imaginary parts of z/sqrt(|z|^2) in
! place (i.e., the input is destroyed, and replaced by the
! output). Note zturn in GHOST is the complex conjugate of the
! quantity called zturn in Brachet's GPE codes (including TYGRES).
!
! Parameters
!     ra : real part of the wavefunction in real space
!     rb : imaginary part of the wavefunction in real space
!     reg: regularization of the density, regu*omegag/beta
!
      USE fprecision
      USE grid
      USE mpivars
      IMPLICIT NONE

      REAL(KIND=GP), INTENT(INOUT), DIMENSION(nx,ny,ksta:kend) :: ra,rb
      REAL(KIND=GP), INTENT(IN) :: reg
      REAL(KIND=GP)    :: rmp,rms
      INTEGER          :: i,j,k

      rms = (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2*reg
!
! Computes zbar/sqrt(|z|^2)
!
!$omp parallel do collapse(2) private (i,rmp)
      DO k = ksta,kend
         DO j = 1,ny
            DO i = 1,nx
               rmp = 1.0_GP/sqrt(ra(i,j,k)**2+rb(i,j,k)**2+rms)
               ra(i,j,k) = ra(i,j,k)*rmp
               rb(i,j,k) = rb(i,j,k)*rmp
            END DO
         END DO
      END DO

      RETURN
      END SUBROUTINE zturn

!*******************************************************************
      SUBROUTINE gpekfield(a,b,beta,omegag,c,d,e,f,g,h)
!-------------------------------------------------------------------
!
! Computes all the components of the compressible and incompressible
! parts of sqrt(rho)*v. These quantities must be computed in real
! space and then transformed to Fourier space.
!
! Parameters
!     a : real part of the wavefunction in Fourier space
!     b : imaginary part of the wavefunction in Fourier space
!     beta,omegag: constants of the GPE
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
      USE mpivars
      USE filefmt
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN),  DIMENSION(nz,ny,ista:iend) :: a,b
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(nz,ny,ista:iend) :: c,d
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(nz,ny,ista:iend) :: e,f
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(nz,ny,ista:iend) :: g,h
      REAL(KIND=GP), INTENT(IN) :: beta,omegag
      COMPLEX(KIND=GP), POINTER, DIMENSION(:,:,:) :: c1,c2,c3,c4
      REAL(KIND=GP), POINTER, DIMENSION(:,:,:) :: r1,r2,r3,kin
      INTEGER          :: i,j,k,m
      LOGICAL          :: bret

      CALL gws%get_complex_htmp(c1,bret)
      CALL gws%get_complex_htmp(c2,bret)
      CALL gws%get_complex_htmp(c3,bret)
      CALL gws%get_complex_htmp(c4,bret)
      CALL gws%get_real_htmp(r1,bret)
      CALL gws%get_real_htmp(r2,bret)
      CALL gws%get_real_htmp(r3,bret)
      CALL gws%get_real_htmp(kin,bret)
!
! Transforms the wavefunction to real space
!
      CALL copy3(a,c1)
      CALL copy3(b,c2)
      CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,c2,r2,MPI_COMM_WORLD)
!
! Computes z/sqrt(|z|^2) inplace
!
      CALL zturn(r1,r2,regu*omegag/beta)
!
! Prepares to compute the kinetic energy spectra
! Ekin ~ (zturnre*grad(zim)-zturnim*grad(zre))^2
!
      DO m = 1,3
         CALL derivk3(a,c1,m)
         CALL fftp3d_complex_to_real(plancr,c1,r3,MPI_COMM_WORLD)
!$omp parallel do collapse(2) private (i)
         DO k = ksta,kend
            DO j = 1,ny
               DO i = 1,nx
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
                  kin(i,j,k) = kin(i,j,k)+r1(i,j,k)*r3(i,j,k)
               END DO
            END DO
         END DO
         IF (m.eq.1) THEN
            CALL fftp3d_real_to_complex(planrc,kin,c2,MPI_COMM_WORLD)
         ELSE IF (m.eq.2) THEN
            CALL fftp3d_real_to_complex(planrc,kin,c3,MPI_COMM_WORLD)
         ELSE
            CALL fftp3d_real_to_complex(planrc,kin,c4,MPI_COMM_WORLD)
         ENDIF
      END DO
!
! Computes the compressible and incompressible parts
!
      CALL gauge3(c2,c3,c4,c1,1)      ! x component
      CALL copy3(c1,c)                ! incompressible
      CALL saxpby_c(d,c2,1.0_GP,c1,-1.0_GP) ! compressible

      CALL gauge3(c2,c3,c4,c1,2)      ! y component
      CALL copy3(c1,e)                ! incompressible
      CALL saxpby_c(f,c3,1.0_GP,c1,-1.0_GP) ! compressible

      CALL gauge3(c2,c3,c4,c1,3)      ! z component
      CALL copy3(c1,g)                ! incompressible
      CALL saxpby_c(h,c4,1.0_GP,c1,-1.0_GP) ! compressible

      CALL gws%free_real_htmp(kin)
      CALL gws%free_real_htmp(r3)
      CALL gws%free_real_htmp(r2)
      CALL gws%free_real_htmp(r1)
      CALL gws%free_complex_htmp(c4)
      CALL gws%free_complex_htmp(c3)
      CALL gws%free_complex_htmp(c2)
      CALL gws%free_complex_htmp(c1)

      RETURN
      END SUBROUTINE gpekfield

!***********************************************************************
      SUBROUTINE gpehelicity(a,b,alpha,beta,omegag,t,dt,path)
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
!     alpha,beta,omegag: constants of the GPE
!     t    : number of time steps made
!     dt   : time step
!     path : path for the output
!
      USE fprecision
      USE commtypes
      USE kes
      USE fft
      USE ali
      USE grid
      USE mpivars
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a,b
      REAL(KIND=GP), INTENT(IN)          :: alpha,beta,omegag,dt
      INTEGER, INTENT(IN)                :: t
      CHARACTER(len=*), INTENT(IN)       :: path
      COMPLEX(KIND=GP), POINTER, DIMENSION(:,:,:) :: c1,c2,c3,c4,c5,c6
      COMPLEX(KIND=GP), POINTER, DIMENSION(:,:,:) :: c7,c8,c9
      DOUBLE PRECISION    :: Htot1,Htot2
      LOGICAL             :: bret

      CALL gws%get_complex_htmp(c1,bret)
      CALL gws%get_complex_htmp(c2,bret)
      CALL gws%get_complex_htmp(c3,bret)
      CALL gws%get_complex_htmp(c4,bret)
      CALL gws%get_complex_htmp(c5,bret)
      CALL gws%get_complex_htmp(c6,bret)
      CALL gws%get_complex_htmp(c7,bret)
      CALL gws%get_complex_htmp(c8,bret)
      CALL gws%get_complex_htmp(c9,bret)
!
! Regularized parallel velocity (c2,c3,c4) and velocity (c5,c6,c1)
!
      CALL gpevfields(a,b,alpha,beta,omegag,c2,c3,c4,c5,c6,c1)
!
! Computes the helicity
!
      CALL helicity(c5,c6,c1,Htot2)
!
! Computes the regularized helicity, v_parallel.curl(v)
!
      CALL rotor3(c6,c1,c7,1)
      CALL rotor3(c5,c1,c8,2)
      CALL rotor3(c5,c6,c9,3)
      CALL cross(c2,c3,c4,c7,c8,c9,Htot1,1)
!
! Writes the result to a file
!
      IF (myrank.eq.0) THEN
         OPEN(1,file=trim(path) // '/helicity.txt',position='append')
         WRITE(1,20) (t-1)*dt,Htot1,Htot2
20    FORMAT( E13.6,E22.14,E22.14 )
         CLOSE(1)
      ENDIF

      CALL gws%free_complex_htmp(c9)
      CALL gws%free_complex_htmp(c8)
      CALL gws%free_complex_htmp(c7)
      CALL gws%free_complex_htmp(c6)
      CALL gws%free_complex_htmp(c5)
      CALL gws%free_complex_htmp(c4)
      CALL gws%free_complex_htmp(c3)
      CALL gws%free_complex_htmp(c2)
      CALL gws%free_complex_htmp(c1)

      RETURN
      END SUBROUTINE gpehelicity

!***********************************************************************
      SUBROUTINE gpehelspec(a,b,alpha,beta,omegag,path,nmb)
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
!     alpha,beta,omegag: constants of the GPE
!     path : path for the output
!     nmb  : the extension used when writting the file
!
      USE fprecision
      USE commtypes
      USE kes
      USE fft
      USE ali
      USE grid
      USE boxsize
      USE mpivars
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a,b
      REAL(KIND=GP), INTENT(IN)          :: alpha,beta,omegag
      CHARACTER(len=*), INTENT(IN)       :: path,nmb
      COMPLEX(KIND=GP), POINTER, DIMENSION(:,:,:) :: c1,c2,c3,c4,c5,c6
      DOUBLE PRECISION, DIMENSION(nmax/2+1) :: Htot3,Htot4
      INTEGER             :: i
      LOGICAL             :: bret

      CALL gws%get_complex_htmp(c1,bret)
      CALL gws%get_complex_htmp(c2,bret)
      CALL gws%get_complex_htmp(c3,bret)
      CALL gws%get_complex_htmp(c4,bret)
      CALL gws%get_complex_htmp(c5,bret)
      CALL gws%get_complex_htmp(c6,bret)
!
! Regularized parallel velocity (c2,c3,c4) and velocity (c5,c6,c1)
!
      CALL gpevfields(a,b,alpha,beta,omegag,c2,c3,c4,c5,c6,c1)
!
! Computes the helicity spectra
!
      CALL gpespectrumc(c5,c6,c1,Htot4)
!
! Computes the regularized helicity spectra
!
      CALL crosspecc(c2,c3,c4,c5,c6,c1,Htot3,1.0_GP)
      IF (myrank.eq.0) THEN
         OPEN(1,file=trim(path) // '/hspectrum.' // nmb // '.txt')
         DO i = 1,nmax/2+1
            WRITE(1,30) Dkk*(i-1),Htot3(i)/Dkk,Htot4(i)/Dkk
         END DO
30       FORMAT( E13.6,E23.15,E23.15 )
         CLOSE(1)
      ENDIF

      CALL gws%free_complex_htmp(c6)
      CALL gws%free_complex_htmp(c5)
      CALL gws%free_complex_htmp(c4)
      CALL gws%free_complex_htmp(c3)
      CALL gws%free_complex_htmp(c2)
      CALL gws%free_complex_htmp(c1)

      RETURN
      END SUBROUTINE gpehelspec

!***********************************************************************
      SUBROUTINE gpevfields(a,b,alpha,beta,omegag,px,py,pz,vx,vy,vz)
!-----------------------------------------------------------------------
!
! Computes the fields used by the helicity diagnostics: the
! regularized parallel velocity v_parallel.e_parallel, with
! e_parallel the unit vector along grad(zre) x grad(zim) and
!    v_parallel = 2.alpha.e_parallel.[(grad(zre).grad)grad(zim)
!                 - (grad(zim).grad)grad(zre)]/(|grad(zre)|^2+|grad(zim)|^2)
! and the velocity v = 2.alpha[zre.grad(zim)-zim.grad(zre)]/|z|^2.
! Both are returned in Fourier space.
!
! Parameters
!     a    : real part of the wavefunction in Fourier space
!     b    : imaginary part of the wavefunction in Fourier space
!     alpha,beta,omegag: constants of the GPE
!     px,py,pz: components of v_parallel.e_parallel [output]
!     vx,vy,vz: components of v [output]
!
      USE fprecision
      USE commtypes
      USE kes
      USE fft
      USE ali
      USE grid
      USE mpivars
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend)  :: a,b
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(nz,ny,ista:iend) :: px,py,pz
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(nz,ny,ista:iend) :: vx,vy,vz
      REAL(KIND=GP), INTENT(IN)          :: alpha,beta,omegag
      COMPLEX(KIND=GP), POINTER, DIMENSION(:,:,:) :: c1,c2,c3
      REAL(KIND=GP), POINTER, DIMENSION(:,:,:)    :: r1,r2,r3,r4
      REAL(KIND=GP), POINTER, DIMENSION(:,:,:)    :: r5,r6,r7,r8
      REAL(KIND=GP), POINTER, DIMENSION(:,:,:)    :: r9,r10,r11,r12
      REAL(KIND=GP)       :: rmp,rms
      REAL(KIND=GP)       :: tmp
      INTEGER             :: i,j,k,m,n
      LOGICAL             :: bret

      CALL gws%get_complex_htmp(c1,bret)
      CALL gws%get_complex_htmp(c2,bret)
      CALL gws%get_complex_htmp(c3,bret)
      CALL gws%get_real_htmp(r1,bret)
      CALL gws%get_real_htmp(r2,bret)
      CALL gws%get_real_htmp(r3,bret)
      CALL gws%get_real_htmp(r4,bret)
      CALL gws%get_real_htmp(r5,bret)
      CALL gws%get_real_htmp(r6,bret)
      CALL gws%get_real_htmp(r7,bret)
      CALL gws%get_real_htmp(r8,bret)
      CALL gws%get_real_htmp(r9,bret)
      CALL gws%get_real_htmp(r10,bret)
      CALL gws%get_real_htmp(r11,bret)
      CALL gws%get_real_htmp(r12,bret)

      rms = (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2* &
            regu*omegag/beta
!
! Transforms the wavefunction to real space
!
      CALL copy3(a,c1)
      CALL copy3(b,c2)
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
!$omp parallel do collapse(2) private (i)
      DO k = ksta,kend
         DO j = 1,ny
            DO i = 1,nx
               r9(i,j,k)  = (r5(i,j,k)*r8(i,j,k)-r7(i,j,k)*r6(i,j,k))*tmp
               r10(i,j,k) = (r7(i,j,k)*r4(i,j,k)-r3(i,j,k)*r8(i,j,k))*tmp
               r11(i,j,k) = (r3(i,j,k)*r6(i,j,k)-r5(i,j,k)*r4(i,j,k))*tmp
               r12(i,j,k) = (real(nx,kind=GP)*real(ny,kind=GP)  *      &
                            real(nz,kind=GP))**2/(r3(i,j,k)**2  +      &
                            r4(i,j,k)**2 + r5(i,j,k)**2         +      &
                            r6(i,j,k)**2 + r7(i,j,k)**2         +      &
                            r8(i,j,k)**2 + rms)
            END DO
         END DO
      END DO
!
! Calculates v_parallel by doing regularizing and the projecting
! The magnitud calculated is:
! 2*alpha e_parallel.((grad(zre).grad)(grad(zim)) -
! (grad(zim).grad)(grad(zre)))/(grad(zre)**2+grad(zim)**2)
! First pass (n=1): +d_m(zre) e_parallel.d_m(grad(zim)), summed
! over m; second pass (n=2): -d_m(zim) e_parallel.d_m(grad(zre)).
!
      DO n = 1,2
         IF (n.eq.1) THEN
            CALL derivk3(b,c1,1)
            CALL derivk3(b,c2,2)
            CALL derivk3(b,c3,3)
         ELSE
            CALL derivk3(a,c1,1)
            CALL derivk3(a,c2,2)
            CALL derivk3(a,c3,3)
         ENDIF
         DO m = 1,3
            IF (n.eq.1) THEN
               CALL derivk3(a,px,m)
            ELSE
               CALL derivk3(b,px,m)
            ENDIF
            CALL derivk3(c1,py,m)
            CALL fftp3d_complex_to_real(plancr,px,r3,MPI_COMM_WORLD)
            CALL fftp3d_complex_to_real(plancr,py,r5,MPI_COMM_WORLD)
            CALL derivk3(c2,py,m)
            CALL fftp3d_complex_to_real(plancr,py,r6,MPI_COMM_WORLD)
            CALL derivk3(c3,py,m)
            CALL fftp3d_complex_to_real(plancr,py,r7,MPI_COMM_WORLD)
            IF ((n.eq.1).and.(m.eq.1)) THEN
!$omp parallel do collapse(2) private (i)
               DO k = ksta,kend
                  DO j = 1,ny
                     DO i = 1,nx
                        r4(i,j,k) = r3(i,j,k)*(r9(i,j,k)*r5(i,j,k) +   &
                                    r10(i,j,k)*r6(i,j,k)           +   &
                                    r11(i,j,k)*r7(i,j,k))*tmp**2
                     END DO
                  END DO
               END DO
            ELSE IF (n.eq.1) THEN
!$omp parallel do collapse(2) private (i)
               DO k = ksta,kend
                  DO j = 1,ny
                     DO i = 1,nx
                        r4(i,j,k) = r4(i,j,k)                      +   &
                                    r3(i,j,k)*(r9(i,j,k)*r5(i,j,k) +   &
                                    r10(i,j,k)*r6(i,j,k)           +   &
                                    r11(i,j,k)*r7(i,j,k))*tmp**2
                     END DO
                  END DO
               END DO
            ELSE
!$omp parallel do collapse(2) private (i)
               DO k = ksta,kend
                  DO j = 1,ny
                     DO i = 1,nx
                        r4(i,j,k) = r4(i,j,k)                      -   &
                                    r3(i,j,k)*(r9(i,j,k)*r5(i,j,k) +   &
                                    r10(i,j,k)*r6(i,j,k)           +   &
                                    r11(i,j,k)*r7(i,j,k))*tmp**2
                     END DO
                  END DO
               END DO
            ENDIF
         END DO
      END DO
!
! Does v_parallel*e_parallel
!
!$omp parallel do collapse(2) private (i,rmp)
      DO k = ksta,kend
         DO j = 1,ny
            DO i = 1,nx
               r4(i,j,k) = 2*alpha*r4(i,j,k)
               rmp = (real(nx,kind=GP)*real(ny,kind=GP)* &
                     real(nz,kind=GP))**2/(r9(i,j,k)**2 + &
                     r10(i,j,k)**2 + r11(i,j,k)**2 + rms)
               r9(i,j,k)  = r4(i,j,k)*r9(i,j,k)*rmp*r12(i,j,k)
               r10(i,j,k) = r4(i,j,k)*r10(i,j,k)*rmp*r12(i,j,k)
               r11(i,j,k) = r4(i,j,k)*r11(i,j,k)*rmp*r12(i,j,k)
            END DO
         END DO
      END DO
      CALL fftp3d_real_to_complex(planrc,r9, px,MPI_COMM_WORLD)
      CALL fftp3d_real_to_complex(planrc,r10,py,MPI_COMM_WORLD)
      CALL fftp3d_real_to_complex(planrc,r11,pz,MPI_COMM_WORLD)
!
! Obtains v by computing 2.alpha[zre*grad(zim)-zim*grad(zre)]
! and dividing by rho.
!
      DO m = 1,3
         CALL derivk3(a,c1,m)
         CALL fftp3d_complex_to_real(plancr,c1,r3,MPI_COMM_WORLD)
         CALL derivk3(b,c1,m)
         CALL fftp3d_complex_to_real(plancr,c1,r4,MPI_COMM_WORLD)
!$omp parallel do collapse(2) private (i,rmp)
         DO k = ksta,kend
            DO j = 1,ny
               DO i = 1,nx
                  rmp = 1.0_GP/(r1(i,j,k)**2+r2(i,j,k)**2+rms)
                  r3(i,j,k) = 2*alpha*(r1(i,j,k)*r4(i,j,k)-      &
                              r2(i,j,k)*r3(i,j,k))*rmp
               END DO
            END DO
         END DO
         IF (m.eq.1) THEN
            CALL fftp3d_real_to_complex(planrc,r3,vx,MPI_COMM_WORLD)
         ELSE IF (m.eq.2) THEN
            CALL fftp3d_real_to_complex(planrc,r3,vy,MPI_COMM_WORLD)
         ELSE
            CALL fftp3d_real_to_complex(planrc,r3,vz,MPI_COMM_WORLD)
         ENDIF
      END DO

      CALL gws%free_real_htmp(r12)
      CALL gws%free_real_htmp(r11)
      CALL gws%free_real_htmp(r10)
      CALL gws%free_real_htmp(r9)
      CALL gws%free_real_htmp(r8)
      CALL gws%free_real_htmp(r7)
      CALL gws%free_real_htmp(r6)
      CALL gws%free_real_htmp(r5)
      CALL gws%free_real_htmp(r4)
      CALL gws%free_real_htmp(r3)
      CALL gws%free_real_htmp(r2)
      CALL gws%free_real_htmp(r1)
      CALL gws%free_complex_htmp(c3)
      CALL gws%free_complex_htmp(c2)
      CALL gws%free_complex_htmp(c1)

      RETURN
      END SUBROUTINE gpevfields

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
      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(nmax/2+1) :: Ek
      DOUBLE PRECISION, INTENT(OUT), DIMENSION(nmax/2+1) :: Hktot
      DOUBLE PRECISION    :: tmq
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a,b,c
      COMPLEX(KIND=GP), POINTER, DIMENSION(:,:,:)              :: c1,c2,c3
      REAL(KIND=GP)       :: tmp
      INTEGER             :: i,j,k
      INTEGER             :: kmn
      LOGICAL             :: bret

      CALL gws%get_complex_htmp(c1,bret)
      CALL gws%get_complex_htmp(c2,bret)
      CALL gws%get_complex_htmp(c3,bret)
!
! Computes the curl of the field
!
      CALL rotor3(b,c,c1,1)
      CALL rotor3(a,c,c2,2)
      CALL rotor3(a,b,c3,3)
!
! Computes the helicity spectrum
!
      tmp = 1.0_GP/ &
            (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2
      DO i = 1,nmax/2+1
         Ek(i) = 0.0D0
         Hktot(i) = 0.0D0
      END DO
      IF (ista.eq.1) THEN
!$omp parallel private (k,kmn,tmq) reduction(+:Ek)
!$omp do
         DO j = 1,ny
            DO k = 1,nz
               kmn = int(sqrt(kk2(k,j,1))/Dkk+1.501)
               IF ((kmn.gt.0).and.(kmn.le.nmax/2+1)) THEN
                  tmq = (real(a(k,j,1)*conjg(c1(k,j,1)))+          &
                         real(b(k,j,1)*conjg(c2(k,j,1)))+          &
                         real(c(k,j,1)*conjg(c3(k,j,1))))*tmp
                  Ek(kmn) = Ek(kmn)+tmq
               ENDIF
            END DO
         END DO
!$omp end do
!$omp do collapse(2)
         DO i = 2,iend
            DO j = 1,ny
               DO k = 1,nz
                  kmn = int(sqrt(kk2(k,j,i))/Dkk+1.501)
                  IF ((kmn.gt.0).and.(kmn.le.nmax/2+1)) THEN
                     tmq = 2*(real(a(k,j,i)*conjg(c1(k,j,i)))+     &
                              real(b(k,j,i)*conjg(c2(k,j,i)))+     &
                              real(c(k,j,i)*conjg(c3(k,j,i))))*tmp
                     Ek(kmn) = Ek(kmn)+tmq
                  ENDIF
              END DO
            END DO
         END DO
!$omp end do
!$omp end parallel
      ELSE
!$omp parallel do collapse(2) private (k,kmn,tmq) reduction(+:Ek)
         DO i = ista,iend
            DO j = 1,ny
               DO k = 1,nz
                  kmn = int(sqrt(kk2(k,j,i))/Dkk+1.501)
                  IF ((kmn.gt.0).and.(kmn.le.nmax/2+1)) THEN
                     tmq = 2*(real(a(k,j,i)*conjg(c1(k,j,i)))+     &
                              real(b(k,j,i)*conjg(c2(k,j,i)))+     &
                              real(c(k,j,i)*conjg(c3(k,j,i))))*tmp
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

      CALL gws%free_complex_htmp(c3)
      CALL gws%free_complex_htmp(c2)
      CALL gws%free_complex_htmp(c1)

      RETURN
      END SUBROUTINE gpespectrumc

!***********************************************************************
      SUBROUTINE gpemomtspec(a,b,alpha,path,nmb)
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
!     alpha: constant of the GPE
!     path: path for the output
!     nmb: the extension used when writting the file
!
      USE fprecision
      USE commtypes
      USE kes
      USE fft
      USE ali
      USE var
      USE grid
      USE mpivars
      USE filefmt
      USE boxsize
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a,b
      REAL(KIND=GP), INTENT(IN)    :: alpha
      CHARACTER(len=*), INTENT(IN) :: path,nmb
      COMPLEX(KIND=GP), POINTER, DIMENSION(:,:,:) :: c1,c2,c3,c4
      REAL(KIND=GP), POINTER, DIMENSION(:,:,:)    :: r1,r2,r3,r4,r5
      DOUBLE PRECISION, DIMENSION(nmax/2+1)       :: Ek,Ektot
      REAL(KIND=GP)    :: rmq
      INTEGER          :: i,j,k,m
      LOGICAL          :: bret

      CALL gws%get_complex_htmp(c1,bret)
      CALL gws%get_complex_htmp(c2,bret)
      CALL gws%get_complex_htmp(c3,bret)
      CALL gws%get_complex_htmp(c4,bret)
      CALL gws%get_real_htmp(r1,bret)
      CALL gws%get_real_htmp(r2,bret)
      CALL gws%get_real_htmp(r3,bret)
      CALL gws%get_real_htmp(r4,bret)
      CALL gws%get_real_htmp(r5,bret)
!
! Transforms the wavefunction to real space
!
      CALL copy3(a,c1)
      CALL copy3(b,c2)
      CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,c2,r2,MPI_COMM_WORLD)
!
! Do derivatives and multiply in real space
!
      rmq = 2*alpha/(real(nx,kind=GP)*real(ny,kind=GP)* &
            real(nz,kind=GP))**2
      DO m = 1,3
         CALL derivk3(a,c1,m)
         CALL fftp3d_complex_to_real(plancr,c1,r3,MPI_COMM_WORLD)
         CALL derivk3(b,c1,m)
         CALL fftp3d_complex_to_real(plancr,c1,r5,MPI_COMM_WORLD)
!$omp parallel do collapse(2) private (i)
         DO k = ksta,kend
            DO j = 1,ny
               DO i = 1,nx
                  r4(i,j,k) = (r2(i,j,k)*r3(i,j,k) - &
                               r1(i,j,k)*r5(i,j,k))*rmq
               END DO
            END DO
         END DO
         IF (m.eq.1) THEN
            CALL fftp3d_real_to_complex(planrc,r4,c2,MPI_COMM_WORLD)
         ELSE IF (m.eq.2) THEN
            CALL fftp3d_real_to_complex(planrc,r4,c3,MPI_COMM_WORLD)
         ELSE
            CALL fftp3d_real_to_complex(planrc,r4,c4,MPI_COMM_WORLD)
         ENDIF
      END DO
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
         OPEN(1,file=trim(path) // '/momtspectrum.' // nmb // '.txt')
         DO i=1,nmax/2+1
            WRITE(1,FMT='(E13.6,E23.15)') Dkk*(i-1),Ektot(i)/Dkk
         END DO
         CLOSE(1)
      ENDIF

      CALL gws%free_real_htmp(r5)
      CALL gws%free_real_htmp(r4)
      CALL gws%free_real_htmp(r3)
      CALL gws%free_real_htmp(r2)
      CALL gws%free_real_htmp(r1)
      CALL gws%free_complex_htmp(c4)
      CALL gws%free_complex_htmp(c3)
      CALL gws%free_complex_htmp(c2)
      CALL gws%free_complex_htmp(c1)

      RETURN
      END SUBROUTINE gpemomtspec

END MODULE pseudospec_gpe
