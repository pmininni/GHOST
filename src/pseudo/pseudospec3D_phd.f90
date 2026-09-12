!=================================================================
! PSEUDOSPECTRAL modules
!
! CONTAINS:
!      MODULE pseudospec_scalar
!      MODULE pseudospec_phd
!
! Extra subroutines to compute the passive/active scalar 
! spectrum, transfer function, and associated global quantities 
! in the HD, MHD, Hall-MHD, and Boussinesq equations when a 
! passive or active scalar is present. You should use the 
! FFTPLANS and MPIVARS modules (see the file 'fftp_mod.f90') in 
! each program that calls any of the subroutines in this file. 
!
! NOTATION: index 'i' is 'x' 
!           index 'j' is 'y'
!           index 'k' is 'z'
!
! 2009 Pablo D. Mininni.
!      Department of Physics, 
!      Facultad de Ciencias Exactas y Naturales.
!      Universidad de Buenos Aires.
!      e-mail: mininni@df.uba.ar 
!=================================================================

MODULE pseudospec_scalar
   USE pseudospec_fluid
   USE class_GWorkspace3D, ONLY: gws
   CONTAINS

!*****************************************************************
      SUBROUTINE advect3(a,b,c,d,e)
!-----------------------------------------------------------------
!
! Three-dimensional inner product -A.grad(B) in 
! real space. The components of the field A are 
! given by the arrays a, b and c, B is a scalar 
! quantity given by d.
!
! Parameters
!     a: input matrix in the x-direction
!     b: input matrix in the y-direction
!     c: input matrix in the z-direction
!     d: input matrix with the scalar
!     e: product (A.grad)B in Fourier space [output]
!
      USE fprecision
      USE commtypes
      USE mpivars
      USE grid
      USE fft
      USE pseudospec_fluid
      USE class_GWorkspace3D, ONLY: gws
      USE gdevice, ONLY: gdev_active
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT (IN), DIMENSION(nz,ny,ista:iend) :: a,b
      COMPLEX(KIND=GP), INTENT (IN), DIMENSION(nz,ny,ista:iend) :: c,d
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(nz,ny,ista:iend) :: e
      COMPLEX(KIND=GP), POINTER, DIMENSION(:,:,:) :: c1,c2
      REAL(KIND=GP), POINTER, DIMENSION(:,:,:)    :: r1,r2,r3
      REAL(KIND=GP)    :: tmp
      INTEGER :: i,j,k,dir
      LOGICAL :: bret

      CALL gws%get_complex_tmp(c1,bret)
      CALL gws%get_complex_tmp(c2,bret)
      CALL gws%get_real_tmp(r1,bret)
      CALL gws%get_real_tmp(r2,bret)
      CALL gws%get_real_tmp(r3,bret)
!
! Computes (A_x.dx)B, (A_y.dy)B and (A_z.dz)B and accumulates
! them in r3. We need -A.grad(B), hence the sign of tmp.
!
      tmp = -1.0_GP/ &
            (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2
      DO dir = 1,3
         IF (dir.eq.1) THEN
            CALL copy3(a,c1)
         ELSE IF (dir.eq.2) THEN
            CALL copy3(b,c1)
         ELSE
            CALL copy3(c,c1)
         ENDIF
         CALL derivk3(d,c2,dir)
         CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
         CALL fftp3d_complex_to_real(plancr,c2,r2,MPI_COMM_WORLD)
         IF (dir.eq.1) THEN
#if defined(GHOST_GPU)
!$omp target teams distribute parallel do collapse(3) if(target: gdev_active)
#else
!$omp parallel do collapse(2) private (i)
#endif
            DO k = ksta,kend
               DO j = 1,ny
                  DO i = 1,nx
                     r3(i,j,k) = r1(i,j,k)*r2(i,j,k)
                  END DO
               END DO
            END DO
         ELSE IF (dir.eq.2) THEN
#if defined(GHOST_GPU)
!$omp target teams distribute parallel do collapse(3) if(target: gdev_active)
#else
!$omp parallel do collapse(2) private (i)
#endif
            DO k = ksta,kend
               DO j = 1,ny
                  DO i = 1,nx
                     r3(i,j,k) = r3(i,j,k)+r1(i,j,k)*r2(i,j,k)
                  END DO
               END DO
            END DO
         ELSE
#if defined(GHOST_GPU)
!$omp target teams distribute parallel do collapse(3) if(target: gdev_active)
#else
!$omp parallel do collapse(2) private (i)
#endif
            DO k = ksta,kend
               DO j = 1,ny
                  DO i = 1,nx
                     r3(i,j,k) = (r3(i,j,k)+r1(i,j,k)*r2(i,j,k))*tmp
                  END DO
               END DO
            END DO
         ENDIF
      END DO
      CALL fftp3d_real_to_complex(planrc,r3,e,MPI_COMM_WORLD)
      CALL gws%free_complex_tmp(c1)
      CALL gws%free_complex_tmp(c2)
      CALL gws%free_real_tmp(r1)
      CALL gws%free_real_tmp(r2)
      CALL gws%free_real_tmp(r3)

      RETURN
      END SUBROUTINE advect3

!*****************************************************************
      SUBROUTINE rhs_scalar3(lapl,adve,f,kappa,out)
!-----------------------------------------------------------------
!
! Assembles the dealiased right-hand side of a scalar equation,
! out = kappa*Laplacian + advection + forcing, inside the
! dealiasing sphere and zero outside.
!
! Parameters
!     lapl : Laplacian of the scalar in Fourier space
!     adve : advection term in Fourier space
!     f    : forcing in Fourier space
!     kappa: diffusivity
!     out  : the right-hand side [output]
!
      USE fprecision
      USE kes
      USE ali
      USE grid
      USE mpivars
      USE gdevice, ONLY: gdev_active
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT (IN), DIMENSION(nz,ny,ista:iend) :: lapl,adve,f
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(nz,ny,ista:iend) :: out
      REAL(KIND=GP), INTENT(IN) :: kappa
      INTEGER :: i,j,k

#if defined(GHOST_GPU)
!$omp target teams distribute parallel do collapse(3) if(target: gdev_active)
#else
!$omp parallel do collapse(2) private (k)
#endif
      DO i = ista,iend
         DO j = 1,ny
            DO k = 1,nz
               IF ((kn2(k,j,i).le.kmax).and.(kn2(k,j,i).ge.tiny)) THEN
                  out(k,j,i) = kappa*lapl(k,j,i) + adve(k,j,i) + f(k,j,i)
               ELSE
                  out(k,j,i) = 0.0_GP
               ENDIF
            END DO
         END DO
      END DO
      RETURN
      END SUBROUTINE rhs_scalar3

!*****************************************************************
      SUBROUTINE variance(a,b,kin)
!-----------------------------------------------------------------
!
! Computes the mean variance of the scalar.
! The output is only valid in the first node.
!
! Parameters
!     a  : input matrix with the scalar
!     d  : at the output contains the variance
!     kin: =0 computes the variance of k^2 times the scalar
!          =1 computes the variance of the scalar
!
      USE fprecision
      USE commtypes
      USE kes
      USE grid
      USE mpivars
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a
      DOUBLE PRECISION, INTENT(OUT) :: b
      DOUBLE PRECISION              :: bloc
      REAL(KIND=GP)                 :: tmp
      INTEGER, INTENT(IN) :: kin
      INTEGER             :: i,j,k

      bloc = 0.0D0
      tmp = 1.0_GP/ &
            (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2
!
! Computes the variance
!
      IF (kin.eq.1) THEN
         IF (ista.eq.1) THEN
!$omp parallel private (k) reduction(+:bloc)
!$omp do
            DO j = 1,ny
               DO CONCURRENT (k=1:nz)
                  bloc = bloc+tmp*abs(a(k,j,1))**2
               END DO
            END DO
!$omp end do
!$omp do collapse(2)
            DO i = 2,iend
               DO j = 1,ny
                  DO CONCURRENT (k=1:nz)
                     bloc = bloc+2*tmp*abs(a(k,j,i))**2
                  END DO
               END DO
            END DO
!$omp end do
!$omp end parallel
         ELSE
!$omp parallel do collapse(2) private (k) reduction(+:bloc)
            DO i = ista,iend
               DO j = 1,ny
                  DO CONCURRENT (k=1:nz)
                     bloc = bloc+2*tmp*abs(a(k,j,i))**2
                  END DO
               END DO
            END DO
         ENDIF
!
! Computes the variance of k^2 times the scalar
!
      ELSE IF (kin.eq.0) THEN
         IF (ista.eq.1) THEN
!$omp parallel private (k) reduction(+:bloc)
!$omp do
            DO j = 1,ny
               DO CONCURRENT (k=1:nz)
                  bloc = bloc+tmp*kk2(k,j,1)*abs(a(k,j,1))**2
               END DO
            END DO
!$omp end do
!$omp do collapse(2)
            DO i = 2,iend
               DO j = 1,ny
                  DO CONCURRENT (k=1:nz)
                     bloc = bloc+2*tmp*kk2(k,j,i)*abs(a(k,j,i))**2
                  END DO
               END DO
            END DO
!$omp end do
!$omp end parallel
         ELSE
!$omp parallel do collapse(2) private (k) reduction(+:bloc)
            DO i = ista,iend
               DO j = 1,ny
                  DO CONCURRENT (k=1:nz)
                     bloc = bloc+2*tmp*kk2(k,j,i)*abs(a(k,j,i))**2
                  END DO
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
      END SUBROUTINE variance

!*****************************************************************
      SUBROUTINE product(a,b,c)
!-----------------------------------------------------------------
!
! Computes the integral of the product of two scalars. 
! The output is only valid in the first node.
!
! Parameters
!     a  : first scalar
!     b  : second scalar
!     c  : at the output contains the product
!
      USE fprecision
      USE commtypes
      USE grid
      USE mpivars
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a,b
      DOUBLE PRECISION, INTENT(OUT) :: c
      DOUBLE PRECISION              :: cloc
      REAL(KIND=GP)                 :: tmp
      INTEGER             :: i,j,k

      cloc = 0.0D0
      tmp = 1.0_GP/ &
            (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2
!
! Computes the averaged inner product between the fields
!
      IF (ista.eq.1) THEN
!$omp parallel private (k) reduction(+:cloc)
!$omp do
         DO j = 1,ny
            DO CONCURRENT (k=1:nz)
               cloc = cloc+tmp*real(a(k,j,1)*conjg(b(k,j,1)))
            END DO
         END DO
!$omp end do
!$omp do collapse(2)
         DO i = 2,iend
            DO j = 1,ny
               DO CONCURRENT (k=1:nz)
                  cloc = cloc+2*tmp*real(a(k,j,i)*conjg(b(k,j,i)))
               END DO
            END DO
         END DO
!$omp end do
!$omp end parallel
      ELSE
!$omp parallel do collapse(2) private (k) reduction(+:cloc)
         DO i = ista,iend
            DO j = 1,ny
               DO CONCURRENT (k=1:nz)
                  cloc = cloc+2*tmp*real(a(k,j,i)*conjg(b(k,j,i)))
               END DO
            END DO
         END DO
      ENDIF
!
! Computes the reduction between nodes
!
      CALL MPI_REDUCE(cloc,c,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                         MPI_COMM_WORLD,ierr)

      RETURN
      END SUBROUTINE product

!*****************************************************************
      SUBROUTINE spectrsc(a,path,nmb,isc,tail)
!-----------------------------------------------------------------
!
! Computes the passive/active scalar power spectrum.
! Normalization of the spectrum is such that E = sum[E(k).Dkk],
! where Dkk is the width of the Fourier shells. The output
! is written to a file by the first node.
!
! Output files contain:
! 'sspectrum.XXX.txt' : k, V(k) (power spectrum of the scalar)
! 'sNspectrum.XXX.txt': k, V(k) (same for the N-th scalar)
!
! Parameters
!     a    : input matrix with the scalar
!     path : path for the output
!     nmb  : the extension used when writting the file
!     isc  : index to specify which scalar the spectrum 
!            represents; modifies output file name. If 
!            isc < 0, then we get the following prefixes:
!            -1 ==> 'rhospect.XXX.txt'
!     tail : Appends tail at the end of the file name [optional]
!
      USE kes
      USE grid
      USE mpivars
      USE filefmt
      USE boxsize
      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(nmax/2+1)                    :: Ek
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a
      INTEGER,          INTENT(IN)                             :: isc
      INTEGER                                :: i
      CHARACTER(len=*), INTENT(IN)           :: path,nmb
      CHARACTER(len=*), INTENT(IN), OPTIONAL :: tail
      CHARACTER(len=1)                       :: si
      CHARACTER(len=128)                     :: fname

!
! Computes the power spectrum
!
      CALL spectrscc(a,Ek,0.0_GP)
!
! Exports the spectrum to a file
!
      IF ( myrank.eq.0 ) THEN
         IF ( isc.ge.0 ) THEN
           IF ( isc.gt.0 ) THEN
              WRITE(si,'(i1.1)') isc
              fname = 's' // si // 'spectrum'
           ELSE
              fname = 'sspectrum'
           ENDIF
         ELSE IF ( isc .eq. -1 ) THEN
           fname = 'rhospectrum'
         ENDIF
         if (present(tail)) then
            fname = trim(adjustl(fname)) // '_' // trim(adjustl(tail))
         endif
         OPEN(1,file= trim(path) // '/' // trim(adjustl(fname)) // '.' &
              // nmb // '.txt')
         DO i=1,nmax/2+1
            WRITE(1,FMT='(E13.6,E23.15)')  Dkk*i,Ek(i)/Dkk
         END DO
         CLOSE(1)
      ENDIF

      RETURN
      END SUBROUTINE spectrsc

!*****************************************************************
      SUBROUTINE spectrscc(a,Ektot,shift)
!-----------------------------------------------------------------
!
! Computes the passive/active scalar power spectrum, returning it.
!
! Parameters
!     a    : input matrix with the scalar
!     Ektot: output power spectrum
!     shift: value that can be used to shift wavenumbers 
!            (usually by 1) and get the spetrum to start at k=0 
!
      USE fprecision
      USE commtypes
      USE kes
      USE grid
      USE mpivars
      USE boxsize
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT (IN), DIMENSION(nz,ny,ista:iend) :: a
      DOUBLE PRECISION, INTENT(OUT), DIMENSION(nmax/2+1)    :: Ektot
      DOUBLE PRECISION,              DIMENSION(nmax/2+1)    :: Ek
      DOUBLE PRECISION :: tmq
      REAL(KIND=GP),    INTENT(IN)                          :: shift
      REAL(KIND=GP)    :: tmp,round
      INTEGER          :: i,j,k
      INTEGER          :: kmn

!
! Sets Ek to zero
!
      DO i = 1,nmax/2+1
         Ek(i) = 0.0D0
      END DO
!
! Sets the zero for the wavenumbers
!
      round = shift+.501_GP
!
! Computes the power spectrum
!
      tmp = 1.0_GP/ &
            (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2
      IF (ista.eq.1) THEN
!$omp parallel private (k,kmn,tmq) reduction(+:Ek)
!$omp do
         DO j = 1,ny
            DO k = 1,nz
               kmn = int(sqrt(kk2(k,j,1))/Dkk+round)
               IF ((kmn.gt.0).and.(kmn.le.nmax/2+1)) THEN
                  tmq = tmp*abs(a(k,j,1))**2
                  Ek(kmn) = Ek(kmn)+tmq
               ENDIF
            END DO
         END DO
!$omp end do
!$omp do collapse(2)
         DO i = 2,iend
            DO j = 1,ny
               DO k = 1,nz
                  kmn = int(sqrt(kk2(k,j,i))/Dkk+round)
                  IF ((kmn.gt.0).and.(kmn.le.nmax/2+1)) THEN
                     tmq = 2*tmp*abs(a(k,j,i))**2
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
                  kmn = int(sqrt(kk2(k,j,i))/Dkk+round)
                  IF ((kmn.gt.0).and.(kmn.le.nmax/2+1)) THEN
                     tmq = 2*tmp*abs(a(k,j,i))**2
                     Ek(kmn) = Ek(kmn)+tmq
                  ENDIF
               END DO
            END DO
         END DO
      ENDIF
!
! Computes the reduction between nodes
!
      CALL MPI_ALLREDUCE(Ek,Ektot,nmax/2+1,MPI_DOUBLE_PRECISION, &
                      MPI_SUM,MPI_COMM_WORLD,ierr)

      RETURN
      END SUBROUTINE spectrscc

END MODULE pseudospec_scalar

!=================================================================

MODULE pseudospec_phd
   USE pseudospec_scalar
   CONTAINS

!*****************************************************************
      SUBROUTINE pscheck(a,b,t,dt,path,ext)
!-----------------------------------------------------------------
!
! Writes a global file with quantities relevant for a passive
! or active scalar.
!
! Output file contains:
! 'scalar.txt':  time, <theta^2>, <|grad(theta)|^2>, injection rate
!
! Parameters
!     a   : scalar concentration
!     b   : source of the scalar
!     t   : number of time steps made
!     dt  : time step
!     path: path for the output
!     ext : file extension [optional]
!
      USE fprecision
      USE grid
      USE mpivars
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a,b
      DOUBLE PRECISION          :: eng,ens,pot
      REAL(KIND=GP), INTENT(IN) :: dt
      INTEGER, INTENT(IN)       :: t
      INTEGER                   :: i,j,k
      CHARACTER(len=*), INTENT(IN)           :: path
      CHARACTER(len=*), INTENT(IN), OPTIONAL :: ext
      CHARACTER(len=128)        :: fname

!
! Computes the variance and the variance of k^2 times the scalar
!
      CALL variance(a,eng,1)
      CALL variance(a,ens,0)
!
! Computes the scalar injection rate
!
      CALL product(a,b,pot)
!
! Creates external files to store the results
!
      IF (myrank.eq.0) THEN
         if (present(ext)) then
            fname = 'scalar_' // trim(adjustl(ext)) // '.txt'
         else
            fname = 'scalar.txt'
         endif
         OPEN(1,file=trim(path) // '/' // trim(adjustl(fname)), &
              position='append')
         WRITE(1,10) (t-1)*dt,eng,ens,pot
   10    FORMAT( E13.6,E22.14,E22.14,E22.14 )
         CLOSE(1)
      ENDIF

      RETURN
    END SUBROUTINE pscheck

END MODULE pseudospec_phd
  
