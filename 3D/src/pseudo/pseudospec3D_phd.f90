!=================================================================
! PSEUDOSPECTRAL subroutines
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
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT (IN), DIMENSION(nz,ny,ista:iend) :: a,b
      COMPLEX(KIND=GP), INTENT (IN), DIMENSION(nz,ny,ista:iend) :: c,d
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(nz,ny,ista:iend) :: e
      COMPLEX(KIND=GP), DIMENSION(nz,ny,ista:iend) :: c1,c2
      REAL(KIND=GP),    DIMENSION(nx,ny,ksta:kend) :: r1,r2
      REAL(KIND=GP),    DIMENSION(nx,ny,ksta:kend) :: r3
      REAL(KIND=GP)    :: tmp
      INTEGER :: i,j,k

!
! Computes (A_x.dx)B
!
      c1 = a
      CALL derivk3(d,c2,1)
      CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,c2,r2,MPI_COMM_WORLD)

!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
         DO j = 1,ny
            DO i = 1,nx
               r3(i,j,k) = r1(i,j,k)*r2(i,j,k)
            END DO
         END DO
      END DO
!
! Computes (A_y.dy)B
!
      c1 = b
      CALL derivk3(d,c2,2)
      CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,c2,r2,MPI_COMM_WORLD)

!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
         DO j = 1,ny
            DO i = 1,nx
               r3(i,j,k) = r3(i,j,k)+r1(i,j,k)*r2(i,j,k)
            END DO
         END DO
      END DO
!
! Computes (A_z.dz)B
!
      c1 = c
      CALL derivk3(d,c2,3)
      CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,c2,r2,MPI_COMM_WORLD)

! We need -A.grad(B)
      tmp = -1.0_GP/ &
            (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
         DO j = 1,ny
            DO i = 1,nx
               r3(i,j,k) = (r3(i,j,k)+r1(i,j,k)*r2(i,j,k))*tmp
            END DO
         END DO
      END DO

      CALL fftp3d_real_to_complex(planrc,r3,e,MPI_COMM_WORLD)

      RETURN
      END SUBROUTINE advect3

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
!$    USE threads
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
!$omp parallel do private (k) reduction(+:bloc)
            DO j = 1,ny
               DO k = 1,nz
                  bloc = bloc+tmp*abs(a(k,j,1))**2
               END DO
            END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k) reduction(+:bloc)
            DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k) reduction(+:bloc)
               DO j = 1,ny
                  DO k = 1,nz
                     bloc = bloc+2*tmp*abs(a(k,j,i))**2
                  END DO
               END DO
            END DO
         ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k) reduction(+:bloc)
            DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k) reduction(+:bloc)
               DO j = 1,ny
                  DO k = 1,nz
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
!$omp parallel do private (k) reduction(+:bloc)
            DO j = 1,ny
               DO k = 1,nz
                  bloc = bloc+tmp*kk2(k,j,1)*abs(a(k,j,1))**2
               END DO
            END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k) reduction(+:bloc)
            DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k) reduction(+:bloc)
               DO j = 1,ny
                  DO k = 1,nz
                     bloc = bloc+2*tmp*kk2(k,j,i)*abs(a(k,j,i))**2
                  END DO
               END DO
            END DO
         ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k) reduction(+:bloc)
            DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k) reduction(+:bloc)
               DO j = 1,ny
                  DO k = 1,nz
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
!$    USE threads
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
!$omp parallel do private (k) reduction(+:cloc)
         DO j = 1,ny
            DO k = 1,nz
               cloc = cloc+tmp*real(a(k,j,1)*conjg(b(k,j,1)))
            END DO
         END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k) reduction(+:cloc)
         DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k) reduction(+:cloc)
            DO j = 1,ny
               DO k = 1,nz
                  cloc = cloc+2*tmp*real(a(k,j,i)*conjg(b(k,j,i)))
               END DO
            END DO
         END DO
      ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k) reduction(+:cloc)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k) reduction(+:cloc)
            DO j = 1,ny
               DO k = 1,nz
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
      SUBROUTINE pscheck(a,b,t,dt)
!-----------------------------------------------------------------
!
! Consistency check for the conservation of energy, 
! helicity, and null divergency of the velocity field
!
! Output file contains:
! 'scalar.txt':  time, <theta^2>, <|grad(theta)|^2>, injection rate
!
! Parameters
!     a : scalar concentration
!     b : source of the scalar
!     t : number of time steps made
!     dt: time step
!
      USE fprecision
      USE grid
      USE mpivars
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a,b
      DOUBLE PRECISION    :: eng,ens,pot
      REAL(KIND=GP), INTENT(IN)    :: dt
      INTEGER, INTENT(IN) :: t
      INTEGER             :: i,j,k

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
         OPEN(1,file='scalar.txt',position='append')
         WRITE(1,10) (t-1)*dt,eng,ens,pot
   10    FORMAT( E13.6,E22.14,E22.14,E22.14 )
         CLOSE(1)
      ENDIF

      RETURN
      END SUBROUTINE pscheck

!*****************************************************************
      SUBROUTINE mpscheck2(a1,b1,a2,b2,t,dt)
!-----------------------------------------------------------------
!
! Consistency check for the conservation of 2 'multiscalar' energies
!
! Output file contains:
! 'mscalar.txt':  time,
!     [FIRST SCALAR:]  <theta^2>, <|grad(theta)|^2>, injection rate,
!     [SECOND SCALAR:] <theta^2>, <|grad(theta)|^2>, injection rate
!
! Parameters
!     a_i : i_th scalar concentration, i=1-2
!     b_i : i_th source/force of the scalar i
!     t : number of time steps made
!     dt: time step
!
      USE fprecision
      USE grid
      USE mpivars
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a1,a2
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: b1,b2
      DOUBLE PRECISION    :: eng(2),ens(2),pot(2)
      REAL(KIND=GP), INTENT(IN)    :: dt
      INTEGER, INTENT(IN) :: t
      INTEGER             :: i,j,k

!
! Computes the variance and the variance of k^2 times the scalar
!
      CALL variance(a1,eng(1),1)
      CALL variance(a1,ens(1),0)
      CALL variance(a2,eng(2),1)
      CALL variance(a2,ens(2),0)
!
! Computes the scalar injection rate
!
      CALL product(a1,b1,pot(1))
      CALL product(a2,b2,pot(2))
!
! Creates external files to store the results
!
      IF (myrank.eq.0) THEN
         OPEN(1,file='mscalar.txt',position='append')
         WRITE(1,20) (t-1)*dt,eng(1),ens(1),pot(1),eng(2),ens(2),pot(2)
   20    FORMAT( E13.6,E22.14,E22.14,E22.14,E22.14,E22.14,E22.14 )
         CLOSE(1)
      ENDIF

      RETURN
      END SUBROUTINE mpscheck2

!*****************************************************************
      SUBROUTINE mpscheck3(a1,b1,a2,b2,a3,b3,t,dt)
!-----------------------------------------------------------------
!
! Consistency check for the conservation of 3 'multiscalar' energies
!
! Output file contains:
! 'mscalar.txt':  time,
!     [FIRST SCALAR:]  <theta^2>, <|grad(theta)|^2>, injection rate,
!     [SECOND SCALAR:] <theta^2>, <|grad(theta)|^2>, injection rate,
!     [THIRD SCALAR:]  <theta^2>, <|grad(theta)|^2>, injection rate
!
! Parameters
!     a_i : i_th scalar concentration, i=1-3
!     b_i : i_th source/force of the scalar i
!     t : number of time steps made
!     dt: time step
!
      USE fprecision
      USE grid
      USE mpivars
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a1,a2
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a3,b1
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: b2,b3
      DOUBLE PRECISION    :: eng(3),ens(3),pot(3)
      REAL(KIND=GP), INTENT(IN)    :: dt
      INTEGER, INTENT(IN) :: t
      INTEGER             :: i,j,k

!
! Computes the variance and the variance of k^2 times the scalar
!
      CALL variance(a1,eng(1),1)
      CALL variance(a1,ens(1),0)
      CALL variance(a2,eng(2),1)
      CALL variance(a2,ens(2),0)
      CALL variance(a3,eng(3),1)
      CALL variance(a3,ens(3),0)
!
! Computes the scalar injection rate
!
      CALL product(a1,b1,pot(1))
      CALL product(a2,b2,pot(2))
      CALL product(a3,b3,pot(3))
!
! Creates external files to store the results
!
      IF (myrank.eq.0) THEN
         OPEN(1,file='mscalar.txt',position='append')
         WRITE(1,30) (t-1)*dt,eng(1),ens(1),pot(1),eng(2),ens(2), &
              pot(2),eng(3),ens(3),pot(3)
   30    FORMAT( E13.6,E22.14,E22.14,E22.14,E22.14,E22.14,        &
              E22.14,E22.14,E22.14,E22.14 )
         CLOSE(1)
      ENDIF

      RETURN
      END SUBROUTINE mpscheck3

!*****************************************************************
      SUBROUTINE spectrsc(a,nmb,isc)
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
!     nmb  : the extension used when writting the file
!     isc  : index to specify which scalar the spectrum 
!            represents; modifies output file name. If 
!            isc < 0, then we get the following prefixes:
!            -1 ==> 'rhospect.XXX.txt'
!
      USE kes
      USE grid
      USE mpivars
      USE filefmt
      USE boxsize
!$    USE threads
      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(nmax/2+1)                   :: Ek
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a
      INTEGER, INTENT(IN)                                    :: isc
      INTEGER                                :: i
      CHARACTER(len=*),           INTENT(IN) :: nmb
      CHARACTER(len=1)                       :: si

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
             OPEN(1,file='s' // si // 'spectrum.' // nmb // '.txt')
           ELSE
             OPEN(1,file='sspectrum.' // nmb // '.txt')
           ENDIF
         ELSE IF ( isc .eq. -1 ) THEN
             OPEN(1,file='rhospectrum.' // nmb // '.txt')
         ELSE
             PRINT*, 'SPECTRSC: Invalid option'
             STOP
         ENDIF
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
!$    USE threads
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
!$omp parallel do private (k,kmn,tmq)
         DO j = 1,ny
            DO k = 1,nz
               kmn = int(sqrt(kk2(k,j,1))/Dkk+round)
               IF ((kmn.gt.0).and.(kmn.le.nmax/2+1)) THEN
                  tmq = tmp*abs(a(k,j,1))**2
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
                  kmn = int(sqrt(kk2(k,j,i))/Dkk+round)
                  IF ((kmn.gt.0).and.(kmn.le.nmax/2+1)) THEN
                     tmq = 2*tmp*abs(a(k,j,i))**2
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
                  kmn = int(sqrt(kk2(k,j,i))/Dkk+round)
                  IF ((kmn.gt.0).and.(kmn.le.nmax/2+1)) THEN
                     tmq = 2*tmp*abs(a(k,j,i))**2
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
      CALL MPI_ALLREDUCE(Ek,Ektot,nmax/2+1,MPI_DOUBLE_PRECISION, &
                      MPI_SUM,MPI_COMM_WORLD,ierr)

      RETURN
      END SUBROUTINE spectrscc

!*****************************************************************
      SUBROUTINE sctrans(a,b,nmb,isc)
!-----------------------------------------------------------------
!
! Computes the scalar transfer in Fourier space in 3D.
! Normalization of the transfer function is such that the
! flux is Pi = -sum[T(k).Dkk], where Dkk is the width of the
! Fourier shells. The output is written to a file by the 
! first node.
!
! Output files contain:
! 'stransfer.XXX.txt' : k, Ts(k) (scalar transfer function)
! 'sNtransfer.XXX.txt': k, Ts(k) (same for the N-th scalar)
!
! Parameters
!     a  : scalar
!     b  : nonlinear term
!     nmb: the extension used when writting the file
!     isc: if doing multi-scalar, gives index of scalar 
!          whose transfer is being computed (1, 2, or 3) and
!          names file as s<isc>transfer.XXX.txt. If isc=0, then
!          filename is stransfer.XXX.txt
!          
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

      COMPLEX (KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a,b
      DOUBLE PRECISION, DIMENSION(nmax/2+1) :: Ek,Ektot
      DOUBLE PRECISION :: tmq
      REAL(KIND=GP)    :: tmp
      INTEGER,           INTENT(IN)                             :: isc
      INTEGER          :: i,j,k
      INTEGER          :: kmn
      CHARACTER(len=*) , INTENT(IN)                             :: nmb
      CHARACTER(len=1) :: si
!
! Sets Ek to zero
!
      DO i = 1,nmax/2+1
         Ek(i) = 0.0D0
      END DO
!
! Computes the scalar transfer
!
      tmp = 1.0_GP/ &
            (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2
      IF (ista.eq.1) THEN
!$omp parallel do private (k,kmn,tmq)
         DO j = 1,ny
            DO k = 1,nz
               kmn = int(sqrt(kk2(k,j,1))/Dkk+.501)
               IF ((kmn.gt.0).and.(kmn.le.nmax/2+1)) THEN
                  tmq = tmp*real(a(k,j,1)*conjg(b(k,j,1)))
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
                  kmn = int(sqrt(kk2(k,j,i))/Dkk+.501)
                  IF ((kmn.gt.0).and.(kmn.le.nmax/2+1)) THEN
                     tmq = 2*tmp*real(a(k,j,i)*conjg(b(k,j,i)))
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
                  kmn = int(sqrt(kk2(k,j,i))/Dkk+.501)
                  IF ((kmn.gt.0).and.(kmn.le.nmax/2+1)) THEN
                     tmq = 2*tmp*real(a(k,j,i)*conjg(b(k,j,i)))
!$omp atomic
                     Ek(kmn) = Ek(kmn)+tmq
                  ENDIF
               END DO
            END DO
         END DO
      ENDIF
!
! Computes the reduction between nodes
! and exports the result to a file
!
      CALL MPI_REDUCE(Ek,Ektot,nmax/2+1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
      IF (myrank.eq.0) THEN
        IF ( isc.GT.0 ) THEN
          WRITE(si,'(i1.1)') isc
          OPEN(1,file='s' // trim(si) // 'transfer.' // nmb // '.txt')
        ELSE
          OPEN(1,file='stransfer.' // nmb // '.txt')
        ENDIF
        DO i=1,nmax/2+1
           WRITE(1,FMT='(E13.6,E23.15)')  Dkk*i,Ektot(i)/Dkk
        END DO
      ENDIF

      RETURN
      END SUBROUTINE sctrans

!*****************************************************************
      SUBROUTINE difucx(a,b,nmb)
!-----------------------------------------------------------------
!
! Computes the mean profiles in x of the velocity, the 
! passive scalar, and their product. The output is 
! written to a file by the first node.
!
! Output file contains:
! 'profilex.txt': x, <v_i>(x), <theta>(x), <v_i.theta>(x)
!
! Parameters
!     a    : vector field component in the x-direction
!     b    : scalar field
!     nmb: the extension used when writting the file
!
      USE fprecision
      USE commtypes
      USE var
      USE kes
      USE fft
      USE grid
      USE mpivars
      USE boxsize
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a,b
      COMPLEX(KIND=GP),             DIMENSION(nz,ny,ista:iend) :: c1,c2
      REAL(KIND=GP),                DIMENSION(nx,ny,ksta:kend) :: r1,r2
      REAL(KIND=GP), DIMENSION(nx) :: meth,mev,methv
      REAL(KIND=GP), DIMENSION(nx) :: mth,mv,mthv
      REAL(KIND=GP)                :: tmp,tmq
      INTEGER                      :: i,j,k
      CHARACTER(len=*), INTENT(IN) :: nmb

!
! Transforms the input arrays to real space
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
! Computes the mean profiles
!    
      DO i = 1,nx
         mev(i) = 0.0_GP
         meth(i) = 0.0_GP
         methv(i) = 0.0_GP
      END DO
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
         DO j = 1,ny
            DO i = 1,nx
!$omp critical
               mev(i) = mev(i)+r1(i,j,k)
               meth(i) = meth(i)+r2(i,j,k)
               methv(i) = methv(i)+r1(i,j,k)*r2(i,j,k)
!$omp end critical
            END DO
         END DO
      END DO
      tmp = 1.0_GP/(real(nx,kind=GP)*      &
                    real(ny,kind=GP)**2*real(nz,kind=GP)**2)
      tmq = 1.0_GP/(real(nx,kind=GP)**2*   &
                    real(ny,kind=GP)**3*real(nz,kind=GP)**3)
      DO i = 1,nx
         mev(i) = mev(i)*tmp
         meth(i) = meth(i)*tmp
         methv(i) = methv(i)*tmq
      END DO
!
! Computes the reduction between nodes
!
      CALL MPI_REDUCE(mev,mv,nx,GC_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(meth,mth,nx,GC_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(methv,mthv,nx,GC_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)

      IF (myrank.eq.0) THEN
         OPEN(1,file='profilex.' // nmb // '.txt')
         DO i = 1,nx
            WRITE(1,40) 2*pi*Lx*(real(i,kind=GP)-1)/real(nx,kind=GP), &
                        mv(i),mth(i),mthv(i)
         END DO
         CLOSE(1) 
   40    FORMAT( E23.15,E23.15,E23.15,E23.15 )
      ENDIF

      RETURN
      END SUBROUTINE difucx

!*****************************************************************
      SUBROUTINE difucz(a,b,nmb)
!-----------------------------------------------------------------
!
! Computes the mean profiles in z of the velocity, the 
! passive scalar, and their product. The output is 
! written to a file by the first node.
!
! Output file contains:
! 'profilez.txt': z, <v_i>(z), <theta>(z), <v_i.theta>(z)
!
! Parameters
!     a    : vector field component in the z-direction
!     b    : scalar field
!     nmb  : number of blocks
!
      USE fprecision
      USE commtypes
      USE var
      USE kes
      USE fft
      USE grid
      USE mpivars
      USE boxsize
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a,b
      COMPLEX(KIND=GP),             DIMENSION(nz,ny,ista:iend) :: c1,c2
      REAL(KIND=GP),                DIMENSION(nx,ny,ksta:kend) :: r1,r2
      REAL(KIND=GP), DIMENSION(nz) :: meth,mev,methv
      REAL(KIND=GP), DIMENSION(nz) :: mth,mv,mthv
      REAL(KIND=GP)                :: tmp,tmq
      INTEGER                      :: i,j,k
      CHARACTER(len=*), INTENT(IN) :: nmb

!
! Transforms the input arrays to real space
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
! Computes the mean profiles
!    
      DO i = 1,nz
         mev(i) = 0.0_GP
         meth(i) = 0.0_GP
         methv(i) = 0.0_GP
      END DO
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
         DO j = 1,ny
            DO i = 1,nx
!$omp critical
               mev(k) =  mev(k)+r1(i,j,k)
               meth(k) = meth(k)+r2(i,j,k)
               methv(k) = methv(k)+r1(i,j,k)*r2(i,j,k)  
!$omp end critical
            END DO
         END DO
      END DO
      tmp = 1.0_GP/(real(nx,kind=GP)**2*   &
                    real(ny,kind=GP)**2*real(nz,kind=GP))
      tmq = 1.0_GP/(real(nx,kind=GP)**3*   &
                    real(ny,kind=GP)**3*real(nz,kind=GP)**2)
      DO k = 1,nz
         mev(k) = mev(k)*tmp
         meth(k) = meth(k)*tmp
         methv(k) = methv(k)*tmq
      END DO
!
! Computes the reduction between nodes
!
      CALL MPI_REDUCE(mev,mv,nz,GC_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(meth,mth,nz,GC_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(methv,mthv,nz,GC_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)

      IF (myrank.eq.0) THEN
         DO k = 1,nz
            WRITE(1,50) 2*pi*Lz*(real(k,kind=GP)-1)/real(nz,kind=GP), &
                        mv(k),mth(k),mthv(k)
         END DO
         CLOSE(1) 
   50    FORMAT( E23.15,E23.15,E23.15,E23.15 ) 
      ENDIF

      RETURN
      END SUBROUTINE difucz
