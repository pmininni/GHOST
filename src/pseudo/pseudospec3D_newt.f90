!=================================================================
! PSEUDOSPECTRAL subroutines
!
! Subroutines to do Newton and stabilized biconjugate gradient
! methods to prepare data for the GPE solver. Can be easily
! modified to prepare initial conditions close to a fixed
! point for any solver. You should use the FFTPLANS and
! MPIVARS modules (see the file 'fftp_mod.f90') in each 
! program that calls any of the subroutines in this file. 
!
! NOTATION: index 'i' is 'x' 
!           index 'j' is 'y'
!           index 'k' is 'z'
!
! 2015 P.D. Mininni and M.E. Brachet
!
! 17 May 2018: Support for elongated box (N.Muller & P.D.Mininni) 
!=================================================================

!*****************************************************************
      SUBROUTINE ThreeTo1D(zre,zim,xvec1d)
!-----------------------------------------------------------------
!
! Copies the data into a 1D real array to do the Newton method.
!
! Parameters
!     zre     : Real part of the wavefunction
!     zim     : Imaginary part of the wavefunction
!     xvec1d  : 1D vector
!
      USE fprecision
      USE mpivars
      USE newtmod
      USE grid
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: zre,zim
      REAL(KIND=GP), INTENT(OUT), DIMENSION(n_dim_1d)          :: xvec1D
      INTEGER             :: i,j,k
      INTEGER             :: offset1,offset2

!$omp parallel do if ((iend-ista).ge.nth) private(j,k,offset1,offset2)
      DO i = ista,iend
         offset1 = 4*(i-ista)*ny*nz
!$omp parallel do if ((iend-ista).lt.nth) private(k,offset2)
         DO j = 1,ny
            offset2 = offset1 + 4*(j-1)*nz
            DO k = 1,nz
               xvec1D(1+4*(k-1)+offset2) = REAL(zre(k,j,i))
               xvec1D(2+4*(k-1)+offset2) = AIMAG(zre(k,j,i))
               xvec1D(3+4*(k-1)+offset2) = REAL(zim(k,j,i))
               xvec1D(4+4*(k-1)+offset2) = AIMAG(zim(k,j,i))
            END DO
         END DO
      END DO

      RETURN
      END SUBROUTINE ThreeTo1D

!*****************************************************************
      SUBROUTINE OneTo3D(xvec1d,zre,zim)
!-----------------------------------------------------------------
!
! Copies the data back into 3D complex arrays, after doing the 
! Newton method.
!
! Parameters
!     xvec1d  : 1D vector
!     zre     : Real part of the wavefunction
!     rim     : Imaginary part of the wavefunction
!
      USE fprecision
      USE mpivars
      USE newtmod
      USE grid
      USE var
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(nz,ny,ista:iend) :: zre,zim
      REAL(KIND=GP), INTENT(IN), DIMENSION(n_dim_1d)            :: xvec1D
      INTEGER             :: i,j,k
      INTEGER             :: offset1,offset2

!$omp parallel do if ((iend-ista).ge.nth) private(j,k,offset1,offset2)
      DO i = ista,iend
         offset1 = 4*(i-ista)*ny*nz
!$omp parallel do if ((iend-ista).lt.nth) private(k,offset2)
         DO j = 1,ny
            offset2 = offset1+4*(j-1)*nz
            DO k = 1,nz
               zre(k,j,i) = xvec1D(1+4*(k-1)+offset2)+ &
                         im*xvec1D(2+4*(k-1)+offset2)
               zim(k,j,i) = xvec1D(3+4*(k-1)+offset2)+ &
                         im*xvec1D(4+4*(k-1)+offset2)
            END DO
         END DO
      END DO

      RETURN
      END SUBROUTINE OneTo3D

!*****************************************************************
      SUBROUTINE gtrhs(xguess,rhs,cflow,dt,vx,vy,vz,vsq,R1,C3,C4,C5,C6)
!-----------------------------------------------------------------
!
! Computes -rhs of the ARGL equation.
!
! Parameters
!     xguess      : 1D vector state
!     rhs         : rhs of the ARGL equation
!     cflow       : =1 if mean constant velocity is present
!     dt          : time step
!     vx,vy,vz,vsq: arrays with the velocity field
!     R1,C1-C6    : auxiliary arrays
!
      USE fprecision
      USE mpivars
      USE newtmod
      USE grid
      USE hbar
      USE kes
      USE ali
!$    USE threads
      IMPLICIT NONE
      
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend)    :: vx,vy,vz
      COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(nz,ny,ista:iend) :: C3,C4
      COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(nz,ny,ista:iend) :: C5,C6
      COMPLEX(KIND=GP), DIMENSION(nz,ny,ista:iend)              :: zre,zim
      REAL(KIND=GP), INTENT(IN), DIMENSION(nx,ny,ksta:kend)     :: vsq
      REAL(KIND=GP), INTENT(INOUT), DIMENSION(nx,ny,ksta:kend)  :: R1
      REAL(KIND=GP), INTENT(IN), DIMENSION(n_dim_1d)  :: xguess
      REAL(KIND=GP), INTENT(OUT), DIMENSION(n_dim_1d) :: rhs
      REAL(KIND=GP), INTENT(IN) :: dt
      REAL(KIND=GP)       :: rmp,rmq
      INTEGER, INTENT(IN) :: cflow
      INTEGER             :: i,j,k

      CALL OneTo3D(xguess,zre,zim)
      
! Here we basically have the contents of "include/argl_rkstep2.f90".
! We use Newton when we don't have counterflow, so we have no forcing
! function. However we still allow for counterflow.neq.0 for the case 
! in which we have a mean constant velocity (e.g., to advect a ring).
! The equations are:
!    rhs.re = dt.(omegag.zre - beta.|z|^2 zre - |v^2|.zre/(4.alpha) +
!             + v.grad(zim) - alpha.k^2.zre)/(1+alpha.k^2.dt)
!    rhs.im = dt.(omegag.zim - beta.|z|^2 zim - |v^2|.zim/(4.alpha) -
!             - v.grad(zre) - alpha.k^2.zim)/(1+alpha.k^2.dt)

      CALL squareabs(zre,zim,R1,1)
      rmq = real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP)*omegag/beta
      IF (cflow.eq.0) THEN ! If not doing counterflow we have |v^2|
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
         DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
            DO j = 1,ny
               DO i = 1,nx
                  R1(i,j,k) = rmq-R1(i,j,k)-vsq(i,j,k)
               END DO
            END DO
         END DO
      ELSE
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
         DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
            DO j = 1,ny
               DO i = 1,nx
                  R1(i,j,k) = rmq-R1(i,j,k)
               END DO
            END DO
         END DO
      ENDIF
      CALL nonlgpe(R1,zre,C3)
      CALL nonlgpe(R1,zim,C4)
      CALL advect3(vx,vy,vz,zre,C5) ! -v.grad(zre)
      CALL advect3(vx,vy,vz,zim,C6) ! -v.grad(zim)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
            DO k = 1,nz
               IF (kn2(k,j,i).le.kmax) THEN
                  rmp = dt/(1.0_GP+alpha*kk2(k,j,i)*dt)
                  zre(k,j,i) = -(beta*C3(k,j,i)-C6(k,j,i)-    &
                               alpha*kk2(k,j,i)*zre(k,j,i))*rmp
                  zim(k,j,i) = -(beta*C4(k,j,i)+C5(k,j,i)-    &
                               alpha*kk2(k,j,i)*zim(k,j,i))*rmp
               ELSE
                  zre(k,j,i) = 0.0_GP
                  zim(k,j,i) = 0.0_GP
               ENDIF
            END DO
         END DO
      END DO
      
      CALL ThreeTo1D(zre,zim,rhs)
      RETURN
      END SUBROUTINE gtrhs

!*****************************************************************
      SUBROUTINE lin(xguess,from,to,cflow,dt,vx,vy,vz,vsq,R1,C3,C4,C5,C6)
!-----------------------------------------------------------------
!
! Computes the linearized rhs of the ARGL equation:
! lin = dt.(omega.dz -|v|^2.dz/(4.alpha) - i v.grad(dz) - beta. 
!       (2 |z|^2.dz + z^2.dz*) - alpha k^2 dz)/(1 + alpha k^2 dt)
!
! Parameters
!     xguess      : 1D vector state
!     from        : perturbation to the 1D vector state
!     to          : linearized rhs
!     cflow       : =1 if mean constant velocity is present
!     dt          : time step
!     vx,vy,vz,vsq: arrays with the velocity field
!     R1,C1-C6    : auxiliary arrays
!
      USE fprecision
      USE commtypes
      USE mpivars
      USE newtmod
      USE grid
      USE hbar
      USE kes
      USE ali
      USE fft
!$    USE threads
      IMPLICIT NONE
      
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend)    :: vx,vy,vz
      COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(nz,ny,ista:iend) :: C3,C4
      COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(nz,ny,ista:iend) :: C5,C6
      COMPLEX(KIND=GP), DIMENSION(nz,ny,ista:iend) :: C1,C2
      COMPLEX(KIND=GP), DIMENSION(nz,ny,ista:iend)             :: zre,zim
      COMPLEX(KIND=GP), DIMENSION(nz,ny,ista:iend)             :: zreg,zimg
      REAL(KIND=GP), INTENT(IN), DIMENSION(nx,ny,ksta:kend)    :: vsq
      REAL(KIND=GP), INTENT(INOUT), DIMENSION(nx,ny,ksta:kend) :: R1
      REAL(KIND=GP), DIMENSION(nx,ny,ksta:kend)       :: R2
      REAL(KIND=GP), INTENT(IN), DIMENSION(n_dim_1d)  :: xguess,from
      REAL(KIND=GP), INTENT(OUT), DIMENSION(n_dim_1d) :: to
      REAL(KIND=GP), INTENT(IN) :: dt
      REAL(KIND=GP)       :: rmp,rmq
      INTEGER, INTENT(IN) :: cflow
      INTEGER             :: i,j,k

      CALL OneTo3D(xguess,zre,zim)
      CALL OneTo3D(from,zreg,zimg)

      CALL squareabs(zre,zim,R1,1)
      rmq = real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP)*omegag/beta
      IF (cflow.eq.0) THEN ! If not doing counterflow we have |v^2|
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
         DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
            DO j = 1,ny
               DO i = 1,nx
                  R1(i,j,k) = rmq-2*R1(i,j,k)-vsq(i,j,k)
               END DO
            END DO
         END DO
      ELSE
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
         DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
            DO j = 1,ny
               DO i = 1,nx
                  R1(i,j,k) = rmq-2*R1(i,j,k)
               END DO
            END DO
         END DO
      ENDIF
      CALL nonlgpe(R1,zreg,C3)
      CALL nonlgpe(R1,zimg,C4)
      CALL advect3(vx,vy,vz,zreg,C5)        ! -v.grad(zreg)
      CALL advect3(vx,vy,vz,zimg,C6)        ! -v.grad(zimg)
      CALL zsquare(zre,zim,R1,R2,1)
      CALL pertgpe(R1,R2,zreg,zimg,zre,zim) ! (z^2).zreg*
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
            DO k = 1,nz
               IF (kn2(k,j,i).le.kmax) THEN
                  rmp = dt/(1.0_GP+alpha*kk2(k,j,i)*dt)
                  zre(k,j,i) = (beta*(C3(k,j,i)-zre(k,j,i))-C6(k,j,i)- &
                               alpha*kk2(k,j,i)*zreg(k,j,i))*rmp
                  zim(k,j,i) = (beta*(C4(k,j,i)-zim(k,j,i))+C5(k,j,i)- &
                               alpha*kk2(k,j,i)*zimg(k,j,i))*rmp
               ELSE
                  zre(k,j,i) = 0.0_GP
                  zim(k,j,i) = 0.0_GP
               ENDIF
            END DO
         END DO
      END DO

      CALL ThreeTo1D(zre,zim,to)

      RETURN
      END SUBROUTINE lin

!*****************************************************************
      SUBROUTINE newton(zre,zim,cflow,dt,vx,vy,vz,vsq,R1,C3,C4,C5,C6)
!-----------------------------------------------------------------
!
! Newton's method to find a steady state by solving:
!   eulstep(x)-x = 0.
! Starting from xguess, one iterates:
!   x_{n+1} = x_{n} + deltax
! with
!   lin deltax = rhs
! where
!   rhs = -(eulstep(x_{n})-x_{n})
! and
!   lin = d/dx (eulstep(x)-x) evaluated at x_{n}.
!
! Parameters
!     zre         : Real part of the wavefunction
!     zim         : Imaginary part of the wavefunction
!     cflow       : =1 if mean constant velocity is present
!     dt          : time step
!     vx,vy,vz,vsq: arrays with the velocity field
!     R1,R2,C3-C6 : auxiliary arrays
!     outs        : controls the amount of output
!
      USE fprecision
      USE mpivars
      USE newtmod
      USE grid
      USE var
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend)    :: vx,vy,vz
      COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(nz,ny,ista:iend) :: zre,zim
      COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(nz,ny,ista:iend) :: C3,C4
      COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(nz,ny,ista:iend) :: C5,C6
      DOUBLE PRECISION    :: errnewt,ss
      REAL(KIND=GP), INTENT(IN), DIMENSION(nx,ny,ksta:kend)    :: vsq
      REAL(KIND=GP), INTENT(INOUT), DIMENSION(nx,ny,ksta:kend) :: R1
      REAL(KIND=GP), ALLOCATABLE, DIMENSION(:) :: xguess,rhs
      REAL(KIND=GP), ALLOCATABLE, DIMENSION(:) :: deltax
      REAL(KIND=GP), INTENT(IN) :: dt
      INTEGER, INTENT(IN) :: cflow
      INTEGER             :: i,j,k,inewt
      
      ALLOCATE( xguess(1:n_dim_1d) )
      ALLOCATE( rhs(1:n_dim_1d) )
      ALLOCATE( deltax(1:n_dim_1d) )
      CALL ThreeTo1D(zre,zim,xguess)

! The Newton loop
      DO inewt=1,iter_max_newt

! Compute actual error
         CALL gtrhs(xguess,rhs,cflow,dt,vx,vy,vz,vsq,R1,C3,C4,C5,C6)
         CALL scal(rhs,rhs,ss)
         errnewt = sqrt(ss)
         IF (myrank.eq.0) &
            PRINT *,' NEWTON LOOP, iter= ',inewt,' err= ',errnewt
         IF (errnewt.lt.tol_newt) then
            IF (myrank.eq.0) print *, 'Newton converged!'
            CALL OneTo3D(xguess,zre,zim)
            DEALLOCATE ( rhs,deltax )
            RETURN
         ELSE
            CALL bicgstab(xguess,rhs,deltax,cflow,dt,errnewt,vx,vy,vz, &
                          vsq,R1,C3,C4,C5,C6)
            DO i = 1,n_dim_1d
               xguess(i)=xguess(i)+deltax(i)
            ENDDO
         ENDIF
         
      ENDDO
      
      IF (myrank.eq.0) &
           PRINT *,'Newton failed to converge in ',iter_max_newt,' steps!'
      CALL OneTo3D(xguess,zre,zim)
      DEALLOCATE ( rhs,deltax )

      RETURN
    END SUBROUTINE newton

!*****************************************************************
      SUBROUTINE scal(u1,u2,s)
!-----------------------------------------------------------------
!
! Routine to compute the reduced scalar product of two 1D
! vectors in double precision (even if GP=SINGLE).
!
! Parameters
!     u1      : First 1D vector
!     u2      : Second 1D vector
!     s       : at the output contais the reduced scalar product
!     n_dim_1d: size of the 1D vectors
!
      USE fprecision
      USE commtypes
      USE newtmod
      USE mpivars
      USE grid
!$    USE threads
      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(OUT) :: s
      DOUBLE PRECISION              :: stemp,tmp
      REAL(KIND=GP), INTENT(IN), DIMENSION(n_dim_1d) :: u1,u2
      INTEGER :: i

      stemp = 0.0D0
      tmp = 1.0_GP/  &
            (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2
!$omp parallel do reduction(+:stemp)
      DO i = 1,n_dim_1d
         stemp = stemp+u1(i)*u2(i)*tmp
      ENDDO
      CALL MPI_ALLREDUCE(stemp,s,1,MPI_DOUBLE_PRECISION,MPI_SUM, &
                         MPI_COMM_WORLD,ierr)

      RETURN
      END SUBROUTINE scal

!*****************************************************************
      SUBROUTINE zsquare(a,b,ra,rb,dealias)
!-----------------------------------------------------------------
!
! Square of the complex wavefunction Z. Note the output is not
! normalized (i.e., not divided by N^3).
!
! Parameters
!     a : real part of the wavefunction in Fourier space
!     b : imaginary part of the wavefunction in Fourier space
!     ra: real part of Z^2 (= a^2-b^2) [output]
!     rb: imag part of Z^2 (= 2 a.b  ) [output]
!     dealias: =0 does not dealias the result
!              =1 dealiases the result

      USE fprecision
      USE commtypes
      USE mpivars
      USE newtmod
      USE grid
      USE kes
      USE ali
      USE fft
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a,b
      COMPLEX(KIND=GP), DIMENSION(nz,ny,ista:iend)             :: c1,c2
      REAL(KIND=GP), INTENT(OUT), DIMENSION(nx,ny,ksta:kend)   :: ra,rb
      REAL(KIND=GP), DIMENSION(nx,ny,ksta:kend)    :: r1
      REAL(KIND=GP)       :: rmp
      INTEGER, INTENT(IN) :: dealias
      INTEGER :: i,j,k

!
! Computes real and imaginary parts
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
      CALL fftp3d_complex_to_real(plancr,c2,rb,MPI_COMM_WORLD)
      rmp = 1.0_GP/  &
            (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
         DO j = 1,ny
            DO i = 1,nx
               ra(i,j,k) = (r1(i,j,k)**2-rb(i,j,k)**2)*rmp
               rb(i,j,k) = 2*r1(i,j,k)*rb(i,j,k)*rmp
            END DO
         END DO
      END DO
!
! Dealiases the result and returns to real space
!
      IF (dealias.eq.1) THEN
         CALL fftp3d_real_to_complex(planrc,ra,c1,MPI_COMM_WORLD)
         CALL fftp3d_real_to_complex(planrc,rb,c2,MPI_COMM_WORLD)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
            DO j = 1,ny
               DO k = 1,nz
                  IF (kn2(k,j,i).gt.kmax) THEN
                     c1(k,j,i) = 0.0_GP
                     c2(k,j,i) = 0.0_GP
                  ENDIF
               END DO
            END DO
         END DO
         CALL fftp3d_complex_to_real(plancr,c1,ra,MPI_COMM_WORLD)
         CALL fftp3d_complex_to_real(plancr,c2,rb,MPI_COMM_WORLD)
      ENDIF

      RETURN
      END SUBROUTINE zsquare

!*****************************************************************
      SUBROUTINE pertgpe(ra,rb,ca,cb,a,b)
!-----------------------------------------------------------------
!
! Computes one of the linearized terms of Z.|Z|^2 in real space:
! pertgpe = Z^2.dZ* = (Zre^2 - Zim^2 + 2i Zre.Zim)*(dZre - i dZim)
!
! Parameters
!     ra : input matrix with real(Z^2) in real space (not normalized)
!     rb : input matrix with imag(Z^2) in real space (not normalized)
!     ca : input matrix with real(dZ) in Fourier space
!     cb : input matrix with imag(dZ) in Fourier space
!     a  : real part of pertgpe in Fourier space [output]
!     b  : imag part of pertgpe in Fourier space [output]
!
      USE fprecision
      USE commtypes
      USE mpivars
      USE newtmod
      USE grid
      USE fft
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN),  DIMENSION(nz,ny,ista:iend) :: ca,cb
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(nz,ny,ista:iend) :: a,b
      COMPLEX(KIND=GP), DIMENSION(nz,ny,ista:iend) :: c1,c2
      REAL(KIND=GP), INTENT(IN), DIMENSION(nx,ny,ksta:kend)     :: ra,rb
      REAL(KIND=GP), DIMENSION(nx,ny,ksta:kend)    :: r1,r2,r3
      REAL(KIND=GP)    :: rmp
      INTEGER :: i,j,k

!
! Computes real and imaginary parts
!
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
            DO k = 1,nz
               c1(k,j,i) = ca(k,j,i)
               c2(k,j,i) = cb(k,j,i)
            END DO
         END DO
      END DO
      CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,c2,r2,MPI_COMM_WORLD)
      rmp = 1.0_GP/ &
            (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
         DO j = 1,ny
            DO i = 1,nx
               r3(i,j,k) = (ra(i,j,k)*r1(i,j,k)+rb(i,j,k)*r2(i,j,k))*rmp
               r1(i,j,k) = (rb(i,j,k)*r1(i,j,k)-ra(i,j,k)*r2(i,j,k))*rmp
            END DO
         END DO
      END DO
      CALL fftp3d_real_to_complex(planrc,r3,a,MPI_COMM_WORLD)
      CALL fftp3d_real_to_complex(planrc,r1,b,MPI_COMM_WORLD)

      RETURN
      END SUBROUTINE pertgpe
      
!*****************************************************************
      SUBROUTINE bicgstab(xguess,rhs,xnew,cflow,dt,errnewt,vx,vy,vz, &
                          vsq,R1,C3,C4,C5,C6)
!-----------------------------------------------------------------
!
! Stabilized biconjugate gradient method.
! Solves lin xnew = rhs
! xnew and rhs are vectors of dimension n
! the action of linear operator lin on foo
! goo = lin(foo) is given by the call
! call lin(xguess,foo,goo,n)
! and the scalar product ss= foo1 . foo2
! is given by the call
! call scal(foo1,foo2,ss,n).
! The program exits when the following condition is met:
! sqrt((lin(xnew)-rhs).(lin(xnew)-rhs)) < tol
!
! Parameters
!     xguess: 1D vector state
!     rhs   : rhs of the ARGL equation
!     xnew  : correction to the 1D vector state
!     vx,vy,vz,vsq : arrays with the velocity field
!     R1,C3-C6     : auxiliary arrays
!
      USE fprecision
      USE commtypes
      USE mpivars
      USE newtmod
      USE grid
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend)    :: vx,vy,vz
      COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(nz,ny,ista:iend) :: C3,C4
      COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(nz,ny,ista:iend) :: C5,C6
      DOUBLE PRECISION, INTENT(IN) :: errnewt
      DOUBLE PRECISION    :: ss,tt,ts,rhonew
      DOUBLE PRECISION    :: omeganew,omegaold,rhoold
      DOUBLE PRECISION    :: alpha,beta,err
      
      REAL(KIND=GP), INTENT(IN), DIMENSION(nx,ny,ksta:kend)    :: vsq
      REAL(KIND=GP), INTENT(INOUT), DIMENSION(nx,ny,ksta:kend) :: R1
      REAL(KIND=GP), INTENT(IN), DIMENSION(n_dim_1d)    :: xguess
      REAL(KIND=GP), INTENT(INOUT), DIMENSION(n_dim_1d) :: rhs
      REAL(KIND=GP), INTENT(OUT), DIMENSION(n_dim_1d)   :: xnew
      REAL(KIND=GP), INTENT(IN) :: dt
      REAL(KIND=GP), ALLOCATABLE, DIMENSION(:) :: pold,pnew,xold,rold,rnew
      REAL(KIND=GP), ALLOCATABLE, DIMENSION(:) :: rhat0,vold,vnew,s,t
      REAL(KIND=GP)       :: tol
      INTEGER, INTENT(IN) :: cflow
      INTEGER             :: i,ilp
     
      ALLOCATE ( pold(1:n_dim_1d) )
      ALLOCATE ( pnew(1:n_dim_1d) )
      ALLOCATE ( xold(1:n_dim_1d) )
      ALLOCATE ( rold(1:n_dim_1d) )
      ALLOCATE ( rnew(1:n_dim_1d) )
      ALLOCATE ( rhat0(1:n_dim_1d) )
      ALLOCATE ( vold(1:n_dim_1d) )
      ALLOCATE ( vnew(1:n_dim_1d) )
      ALLOCATE ( s(1:n_dim_1d) )
      ALLOCATE ( t(1:n_dim_1d) )
   
! Set absolute error goal
      tol = tolbicg_rel*errnewt

! Set xold to zero
!$omp parallel do 
      DO i = 1,n_dim_1d
         xold(i) = 0.0_GP
      ENDDO
      
      CALL lin(xguess,xold,rold,cflow,dt,vx,vy,vz,vsq,R1,C3,C4,C5,C6)
      CALL gtrhs(xguess,rhs,cflow,dt,vx,vy,vz,vsq,R1,C3,C4,C5,C6)
!$omp parallel do
      DO i = 1,n_dim_1d
         rold(i) = rhs(i)-rold(i)
      ENDDO
!$omp parallel do 
      DO i = 1,n_dim_1d
         rhat0(i) = rold(i)
      ENDDO
      rhoold = 1.0D0
      alpha = 1.0D0
      omegaold = 1.0D0
!$omp parallel do
      DO i = 1,n_dim_1d
         vold(i) = 0.0_GP
         pold(i) =0.0_GP
      ENDDO
      
! This is the main bicgstab loop
      DO ilp = 1,iter_max_bicg
         CALL scal(rhat0,rold,rhonew)
         beta = (rhonew/rhoold)*(alpha/omegaold)
         
!$omp parallel do 
         DO i = 1,n_dim_1d
            pnew(i) = rold(i) + beta*(pold(i) - omegaold*vold(i))
         ENDDO

         CALL lin(xguess,pnew,vnew,cflow,dt,vx,vy,vz,vsq,R1,C3,C4,C5,C6)
         CALL scal(rhat0,vnew,ss)
         alpha = rhonew/ss

!$omp parallel do
         DO i = 1,n_dim_1d
            s(i) = rold(i) - alpha*vnew(i)
         ENDDO

         CALL lin(xguess,s,t,cflow,dt,vx,vy,vz,vsq,R1,C3,C4,C5,C6)
         call scal(s,s,ss)
         err = sqrt(ss)
         IF (myrank.eq.0) PRINT *, ' BiCgStab LOOP, iter= ',ilp,' err= ',err

         IF (err.lt.tol) THEN
!$omp parallel do
            DO i=1,n_dim_1d
               xnew(i) = xold(i) + alpha*pnew(i)
            ENDDO
            DEALLOCATE ( pold,pnew,xold,rold,rnew,rhat0,vold,vnew,s,t )
            RETURN
         ELSE
            CALL scal(t,s,ts)
            CALL scal(t,t,tt)
            omeganew = ts/tt
!$omp parallel do
            DO i=1,n_dim_1d
               xnew(i) = xold(i) + alpha*pnew(i) + omeganew*s(i)
               rnew(i) = s(i) - omeganew*t(i)
               xold(i) = xnew(i)
               rold(i) = rnew(i)
               vold(i) = vnew(i)
               pold(i) = pnew(i)
            ENDDO
            rhoold = rhonew
            omegaold = omeganew
         ENDIF

      ENDDO

! We're out of the loop without converging
      IF (myrank.eq.0) PRINT *, ' BiCgStab FAILED TO CONVERGE IN ',iter_max_bicg,' ITERATIONS'
      DEALLOCATE ( pold,pnew,xold,rold,rnew,rhat0,vold,vnew,s,t )

      RETURN
      END SUBROUTINE bicgstab
