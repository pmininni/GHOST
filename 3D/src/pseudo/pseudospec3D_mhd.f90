!=================================================================
! PSEUDOSPECTRAL subroutines
!
! Extra subroutines to compute spatial derivatives and 
! nonlinear terms in MHD and Hall-MHD equations in 3D using 
! a pseudo-spectral method. You should use the FFTPLANS and 
! MPIVARS modules (see the file 'fftp_mod.f90') in each 
! program that calls any of the subroutines in this file. 
!
! NOTATION: index 'i' is 'x' 
!           index 'j' is 'y'
!           index 'k' is 'z'
!
! 2003 Pablo D. Mininni.
!      Department of Physics, 
!      Facultad de Ciencias Exactas y Naturales.
!      Universidad de Buenos Aires.
!      e-mail: mininni@df.uba.ar
!=================================================================

!*****************************************************************
      SUBROUTINE vector3(a,b,c,d,e,f,x,y,z)
!-----------------------------------------------------------------
!
! Computes the product AxB in real space. The components 
! of the vector fields A and B are given by the matrixes 
! a,b,c,d,e and f, following the right hand convention.
!
! Parameters
!     a: input matrix with A_x
!     b: input matrix with A_y
!     c: input matrix with A_z
!     d: input matrix with B_x
!     e: input matrix with B_y
!     f: input matrix with B_z
!     x: at the output contains (AxB)_x in Fourier space
!     y: at the output contains (AxB)_y in Fourier space
!     z: at the output contains (AxB)_z in Fourier space
!
      USE fprecision
      USE commtypes
      USE mpivars
      USE grid
      USE fft
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT (IN), DIMENSION(nz,ny,ista:iend) :: a,b,c
      COMPLEX(KIND=GP), INTENT (IN), DIMENSION(nz,ny,ista:iend) :: d,e,f
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(nz,ny,ista:iend) :: x,y,z
      REAL(KIND=GP), DIMENSION(nx,ny,ksta:kend) :: r1,r2
      REAL(KIND=GP), DIMENSION(nx,ny,ksta:kend) :: r3,r4
      REAL(KIND=GP), DIMENSION(nx,ny,ksta:kend) :: r5,r6
      REAL(KIND=GP), DIMENSION(nx,ny,ksta:kend) :: r7
      REAL(KIND=GP)    :: tmp
      INTEGER :: i,j,k

!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
            DO k = 1,nz
               x(k,j,i) = a(k,j,i)
               y(k,j,i) = b(k,j,i)
               z(k,j,i) = c(k,j,i)
            END DO
         END DO
      END DO
      CALL fftp3d_complex_to_real(plancr,x,r1,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,y,r2,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,z,r3,MPI_COMM_WORLD)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
            DO k = 1,nz
               x(k,j,i) = d(k,j,i)
               y(k,j,i) = e(k,j,i)
               z(k,j,i) = f(k,j,i)
            END DO
         END DO
      END DO
      CALL fftp3d_complex_to_real(plancr,x,r4,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,y,r5,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,z,r6,MPI_COMM_WORLD)

      tmp = 1.0_GP/(real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2
!$omp parallel do if (iend-ista.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (iend-ista.lt.nth) private (i)
         DO j = 1,ny
            DO i = 1,nx
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
      END SUBROUTINE vector3

!*****************************************************************
      SUBROUTINE nonlin3(a,b,c,d,e,f,g,dir)
!-----------------------------------------------------------------
!
! Computes the nonlinear terms in 3D Navier-Stokes 
! equation with the contribution of the Lorentz 
! force. It takes the components of (v.grad)v and 
! (b.grad)b [or curl(v)xv and curl(b)xb] as input 
! matrixes and computes -(v.grad)v+(b.grad)b-grad(p) 
! [or -curl(v)xv+curl(b)xb-grad(p)] in Fourier 
! space, with the pressure chosen to satisfy the 
! incompressibility condition.
!
! Parameters
!     a  : input matrix with (v.grad)v_x or [curl(v)xv]_x
!     b  : input matrix with (v.grad)v_y or [curl(v)xv]_y
!     c  : input matrix with (v.grad)v_z or [curl(v)xv]_z
!     d  : input matrix with (b.grad)b_x or [curl(b)xb]_x
!     e  : input matrix with (b.grad)b_y or [curl(b)xb]_y
!     f  : input matrix with (b.grad)b_z or [curl(b)xb]_z
!     g  : at the output contains the result in Fourier space
!     dir: =1 computes the x-component
!          =2 computes the y-component
!          =3 computes the z-component
!
      USE fprecision
      USE kes
      USE grid
      USE mpivars
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT (IN), DIMENSION(nz,ny,ista:iend) :: a,b,c
      COMPLEX(KIND=GP), INTENT (IN), DIMENSION(nz,ny,ista:iend) :: d,e,f
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(nz,ny,ista:iend) :: g
      INTEGER, INTENT(IN) :: dir
      INTEGER             :: i,j,k

!
! Computes the x-component
!
      IF (dir.eq.1) THEN
         IF (ista.eq.1) THEN
!$omp parallel do private (k)
            DO j = 1,ny
               DO k = 1,nz
                  g(k,j,1) = -a(k,j,1)+d(k,j,1)
               END DO
            END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k)
            DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k)
               DO j = 1,ny
                  DO k = 1,nz
                     g(k,j,i) = -a(k,j,i)+d(k,j,i)+kx(i)*(kx(i)*(a(k,j,i) &
                     -d(k,j,i))+ky(j)*(b(k,j,i)-e(k,j,i))+kz(k)*(c(k,j,i) &
                     -f(k,j,i)))/kk2(k,j,i)
                  END DO
               END DO
            END DO
         ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
            DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
               DO j = 1,ny
                  DO k = 1,nz
                     g(k,j,i) = -a(k,j,i)+d(k,j,i)+kx(i)*(kx(i)*(a(k,j,i) &
                     -d(k,j,i))+ky(j)*(b(k,j,i)-e(k,j,i))+kz(k)*(c(k,j,i) &
                     -f(k,j,i)))/kk2(k,j,i)
                  END DO
               END DO
            END DO
         ENDIF
!
! Computes the y-component
!
      ELSE IF (dir.eq.2) THEN
!$omp parallel do if (iend-ista.ge.nth) private (k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth)
            DO k = 1,nz
               g(k,1,i) = -b(k,1,i)+e(k,1,i)
            END DO
         END DO
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
            DO j = 2,ny
               DO k = 1,nz
                  g(k,j,i) = -b(k,j,i)+e(k,j,i)+ky(j)*(kx(i)*(a(k,j,i) &
                  -d(k,j,i))+ky(j)*(b(k,j,i)-e(k,j,i))+kz(k)*(c(k,j,i) &
                  -f(k,j,i)))/kk2(k,j,i)
               END DO
            END DO
         END DO
!
! Computes the z-component
!
      ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth)
            DO j = 1,ny
               g(1,j,i) = -c(1,j,i)+f(1,j,i)
            END DO
         END DO
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
            DO j = 1,ny
               DO k = 2,nz
                  g(k,j,i) = -c(k,j,i)+f(k,j,i)+kz(k)*(kx(i)*(a(k,j,i) &
                  -d(k,j,i))+ky(j)*(b(k,j,i)-e(k,j,i))+kz(k)*(c(k,j,i) &
                  -f(k,j,i)))/kk2(k,j,i)
               END DO
            END DO
         END DO
      ENDIF

      RETURN
      END SUBROUTINE nonlin3

!*****************************************************************
      SUBROUTINE gauge3(a,b,c,g,dir)
!-----------------------------------------------------------------
!
! Computes the nonlinear terms in the induction 
! equation for the vector potential, imposing a 
! gauge that satisfies the condition div(A)=0.
!
! Parameters
!     a  : input matrix with (vxb)_x
!     b  : input matrix with (vxb)_y
!     c  : input matrix with (vxb)_z
!     g  : at the output contains the result in Fourier space
!     dir: =1 computes the x-component
!          =2 computes the y-component
!          =3 computes the z-component
!
      USE fprecision
      USE kes
      USE grid
      USE mpivars
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT (IN), DIMENSION(nz,ny,ista:iend) :: a,b,c
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(nz,ny,ista:iend) :: g
      INTEGER, INTENT(IN) :: dir
      INTEGER             :: i,j,k

!
! Computes the x-component
!
      IF (dir.eq.1) THEN
         IF (ista.eq.1) THEN
!$omp parallel do private (k)
            DO j = 1,ny
               DO k = 1,nz
                  g(k,j,1) = a(k,j,1)
               END DO
            END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k)
            DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k)
               DO j = 1,ny
                  DO k = 1,nz
                     g(k,j,i) = a(k,j,i)-kx(i)*(kx(i)*a(k,j,i) &
                      +ky(j)*b(k,j,i)+kz(k)*c(k,j,i))/kk2(k,j,i)
                  END DO
               END DO
            END DO
         ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
            DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
               DO j = 1,ny
                  DO k = 1,nz
                     g(k,j,i) = a(k,j,i)-kx(i)*(kx(i)*a(k,j,i) &
                      +ky(j)*b(k,j,i)+kz(k)*c(k,j,i))/kk2(k,j,i)
                  END DO
               END DO
            END DO
         ENDIF
!
! Computes the y-component
!
      ELSE IF (dir.eq.2) THEN
!$omp parallel do if (iend-ista.ge.nth) private (k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth)
            DO k = 1,nz
               g(k,1,i) = b(k,1,i)
            END DO
         END DO
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
            DO j = 2,ny
               DO k = 1,nz
                  g(k,j,i) = b(k,j,i)-ky(j)*(kx(i)*a(k,j,i) &
                   +ky(j)*b(k,j,i)+kz(k)*c(k,j,i))/kk2(k,j,i)
               END DO
            END DO
         END DO
!
! Computes the z-component
!
      ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth)
            DO j = 1,ny
               g(1,j,i) = c(1,j,i)
            END DO
         END DO
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
            DO j = 1,ny
               DO k = 2,nz
                  g(k,j,i) = c(k,j,i)-kz(k)*(kx(i)*a(k,j,i) &
                   +ky(j)*b(k,j,i)+kz(k)*c(k,j,i))/kk2(k,j,i)
                  END DO
            END DO
         END DO
      ENDIF

      RETURN
      END SUBROUTINE gauge3

!*****************************************************************
      SUBROUTINE mhdcheck(a,b,c,ma,mb,mc,t,dt,hel,crs,chk)
!-----------------------------------------------------------------
!
! Consistency check for the conservation of the total
! energy, helicity, and null divergency of the fields
!
! Output files contain:
! 'balance.txt':  time, <v^2>+<b^2>, <omega^2>, <j^2>
! 'energy.txt':   time, <v^2>, <b^2>
! 'helicity.txt': time, kinetic helicity, magnetic helicity
! 'cross.txt':    time, <v.b>, <a^2>, generalized helicity [OPT]
! 'divergence.txt' [OPTIONAL]: time, <(div.v)^2>, <(div.a)^2>
!
! Parameters
!     a  : velocity field in the x-direction
!     b  : velocity field in the y-direction
!     c  : velocity field in the z-direction
!     ma : vector potential in the x-direction
!     mb : vector potential in the y-direction
!     mc : vector potential in the z-direction
!     t  : number of time steps made
!     dt : time step
!     hel: =0 skips helicity computation
!          =1 computes the helicity
!     crs: =0 skips cross helicity computation
!          =1 computes cross helicity
!          =2 computes cross and generalized helicity
!     chk: =0 skips divergency check
!          =1 performs divergency check
!
      USE fprecision
      USE commtypes
      USE grid
      USE hall
      USE mpivars
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a,b,c
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: ma,mb,mc
      COMPLEX(KIND=GP), DIMENSION(nz,ny,ista:iend) :: c1,c2,c3
      DOUBLE PRECISION    :: engk,engm,eng,ens
      DOUBLE PRECISION    :: divk,divm,asq,crh
      DOUBLE PRECISION    :: helk,helm,cur,tmp
      DOUBLE PRECISION    :: helg
      REAL(KIND=GP)                :: tmq
      REAL(KIND=GP), INTENT(IN)    :: dt
      INTEGER, INTENT(IN) :: hel,crs,chk
      INTEGER, INTENT(IN) :: t
      INTEGER             :: i,j,k

      divk = 0.0D0
      divm = 0.0D0
      tmp = 0.0D0
      tmq = 1.0_GP/(real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2

!
! Computes the mean square value of
! the divergence of the velocity and
! vector potential field
!
      IF (chk.eq.1) THEN

      CALL derivk3(a,c1,1)
      CALL derivk3(b,c2,2)
      CALL derivk3(c,c3,3)
      IF (ista.eq.1) THEN
!$omp parallel do private (k) reduction(+:tmp)
         DO j = 1,ny
            DO k = 1,nz
               tmp = tmp+abs(c1(k,j,1)+c2(k,j,1)+c3(k,j,1))**2*tmq
            END DO
         END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k) reduction(+:tmp)
         DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k) reduction(+:tmp)
            DO j = 1,ny
               DO k = 1,nz
                  tmp = tmp+2*abs(c1(k,j,i)+c2(k,j,i)+c3(k,j,i))**2*tmq
               END DO
            END DO
         END DO
      ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k) reduction(+:tmp)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k) reduction(+:tmp)
            DO j = 1,ny
               DO k = 1,nz
                  tmp = tmp+2*abs(c1(k,j,i)+c2(k,j,i)+c3(k,j,i))**2*tmq
               END DO
            END DO
         END DO
      ENDIF
      CALL MPI_REDUCE(tmp,divk,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)

      tmp = 0.0D0
      CALL derivk3(ma,c1,1)
      CALL derivk3(mb,c2,2)
      CALL derivk3(mc,c3,3)
      IF (ista.eq.1) THEN
!$omp parallel do private (k) reduction(+:tmp)
         DO j = 1,ny
            DO k = 1,nz
               tmp = tmp+abs(c1(k,j,1)+c2(k,j,1)+c3(k,j,1))**2*tmq
            END DO
         END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k) reduction(+:tmp)
         DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k) reduction(+:tmp)
            DO j = 1,ny
               DO k = 1,nz
                  tmp = tmp+2*abs(c1(k,j,i)+c2(k,j,i)+c3(k,j,i))**2*tmq
               END DO
            END DO
         END DO
      ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k) reduction(+:tmp)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k) reduction(+:tmp)
            DO j = 1,ny
               DO k = 1,nz
                  tmp = tmp+2*abs(c1(k,j,i)+c2(k,j,i)+c3(k,j,i))**2*tmq
               END DO
            END DO
         END DO
      ENDIF
      CALL MPI_REDUCE(tmp,divm,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)

      ENDIF
!
! Computes the mean energy, enstrophy 
! and square current, and the kinetic 
! and magnetic helicity
!
      CALL energy(a,b,c,engk,1)
      CALL energy(a,b,c,ens,0)
      CALL rotor3(mb,mc,c1,1)
      CALL rotor3(ma,mc,c2,2)
      CALL rotor3(ma,mb,c3,3)
      CALL energy(c1,c2,c3,engm,1)
      CALL energy(c1,c2,c3,cur,0)
      eng = engk+engm
      IF (hel.eq.1) THEN
         CALL helicity(a,b,c,helk)
         CALL helicity(ma,mb,mc,helm)
      ENDIF
!
! Computes the square vector potential, the 
! cross helicity, and the generalized helicity 
! in Hall-MHD
!
      IF (crs.ge.1) THEN
         CALL energy(ma,mb,mc,asq,1)
         CALL cross(a,b,c,c1,c2,c3,crh,1)
      ENDIF
      IF (crs.eq.2) THEN
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
            DO j = 1,ny
               DO k = 1,nz
                  c1(k,j,i) = ma(k,j,i)+ep*a(k,j,i)
                  c2(k,j,i) = mb(k,j,i)+ep*b(k,j,i)
                  c3(k,j,i) = mc(k,j,i)+ep*c(k,j,i)
               END DO
            END DO
         END DO
         CALL helicity(c1,c2,c3,helg)
      ENDIF
!
! Creates external files to store the results
!
      IF (myrank.eq.0) THEN
         OPEN(1,file='balance.txt',position='append')
         WRITE(1,10) (t-1)*dt,eng,ens,cur
   10    FORMAT( E13.6,E22.14,E22.14,E22.14 )
         CLOSE(1)
         OPEN(1,file='energy.txt',position='append')
         WRITE(1,10) (t-1)*dt,engk,engm
         CLOSE(1)
         IF (hel.eq.1) THEN
            OPEN(1,file='helicity.txt',position='append')
            WRITE(1,10) (t-1)*dt,helk,helm
            CLOSE(1)
         ENDIF
         IF (crs.eq.1) THEN
            OPEN(1,file='cross.txt',position='append')
            WRITE(1,10) (t-1)*dt,crh,asq
            CLOSE(1)
         ELSE IF (crs.eq.2) THEN
            OPEN(1,file='cross.txt',position='append')
            WRITE(1,10) (t-1)*dt,crh,asq,helg
            CLOSE(1)
         ENDIF
         IF (chk.eq.1) THEN
            OPEN(1,file='divergence.txt',position='append')
            WRITE(1,10) (t-1)*dt,divk,divm
            CLOSE(1)
         ENDIF
      ENDIF

      RETURN
      END SUBROUTINE mhdcheck

!*****************************************************************
      SUBROUTINE crosspec(a,b,c,d,e,f,nmb)
!-----------------------------------------------------------------
!
! Computes the cross-helicity spectrum. Normalization of the
! spectrum is such that E = sum[E(k).Dkk], where Dkk is the
! width of the Fourier shells. The output is written to a
! file by the first node.
!
! Output file contains:
! 'cspectrum.XXX.txt': k, C(k)
!
! Parameters
!     a  : velocity in the x-direction
!     b  : velocity in the y-direction
!     c  : velocity in the z-direction
!     d  : magnetic potential in the x-direction
!     e  : magnetic potential in the y-direction
!     f  : magnetic potential in the z-direction
!     nmb: the extension used when writting the file
!
      USE kes
      USE grid
      USE mpivars
      USE filefmt
      USE boxsize
      USE fprecision
      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(nmax/2+1) :: Cktot
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a,b,c
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: d,e,f
      INTEGER                      :: i
      CHARACTER(len=*), INTENT(IN) :: nmb

!
! Computes the cross-helicity spectrum
!
      CALL crosspecc(a,b,c,d,e,f,Cktot,0.0_GP)
!
! Exports the spectrum to a file
!
      IF (myrank.eq.0) THEN
         OPEN(1,file='cspectrum.' // nmb // '.txt')
         DO i=1,nmax/2+1
            WRITE(1,FMT='(E13.6,E23.15)')  Dkk*i,Cktot(i)/Dkk
         END DO
         CLOSE(1)
      ENDIF

      RETURN
      END SUBROUTINE crosspec

!*****************************************************************
      SUBROUTINE crosspecc(a,b,c,d,e,f,Cktot,shift)
!-----------------------------------------------------------------
!
! Computes the cross-helicity spectrum, returning it.
!
! Parameters
!     a  : velocity in the x-direction
!     b  : velocity in the y-direction
!     c  : velocity in the z-direction
!     d  : magnetic potential in the x-direction
!     e  : magnetic potential in the y-direction
!     f  : magnetic potential in the z-direction
!     Cktot: output cross-helicity spectrum
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

      DOUBLE PRECISION, DIMENSION(nmax/2+1) :: Ck
      DOUBLE PRECISION, INTENT(OUT), DIMENSION(nmax/2+1) :: Cktot
      DOUBLE PRECISION    :: tmq
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a,b,c
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: d,e,f
      COMPLEX(KIND=GP), DIMENSION(nz,ny,ista:iend)             :: c1,c2,c3
      REAL(KIND=GP), INTENT(IN)                              :: shift
      REAL(KIND=GP)       :: tmp,rmp
      INTEGER             :: i,j,k
      INTEGER             :: kmn

!
! Computes the curl of the field if needed
!
      CALL rotor3(e,f,c1,1)
      CALL rotor3(d,f,c2,2)
      CALL rotor3(d,e,c3,3)
!
! Sets the zero for the wavenumbers
!
      rmp = shift+.501_GP
!
! Computes the cross-helicity spectrum
!
      tmp = 1.0_GP/(real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2
      DO i = 1,nmax/2+1
         Ck   (i) = 0.0D0
         Cktot(i) = 0.0D0
      END DO
      IF (ista.eq.1) THEN
!$omp parallel do private (k,kmn,tmq)
         DO j = 1,ny
            DO k = 1,nz
               kmn = int(sqrt(kk2(k,j,1))/Dkk+rmp)
               IF ((kmn.gt.0).and.(kmn.le.nmax/2+1)) THEN
                     tmq = (real(a(k,j,1)*conjg(c1(k,j,1)))+          &
                            real(b(k,j,1)*conjg(c2(k,j,1)))+          &
                            real(c(k,j,1)*conjg(c3(k,j,1))))*tmp
!$omp atomic
                     Ck(kmn) = Ck(kmn)+tmq
               ENDIF
            END DO
         END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k,kmn,tmq)
         DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k,kmn,tmq)
            DO j = 1,ny
               DO k = 1,nz
                  kmn = int(sqrt(kk2(k,j,i))/Dkk+rmp)
                  IF ((kmn.gt.0).and.(kmn.le.nmax/2+1)) THEN
                     tmq = 2*(real(a(k,j,i)*conjg(c1(k,j,i)))+     &
                              real(b(k,j,i)*conjg(c2(k,j,i)))+     &
                              real(c(k,j,i)*conjg(c3(k,j,i))))*tmp
!$omp atomic
                     Ck(kmn) = Ck(kmn)+tmq
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
                  kmn = int(sqrt(kk2(k,j,i))/Dkk+rmp)
                  IF ((kmn.gt.0).and.(kmn.le.nmax/2+1)) THEN
                     tmq = 2*(real(a(k,j,i)*conjg(c1(k,j,i)))+     &
                              real(b(k,j,i)*conjg(c2(k,j,i)))+     &
                              real(c(k,j,i)*conjg(c3(k,j,i))))*tmp
!$omp atomic
                     Ck(kmn) = Ck(kmn)+tmq
                  ENDIF
               END DO
            END DO
         END DO
      ENDIF
!
! Computes the reduction between nodes
!
      CALL MPI_ALLREDUCE(Ck,Cktot,nmax/2+1,MPI_DOUBLE_PRECISION,      &
                         MPI_SUM,MPI_COMM_WORLD,ierr)

      RETURN
      END SUBROUTINE crosspecc

