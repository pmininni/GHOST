!=================================================================
! PSEUDOSPECTRAL subroutines
!
! Extra subroutines to compute diagnostic quantities
! for the Boussinesq solvers. You should use the FFTPLANS 
! and MPIVARS modules (see the file 'fftp_mod.f90') in each 
! program that calls any of the subroutines in this file. 
!
! NOTATION: index 'i' is 'x' 
!           index 'j' is 'y'
!           index 'k' is 'z'
!
! 2011 D. Rosenberg
!      NCAR
!
!=================================================================

!*****************************************************************
      SUBROUTINE havgcomp(gsh, u, v, w, s, fo, bv, itype)
!-----------------------------------------------------------------
!
! Computes a variety of 'horizontally-averaged' quantities, which
! are functions of z. The results are stored in the return array,
! gsh, which must have a size >= NZ. The quantity to be computed
! is specified by flag itype, which can take the following 
! values:
!    itype == 0 : hor. avg. (H.A.) of shear == <(du_perp/dz)^2>_perp
!    itype == 1 : H.A. of vert. temp. gradient == <(d\theta/dz)^2>_perp
!    itype == 2 : H.A. of correlation == <u_z d\theta/dz>_perp
!    itype == 3 : H.A. of hor. kinetic energy == <u^2 + v^2>_perp
!    itype == 4 : H.A. of vert. kinetic energy == <w^2>_perp
!    itype == 5 : H.A. of perp. helicity == <u_perp.(curl u)_perp>_perp
!    itype == 6 : H.A. of correlation: <\omega_z \theta>_perp
!    itype == 7 : H.A. of PV^2 == <(f d\theta/dz - N \omega_z +
!                                 omega.grad\theta -f N)^2>_perp
!    itype == 8 : H.A. of 'super-helicity' == <\omega.curl \omega>_perp
!    itype == 9:  H.A. of gradient Richardson no. ==
!                                 N(N-d\theta/dz)/(du_\perp/dz)^2
!    itype == 10: H.A. of correlation: <u_z \theta>_perp
!
! Parameters
!     gsh: return array, funcion of z of size >= N
!     u  : input x-velocity, Fourier coeffs
!     v  : input y-velocity, Fourier coeffs
!     w  : input z-velocity, Fourier coeffs
!     s  : input scalar density/temperature, Fourier coeffs
!     fo : Coriolis parameter (usually == 2\Omega)
!     bv : Brunt-Vaisala frequency
!     itype: which quantity to compute
!
      USE fprecision
      USE commtypes
      USE kes
      USE grid
      USE mpivars
      USE filefmt
      USE fft
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: u,v
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: w,s
      COMPLEX(KIND=GP), DIMENSION(nz,ny,ista:iend)     :: c1,c2,c3,c4
      REAL(KIND=GP),    DIMENSION(nx,ny,ksta:kend)     :: r1,r2,r3
      REAL(KIND=GP),    INTENT(IN)                     :: fo, bv
      DOUBLE PRECISION, DIMENSION(nz)                  :: sh
      DOUBLE PRECISION, INTENT(INOUT), DIMENSION(nz)   :: gsh
      DOUBLE PRECISION                                 :: tmp
      INTEGER,INTENT(IN)                               :: itype
      INTEGER                                          :: i,j,k


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IF ( itype .eq. 0 .or. itype .eq. 9 ) THEN  ! <(du_perp/dz)^2>_perp
!
! Find z-derivative of u, v:
        CALL derivk3(u,c1,3)
        CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
        CALL derivk3(v,c1,3)
        CALL fftp3d_complex_to_real(plancr,c1,r2,MPI_COMM_WORLD)
!
! Do hor. average of total shear:
        sh  = 0.0D0
        gsh = 0.0D0
        ! fact 1/(nx*ny*nz)**2 for each dir * 1/(nx*ny) for horiz. avg
        tmp = 1.0D0/(dble(nx)**3*dble(ny)**3*dble(nz)**2)
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
        DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
           DO j = 1,ny
              DO i = 1,nx
!$omp atomic
                 sh(k) = sh(k)+( r1(i,j,k)**2 + r2(i,j,k)**2 )
              END DO
           END DO
           sh(k) = sh(k) * tmp
        END DO
!
! Collect as a fcn of z:
        CALL MPI_ALLREDUCE(sh,gsh,nz,MPI_DOUBLE_PRECISION,      &
                           MPI_SUM,MPI_COMM_WORLD,ierr)
        RETURN

     ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IF ( itype .eq. 1 .or. itype .eq. 2  .or. itype .eq. 9 ) THEN
      ! <(d\theta/dz)^2>_perp
!
! Find z-derivative of s:
        CALL derivk3(s,c1,3)
        CALL fftp3d_complex_to_real(plancr,c1,r3,MPI_COMM_WORLD)
!
! Do hor. average of vert. temp. gradient:
        IF ( itype .eq. 1 ) THEN
          sh  = 0.0D0
          gsh = 0.0D0
          tmp = 1.0D0/(dble(nx)**3*dble(ny)**3*dble(nz)**2)
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
          DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
             DO j = 1,ny
                DO i = 1,nx
!$omp atomic
                   sh(k) = sh(k) + r3(i,j,k)**2 
                END DO
             END DO
             sh(k) = sh(k) * tmp
          END DO
!
! Collect as a fcn of z:
          CALL MPI_ALLREDUCE(sh,gsh,nz,MPI_DOUBLE_PRECISION,    &
                             MPI_SUM,MPI_COMM_WORLD,ierr)
          RETURN
        ENDIF

     ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IF ( itype .eq. 2 ) THEN  ! <u_z d\theta/dz>_perp
! Do correlation <u_z d\theta/dz>; d\theta already computed above:
        c1 = w
        CALL fftp3d_complex_to_real(plancr,c1,r2,MPI_COMM_WORLD)
!
! Do hor. average of vert. temp. gradient:
        sh  = 0.0D0
        gsh = 0.0D0
        tmp = 1.0D0/(dble(nx)**3*dble(ny)**3*dble(nz)**2)
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
        DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
           DO j = 1,ny
              DO i = 1,nx
!$omp atomic
                 sh(k) = sh(k) + ( r1(i,j,k)*r2(i,j,k) )
              END DO
           END DO
           sh(k) = sh(k) * tmp
        END DO
!
! Collect as a fcn of z:
        CALL MPI_ALLREDUCE(sh,gsh,nz,MPI_DOUBLE_PRECISION,      &
                           MPI_SUM,MPI_COMM_WORLD,ierr)
        RETURN

      ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IF ( itype .eq. 3 ) THEN  ! <u^2 + v^2>_perp
!
! Find spatial u, v:
        c1 = u
        CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
        c1 = v
        CALL fftp3d_complex_to_real(plancr,c1,r2,MPI_COMM_WORLD)
!
! Do volume average of perp kinetic energy:
        sh  = 0.0D0
        tmp = 1.0D0/(dble(nx)**3*dble(ny)**3*dble(nz)**2)
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
        DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
           DO j = 1,ny
              DO i = 1,nx
!$omp atomic
                 sh(k) = sh(k)+( r1(i,j,k)**2 + r2(i,j,k)**2 )
              END DO
           END DO
           sh(k) = sh(k) * tmp
        END DO
!
! Collect as a fcn of z:
        CALL MPI_ALLREDUCE(sh,gsh,nz,MPI_DOUBLE_PRECISION,      &
                           MPI_SUM,MPI_COMM_WORLD,ierr)
        RETURN

     ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IF ( itype .eq. 4 ) THEN  ! <w^2 >_perp
!
! Find spatial u_z:
        c1 = w
        CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
!
! Do volume average of perp kinetic energy:
        sh  = 0.0D0
        tmp = 1.0D0/(dble(nx)**3*dble(ny)**3*dble(nz)**2)
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
        DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
           DO j = 1,ny
              DO i = 1,nx
!$omp atomic
                 sh(k) = sh(k)+( r1(i,j,k)**2 )
              END DO
           END DO
           sh(k) = sh(k) * tmp
        END DO
!
! Collect as a fcn of z:
        CALL MPI_ALLREDUCE(sh,gsh,nz,MPI_DOUBLE_PRECISION,      &
                           MPI_SUM,MPI_COMM_WORLD,ierr)
        RETURN

      ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IF ( itype .eq. 5 ) THEN
      ! <u_perp . (curl u)_perp >_perp = <-u dv/dz + v du/dz>_perp
!
! Find z-derivative of v, IFFT of u:
        CALL derivk3(v,c1,3)
        CALL fftp3d_complex_to_real(plancr,c1,r2,MPI_COMM_WORLD)
        c1 = u
        CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
!
! Compute -u dv/dz contrib:
        sh  = 0.0D0
        gsh = 0.0D0
        tmp = 1.0D0/(dble(nx)**3*dble(ny)**3*dble(nz)**2)
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
        DO k = ksta,kend 
!$omp parallel do if (kend-ksta.lt.nth) private (i)
           DO j = 1,ny
              DO i = 1,nx
!$omp atomic
                 sh(k) = sh(k)-( r1(i,j,k)*r2(i,j,k) )
              END DO
           END DO
        END DO
!
! Compute v du/dz contrib:
        CALL derivk3(u,c1,3)
        CALL fftp3d_complex_to_real(plancr,c1,r2,MPI_COMM_WORLD)
        c1 = v
        CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
        DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
           DO j = 1,ny
              DO i = 1,nx
!$omp atomic
                 sh(k) = sh(k)+( r1(i,j,k)*r2(i,j,k) )
              END DO
           END DO
           sh(k) = sh(k) * tmp
        END DO
!
! Collect as a fcn of z:
        CALL MPI_ALLREDUCE(sh,gsh,nz,MPI_DOUBLE_PRECISION,      &
                           MPI_SUM,MPI_COMM_WORLD,ierr)
        RETURN

      ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IF ( itype .eq. 6 ) THEN  ! <\omega_z \theta>_perp
        sh  = 0.0D0
        gsh = 0.0D0
        tmp = 1.0D0/(dble(nx)**3*dble(ny)**3*dble(nz)**2)
        CALL rotor3(u,v,c1,3)
        CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
        c1 = s
        CALL fftp3d_complex_to_real(plancr,c1,r2,MPI_COMM_WORLD)

!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
        DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
           DO j = 1,ny
              DO i = 1,nx
!$omp atomic
                 sh(k) = sh(k)+r1(i,j,k)*r2(i,j,k)
              END DO
           END DO
           sh(k) = sh(k) * tmp
        END DO
!
! Collect as a fcn of z:
        CALL MPI_ALLREDUCE(sh,gsh,nz,MPI_DOUBLE_PRECISION,      &
                           MPI_SUM,MPI_COMM_WORLD,ierr)
        RETURN

      ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IF ( itype .eq. 7 ) THEN  ! <PV^2>_perp
        r2  = 0.0_GP
        sh  = 0.0D0
        gsh = 0.0D0
!
! Find z-derivative of theta contrib:
        CALL derivk3(s,c1,3)
        CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
        ! fac 1/(nx*ny*nz)**2 for each dir
        tmp = 1.0D0/(dble(nx)*dble(ny)*dble(nz))
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
        DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
           DO j = 1,ny
              DO i = 1,nx
!$omp atomic
                 r2(i,j,k) = r2(i,j,k) + fo*r1(i,j,k)*tmp
              END DO
           END DO
        END DO
!
! Find omega_z contrib:
        CALL rotor3(u,v,c1,3) ! omega_z
        CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
        tmp = 1.0D0/(dble(nx)*dble(ny)*dble(nz))
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
        DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
           DO j = 1,ny
              DO i = 1,nx
!$omp atomic
                 r2(i,j,k) = r2(i,j,k) - bv*r1(i,j,k)*tmp
              END DO
           END DO
        END DO
!
!  Find omega.Grad theta contrib:
        CALL rotor3(v,w,c1,1)
        CALL rotor3(u,w,c2,2)
        CALL rotor3(u,v,c3,3)
        CALL advect3(c1,c2,c3,s,c4)
        ! what comes out is -omega.Grad theta, so tmp<0:
        CALL fftp3d_complex_to_real(plancr,c4,r1,MPI_COMM_WORLD)
        tmp = -1.0D0/(dble(nx)*dble(ny)*dble(nz))
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
        DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
           DO j = 1,ny
              DO i = 1,nx
!$omp atomic
                 r2(i,j,k) = r2(i,j,k) + r1(i,j,k)*tmp
              END DO
           END DO
        END DO
!
! Find fN contrib:
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
        DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
           DO j = 1,ny
              DO i = 1,nx
!$omp atomic
                 r2(i,j,k) = r2(i,j,k) - fo*bv 
              END DO
           END DO
        END DO
        ! fac 1/(nx*ny) for horiz. avg
        tmp = 1.0D0/(dble(nx)*dble(ny))
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
        DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
           DO j = 1,ny
              DO i = 1,nx
!$omp atomic
                 sh(k) = sh(k) + r2(i,j,k)**2
              END DO
           END DO
           sh(k) = sh(k)*tmp
        END DO
!
! Collect as a fcn of z:
        CALL MPI_ALLREDUCE(sh,gsh,nz,MPI_DOUBLE_PRECISION,      &
                           MPI_SUM,MPI_COMM_WORLD,ierr)
        RETURN
        
      ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IF ( itype .eq. 8 ) THEN
      ! <\omega . curl \omega>_perp == <super-helicity>_perp
        sh  = 0.0D0
        gsh = 0.0D0
        tmp = 1.0D0/(dble(nx)**3*dble(ny)**3*dble(nz)**2)

!       ...do omega_x * ( d\omega_z/dy - d\omega_y/dz) contrib:
        CALL rotor3(u,v,c1,3)
        CALL derivk3(c1,c2,2)
        CALL fftp3d_complex_to_real(plancr,c2,r1,MPI_COMM_WORLD)
        CALL rotor3(u,w,c1,2)
        CALL derivk3(c1,c2,3)
        CALL fftp3d_complex_to_real(plancr,c2,r2,MPI_COMM_WORLD)

        CALL rotor3(v,w,c1,1) ! omega_x
        CALL fftp3d_complex_to_real(plancr,c1,r3,MPI_COMM_WORLD)
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
        DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
           DO j = 1,ny
              DO i = 1,nx
!$omp atomic
                 sh(k) = sh(k) + r3(i,j,k)*( r1(i,j,k)-r2(i,j,k) )
              END DO
           END DO
        END DO

!       ...do omega_y * ( d\omega_x/dz - d\omega_z/dx) contrib:
        CALL rotor3(v,w,c1,1) ! omega_x
        CALL derivk3(c1,c2,3)
        CALL fftp3d_complex_to_real(plancr,c2,r1,MPI_COMM_WORLD)
        CALL rotor3(u,v,c1,3) ! omega_z
        CALL derivk3(c1,c2,1)
        CALL fftp3d_complex_to_real(plancr,c2,r2,MPI_COMM_WORLD)

        CALL rotor3(u,w,c1,2) ! omega_y
        CALL fftp3d_complex_to_real(plancr,c1,r3,MPI_COMM_WORLD)
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
        DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
           DO j = 1,ny
              DO i = 1,nx
!$omp atomic
                 sh(k) = sh(k) + r3(i,j,k)*( r1(i,j,k)-r2(i,j,k) )
              END DO
           END DO
        END DO

!       ...do omega_z * ( d\omega_y/dx - d\omega_x/dy) contrib:
        CALL rotor3(u,w,c1,2) ! omega_y
        CALL derivk3(c1,c2,1)
        CALL fftp3d_complex_to_real(plancr,c2,r1,MPI_COMM_WORLD)
        CALL rotor3(v,w,c1,1) ! omega_x
        CALL derivk3(c1,c2,2)
        CALL fftp3d_complex_to_real(plancr,c2,r2,MPI_COMM_WORLD)

        CALL rotor3(u,v,c1,3) ! omega_z
        CALL fftp3d_complex_to_real(plancr,c1,r3,MPI_COMM_WORLD)
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
        DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
           DO j = 1,ny
              DO i = 1,nx
!$omp atomic
                 sh(k) = sh(k) + r3(i,j,k)*( r1(i,j,k)-r2(i,j,k) )
              END DO
           END DO
           sh(k) = sh(k) * tmp
        END DO
!
! Collect as a fcn of z:
        CALL MPI_ALLREDUCE(sh,gsh,nz,MPI_DOUBLE_PRECISION,      &
                           MPI_SUM,MPI_COMM_WORLD,ierr)
        RETURN
        
      ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IF ( itype .eq. 9 ) THEN  ! Ri = N (N-d\theta/dz)/(du_\perp/dz) 
!
!      Note, d\theta/dz should be in r3 already
!      du_x/dz should be in r1 already, and du_y/dz in r2 from above
!
       sh  = 0.0D0
       gsh = 0.0D0
       tmp = 1.0D0/(dble(nx)**3*dble(ny)**3*dble(nz)**2)

!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
        DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
           DO j = 1,ny
              DO i = 1,nx
!$omp atomic
                 sh(k) = sh(k)+bv*(bv-r3(i,j,k))/(r1(i,j,k)**2+r2(i,j,k)**2) 
              END DO
           END DO
           sh(k) = sh(k) * tmp
        END DO
!
! Collect as a fcn of z:
        CALL MPI_ALLREDUCE(sh,gsh,nz,MPI_DOUBLE_PRECISION,      &
                           MPI_SUM,MPI_COMM_WORLD,ierr)
        RETURN

      ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IF ( itype .eq. 10 ) THEN  ! <v_z \theta>_perp
        sh  = 0.0D0
        gsh = 0.0D0
        tmp = 1.0D0/(dble(nx)**3*dble(ny)**3*dble(nz)**2)
        c1 = w
        CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
        c1 = s
        CALL fftp3d_complex_to_real(plancr,c1,r2,MPI_COMM_WORLD)

!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
        DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
           DO j = 1,ny
              DO i = 1,nx
!$omp atomic
                 sh(k) = sh(k)+r1(i,j,k)*r2(i,j,k)
              END DO
           END DO
           sh(k) = sh(k) * tmp
        END DO
!
! Collect as a fcn of z:
        CALL MPI_ALLREDUCE(sh,gsh,nz,MPI_DOUBLE_PRECISION,      &
                           MPI_SUM,MPI_COMM_WORLD,ierr)

        RETURN

      ENDIF

      RETURN
      END SUBROUTINE havgcomp

!*****************************************************************
      SUBROUTINE havgwrite(itype,spref,nmb,u,v,w,s,fo,bv)
!-----------------------------------------------------------------
!
! Computes horizontal average of quantity itype (specified in 
! method 'havgcomp'), and outputs it as a function of z, to the
! file whose prefix is 'spref'. Filename will be of the form
! spref.XXX.txt, where XXX is computed from output interval id 'nmb'.
!
! Output file contains:
! 'spref.XXX.txt': z, mean profile of quantity 'itpype'
!         
! Parameters
!     itype : quantity id
!     spref : filename prefix
!     nmb   : output interval id extension
!     u     : input x-velocity, Fourier coeffs
!     v     : input y-velocity, Fourier coeffs
!     w     : input z-velocity, Fourier coeffs
!     s     : input scalar density/temperature, Fourier coeffs
!     fo    : Coriolis parameter (usually == 2\Omega)
!     bv    : Brunt-Vaisala frequency
!
      USE fprecision
      USE commtypes
      USE var
      USE kes
      USE grid
      USE mpivars
      USE filefmt
      USE boxsize
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: u,v
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: w,s
      REAL(KIND=GP),    INTENT(IN)      :: fo, bv
      DOUBLE PRECISION, DIMENSION(nz)   :: gsh
      CHARACTER(len=*), INTENT(IN)      :: nmb,spref
      INTEGER,INTENT(IN)                :: itype
      INTEGER                           :: k

      CALL havgcomp(gsh,u,v,w,s,fo,bv,itype) 
      IF (myrank.eq.0) THEN
         OPEN(1,file=trim(spref) // '.' // trim(nmb) // '.txt')
         DO k = 1,nz
         WRITE(1,10) 2*pi*Lz*(real(k,kind=GP)-1)/real(nz,kind=GP), gsh(k)
         END DO
         CLOSE(1)
   10    FORMAT( E23.15,E23.15 )
      ENDIF

      RETURN
      END SUBROUTINE havgwrite

!*****************************************************************
      SUBROUTINE tbouss(u, v, w, s, t, dt, fo, bv)
!-----------------------------------------------------------------
!
! Computes the volume-average and max/min of quantities computed in
! havgcomp, and outputs them to disk as a function of time.
!
! Output files contain:
! 'tboussavg.txt': time, vol. average of quantities described below
! 'tboussmax.txt': time, max. val. of quantities described below
! 'tboussmin.txt': time, min. val. of quantities described below
!    col== 2 : hor. avg. of shear == <du_perp/dz>^2
!    col== 3 : hor. avg. of vert. temp. gradient == <d\theta/dz>
!    col== 4 : hor. avg. of correlation == <u_z d\theta/dz>
!    col== 5 : hor. avg. of hor. kinetic energy == <u^2 + v^2>
!    col== 6 : hor. avg. of vert. kinetic energy == <w^2>
!    col== 7 : hor. avg. of perp. helicity == <u_perp.(curl u)_perp>_perp
!    col== 8 : hor. avg. of correlation: <\omega_z \theta>
!    col== 9 : hor. avg. of PV^2 (see 'havgcomp')
!    col== 10: hor. avg. of 'super-helicity': <\omega . curl \omega>_perp
!    col== 11: hor. avg. of gradient Richardson no. (see 'havgcomp')
!    col== 12: hor. avg. of correlation: <u_z \theta>
!        
! Parameters
!     u  : input x-velocity, Fourier coeffs
!     v  : input y-velocity, Fourier coeffs
!     w  : input z-velocity, Fourier coeffs
!     s  : input scalar density/temperature, Fourier coeffs
!     t  : number of time steps made
!     dt : time step
!     fo : Coriolis parameter (usually == 2\Omega)
!     bv : Brunt-Vaisala frequency
!
      USE fprecision
      USE commtypes
      USE kes
      USE ali
      USE grid
      USE mpivars
      USE filefmt
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: u,v
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: w,s
      REAL(KIND=GP),    INTENT(IN)    :: dt
      REAL(KIND=GP),    INTENT(IN)    :: fo, bv
      DOUBLE PRECISION                :: gxavg(11),gxmax(11),gxmin(11)
      DOUBLE PRECISION                :: tmp
      DOUBLE PRECISION, DIMENSION(nz) :: gsh
      INTEGER, INTENT(IN)             :: t
      INTEGER                         :: i

      gxavg  = 0.0D0
      gxmax  = 0.0D0
      gxmin  = 1.0D0/tinyf

! Compute volume averages from horizontal averages, and
! also extrema from horiz. averages:
      DO i = 1, 11
        CALL havgcomp(gsh,u,v,w,s,fo,bv,i-1)
        gxavg(i) = sum   (gsh)/dble(nz);
        gxmax(i) = maxval(gsh);
        gxmin(i) = minval(gsh);
      ENDDO
!
! Don't need reductions bec the results from havgcomp are globalized,
! so just do the outputs:
!
! Output averaged quantities as a fcn of t:
      IF (myrank.eq.0) THEN
         OPEN(1,file='tboussavg.txt',position='append')
         WRITE(1,20) (t-1)*dt,gxavg(1),gxavg (2),gxavg (3),gxavg(4), &
                              gxavg(5),gxavg (6),gxavg (7),gxavg(8), &
                              gxavg(9),gxavg(10),gxavg(11)
         CLOSE(1)
      ENDIF
!
! Output max quantities as a fcn of t:
      IF (myrank.eq.0) THEN
         OPEN(1,file='tboussmax.txt',position='append')
         WRITE(1,20) (t-1)*dt,gxmax(1),gxmax (2),gxmax (3),gxmax(4), &
                              gxmax(5),gxmax (6),gxmax (7),gxmax(8), &
                              gxmax(9),gxmax(10),gxmax(11)
         CLOSE(1)
      ENDIF
!
! Output min quantities as a fcn of t:
      IF (myrank.eq.0) THEN
         OPEN(1,file='tboussmin.txt',position='append')
         WRITE(1,20) (t-1)*dt,gxmin(1),gxmin (2),gxmin (3),gxmin(4), &
                              gxmin(5),gxmin (6),gxmin (7),gxmin(8), &
                              gxmin(9),gxmin(10),gxmin(11)
         CLOSE(1)
      ENDIF

      RETURN
   20 FORMAT( E13.6,1x,11(E26.18,1x) )
      END SUBROUTINE tbouss

!*****************************************************************
      SUBROUTINE spectpv(u,v,w,s,nmb)
!-----------------------------------------------------------------
!
! Computes the potential vorticity power spectrum as a function
! of kz, wigh kz = Dkz*(0,...,nz/2). The pot'l vorticity is
! computed as:
!  PV = curl(v) . Grad(theta)
! Output is written to a file by the first node.
!
! Output file contains:
! 'pvspectrum.XXX.txt': kz, PV(kz) (power spectrum of pot. vort.)
!
! Parameters
!     u,v,w  : input matrix with the 1,2,3 velocity components;
!     s      : input matrix with the scalar density/temperature
!     nmb    : the extension used when writting the file
!
      USE fprecision
      USE commtypes
      USE kes
      USE grid
      USE mpivars
      USE filefmt
      USE boxsize
      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(nz/2+1)                 :: Ek,Ektot
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: u,v
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: w,s
      COMPLEX(KIND=GP),             DIMENSION(nz,ny,ista:iend) :: c1,c2
      COMPLEX(KIND=GP),             DIMENSION(nz,ny,ista:iend) :: c3,a
      REAL(KIND=GP)     :: tmp
      INTEGER           :: i,j,k
      INTEGER           :: kmn
      CHARACTER(len=*), INTENT(IN) :: nmb
!
! Compute vorticity:
      CALL rotor3(v,w,c1,1)
      CALL rotor3(u,w,c2,2)
      CALL rotor3(u,v,c3,3)
!
! Compute curl a = v . Grad(s):
      CALL advect3(c1,c2,c3,s,a)
!
! Set Ek to zero
      DO i = 1,nz/2+1
         Ek(i) = 0.
      END DO
!
! Computes the power spectrum
!
      tmp = 1.0_GP/ &
            (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2
      IF (ista.eq.1) THEN
         DO j = 1,ny
            DO k = 1,nz
               kmn = int(abs(kz(k))*Lz+1)
               IF ((kmn.gt.0).and.(kmn.le.nz/2+1)) THEN
                  Ek(kmn) = Ek(kmn)+tmp*abs(a(k,j,1))**2
               ENDIF
            END DO
         END DO
         DO i = 2,iend
            DO j = 1,ny
               DO k = 1,nz
                  kmn = int(abs(kz(k))*Lz+1)
                  IF ((kmn.gt.0).and.(kmn.le.nz/2+1)) THEN
                     Ek(kmn) = Ek(kmn)+2*tmp*abs(a(k,j,i))**2
                  ENDIF
               END DO
            END DO
         END DO
      ELSE
         DO i = ista,iend
            DO j = 1,ny
               DO k = 1,nz
                  kmn = int(abs(kz(k))*Lz+1)
                  IF ((kmn.gt.0).and.(kmn.le.nz/2+1)) THEN
                     Ek(kmn) = Ek(kmn)+2*tmp*abs(a(k,j,i))**2
                  ENDIF
               END DO
            END DO
         END DO
      ENDIF
!
! Computes the reduction between nodes
!
      CALL MPI_REDUCE(Ek,Ektot,nz/2+1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
!
! Exports the spectrum to a file
!
      IF (myrank.eq.0) THEN
         OPEN(1,file='pvspectrum.' // nmb // '.txt')
         DO k = 1,nz/2+1
            WRITE(1,FMT='(E13.6,E23.15)') Dkz*(k-1),Ektot(k)
         END DO
         CLOSE(1)
      ENDIF

      RETURN
      END SUBROUTINE spectpv

