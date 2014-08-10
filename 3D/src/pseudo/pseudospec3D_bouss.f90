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
! gsh, which must have a size >= N. The quantity to be computed
! is specified by flag itype, which can take the following 
! values:
!    itype == 0 : hor. avg. of shear == <du_perp/dz>_perp ^2
!    itype == 1 : hor. avg. of vert. temp. gradient == <d\theta/dz>_perp
!    itype == 2 : hor. avg. of correlation == <u_z d\theta/dz>_perp
!    itype == 3 : hor. avg. of hor. kinetic energy == <u^2 + v^2>_perp
!    itype == 4 : hor. avg. of vert. kinetic energy == <w^2>_perp
!    itype == 5 : hor. avg. of perp. helicity == <u_perp . curl u_perp = -u dv/dz + v du/dz>_perp
!    itype == 6 : hor. avg. of correlation: <\omega_z \theta>_perp
!    itype == 7 : hor. avg. of PV^2: <(fd\theta/dz - N \omega_z + omega.Grad\theta -fN)^2>_perp
!    itype == 8 : hor. avg. of 'super-helicity': <\omega . curl \omega>_perp 
!    itype == 9: hor. avg. of gradient Richardson no: !    N(N-d\theta/dz)/(du_\perp/dz)^2
!    itype == 10: hor. avg. of correlation: <u_z \theta>_perp
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
      USE fprecision
      USE commtypes
      USE kes
      USE grid
      USE mpivars
      USE filefmt
      USE fft
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,n,ista:iend) :: u,v,w,s
      COMPLEX(KIND=GP), DIMENSION(n,n,ista:iend) :: c1,c2,c3,c4
      REAL(KIND=GP), DIMENSION(n,n,ksta:kend)    :: r1,r2,r3
      REAL(KIND=GP),INTENT(IN)                   :: fo, bv
      DOUBLE PRECISION, DIMENSION(n)             :: sh
      DOUBLE PRECISION, INTENT(INOUT), DIMENSION(n):: gsh
      DOUBLE PRECISION                           :: tmp
      INTEGER,INTENT(IN)                         :: itype
      INTEGER                                    :: i,j,k


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1111!1111!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1111!1111!!!

      IF ( itype .eq. 0 .or. itype .eq. 9 ) THEN  ! <du_perp/dz>_perp^2
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
        tmp = 1.0D0/dble(n)**8 ! fact of 1/n^3 for each factor * 1/n^2 for horiz. avg
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
        DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
           DO j = 1,n
              DO i = 1,n
!$omp atomic
                 sh(k) = sh(k)+( r1(i,j,k)**2 + r2(i,j,k)**2 )
              END DO
           END DO
           sh(k) = sh(k) * tmp
        END DO

!
! Collect as a fcn of z:
        CALL MPI_ALLREDUCE(sh,gsh,n,MPI_DOUBLE_PRECISION,      &
                           MPI_SUM,MPI_COMM_WORLD,ierr)

        RETURN
      ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1111!1111!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1111!1111!!!

      IF ( itype .eq. 1 .or. itype .eq. 2  .or. itype .eq. 9 ) THEN  ! <d\theta/dz>_perp
!
! Find z-derivative of s:
        CALL derivk3(s,c1,3)
        CALL fftp3d_complex_to_real(plancr,c1,r3,MPI_COMM_WORLD)
!
! Do hor. average of vert. temp. gradient:
        IF ( itype .eq. 1 ) THEN
          sh  = 0.0D0
          gsh = 0.0D0
          tmp = 1.0D0/dble(n)**8 ! fact of 1/n^3 for r1  * 1/n^2 for horiz. avg
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
          DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
             DO j = 1,n
              DO i = 1,n
!$omp atomic
                   sh(k) = sh(k) + r3(i,j,k)**2 
                END DO
             END DO
             sh(k) = sh(k) * tmp
          END DO
!
! Collect as a fcn of z:
          CALL MPI_ALLREDUCE(sh,gsh,n,MPI_DOUBLE_PRECISION,      &
                             MPI_SUM,MPI_COMM_WORLD,ierr)
          RETURN

        ENDIF
      ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1111!1111!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1111!1111!!!

      IF ( itype .eq. 2 ) THEN  ! <u_z d\theta/dz>_perp
! Do correlation <u_z d\theta/dz>; d\theta already computed above:
        c1 = w
        CALL fftp3d_complex_to_real(plancr,c1,r2,MPI_COMM_WORLD)
!
! Do hor. average of vert. temp. gradient:
        sh  = 0.0D0
        gsh = 0.0D0
        tmp = 1.0D0/dble(n)**8 ! fact of 1/n^3 for each factor * 1/n^2 for horiz. avg
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
        DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
           DO j = 1,n
              DO i = 1,n
!$omp atomic
                 sh(k) = sh(k) + ( r1(i,j,k)*r2(i,j,k) )
              END DO
           END DO
           sh(k) = sh(k) * tmp
        END DO
!
! Collect as a fcn of z:
        CALL MPI_ALLREDUCE(sh,gsh,n,MPI_DOUBLE_PRECISION,      &
                           MPI_SUM,MPI_COMM_WORLD,ierr)

        RETURN
      ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1111!1111!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1111!1111!!!

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
        tmp = 1.0D0/dble(n)**8
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
        DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
           DO j = 1,n
              DO i = 1,n
!$omp atomic
                 sh(k) = sh(k)+( r1(i,j,k)**2 + r2(i,j,k)**2 )
              END DO
           END DO
           sh(k) = sh(k) * tmp
        END DO

!
! Collect as a fcn of z:
        CALL MPI_ALLREDUCE(sh,gsh,n,MPI_DOUBLE_PRECISION,      &
                           MPI_SUM,MPI_COMM_WORLD,ierr)
        RETURN
      ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1111!1111!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1111!1111!!!

      IF ( itype .eq. 4 ) THEN  ! <w^2 >_perp
!
! Find spatial u_z:
        c1 = w
        CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)

!
! Do volume average of perp kinetic energy:
        sh  = 0.0D0
        tmp = 1.0D0/dble(n)**8
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
        DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
           DO j = 1,n
              DO i = 1,n
!$omp atomic
                 sh(k) = sh(k)+( r1(i,j,k)**2 )
              END DO
           END DO
           sh(k) = sh(k) * tmp
        END DO

!
! Collect as a fcn of z:
        CALL MPI_ALLREDUCE(sh,gsh,n,MPI_DOUBLE_PRECISION,      &
                           MPI_SUM,MPI_COMM_WORLD,ierr)
        RETURN
      ENDIF


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1111!1111!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1111!1111!!!

      IF ( itype .eq. 5 ) THEN  ! <u_perp . curl u_perp = -u dv/dz + v du/dz>_perp
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
        tmp = 1.0D0/dble(n)**8
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
        DO k = ksta,kend 
!$omp parallel do if (kend-ksta.lt.nth) private (i)
           DO j = 1,n
              DO i = 1,n
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
           DO j = 1,n
              DO i = 1,n
!$omp atomic
                 sh(k) = sh(k)+( r1(i,j,k)*r2(i,j,k) )
              END DO
           END DO
           sh(k) = sh(k) * tmp
        END DO

!
! Collect as a fcn of z:
        CALL MPI_ALLREDUCE(sh,gsh,n,MPI_DOUBLE_PRECISION,      &
                           MPI_SUM,MPI_COMM_WORLD,ierr)

        RETURN
      ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1111!1111!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1111!1111!!!

      IF ( itype .eq. 6 ) THEN  ! <\omega_z \theta>_perp
        sh  = 0.0D0
        gsh = 0.0D0
        tmp = 1.0D0/dble(n)**8
        CALL rotor3(u,v,c1,3)
        CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
        c1 = s
        CALL fftp3d_complex_to_real(plancr,c1,r2,MPI_COMM_WORLD)

!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
        DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
           DO j = 1,n
              DO i = 1,n
!$omp atomic
                 sh(k) = sh(k)+r1(i,j,k)*r2(i,j,k)
              END DO
           END DO
           sh(k) = sh(k) * tmp
        END DO
!
! Collect as a fcn of z:
        CALL MPI_ALLREDUCE(sh,gsh,n,MPI_DOUBLE_PRECISION,      &
                           MPI_SUM,MPI_COMM_WORLD,ierr)

        RETURN
      ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1111!1111!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1111!1111!!!

      IF ( itype .eq. 7 ) THEN  ! <PV^2>_perp
        r2  = 0.0_GP
        sh  = 0.0D0
        gsh = 0.0D0
!
! Find z-derivative of theta contrib:
        CALL derivk3(s,c1,3)
        CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)

        tmp = 1.0D0/dble(n)**3
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
        DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
           DO j = 1,n
              DO i = 1,n
!$omp atomic
                 r2(i,j,k) = r2(i,j,k) + fo*r1(i,j,k)*tmp
              END DO
           END DO
        END DO

! Find omega_z contrib:
        CALL rotor3(u,v,c1,3) ! omega_z
        CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
        tmp = 1.0D0/dble(n)**3
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
        DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
           DO j = 1,n
              DO i = 1,n
!$omp atomic
                 r2(i,j,k) = r2(i,j,k) - bv*r1(i,j,k)*tmp
              END DO
           END DO
        END DO

!  Find omega.Grad theta contrib:
      CALL rotor3(v,w,c1,1)
      CALL rotor3(u,w,c2,2)
      CALL rotor3(u,v,c3,3)
      CALL advect3(c1,c2,c3,s,c4) ! what comes out is -omega.Grad theta, so tmp<0:
      CALL fftp3d_complex_to_real(plancr,c4,r1,MPI_COMM_WORLD)
      tmp = -1.0D0/dble(n)**3
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
        DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
           DO j = 1,n
              DO i = 1,n
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
           DO j = 1,n
              DO i = 1,n
!$omp atomic
                 r2(i,j,k) = r2(i,j,k) - fo*bv 
              END DO
           END DO
        END DO

      tmp = 1.0D0/dble(n)**2
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
        DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
           DO j = 1,n
              DO i = 1,n
!$omp atomic
                 sh(k) = sh(k) + r2(i,j,k)**2
              END DO
           END DO
           sh(k) = sh(k)*tmp
        END DO
!
! Collect as a fcn of z:
        CALL MPI_ALLREDUCE(sh,gsh,n,MPI_DOUBLE_PRECISION,      &
                           MPI_SUM,MPI_COMM_WORLD,ierr)

        RETURN
      ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1111!1111!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1111!1111!!!

      IF ( itype .eq. 8 ) THEN  ! <\omega . curl \omega>_perp == <super-helicity>_perp
        sh  = 0.0D0
        gsh = 0.0D0
        tmp = 1.0D0/dble(n)**8

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
           DO j = 1,n
              DO i = 1,n
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
           DO j = 1,n
              DO i = 1,n
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
           DO j = 1,n
              DO i = 1,n
!$omp atomic
                 sh(k) = sh(k) + r3(i,j,k)*( r1(i,j,k)-r2(i,j,k) )
              END DO
           END DO
           sh(k) = sh(k) * tmp
        END DO
!
! Collect as a fcn of z:
        CALL MPI_ALLREDUCE(sh,gsh,n,MPI_DOUBLE_PRECISION,      &
                           MPI_SUM,MPI_COMM_WORLD,ierr)

        RETURN
      ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1111!1111!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1111!1111!!!

      IF ( itype .eq. 9 ) THEN  ! Ri = N (N-d\theta/dz)/(du_\perp/dz) 
!
!      Note, d\theta/dz should be in r3 already
!      du_x/dz should be in r1 already, and du_y/dz in r2 from above
!
        sh  = 0.0D0
        gsh = 0.0D0
        tmp = 1.0D0/dble(n)**8 ! fact of 1/n^3 for r1  * 1/n^2 for horiz. avg
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
        DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
           DO j = 1,n
            DO i = 1,n
!$omp atomic
                 sh(k) = sh(k) + bv*(bv-r3(i,j,k))/ (r1(i,j,k)**2+r2(i,j,k)**2) 
              END DO
           END DO
           sh(k) = sh(k) * tmp
        END DO
!
! Collect as a fcn of z:
        CALL MPI_ALLREDUCE(sh,gsh,n,MPI_DOUBLE_PRECISION,      &
                           MPI_SUM,MPI_COMM_WORLD,ierr)
        RETURN

      ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1111!1111!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1111!1111!!!

      IF ( itype .eq. 10 ) THEN  ! <v_z \theta>_perp
        sh  = 0.0D0
        gsh = 0.0D0
        tmp = 1.0D0/dble(n)**8
        c1 = w
        CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
        c1 = s
        CALL fftp3d_complex_to_real(plancr,c1,r2,MPI_COMM_WORLD)

!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
        DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
           DO j = 1,n
              DO i = 1,n
!$omp atomic
                 sh(k) = sh(k)+r1(i,j,k)*r2(i,j,k)
              END DO
           END DO
           sh(k) = sh(k) * tmp
        END DO
!
! Collect as a fcn of z:
        CALL MPI_ALLREDUCE(sh,gsh,n,MPI_DOUBLE_PRECISION,      &
                           MPI_SUM,MPI_COMM_WORLD,ierr)

        RETURN
      ENDIF

      RETURN
      END SUBROUTINE havgcomp


!*****************************************************************
      SUBROUTINE havgwrite(itype,spref,nmb,u,v,w,s,fo,bv)
!-----------------------------------------------------------------
!
! Computes horizontal average of quantity, itype (specified in 
! method 'havgcomp', and outputs it as a function of z, to the
! file whose prefix is 'spref'. Filename will be of the form
! spref.XXX.txt, where XXX is computed from output interval id, 'nmb'
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


      USE fprecision
      USE commtypes
      USE kes
      USE grid
      USE mpivars
      USE filefmt
      USE fft
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,n,ista:iend) :: u,v,w,s
      REAL(KIND=GP),INTENT(IN)                               :: fo, bv
      DOUBLE PRECISION, DIMENSION(n)                         :: gsh
      CHARACTER(len=*), INTENT(IN)                           :: nmb,spref
      INTEGER,INTENT(IN)                                     :: itype

!
      CALL havgcomp(gsh,u,v,w,s,fo,bv,itype) 
      IF (myrank.eq.0) THEN
         OPEN(1,file=trim(spref) // '.' // trim(nmb) // '.txt')
         WRITE(1,10) gsh
   10    FORMAT( E23.15 )
         CLOSE(1)
      ENDIF

      RETURN
      END SUBROUTINE havgwrite


!*****************************************************************
      SUBROUTINE tbouss(u, v, w, s, t, dt, fo, bv)
!-----------------------------------------------------------------
!
! Computes the volume-average and max/min of quantities computed in
! havgcomp, and output them to disk as a function of time.
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


      USE fprecision
      USE commtypes
      USE kes
      USE grid
      USE mpivars
      USE filefmt
      USE fft
      USE ali
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,n,ista:iend) :: u,v,w,s
      REAL(KIND=GP),INTENT(IN)                   :: dt
      REAL(KIND=GP),INTENT(IN)                   :: fo, bv
      DOUBLE PRECISION                           :: gxavg(11),gxmax(11),gxmin(11),tmp
      DOUBLE PRECISION,  DIMENSION(n)            :: gsh
      INTEGER, INTENT(IN)                        :: t
      INTEGER                                    :: i

      gxavg  = 0.0D0
      gxmax  = 0.0D0
      gxmin  = 1.0D0/tinyf

! Compute volume averages from horizontal averages, and
! also extrema from horiz. averages:
      DO i = 1, 11
        CALL havgcomp(gsh,u,v,w,s,fo,bv,i-1)
        gxavg(i) = sum   (gsh)/dble(n);
        gxmax(i) = maxval(gsh);
        gxmin(i) = minval(gsh);
      ENDDO
!
! Don't need reductions bec the results from havgcomp are globalized,
! so just do the outputs:

! NOTE:
!    col== 2 : hor. avg. of shear == <du_perp/dz>^2
!    col== 3 : hor. avg. of vert. temp. gradient == <d\theta/dz>
!    col== 4 : hor. avg. of correlation == <u_z d\theta/dz>
!    col== 5 : hor. avg. of hor. kinetic energy == <u^2 + v^2>
!    col== 6 : hor. avg. of vert. kinetic energy == <w^2>
!    col== 7 : hor. avg. of perp. helicity == <u_perp . curl u_perp = -u dv/dz + v du/dz>
!    col== 8 : hor. avg. of correlation: <\omega_z \theta>
!    col== 9 : hor. avg. of PV^2: <(fd\theta/dz - N \omega_z + omega.Grad\theta -fN)^2>
!    col== 10: hor. avg. of 'super-helicity': <\omega . curl \omega>
!    col== 11: hor. avg. of gradient Richardson no: !    N(N-d\theta/dz)/(du_\perp/dz)^2
!    col== 12: hor. avg. of correlation: <u_z \theta>
!
! Output quantities as a fcn of t:
      IF (myrank.eq.0) THEN
         OPEN(1,file='tboussavg.txt',position='append')
         WRITE(1,10) (t-1)*dt,gxavg (1),gxavg(2),gxavg(3),gxavg(4),gxavg (5), &
                              gxavg (6),gxavg(7),gxavg(8),gxavg(9),gxavg(10), &
                              gxavg(11)
         CLOSE(1)
      ENDIF
!
! Output max quantities as a fcn of t:
      IF (myrank.eq.0) THEN
         OPEN(1,file='tboussmax.txt',position='append')
         WRITE(1,10) (t-1)*dt,gxmax (1),gxmax(2),gxmax(3),gxmax(4),gxmax (5), &
                              gxmax (6),gxmax(7),gxmax(8),gxmax(9),gxmax(10), &
                              gxmax(11)
         CLOSE(1)
      ENDIF
!
! Output min quantities as a fcn of t:
      IF (myrank.eq.0) THEN
         OPEN(1,file='tboussmin.txt',position='append')
         WRITE(1,10) (t-1)*dt,gxmin(1),gxmin(2),gxmin(3),gxmin(4),gxmin (5), &
                              gxmin(6),gxmin(7),gxmin(8),gxmin(9),gxmin(10), &
                              gxmin(11)
         CLOSE(1)
      ENDIF

      RETURN
   10 FORMAT( E13.6,1x,11(E26.18,1x) )
      END SUBROUTINE tbouss


!*****************************************************************
      SUBROUTINE spectpv(u,v,w,s,nmb)
!-----------------------------------------------------------------
!
! Computes the potential vorticity power spectrum. 
! The pot'l vorticity is computed as:
!  PV = curl(v) . Grad(theta)
! output is written to a file by the first node.
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
      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(n/2+1)            :: Ek,Ektot
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,n,ista:iend) :: u,v,w,s
      COMPLEX(KIND=GP), DIMENSION(n,n,ista:iend)             :: c1,c2,c3,a
      REAL(KIND=GP)    :: tmp
      INTEGER :: i,j,k
      INTEGER :: kmn
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
      DO i = 1,n/2+1
         Ek(i) = 0.
      END DO
!
! Computes the power spectrum
!
      tmp = 1./real(n,kind=GP)**6
      IF (ista.eq.1) THEN
         DO j = 1,n
            DO k = 1,n
               kmn = int(abs(ka(k))+1)
               IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                  Ek(kmn) = Ek(kmn)+tmp*abs(a(k,j,1))**2
               ENDIF
            END DO
         END DO
         DO i = 2,iend
            DO j = 1,n
               DO k = 1,n
                  kmn = int(abs(ka(k))+1)
                  IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                     Ek(kmn) = Ek(kmn)+2*tmp*abs(a(k,j,i))**2
                  ENDIF
               END DO
            END DO
         END DO
      ELSE
         DO i = ista,iend
            DO j = 1,n
               DO k = 1,n
                  kmn = int(abs(ka(k))+1)
                  IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                     Ek(kmn) = Ek(kmn)+2*tmp*abs(a(k,j,i))**2
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
         OPEN(1,file='pvspectrum.' // nmb // '.txt')
         WRITE(1,20) Ektot
   20    FORMAT( E23.15 ) 
         CLOSE(1)
      ENDIF

      RETURN
      END SUBROUTINE spectpv

