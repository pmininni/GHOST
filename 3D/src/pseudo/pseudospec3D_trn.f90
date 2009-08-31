!=================================================================
! PSEUDOSPECTRAL subroutines
!
! Subroutines to compute shell-to-shell and triadic transfer
! in the HD, MHD, and Hall-MHD equations in 3D using a 
! pseudo-spectral method. Extensions for the anisotropic case 
! (e.g. in the rotating frame or with an imposed external 
! magnetic field) are also provided.
! You should use the FFTPLANS and MPIVARS modules (see the 
! file 'fftp_mod.f90') in each program that calls any of the 
! subroutines in this file. 
!
! NOTATION: index 'i' is 'x' 
!           index 'j' is 'y'
!           index 'k' is 'z'
! 
! 2005 Alexandros Alexakis and Pablo D. Mininni.
!      National Center for Atmospheric Research.
!=================================================================

!*****************************************************************
      SUBROUTINE shelltran(ax,ay,az,rbx,rby,rbz,cx,cy,cz,k1,q1,out)
!-----------------------------------------------------------------
!
! Computes the shell-to-shell energy transfer A_k1.[(B.grad)C_q1], 
! from the spherical shell at k1 to the spherical shell at q1. The 
! output is only valid in the first node.
!
! Parameters
!     ax : x-component of the array A in Fourier space
!     ay : y-component of the array A in Fourier space
!     az : z-component of the array A in Fourier space
!     rby: x-component of the array B in real space
!     rby: y-component of the array B in real space
!     rby: z-component of the array B in real space
!     cx : x-component of the array C in Fourier space
!     cy : y-component of the array C in Fourier space
!     cz : z-component of the array C in Fourier space
!     k1 : inner radius of the k-shell
!     q1 : inner radius of the q-shell
!     out: at the output contains the transfer
!
      USE fprecision
      USE commtypes
      USE kes
      USE fft
      USE grid
      USE mpivars
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,n,ista:iend) :: ax,ay,az
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,n,ista:iend) :: cx,cy,cz
      COMPLEX(KIND=GP), DIMENSION(n,n,ista:iend)             :: c1
      REAL(KIND=GP), INTENT(IN), DIMENSION(n,n,ksta:kend)    :: rbx,rby,rbz
      REAL(KIND=GP), DIMENSION(n,n,ksta:kend)                :: r1,r2
      REAL(KIND=GP), INTENT(OUT)   :: out
      REAL(KIND=GP)                :: tmp
      INTEGER, INTENT(IN) :: k1,q1
      INTEGER             :: i,j,k

      tmp = 0.

! Step 1 - s stands for filtered
! Computes As_x(k).(B_x.dx)Cs_x(q)
!
      CALL dshell(q1,cx,c1,1)
      CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      DO k = ksta,kend
         DO j = 1,n
            DO i = 1,n
               r1(i,j,k) = rbx(i,j,k)*r1(i,j,k)
            END DO
         END DO
      END DO
      CALL shell(k1,ax,c1)
      CALL fftp3d_complex_to_real(plancr,c1,r2,MPI_COMM_WORLD)
      DO k = ksta,kend
         DO j = 1,n
            DO i = 1,n
               tmp = tmp+r1(i,j,k)*r2(i,j,k)
            END DO
         END DO
      END DO

! Step 2
! Computes As_x(k).(B_y.dy)Cs_x(q)
!
      CALL dshell(q1,cx,c1,2)
      CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      DO k = ksta,kend
         DO j = 1,n
            DO i = 1,n
               r1(i,j,k) = rby(i,j,k)*r1(i,j,k)
            END DO
         END DO
      END DO
      CALL shell(k1,ax,c1)
      CALL fftp3d_complex_to_real(plancr,c1,r2,MPI_COMM_WORLD)
      DO k = ksta,kend
         DO j = 1,n
            DO i = 1,n
               tmp = tmp+r1(i,j,k)*r2(i,j,k)
            END DO
         END DO
      END DO

! Step 3
! Computes As_x(k).(B_z.dz)Cs_x(q)
!
      CALL dshell(q1,cx,c1,3)
      CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      DO k = ksta,kend
         DO j = 1,n
            DO i = 1,n
               r1(i,j,k) = rbz(i,j,k)*r1(i,j,k)
            END DO
         END DO
      END DO
      CALL shell(k1,ax,c1)
      CALL fftp3d_complex_to_real(plancr,c1,r2,MPI_COMM_WORLD)
      DO k = ksta,kend
         DO j = 1,n
            DO i = 1,n
               tmp = tmp+r1(i,j,k)*r2(i,j,k)
            END DO
         END DO
      END DO

! Step 4
! Computes As_y(k).(B_x.dx)Cs_y(q)
!
      CALL dshell(q1,cy,c1,1)
      CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      DO k = ksta,kend
         DO j = 1,n
            DO i = 1,n
               r1(i,j,k) = rbx(i,j,k)*r1(i,j,k)
            END DO
         END DO
      END DO
      CALL shell(k1,ay,c1)
      CALL fftp3d_complex_to_real(plancr,c1,r2,MPI_COMM_WORLD)
      DO k = ksta,kend
         DO j = 1,n
            DO i = 1,n
               tmp = tmp+r1(i,j,k)*r2(i,j,k)
            END DO
         END DO
      END DO

! Step 5
! Computes As_y(k).(B_y.dy)Cs_y(q)
!
      CALL dshell(q1,cy,c1,2)
      CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      DO k = ksta,kend
         DO j = 1,n
            DO i = 1,n
               r1(i,j,k) = rby(i,j,k)*r1(i,j,k)
            END DO
         END DO
      END DO
      CALL shell(k1,ay,c1)
      CALL fftp3d_complex_to_real(plancr,c1,r2,MPI_COMM_WORLD)
      DO k = ksta,kend
         DO j = 1,n
            DO i = 1,n
               tmp = tmp+r1(i,j,k)*r2(i,j,k)
            END DO
         END DO
      END DO

! Step 6
! Computes As_y(k).(B_z.dz)Cs_y(q)
!
      CALL dshell(q1,cy,c1,3)
      CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      DO k = ksta,kend
         DO j = 1,n
            DO i = 1,n
               r1(i,j,k) = rbz(i,j,k)*r1(i,j,k)
            END DO
         END DO
      END DO
      CALL shell(k1,ay,c1)
      CALL fftp3d_complex_to_real(plancr,c1,r2,MPI_COMM_WORLD)
      DO k = ksta,kend
         DO j = 1,n
            DO i = 1,n
               tmp = tmp+r1(i,j,k)*r2(i,j,k)
            END DO
         END DO
      END DO

! Step 7
! Computes As_z(k).(B_x.dx)Cs_z(q)
!
      CALL dshell(q1,cz,c1,1)
      CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      DO k = ksta,kend
         DO j = 1,n
            DO i = 1,n
               r1(i,j,k) = rbx(i,j,k)*r1(i,j,k)
            END DO
         END DO
      END DO
      CALL shell(k1,az,c1)
      CALL fftp3d_complex_to_real(plancr,c1,r2,MPI_COMM_WORLD)
      DO k = ksta,kend
         DO j = 1,n
            DO i = 1,n
               tmp = tmp+r1(i,j,k)*r2(i,j,k)
            END DO
         END DO
      END DO

! Step 8
! Computes As_z(k).(B_y.dy)Cs_z(q)
!
      CALL dshell(q1,cz,c1,2)
      CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      DO k = ksta,kend
         DO j = 1,n
            DO i = 1,n
               r1(i,j,k) = rby(i,j,k)*r1(i,j,k)
            END DO
         END DO
      END DO
      CALL shell(k1,az,c1)
      CALL fftp3d_complex_to_real(plancr,c1,r2,MPI_COMM_WORLD)
      DO k = ksta,kend
         DO j = 1,n
            DO i = 1,n
               tmp = tmp+r1(i,j,k)*r2(i,j,k)
            END DO
         END DO
      END DO

! Step 9
! Computes As_z(k).(B_z.dz)Cs_z(q)
!
      CALL dshell(q1,cz,c1,3)
      CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      DO k = ksta,kend
         DO j = 1,n
            DO i = 1,n
               r1(i,j,k) = rbz(i,j,k)*r1(i,j,k)
            END DO
         END DO
      END DO
      CALL shell(k1,az,c1)
      CALL fftp3d_complex_to_real(plancr,c1,r2,MPI_COMM_WORLD)
      DO k = ksta,kend
         DO j = 1,n
            DO i = 1,n
               tmp = tmp+r1(i,j,k)*r2(i,j,k)
            END DO
         END DO
      END DO

      tmp = tmp/real(n,kind=GP)**9
      CALL MPI_REDUCE(tmp,out,1,GC_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)

      RETURN
      END SUBROUTINE shelltran

!*****************************************************************
      SUBROUTINE paratran(ax,ay,az,rbx,rby,rbz,cx,cy,cz,k1,q1,out)
!-----------------------------------------------------------------
!
! Computes the shell-to-shell energy transfer A_k1.[(B.grad)C_q1], 
! from the shell k1 to the shell q1 in the direction parallel to 
! the preferred direction (rotation or uniform magnetic field). 
! As a result, the k-shells are planes with normal (0,0,k). The 
! output is only valid in the first node
!
! Parameters
!     ax : x-component of the array A in Fourier space
!     ay : y-component of the array A in Fourier space
!     az : z-component of the array A in Fourier space
!     rby: x-component of the array B in real space
!     rby: y-component of the array B in real space
!     rby: z-component of the array B in real space
!     cx : x-component of the array C in Fourier space
!     cy : y-component of the array C in Fourier space
!     cz : z-component of the array C in Fourier space
!     k1 : height of the k-shell
!     q1 : height of the q-shell
!     out: at the output contains the transfer
!
      USE fprecision
      USE commtypes
      USE kes
      USE fft
      USE grid
      USE mpivars
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,n,ista:iend) :: ax,ay,az
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,n,ista:iend) :: cx,cy,cz
      COMPLEX(KIND=GP), DIMENSION(n,n,ista:iend)             :: c1
      REAL(KIND=GP), INTENT(IN), DIMENSION(n,n,ksta:kend)    :: rbx,rby,rbz
      REAL(KIND=GP), DIMENSION(n,n,ksta:kend)                :: r1,r2
      REAL(KIND=GP), INTENT(OUT)   :: out
      REAL(KIND=GP)                :: tmp
      INTEGER, INTENT(IN) :: k1,q1
      INTEGER             :: i,j,k

      tmp = 0.

! Step 1 - s stands for filtered
! Computes As_x(k).(B_x.dx)Cs_x(q)
!
      CALL dshellpara(q1,cx,c1,1)
      CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      DO k = ksta,kend
         DO j = 1,n
            DO i = 1,n
               r1(i,j,k) = rbx(i,j,k)*r1(i,j,k)
            END DO
         END DO
      END DO
      CALL shellpara(k1,ax,c1)
      CALL fftp3d_complex_to_real(plancr,c1,r2,MPI_COMM_WORLD)
      DO k = ksta,kend
         DO j = 1,n
            DO i = 1,n
               tmp = tmp+r1(i,j,k)*r2(i,j,k)
            END DO
         END DO
      END DO

! Step 2
! Computes As_x(k).(B_y.dy)Cs_x(q)
!
      CALL dshellpara(q1,cx,c1,2)
      CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      DO k = ksta,kend
         DO j = 1,n
            DO i = 1,n
               r1(i,j,k) = rby(i,j,k)*r1(i,j,k)
            END DO
         END DO
      END DO
      CALL shellpara(k1,ax,c1)
      CALL fftp3d_complex_to_real(plancr,c1,r2,MPI_COMM_WORLD)
      DO k = ksta,kend
         DO j = 1,n
            DO i = 1,n
               tmp = tmp+r1(i,j,k)*r2(i,j,k)
            END DO
         END DO
      END DO

! Step 3
! Computes As_x(k).(B_z.dz)Cs_x(q)
!
      CALL dshellpara(q1,cx,c1,3)
      CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      DO k = ksta,kend
         DO j = 1,n
            DO i = 1,n
               r1(i,j,k) = rbz(i,j,k)*r1(i,j,k)
            END DO
         END DO
      END DO
      CALL shellpara(k1,ax,c1)
      CALL fftp3d_complex_to_real(plancr,c1,r2,MPI_COMM_WORLD)
      DO k = ksta,kend
         DO j = 1,n
            DO i = 1,n
               tmp = tmp+r1(i,j,k)*r2(i,j,k)
            END DO
         END DO
      END DO

! Step 4
! Computes As_y(k).(B_x.dx)Cs_y(q)
!
      CALL dshellpara(q1,cy,c1,1)
      CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      DO k = ksta,kend
         DO j = 1,n
            DO i = 1,n
               r1(i,j,k) = rbx(i,j,k)*r1(i,j,k)
            END DO
         END DO
      END DO
      CALL shellpara(k1,ay,c1)
      CALL fftp3d_complex_to_real(plancr,c1,r2,MPI_COMM_WORLD)
      DO k = ksta,kend
         DO j = 1,n
            DO i = 1,n
               tmp = tmp+r1(i,j,k)*r2(i,j,k)
            END DO
         END DO
      END DO

! Step 5
! Computes As_y(k).(B_y.dy)Cs_y(q)
!
      CALL dshellpara(q1,cy,c1,2)
      CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      DO k = ksta,kend
         DO j = 1,n
            DO i = 1,n
               r1(i,j,k) = rby(i,j,k)*r1(i,j,k)
            END DO
         END DO
      END DO
      CALL shellpara(k1,ay,c1)
      CALL fftp3d_complex_to_real(plancr,c1,r2,MPI_COMM_WORLD)
      DO k = ksta,kend
         DO j = 1,n
            DO i = 1,n
               tmp = tmp+r1(i,j,k)*r2(i,j,k)
            END DO
         END DO
      END DO

! Step 6
! Computes As_y(k).(B_z.dz)Cs_y(q)
!
      CALL dshellpara(q1,cy,c1,3)
      CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      DO k = ksta,kend
         DO j = 1,n
            DO i = 1,n
               r1(i,j,k) = rbz(i,j,k)*r1(i,j,k)
            END DO
         END DO
      END DO
      CALL shellpara(k1,ay,c1)
      CALL fftp3d_complex_to_real(plancr,c1,r2,MPI_COMM_WORLD)
      DO k = ksta,kend
         DO j = 1,n
            DO i = 1,n
               tmp = tmp+r1(i,j,k)*r2(i,j,k)
            END DO
         END DO
      END DO

! Step 7
! Computes As_z(k).(B_x.dx)Cs_z(q)
!
      CALL dshellpara(q1,cz,c1,1)
      CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      DO k = ksta,kend
         DO j = 1,n
            DO i = 1,n
               r1(i,j,k) = rbx(i,j,k)*r1(i,j,k)
            END DO
         END DO
      END DO
      CALL shellpara(k1,az,c1)
      CALL fftp3d_complex_to_real(plancr,c1,r2,MPI_COMM_WORLD)
      DO k = ksta,kend
         DO j = 1,n
            DO i = 1,n
               tmp = tmp+r1(i,j,k)*r2(i,j,k)
            END DO
         END DO
      END DO

! Step 8
! Computes As_z(k).(B_y.dy)Cs_z(q)
!
      CALL dshellpara(q1,cz,c1,2)
      CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      DO k = ksta,kend
         DO j = 1,n
            DO i = 1,n
               r1(i,j,k) = rby(i,j,k)*r1(i,j,k)
            END DO
         END DO
      END DO
      CALL shellpara(k1,az,c1)
      CALL fftp3d_complex_to_real(plancr,c1,r2,MPI_COMM_WORLD)
      DO k = ksta,kend
         DO j = 1,n
            DO i = 1,n
               tmp = tmp+r1(i,j,k)*r2(i,j,k)
            END DO
         END DO
      END DO

! Step 9
! Computes As_z(k).(B_z.dz)Cs_z(q)
!
      CALL dshellpara(q1,cz,c1,3)
      CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      DO k = ksta,kend
         DO j = 1,n
            DO i = 1,n
               r1(i,j,k) = rbz(i,j,k)*r1(i,j,k)
            END DO
         END DO
      END DO
      CALL shellpara(k1,az,c1)
      CALL fftp3d_complex_to_real(plancr,c1,r2,MPI_COMM_WORLD)
      DO k = ksta,kend
         DO j = 1,n
            DO i = 1,n
               tmp = tmp+r1(i,j,k)*r2(i,j,k)
            END DO
         END DO
      END DO

      tmp = tmp/real(n,kind=GP)**9
      CALL MPI_REDUCE(tmp,out,1,GC_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)

      RETURN
      END SUBROUTINE paratran

!*****************************************************************
      SUBROUTINE perptran(ax,ay,az,rbx,rby,rbz,cx,cy,cz,k1,q1,out)
!-----------------------------------------------------------------
!
! Computes the shell-to-shell energy transfer A_k1.[(B.grad)C_q1], 
! from the cylindrical shell k1 to the cylindrical shell q1 (i.e. 
! in the direction perpendicular to the preferred direction). The 
! output is only valid in the first node
!
! Parameters
!     ax : x-component of the array A in Fourier space
!     ay : y-component of the array A in Fourier space
!     az : z-component of the array A in Fourier space
!     rby: x-component of the array B in real space
!     rby: y-component of the array B in real space
!     rby: z-component of the array B in real space
!     cx : x-component of the array C in Fourier space
!     cy : y-component of the array C in Fourier space
!     cz : z-component of the array C in Fourier space
!     k1 : inner radius of the k-shell
!     q1 : inner radius of the q-shell
!     out: at the output contains the transfer
!
      USE fprecision
      USE commtypes
      USE kes
      USE fft
      USE grid
      USE mpivars
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,n,ista:iend) :: ax,ay,az
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,n,ista:iend) :: cx,cy,cz
      COMPLEX(KIND=GP), DIMENSION(n,n,ista:iend)             :: c1
      REAL(KIND=GP), INTENT(IN), DIMENSION(n,n,ksta:kend)    :: rbx,rby,rbz
      REAL(KIND=GP), DIMENSION(n,n,ksta:kend)                :: r1,r2
      REAL(KIND=GP),INTENT(OUT)    :: out
      REAL(KIND=GP)                :: tmp
      INTEGER, INTENT(IN) :: k1,q1
      INTEGER             :: i,j,k

      tmp = 0.

! Step 1 - s stands for filtered
! Computes As_x(k).(B_x.dx)Cs_x(q)
!
      CALL dshellperp(q1,cx,c1,1)
      CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      DO k = ksta,kend
         DO j = 1,n
            DO i = 1,n
               r1(i,j,k) = rbx(i,j,k)*r1(i,j,k)
            END DO
         END DO
      END DO
      CALL shellperp(k1,ax,c1)
      CALL fftp3d_complex_to_real(plancr,c1,r2,MPI_COMM_WORLD)
      DO k = ksta,kend
         DO j = 1,n
            DO i = 1,n
               tmp = tmp+r1(i,j,k)*r2(i,j,k)
            END DO
         END DO
      END DO

! Step 2
! Computes As_x(k).(B_y.dy)Cs_x(q)
!
      CALL dshellperp(q1,cx,c1,2)
      CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      DO k = ksta,kend
         DO j = 1,n
            DO i = 1,n
               r1(i,j,k) = rby(i,j,k)*r1(i,j,k)
            END DO
         END DO
      END DO
      CALL shellperp(k1,ax,c1)
      CALL fftp3d_complex_to_real(plancr,c1,r2,MPI_COMM_WORLD)
      DO k = ksta,kend
         DO j = 1,n
            DO i = 1,n
               tmp = tmp+r1(i,j,k)*r2(i,j,k)
            END DO
         END DO
      END DO

! Step 3
! Computes As_x(k).(B_z.dz)Cs_x(q)
!
      CALL dshellperp(q1,cx,c1,3)
      CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      DO k = ksta,kend
         DO j = 1,n
            DO i = 1,n
               r1(i,j,k) = rbz(i,j,k)*r1(i,j,k)
            END DO
         END DO
      END DO
      CALL shellperp(k1,ax,c1)
      CALL fftp3d_complex_to_real(plancr,c1,r2,MPI_COMM_WORLD)
      DO k = ksta,kend
         DO j = 1,n
            DO i = 1,n
               tmp = tmp+r1(i,j,k)*r2(i,j,k)
            END DO
         END DO
      END DO

! Step 4
! Computes As_y(k).(B_x.dx)Cs_y(q)
!
      CALL dshellperp(q1,cy,c1,1)
      CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      DO k = ksta,kend
         DO j = 1,n
            DO i = 1,n
               r1(i,j,k) = rbx(i,j,k)*r1(i,j,k)
            END DO
         END DO
      END DO
      CALL shellperp(k1,ay,c1)
      CALL fftp3d_complex_to_real(plancr,c1,r2,MPI_COMM_WORLD)
      DO k = ksta,kend
         DO j = 1,n
            DO i = 1,n
               tmp = tmp+r1(i,j,k)*r2(i,j,k)
            END DO
         END DO
      END DO

! Step 5
! Computes As_y(k).(B_y.dy)Cs_y(q)
!
      CALL dshellperp(q1,cy,c1,2)
      CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      DO k = ksta,kend
         DO j = 1,n
            DO i = 1,n
               r1(i,j,k) = rby(i,j,k)*r1(i,j,k)
            END DO
         END DO
      END DO
      CALL shellperp(k1,ay,c1)
      CALL fftp3d_complex_to_real(plancr,c1,r2,MPI_COMM_WORLD)
      DO k = ksta,kend
         DO j = 1,n
            DO i = 1,n
               tmp = tmp+r1(i,j,k)*r2(i,j,k)
            END DO
         END DO
      END DO

! Step 6
! Computes As_y(k).(B_z.dz)Cs_y(q)
!
      CALL dshellperp(q1,cy,c1,3)
      CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      DO k = ksta,kend
         DO j = 1,n
            DO i = 1,n
               r1(i,j,k) = rbz(i,j,k)*r1(i,j,k)
            END DO
         END DO
      END DO
      CALL shellperp(k1,ay,c1)
      CALL fftp3d_complex_to_real(plancr,c1,r2,MPI_COMM_WORLD)
      DO k = ksta,kend
         DO j = 1,n
            DO i = 1,n
               tmp = tmp+r1(i,j,k)*r2(i,j,k)
            END DO
         END DO
      END DO

! Step 7
! Computes As_z(k).(B_x.dx)Cs_z(q)
!
      CALL dshellperp(q1,cz,c1,1)
      CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      DO k = ksta,kend
         DO j = 1,n
            DO i = 1,n
               r1(i,j,k) = rbx(i,j,k)*r1(i,j,k)
            END DO
         END DO
      END DO
      CALL shellperp(k1,az,c1)
      CALL fftp3d_complex_to_real(plancr,c1,r2,MPI_COMM_WORLD)
      DO k = ksta,kend
         DO j = 1,n
            DO i = 1,n
               tmp = tmp+r1(i,j,k)*r2(i,j,k)
            END DO
         END DO
      END DO

! Step 8
! Computes As_z(k).(B_y.dy)Cs_z(q)
!
      CALL dshellperp(q1,cz,c1,2)
      CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      DO k = ksta,kend
         DO j = 1,n
            DO i = 1,n
               r1(i,j,k) = rby(i,j,k)*r1(i,j,k)
            END DO
         END DO
      END DO
      CALL shellperp(k1,az,c1)
      CALL fftp3d_complex_to_real(plancr,c1,r2,MPI_COMM_WORLD)
      DO k = ksta,kend
         DO j = 1,n
            DO i = 1,n
               tmp = tmp+r1(i,j,k)*r2(i,j,k)
            END DO
         END DO
      END DO

! Step 9
! Computes As_z(k).(B_z.dz)Cs_z(q)
!
      CALL dshellperp(q1,cz,c1,3)
      CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      DO k = ksta,kend
         DO j = 1,n
            DO i = 1,n
               r1(i,j,k) = rbz(i,j,k)*r1(i,j,k)
            END DO
         END DO
      END DO
      CALL shellperp(k1,az,c1)
      CALL fftp3d_complex_to_real(plancr,c1,r2,MPI_COMM_WORLD)
      DO k = ksta,kend
         DO j = 1,n
            DO i = 1,n
               tmp = tmp+r1(i,j,k)*r2(i,j,k)
            END DO
         END DO
      END DO

      tmp = tmp/real(n,kind=GP)**9
      CALL MPI_REDUCE(tmp,out,1,GC_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)

      RETURN
      END SUBROUTINE perptran

!*****************************************************************
      SUBROUTINE crosstran(ax,ay,az,rbx,rby,rbz,k1,q1,out)
!-----------------------------------------------------------------
!
! Computes the energy transfer A_k1.(B x A_q1) from the 
! spherical shell at k1 to the spherical shell at q1. The 
! output is only valid in the first node
!
! Parameters
!     ax : x-component of the array A in Fourier space
!     ay : y-component of the array A in Fourier space
!     az : z-component of the array A in Fourier space
!     rbx: x-component of the array B in real space
!     rby: y-component of the array B in real space
!     rbz: z-component of the array B in real space
!     k1 : inner radius of the k-shell
!     q1 : inner radius of the q-shell
!     out: at the output contains the transfer
!
      USE fprecision
      USE commtypes
      USE kes
      USE fft
      USE grid
      USE mpivars
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,n,ista:iend) :: ax,ay,az
      COMPLEX(KIND=GP), DIMENSION(n,n,ista:iend)             :: c1
      REAL(KIND=GP), INTENT(IN), DIMENSION(n,n,ksta:kend)    :: rbx,rby,rbz
      REAL(KIND=GP), DIMENSION(n,n,ksta:kend)                :: r1,r2,r3,r4
      REAL(KIND=GP),INTENT(OUT)    :: out
      REAL(KIND=GP)                :: tmp
      INTEGER, INTENT(IN) :: k1,q1
      INTEGER             :: i,j,k

      tmp = 0.

! Step 1 - s stands for filtered
! Computes As_x(k).[B x curl(As(q))_x]
!
      CALL shell(k1,ax,c1)
      CALL fftp3d_complex_to_real(plancr,c1,r4,MPI_COMM_WORLD)
      CALL shell(q1,ay,c1)
      CALL fftp3d_complex_to_real(plancr,c1,r2,MPI_COMM_WORLD)
      CALL shell(q1,az,c1)
      CALL fftp3d_complex_to_real(plancr,c1,r3,MPI_COMM_WORLD)
      DO k = ksta,kend
         DO j = 1,n
            DO i = 1,n
               tmp = tmp+r4(i,j,k)*(rby(i,j,k)*r3(i,j,k)      &
                     -rbz(i,j,k)*r2(i,j,k))
            END DO
         END DO
      END DO

! Step 2 - s stands for filtered
! Computes As_y(k).[B x curl(As(q))_y]
!
      CALL shell(k1,ay,c1)
      CALL fftp3d_complex_to_real(plancr,c1,r4,MPI_COMM_WORLD)
      CALL shell(q1,ax,c1)
      CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      DO k = ksta,kend
         DO j = 1,n
            DO i = 1,n
               tmp = tmp+r4(i,j,k)*(rbz(i,j,k)*r1(i,j,k)      &
                     -rbx(i,j,k)*r3(i,j,k))
            END DO
         END DO
      END DO

! Step 2 - s stands for filtered
! Computes As_z(k).[B x curl(As(q))_z]
!
      CALL shell(k1,az,c1)
      CALL fftp3d_complex_to_real(plancr,c1,r4,MPI_COMM_WORLD)
      DO k = ksta,kend
         DO j = 1,n
            DO i = 1,n
               tmp = tmp+r4(i,j,k)*(rbx(i,j,k)*r2(i,j,k)      &
                     -rby(i,j,k)*r1(i,j,k))
            END DO
         END DO
      END DO

      tmp = tmp/real(n,kind=GP)**9
      CALL MPI_REDUCE(tmp,out,1,GC_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)

      RETURN
      END SUBROUTINE crosstran

!*****************************************************************
      SUBROUTINE triagtran(ax,ay,az,rxx,rxy,rxz,ryx,ryy,ryz, &
                 rzx,rzy,rzz,k1,p1,out)
!-----------------------------------------------------------------
!
! Computes the shell-to-shell energy transfer 
! A_k1.[(A_p1.grad)B_q1] with q1 fixed, from the spherical 
! shell k1 to the spherical shell at q1 mediated by the modes 
! in the spherical shell p1. The output is only valid in the 
! first node
!
! Parameters
!     ax : x-component of the array A in Fourier space
!     ay : y-component of the array A in Fourier space
!     az : z-component of the array A in Fourier space
!     rxx: x-component of the x-derivative of B in real space
!     rxy: x-component of the y-derivative of B in real space
!     rxz: x-component of the z-derivative of B in real space
!     ryx: y-component of the x-derivative of B in real space
!     ryy: y-component of the y-derivative of B in real space
!     ryz: y-component of the z-derivative of B in real space
!     rzx: z-component of the x-derivative of B in real space
!     rzy: z-component of the y-derivative of B in real space
!     rzz: z-component of the z-derivative of B in real space
!     k1 : inner radius of the k-shell
!     p1 : inner radius of the p-shell
!     out: at the output contains the transfer
!
      USE fprecision
      USE commtypes
      USE kes
      USE fft
      USE grid
      USE mpivars
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,n,ista:iend) :: ax,ay,az
      COMPLEX(KIND=GP), DIMENSION(n,n,ista:iend)             :: c1
      REAL(KIND=GP), INTENT(IN), DIMENSION(n,n,ksta:kend)    :: rxx,rxy,rxz
      REAL(KIND=GP), INTENT(IN), DIMENSION(n,n,ksta:kend)    :: ryx,ryy,ryz
      REAL(KIND=GP), INTENT(IN), DIMENSION(n,n,ksta:kend)    :: rzx,rzy,rzz
      REAL(KIND=GP), DIMENSION(n,n,ksta:kend)                :: r1,r2
      REAL(KIND=GP), INTENT(OUT)   :: out
      REAL(KIND=GP)                :: tmp
      INTEGER, INTENT(IN) :: k1,p1
      INTEGER             :: i,j,k

      tmp = 0.

! Step 1 - s stands for filtered
! Computes As_x(k).(As_x(p).dx)Bs_x(q)
!
      CALL shell(p1,ax,c1)
      CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      DO k = ksta,kend
         DO j = 1,n
            DO i = 1,n
               r1(i,j,k) = rxx(i,j,k)*r1(i,j,k)
            END DO
         END DO
      END DO
      CALL shell(k1,ax,c1)
      CALL fftp3d_complex_to_real(plancr,c1,r2,MPI_COMM_WORLD)
      DO k = ksta,kend
         DO j = 1,n
            DO i = 1,n
               tmp = tmp+r1(i,j,k)*r2(i,j,k)
            END DO
         END DO
      END DO

! Step 2
! Computes As_x(k).(As_y(p).dy)Bs_x(q)
!
      CALL shell(p1,ay,c1)
      CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      DO k = ksta,kend
         DO j = 1,n
            DO i = 1,n
               r1(i,j,k) = rxy(i,j,k)*r1(i,j,k)
            END DO
         END DO
      END DO
      CALL shell(k1,ax,c1)
      CALL fftp3d_complex_to_real(plancr,c1,r2,MPI_COMM_WORLD)
      DO k = ksta,kend
         DO j = 1,n
            DO i = 1,n
               tmp = tmp+r1(i,j,k)*r2(i,j,k)
            END DO
         END DO
      END DO

! Step 3
! Computes As_x(k).(As_z(p).dz)Bs_x(q)
!
      CALL shell(p1,az,c1)
      CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      DO k = ksta,kend
         DO j = 1,n
            DO i = 1,n
               r1(i,j,k) = rxz(i,j,k)*r1(i,j,k)
            END DO
         END DO
      END DO
      CALL shell(k1,ax,c1)
      CALL fftp3d_complex_to_real(plancr,c1,r2,MPI_COMM_WORLD)
      DO k = ksta,kend
         DO j = 1,n
            DO i = 1,n
               tmp = tmp+r1(i,j,k)*r2(i,j,k)
            END DO
         END DO
      END DO

! Step 4
! Computes As_y(k).(As_x(p).dx)Bs_y(q)
!
      CALL shell(p1,ax,c1)
      CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      DO k = ksta,kend
         DO j = 1,n
            DO i = 1,n
               r1(i,j,k) = ryx(i,j,k)*r1(i,j,k)
            END DO
         END DO
      END DO
      CALL shell(k1,ay,c1)
      CALL fftp3d_complex_to_real(plancr,c1,r2,MPI_COMM_WORLD)
      DO k = ksta,kend
         DO j = 1,n
            DO i = 1,n
               tmp = tmp+r1(i,j,k)*r2(i,j,k)
            END DO
         END DO
      END DO

! Step 5
! Computes As_y(k).(As_y(p).dy)Bs_y(q)
!
      CALL shell(p1,ay,c1)
      CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      DO k = ksta,kend
         DO j = 1,n
            DO i = 1,n
               r1(i,j,k) = ryy(i,j,k)*r1(i,j,k)
            END DO
         END DO
      END DO
      CALL shell(k1,ay,c1)
      CALL fftp3d_complex_to_real(plancr,c1,r2,MPI_COMM_WORLD)
      DO k = ksta,kend
         DO j = 1,n
            DO i = 1,n
               tmp = tmp+r1(i,j,k)*r2(i,j,k)
            END DO
         END DO
      END DO

! Step 6
! Computes As_y(k).(As_z(p).dz)Bs_y(q)
!
      CALL shell(p1,az,c1)
      CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      DO k = ksta,kend
         DO j = 1,n
            DO i = 1,n
               r1(i,j,k) = ryz(i,j,k)*r1(i,j,k)
            END DO
         END DO
      END DO
      CALL shell(k1,ay,c1)
      CALL fftp3d_complex_to_real(plancr,c1,r2,MPI_COMM_WORLD)
      DO k = ksta,kend
         DO j = 1,n
            DO i = 1,n
               tmp = tmp+r1(i,j,k)*r2(i,j,k)
            END DO
         END DO
      END DO

! Step 7
! Computes As_z(k).(As_x(p).dx)Bs_z(q)
!
      CALL shell(p1,ax,c1)
      CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      DO k = ksta,kend
         DO j = 1,n
            DO i = 1,n
               r1(i,j,k) = rzx(i,j,k)*r1(i,j,k)
            END DO
         END DO
      END DO
      CALL shell(k1,az,c1)
      CALL fftp3d_complex_to_real(plancr,c1,r2,MPI_COMM_WORLD)
      DO k = ksta,kend
         DO j = 1,n
            DO i = 1,n
               tmp = tmp+r1(i,j,k)*r2(i,j,k)
            END DO
         END DO
      END DO

! Step 8
! Computes As_z(k).(As_y(p).dy)Bs_z(q)
!
      CALL shell(p1,ay,c1)
      CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      DO k = ksta,kend
         DO j = 1,n
            DO i = 1,n
               r1(i,j,k) = rzy(i,j,k)*r1(i,j,k)
            END DO
         END DO
      END DO
      CALL shell(k1,az,c1)
      CALL fftp3d_complex_to_real(plancr,c1,r2,MPI_COMM_WORLD)
      DO k = ksta,kend
         DO j = 1,n
            DO i = 1,n
               tmp = tmp+r1(i,j,k)*r2(i,j,k)
            END DO
         END DO
      END DO

! Step 9
! Computes As_z(k).(As_z(p).dz)Bs_z(q)
!
      CALL shell(p1,az,c1)
      CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      DO k = ksta,kend
         DO j = 1,n
            DO i = 1,n
               r1(i,j,k) = rzz(i,j,k)*r1(i,j,k)
            END DO
         END DO
      END DO
      CALL shell(k1,az,c1)
      CALL fftp3d_complex_to_real(plancr,c1,r2,MPI_COMM_WORLD)
      DO k = ksta,kend
         DO j = 1,n
            DO i = 1,n
               tmp = tmp+r1(i,j,k)*r2(i,j,k)
            END DO
         END DO
      END DO

      tmp = tmp/real(n,kind=GP)**9
      tmp = tmp/real(n,kind=GP)**3
      CALL MPI_REDUCE(tmp,out,1,GC_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)

      RETURN
      END SUBROUTINE triagtran

!*****************************************************************
      SUBROUTINE dshell(k1,a,b,dir)
!-----------------------------------------------------------------
!
! Computes the filtered derivative at the spherical shell k1
!
! Parameters
!     k1: inner radius of the k-shell
!     a : field component in Fourier space [in]
!     b : fieltered derivative of a [out]
!     dir: =1 derivative in the x-direction
!          =2 derivative in the y-direction
!          =3 derivative in the z-direction
!
      USE fprecision
      USE kes
      USE var
      USE grid
      USE mpivars
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,n,ista:iend)  :: a
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(n,n,ista:iend) :: b
      INTEGER, INTENT(IN) :: k1
      INTEGER, INTENT(IN) :: dir
      INTEGER             :: i,j,k

      DO i = ista,iend
         DO j = 1,n
            DO k = 1,n

            IF ((ka2(k,j,i).ge.k1**2).and.(ka2(k,j,i).lt.(k1+1)**2)) THEN
               IF (dir.eq.1) THEN
                  b(k,j,i) = im*ka(i)*a(k,j,i)
               ELSE IF (dir.eq.2) THEN
                  b(k,j,i) = im*ka(j)*a(k,j,i)
               ELSE
                  b(k,j,i) = im*ka(k)*a(k,j,i)
               ENDIF
            ELSE
               b(k,j,i) = 0.
            ENDIF

            END DO
         END DO
      END DO

      RETURN
      END SUBROUTINE dshell

!*****************************************************************
      SUBROUTINE dshellpara(k1,a,b,dir)
!-----------------------------------------------------------------
!
! Computes the filtered derivative at the shell k1. The k-shells 
! are planes with normal (0,0,k)
!
! Parameters
!     k1: height of the k-shell
!     a : field component in Fourier space [in]
!     b : fieltered derivative of a [out]
!     dir: =1 derivative in the x-direction
!          =2 derivative in the y-direction
!          =3 derivative in the z-direction
!
      USE fprecision
      USE kes
      USE var
      USE grid
      USE mpivars
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,n,ista:iend)  :: a
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(n,n,ista:iend) :: b
      INTEGER, INTENT(IN) :: k1
      INTEGER, INTENT(IN) :: dir
      INTEGER             :: i,j,k

      DO i = ista,iend
         DO j = 1,n
            DO k = 1,n

            IF ((ka(k)**2.ge.k1**2).and.(ka(k)**2.lt.(k1+1)**2)) THEN
               IF (dir.eq.1) THEN
                  b(k,j,i) = im*ka(i)*a(k,j,i)
               ELSE IF (dir.eq.2) THEN
                  b(k,j,i) = im*ka(j)*a(k,j,i)
               ELSE
                  b(k,j,i) = im*ka(k)*a(k,j,i)
               ENDIF
            ELSE
               b(k,j,i) = 0.
            ENDIF

            END DO
         END DO
      END DO

      RETURN
      END SUBROUTINE dshellpara

!*****************************************************************
      SUBROUTINE dshellperp(k1,a,b,dir)
!-----------------------------------------------------------------
!
! Computes the filtered derivative at the cylindrical shell k1
!
! Parameters
!     k1: inner radius of the k-shell
!     a : field component in Fourier space [in]
!     b : fieltered derivative of a [out]
!     dir: =1 derivative in the x-direction
!          =2 derivative in the y-direction
!          =3 derivative in the z-direction
!
      USE fprecision
      USE kes
      USE var
      USE grid
      USE mpivars
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,n,ista:iend)  :: a
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(n,n,ista:iend) :: b
      REAL(KIND=GP)                :: kperp2
      INTEGER, INTENT(IN) :: k1
      INTEGER, INTENT(IN) :: dir
      INTEGER             :: i,j,k

      DO i = ista,iend
         DO j = 1,n
            kperp2 = ka(i)**2+ka(j)**2
            DO k = 1,n

            IF ((kperp2.ge.k1**2).and.(kperp2.lt.(k1+1)**2)) THEN
               IF (dir.eq.1) THEN
                  b(k,j,i) = im*ka(i)*a(k,j,i)
               ELSE IF (dir.eq.2) THEN
                  b(k,j,i) = im*ka(j)*a(k,j,i)
               ELSE
                  b(k,j,i) = im*ka(k)*a(k,j,i)
               ENDIF
            ELSE
               b(k,j,i) = 0.
            ENDIF

            END DO
         END DO
      END DO

      RETURN
      END SUBROUTINE dshellperp

!*****************************************************************
      SUBROUTINE shell(k1,a,b)
!-----------------------------------------------------------------
!
! Computes the filtered field at the spherical shell k1
!
! Parameters
!     k1: inner radius of the k-shell
!     a : field component in Fourier space [in]
!     b : fieltered field [out]
!
      USE fprecision
      USE kes
      USE var
      USE grid
      USE mpivars
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,n,ista:iend)  :: a
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(n,n,ista:iend) :: b
      INTEGER, INTENT(IN) :: k1
      INTEGER             :: i,j,k

      DO i = ista,iend
         DO j = 1,n
            DO k = 1,n

            IF ((ka2(k,j,i).ge.k1**2).and.(ka2(k,j,i).lt.(k1+1)**2)) THEN
               b(k,j,i) = a(k,j,i)
            ELSE
               b(k,j,i) = 0.
            ENDIF

            END DO
         END DO
      END DO

      RETURN
      END SUBROUTINE shell

!*****************************************************************
      SUBROUTINE shellpara(k1,a,b)
!-----------------------------------------------------------------
!
! Computes the filtered field at the shell k1. The k-shells are 
! planes with normal (0,0,k)
!
! Parameters
!     k1: height of the k-shell
!     a : field component in Fourier space [in]
!     b : fieltered field [out]
!
      USE fprecision
      USE kes
      USE var
      USE grid
      USE mpivars
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,n,ista:iend)  :: a
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(n,n,ista:iend) :: b
      INTEGER, INTENT(IN) :: k1
      INTEGER             :: i,j,k

      DO i = ista,iend
         DO j = 1,n
            DO k = 1,n

            IF ((ka(k)**2.ge.k1**2).and.(ka(k)**2.lt.(k1+1)**2)) THEN
               b(k,j,i) = a(k,j,i)
            ELSE
               b(k,j,i) = 0.
            ENDIF

            END DO
         END DO
      END DO

      RETURN
      END SUBROUTINE shellpara

!*****************************************************************
      SUBROUTINE shellperp(k1,a,b)
!-----------------------------------------------------------------
!
! Computes the filtered field at the cylindrical shell k1
!
! Parameters
!     k1: inner radius of the k-shell
!     a : field component in Fourier space [in]
!     b : fieltered field [out]
!
      USE fprecision
      USE kes
      USE var
      USE grid
      USE mpivars
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,n,ista:iend)  :: a
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(n,n,ista:iend) :: b
      REAL(KIND=GP)                :: kperp2
      INTEGER, INTENT(IN) :: k1
      INTEGER             :: i,j,k

      DO i = ista,iend
         DO j = 1,n
            kperp2 = ka(i)**2+ka(j)**2
            DO k = 1,n

            IF ((kperp2.ge.k1**2).and.(kperp2.lt.(k1+1)**2)) THEN
               b(k,j,i) = a(k,j,i)
            ELSE
               b(k,j,i) = 0.
            ENDIF

            END DO
         END DO
      END DO

      RETURN
      END SUBROUTINE shellperp

