!=================================================================
! PSEUDOSPECTRAL subroutines
!
! Extra subroutines to compute spatial derivatives and 
! nonlinear terms in Hybrid-PIC solvers using a
! pseudo-spectral method. You should use the FFTPLANS and 
! MPIVARS modules (see the file 'fftp_mod.f90') in each 
! program that calls any of the subroutines in this file. 
!
! NOTATION: index 'i' is 'x' 
!           index 'j' is 'y'
!           index 'k' is 'z'
!
! 2023 Facundo Pugliese
!      Department of Physics, 
!      Facultad de Ciencias Exactas y Naturales.
!      Universidad de Buenos Aires.
!      e-mail: fpugliese@df.uba.ar 
!=================================================================

!*****************************************************************
      SUBROUTINE poisson_elecstat(rhoc,kde2,phi)
!-----------------------------------------------------------------
! 
! Solves poisson problem lap(phi) = -rhoc + kde2*phi wherw kde2 is
! the square inverse Debye length.
!
! Parameters
!     rhoc : input complex density
!     kde2 : square inverse Debye length
!     phi  : output complex potential
!
      USE fprecision
      USE kes
      USE grid
      USE mpivars
!$    USE threads

      COMPLEX(KIND=GP), INTENT (IN), DIMENSION(nz,ny,ista:iend) :: rhoc
      REAL(KIND=GP)   , INTENT (IN)                             :: kde2
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(nz,ny,ista:iend) :: phi

      IF (myrank.eq.0) THEN
         phi(1,1,1) = 0.0_GP
         DO k = 2,nz
            phi(k,1,1) = rhoc(k,1,1)/(kk2(k,1,1)+kde2)
         END DO
!$omp parallel do if (iend-ista.ge.nth) private (k)
         DO j = 2,ny
            DO k = 1,nz
               phi(k,j,1) = rhoc(k,j,1)/(kk2(k,j,1)+kde2)
            END DO
         END DO
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = 2,iend
!$omp parallel do if (iend-ista.ge.nth) private (k)
            DO j = 1,ny
               DO k = 1,nz
                  phi(k,j,i) = rhoc(k,j,i)/(kk2(k,j,i)+kde2)
               END DO
            END DO
         END DO
      ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.ge.nth) private (k)
            DO j = 1,ny
               DO k = 1,nz
                  phi(k,j,i) = rhoc(k,j,i)/(kk2(k,j,i)+kde2)
               END DO
            END DO
         END DO
      END IF

      RETURN
      END SUBROUTINE poisson_elecstat

!*****************************************************************
      SUBROUTINE hpiccheck(Tem,ux,uy,uz,ni,ax,ay,az,gamm,beta,t,dt)
!-----------------------------------------------------------------
!
! Consistency check for the conservation of the total
! energy of electrostatic hybrid-pic solver.
!
! Output files contain:
! 'energy.txt':   time, <u^2>, <(v-u)^2>, <B^2>, <Pe*beta/(gamma-1)>
!
! Parameters
!       Tem: ion kinetic energy density in real space
!  ux,uy,uz: ion bulk velocity components in real space
!       ni : ion density in real space
!  ax,ay,az: vector potential components in Fourier space
!      gamm: barotropic exponent (fluid electrons)
!      beta: electronic plasma beta
!       t  : number of time steps made
!       dt : time step
!
      USE fprecision
      USE commtypes
      USE fft
      USE grid
      USE mpivars
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP),INTENT(IN),DIMENSION(nz,ny,ista:iend) :: ax,ay,az
      REAL(KIND=GP)   ,INTENT(IN),DIMENSION(ksta:kend,ny,nz) :: Tem,ni
      REAL(KIND=GP)   ,INTENT(IN),DIMENSION(ksta:kend,ny,nz) :: ux,uy,uz
      REAL(KIND=GP), INTENT(IN)    :: dt,gamm,beta
      INTEGER      , INTENT(IN)    :: t
      COMPLEX(KIND=GP),DIMENSION(nz,ny,ista:iend) :: c1,c2,c3,c4,c5,c6
      DOUBLE PRECISION    :: et,etloc,ei,eiloc,em,ek,ekloc,div
      INTEGER             :: i,j,k
      REAL(KIND=GP)       :: rmp,tmp

      CALL energy(ax,ay,az,em,0) 
      
      
      rmp = 1/(real(nx,KIND=GP)*real(ny,KIND=GP)*real(nz,KIND=GP))
      ekloc = 0.0D0
      eiloc = 0.0D0
      IF (gamm.EQ.1) THEN
      DO k = ksta,kend
!$omp parallel do if (kend-2.ge.nth) private (j,i) reduction(+:ekloc)
         DO j = 1,ny
!$omp parallel do if (kend-2.ge.nth) private (i) reduction(+:ekloc)
            DO i = 1,nx
               etloc = etloc + Tem(k,j,i)*rmp
               ekloc = ekloc + (ux(k,j,i)**2 +uy(k,j,i)**2   &
                               +uz(k,j,i)**2)*rmp*ni(k,j,i)
               eiloc = eiloc + ni(k,j,i)*LOG(ni(k,j,i))*rmp
            ENDDO
         ENDDO
      ENDDO
      tmp = beta
      ELSE
      DO k = ksta,kend
!$omp parallel do if (kend-2.ge.nth) private (j,i) reduction(+:ekloc)
         DO j = 1,ny
!$omp parallel do if (kend-2.ge.nth) private (i) reduction(+:ekloc)
            DO i = 1,nx
               etloc = etloc + Tem(k,j,i)*rmp
               ekloc = ekloc + (ux(k,j,i)**2 +uy(k,j,i)**2   &
                               +uz(k,j,i)**2)*rmp
               eiloc = eiloc + (ni(k,j,i)**gamm)*rmp 
            ENDDO
         ENDDO
      ENDDO
      tmp = beta/(gamm-1)
      END IF

      CALL MPI_REDUCE(ekloc,ek,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(etloc,et,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(eiloc,ei,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
!
! Creates external files to store the results
!
      IF (myrank.eq.0) THEN
         ei = ei*tmp
         et = et - ek
         OPEN(1,file='energy.txt',position='append')
         WRITE(1,10) (t-1)*dt,ek,et,em,ei
         CLOSE(1)
         et = ek + et + em + ei
      ENDIF

      CALL fftp3d_real_to_complex(planrc,ux,c1,MPI_COMM_WORLD)
      CALL fftp3d_real_to_complex(planrc,uy,c2,MPI_COMM_WORLD)
      CALL fftp3d_real_to_complex(planrc,uz,c3,MPI_COMM_WORLD)

      CALL maxabs(c1,c2,c3,rmp,0)  ! maximum vorticity
      CALL maxabs(ax,ay,az,tmp,1)  ! maximum current density

      CALL derivk3(c1,c4,1)
      CALL derivk3(c2,c5,2)
      CALL derivk3(c3,c6,3)
      CALL energy(c4,c5,c6,div,1)  ! <div(u)^2>

      CALL rotor3(c2,c3,c4,1)
      CALL rotor3(c1,c3,c5,2)
      CALL rotor3(c1,c2,c6,3)
      CALL energy(c4,c5,c6,ek,1)   ! <w^2>

      CALL laplak3(ax,c1)
      CALL laplak3(ay,c2) 
      CALL laplak3(az,c3) 
      CALL energy(c1,c2,c3,em,1)   ! <j^2>

      IF (myrank.eq.0) THEN
         OPEN(1,file='balance.txt',position='append')
         WRITE(1,10) (t-1)*dt,et,ek,em,div
   10    FORMAT( E13.6,E22.14,E22.14,E22.14,E22.14 )
         CLOSE(1)
         OPEN(1,file='maximum.txt',position='append')
         WRITE(1,FMT='(E13.6,E13.6,E13.6)') (t-1)*dt,rmp,tmp
         CLOSE(1)
      ENDIF

      RETURN
      END SUBROUTINE hpiccheck

!*****************************************************************
      SUBROUTINE ehpiccheck(Tem, rho, phi, t, dt)
!-----------------------------------------------------------------
!
! Consistency check for the conservation of the total
! energy of electrostatic hybrid-pic solver.
!
! Output files contain:
! 'energy.txt':   time, <v^2>, <E^2>
!
! Parameters
!     Tem: kinetic energy (of particles) in real space
!     rho: charge density in Fourier space
!     phi: electric potential in Fourier space
!     t  : number of time steps made
!     dt : time step
!
      USE fprecision
      USE commtypes
      USE grid
      USE mpivars
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP),INTENT(IN),DIMENSION(nz,ny,ista:iend) :: rho,phi
      REAL(KIND=GP)   ,INTENT(IN),DIMENSION(ksta:kend,ny,nz) :: Tem
      REAL(KIND=GP), INTENT(IN)    :: dt
      INTEGER      , INTENT(IN)    :: t
      DOUBLE PRECISION    :: ek,ekloc,ep
      INTEGER             :: i,j,k
      REAL(KIND=GP)       :: rmp

      CALL product(rho,phi,ep)
      
      rmp = 1/(real(nx,KIND=GP)*real(ny,KIND=GP)*real(nz,KIND=GP))
      ekloc = 0.0D0
      DO k = ksta,kend
!$omp parallel do if (kend-2.ge.nth) private (j,i) reduction(+:ekloc)
         DO j = 1,ny
!$omp parallel do if (kend-2.ge.nth) private (i) reduction(+:ekloc)
            DO i = 1,nx
               ekloc = ekloc + Tem(k,j,i)*rmp
            ENDDO
         ENDDO
      ENDDO

      CALL MPI_REDUCE(ekloc,ek,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
!
! Creates external files to store the results
!
      IF (myrank.eq.0) THEN
         OPEN(1,file='energy.txt',position='append')
         WRITE(1,10) (t-1)*dt,ek,ep
   10    FORMAT( E13.6,E22.14,E22.14)
         CLOSE(1)
      ENDIF

      RETURN
      END SUBROUTINE ehpiccheck

!*****************************************************************
      SUBROUTINE dealias(a)
!-----------------------------------------------------------------
!
! Applies 2/3 rule to dealias field a in Fourier space
!
! Parameters
!   a  : field in Fourier space
!

      USE fprecision
      USE kes
      USE grid
      USE mpivars
      USE ali
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP),INTENT(INOUT),DIMENSION(nz,ny,ista:iend):: a
      INTEGER                                                  :: i,j,k

!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.ge.nth) private (k)
         DO j = 1,ny
            DO k = 1,nz
               IF (kn2(k,j,i).GT.kmax) THEN
                  a(k,j,i) = 0.0_GP
               END IF
            END DO
         END DO
      END DO
 
      RETURN
      END SUBROUTINE dealias
