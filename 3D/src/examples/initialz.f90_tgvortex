! Initial condition for the wavefunction.
! This file contains the expression used for the initial
! wavefunction. You can use temporary real arrays R1-R3
! of size (1:n,1:n,ksta:kend) and temporary complex arrays
! C1-C8 of size (n,n,ista:iend) to do intermediate
! computations. The variables rho0 and zparam0-9 can be used 
! to control properties of the initial wavefunction. At the
! end, the real and imaginary parts of the wavefunction in 
! spectral space should be stored in the arrays zre and zim.

! Array of vortices for the TG flow
!     zparam0 : amplitude of the perturbation
!     zparam1 : length of the perturbation (in 2.pi units)

! We generate the functions lambda (R1) and mu (R2)
!$omp parallel do if (kend-ksta.ge.nth) private (i,j,rmp,rmq)
      DO k = ksta,kend
         rmp = sqrt(2*abs(cos(2*pi*(real(k,kind=GP)-1)/real(n,kind=GP))))
         rmq = rmp*sign(1.,cos(2*pi*(real(k,kind=GP)-1)/real(n,kind=GP)))
!$omp parallel do if (kend-ksta.lt.nth) private (i)
         DO j = 1,n
            DO i = 1,n
               R1(i,j,k) = cos(2*pi*(real(i,kind=GP)-1)/real(n,kind=GP))*rmp
               R2(i,j,k) = cos(2*pi*(real(j,kind=GP)-1)/real(n,kind=GP))*rmq
            END DO
         END DO
      END DO

! We generate phi_e(R1-1./sqrt(2),R2)

      dump = 1.0_GP/sqrt(2.0_GP)
!$omp parallel do if (kend-ksta.ge.nth) private (i,j,rmp,rmq,cdump,cdumq)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i,rmp,rmq,cdump,cdumq)
         DO j = 1,n
            DO i = 1,n
               rmp = 1.0_GP/sqrt(((R1(i,j,k)-dump)**2+R2(i,j,k)**2) &
                               *((R2(i,j,k)-dump)**2+R1(i,j,k)**2)  &
                               *((R1(i,j,k)+dump)**2+R2(i,j,k)**2)  &
                               *((R2(i,j,k)+dump)**2+R1(i,j,k)**2))
               rmq = tanh(dump*sqrt((R1(i,j,k)-dump)**2+R2(i,j,k)**2)/lambda) &
                    *tanh(dump*sqrt((R2(i,j,k)-dump)**2+R1(i,j,k)**2)/lambda) &
                    *tanh(dump*sqrt((R1(i,j,k)+dump)**2+R2(i,j,k)**2)/lambda) &
                    *tanh(dump*sqrt((R2(i,j,k)+dump)**2+R1(i,j,k)**2)/lambda)
               cdump = ((R1(i,j,k)-dump)+im*R2(i,j,k)) &
                      *(R1(i,j,k)+im*(R2(i,j,k)-dump)) &
                      *((R1(i,j,k)+dump)+im*R2(i,j,k)) &
                      *(R1(i,j,k)+im*(R2(i,j,k)+dump))
               cdumq = sqrt(omegag/beta)*(cdump*rmp*rmq)**int(1.0_GP/(2*pi*alpha))
               R1(i,j,k) =  real(cdumq)
               R2(i,j,k) = aimag(cdumq)
            END DO
         END DO
      END DO

      CALL fftp3d_real_to_complex(planrc,R1,zre,MPI_COMM_WORLD)
      CALL fftp3d_real_to_complex(planrc,R2,zim,MPI_COMM_WORLD)
