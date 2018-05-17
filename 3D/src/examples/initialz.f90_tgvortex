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
! Use with initialv.f90_tg without normalization (i.e., the 
! call to 'normalize' in initialv should be commented out).
!     kdn    : minimum wave number
!     kup    : maximum wave number

      dump = 1.0_GP/sqrt(2.0_GP)

! Sets the normalization for the entire wavefunction
!$omp parallel do if (kend-ksta.ge.nth) private (i,j,rmp,rmq,cdump,cdumq)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i,rmp,rmq,cdump,cdumq)
         DO j = 1,ny
            DO i = 1,nx
               R1(i,j,k) = sqrt(omegag/beta)
               R2(i,j,k) = sqrt(omegag/beta)
            END DO
         END DO
      END DO

! Computes a superposition of vortex filaments at different wavenumbers
      DO ki = INT(kdn),INT(kup)

! We generate the functions lambda (rm1) and mu (rm2)
!$omp parallel do if (kend-ksta.ge.nth) private (i,j,rmp,rmq)
      DO k = ksta,kend
         rmp = sqrt(2*abs(cos(2*pi*ki*(real(k,kind=GP)-1)/real(nz,kind=GP))))
         rmq = rmp*sign(1.0_GP,cos(2*pi*ki*(real(k,kind=GP)-1)/real(nz,kind=GP)))
!$omp parallel do if (kend-ksta.lt.nth) private (i)
         DO j = 1,ny
            DO i = 1,nx

               rm1 = cos(2*pi*ki*(real(i,kind=GP)-1)/real(nx,kind=GP))*rmp
               rm2 = cos(2*pi*ki*(real(j,kind=GP)-1)/real(ny,kind=GP))*rmq

! We generate phi_e(rm1-1./sqrt(2),rm2)

               rms = 1.0_GP/sqrt(((rm1-dump)**2+rm2**2)  &
                                *((rm2-dump)**2+rm1**2)  &
                                *((rm1+dump)**2+rm2**2)  &
                                *((rm2+dump)**2+rm1**2))
               rmt = tanh(dump*sqrt((rm1-dump)**2+rm2**2)/lambda) &
                    *tanh(dump*sqrt((rm2-dump)**2+rm1**2)/lambda) &
                    *tanh(dump*sqrt((rm1+dump)**2+rm2**2)/lambda) &
                    *tanh(dump*sqrt((rm2+dump)**2+rm1**2)/lambda)
               cdump = ((rm1-dump)+im*rm2)  &
                       *(rm1+im*(rm2-dump)) &
                      *((rm1+dump)+im*rm2)  &
                       *(rm1+im*(rm2+dump))               
               cdumq = (cdump*rms*rmt)**int(1.0_GP/(2*pi*alpha*ki))

               R1(i,j,k) = R1(i,j,k)*real(cdumq)
               R2(i,j,k) = R2(i,j,k)*aimag(cdumq)

            END DO
         END DO
      END DO

      END DO

      CALL fftp3d_real_to_complex(planrc,R1,zre,MPI_COMM_WORLD)
      CALL fftp3d_real_to_complex(planrc,R2,zim,MPI_COMM_WORLD)
