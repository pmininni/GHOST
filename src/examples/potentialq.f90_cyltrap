! Trapping potential and linear ramp functions for quantum solvers.
! This file contains the expression used for the trapping potential
! in quantum fluid solvers with external potentials, as well as
! linear ramp functions used to compute angular momentum operators
! in the presence of rotation. Note the z angular momentum operator
! is Lz = -i hbar (x d/dy - y d/dx). You can use temporary real arrays
! R1-R3 of size (1:nx,1:ny,ksta:kend) and temporary complex arrays
! C1-C8 of size (nz,ny,ista:iend) to do intermediate computations. At
! the end of the computation, the potential divided by beta in real
! space should be stored in Vtrap, and the linear ramps in x and y in
! Vlinx and Vliny, all of size (1:nx,1:ny,ksta:kend), and not normalized
! (i.e., of order nx*ny*nz). The parameter V0 = m.w0^2/(2.hbar), where
! w0 is the trapping frequency, can be used to control the amplitude of
! the potential. This results in a characteristic lengthscale for the
! trap of a0 = sqrt(hbar/(m.w0)) = (cspeed.lambda/(V0.sqrt(2)))^(1/4).
! The linear ramps should also be multiplied by the rotation rate if
! solvers in a rotating frame (RGPE and RARGL) are used.

! Harmonic cylindrical (cigar) trap in x and y, V = V0 (x^2 + y^2), and
! linear functions of the x and y coordinates multiplied by omegaz:

! We compute the linear functions and the parabola
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
         DO j = 1,ny
            rmt = pi*(2*real(j-1,kind=GP)/real(ny,kind=GP)-1.0_GP)
            DO i = 1,nx
               rms = pi*(2*real(i-1,kind=GP)/real(nx,kind=GP)-1.0_GP)
               Vlinx(i,j,k) = rms              ! omegaz.x
               Vliny(i,j,k) = rmt              ! omegaz.y
               Vtrap(i,j,k) = rms**2 + rmt**2  ! V0.(x^2 + y^2)
            END DO
         END DO
      END DO

! We filter these functions to satisfy the periodic boundary conditions
      CALL fftp3d_real_to_complex(planrc,Vlinx,C1,MPI_COMM_WORLD)
      CALL fftp3d_real_to_complex(planrc,Vliny,C2,MPI_COMM_WORLD)
      CALL fftp3d_real_to_complex(planrc,Vtrap,C3,MPI_COMM_WORLD)
      rmp  = nx*Dkx/17.0 ! sigma_x: width of the filter in x
      rmq  = ny*Dky/17.0 ! sigma_y: width of the filter in y
      dump = 1.0_GP/(real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))
!$omp parallel do if (iend-ista.ge.nth) private (j,k,rm1,rm2)
      DO i = ista,iend
         rms = exp(-(kx(i)/rmp)**2/2)
!$omp parallel do if (iend-ista.lt.nth) private (k,rm2)
         DO j = 1,ny
            rmt = exp(-(ky(j)/rmq)**2/2)
            DO k = 1,nz
               IF (kn2(k,j,i).le.kmax) THEN
                  C1(k,j,i) = C1(k,j,i) * rms       * dump
                  C2(k,j,i) = C2(k,j,i)       * rmt * dump
                  C3(k,j,i) = C3(k,j,i) * rms * rmt * dump
               ELSE
                  C1(k,j,i) = 0.0_GP
                  C2(k,j,i) = 0.0_GP
                  C3(k,j,i) = 0.0_GP
               END IF
            END DO
         END DO
      END DO
      CALL fftp3d_complex_to_real(plancr,C1,Vlinx,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,C2,Vliny,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,C3,Vtrap,MPI_COMM_WORLD)

! After filtering, the amplitudes may be off. We correct them to have the
! right values of the derivatives in the center of the domain for the
! trapping potential and for the linear ramps.

      rmp = Vlinx(nx/2+nx/4,ny/2,ksta) - Vlinx(nx/2-nx/4,ny/2,ksta)
      rmq = Vliny(nx/2,ny/2+ny/4,ksta) - Vliny(nx/2,ny/2-ny/4,ksta)
      rmp = real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP)* &
            omegaz*pi*Lx/rmp ! unnormalized omegaz/slope_x
      rmq = real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP)* &
            omegaz*pi*Ly/rmq ! unnormalized omegaz/slope_y
      rms = Vtrap(nx/2-nx/4-1,ny/2,ksta) - Vtrap(nx/2-nx/4,ny/2,ksta) &
          - Vtrap(nx/2+nx/4-1,ny/2,ksta) + Vtrap(nx/2+nx/4,ny/2,ksta)
      rmt = Vtrap(nx/2,ny/2-ny/4-1,ksta) - Vtrap(nx/2,ny/2-ny/4,ksta) &
          - Vtrap(nx/2,ny/2+ny/4-1,ksta) + Vtrap(nx/2,ny/2+ny/4,ksta)
      rms = real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP)* &
            8*V0*pi**2/((rms*real(nx,kind=GP)/Lx**2 +           &
            rmt*real(ny,kind=GP)/Ly**2)*beta) ! 2V0/(second_deriv*beta)
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
         DO j = 1,ny
            DO i = 1,nx
               Vlinx(i,j,k) = Vlinx(i,j,k)*rmp
               Vliny(i,j,k) = Vliny(i,j,k)*rmq
               Vtrap(i,j,k) = Vtrap(i,j,k)*rms
            END DO
         END DO
      END DO

! We write the normalized potentials and linear ramps to a file
      rmp = 1.0_GP/(real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))
      rmq = beta  /(real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
         DO j = 1,ny
            DO i = 1,nx
               R1(i,j,k) = Vlinx(i,j,k)*rmp
               R2(i,j,k) = Vliny(i,j,k)*rmp
               R3(i,j,k) = Vtrap(i,j,k)*rmq
            END DO
         END DO
      END DO
      CALL io_write(1,odir,'Vlinx','init',planio,R1)
      CALL io_write(1,odir,'Vliny','init',planio,R2)
      CALL io_write(1,odir,'Vtrap','init',planio,R3)
