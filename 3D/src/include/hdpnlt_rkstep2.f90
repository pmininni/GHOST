! Step 2 of Runge-Kutta for the HDPNLT equations
! This solver integrates the Navier-Stokes equations with the penalty
! method, to model an obstacle submerged in the flow. A mean flow is
! imposed in the flow in the x direction, initially started from 0 and
! slowly increased in time. The shape of the obstacle should be set
! in initialv.f90 (see the file initialv.f90_penalty in the examples
! folder). No external forcing is used except for the mean flow and the
! penalization. When solving with particles, collisions between particles
! and the obstacle are computed if heavy inertial particles are used.
!
! Jan 2024: New solver to do penalty in GHOST by Kaustuv Lahiri 

! Computes the nonlinear terms and evolves the equations in dt/o
! We first compute the normalization factor for the FFTs
     tmp = 1.0_GP/ &
            (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))

! Ramp up the mean velocity
! Note fvparam8 and fparam9 are used to control the mean velocity
!    u0     : mean velocity in the steady state
!    vparam0: rate at which vparam8 is increased
!    vparam8: instantaneous mean velocity (used internally)

! Creates the penalty function everytime. If the body moves, this has to
! be computed every time after computing forces over the body. In the case
! the body is fixed (this example), this is required only after a fresh
! start or a restart. Otherwise we resuse the penalty function already
! stored in memory. The behavior is controlled by the status of 'ischi'.

     IF (ischi.eq.0) THEN
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
       DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
         DO j = 1,ny
           DO i = 1,nx	    
	   ! Creates the sphere
	   vparam3 = sqrt((2*pi*Lx*(i-nx*x0)/nx)**2 + (2*pi*Ly*(j-ny*y0)/ny)**2 + (2*pi*Lz*(k-nz*z0)/nz)**2)
	   chi(i,j,k) = (0.5)*(1 + tanh(75.0*(1-(vparam3/radius))))
	   ! Creates a grid
	   !IF ((i == 1 .OR. i == 2) .AND. &
    	   !	(((j >= (ny/6-ny/50) .AND. j <= (ny/6+ny/50)) .OR. (j >= (ny/2-ny/50) .AND. j <= (ny/2+ny/50)) &
	   !   .OR. (j >= (5*ny/6-ny/50) .AND. j <= (5*ny/6+ny/50)) .OR. &
    	   !	((k >= (nz/6-nz/50) .AND. k <= (nz/6+nz/50))) .OR. (k >= (nz/2-nz/50) .AND. k <= (nz/2+nz/50)) &
	   !   .OR. (k >= (5*ny/6-ny/50) .AND. k <= (5*ny/6+ny/50))))) THEN
	   !   chi(i,j,k) = 1.0_GP
	   !ENDIF
           END DO
         END DO
       END DO
     ENDIF
     ischi = 1

! We add the mean velocity temporarily to compute advection
     IF (myrank.eq.0) THEN          ! ux = ux + U_0
        vx(1,1,1) = vx(1,1,1) + vparam8*real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP)
     ENDIF
     CALL gradre3(vx,vy,vz,C4,C5,C6)
     IF (myrank.eq.0) THEN          ! ux = ux - U_0
        vx(1,1,1) = vx(1,1,1) - vparam8*real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP)
     ENDIF

! Volume penalization 
     C34 = vx
	 CALL fftp3d_complex_to_real(plancr,C34,R1,MPI_COMM_WORLD)
	 DO k = ksta, kend
	 DO j = 1,ny
     	 DO i = 1,nx 
		R2(i,j,k) = inveta*R1(i,j,k)*tmp*chi(i,j,k)
	 END DO
	 END DO
	 END DO 
	 CALL fftp3d_real_to_complex(planrc, R2, C31, MPI_COMM_WORLD)
	 C4 = C4 + C31 

     C34 = vy
	 CALL fftp3d_complex_to_real(plancr,C34,R1,MPI_COMM_WORLD)
	 DO k = ksta, kend
	 DO j = 1,ny
         DO i = 1,nx 
		R2(i,j,k) = inveta*R1(i,j,k)*tmp*chi(i,j,k)
	 END DO
	 END DO
	 END DO 
	 CALL fftp3d_real_to_complex(planrc, R2, C32, MPI_COMM_WORLD)
	 C5 = C5 + C32

     C34 = vz
	 CALL fftp3d_complex_to_real(plancr,C34,R1,MPI_COMM_WORLD)
	 DO k = ksta,kend
	 DO j = 1,ny
         DO i = 1,nx 
		R2(i,j,k) = inveta*R1(i,j,k)*tmp*chi(i,j,k)
	 END DO
	 END DO
	 END DO 
	 CALL fftp3d_real_to_complex(planrc, R2, C33, MPI_COMM_WORLD)
	 C6 = C6 + C33

! Calculate the U0.del(U) terms for the time derivative
     CALL derivk3(vx, C31, 1)
     CALL derivk3(vy, C32, 1)
     CALL derivk3(vz, C33, 1)
     C4 = C4 - vparam8*C31
     C5 = C5 - vparam8*C32
     C6 = C6 - vparam8*C33

! Evaluate the individual terms and solve the NS equations in Fourier
! space
    	 CALL nonlhd3(C4,C5,C6,C7,1)
    	 CALL nonlhd3(C4,C5,C6,C8,2)
    	 CALL nonlhd3(C4,C5,C6,C4,3)	

		 
         CALL laplak3(vx,vx)
         CALL laplak3(vy,vy)
         CALL laplak3(vz,vz)

         IF ((trans.eq.1).and.(times.eq.0).and.(bench.eq.0).and.(o.eq.ord)) &
            CALL entrans(C1,C2,C3,C7,C8,C4,ext,1)
         rmp = 1./real(o,kind=GP)

!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
         DO k = 1,nz

            ! We don't dealias the 0 frequency mode
            IF ((kn2(k,j,i).le.kmax)) THEN
               vx(k,j,i) = C1(k,j,i) + dt*(nu*vx(k,j,i) + C7(k,j,i))*rmp
               vy(k,j,i) = C2(k,j,i) + dt*(nu*vy(k,j,i) + C8(k,j,i))*rmp
               vz(k,j,i) = C3(k,j,i) + dt*(nu*vz(k,j,i) + C4(k,j,i))*rmp
                        ! Semi-implicit time integration
	                !tmp = 1/(1 + 0.5*dt*nu*kk2(k,j,i))
			!vx(k,j,i) = C1(k,j,i)*(1 - 0.5*dt*nu*kk2(k,j,i))*tmp + dt*C7(k,j,i)*tmp*rmp
			!vy(k,j,i) = C2(k,j,i)*(1 - 0.5*dt*nu*kk2(k,j,i))*tmp + dt*C8(k,j,i)*tmp*rmp
			!vz(k,j,i) = C3(k,j,i)*(1 - 0.5*dt*nu*kk2(k,j,i))*tmp + dt*C4(k,j,i)*tmp*rmp
            ELSE
               vx(k,j,i) = 0.0_GP
               vy(k,j,i) = 0.0_GP
               vz(k,j,i) = 0.0_GP
            ENDIF

	 END DO	
 	 END DO	 
	 END DO
