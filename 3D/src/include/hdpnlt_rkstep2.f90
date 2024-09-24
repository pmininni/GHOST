! Step 2 of Runge-Kutta for the HDPNLT equations
! This solver integrates the Navier-Stokes equations with the penalty
! method, to model an obstacle submerged in the flow. A mean flow is
! imposed in the flow in the x direction, initially started from 0 and
! slowly increased in time. The shape of the obstacle should be set
! in initialfv.f90 (see the file initialfv.f90_penalty in the examples
! folder). No external forcing is used except for the mean flow and the
! penalization. When solving with particles, collisions between particles
! and the obstacle are computed if heavy inertial particles are used.

! Computes the nonlinear terms and evolves the equations in dt/o
     tmp = 1.0_GP/ &
            (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))
! Ramp up the mean velocity
! Note fvparam8 and fparam9 are used to control the mean velocity
!    u0     : mean velocity in the steady state
!    vparam0: rate at which vparam8 is increased
!    vparam8: instantaneous mean velocity (used internally)
	 IF (vparam8.ge.u0) THEN
	 	vparam8 = vparam8
	 ELSE
	 	vparam8 = vparam8 + vparam0
	 ENDIF

! Add the constant x velocity
     C34 = vx
         IF (myrank.eq.0) C34(1,1,1) = vparam8*real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP)
         CALL prodre3(C34,vy,vz,C4,C5,C6)

! Volume penalization 
     C34 = vx
	 CALL fftp3d_complex_to_real(plancr,C34,R1,MPI_COMM_WORLD)
	 DO k = ksta, kend
	 DO j = 1,ny
         DO i = 1,nx 
			R2(i,j,k) = inveta*R1(i,j,k)*chi(i,j,k)*tmp
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
			R2(i,j,k) = inveta*R1(i,j,k)*chi(i,j,k)*tmp
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
			R2(i,j,k) = inveta*R1(i,j,k)*chi(i,j,k)*tmp
	 END DO
	 END DO
	 END DO 
	 CALL fftp3d_real_to_complex(planrc, R2, C33, MPI_COMM_WORLD)
	 C6 = C6 + C33

! Evaluate the individual terms and solve the NS eqautions in Fourier
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
               vx(k,j,i) = C1(k,j,i)+dt*(nu*vx(k,j,i)+C7(k,j,i))*rmp
               vy(k,j,i) = C2(k,j,i)+dt*(nu*vy(k,j,i)+C8(k,j,i))*rmp
               vz(k,j,i) = C3(k,j,i)+dt*(nu*vz(k,j,i)+C4(k,j,i))*rmp
            ELSE
               vx(k,j,i) = 0.0_GP
               vy(k,j,i) = 0.0_GP
               vz(k,j,i) = 0.0_GP
            ENDIF

	 END DO	
 	 END DO	 
	 END DO
