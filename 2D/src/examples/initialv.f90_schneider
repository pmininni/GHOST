! Initial condition for the velocity.
! This file contains the expression used for the initial
! streamfunction or velocity field (depending on the solver). 
! You can use the temporary real array R1 of size 
! (n,jsta:jend), and temporary complex arrays C1, C2 of size 
! (n,ista:iend) to do intermediate computations. The variable 
! u0 should control the global amplitude of the field, and 
! variables vparam0-9 can be used to control the amplitudes 
! of individual terms. At the end, the streamfunction should 
! be stored in the array ps, or the velocity field components 
! in the arrays vx and vy (plus vz in 2.5D solvers).

! Merging of two positive vortices and a negative vortex 
! (Schneider et al., Teor. Comp. Fluid Dyn. 9, 191).
! Vorticity-stream function formulation is used (2D).
!     u0     : field amplitude
!     vparam0: size of the vortices
      DO j = jsta,jend
         DO i = 1,n
            R1(i,j) = u0*                                                &
          (exp(-((2*pi*(real(i,kind=GP)-1)/real(n,kind=GP)-3*pi/4)**2    &
          +(2*pi*(real(j,kind=GP)-1)/real(n,kind=GP)-pi)**2)/vparam0**2) &
          +exp(-((2*pi*(real(i,kind=GP)-1)/real(n,kind=GP)-5*pi/4)**2    &
          +(2*pi*(real(j,kind=GP)-1)/real(n,kind=GP)-pi)**2)/vparam0**2) &
          -.5*exp(-((2*pi*(real(i,kind=GP)-1)/real(n,kind=GP)-5*pi/4)**2 &
          +(2*pi*(real(j,kind=GP)-1)/real(n,kind=GP)                     &
          -pi*(1+1./(2*sqrt(2.))))**2)/vparam0**2))
         END DO
      END DO
      CALL fftp2d_real_to_complex(planrc,R1,ps,MPI_COMM_WORLD)
      DO i = ista,iend
         DO j = 1,n
            IF ((ka2(j,i).le.kmax).and.(ka2(j,i).ge.tiny)) THEN
               ps(j,i) = ps(j,i)/ka2(j,i)
            ELSE
               ps(j,i) = 0.0_GP
            ENDIF
         END DO
      END DO
