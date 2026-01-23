! Step 2 of Runge-Kutta for the compressible HD equations
! Computes the nonlinear terms and evolves the equations in dt/o.
! Remember: we evolve the momentum *density* (rho v) here:

       if ( useRK3 .eq. 0 ) THEN

         CALL mom2vel(rho,sx,sy,sz,0,vx,vy,vz)    ! compute velocity
         CALL divrhov(rho,vx,vy,vz,0,C7)          ! div(rho.v)

         CALL divrhov(th,vx,vy,vz,0,C8)           ! div(e.v)
         CALL pdVwork(gam1,th,vx,vy,vz,0,C34)     ! p.div(v) = (gamma-1) e.div(v)
     !   CALL viscHeatRayleigh(vx,vy,vz,C36)      ! phi/mu, visc. heat

         CALL divrhov(sx,vx,vy,vz,0,C4)           ! div(rho vx.v)
         CALL divrhov(sy,vx,vy,vz,0,C5)           ! div(rho vy.v)
         CALL divrhov(sz,vx,vy,vz,0,C6)           ! div(rho vz.v)
     !   CALL gradre3(vx,vy,vz,C4,C5,C6)          ! v.Grad v
         CALL gradpressi(gam1,th,C31,C32,C33)     ! Grad p term
     !   CALL divide(rho,C31,C32,C33)             ! divide Grad p by rho
         ! Note: th, si overwritten, so make these calls last:
         CALL vdiss(nu,nu2,sx,sy,sz)              ! viscous term
     !   CALL divide(rho,sx,sy,sz)                ! divide viscous term by rho
         CALL laplak3(th,th)                      ! laplacian(e)
        
         rmp = 1./real(o,kind=GP)

         if ( .NOT. use_voigt ) THEN

!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
         DO k = 1,nz

            IF (kn2(k,j,i).le.kmax) THEN
               sx (k,j,i) = C1(k,j,i)+dt*(sx(k,j,i)-C4(k,j,i)-C31(k,j,i) &
              +fx(k,j,i))*rmp
               sy (k,j,i) = C2(k,j,i)+dt*(sy(k,j,i)-C5(k,j,i)-C32(k,j,i) &
              +fy(k,j,i))*rmp
               sz (k,j,i) = C3(k,j,i)+dt*(sz(k,j,i)-C6(k,j,i)-C33(k,j,i) &
              +fz(k,j,i))*rmp
               rho(k,j,i) = C20(k,j,i)-dt*C7(k,j,i)*rmp
!              th (k,j,i) = C35(k,j,i)+dt*(nu*C36(k,j,i)-C8(k,j,i)-C34(k,j,i) &
               th (k,j,i) = C35(k,j,i)+dt*(kappa*th(k,j,i)-C8(k,j,i)-C34(k,j,i) &
              +fs(k,j,i))*rmp
!           ELSE IF (kn2(k,j,i).gt.kmax) THEN
            ELSE 
               sx (k,j,i) = 0.0_GP
               sy (k,j,i) = 0.0_GP
               sz (k,j,i) = 0.0_GP
               rho(k,j,i) = 0.0_GP
               th (k,j,i) = 0.0_GP
!           ELSE IF (kn2(k,j,i).lt.tinyd) THEN
!              sx (k,j,i) = 0.0_GP
!              sy (k,j,i) = 0.0_GP
!              sz (k,j,i) = 0.0_GP
!              rho(k,j,i) = C20(k,j,i)
!              th (k,j,i) = C35(k,j,i)
            ENDIF

         END DO
         END DO
         END DO
     
         ELSE ! using Voigt
!
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
         DO k = 1,nz

            IF (kn2(k,j,i).le.kmax) THEN
               sx (k,j,i) = C1(k,j,i)+dt*(sx(k,j,i)-C4(k,j,i)-C31(k,j,i) &
              +fx(k,j,i))*rmp * Hvinv(k,j,i)
               sy (k,j,i) = C2(k,j,i)+dt*(sy(k,j,i)-C5(k,j,i)-C32(k,j,i) &
              +fy(k,j,i))*rmp * Hvinv(k,j,i)
               sz (k,j,i) = C3(k,j,i)+dt*(sz(k,j,i)-C6(k,j,i)-C33(k,j,i) &
              +fz(k,j,i))*rmp * Hvinv(k,j,i)
               rho(k,j,i) = C20(k,j,i)-dt*C7(k,j,i)*rmp
               th (k,j,i) = C35(k,j,i)+dt*(kappa*th(k,j,i)-C8(k,j,i)-C34(k,j,i) &
              +fs(k,j,i))*rmp * Hsinv(k,j,i)
            ELSE 
               sx (k,j,i) = 0.0_GP
               sy (k,j,i) = 0.0_GP
               sz (k,j,i) = 0.0_GP
               rho(k,j,i) = 0.0_GP
               th (k,j,i) = 0.0_GP
            ENDIF

         END DO
         END DO
         END DO

         ENDIF ! Voigt check
!
         CALL mom2vel(rho,sx,sy,sz,0,vx,vy,vz)    ! compute velocity update

      else ! useRK3 != 0: 

! Using COMPIHD solver: do RK3 step:
         IF ( o .NE. 1 ) THEN
           STOP 'ord must equal 1'
         ENDIF
 

         CALL COMPIRHS(sx,sy,sz,vx,vy,vz,rho,th, &
                       fx,fy,fz,fs,ttime,nu,nu2,kappa,gam1, &
                       c4,c5,c6,c7,c8,c20,c31,c32,c33,c34,c35,&
                       KVX1,KVY1,KVZ1,KD1,KE1)
        
         sx  = C1  + ( KVX1 * dt * 0.5 ) 
         sy  = C2  + ( KVY1 * dt * 0.5 ) 
         sz  = C3  + ( KVZ1 * dt * 0.5 ) 
         rho = C20 + ( KD1  * dt * 0.5 ) 
         th  = C35 + ( KE1  * dt * 0.5 ) 
         CALL COMPIRHS(sx,sy,sz,vx,vy,vz,rho,th, &
                       fx,fy,fz,fs,ttime+0.5*dt,nu,nu2,kappa,gam1, &
                       c4,c5,c6,c7,c8,c20,c31,c32,c33,c34,c35,&
                       KVX2,KVY2,KVZ2,KD2,KE2)

         sx  = C1  + ( KVX2 * 1.0 * dt )
         sy  = C2  + ( KVY2 * 1.0 * dt )
         sz  = C3  + ( KVZ2 * 1.0 * dt )
         rho = C20 + ( KD2  * 1.0 * dt )
         th  = C35 + ( KE2  * 1.0 * dt )
         CALL COMPIRHS(sx,sy,sz,vx,vy,vz,rho,th, &
                       fx,fy,fz,fs,ttime+dt,nu,nu2,kappa,gam1, &
                       c4,c5,c6,c7,c8,c20,c31,c32,c33,c34,c35,&
                       KVX3,KVY3,KVZ3,KD3,KE3)

         sx  = C1  + ( KVX1 + 4.0* KVX2 + KVX3 ) * ( dt / 6.0 ) 
         sy  = C2  + ( KVY1 + 4.0* KVY2 + KVX3 ) * ( dt / 6.0 ) 
         sz  = C3  + ( KVZ1 + 4.0* KVZ2 + KVX3 ) * ( dt / 6.0 ) 
         rho = C20 + ( KD1  + 4.0* KD2  + KD3  ) * ( dt / 6.0 ) 
         th  = C35 + ( KE1  + 4.0* KE2  + KE3  ) * ( dt / 6.0 ) 

!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
         DO k = 1,nz

            IF (kn2(k,j,i).gt.kmax) THEN
               sx (k,j,i) = 0.0_GP
               sy (k,j,i) = 0.0_GP
               sz (k,j,i) = 0.0_GP
               rho(k,j,i) = 0.0_GP
               th (k,j,i) = 0.0_GP
            ENDIF

         END DO
         END DO
         END DO

         CALL mom2vel(rho,sx,sy,sz,0,vx,vy,vz)    ! compute velocity update


      endif ! end, useRK3 check
