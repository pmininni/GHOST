! Step 1 of Runge-Kutta for the MOIST equations
! Copies v_x,y,z into the auxiliary arrays C1-C3
! Copies th_i into the auxiliary arrays C20-C23

         C1 (k,j,i) = vx (k,j,i)
         C2 (k,j,i) = vy (k,j,i)
         C3 (k,j,i) = vz (k,j,i)
         C20(k,j,i) = th(k,j,i)   ! active scalar 1, zeta_u
         C21(k,j,i) = th1(k,j,i)  ! active scalar 2, zeta_s
!!       C22(k,j,i) = th2(k,j,i)  ! passive scalar 1
!!       C23(k,j,i) = th3(k,j,i)  ! passive scalar 2
