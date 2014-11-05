! Step 1 of Runge-Kutta for the PHD equations
! Copies v_x,y,z into the auxiliary arrays C1-C3
! Copies th into the auxiliary array C20

         C1 (k,j,i) = vx (k,j,i)
         C2 (k,j,i) = vy (k,j,i)
         C3 (k,j,i) = vz (k,j,i)
         C21(k,j,i) = th1(k,j,i)
         C22(k,j,i) = th2(k,j,i)
!!       C23(k,j,i) = th3(k,j,i)
