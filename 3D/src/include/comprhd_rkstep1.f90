! Step 1 of Runge-Kutta for the compressible HD equations
! Copies v_x,y,z into the auxiliary arrays C1-3
! and th into the auxiliary array C20

         C1(k,j,i) = vx(k,j,i)
         C2(k,j,i) = vy(k,j,i)
         C3(k,j,i) = vz(k,j,i)
         C20(k,j,i) = th(k,j,i)
