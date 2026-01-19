! Step 1 of Runge-Kutta for the CHMHDB equations
! Copies v_x,y,z into the auxiliary arrays C1-3
! and a_x,y,z into the auxiliary arrays C4-6

         C1(k,j,i) = vx(k,j,i)
         C2(k,j,i) = vy(k,j,i)
         C3(k,j,i) = vz(k,j,i)
         C4(k,j,i) = ax(k,j,i)
         C5(k,j,i) = ay(k,j,i)
         C6(k,j,i) = az(k,j,i)
         C20(k,j,i) = th(k,j,i)
