! Step 1 of Runge-Kutta for the MHD equations
! Copies v_x,y,z into the auxiliary arrays C1-3
! and a_x,y,z into the auxiliary arrays C4-6

         C1(k,j,i) = ax(k,j,i)
         C2(k,j,i) = ay(k,j,i)
         C3(k,j,i) = az(k,j,i)
