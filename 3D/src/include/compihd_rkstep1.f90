! Step 1 of Runge-Kutta for the compressible HD 
! equations including mass & int. energy
! Copies v_x,y,z into the auxiliary arrays C1-3,
! th into the auxiliary array C20, rho into C21

         C1 (k,j,i) = sx(k,j,i)
         C2 (k,j,i) = sy(k,j,i)
         C3 (k,j,i) = sz(k,j,i)
         C20(k,j,i) = rho(k,j,i) ! mass
         C35(k,j,i) = th(k,j,i)  ! energy
