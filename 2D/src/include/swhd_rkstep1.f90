! Step 1 of Runge-Kutta for the SWHD equations
! Copies vx and vy into the auxiliary arrays C1 and C2
! Copies th-fs into the auxiliary arrays C5 and C12

         C1(j,i)  = vx(j,i)
         C2(j,i)  = vy(j,i)
         C12(j,i) = th(j,i)
