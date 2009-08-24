! Adjust the amplitude of the mechanical forcing to keep
! the energy constant in an EDQNM-LES hydrodynamic run

               CALL energy(vx,vy,vz,tmp,0)
               CALL cross(vx,vy,vz,fx,fy,fz,eps,1)
               CALL spectrumc(vx,vy,vz,1,heli,Eold,Hold)
               rmp = 0.
               DO i = 1,n/2+1
                   rmp = rmp+tve(i)*Eold(i)*ka(i)**2-tvh(i)*Hold(i)
               END DO
               CALL MPI_BCAST(rmp,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
               CALL MPI_BCAST(tmp,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
               CALL MPI_BCAST(eps,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
               Faux1(1:9) = Faux1(2:10)
               Faux2(1:9) = Faux2(2:10)
               ampl = ampl*(nu*tmp+rmp)/eps
               Faux1(10) = ampl
               Faux2(10) = nu*tmp+rmp-eps
               ampl = .9*ampl
               DO i = 1,9
                  ampl = ampl+(.1*Faux1(i)+.01*Faux2(i))/9.
               END DO
               fx = fx*ampl
               fy = fy*ampl
               fz = fz*ampl
