! Adjust the amplitude of the mechanical forcing
! to keep the energy constant in a hydrodynamic run

               CALL aenergy(vx,vy,vz,tmp,alpk,0)
               CALL across(vx,vy,vz,fx,fy,fz,eps,alpk,1)
               CALL MPI_BCAST(tmp,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
               CALL MPI_BCAST(eps,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
               Faux1(1:9) = Faux1(2:10)
               Faux2(1:9) = Faux2(2:10)
               ampl = ampl*nu*tmp/eps
               Faux1(10) = ampl
               Faux2(10) = nu*tmp-eps
               ampl = .9*ampl
               DO i = 1,9
                  ampl = ampl+(.1*Faux1(i)+.01*Faux2(i))/9.
               END DO
               fx = fx*ampl
               fy = fy*ampl
               fz = fz*ampl
