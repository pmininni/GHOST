! Global quantities computed in BOUSS runs

            CALL hdcheck(vx,vy,vz,fx,fy,fz,t,dt,1,0)
            CALL pscheck(th,fs,t,dt)
            CALL maxabs(vx,vy,vz,rmp,0)
            CALL derivk3(th,C1,1)
            CALL derivk3(th,C2,2)
            CALL derivk3(th,C3,3)
            CALL maxabs(C1,C2,C3,rmq,2)
            IF (myrank.eq.0) THEN
               OPEN(1,file='maximum.txt',position='append')
               WRITE(1,FMT='(E13.6,E13.6,E13.6)') (t-1)*dt,rmp,rmq
               CLOSE(1)
            ENDIF
            CALL tbouss(vx,vy,vz,t,dt)
