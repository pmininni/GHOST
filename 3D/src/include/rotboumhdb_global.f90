! Global quantities computed in BOUSS runs

            CALL mhdcheck(vx,vy,vz,ax,ay,az,t,dt,1,1,0)
            CALL pscheck(th,fs,t,dt)
            CALL cross(vx,vy,vz,fx,fy,fz,eps,1)
            CALL cross(ax,ay,az,mx,my,mz,epm,0)
            IF (myrank.eq.0) THEN
               OPEN(1,file='injection.txt',position='append')
               WRITE(1,FMT='(E13.6,E22.14,E22.14)') (t-1)*dt,eps,epm
               CLOSE(1)
            ENDIF
