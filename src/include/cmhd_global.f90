! Global quantities computed in compressible MHD runs

            CALL energycompr (gam1,cp1, th,vx,vy,vz,t,dt)
            CALL mhdcheck(vx,vy,vz,ax,ay,az,t,dt,1,1,1)
            CALL cross(vx,vy,vz,fx,fy,fz,eps,1)
            CALL cross(ax,ay,az,mx,my,mz,epm,0)
            CALL maxabs(vx,vy,vz,rmp,0)
            CALL maxabs(ax,ay,az,rmq,1)
            IF (myrank.eq.0) THEN
               OPEN(1,file='injection.txt',position='append')
               WRITE(1,FMT='(E13.6,E22.14,E22.14)') (t-1)*dt,eps,epm
               CLOSE(1)
               OPEN(1,file='maximum.txt',position='append')
               WRITE(1,FMT='(E13.6,E13.6,E13.6)') (t-1)*dt,rmp,rmq
               CLOSE(1)
            ENDIF
