! Global quantities computed in Lagrangian averaged MHD runs

            CALL amhdcheck(vx,vy,vz,ax,ay,az,t,dt,alpk,alpm,1,1)
            CALL mhdcheck(vx,vy,vz,ax,ay,az,t,dt,1,1,0)
            CALL across(vx,vy,vz,fx,fy,fz,eps,alpk,1)
            CALL across(ax,ay,az,mx,my,mz,epm,alpm,0)
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

