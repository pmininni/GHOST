! Global quantities computed in BOUSS runs

            CALL hdcheck(vx,vy,vz,fx,fy,fz,t,dt,1,1)
            CALL pscheck(th,fs,t,dt)
            CALL tbouss(vx,vy,vz,th,t,dt,0.0_GP,bvfreq)
            CALL maxabs(vx,vy,vz,rmp,0)
            CALL mpscheck2(th1,fs1,th2,fs2,t,dt)
!!          CALL mpscheck3(th1,fs1,th2,fs2,th3,fs3,t,dt)
            CALL derivk3(th,C1,1)
            CALL derivk3(th,C2,2)
            CALL derivk3(th,C3,3)
            CALL maxabs(C1,C2,C3,rmq,2)
            CALL derivk3(th1,C1,1)
            CALL derivk3(th1,C2,2)
            CALL derivk3(th1,C3,3)
            CALL maxabs(C1,C2,C3,rmq1,2)
            CALL derivk3(th2,C1,1)
            CALL derivk3(th2,C2,2)
            CALL derivk3(th2,C3,3)
            CALL maxabs(C1,C2,C3,rmq2,2)
!!          CALL derivk3(th3,C1,1)
!!          CALL derivk3(th3,C2,2)
!!          CALL derivk3(th3,C3,3)
!!          CALL maxabs(C1,C2,C3,rmq3,2)
            IF (myrank.eq.0) THEN
               OPEN(1,file='maximum.txt',position='append')
               WRITE(1,FMT='(E13.6,E13.6,E13.6,E13.6,E13.6)')       &
                            (t-1)*dt,rmp,rmq,rmq1,rmq2
!!             WRITE(1,FMT='(E13.6,E13.6,E13.6,E13.6,E13.6,E13.6)') &
!!                          (t-1)*dt,rmp,rmq,rmq1,rmq2,rmq3
               CLOSE(1)
            ENDIF
