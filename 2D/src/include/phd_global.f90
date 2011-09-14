! Global quantities computed in HD runs
            CALL hdcheck(ps,fk,(t-1)*dt)
            CALL pscheck(th,fs,(t-1)*dt,1)
            CALL derivk2(th,C1,1)
            CALL derivk2(th,C2,2)
            CALL maxabs2C(C1,C2,rmp,ki,kj)
            CALL maxabs(th,rmq,ki,kj)
            IF (myrank.eq.0) THEN
               OPEN(1,file='maximum.txt',position='append')
               WRITE(1,'(E13.6,E13.6,E13.6)') (t-1)*dt,rmp,rmq
               CLOSE(1)
            ENDIF
            CALL energy25(ps,th,tmp,1) ! energy
            CALL energy25(ps,th,tmq,0) ! vorticity
            CALL helicity25(ps,th,tmpp) ! helicity
            IF (myrank.eq.0) THEN
               OPEN(1,file='global25.txt',position='append')
               WRITE(1,'(2(E13.6,1X,E13.6,1X))') (t-1)*dt,tmp,tmq,tmpp
               CLOSE(1)
            ENDIF
