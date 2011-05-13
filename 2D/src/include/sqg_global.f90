! Global quantities computed in SQG runs

            CALL sqgcheck(ps,fk,(t-1)*dt)
            CALL maxabs(ps,rmp,ki,kj)
            IF (myrank.eq.0) THEN
               OPEN(1,file='maximum.txt',position='append')
               WRITE(1,'(E13.6,E13.6)') (t-1)*dt,rmp
               CLOSE(1)
            ENDIF
