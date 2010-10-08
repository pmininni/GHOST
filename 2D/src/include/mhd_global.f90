! Global quantities computed in MHD runs

            CALL mhdcheck(ps,az,fk,mk,(t-1)*dt)
            CALL maxabs(ps,rmp,ki,kj)
            CALL maxabs(az,rmq,ki,kj)
            IF (myrank.eq.0) THEN
               OPEN(1,file='maximum.txt',position='append')
               WRITE(1,'(E13.6,E13.6,E13.6)') (t-1)*dt,rmp,rmq
               CLOSE(1)
            ENDIF
