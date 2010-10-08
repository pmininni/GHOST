! Global quantities computed in HMHD runs

            CALL mhdcheck25(ps,az,vz,bz,fk,mk,fz,mz,(t-1)*dt,1)
            CALL maxabs(ps,rmp,ki,kj)
            CALL maxabs(az,rmq,ki,kj)
            IF (myrank.eq.0) THEN
               OPEN(1,file='maximum.txt',position='append')
               WRITE(1,'(E13.6,E13.6,E13.6)') (t-1)*dt,rmp,rmq
               CLOSE(1)
            ENDIF
