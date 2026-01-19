! Global quantities computed in ARGL runs

            CALL gpecheck(zre,zim,t,dt)
            CALL momentum(zre,zim,t,dt)
            CALL gpehelicity(zre,zim,t,dt)
