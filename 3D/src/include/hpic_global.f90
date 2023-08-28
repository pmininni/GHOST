            CALL picpart%GetDensity(R1)
            CALL picpart%GetTemperature(R2)
            CALL hpiccheck(R2,R1,ax,ay,az,gammae,betae,t,dt)
