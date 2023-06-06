            CALL picpart%GetTemperature(R1)

            CALL derivk3(phi,C1,1)
            CALL derivk3(phi,C2,2)
            CALL derivk3(phi,C3,3)

            CALL ehpiccheck(R1,C1,C2,C3,t,dt)
