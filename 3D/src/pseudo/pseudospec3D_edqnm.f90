!=================================================================
! PSEUDOSPECTRAL subroutines
!
! Extra subroutines to compute nonlinear and other terms in the 
! EDQNM-based LES model in 3D using a pseudo-spectral method.
! You should use the FFTPLANS and MPIVARS modules (see the 
! file 'fftp_mod.f90') in each program that calls any of the 
! subroutines in this file. 
!
! This EDQNM model was developed by Julien Baerenzung
! (JAS 2008); most of the code was adapted from Julien's 
! C++ code. 
!
! NOTATION: index 'i' is 'x' 
!           index 'j' is 'y'
!           index 'k' is 'z'
!
! 2009 D. Rosenberg and Pablo D. Mininni.
!      National Center for Atmospheric Research.
!      e-mail: mininni@ucar.edu
!=================================================================

!*****************************************************************
      SUBROUTINE evedqnm(Eold,Hold,nu,time,hel,tepq,thpq,tve,tvh,Ext,Hxt)
!-----------------------------------------------------------------
!
! Compute the turbulent eddy viscosity and helicity eddy diffusivity 
! as well as the kinetic and helicity modal transfer function for EDQNM 
! closure in spectral space. 
!
! Parameters
!     Eold: resolved energy spectum (in)
!     Hold: resolved helicity spectrum (in)
!     nu  : kinematic (real) viscosity (in)
!     time: time (int)
!     hel : flag (int) s.t.:
!            1: compute helicity quantities
!            0: skip computation of helicity quantities
!     tepq: modal energy transfer function (out)
!     thpq: modal helicity transfer function (out)
!     tve : turbulent eddy viscosity (out)
!     tvh : turbulent eddy helicity diffusivity (out)
!     Ext : enery spectrum extrapolated to 3kc grid (out)
!     Hxt : enery spectrum extrapolated to 3kc grid (out)
!
      USE ali
      USE mpivars
      USE grid
      USE fft
      IMPLICIT NONE

      DOUBLE PRECISION, INTENT (IN), DIMENSION(n/2+1)     :: Eold,Hold
      DOUBLE PRECISION, INTENT(OUT), DIMENSION(n/2+1)     :: tepq,thpq
      DOUBLE PRECISION, INTENT(OUT), DIMENSION(n/2+1)     :: tve,tvh
      DOUBLE PRECISION, INTENT(OUT), DIMENSION(3*(n/2+1)) :: Ext,Hxt
      REAL   , INTENT(IN)   :: nu,time
      REAL                  :: ddelta,delk,power,hddelta,hpower
      INTEGER, INTENT(IN)   :: hel
      INTEGER               :: k,kc

      kc = int(sqrt(kmax))
      delk = 3.*kc/4.
      DO k = 1,n/2+1
         tve (k) = 0.
         tvh (k) = 0.
         tepq(k) = 0.
         thpq(k) = 0.
      ENDDO
      IF ( hel.eq.1 ) THEN
         call spec_fit(Eold,delk-2,Ext,power,ddelta)
         call spec_fit(Hold,delk,Hxt,hpower,hddelta)
         IF ( power.gt.0 ) THEN
            IF ( hpower.gt.0. ) THEN
              call full_cusp(Ext,Hxt,nu,time,tepq,thpq,tve,tvh)
            ELSE
              call cusp(Ext,nu,time,tepq,tve)
            ENDIF
         ENDIF
      ELSE
         call spec_fit(Eold,delk,Ext,power,ddelta)
         IF ( power.gt.0. ) THEN
            call cusp(Ext,nu,time,tepq,tve)
         ENDIF     
      ENDIF

      RETURN
      END SUBROUTINE evedqnm

!*****************************************************************
      SUBROUTINE cusp(Ext,nu,t0,tepq,nut)
!-----------------------------------------------------------------
!
! Compute the turbulent eddy viscosity and kinetic transfer 
! function for EDQNM closure in spectral space, without helicity. 
!
! Parameters
!     Ext : energy spectrum extrapolated to 3k_c (in)
!     nu  : kinematic viscosity (in)
!     t0  : time (in)
!     tepq: modal transfer function (out)
!     nut : eddy viscosity (out)
!
      USE ali
      USE mpivars
      USE grid
      USE fft
      USE edqnm
      IMPLICIT NONE

      DOUBLE PRECISION, INTENT (IN), DIMENSION(3*(n/2+1)) :: Ext
      DOUBLE PRECISION, INTENT(OUT), DIMENSION(n/2+1)     :: tepq,nut
      REAL, INTENT(IN)                                    :: nu,t0
      REAL, DIMENSION(n,n,ksta:kend) :: r1,r2
      REAL, DIMENSION(3*(n/2+1))     :: mu
      REAL    :: coeff,en,frac,ppi,qqi,kk,lambda,pp,qq
      REAL    :: theta,tmp,tvvkp,tvvkpsym,tvvpq,tvvpqsym,x,y,z
      INTEGER :: i,j,k,kc,ktot,p,q

      lambda = 0.218*(kolmo**1.5)
      kc     = int(sqrt(kmax))
      ktot   = 3*kc

      IF ( (Ext(kc+1).eq.0.) ) RETURN

      en = 0.
      DO k = 1,ktot
        tmp   = real(k*k)
        en    = en + tmp*Ext(k)
        mu(k) = lambda*sqrt(en)
      ENDDO

      DO k = 1,kc
         kk = real(k)
         DO p = kc+1, ktot
             pp  = real(p)
             ppi = 1.0/pp
             DO q = p-k, p
                qq  = real(q)
                qqi = 1.0/qq
                IF ( q.eq.0. ) CYCLE
                coeff = 0.5
                IF ( ((qq-0.5).ge.(pp-kk)) .and. ((qq+0.5).le.pp) ) THEN
                  coeff = 1.0
                ENDIF
                x = ( pp**2 + qq**2 - kk**2 ) / (2.0*pp*qq)
                y = ( kk**2 + qq**2 - pp**2 ) / (2.0*kk*qq)
                z = ( pp**2 + kk**2 - qq**2 ) / (2.0*pp*kk)
                tvvkp    = -qqi*(x*y + z**3) * pp**2 * Ext(q)*Ext(k)
                tvvpq    =  qqi*(x*y + z**3) * kk**2 * Ext(q)*Ext(p)
                tvvkpsym = -ppi*(x*z + y**3) * qq**2 * Ext(p)*Ext(k) 
                tvvpqsym =  ppi*(x*z + y**3) * kk**2 * Ext(q)*Ext(p) 
                tmp      = mu(k)+mu(p)+mu(q) +nu*(kk**2 + pp**2 + qq**2)
                frac     = 1.0/tmp
                theta    = -exp(-tmp*t0) / tmp 
                nut(k)   = nut(k) - coeff*(frac+theta) &
                                  *(tvvkp+tvvkpsym)/(2.*kk*kk*Ext(k))
                tepq(k)  = tepq(k)+ coeff*(frac+theta) &
                                  *(tvvkp+tvvkpsym)
             ENDDO
         ENDDO
      ENDDO

      RETURN
      END SUBROUTINE cusp

!*****************************************************************
      SUBROUTINE full_cusp(Ext,Hxt,nu,t0,tepq,thpq,nut,tvh)
!-----------------------------------------------------------------
!
! Compute the turbulent eddy viscosity and kinetic transfer 
! function for EDQNM closure in spectral space with helicity. 
!
! Parameters
!     Ext : energy spectrum extrapolated to 3k_c (in)
!     Hxt : helicity spectrum extrapolated to 3k_c (in)
!     nu  : kinematic viscosity (in)
!     t0  : time (in)
!     tepq: modal energy transfer function (out)
!     thpq: modal helicity transfer function (out)
!     nut : eddy viscosity (out)
!     tvh : helicity eddy diffusivity (out)
!
      USE ali
      USE mpivars
      USE grid
      USE fft
      USE edqnm
      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(OUT), DIMENSION(3*(n/2+1)) :: Ext,Hxt
      DOUBLE PRECISION, INTENT(OUT), DIMENSION(n/2+1)     :: tepq,thpq,nut,tvh
      REAL, INTENT(IN)                                    :: nu,t0
      REAL, DIMENSION(n,n,ksta:kend) :: r1,r2
      REAL, DIMENSION(3*(n/2+1))     :: mu
      REAL    :: coeff,en,ppi,qqi,kk,lambda,pp,qq
      REAL    :: tmp,theta,tvvkp,tvvkpsym,tvvpq,tvvpqsym,tvw,tvw2
      REAL    :: twwkp,twwkpsym,twwpq,twwpqsym,x,y,z,xyz,xzy,ym1,zm1
      INTEGER :: i,j,k,kc,ktot,p,q

      lambda = 0.218*(kolmo**1.5)
      kc     = int(sqrt(kmax))
      ktot   = 3*kc

      IF ( (Ext(kc+1).eq.0.) ) RETURN

      DO k = 1,ktot
         Ext(k) = 0.5*Ext(k)
         Hxt(k) = 0.5*Hxt(k)
      ENDDO

      en = 0.
      DO k = 1,ktot
        tmp   = real(k*k)
        en    = en + tmp*Ext(k)
        mu(k) = lambda*sqrt(en)
      ENDDO

      DO k = 1,kc
         kk = real(k)
         DO p = kc+1, ktot
             pp  = real(p)
             ppi = 1.0/pp
             DO q = p-k, p
                qq  = real(q)
                qqi = 1.0/qq
                IF ( q.eq.0 ) CYCLE
                coeff = 0.5
                IF ( ((qq-0.5).ge.(pp-kk)) .and. ((qq+0.5).le.pp) ) THEN
                  coeff = 1.0
                ENDIF
                x = ( pp**2 + qq**2 - kk**2 ) / (2.0*pp*qq)
                y = ( kk**2 + qq**2 - pp**2 ) / (2.0*kk*qq)
                z = ( pp**2 + kk**2 - qq**2 ) / (2.0*pp*kk)
                xyz      = x*y + z**3
                xzy      = x*z + y**3
                ym1      = 1.-y**2
                zm1      = 1.-z**2
                tvvkp    = -qqi*xyz* pp**2 * Ext(q)*Ext(k)
                tvvpq    =  qqi*xyz* kk**2 * Ext(q)*Ext(p)
                tvvkpsym = -ppi*xzy* qq**2 * Ext(p)*Ext(k) 
                tvvpqsym =  ppi*xzy* kk**2 * Ext(q)*Ext(p) 
              
                twwpq    = -z*ym1*kk*kk*Hxt(p)*Hxt(q)/(pp*pp*qq)
                twwkp    =  z*ym1*      Hxt(k)*Hxt(q)*qqi
                twwpqsym = -y*zm1*kk*kk*Hxt(p)*Hxt(q)/(pp*qq*qq)
                twwkpsym =  y*zm1*      Hxt(k)*Hxt(p)*ppi

                tvw      =  qqi*xyz*kk*kk*Ext(q)*Hxt(p)       &
                         +  ppi*xzy*kk*kk*Hxt(q)*Ext(p)       &
                         -  qqi*z*ym1*kk*kk*Ext(p)*Hxt(q)     &
                         -  ppi*y*zm1*kk*kk*Ext(q)*Hxt(p) 

                tmp      = mu(k)+mu(p)+mu(q) + nu*(kk**2 + pp**2 + qq**2)
                theta    = (1.0-exp(-tmp*t0))/tmp

                nut (k)  = nut (k)-coeff*theta*(tvvkp+tvvkpsym)/ &
                           (2.*kk*kk*Ext(k))
                tepq(k)  = tepq(k)+coeff*theta*(tvvpq+tvvpqsym+twwpq+twwpqsym)
                tvh (k)  = tvh (k)+coeff*theta*(twwkp+twwkpsym)/(2.*Hxt(k))
                thpq(k)  = thpq(k)+coeff*theta*tvw
             ENDDO
         ENDDO
      ENDDO

      RETURN
      END SUBROUTINE full_cusp

!*****************************************************************
      SUBROUTINE spec_fit(Ein,delk,Efit,power,ddelta)
!-----------------------------------------------------------------
!
! Compute the fit to the end of the input spectrum, using logarithmic 
! decrement, and extrapolate it to create extended spectrum
!
! Parameters
!     Ein   : input spectrum (of the resolved scales; in)
!     delk  : number of resolved wave numbers to use in fit (in)
!     Efit  : spectrum including resolved to fit scale in buffer (out)
!     power : spectral index of fit (out)
!     ddelta: logarithmic decrement of spectral fit (out)
!
      USE ali
      USE mpivars
      USE grid
      USE fft
      IMPLICIT NONE

      DOUBLE PRECISION, INTENT (IN), DIMENSION(n/2+1)     :: Ein
      DOUBLE PRECISION, INTENT(OUT), DIMENSION(3*(n/2+1)) :: Efit
      REAL, INTENT (IN) :: delk
      REAL, INTENT(OUT) :: power, ddelta
      REAL      :: x,y,z,v,w,a,b,m,o,p,nn
      REAL      :: etot,fksup,ii,kk,expo,pnum,dden,pden
      INTEGER   :: i,j,k,kc,kinf,ksup,ktot

      kc     = int(sqrt(kmax))
      ktot   = 3*kc
      ksup   = kc
      kinf   = kc-delk

      power  = -1.
      ddelta = -1.
      DO k = 1,3*(n/2+1)
         Efit(k) = 0.
      ENDDO

      IF ( Ein(ksup).le.tinyf ) RETURN

      x = 0.
      y = 0.    
      z = 0.    
      v = 0.    
      w = 0.    
      m = 0.    
      nn= 0.    
      o = 0.    
      p = 0.    
      DO i = kinf,ksup
         ii = real(i)
         x  = x + 1.
         y  = y + log(ii)
         z  = z + log(Ein(i)+tinyf)
         v  = v + log(ii)*log(Ein(i)+tinyf)
         w  = w + log(ii)*log(ii)
         m  = m + ii
         nn = nn+ ii*log(ii)
         o  = o + ii*log(Ein(i)+tinyf)
         p  = p + ii*ii
      ENDDO
      pden  = (w*p-nn*nn)*(x*nn-y*m)-(y*nn-w*m)*(y*p-m*nn)
      pnum  = ((z*nn-v*m)*(y*p-m*nn)-(v*p-o*nn)*(x*nn-y*m)) 
      dden  = x*nn-y*m
      IF ( pden.ne.0. .and. dden.ne.0. ) THEN
         power  = pnum/pden
         expo   = ((z*nn-v*m)-power*(w*m-y*nn))/dden
         ddelta = expo*x/m - power*y/m - z/m
      ENDIF
      IF ( power.ge.0. .and. ddelta.ge.0. ) THEN
         expo   = exp(expo)
         DO k = 1,3*(n/2+1)
            Efit(k) = 0.0
         ENDDO
         DO k = 1,ksup
            Efit(k) = Ein(k)
         ENDDO
         DO k = ksup+1,ktot
            kk = real(k)
            Efit(k) = expo*(kk**(-power))*exp(-ddelta*kk)
         ENDDO
      ELSE
         x = 0.
         y = 0.    
         z = 0.    
         v = 0.    
         w = 0.    
         DO i = kinf,ksup
            ii = real(i)
            x  = x + 1.
            y  = y + log(ii)
            z  = z + log(Ein(i)+tinyf)
            v  = v + log(ii)*log(Ein(i)+tinyf)
            w  = w + log(ii)*log(ii)
         ENDDO
         pden  = y*y - w*x
         pnum  = (v*x - y*z)
         dden  = x*y
         IF ( pden.ne.0. .and. dden.ne.0. ) THEN
            power = pnum/pden
            expo  = exp((y*z + power*y*y)/dden)
            DO k=1,3*(n/2+1)
               Efit(k) = 0.
            ENDDO
            DO k = 1,ksup
               Efit(k) = Ein(k)
            ENDDO
            DO k = ksup+1,ktot
               kk = real(k)
               Efit(k) = expo*(kk**(-power)) + tinyf
            ENDDO
         ELSE
            DO k = 1,ksup
               Efit(k) = Ein(k)
            ENDDO
            DO k = ksup+1,3*(n/2+1)
               Efit(k) = 0.
            ENDDO
         ENDIF
      ENDIF

      RETURN
      END SUBROUTINE spec_fit

!*****************************************************************
      SUBROUTINE vcorrect(vx,vy,vz,ef,hf,tepq,thpq,eold,hold,dt,hel)
!-----------------------------------------------------------------
!
! After correcting spectra for subgrid effects, reconstruct the velocity
! corrections.
!
! Parameters
!     vx  : x-velocity (in/out)
!     vy  : y-velocity (in/out)
!     vz  : z-velocity (in/out)
!     ef  : array with the previous energy field
!     hf  : array with the previous helicity field
!     tepq: modal energy transfer function (in)
!     thpq: modal helicity transfer function (in)
!     eold: old energy spectrum
!     hold: old helicity spectrum
!     dt  : delta_t/ord
!     hel : flag (in) s.t.:
!            1: compute helicity quantities
!            0: skip computation of helicity quantities
!
      USE ali
      USE mpivars
      USE grid
      USE kes
      IMPLICIT NONE

      COMPLEX, INTENT(INOUT), DIMENSION(n,n,ista:iend) :: vx,vy,vz
      DOUBLE PRECISION, INTENT (IN), DIMENSION(n/2+1)  :: tepq,thpq
      DOUBLE PRECISION, INTENT (IN), DIMENSION(n/2+1)  :: eold,hold
      REAL, INTENT(IN), DIMENSION(n,n,ista:iend)       :: ef,hf
      REAL, INTENT(IN)                                 :: dt
      INTEGER, INTENT(IN)                              :: hel

      IF ( ista.eq.1 ) THEN 
! kx==0, ky,kz!=0
         call vcorr_x0ynzn(vx,vy,vz,ef,hf,tepq,thpq,eold,hold,dt,hel,0)
! kx!=0, ky,kz==0
         call vcorr_xny0z0(vx,vy,vz,ef,hf,tepq,thpq,eold,hold,dt,hel,1)
! kx==0, ky==0, kz!=0
         call vcorr_x0y0zn(vx,vy,vz,ef,hf,tepq,thpq,eold,hold,dt,hel,0)
! kx!=0, ky==0, kz!=0
         call vcorr_xny0zn(vx,vy,vz,ef,hf,tepq,thpq,eold,hold,dt,hel,1)
! kx==0, ky!=0, kz==0
         call vcorr_x0ynz0(vx,vy,vz,ef,hf,tepq,thpq,eold,hold,dt,hel,0)
! kx!=0, ky!=0, kz==0
         call vcorr_xnynz0(vx,vy,vz,ef,hf,tepq,thpq,eold,hold,dt,hel,1)
! kx!=0, ky!=0, kz!=0
         call vcorr_xnynzn(vx,vy,vz,ef,hf,tepq,thpq,eold,hold,dt,hel,1)
      ELSE
! kx!=0, ky,kz==0
         call vcorr_xny0z0(vx,vy,vz,ef,hf,tepq,thpq,eold,hold,dt,hel,0)
! kx!=0, ky==0, kz!=0
         call vcorr_xny0zn(vx,vy,vz,ef,hf,tepq,thpq,eold,hold,dt,hel,0)
! kx!=0, ky!=0, kz==0
         call vcorr_xnynz0(vx,vy,vz,ef,hf,tepq,thpq,eold,hold,dt,hel,0)
! kx!=0, ky!=0, kz!=0
         call vcorr_xnynzn(vx,vy,vz,ef,hf,tepq,thpq,eold,hold,dt,hel,0)
      ENDIF

      RETURN
      END SUBROUTINE vcorrect

!*****************************************************************
      SUBROUTINE vcorr_x0y0zn(vx,vy,vz,ef,hf,tepq,thpq,eold,hold,dt,hel,ip)
!-----------------------------------------------------------------
!
! Correction for case where kz != 0,  kx,ky == 0
!
! Parameters
!     vx  : x-velocity (in/out)
!     vy  : y-velocity (in/out)
!     vz  : z-velocity (in/out)
!     ef  : array with the previous energy field
!     hf  : array with the previous helicity field
!     tepq: modal energy transfer function (in)
!     thpq: modal helicity transfer function (in)
!     eold: old energy spectrum
!     hold: old helicity spectrum
!     dt  : delta_t/ord
!     hel : flag (in) s.t.:
!            1: compute helicity quantities
!            0: skip computation of helicity quantities
!
      USE ali
      USE mpivars
      USE grid
      USE kes
      USE var
      IMPLICIT NONE

      COMPLEX, INTENT(INOUT), DIMENSION(n,n,ista:iend) :: vx,vy,vz
      DOUBLE PRECISION, INTENT (IN), DIMENSION(n/2+1)  :: tepq,thpq
      DOUBLE PRECISION, INTENT (IN), DIMENSION(n/2+1)  :: eold,hold
      DOUBLE PRECISION    :: E,H,mx,my,v2,vx2,vy2,vz2,w2
      DOUBLE PRECISION    :: sxy,cxy,dsxy
      DOUBLE PRECISION    :: sx,sy,cx,cy,tx,ty,txy
      DOUBLE PRECISION    :: tmp,tmp1,tmp2,inorm
      DOUBLE PRECISION    :: vxr,vxi,vyr,vyi,vzr,vzi
      DOUBLE PRECISION    :: kx2,ky2,kz2
      REAL, INTENT(IN), DIMENSION(n,n,ista:iend) :: ef,hf
      REAL, INTENT(IN)    :: dt
      REAL                :: norm
      INTEGER, INTENT(IN) :: hel,ip
      INTEGER             :: i,j,k,km

      i = 1
      j = 1
      kx2= ka(i)**2
      ky2= ka(j)**2
      norm = float(n)**3
      inorm = 1./dble(norm)
      DO k = 2,n
         IF (ka2(k,j,i).gt.kmax) CYCLE
         kz2= ka(k)**2
         km = int(sqrt(ka2(k,j,i))+.501)        
         vxr = dble (vx(k,j,i))*inorm
         vxi = aimag(vx(k,j,i))*inorm
         vyr = dble (vy(k,j,i))*inorm
         vyi = aimag(vy(k,j,i))*inorm
         vzr = 0.
         vzi = 0.
         v2 = vxr**2+vxi**2+vyr**2+vyi**2+vzr**2+vzi**2
         E  = v2 + 2.0*tepq(km)*dt*ef(k,j,i)/eold(km)
         w2 = ka(k)*( vxr*vyi - vxi*vyr ) &
            + ka(j)*( vzr*vxi - vzi*vxr ) &
            + ka(i)*( vyr*vzi - vyi*vzr )
         IF ( hel.eq.1 ) THEN
            H = w2 + thpq(km)*dt*hf(k,j,i)/hold(km)
         ELSE
            H = w2
         ENDIF

         cx = vxr/sqrt(vxr**2+vxi**2)
         sx = vxi/sqrt(vxr**2+vxi**2)
         cy = vyr/sqrt(vyr**2+vyi**2)
         sy = vyi/sqrt(vyr**2+vyi**2)

         IF ( sx.gt.0. ) THEN
           tx = acos(cx)
         ELSE
           tx = 2.0*pi - acos(cx)
         ENDIF
         IF ( sy.gt.0. ) THEN
           ty = acos(cy)
         ELSE
           ty = 2.0*pi - acos(cy)
         ENDIF
         txy = ty-tx

         tmp1 = 1.
         tmp2 = 1.
         tmp = (E*sqrt(1.-H**2/(0.25*ka2(k,j,i)*E**2)))/   &
               (v2*sqrt(1.-w2**2/(0.25*ka2(k,j,i)*v2**2)))
         vx2 = tmp*(vxr**2+vxi**2-0.5*v2) + 0.5*E
         vy2 = tmp*(vyr**2+vyi**2-0.5*v2) + 0.5*E
         tmp1 = sign(tmp1,vx2)
         tmp2 = sign(tmp2,vy2)
         mx  = sqrt(abs(vx2))
         my  = sqrt(abs(vy2))
         sxy = H
         dsxy= ka(k)*mx*my

         IF ( tmp1.ge.0. .and. tmp2.ge.0. .and. abs(sxy).le.abs(dsxy) & 
            .and. dsxy.ne.0. ) THEN
            sxy = sxy/dsxy
            cxy = sqrt(1.-sxy*sxy)
            IF ( cos(txy).lt.0. ) THEN
            cxy = -cxy
            ENDIF
            IF ( mx.ge.my ) THEN
            vxr = mx*cx
            vxi = mx*sx
            tmp = my/mx
            vyr = tmp*(cxy*vxr-sxy*vxi)
            vyi = tmp*(sxy*vxr+cxy*vxi)
             ELSE
            vyr = my*cy
            vyi = my*sy
            tmp = mx/my
            vxr = tmp*(cxy*vyr+sxy*vyi)
            vxi = tmp*(-sxy*vyr+cxy*vyi)
            ENDIF
         ENDIF
         vx(k,j,i) = cmplx(real(vxr),real(vxi))*norm
         vy(k,j,i) = cmplx(real(vyr),real(vyi))*norm
         vz(k,j,i) = cmplx(real(vzr),real(vzi))*norm
      END DO

      RETURN
      END SUBROUTINE vcorr_x0y0zn 

!*****************************************************************
      SUBROUTINE vcorr_x0ynz0(vx,vy,vz,ef,hf,tepq,thpq,eold,hold,dt,hel,ip)
!-----------------------------------------------------------------
!
! Correction for case where ky != 0,  kx,kz == 0
!
! Parameters
!     vx  : x-velocity (in/out)
!     vy  : y-velocity (in/out)
!     vz  : z-velocity (in/out)
!     ef  : array with the previous energy field
!     hf  : array with the previous helicity field
!     tepq: modal energy transfer function (in)
!     thpq: modal helicity transfer function (in)
!     eold: old energy spectrum
!     hold: old helicity spectrum
!     dt  : delta_t/ord
!     hel : flag (in) s.t.:
!            1: compute helicity quantities
!            0: skip computation of helicity quantities
!
      USE ali
      USE mpivars
      USE grid
      USE kes
      USE var
      IMPLICIT NONE

      COMPLEX, INTENT(INOUT), DIMENSION(n,n,ista:iend) :: vx,vy,vz
      DOUBLE PRECISION, INTENT (IN), DIMENSION(n/2+1)  :: tepq,thpq
      DOUBLE PRECISION, INTENT (IN), DIMENSION(n/2+1)  :: eold,hold
      DOUBLE PRECISION    :: E,H,mx,mz,v2,vx2,vy2,vz2,w2
      DOUBLE PRECISION    :: szx,czx,dszx
      DOUBLE PRECISION    :: sx,sz,cx,cz,tx,tz,tzx
      DOUBLE PRECISION    :: tmp,tmp1,tmp2,inorm
      DOUBLE PRECISION    :: vxr,vxi,vyr,vyi,vzr,vzi
      DOUBLE PRECISION    :: kx2,ky2,kz2
      REAL, INTENT(IN), DIMENSION(n,n,ista:iend) :: ef,hf
      REAL, INTENT(IN)    :: dt
      REAL                :: norm
      INTEGER, INTENT(IN) :: hel,ip
      INTEGER             :: i,j,k,km

      i = 1
      k = 1
      kx2= ka(i)**2
      kz2= ka(k)**2
      norm = float(n)**3
      inorm = 1./dble(norm)
      DO j = 2,n
         IF (ka2(k,j,i).gt.kmax) CYCLE
         ky2= ka(j)**2
         km = int(sqrt(ka2(k,j,i))+.501)        
         vxr = dble (vx(k,j,i))*inorm
         vxi = aimag(vx(k,j,i))*inorm
         vyr = 0.
         vyi = 0.
         vzr = dble (vz(k,j,i))*inorm
         vzi = aimag(vz(k,j,i))*inorm
         v2 = vxr**2+vxi**2+vyr**2+vyi**2+vzr**2+vzi**2
         E  = v2 + 2.0*tepq(km)*dt*ef(k,j,i)/eold(km)
         w2 = ka(k)*( vxr*vyi - vxi*vyr ) &
            + ka(j)*( vzr*vxi - vzi*vxr ) &
            + ka(i)*( vyr*vzi - vyi*vzr )  
         IF ( hel.eq.1 ) THEN
            H = w2 + thpq(km)*dt*hf(k,j,i)/hold(km)
         ELSE
            H = w2
         ENDIF

         cx = vxr/sqrt(vxr**2+vxi**2)
         sx = vxi/sqrt(vxr**2+vxi**2)
         cz = vzr/sqrt(vzr**2+vzi**2)
         sz = vzi/sqrt(vzr**2+vzi**2)

         IF ( sx.gt.0. ) THEN
           tx = acos(cx)
         ELSE
           tx = 2.0*pi - acos(cx)
         ENDIF
         IF ( sz.gt.0. ) THEN
           tz = acos(cz)
         ELSE
           tz = 2.0*pi - acos(cz)
         ENDIF
         tzx = tx-tz

         tmp1 = 1.
         tmp2 = 1.
         tmp = (E*sqrt(1.-H**2/(0.25*ka2(k,j,i)*E**2)))/   &
               (v2*sqrt(1.-w2**2/(0.25*ka2(k,j,i)*v2**2)))
         vx2 = tmp*(vxr**2+vxi**2-0.5*v2) + 0.5*E
         vz2 = tmp*(vzr**2+vzi**2-0.5*v2) + 0.5*E
         tmp1 = sign(tmp1,vx2)
         tmp2 = sign(tmp2,vz2)
         mx  = sqrt(abs(vx2))
         mz  = sqrt(abs(vz2))
         szx = H
         dszx= ka(j)*mx*mz

         IF ( tmp1.ge.0. .and. tmp2.ge.0. .and. abs(szx).le.abs(dszx) &
            .and. dszx.ne.0. ) THEN
            szx = szx/dszx
            czx = sqrt(1.-szx*szx)
            IF ( cos(tzx).lt.0. ) THEN
            czx = -czx
            ENDIF
            IF ( mx.ge.mz ) THEN
            vxr = mx*cx
            vxi = mx*sx
            tmp = mz/mx
            vzr = tmp*( szx*vxi+czx*vxr)
            vzi = tmp*( czx*vxi-szx*vxr)
             ELSE
            vzr = mz*cz
            vzi = mz*sz
            tmp = mx/mz
            vxr = tmp*(-szx*vzi+czx*vzr)
            vxi = tmp*( szx*vzr+czx*vzi)
            ENDIF
         ENDIF
         vx(k,j,i) = cmplx(real(vxr),real(vxi))*norm
         vy(k,j,i) = cmplx(real(vyr),real(vyi))*norm
         vz(k,j,i) = cmplx(real(vzr),real(vzi))*norm
      END DO

      RETURN
      END SUBROUTINE vcorr_x0ynz0

!*****************************************************************
      SUBROUTINE vcorr_xny0z0(vx,vy,vz,ef,hf,tepq,thpq,eold,hold,dt,hel,ip)
!-----------------------------------------------------------------
!
! Correction for case where kx != 0,  ky,kz == 0
!
! Parameters
!     vx  : x-velocity (in/out)
!     vy  : y-velocity (in/out)
!     vz  : z-velocity (in/out)
!     ef  : array with the previous energy field
!     hf  : array with the previous helicity field
!     tepq: modal energy transfer function (in)
!     thpq: modal helicity transfer function (in)
!     eold: old energy spectrum
!     hold: old helicity spectrum
!     dt  : delta_t/ord
!     hel : flag (in) s.t.:
!            1: compute helicity quantities
!            0: skip computation of helicity quantities
!
      USE ali
      USE mpivars
      USE grid
      USE kes
      USE var
      IMPLICIT NONE

      COMPLEX, INTENT(INOUT), DIMENSION(n,n,ista:iend) :: vx,vy,vz
      DOUBLE PRECISION, INTENT (IN), DIMENSION(n/2+1)  :: tepq,thpq
      DOUBLE PRECISION, INTENT (IN), DIMENSION(n/2+1)  :: eold,hold
      DOUBLE PRECISION    :: E,H,my,mz,v2,vx2,vy2,vz2,w2
      DOUBLE PRECISION    :: syz,cyz,dsyz
      DOUBLE PRECISION    :: sy,sz,cy,cz,ty,tz,tyz
      DOUBLE PRECISION    :: tmp,tmp1,tmp2,inorm
      DOUBLE PRECISION    :: vxr,vxi,vyr,vyi,vzr,vzi
      DOUBLE PRECISION    :: kx2,ky2,kz2
      REAL, INTENT(IN), DIMENSION(n,n,ista:iend) :: ef,hf
      REAL, INTENT(IN)    :: dt
      REAL                :: norm
      INTEGER, INTENT(IN) :: hel,ip
      INTEGER             :: i,j,k,km

      j = 1
      k = 1
      ky2= ka(j)**2
      kz2= ka(k)**2
      norm = float(n)**3
      inorm = 1./dble(norm)
      DO i = ista+ip,iend
         IF (ka2(k,j,i).gt.kmax) CYCLE
         kx2= ka(i)**2
         km = int(sqrt(ka2(k,j,i))+.501)        
         vxr = 0.
         vxi = 0.
         vyr = dble (vy(k,j,i))*inorm
         vyi = aimag(vy(k,j,i))*inorm
         vzr = dble (vz(k,j,i))*inorm
         vzi = aimag(vz(k,j,i))*inorm
         v2 = vxr**2+vxi**2+vyr**2+vyi**2+vzr**2+vzi**2
         E  = v2 + 2.0*tepq(km)*dt*ef(k,j,i)/eold(km)
         w2 = ka(k)*( vxr*vyi - vxi*vyr ) &
            + ka(j)*( vzr*vxi - vzi*vxr ) &
            + ka(i)*( vyr*vzi - vyi*vzr )  
         IF ( hel.eq.1 ) THEN
            H = w2 + thpq(km)*dt*hf(k,j,i)/hold(km)
         ELSE
            H = w2
         ENDIF

         cy = vyr/sqrt(vyr**2+vyi**2)
         sy = vyi/sqrt(vyr**2+vyi**2)
         cz = vzr/sqrt(vzr**2+vzi**2)
         sz = vzi/sqrt(vzr**2+vzi**2)

         IF ( sy.gt.0. ) THEN
           ty = acos(cy)
         ELSE
           ty = 2.0*pi - acos(cy)
         ENDIF
         IF ( sz.gt.0. ) THEN
           tz = acos(cz)
         ELSE
           tz = 2.0*pi - acos(cz)
         ENDIF
         tyz = tz-ty

         tmp1 = 1.
         tmp2 = 1.
         tmp = (E*sqrt(1.-H**2/(0.25*ka2(k,j,i)*E**2)))/   &
               (v2*sqrt(1.-w2**2/(0.25*ka2(k,j,i)*v2**2)))
         vy2 = tmp*(vyr**2+vyi**2-0.5*v2) + 0.5*E
         vz2 = tmp*(vzr**2+vzi**2-0.5*v2) + 0.5*E
         tmp1 = sign(tmp1,vy2)
         tmp2 = sign(tmp2,vz2)
         my  = sqrt(abs(vy2))
         mz  = sqrt(abs(vz2))

         syz = H
         dsyz= ka(i)*mz*my

         IF ( tmp1.ge.0. .and. tmp2.ge.0. .and. abs(syz).le.abs(dsyz) &
            .and. dsyz.ne.0. ) THEN
            syz = syz/dsyz
            cyz = sqrt(1.-syz*syz)
            IF ( cos(tyz).lt.0. ) THEN
            cyz = -cyz
            ENDIF
            IF ( my.ge.mz ) THEN
            vyr = my*cy
            vyi = my*sy
            tmp = mz/my
            vzr = tmp*(-syz*vyi+cyz*vyr)
            vzi = tmp*( cyz*vyi+syz*vyr)
             ELSE
            vzr = mz*cz
            vzi = mz*sz
            tmp = my/mz
            vyr = tmp*( cyz*vzr+syz*vzi)
            vyi = tmp*(-syz*vzr+cyz*vzi)
            ENDIF
         ENDIF
         vx(k,j,i) = cmplx(real(vxr),real(vxi))*norm
         vy(k,j,i) = cmplx(real(vyr),real(vyi))*norm
         vz(k,j,i) = cmplx(real(vzr),real(vzi))*norm
      END DO

      RETURN
      END SUBROUTINE vcorr_xny0z0

!*****************************************************************
      SUBROUTINE vcorr_x0ynzn(vx,vy,vz,ef,hf,tepq,thpq,eold,hold,dt,hel,ip)
!-----------------------------------------------------------------
!
! Correction for case where ky,kz != 0,  kx == 0
!
!     vx  : x-velocity (in/out)
!     vy  : y-velocity (in/out)
!     vz  : z-velocity (in/out)
!     ef  : array with the previous energy field
!     hf  : array with the previous helicity field
!     tepq: modal energy transfer function (in)
!     thpq: modal helicity transfer function (in)
!     eold: old energy spectrum
!     hold: old helicity spectrum
!     dt  : delta_t/ord
!     hel : flag (in) s.t.:
!            1: compute helicity quantities
!            0: skip computation of helicity quantities
!
      USE ali
      USE mpivars
      USE grid
      USE kes
      USE var
      IMPLICIT NONE

      COMPLEX, INTENT(INOUT), DIMENSION(n,n,ista:iend) :: vx,vy,vz
      DOUBLE PRECISION, INTENT (IN), DIMENSION(n/2+1)  :: tepq,thpq
      DOUBLE PRECISION, INTENT (IN), DIMENSION(n/2+1)  :: eold,hold
      DOUBLE PRECISION    :: E,H,mx,my,mz,v2,vx2,vy2,vz2,w2
      DOUBLE PRECISION    :: szx,sxy,syz,czx,cxy,cyz
      DOUBLE PRECISION    :: dszx,dsxy,dsyz,dczx,dcxy,dcyz
      DOUBLE PRECISION    :: sx,sy,sz,cx,cy,cz,tx,ty,tz,tzx,txy,tyz
      DOUBLE PRECISION    :: tmp,tmp1,tmp2,tmp3,inorm
      DOUBLE PRECISION    :: vxr,vxi,vyr,vyi,vzr,vzi
      DOUBLE PRECISION    :: kx2,ky2,kz2
      REAL, INTENT(IN), DIMENSION(n,n,ista:iend) :: ef,hf
      REAL, INTENT(IN)    :: dt
      REAL                :: norm
      INTEGER, INTENT(IN) :: hel,ip
      INTEGER             :: i,j,k,km

      i = 1
      kx2= ka(i)**2
      norm = float(n)**3
      inorm = 1./dble(norm)
      DO j = 2,n
      ky2= ka(j)**2
      DO k = 2,n
         IF (ka2(k,j,i).gt.kmax) CYCLE
         kz2= ka(k)**2
         km = int(sqrt(ka2(k,j,i))+.501)        
         vxr = dble (vx(k,j,i))*inorm
         vxi = aimag(vx(k,j,i))*inorm
         vyr = dble (vy(k,j,i))*inorm
         vyi = aimag(vy(k,j,i))*inorm
         vzr = dble (vz(k,j,i))*inorm
         vzi = aimag(vz(k,j,i))*inorm
         v2 = vxr**2+vxi**2+vyr**2+vyi**2+vzr**2+vzi**2
         w2 = ka(k)*( vxr*vyi - vxi*vyr ) &
            + ka(j)*( vzr*vxi - vzi*vxr ) &
            + ka(i)*( vyr*vzi - vyi*vzr )  
         E  = v2 + 2.0*tepq(km)*dt*ef(k,j,i)/eold(km)
         IF ( hel.eq.1 ) THEN
            H = w2 + thpq(km)*dt*hf(k,j,i)/hold(km)
         ELSE
            H = w2
         ENDIF

         cx = vxr/sqrt(vxr**2+vxi**2)
         sx = vxi/sqrt(vxr**2+vxi**2)
         cy = vyr/sqrt(vyr**2+vyi**2)
         sy = vyi/sqrt(vyr**2+vyi**2)
         cz = vzr/sqrt(vzr**2+vzi**2)
         sz = vzi/sqrt(vzr**2+vzi**2)

         IF ( sx.gt.0. ) THEN
           tx = acos(cx)
         ELSE
           tx = 2.0*pi - acos(cx)
         ENDIF
         IF ( sy.gt.0. ) THEN
           ty = acos(cy)
         ELSE
           ty = 2.0*pi - acos(cy)
         ENDIF
         IF ( sz.gt.0. ) THEN
           tz = acos(cz)
         ELSE
           tz = 2.0*pi - acos(cz)
         ENDIF
         txy = ty-tx
         tzx = tx-tz
         tyz = tz-ty

         tmp1 = 1.
         tmp2 = 1.
         vx2 = ((E*sqrt(1.-H**2/(0.25*ka2(k,j,i)*E**2)))/ &
             (v2*sqrt(1.-w2**2/(0.25*ka2(k,j,i)*v2**2)))) &
             * (vxr**2+vxi**2-0.5*v2) + 0.5*E
         tmp1= sign(tmp1,vx2)
         mx  = sqrt(abs(vx2))
         tmp3= E-vx2
         tmp2= sign(tmp2,tmp3)
         tmp = 1./(ky2+kz2)
         my  = sqrt(abs(tmp3)*kz2*tmp)
         mz  = sqrt(abs(tmp3)*ky2*tmp)

         sxy = ka(k)*H
         dsxy= (ky2+kz2)*mx*my
         szx = ka(j)*H
         dszx= (ky2+kz2)*mx*mz

         IF ( tmp1.ge.0. .and. tmp2.ge.0. .and. abs(sxy).le.abs(dsxy) &
         .and. dsxy.ne.0. .and. abs(szx).le.abs(dszx) .and. dszx.ne.0. ) THEN
            sxy = sxy/dsxy
            szx = szx/dszx
            cxy = sqrt(1.-sxy*sxy)
            czx = sqrt(1.-szx*szx)
            IF ( cos(txy).lt.0. ) THEN
            cxy = -cxy
            ENDIF
            IF ( cos(tzx).lt.0. ) THEN
            czx = -czx
            ENDIF
            syz = -szx*cxy - sxy*czx;
            cyz =  czx*cxy - szx*sxy;
            IF ( mx.ge.my .and. mx.ge.mz ) THEN
            vxr = mx*cx
            vxi = mx*sx
            tmp = my/mx
            vyr = tmp*(cxy*vxr-sxy*vxi)
            vyi = tmp*(sxy*vxr+cxy*vxi)
            tmp = mz/mx
            vzr = tmp*(szx*vxi+czx*vxr)
            vzi = tmp*(czx*vxi-szx*vxr)
            ENDIF
            IF ( my.ge.mx .and. my.ge.mz ) THEN
            vyr = my*cy
            vyi = my*sy
            tmp = mx/my
            vxr = tmp*( cxy*vyr+sxy*vyi)
            vxi = tmp*(-sxy*vyr+cxy*vyi)
            tmp = mz/my
            vzr = tmp*(-syz*vyi+cyz*vyr)
            vzi = tmp*( cyz*vyi+syz*vyr)
            ENDIF
            IF ( mz.ge.mx .and. mz.ge.my ) THEN
            vzr = mz*cz
            vzi = mz*sz
            tmp = mx/mz
            vxr = tmp*(-szx*vzi+czx*vzr)
            vxi = tmp*( szx*vzr+czx*vzi)
            tmp = my/mz
            vyr = tmp*( cyz*vzr+syz*vzi)
            vyi = tmp*(-syz*vzr+cyz*vzi)
            ENDIF
         ENDIF
         vx(k,j,i) = cmplx(real(vxr),real(vxi))*norm
         vy(k,j,i) = cmplx(real(vyr),real(vyi))*norm
         vz(k,j,i) = cmplx(real(vzr),real(vzi))*norm
      END DO
      END DO

      RETURN
      END SUBROUTINE vcorr_x0ynzn

!*****************************************************************
      SUBROUTINE vcorr_xny0zn(vx,vy,vz,ef,hf,tepq,thpq,eold,hold,dt,hel,ip)
!-----------------------------------------------------------------
!
! Correction for case where kx,kz != 0,  ky == 0
!
! Parameters
!     vx  : x-velocity (in/out)
!     vy  : y-velocity (in/out)
!     vz  : z-velocity (in/out)
!     ef  : array with the previous energy field
!     hf  : array with the previous helicity field
!     tepq: modal energy transfer function (in)
!     thpq: modal helicity transfer function (in)
!     eold: old energy spectrum
!     hold: old helicity spectrum
!     dt  : delta_t/ord
!     hel : flag (in) s.t.:
!            1: compute helicity quantities
!            0: skip computation of helicity quantities
!
      USE ali
      USE mpivars
      USE grid
      USE kes
      USE var
      IMPLICIT NONE

      COMPLEX, INTENT(INOUT), DIMENSION(n,n,ista:iend) :: vx,vy,vz
      DOUBLE PRECISION, INTENT (IN), DIMENSION(n/2+1)  :: tepq,thpq
      DOUBLE PRECISION, INTENT (IN), DIMENSION(n/2+1)  :: eold,hold
      DOUBLE PRECISION    :: E,H,mx,my,mz,v2,vx2,vy2,vz2,w2
      DOUBLE PRECISION    :: szx,sxy,syz,czx,cxy,cyz
      DOUBLE PRECISION    :: dszx,dsxy,dsyz,dczx,dcxy,dcyz
      DOUBLE PRECISION    :: sx,sy,sz,cx,cy,cz,tx,ty,tz,tzx,txy,tyz
      DOUBLE PRECISION    :: tmp,tmp1,tmp2,tmp3,inorm
      DOUBLE PRECISION    :: vxr,vxi,vyr,vyi,vzr,vzi
      DOUBLE PRECISION    :: kx2,ky2,kz2
      REAL, INTENT(IN), DIMENSION(n,n,ista:iend) :: ef,hf
      REAL, INTENT(IN)    :: dt
      REAL                :: norm
      INTEGER, INTENT(IN) :: hel,ip
      INTEGER             :: i,j,k,km

      j = 1
      ky2= ka(j)**2
      norm = float(n)**3
      inorm = 1./dble(norm)
      DO i = ista+ip,iend
      kx2= ka(i)**2
      DO k = 2,n     
         IF (ka2(k,j,i).gt.kmax) CYCLE
         kz2= ka(k)**2
         km = int(sqrt(ka2(k,j,i))+.501)
         vxr = real (vx(k,j,i))*inorm
         vxi = aimag(vx(k,j,i))*inorm
         vyr = real (vy(k,j,i))*inorm
         vyi = aimag(vy(k,j,i))*inorm
         vzr = real (vz(k,j,i))*inorm
         vzi = aimag(vz(k,j,i))*inorm
         v2 = vxr**2+vxi**2+vyr**2+vyi**2+vzr**2+vzi**2
         E  = v2 + 2.0*tepq(km)*dt*ef(k,j,i)/eold(km)
         w2 = ka(k)*( vxr*vyi - vxi*vyr ) &
            + ka(j)*( vzr*vxi - vzi*vxr ) &
            + ka(i)*( vyr*vzi - vyi*vzr )  
         IF ( hel.eq.1 ) THEN
            H = w2 + thpq(km)*dt*hf(k,j,i)/hold(km)
         ELSE
            H = w2
         ENDIF

         cx = vxr/sqrt(vxr**2+vxi**2)
         sx = vxi/sqrt(vxr**2+vxi**2)
         cy = vyr/sqrt(vyr**2+vyi**2)
         sy = vyi/sqrt(vyr**2+vyi**2)
         cz = vzr/sqrt(vzr**2+vzi**2)
         sz = vzi/sqrt(vzr**2+vzi**2)

         IF ( sx.gt.0. ) THEN
           tx = acos(cx)
         ELSE
           tx = 2.0*pi - acos(cx)
         ENDIF
         IF ( sy.gt.0. ) THEN
           ty = acos(cy)
         ELSE
           ty = 2.0*pi - acos(cy)
         ENDIF
         IF ( sz.gt.0. ) THEN
           tz = acos(cz)
         ELSE
           tz = 2.0*pi - acos(cz)
         ENDIF
         txy = ty-tx
         tzx = tx-tz
         tyz = tz-ty

         tmp1 = 1.
         tmp2 = 1.
         vy2 = ((E*sqrt(1.-H*H/(0.25*ka2(k,j,i)*E*E)))/   &
             (v2*sqrt(1.-w2*w2/(0.25*ka2(k,j,i)*v2*v2)))) &
             * (vyr**2+vyi**2-0.5*v2) + 0.5*E
         tmp1 = sign(tmp1,vy2)
         my  = sqrt(vy2)
         tmp3= E-vy2
         tmp2= sign(tmp2,tmp3)
         tmp = 1./(kx2+kz2)
         mx  = sqrt(abs(tmp3)*kz2*tmp)
         mz  = sqrt(abs(tmp3)*kx2*tmp)

         sxy = ka(k)*H
         dsxy= (kx2+kz2)*mx*my
         syz = ka(i)*H
         dsyz= (kx2+kz2)*my*mz

         IF ( tmp1.ge.0. .and. tmp2.ge.0. .and. abs(sxy).le.abs(dsxy) &
         .and. dsxy.ne.0. .and. abs(syz).le.abs(dsyz) .and. dsyz.ne.0. ) THEN
            sxy = sxy/dsxy
            syz = syz/dsyz
            cxy = sqrt(1.-sxy*sxy)
            cyz = sqrt(1.-syz*syz)
            IF ( cos(txy).lt.0. ) THEN
            cxy = -cxy
            ENDIF
            IF ( cos(tyz).lt.0. ) THEN
            cyz = -cyz
            ENDIF
            szx = -sxy*cyz - syz*cxy
            czx =  cxy*cyz - sxy*syz
            IF ( mx.ge.my .and. mx.ge.mz ) THEN
            vxr = mx*cx
            vxi = mx*sx
            tmp = my/mx
            vyr = tmp*( cxy*vxr-sxy*vxi)
            vyi = tmp*( sxy*vxr+cxy*vxi)
            tmp = mz/mx
            vzr = tmp*( szx*vxi+czx*vxr)
            vzi = tmp*( czx*vxi-szx*vxr)
            ENDIF
            IF ( my.ge.mx .and. my.ge.mz ) THEN
            vyr = my*cy
            vyi = my*sy
            tmp = mx/my
            vxr = tmp*( cxy*vyr+sxy*vyi)
            vxi = tmp*(-sxy*vyr+cxy*vyi)
            tmp = mz/my
            vzr = tmp*(-syz*vyi+cyz*vyr)
            vzi = tmp*( cyz*vyi+syz*vyr)
            ENDIF
            IF ( mz.ge.mx .and. mz.ge.my ) THEN
            vzr = mz*cz
            vzi = mz*sz
            tmp = mx/mz
            vxr = tmp*(-szx*vzi+czx*vzr)
            vxi = tmp*( szx*vzr+czx*vzi)
            tmp = my/mz
            vyr = tmp*( cyz*vzr+syz*vzi)
            vyi = tmp*(-syz*vzr+cyz*vzi)
            ENDIF
         ENDIF
         vx(k,j,i) = cmplx(real(vxr),real(vxi))*norm
         vy(k,j,i) = cmplx(real(vyr),real(vyi))*norm
         vz(k,j,i) = cmplx(real(vzr),real(vzi))*norm
      END DO
      END DO

      RETURN
      END SUBROUTINE vcorr_xny0zn 

!*****************************************************************
      SUBROUTINE vcorr_xnynz0(vx,vy,vz,ef,hf,tepq,thpq,eold,hold,dt,hel,ip)
!-----------------------------------------------------------------
!
! Correction for case where kx, ky, !=0,  kz == 0
!
! Parameters
!     vx  : x-velocity (in/out)
!     vy  : y-velocity (in/out)
!     vz  : z-velocity (in/out)
!     ef  : array with the previous energy field
!     hf  : array with the previous helicity field
!     tepq: modal energy transfer function (in)
!     thpq: modal helicity transfer function (in)
!     eold: old energy spectrum
!     hold: old helicity spectrum
!     dt  : delta_t/ord
!     hel : flag (in) s.t.:
!            1: compute helicity quantities
!            0: skip computation of helicity quantities
!
      USE ali
      USE mpivars
      USE grid
      USE kes
      USE var
      IMPLICIT NONE

      COMPLEX, INTENT(INOUT), DIMENSION(n,n,ista:iend) :: vx,vy,vz
      DOUBLE PRECISION, INTENT (IN), DIMENSION(n/2+1)  :: tepq,thpq
      DOUBLE PRECISION, INTENT (IN), DIMENSION(n/2+1)  :: eold,hold
      DOUBLE PRECISION    :: E,H,mx,my,mz,v2,vx2,vy2,vz2,w2
      DOUBLE PRECISION    :: szx,sxy,syz,czx,cxy,cyz
      DOUBLE PRECISION    :: dszx,dsxy,dsyz,dczx,dcxy,dcyz
      DOUBLE PRECISION    :: sx,sy,sz,cx,cy,cz,tx,ty,tz,tzx,txy,tyz
      DOUBLE PRECISION    :: tmp,tmp1,tmp2,tmp3,inorm
      DOUBLE PRECISION    :: vxr,vxi,vyr,vyi,vzr,vzi
      DOUBLE PRECISION    :: kx2,ky2,kz2
      REAL, INTENT(IN), DIMENSION(n,n,ista:iend) :: ef,hf
      REAL, INTENT(IN)    :: dt
      REAL                :: norm
      INTEGER, INTENT(IN) :: hel,ip
      INTEGER             :: i,j,k,km

      k = 1
      kz2= ka(k)**2
      norm = float(n)**3
      inorm = 1./dble(norm)
      DO i = ista+ip,iend
      kx2= ka(i)**2
      DO j = 2,n     
         IF (ka2(k,j,i).gt.kmax) CYCLE
         ky2= ka(j)**2
         km = int(sqrt(ka2(k,j,i))+.501)        
         vxr = real (vx(k,j,i))*inorm
         vxi = aimag(vx(k,j,i))*inorm
         vyr = real (vy(k,j,i))*inorm
         vyi = aimag(vy(k,j,i))*inorm
         vzr = real (vz(k,j,i))*inorm
         vzi = aimag(vz(k,j,i))*inorm
         v2 = vxr**2+vxi**2+vyr**2+vyi**2+vzr**2+vzi**2
         E  = v2 + 2.0*tepq(km)*dt*ef(k,j,i)/eold(km)
         w2 = ka(k)*( vxr*vyi - vxi*vyr ) &
            + ka(j)*( vzr*vxi - vzi*vxr ) &
            + ka(i)*( vyr*vzi - vyi*vzr )  
         IF ( hel.eq.1 ) THEN
            H = w2 + thpq(km)*dt*hf(k,j,i)/hold(km)
         ELSE
            H = w2
         ENDIF

         cx = vxr/(sqrt(vxr**2+vxi**2) + tinyf)
         sx = vxi/(sqrt(vxr**2+vxi**2) + tinyf)
         cy = vyr/(sqrt(vyr**2+vyi**2) + tinyf)
         sy = vyi/(sqrt(vyr**2+vyi**2) + tinyf)
         cz = vzr/(sqrt(vzr**2+vzi**2) + tinyf)
         sz = vzi/(sqrt(vzr**2+vzi**2) + tinyf)

         IF ( sx.gt.0. ) THEN
           tx = acos(cx)
         ELSE
           tx = 2.0*pi - acos(cx)
         ENDIF
         IF ( sy.gt.0. ) THEN
           ty = acos(cy)
         ELSE
           ty = 2.0*pi - acos(cy)
         ENDIF
         IF ( sz.gt.0. ) THEN
           tz = acos(cz)
         ELSE
           tz = 2.0*pi - acos(cz)
         ENDIF
         txy = ty-tx
         tzx = tx-tz
         tyz = tz-ty

         tmp1 = 1.
         tmp2 = 1.
         vz2 = ((E*sqrt(1.-H*H/(0.25*ka2(k,j,i)*E*E)))/   &
             (v2*sqrt(1.-w2*w2/(0.25*ka2(k,j,i)*v2*v2)))) &
             * ((vzr**2+vzi**2)-0.5*v2) + 0.5*E
         tmp1= sign(tmp1,vz2)
         mz  = sqrt(abs(vz2))
         tmp3= E-vz2
         tmp2= sign(tmp2,tmp3)
         tmp = 1./(kx2+ky2)
         mx  = sqrt(abs(tmp3)*ky2*tmp)
         my  = sqrt(abs(tmp3)*kx2*tmp)

         szx = ka(j)*H
         dszx= (kx2+ky2)*mx*mz
         syz = ka(i)*H
         dsyz= (kx2+ky2)*my*mz

         IF ( tmp1.ge.0. .and. tmp2.ge.0. .and. abs(szx).le.abs(dszx) &
         .and. dszx.ne.0. .and. abs(syz).le.abs(dsyz) .and. dsyz.ne.0. ) THEN
            szx = szx/dszx
            syz = syz/dsyz
            czx = sqrt(1.-szx*szx)
            cyz = sqrt(1.-syz*syz)
            IF ( cos(tzx).lt.0. ) THEN
            czx = -czx
            ENDIF
            IF ( cos(tyz).lt.0. ) THEN
            cyz = -cyz
            ENDIF
            sxy = -syz*czx - szx*cyz
            cxy =  cyz*czx - syz*szx
            IF ( mx.ge.my .and. mx.ge.mz ) THEN
            vxr = mx*cx
            vxi = mx*sx
            tmp = my/mx
            vyr = tmp*(cxy*vxr-sxy*vxi)
            vyi = tmp*(sxy*vxr+cxy*vxi)
            tmp = mz/mx
            vzr = tmp*(szx*vxi+czx*vxr)
            vzi = tmp*(czx*vxi-szx*vxr)
            ENDIF
            IF ( my.ge.mx .and. my.ge.mz ) THEN
            vyr = my*cy
            vyi = my*sy
            tmp = mx/my
            vxr = tmp*( cxy*vyr+sxy*vyi)
            vxi = tmp*(-sxy*vyr+cxy*vyi)
            tmp = mz/my
            vzr = tmp*(-syz*vyi+cyz*vyr)
            vzi = tmp*( cyz*vyi+syz*vyr)
            ENDIF
            IF ( mz.ge.mx .and. mz.ge.my ) THEN
            vzr = mz*cz
            vzi = mz*sz
            tmp = mx/mz
            vxr = tmp*(-szx*vzi+czx*vzr)
            vxi = tmp*( szx*vzr+czx*vzi)
            tmp = my/mz
            vyr = tmp*( cyz*vzr+syz*vzi)
            vyi = tmp*(-syz*vzr+cyz*vzi)
            ENDIF
         ENDIF
         vx(k,j,i) = cmplx(real(vxr),real(vxi))*norm
         vy(k,j,i) = cmplx(real(vyr),real(vyi))*norm
         vz(k,j,i) = cmplx(real(vzr),real(vzi))*norm
      END DO
      END DO

      RETURN
      END SUBROUTINE vcorr_xnynz0

!*****************************************************************
      SUBROUTINE vcorr_xnynzn(vx,vy,vz,ef,hf,tepq,thpq,eold,hold,dt,hel,ip)
!-----------------------------------------------------------------
!
! Correction for case where kx, ky, kz != 0
!
! Parameters
!     vx  : x-velocity (in/out)
!     vy  : y-velocity (in/out)
!     vz  : z-velocity (in/out)
!     ef  : array with the previous energy field
!     hf  : array with the previous helicity field
!     tepq: modal energy transfer function (in)
!     thpq: modal helicity transfer function (in)
!     eold: old energy spectrum
!     hold: old helicity spectrum
!     dt  : delta_t/ord
!     hel : flag (in) s.t.:
!            1: compute helicity quantities
!            0: skip computation of helicity quantities
!
      USE ali
      USE mpivars
      USE grid
      USE kes
      USE var
      IMPLICIT NONE

      COMPLEX, INTENT(INOUT), DIMENSION(n,n,ista:iend) :: vx,vy,vz
      DOUBLE PRECISION, INTENT (IN), DIMENSION(n/2+1)  :: tepq,thpq
      DOUBLE PRECISION, INTENT (IN), DIMENSION(n/2+1)  :: eold,hold
      DOUBLE PRECISION    :: E,H,mx,my,mz,v2,vx2,vy2,vz2,w2
      DOUBLE PRECISION    :: szx,sxy,syz,czx,cxy,cyz
      DOUBLE PRECISION    :: dszx,dsxy,dsyz,dczx,dcxy,dcyz
      DOUBLE PRECISION    :: cx,sx,cy,sy,cz,sz
      DOUBLE PRECISION    :: a2,delta,tmp,vxr,vxi,vyr,vyi,vzr,vzi
      DOUBLE PRECISION    :: kx2,ky2,kz2,tmp1,tmp2,tmp3,tmp4,tmq,inorm
      REAL, INTENT(IN), DIMENSION(n,n,ista:iend) :: ef,hf
      REAL, INTENT(IN)    :: dt
      REAL                :: norm
      INTEGER, INTENT(IN) :: hel,ip
      INTEGER             :: i,j,k,km

      norm = float(n)**3
      inorm = 1./dble(norm)
      DO i = ista+ip,iend
      kx2= ka(i)**2
      DO j = 2,n     
      ky2= ka(j)**2
      DO k = 2,n
         IF (ka2(k,j,i).gt.kmax) CYCLE
         kz2= ka(k)**2
         km = int(sqrt(ka2(k,j,i))+.501)        
         vxr = real (vx(k,j,i))*inorm
         vxi = aimag(vx(k,j,i))*inorm
         vyr = real (vy(k,j,i))*inorm
         vyi = aimag(vy(k,j,i))*inorm
         vzr = (-ka(i)*vxr-ka(j)*vyr)/ka(k)    !Divergence-free
         vzi = (-ka(i)*vxi-ka(j)*vyi)/ka(k)    !Divergence-free
         v2 = vxr**2+vxi**2+vyr**2+vyi**2+vzr**2+vzi**2
         E  = v2 + 2.0*tepq(km)*dt*ef(k,j,i)/eold(km)
         w2 = ka(k)*( vxr*vyi - vxi*vyr ) &
            + ka(j)*( vzr*vxi - vzi*vxr ) &
            + ka(i)*( vyr*vzi - vyi*vzr )  
         IF ( hel.eq.1 ) THEN
            H = w2 + thpq(km)*dt*hf(k,j,i)/hold(km)
         ELSE
            H = w2
         ENDIF

         cx = vxr/sqrt(vxr**2+vxi**2)
         sx = vxi/sqrt(vxr**2+vxi**2)
         cy = vyr/sqrt(vyr**2+vyi**2)
         sy = vyi/sqrt(vyr**2+vyi**2)
         cz = vzr/sqrt(vzr**2+vzi**2)
         sz = vzi/sqrt(vzr**2+vzi**2)

         tmp1= 1.
         tmp2= 1.
         tmp3= 1.
         tmq = (E*sqrt(1.-H**2/(0.25*ka2(k,j,i)*E**2)))/        &
               (v2*sqrt(1-w2**2/(0.25*ka2(k,j,i)*v2**2)))
         vx2 = tmq*(vxr**2+vxi**2-(1.-kx2/ka2(k,j,i))*0.5*v2) + &
             (1.-kx2/ka2(k,j,i))*0.5*E
         vy2 = tmq*(vyr**2+vyi**2-(1.-ky2/ka2(k,j,i))*0.5*v2) + &
             (1.-ky2/ka2(k,j,i))*0.5*E;
         vz2 = tmq*(vzr**2+vzi**2-(1.-kz2/ka2(k,j,i))*0.5*v2) + &
             (1.-kz2/ka2(k,j,i))*0.5*E
         tmp1= sign(tmp1,vx2)
         tmp2= sign(tmp2,vy2)
         tmp3= sign(tmp3,vz2)
         mx  = sqrt(abs(vx2))
         my  = sqrt(abs(vy2))
         mz  = sqrt(abs(vz2))
         szx = H*ka(j)
         dszx= ka2(k,j,i)*mx*mz
         sxy = H*ka(k)
         dsxy= ka2(k,j,i)*my*mx
         syz = H*ka(i)
         dsyz= ka2(k,j,i)*my*mz

         czx = (ky2*E - mx*mx*(kx2+ky2) - mz*mz*(ky2+kz2))
         dczx= 2.*ka(i)*ka(k)*mx*mz

         cxy = (kz2*E - mx*mx*(kx2+kz2) - my*my*(ky2+kz2))
         dcxy= 2.*ka(i)*ka(j)*mx*my

         cyz = (kx2*E - my*my*(kx2+ky2) - mz*mz*(kx2+kz2))
         dcyz= 2.*ka(j)*ka(k)*my*mz

         IF ( tmp1.ge.0. .and. tmp2.ge.0. .and. tmp3.ge.0. &
        .and. abs(sxy).le.abs(dsxy) .and. dsxy.ne.0. .and. &
         abs(szx).le.abs(dszx) .and. dszx.ne.0. &
        .and. abs(syz).le.abs(dsyz) .and. dsyz.ne.0. .and. &
         abs(cxy).le.abs(dcxy) .and. dcxy.ne.0. &
        .and. abs(cyz).le.abs(dcyz) .and. dcyz.ne.0. .and. &
         abs(czx).le.abs(dczx) .and. dczx.ne.0. ) THEN
           szx = szx / dszx
           sxy = sxy / dsxy
           syz = syz / dsyz
           czx = czx / dczx
           cxy = cxy / dcxy
           cyz = cyz / dcyz
           IF ( mx.ge.my .and. mx.ge.mz ) THEN
              vxr = mx*cx
              vxi = mx*sx
              tmq = my/mx
              vyr = tmq*(cxy*vxr-sxy*vxi)
              vyi = tmq*(sxy*vxr+cxy*vxi)
              tmq = mz/mx
              vzr = tmq*(szx*vxi+czx*vxr)
              vzi = tmq*(czx*vxi-szx*vxr)
           ENDIF
           IF ( my.gt.mx .and. my.ge.mz ) THEN
              vyr = my*cy
              vyi = my*sy
              tmq = mx/my
              vxr = tmq*( cxy*vyr+sxy*vyi)
              vxi = tmq*(-sxy*vyr+cxy*vyi)
              tmq = mz/my
              vzr = tmq*(-syz*vyi+cyz*vyr)
              vzi = tmq*( cyz*vyi+syz*vyr)
           ENDIF
           IF ( mz.gt.mx .and. mz.gt.my ) THEN
              vzr = mz*cz
              vzi = mz*sz
              tmq = mx/mz
              vxr = tmq*(-szx*vzi+czx*vzr)
              vxi = tmq*( szx*vzr+czx*vzi)
              tmq = my/mz
              vyr = tmq*( cyz*vzr+syz*vzi)
              vyi = tmq*(-syz*vzr+cyz*vzi)
           ENDIF
         ELSE
           tmp = (1.-kx2/ka2(k,j,i))*0.5*E*(1.+sqrt(1.-H**2/(0.25*E**2* &
                  ka2(k,j,i))))
           tmq = (1.-kx2/ka2(k,j,i))*0.5*E*(1.-sqrt(1.-H**2/(0.25*E**2* &
                  ka2(k,j,i))))
           IF ( vx2.gt.tmp) vx2 = tmp
           IF ( vx2.lt.tmq) vx2 = tmq
           tmp1  = 1.
           tmp1  = sign(tmp1,vx2)
           mx    = sqrt(abs(vx2))
           a2    = 2.*(ky2+kz2)**2
           tmq   = 1./a2
           tmp   = 1.
           tmp2  = 1.
           tmp3  = 1.
           tmp4  = kx2*ky2*kz2*(-mx**4*ka2(k,j,i) + mx*mx*(ky2 + kz2)*E &
              - H*H*(ky2+kz2)**2/(ka2(k,j,i)**2))
           tmp   = sign(tmp ,tmp4)
           delta = 4.*sqrt(abs(tmp4))*tmq
           tmp4  = (-vx2*(2.*(kx2+kz2)*(ky2+kz2)-4.*kx2*ky2) &
              + 2.*kz2*E*(ky2+kz2))*tmq + delta
           tmp2  = sign(tmp2,tmp4)
           my    = sqrt(abs(tmp4))
           tmp4  = (-vx2*(2.*(kx2+ky2)*(ky2+kz2)-4.*kx2*kz2) &
              + 2.*ky2*E*(ky2+kz2))*tmq - delta
           tmp3  = sign(tmp3,tmp4)
           mz    = sqrt(abs(tmp4))

           szx = H*ka(j)
           dszx= ka2(k,j,i)*mx*mz

           sxy = H*ka(k)
           dsxy= ka2(k,j,i)*my*mx

           syz = H*ka(i)
           dsyz= ka2(k,j,i)*my*mz

           czx = (ky2*E - mx*mx*(kx2+ky2) - mz*mz*(ky2+kz2))
           dczx= 2.*ka(i)*ka(k)*mx*mz

           cxy = (kz2*E - mx*mx*(kx2+kz2) - my*my*(ky2+kz2))
           dcxy= 2.*ka(i)*ka(j)*mx*my

           cyz = (kx2*E - my*my*(kx2+ky2) - mz*mz*(kx2+kz2))
           dcyz= 2.*ka(j)*ka(k)*my*mz
           IF ( tmp.ge.0. .and. tmp1.ge.0. .and. tmp2.ge.0. .and. tmp3.ge.0. &
           .and.abs(sxy).le.abs(dsxy) .and. dsxy.ne.0. .and. &
            abs(szx).le.abs(dszx) .and. dszx.ne.0. &
           .and.abs(syz).le.abs(dsyz) .and. dsyz.ne.0. .and. &
            abs(czx).le.abs(dczx) .and. dczx.ne.0. &
           .and.abs(cxy).le.abs(dcxy) .and. dcxy.ne.0. .and. &
            abs(cyz).le.abs(dcyz) .and. dcyz.ne.0. ) THEN
              szx = szx / dszx
              sxy = sxy / dsxy
              syz = syz / dsyz
              czx = czx / dczx
              cxy = cxy / dcxy
              cyz = cyz / dcyz
              IF ( mx.ge.my .and. mx.ge.mz ) THEN
              vxr = mx*cx
              vxi = mx*sx
              tmq = my/mx
              vyr = tmq*(cxy*vxr-sxy*vxi)
              vyi = tmq*(sxy*vxr+cxy*vxi)
              tmq = mz/mx
              vzr = tmq*(szx*vxi+czx*vxr)
              vzi = tmq*(czx*vxi-szx*vxr)
              ENDIF
              IF ( my.ge.mx .and. my.ge.mz ) THEN
              vyr = my*cy
              vyi = my*sy
              tmq = mx/my
              vxr = tmq*( cxy*vyr+sxy*vyi)
              vxi = tmq*(-sxy*vyr+cxy*vyi)
              tmq = mz/my
              vzr = tmq*(-syz*vyi+cyz*vyr)
              vzi = tmq*( cyz*vyi+syz*vyr)
              ENDIF
              IF ( mz.ge.mx .and. mz.ge.my ) THEN
              vzr = mz*cz
              vzi = mz*sz
              tmq = mx/mz
              vxr = tmq*(-szx*vzi+czx*vzr)
              vxi = tmq*( szx*vzr+czx*vzi)
              tmq = my/mz
              vyr = tmq*( cyz*vzr+syz*vzi)
              vyi = tmq*(-syz*vzr+cyz*vzi)
              ENDIF
           ENDIF
         ENDIF
         vx(k,j,i) = cmplx(real(vxr),real(vxi))*norm
         vy(k,j,i) = cmplx(real(vyr),real(vyi))*norm
         vz(k,j,i) = cmplx(real(vzr),real(vzi))*norm
      END DO
      END DO
      END DO

      RETURN
      END SUBROUTINE vcorr_xnynzn

!*****************************************************************
      SUBROUTINE edqnmout(tv,ntv,time,fn)
!-----------------------------------------------------------------
!
! Print double prec. array tv to file
!
! Parameters
!     tv    : EDQNM spectrum or visc 
!     ntv   : no. tv array members
!     time  : time
!     fn    : filename
!
      USE ali
      USE mpivars
      USE grid
      USE kes
      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN), DIMENSION(*) :: tv
      REAL   , INTENT(IN)      :: time
      INTEGER, INTENT(IN)      :: ntv
      INTEGER                  :: j
      CHARACTER(*), INTENT(IN) :: fn

      IF (myrank.eq.0) THEN
         OPEN(1,file=fn)
         DO j = 1,ntv
            WRITE(1,*) tv(j)
         END DO
         CLOSE(1)
      ENDIF

      RETURN
      END SUBROUTINE edqnmout
