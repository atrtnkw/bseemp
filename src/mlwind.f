***
      real*8 FUNCTION mlwind(kw,lum,r,mt,mc,rl,z)
      implicit none
      integer kw,mdflag
      real*8 lum,r,mt,mc,rl,z
      real*8 teff,teff1,dml,dms,dmt,p0,x,mew,vw,t40
      real*8 lum0,kap,flbv
      parameter(lum0=7.0d+04,kap=-0.5d0)
      parameter(flbv=1.5d0)
* A. Tanikawa modifies here 20/02/17
!      real*8 neta,bwind
!      parameter(neta=0.5d0)
!      parameter(bwind=0.d0)
      real*8 neta,bwind,hewind,mxns,betaacc
      common /value1/ neta,bwind,hewind,mxns,betaacc
* Tanikawa adds Pulsation-driven wind (Nakauchi et al. 2020, ApJ, 902, 81)
      logical switchMlwindPulsationDriven
      parameter(switchMlwindPulsationDriven=.false.)
      real*8 mlwindPulsationDriven
*
* Calculate stellar wind mass loss.
*
* Set mdflag to determine the mass-loss you want: 
*  = 1 - Hurley, Pols & Tout (2000) SSE basic rates; 
*  = 2 - SSE + LBV added; 
*  = 3 - Belczynski et al. (2010, ApJ, 714, 1217) updates. 
*  = 4 - Belczynski without bi-stability jump. 
* [Currently mdflag = 1 is recommended]
      mdflag = 3
********* Tanikawa found this is not consistent with Belczynski+10 *********
*      mdflag = 4
****************************************************************************
*
      dms = 0.d0
*
      teff = 3.762d0 + 0.25d0*log10(lum) - 0.5d0*log10(r)
      teff = 10.d0**teff
      teff = MIN(teff,50000.d0)
*
      if(lum.gt.4000.d0.and.(mdflag.le.2.or.
     &                      (mdflag.ge.3.and.kw.le.6)))then
* Apply mass loss of Nieuwenhuijzen & de Jager (1990, A&A, 231, 134)
* for massive stars over the entire HRD with a metallicity factor 
* from Kudritzki et al. (1989, A&A, 219, 205). 
         x = MIN(1.d0,(lum-4000.d0)/500.d0)
         dms = 9.6d-15*x*(r**0.81d0)*(lum**1.24d0)*(mt**0.16d0)
         dms = dms*(z/0.02d0)**(1.d0/2.d0)
* Belczynski+2010's eq. (3) noted by A. Tanikawa.
      endif
*
      if(mdflag.ge.3.and.kw.le.6)then
* Apply mass loss for hot, massive H-rich O/B stars following 
* Vink et al. (2001, A&A, 369,574). 
         if(teff.ge.12500.0.and.teff.le.25000.0.and.mdflag.eq.3)then
            teff1 = MIN(teff,22500.d0)
            vw = 1.3d0/2.d0
            dml = -6.688d0 + 2.21d0*log10(lum/1.0d+05) -  
     &             1.339d0*log10(mt/30.d0) - 1.601d0*log10(vw) + 
     &             0.85d0*log10(z/0.02d0) + 1.07d0*log10(teff1/20000.d0)
            dml = 10.d0**dml
            dms = MAX(dms,dml)
* Belczynski+2010's eq. (6) noted by A. Tanikawa.
         elseif(teff.gt.12500.0.and.teff.le.50000.1)then
            teff1 = teff
            if(mdflag.eq.3) teff1 = MAX(teff,27500.d0)
            vw = 2.6d0/2.d0
            t40 = log10(teff1/40000.d0)
            dml = -6.697d0 + 2.194d0*log10(lum/1.0d+05) -  
     &             1.313d0*log10(mt/30.d0) - 1.226d0*log10(vw) + 
     &             0.85d0*log10(z/0.02d0) + 
     &             0.933d0*t40*(1.d0 - 11.704d0*t40)
            dml = 10.d0**dml
            dms = MAX(dms,dml)
* Belczynski+2010's eq. (7) noted by A. Tanikawa.
         endif
      endif
*
* A. Tanikawa modifies here 20/05/07
!      if(kw.ge.2.and.kw.le.9)then
      if(kw.ge.1.and.kw.le.9)then
*
* Standard 'Reimers' mass loss for giants from Kudritzki & Reimers 
* (1978, A&A, 70, 227). 
         dml = neta*4.0d-13*r*lum/mt
* Belczynski+2010's eq. (1) noted by A. Tanikawa.
* Check for any tidally enhanced mass loss in binary systems (optional): 
* see Tout & Eggleton (1988, MNRAS, 231, 823).  
         if(rl.gt.0.d0) dml = dml*(1.d0 + bwind*(MIN(0.5d0,(r/rl)))**6)
*
         if(kw.eq.5.or.kw.eq.6)then
* Apply mass loss of Vassiliadis & Wood (1993, ApJ, 413, 641) 
* for high pulsation periods on AGB.
            p0 = -2.07d0 - 0.9d0*log10(mt) + 1.94d0*log10(r)
            p0 = 10.d0**p0
            p0 = MIN(p0,2000.d0)
            dmt = -11.4d0+0.0125d0*(p0-100.d0*MAX(mt-2.5d0,0.d0))
            dmt = 10.d0**dmt
            dmt = 1.d0*MIN(dmt,1.36d-09*lum)
            dml = MAX(dml,dmt)
* Belczynski+2010's eq. (2) noted by A. Tanikawa.
         endif
*
         dms = MAX(dml,dms)
         if(kw.gt.6)then
* Apply mass loss of Hamann & Koesterke (1998, A&A, 335, 1003) 
* for WR (naked helium) stars. 
* A. Tanikawa modifies here 20/02/17
!            dml = 1.0d-13*lum**(3.d0/2.d0)
            dml = 1.0d-13*hewind*lum**(3.d0/2.d0)
* Belczynski+2010's eq. (9) noted by A. Tanikawa. Z-dependence added later.
*
            if(mdflag.eq.3)then
* Add metallicity factor from Vink & de Koter (2005, A&A, 442, 587). 
               dml = dml*(z/0.02d0)**0.86d0
            endif
            dms = MAX(dml,dms)
         else
            mew = ((mt-mc)/mt)*MIN(5.d0,MAX(1.2d0,(lum/lum0)**kap))
* Apply the reduced WR-like mass loss for small H-envelope mass 
* as described in the Hurley, Pols & Tout (200) SSE paper. 
* A. Tanikawa adds here 20/05/07
            if(kw.le.1)mew=1.
*
            if(mew.lt.1.d0)then
* A. Tanikawa modifies here 20/02/17
!               dml = 1.0d-13*lum**(3.d0/2.d0)*(1.d0 - mew)
               dml = 1.0d-13*hewind*lum**(3.d0/2.d0)
     &              *(z/0.02d0)**(0.86d0)*(1.d0 - mew)
* Belczynski+2010's eq. (9) noted by A. Tanikawa. Modified by H envelope
*
               dms = MAX(dml,dms)
            endif
* LBV-like mass loss beyond the Humphreys-Davidson limit
* (see Humphreys & Davidson 1994 and Belczynski et al. 2010). 
            x = 1.0d-5*r*SQRT(lum)
            if(mdflag.gt.1.and.lum.gt.6.0d+05.and.x.gt.1.d0)then
               if(mdflag.eq.2)then
                  dml = 0.1d0*(x-1.d0)**3*(lum/6.0d+05-1.d0)
* Belczynski+2010's eq. (5) noted by A. Tanikawa. Not adopted.
                  dms = dms + dml
               else
* A. Tanikawa modifies here to MOBSE1 20/03/03
                  dml = 1.0d-04*flbv
!                  dml = 1.0d-04*flbv*(z/0.02d0)**(0.86d0)
*
* Belczynski+2010's eq. (8) noted by A. Tanikawa.
                  dms = MAX(dml,dms) ! original prescription
!                  dms = dms + dml ! kinugawa prescription
*
               endif
            endif
         endif
*
      endif
      if(switchMlwindPulsationDriven .and. kw.lt.7)then
         dml = mlwindPulsationDriven(mt,teff,z)
         dms = max(dml,dms)
      endif
*
      mlwind = dms
*
      return
      end
***

      real*8 function mlwindenhanced(kw,lum,r,mt,mc,ospin,ml0)
      implicit none
      integer kw
      real*8 lum,r,mt,mc,ospin,ml0
      real*8 eddlum,gamma,corlum
      real*8 vrot,vcrit,tkh
      real*8 ml1,ml2
      real*8 gammac
      parameter(gammac=0.99)
!      parameter(gammac=1d2)
      real*8 xfrac
      parameter(xfrac=0.76)

      vrot   = ospin * r
!      vrot   = (1./(1d4/365.))*r
      eddlum = 10**4.813 * mt / (1.+xfrac)
      gamma  = lum / eddlum
* A. Tanikawa modifies here 20/06/01
      if(7.le.kw .and. kw.le.9 .and. gamma .gt. gammac) then
         gamma  = gammac
         corlum = gammac * eddlum
! Grafener et al. (2011), Maeder et al. (2012)
!      if(7.le.kw .and. kw.le.9) then
!         eddlum = 10**6.8*(mt/1e2)
!         corlum = 10**(3.017+2.446*log10(mt)-0.306*(log10(mt))**2)
!         gamma  = corlum / eddlum
*
      else
         corlum = lum
      endif
      vcrit = (3.916d8 * mt * (1. - gamma) / r)**0.5
      ml1   = ml0 * (1. - vrot/vcrit)**(-0.43)

      tkh   = 3.12d+07 * mt / (r*corlum)
      if(kw.le.1 .or. kw.eq.7 .or. kw.eq.10) then
         tkh = tkh * mt
      else
         tkh = tkh * (mt - mc)
      endif
      ml2   = 0.1 * mt / tkh

      if(kw.le.9) then
         mlwindenhanced = min(ml1,ml2)
      else
         mlwindenhanced = ml0
      endif
!      write(*,*)kw,mt,mlwindenhanced,lum/eddlum,vrot/vcrit

      return
      end

***

      real*8 function mlwindPulsationDriven(mt,teff,z)
      implicit none
      real*8 mt,teff,z
      real*8 ialph1,ialph2,ibeta1,ibeta2,igamm
      real*8 ilogdms1,ilogdms2

      if(z.le.6d-7) then ! Z=0Zsun
         ialph1 = 1.6
         ialph2 = 4.03
         ibeta1 = 1.0d4
         ibeta2 = 5.04
         igamm  = 4.
      else if(z.le.6d-6) then ! Z=10^-4Zsun (2e-6)
         ialph1 = 1.4
         ialph2 = 3.65
         ibeta1 = 150.
         ibeta2 = 4.93
         igamm  = 2.
      else if(z.le.6d-5) then ! Z=10^-3Zsun (2e-5)
         ialph1 = 1.4
         ialph2 = 3.55
         ibeta1 = 200.
         ibeta2 = 4.89
         igamm  = 2.
      else if(z.le.6d-4) then ! Z=10^-2Zsun (2e-4)
         ialph1 = 1.65
         ialph2 = 3.4
         ibeta1 = 200.
         ibeta2 = 4.85
         igamm  = 2.
      else if(z.le.6d-3) then ! Z=10^-1Zsun (2e-3)
         ialph1 = 0.0
         ialph2 = 3.5
         ibeta1 = 300.
         ibeta2 = 4.82
         igamm  = 2.
      else
         write(0,*) 'Error: Not support z=', z
         stop
      endif
      ilogdms1 = ialph1*log10(mt*1d-3)-ialph2
     &     -ibeta1*(log10(teff)-ibeta2)**igamm
      ilogdms2 = -2.88+log10(mt*1d-3)-15.6*(log10(teff)-3.7)
!      write(*,*) mt,teff,10**ilogdms1+10**ilogdms2,'hoge'
      mlwindPulsationDriven = 10**ilogdms1+10**ilogdms2
      return
      end
