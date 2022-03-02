***
      SUBROUTINE evolv1(kw,mass,mt,r,lum,mc,rc,menv,renv,ospin,
     &                  epoch,tm,tphys,tphysf,dtp,z,zpars,
* A. Tanikawa adds here for merger inflation 21/10/14
     &                  tmrg,tcon,rmrg,emrg,
*
* A. Tanikawa corrects reinitialization 21/08/27
     &     reinit,ip,jp,ffb,typesn,arem,
     &     catcherr,catchsn,catchnr,
     &     msec,dmsec,delv)
*
c-------------------------------------------------------------c
c
c     Evolves a single star.
c     Mass loss is an option.
c     The timestep is not constant but determined by certain criteria.
c     Plots the HRD and variables as a function of time.
c
c     Written by Jarrod Hurley 26/08/97 at the Institute of
c     Astronomy, Cambridge.
c
c-------------------------------------------------------------c
c
c     STELLAR TYPES - KW
c
c        0 - deeply or fully convective low mass MS star
c        1 - Main Sequence star
c        2 - Hertzsprung Gap
c        3 - First Giant Branch
c        4 - Core Helium Burning
c        5 - First Asymptotic Giant Branch
c        6 - Second Asymptotic Giant Branch
c        7 - Main Sequence Naked Helium star
c        8 - Hertzsprung Gap Naked Helium star
c        9 - Giant Branch Naked Helium star
c       10 - Helium White Dwarf
c       11 - Carbon/Oxygen White Dwarf
c       12 - Oxygen/Neon White Dwarf
c       13 - Neutron Star
c       14 - Black Hole
c       15 - Massless Supernova
c
c-------------------------------------------------------------c
      use iso_c_binding
*
      implicit none
*
      integer kw,it,ip,jp,j,kwold,rflag
      integer nv
      parameter(nv=50000)
*
      real*8 mass,z,aj
      real*8 epoch,tphys,tphys2,tmold,tbgold
      real*8 mt,tm,tn,tphysf,dtp,tsave
      real*8 tscls(20),lums(10),GB(10),zpars(20)
      real*8 r,lum,mc,teff,rc,menv,renv,vs(3)
      real*8 ospin,jspin,djt,djmb,k2,k3
      parameter(k3=0.21d0)
      real*8 m0,r1,lum1,mc1,rc1,menv1,renv1,k21
      real*8 dt,dtm,dtr,dr,dtdr,dms,dml,mt2,rl
      real*8 tol,tiny
      parameter(tol=1.0d-10,tiny=1.0d-14)
      real*8 ajhold,rm0,eps,alpha2
      parameter(eps=1.0d-06,alpha2=0.09d0)
      real*8 mlwind,vrotf,mlwindenhanced
      external mlwind,vrotf,mlwindenhanced
      logical iplot,isave
      REAL*8 neta,bwind,hewind,mxns,betaacc
      COMMON /VALUE1/ neta,bwind,hewind,mxns,betaacc
      REAL*8 pts1,pts2,pts3
      COMMON /POINTS/ pts1,pts2,pts3
      REAL scm(50000,15),spp(20,3)
      COMMON /SINGLE/ scm,spp
      include 'cppinterface.h'
      real*8 zwind
* Tanikawa adds here 21/02/19
      logical NewStarModel,WindEnhanced,RadiusShrinkage
      common /SingleFlags/ NewStarModel,WindEnhanced,RadiusShrinkage
*
* Tanikawa adds here 20/05/07
      real*8,dimension(2) :: cinc
      real*8,dimension(3) :: nov1
      real*8 kwoldpatch
      real*8 dmstmp
      real*8 ffb,typesn
      real*8 arem
      real*8 tmrg,tcon,rmrg,emrg
*
* A. Tanikawa corrects reinitialization 21/08/27
      logical reinit
      real*8 trsrt
      integer catcherr,catchsn,catchnr
      real*8 msec,dmsec
      real*8,dimension(3) :: delv,delj,dele
*
      dtm = 0.d0
      r = 0.d0
      lum = 0.d0
      mc = 0.d0
      mc1 = 0.d0
      rc = 0.d0
      rl = 0.d0
      if(ospin.le.0.d0)then
         ospin = 1.0d-10
         jspin = 1.0d-10
      endif
      k2 = 0.15d0
      rflag = 0
*
* Setup variables which control the output (if it is required).
*
* A. Tanikawa corrects reinitialization 21/08/27
      if(reinit) then
         trsrt = tphys
         ip  = 0
         ffb = 0.
         catcherr = 0
         catchsn  = 0
         catchnr  = 0
         msec     = mass
         dmsec    = 0.
         delv(1)  = 0.
         delv(2)  = 0.
         delv(3)  = 0.
      else
         ip = 0
         jp = 0
      endif
*
      tsave = tphys
      isave = .true.
      iplot = .false.
      if(dtp.le.0.d0)then
         iplot = .true.
         isave = .false.
         tsave = tphysf
      elseif(dtp.gt.tphysf)then
         isave = .false.
         tsave = tphysf
      endif
* 
      do 10 , j = 1,nv
*
         if(neta.gt.tiny.and.j.gt.1)then
*
* Calculate mass loss from the previous timestep.
*
            dt = 1.0d+06*dtm
            if(askInUseOrNot()) then
!               zwind = getMetallicity()
               zwind = getWindMetallicity()
               dms = mlwind(kw,lum,r,mt,mc,rl,zwind)
* Tanikawa modifies here for wind 21/02/19
!               dms = mlwindenhanced(kw,lum,r, mt,mc,ospin,dms)*dt
               if(WindEnhanced) then
                  dms = mlwindenhanced(kw,lum,r, mt,mc,ospin,dms)*dt
               else
                  dms = dms*dt
               endif
*
* A. Tanikawa observe wind mass loss temporarily 20/02/18
!               dms = mlwind(kw,lum,r,mt,mc,rl,zwind)
!               write(*,*)'wind',tphys,dms,kw
!               dms = dms*dt
*
            else
               dms = mlwind(kw,lum,r,mt,mc,rl,z)*dt
* A. Tanikawa observe wind mass loss temporarily 20/02/18
!               dms = mlwind(kw,lum,r,mt,mc,rl,zwind)
!               write(*,*)'wind',tphys,dms,kw
!               dms = dms*dt
*
            endif
            if(kw.lt.10)then
               dml = mt - mc
               if(dml.lt.dms)then
                  dtm = (dml/dms)*dtm
                  dms = dml
               endif
            endif
         else
            dms = 0.d0
         endif
*
* Limit to 1% mass loss.
*
         if(dms.gt.0.01d0*mt)then
            dtm = 0.01d0*mt*dtm/dms
            dms = 0.01d0*mt
         endif
*
* Calculate the rate of angular momentum loss due to magnetic braking 
* and/or mass loss.
*
         if(j.gt.1)then
            djt = (2.d0/3.d0)*(dms/(1.0d+06*dtm))*r*r*ospin
            if(mt.gt.0.35d0.and.kw.lt.10)then
               djmb = 5.83d-16*menv*(r*ospin)**3/mt
               djt = djt + djmb
            endif
         endif
*
* Update mass and time and reset epoch for a MS (and possibly a HG) star.
*
         if(dms.gt.0.d0)then
            mt = mt - dms
            if(kw.le.2.or.kw.eq.7)then
               m0 = mass
               mc1 = mc
               mass = mt
               tmold = tm
               tbgold = tscls(1)
               CALL star(kw,mass,mt,tm,tn,tscls,lums,GB,zpars)
               if(kw.eq.2)then
                  if(GB(9).lt.mc1.or.m0.gt.zpars(3))then
                     mass = m0
                  else
                     epoch = tm + (tscls(1) - tm)*(ajhold-tmold)/
     &                            (tbgold - tmold)
                     epoch = tphys - epoch
                  endif
               else
                  epoch = tphys - ajhold*tm/tmold
               endif
            endif
         endif
         tphys2 = tphys
         tphys = tphys + dtm
*
* Find the landmark luminosities and timescales as well as setting
* the GB parameters.
*
         aj = tphys - epoch
         CALL star(kw,mass,mt,tm,tn,tscls,lums,GB,zpars)
*
* Find the current radius, luminosity, core mass and stellar type
* given the initial mass, current mass, metallicity and age
*
         CALL hrdiag(mass,aj,mt,tm,tn,tscls,lums,GB,zpars,
     &        r,lum,kw,mc,rc,menv,renv,k2,ffb,typesn,
     &        jspin,arem,
* A. Tanikawa adds here for merger inflation 21/10/14
     &        tphys,tmrg,tcon,rmrg)
*
*
* If mass loss has occurred and no type change then check that we
* have indeed limited the radius change to 10%.
*
         if(kw.eq.kwold.and.dms.gt.0.d0.and.rflag.ne.0)then
            mt2 = mt + dms
            dml = dms/dtm
            it = 0
 20         dr = r - rm0
            if(ABS(dr).gt.0.1d0*rm0)then
               it = it + 1
               if(it.eq.20.and.kw.eq.4) goto 30
               if(it.gt.30)then
                  WRITE(99,*)' DANGER1! ',it,kw,mass,dr,rm0
                  WRITE(*,*)' STOP: EVOLV1 FATAL ERROR '
* A. Tanikawa corrects reinitialization 21/08/27
                  if(reinit) then
                     catcherr = 1
                  else
                     CALL exit(0)
                     STOP 
                  endif
*
               endif
               dtdr = dtm/ABS(dr)
               dtm = alpha2*MAX(rm0,r)*dtdr
               if(it.ge.20) dtm = 0.5d0*dtm
               if(dtm.lt.1.0d-07*aj) goto 30
               dms = dtm*dml
               mt = mt2 - dms
               if(kw.le.2.or.kw.eq.7)then
                  mass = mt
                  CALL star(kw,mass,mt,tm,tn,tscls,lums,GB,zpars)
                  if(kw.eq.2)then
                     if(GB(9).lt.mc1.or.m0.gt.zpars(3))then
                        mass = m0
                     else
                        epoch = tm + (tscls(1) - tm)*(ajhold-tmold)/
     &                               (tbgold - tmold)
                        epoch = tphys2 - epoch
                     endif
                  else
                     epoch = tphys2 - ajhold*tm/tmold
                  endif
               endif
               tphys = tphys2 + dtm
               aj = tphys - epoch
               mc = mc1
               CALL star(kw,mass,mt,tm,tn,tscls,lums,GB,zpars)
               CALL hrdiag(mass,aj,mt,tm,tn,tscls,lums,GB,zpars,
     &              r,lum,kw,mc,rc,menv,renv,k2,ffb,typesn,
     &              jspin,arem,
* A. Tanikawa adds here for merger inflation 21/10/14
     &              tphys,tmrg,tcon,rmrg)
*
               goto 20
            endif
 30         continue
         endif
*
* Initialize or adjust the spin of the star.
*
         if(j.eq.1)then
            if(tphys.lt.tiny.and.ospin.lt.0.001d0)then
               ospin = 45.35d0*vrotf(mt)/r
            endif
            jspin = ospin*(k2*r*r*(mt-mc)+k3*rc*rc*mc)
         else
            jspin = MAX(1.0d-10,jspin - djt*1.0d+06*dtm)
            ospin = jspin/(k2*r*r*(mt-mc)+k3*rc*rc*mc)
         endif
*
* Test for changes in evolution type.
*
         if(j.eq.1.or.kw.ne.kwold)then
*
* Force new NS or BH to have a one second period. 
* 
            if(kw.eq.13.or.kw.eq.14)then
* A. Tanikawa corrects reinitialization 21/08/27
               if(reinit .and. j.ne.1) then 
                  catchsn = 1
                  tphysf  = tphys
               endif
*
               ospin = 2.0d+08
               jspin = k3*rc*rc*mc*ospin
               CALL kick(kw,mass,mt,0.d0,0.d0,-1.d0,0.d0,
     &              ffb,vs,cinc,nov1,
* A. Tanikawa corrects reinitialization 21/08/27
     &              delv,delj,dele)
*
            endif
            jp = jp + 1
            spp(jp,1) = tphys
            spp(jp,2) = float(kw)
            if(kw.eq.15)then
* A. Tanikawa corrects reinitialization 21/08/27
               if(reinit) catchnr = 1
*
               spp(jp,3) = mass 
               goto 90
            else
               spp(jp,3) = mt
            endif
* A. Tanikawa corrects reinitialization 21/08/27
            if(reinit .and. jp.gt.1 .and. spp(jp,2).eq.spp(jp-1,2)) then
               jp = jp - 1
            endif
*
         endif
* A. Tanikawa corrects reinitialization 21/08/27
         if(reinit .and. catchsn+catchnr.eq.0) then
            msec = mt
            if(dt.eq.0.) then
               dmsec = 0.
            else
               dmsec = - dms/dt
            endif
            delv(1) = 0.
            delv(2) = 0.
            delv(3) = 0.
         endif
*
*
* Record values for plotting and reset epoch.
*
         epoch = tphys - aj
         if((isave.and.tphys.ge.tsave).or.iplot)then
*
            ip = ip + 1
            scm(ip,1) = tphys
            scm(ip,2) = float(kw)
            scm(ip,3) = mass
            scm(ip,4) = mt
            scm(ip,5) = log10(lum)
            scm(ip,6) = log10(r)
            teff = 1000.d0*((1130.d0*lum/(r**2.d0))**(1.d0/4.d0))
            scm(ip,7) = log10(teff)
            scm(ip,8) = mc
            scm(ip,9) = rc
            scm(ip,10) = menv
            scm(ip,11) = renv
            scm(ip,12) = epoch
            scm(ip,13) = ospin
            scm(ip,14) = k2
            if(isave) tsave = tsave + dtp
            if(tphysf.lt.tiny)then
               ip = ip + 1
               do 35 , it = 1,13
                  scm(ip,it) = scm(ip-1,it)
 35            continue
            endif
         endif
*
         if(tphys.ge.tphysf)then
            jp = jp + 1
            spp(jp,1) = tphys
            spp(jp,2) = float(kw)
            spp(jp,3) = mt
* A. Tanikawa corrects reinitialization 21/08/27
            if(reinit .and. jp.gt.1 .and. spp(jp,2).eq.spp(jp-1,2)) then
               jp = jp - 1
            endif
*
            goto 90
         endif
*
* Record radius and current age.
*
         rm0 = r
         ajhold = aj
         if(kw.ne.kwold) kwold = kw
         CALL deltat(kw,aj,tm,tn,tscls,dtm,dtr)
* A. Tanikawa modifies here 21/03/24
!         if(askInUseOrNot()) then
         if(askInUseOrNot() .and. askInScopeOfApplication(mass)) then ! mt-> mass (21/09/24)
*
            call calcTimestepAGBPhase(kw, aj, mass, tn, pts3, dtm, dtr)
         endif
*
* Check for type change.
*
         it = 0
         m0 = mass
         if((dtr-dtm).le.tol.and.kw.le.9)then
*
* Check final radius for too large a jump.
*
            aj = MAX(aj,aj*(1.d0-eps)+dtr)
            mc1 = mc 
* Tanikawa adds here 20/05702
            kwoldpatch = kw
*
            CALL hrdiag(mass,aj,mt,tm,tn,tscls,lums,GB,zpars,
     &           r1,lum1,kw,mc1,rc1,menv1,renv1,k21,ffb,typesn,
     &           jspin,arem,
* A. Tanikawa adds here for merger inflation 21/10/14
     &           tphys,tmrg,tcon,rmrg)
*
            dr = r1 - rm0
* Tanikawa change here 20/05702
!            if(ABS(dr).gt.0.1d0*rm0)then
            if(ABS(dr).gt.0.1d0*rm0 .and. kwoldpatch.eq.kw)then
*
               dtm = dtr - ajhold*eps
               dtdr = dtm/ABS(dr)
               dtm = alpha2*MAX(r1,rm0)*dtdr
               goto 40
            else
               dtm = dtr
               goto 50
            endif
         endif
*
* Limit to a 10% increase in radius assuming no further mass loss
* and thus that the pertubation functions due to small envelope mass
* will not change the radius.
*
 40      aj = ajhold + dtm
         mc1 = mc 
         CALL hrdiag(mass,aj,mt,tm,tn,tscls,lums,GB,zpars,
     &        r1,lum1,kw,mc1,rc1,menv1,renv1,k21,ffb,typesn,
     &        jspin,arem,
* A. Tanikawa adds here for merger inflation 21/10/14
     &        tphys,tmrg,tcon,rmrg)
*
         dr = r1 - rm0
         it = it + 1
         if(it.eq.20.and.kw.eq.4) goto 50
         if(it.gt.30)then
            WRITE(99,*)' DANGER2! ',it,kw,mass,dr,rm0
            WRITE(*,*)' STOP: EVOLV1 FATAL ERROR '
* A. Tanikawa corrects reinitialization 21/08/27
            if(reinit) then
               catcherr = 1
            else
               CALL exit(0)
               STOP 
            endif
*
         endif
         if(ABS(dr).gt.0.1d0*rm0)then
            dtdr = dtm/ABS(dr)
            dtm = alpha2*MAX(rm0,r1)*dtdr
            if(it.ge.20) dtm = 0.5d0*dtm
            goto 40
         endif
*
 50      continue
*
* Ensure that change of type has not occurred during radius check. 
* This is rare but may occur for HG stars of ZAMS mass > 50 Msun. 
*
         if(kw.ne.kwold)then
            kw = kwold
            mass = m0
            CALL star(kw,mass,mt,tm,tn,tscls,lums,GB,zpars)
         endif
*
* Choose minimum of time-scale and remaining interval (> 100 yrs).
*
         dtm = MAX(dtm,1.0d-07*aj)
         dtm = MIN(dtm,tsave-tphys)
*
 10   continue
*
 90   continue
*
      tphysf = tphys
* Tanikawa modifies here 20/05/12
      scm(ip+1,1) = tphys
      scm(ip+1,2) = 15.
      scm(ip+1,3) = mass
      scm(ip+1,4) = 0.
      scm(ip+1,5) = -10.
      scm(ip+1,6) = -4.
      scm(ip+2,1) = -1.0
!      scm(ip+1,1) = -1.0
*
      spp(jp+1,1) = -1.0
      if(ip.ge.nv)then
* A. Tanikawa corrects reinitialization 21/08/27
         WRITE(99,*)' EVOLV1 ARRAY ERROR ',mass
         WRITE(*,*)' STOP: EVOLV1 ARRAY ERROR '
         if(reinit) then
            catcherr = 1
         else
            CALL exit(0)
            STOP
         endif
*
      endif
*
      RETURN
      END
***
