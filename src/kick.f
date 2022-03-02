***
      SUBROUTINE kick(kw,m1,m1n,m2,ecc,sep,jorb,ffb,vs,cinc,nov1,
* A. Tanikawa corrects reinitialization 21/08/27
     &     delv,delj,dele)
*
      implicit none
*
      integer kw,k
      INTEGER idum
      COMMON /VALUE3/ idum
      INTEGER idum2,iy,ir(32)
      COMMON /RAND3/ idum2,iy,ir
      integer bhflag
      real*8 m1,m2,m1n,ecc,sep,jorb,ecc2
      real*8 pi,twopi,gmrkm,yearsc,rsunkm
      parameter(yearsc=3.1557d+07,rsunkm=6.96d+05)
      real*8 mm,em,dif,der,del,r
      real*8 u1,u2,vk,v(4),s,theta,phi
      real*8 sphi,cphi,stheta,ctheta,salpha,calpha
      real*8 vr,vr2,vk2,vn2,hn2
      real*8 mu,cmu,vs(3),v1,v2,mx1,mx2
      real*8 sigma
      COMMON /VALUE4/ sigma,bhflag
      real ran3,xx
      external ran3
*     A. Tanikawa adds here 20/05/05
      real*8,dimension(2) :: cinc
      real*8,dimension(3) :: nov1,xorb,vorb
      real*8 jold
      real*8 vabs,angl
      real*8 ffb
      real*8 sqdisc
*
* A. Tanikawa corrects reinitialization 21/08/27
      real*8,dimension(3) :: delv,delj,dele
      real*8,dimension(3,2) :: vloc0,vloc1
      real*8,dimension(3) :: jmom0,jmom1,lrl0,lrl1,vtmp
      real*8 minv,mred,gminv,rinv,lrlinv,lrlabs,vunit
      real*8 srotlrl,crotlrl
      real*8 GravConst
      parameter (GravConst=3.916d8)
      real*8 tolecc,tolvel,errvel,delvabs,vsabs
      parameter (tolecc=0.01,tolvel=0.01)
*
*
      do k = 1,3
         vs(k) = 0.d0
      enddo
*     if(kw.eq.14.and.bhflag.eq.0) goto 95
*
      pi = ACOS(-1.d0)
      twopi = 2.d0*pi
* Conversion factor to ensure velocities are in km/s using mass and
* radius in solar units.
      gmrkm = 1.906125d+05
*
* Find the initial separation by randomly choosing a mean anomaly.
      if(sep.gt.0.d0.and.ecc.ge.0.d0)then
         xx = RAN3(idum)
         mm = xx*twopi
         em = mm
 2       dif = em - ecc*SIN(em) - mm
         if(ABS(dif/mm).le.1.0d-04) goto 3
         der = 1.d0 - ecc*COS(em)
         del = dif/der
         em = em - del
         goto 2
 3       continue
         r = sep*(1.d0 - ecc*COS(em)) ! mean anomaly: u
*
* Find the initial relative velocity vector.
         salpha = SQRT((sep*sep*(1.d0-ecc*ecc))/(r*(2.d0*sep-r)))
         calpha = (-1.d0*ecc*SIN(em))/SQRT(1.d0-ecc*ecc*COS(em)*COS(em))
         vr2 = gmrkm*(m1+m2)*(2.d0/r - 1.d0/sep)
         vr = SQRT(vr2)
      else
         vr = 0.d0
         vr2 = 0.d0
         salpha = 0.d0
         calpha = 0.d0
      endif
*
* Generate Kick Velocity using Maxwellian Distribution (Phinney 1992).
* Use Henon's method for pairwise components (Douglas Heggie 22/5/97).
      do 20 k = 1,2
         u1 = RAN3(idum)
         u2 = RAN3(idum)
* Generate two velocities from polar coordinates S & THETA.
         s = sigma*SQRT(-2.d0*LOG(1.d0 - u1))
         theta = twopi*u2
         v(2*k-1) = s*COS(theta)
         v(2*k) = s*SIN(theta)
 20   continue
      vk2 = v(1)**2 + v(2)**2 + v(3)**2
      vk = SQRT(vk2)
      if((kw.eq.14.and.bhflag.eq.0).or.kw.lt.0)then
         vk2 = 0.d0
         vk = 0.d0
         if(kw.lt.0) kw = 13
*     A. Tanikawa adds here 20/01/12
      else if(kw.eq.14.and.bhflag.eq.1) then
!         vk  = vk * 1.4 / m1
         vk  = vk * (1. - ffb)
         vk2 = vk * vk
*
*     A. Tanikawa adds here 20/03/03
      else if(kw.eq.14.and.bhflag.eq.2) then
         vk  = vk
         vk2 = vk * vk
*
      endif
      sphi = -1.d0 + 2.d0*u1
      phi = ASIN(sphi)
      cphi = COS(phi)
      stheta = SIN(theta)
      ctheta = COS(theta)
*     WRITE(66,*)' KICK VK PHI THETA ',vk,phi,theta
* A. Tanikawa corrects reinitialization 21/08/27
      delv(1) = vk*ctheta*cphi
      delv(2) = vk*stheta*cphi
      delv(3) = vk       *sphi
*
      if(sep.le.0.d0.or.ecc.lt.0.d0) goto 90
*     A. Tanikawa adds here 20/05/05  ! vr(1)=-vr*salpha, vr(2)=-vr*calpha, calpha=-rxv/(rv)
      jold    = (sep*(1.-ecc))*
     &     sqrt(gmrkm*(m1+m2)*(2.d0/(sep*(1.-ecc))-1.d0/sep))
      vorb(1) = -vr*salpha !-jold/(r*sqrt(1.-ecc**2))*sin(em)
      vorb(2) = -vr*calpha ! jold/ r                 *cos(em)
      vorb(3) = 0.
      sqdisc   = sqrt(1./(salpha)**2-1.) !disc = vorb(2)**2*jold**2/vr**4-(jold**2-(r*vorb(1))**2)/vr**2
      xorb(1) = (+vorb(2)+vorb(1)*sqdisc)*jold/vr**2 !vorb(2)*jold/vr**2+sqrt(disc)
      xorb(2) = (-vorb(1)+vorb(2)*sqdisc)*jold/vr**2 !(xorb(1)*vorb(2)-jold)/vorb(1)
      xorb(3) = 0.
      if((em-pi)*(xorb(1)*vorb(1)+xorb(2)*vorb(2))>0.) then
         xorb(1) = (+vorb(2)-vorb(1)*sqdisc)*jold/vr**2 !vorb(2)*jold/vr**2-sqrt(disc)
         xorb(2) = (-vorb(1)-vorb(2)*sqdisc)*jold/vr**2 !(xorb(1)*vorb(2)-jold)/vorb(1)
         xorb(3) = 0.
      endif
* A. Tanikawa corrects reinitialization 21/08/27
      minv       = 1./(m1+m2)
      vloc0(1,1) = vorb(1)*(+m2)*minv
      vloc0(2,1) = vorb(2)*(+m2)*minv
      vloc0(3,1) = vorb(3)*(+m2)*minv
      vloc0(1,2) = vorb(1)*(-m1)*minv
      vloc0(2,2) = vorb(2)*(-m1)*minv
      vloc0(3,2) = vorb(3)*(-m1)*minv
      mred       = m1*m2*minv
      vunit      = yearsc/rsunkm
      jmom0(1)   = mred*(xorb(2)*vorb(3)-xorb(3)*vorb(2))*vunit
      jmom0(2)   = mred*(xorb(3)*vorb(1)-xorb(1)*vorb(3))*vunit
      jmom0(3)   = mred*(xorb(1)*vorb(2)-xorb(2)*vorb(1))*vunit
      gminv      = 1./(GravConst*(m1+m2))*vunit/mred
      rinv       = 1./sqrt(xorb(1)**2+xorb(2)**2+xorb(3)**2)
      lrl0(1) = gminv*(vorb(2)*jmom0(3)-vorb(3)*jmom0(2))-xorb(1)*rinv
      lrl0(2) = gminv*(vorb(3)*jmom0(1)-vorb(1)*jmom0(3))-xorb(2)*rinv
      lrl0(3) = gminv*(vorb(1)*jmom0(2)-vorb(2)*jmom0(1))-xorb(3)*rinv
      lrlabs = sqrt(lrl0(1)**2+lrl0(2)**2+lrl0(3)**2)
      lrlinv = 1./lrlabs
      if(abs(lrlabs-ecc).gt.tolecc) then
         write(*,*)'Error: LR vec before SN: ',abs(lrlabs-ecc),ecc
         stop
      endif
      crotlrl =  lrl0(1)/(lrl0(1)**2+lrl0(2)**2)*lrlabs
      srotlrl = -lrl0(2)/(lrl0(1)**2+lrl0(2)**2)*lrlabs
*
* A. Tanikawa adds here for recording kicks 20/07/07 ! vk(1)=-vk*ctheta*cphi, vk(2)=-vk*stheta*cphi, vk(3)=-vk*sphi
      vabs = sqrt(vorb(1)**2+vorb(2)**2+vorb(3)**2)
      angl = (vorb(1)*(-vk*ctheta*cphi)   !(vorb(1)*(vk*stheta*cphi)
     &     +vorb(2)*(-vk*stheta*cphi)     !+vorb(2)*(-vk*ctheta*cphi)
     &     +vorb(3)*(-vk*sphi))/(vabs*vk) !+vorb(3)*(vk*sphi))/(vabs*vk)
!      if(nov1(1).eq.-1..and.nov1(2).eq.-1..and.nov1(3).eq.-1.)then
!         write(*,*)'kick1',vabs,vk,angl
!      else
!         write(*,*)'kick2',vabs,vk,angl
!      endif
*
      vorb(1) = vorb(1)+vk*ctheta*cphi !vorb(1)+vk*stheta*cphi
      vorb(2) = vorb(2)+vk*stheta*cphi !vorb(2)-vk*ctheta*cphi
      vorb(3) = vorb(3)+vk       *sphi !vorb(3)+vk       *sphi
* A. Tanikawa corrects reinitialization 21/08/27
      vloc1(1,1) = vloc0(1,1)+vk*ctheta*cphi
      vloc1(2,1) = vloc0(2,1)+vk*stheta*cphi
      vloc1(3,1) = vloc0(3,1)+vk       *sphi
      vloc1(1,2) = vloc0(1,2)
      vloc1(2,2) = vloc0(2,2)
      vloc1(3,2) = vloc0(3,2)
      minv       = 1/(m1n+m2)
      vtmp(1) = (m1n*vloc1(1,1)+m2*vloc1(1,2))*minv
      vtmp(2) = (m1n*vloc1(2,1)+m2*vloc1(2,2))*minv
      vtmp(3) = (m1n*vloc1(3,1)+m2*vloc1(3,2))*minv
      delv(1) = vtmp(1)*crotlrl-vtmp(2)*srotlrl
      delv(2) = vtmp(1)*srotlrl+vtmp(2)*crotlrl
      delv(3) = vtmp(3)
      mred       = m1n*m2*minv
      vunit      = yearsc/rsunkm
      jmom1(1)   = mred*(xorb(2)*vorb(3)-xorb(3)*vorb(2))*vunit
      jmom1(2)   = mred*(xorb(3)*vorb(1)-xorb(1)*vorb(3))*vunit
      jmom1(3)   = mred*(xorb(1)*vorb(2)-xorb(2)*vorb(1))*vunit
      delj(1)    = jmom1(1)*crotlrl-jmom1(2)*srotlrl !jmom1(1)!-jmom0(1)
      delj(2)    = jmom1(1)*srotlrl+jmom1(2)*crotlrl !jmom1(2)!-jmom0(2)
      delj(3)    = jmom1(3) !jmom1(3)!-jmom0(3)
      gminv      = 1./(GravConst*(m1n+m2))*vunit/mred
      rinv       = 1./sqrt(xorb(1)**2+xorb(2)**2+xorb(3)**2)
      lrl1(1) = gminv*(vorb(2)*jmom1(3)-vorb(3)*jmom1(2))-xorb(1)*rinv
      lrl1(2) = gminv*(vorb(3)*jmom1(1)-vorb(1)*jmom1(3))-xorb(2)*rinv
      lrl1(3) = gminv*(vorb(1)*jmom1(2)-vorb(2)*jmom1(1))-xorb(3)*rinv
      lrlabs  = sqrt(lrl1(1)**2+lrl1(2)**2+lrl1(3)**2)
      lrlinv  = 1./lrlabs
      dele(1) = lrl1(1)*crotlrl-lrl1(2)*srotlrl
      dele(2) = lrl1(1)*srotlrl+lrl1(2)*crotlrl
      dele(3) = lrl1(3)
*
*
*
* Determine the magnitude of the new relative velocity.
!      vn2 = vk2+vr2-2.d0*vk*vr*(ctheta*cphi*salpha-stheta*cphi*calpha)
      vn2 = vk2+vr2-2.d0*vk*vr*(ctheta*cphi*salpha+stheta*cphi*calpha)
* Calculate the new semi-major axis.
      sep = 2.d0/r - vn2/(gmrkm*(m1n+m2))
      sep = 1.d0/sep
*     if(sep.le.0.d0)then
*        ecc = 1.1d0
*        goto 90
*     endif
* Determine the magnitude of the cross product of the separation vector
* and the new relative velocity.
      v1 = vk2*sphi*sphi
      v2 = (vk*ctheta*cphi-vr*salpha)**2
      hn2 = r*r*(v1 + v2)
* Calculate the new eccentricity.
      ecc2 = 1.d0 - hn2/(gmrkm*sep*(m1n+m2))
      ecc2 = MAX(ecc2,0.d0)
      ecc = SQRT(ecc2)
* Calculate the new orbital angular momentum taking care to convert
* hn to units of Rsun^2/yr.
      jorb = (m1n*m2/(m1n+m2))*SQRT(hn2)*(yearsc/rsunkm)
* A. Tanikawa corrects reinitialization 21/08/27
      if(abs(lrlabs-ecc).gt.tolecc .and. sep.gt.0.) then
         write(*,*)'Error: LR vec after SN: ',abs(lrlabs-ecc)
         write(*,*)1./lrlinv,ecc
         stop
      else
         dele(1) = ecc*lrlinv*dele(1)
         dele(2) = ecc*lrlinv*dele(2)
         dele(3) = ecc*lrlinv*dele(3)
      endif
*
* Determine the angle between the new and old orbital angular
* momentum vectors.
      cmu = (vr*salpha-vk*ctheta*cphi)/SQRT(v1 + v2)
      mu = ACOS(cmu)
* Calculate the components of the velocity of the new centre-of-mass.
 90   continue
      if(ecc.le.1.0)then
* Calculate the components of the velocity of the new centre-of-mass.
         mx1 = vk*m1n/(m1n+m2)
         mx2 = vr*(m1-m1n)*m2/((m1n+m2)*(m1+m2))
         vs(1) = mx1*ctheta*cphi + mx2*salpha
         vs(2) = mx1*stheta*cphi + mx2*calpha
         vs(3) = mx1*sphi
* A. Tanikawa corrects reinitialization 21/08/27
         delvabs = sqrt(delv(1)**2+delv(2)**2+delv(3)**2)
         vsabs   = sqrt(vs(1)**2+vs(2)**2+vs(3)**2)
         errvel = (delvabs-vsabs)/vsabs
         if(errvel.gt.tolvel .and. vsabs.ne.0.) then
            write(*,*)'Error: v vec after SN: '
            stop
         endif
*
*     A. Tanikawa adds here 20/05/05
         call inclinationByKick(cinc,nov1,xorb,vorb)
*
      else
* Calculate the relative hyperbolic velocity at infinity (simple method).
         sep = r/(ecc-1.d0)
*        cmu = SQRT(ecc-1.d0)
*        mu = ATAN(cmu)
         mu = ACOS(1.d0/ecc)
         vr2 = gmrkm*(m1n+m2)/sep
         vr = SQRT(vr2)
         vs(1) = vr*SIN(mu)
         vs(2) = vr*COS(mu)
         vs(3) = 0.d0
         ecc = MIN(ecc,99.99d0)
      endif
*
 95   continue
*
      RETURN
      END
***

      subroutine inclinationByKick(cinc,nov1,xorb,vorb)
      implicit none
*
      real*8 pi
      integer idum
      common /value3/ idum
      real*8,dimension(2) :: cinc
      real*8,dimension(3) :: nov0,nov1,nov2
      real*8,dimension(3) :: xorb,vorb,jorb
      integer k,k1,k2
      real*8 jorb2,jorbinv
      real*8 ctheta,stheta,cphi,sphi,phi
      real ran3
      external ran3

      jorb2 = 0.
      do k = 1,3
         k1 = mod(k+1,3)
         if(k1.eq.0)k1=3
         k2 = mod(k+2,3)
         if(k2.eq.0)k2=3
         jorb(k) = xorb(k1)*vorb(k2) - xorb(k2)*vorb(k1)
         jorb2   = jorb2 + jorb(k)*jorb(k)
      enddo

      jorbinv = 1./sqrt(jorb2)
      do k = 1,3
         jorb(k) = jorb(k)*jorbinv
      enddo

      if(nov1(1).eq.-1..and.nov1(2).eq.-1..and.nov1(3).eq.-1.)then
         do k = 1,3
            nov1(k) = jorb(k)
         enddo
         cinc(1) = nov1(3)
      else
         pi      = ACOS(-1.d0)
         phi     = 2.*pi*ran3(idum)
         cphi    = cos(phi)
         sphi    = sin(phi)
         ctheta  = nov1(3)
         stheta  = sqrt(1.-ctheta**2)
         nov0(1) = stheta*cphi
         nov0(2) = stheta*sphi
         nov0(3) = ctheta
         nov1(1) = 0.
         nov1(2) = 0.
         nov1(3) = 1.
         nov2(1) = jorb(1)
         nov2(2) = jorb(2)
         nov2(3) = jorb(3)
         cinc(1) = nov0(1)*nov2(1)+nov0(2)*nov2(2)+nov0(3)*nov2(3)
         cinc(2) = nov1(1)*nov2(1)+nov1(2)*nov2(2)+nov1(3)*nov2(3)
      endif

      return
      end

