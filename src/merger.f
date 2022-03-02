      subroutine bhmerger(mass,arem,cinc,
     &     kstar1,kstar2,mass1,mass2,
     &     arem1,arem2,cinc1,cinc2)
      implicit none
*
      real*8 p0,p1
      parameter(p0=0.04827d0,p1=0.01707d0)
      real*8 t0,t2,t3,s4,s5
      parameter(t0=-2.8904d0,t2=-3.51712d0,t3=2.5763d0,
     &     s4=-0.1229d0,s5=0.4537d0)
      real*8 xideg
      parameter(xideg=145.d0)
      real*8 am,bm
      parameter(am=1.2d4,bm=-0.93)
      real*8 hs
      parameter(hs=6.9d3)
      real*8 v11,va,vb,vc
      parameter(v11=3.678d3,va=2.481d3,vb=1.792d3,vc=1.507d3)
      real*8 c2,c3
      parameter(c2=1.507d3,c3=2.481d3)
*
      real*8 mass(2),arem(2),cinc(2)
*
      integer idum
      common /value3/ idum
      integer idum2,iy,ir(32)
      common /rand3/ idum2,iy,ir
      real*8 ran3,xx
*
      integer k,i1,i2
      real*8 pi,phi,cphi,sphi,sinc
      real*8 q,eta
      real*8 chi(3,2),delta(2),xieff(2),chiabs
      real*8 delta2(3),xieff2(3)
      real*8 eisco,risco,z1,z2
      real*8 massf,chif(3)
      real*8 llarge,chifabs
      real*8 vm,vs1,vs2
* Output
      integer kstar1,kstar2
      real*8 mass1,mass2
      real*8 arem1,arem2,cinc1,cinc2
      real*8 vrec(3)
*
      pi = ACOS(-1.d0)
*
      do 10 , k = 1,2
         phi = 2.d0*pi*ran3(idum)
         cphi = cos(phi)
         sphi = sin(phi)
         sinc = sqrt(1.d0-cinc(k)**2)
         chiabs   = min(1.,arem(k))
         chi(1,k) = chiabs*sinc*cphi
         chi(2,k) = chiabs*sinc*sphi
         chi(3,k) = chiabs*cinc(k)
 10   continue
*
      if(mass(1).ge.mass(2)) then
         i1 = 1
         i2 = 2
      else
         i1 = 2
         i2 = 1
      endif   
      q   = mass(i2)/mass(i1)
      eta = q/(1.d0+q)**2
*
      do 20 , k = 1,3
         delta2(k) = (chi(k,i1)-q   *chi(k,i2))/(1.d0+q)
         xieff2(k) = (chi(k,i1)+q**2*chi(k,i2))/(1.d0+q)**2
 20   continue
      delta(1) = sqrt(delta2(1)**2+delta2(2)**2)
      delta(2) = delta2(3)
      xieff(1) = sqrt(xieff2(1)**2+xieff2(2)**2)
      xieff(2) = xieff2(3)
* Calculate final mass
      z1 = 1.d0+(1.d0-xieff(2)**2)**(1/3.)
     &     *((1.d0+xieff(2))**(1/3.)+(1.d0-xieff(2))**(1/3.))
      z2 = sqrt(3.d0*xieff(2)**2+z1**2)
      risco = 3.d0+z2-sign(1.d0,xieff(2))
     &     *sqrt((3.d0-z1)*(3.d0+z1+2.d0*z2))
      eisco = sqrt(1.d0-2.d0/(3.d0*risco))
      massf = 1.d0-eta*(1.d0+4.d0*eta)*(1.d0-eisco)
     &     -16.d0*eta**2*(p0+4.d0*p1*xieff(2)*(xieff(2)+1.d0))
      massf = massf*(mass(i1)+mass(i2))
* Calculate final spin
      llarge = 2.d0*sqrt(3.d0)+t2*eta+t3*eta**2
     &     +s4*(1.d0+q)**4/(1.d0+q**2)**2*(xieff(1)**2+xieff(2)**2)
     &     +(s5*eta+t0+2.d0)*(1.d0+q)**2/(1.d0+q**2)*abs(xieff(2))
      chifabs = 0.d0
      do 30 , k = 1,3
         chif(k) = xieff2(k)
         if(k.eq.3) then
            chif(k) = chif(k)+q/(1.d0+q)**2*llarge
         endif
         chifabs = chifabs+chif(k)**2
 30   continue
      chifabs = sqrt(chifabs)
      if(chifabs.gt.1) then
         do 40 , k = 1,3
            chif(k) = chif(k)/chifabs
 40      continue
      endif
* Calculate recoil kick
      vm  = am*eta**2*(1.d0-q)/(1.d0+q)*(1.d0+bm*eta)
      vs1 = hs*eta**2*delta(2)
      vs2 = 16.d0*eta**2*cos(pi*ran3(idum))
     &     *(delta(1)*(v11+2.d0*va*xieff(2)
     &     +4.d0*vb*xieff(2)**2+8.d0*vc*xieff(2)**3)
     &     +2.d0*xieff(1)*delta(2)*(c2+2.d0*c3*xieff(2)))
* Interface
      kstar1 = 15
      kstar2 = 14
      mass1  = 0.
      mass2  = massf
      arem1  = 0.
      arem2  = chifabs
      cinc1  = 0.
      cinc2  = chif(3)/chifabs
      vrec(1) = vm+vs1*cos(xideg/180.*pi)
      vrec(2) =    vs1*sin(xideg/180.*pi)
      vrec(3) = vs2
*
      return
      end
