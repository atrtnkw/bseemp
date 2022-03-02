***
      PROGRAM popbin
***
*
* Evolves a population of binaries using input parameters 
* read from input file binaries.in (M1, M2, P, e, Z, Tmax). 
*
***
*
* Tanikawa adds this announcement.
*
      use iso_c_binding
*
      implicit none
*
      INCLUDE 'const_bse.h'
*
      integer i,j,k,jj,nm1
      integer kw,kw2,kwx,kwx2,kstar(2)
      integer i1,i2,kdum
*
      real*8 m1,m2,tmax
      real*8 mass0(2),mass(2),z,zpars(20)
      real*8 epoch(2),tms(2),tphys,tphysf,dtp
      real*8 rad(2),lum(2),ospin(2)
      real*8 massc(2),radc(2),menv(2),renv(2)
      real*8 sep0,tb0,tb,ecc0,ecc,aursun,yeardy,yearsc,tol
      PARAMETER(aursun=214.95d0,yeardy=365.25d0,yearsc=3.1557d+07)
      PARAMETER(tol=1.d-07)
      real*8 t1,t2,mx,mx2,tbx,eccx
*
* Tanikawa adds variables and header
*
      real*8 zeta
      include 'cppinterface.h'
      character(64) :: filename2
      CHARACTER*8 label(14)
*
*
************************************************************************
* BSE parameters:
*
* neta is the Reimers mass-loss coefficent (neta*4x10^-13: 0.5 normally). 
* bwind is the binary enhanced mass loss parameter (inactive for single).
* hewind is a helium star mass loss factor (1.0 normally).
* alpha1 is the common-envelope efficiency parameter (1.0).  
* lambda is the binding energy factor for common envelope evolution (0.5).
*
* ceflag > 0 activates spin-energy correction in common-envelope (0). 
* tflag > 0 activates tidal circularisation (1).
* ifflag > 0 uses WD IFMR of HPE, 1995, MNRAS, 272, 800 (0). 
* wdflag > 0 uses modified-Mestel cooling for WDs (0). 
* bhflag > 0 allows velocity kick at BH formation (0). 
* nsflag > 0 takes NS/BH mass from Belczynski et al. 2002, ApJ, 572, 407 (1). 
* mxns is the maximum NS mass (1.8, nsflag=0; 3.0, nsflag=1). 
* idum is the random number seed used by the kick routine. 
*
* Next come the parameters that determine the timesteps chosen in each
* evolution phase:
*                 pts1 - MS                  (0.05) 
*                 pts2 - GB, CHeB, AGB, HeGB (0.01)
*                 pts3 - HG, HeMS            (0.02)
* as decimal fractions of the time taken in that phase.
*
* sigma is the dispersion in the Maxwellian for the SN kick speed (190 km/s). 
* beta is wind velocity factor: proportional to vwind**2 (1/8). 
* xi is the wind accretion efficiency factor (1.0). 
* acc2 is the Bondi-Hoyle wind accretion factor (3/2). 
* epsnov is the fraction of accreted matter retained in nova eruption (0.001). 
* eddfac is Eddington limit factor for mass transfer (1.0).
* gamma is the angular momentum factor for mass lost during Roche (-1.0). 
*
      neta = 0.5 !0.0
      bwind = 0.0
      hewind = 1.0
      alpha1 = 1.0 !3.0
      lambda = 1.0 !0.5 
      ceflag = 0
      tflag = 1
      ifflag = 0 
      wdflag = 1 
      bhflag = 0
      nsflag = 1
      mxns = 3.0
      pts1 = 0.05
      pts2 = 0.01
      pts3 = 0.02
      sigma = 0.0 !190.0
      beta = 0.125
      xi = 1.0 
      acc2 = 1.5
      epsnov = 0.001
      eddfac = 10.0
      gamma = -1.0
*
* Set the seed for the random number generator. 
*
      idum = 3234
      if(idum.gt.0) idum = -idum
*
* Set the collision matrix.
*
      CALL instar
*
* Tanikawa adds here
*
      label(1) = 'INITIAL '
      label(2) = 'KW CHNGE'
      label(3) = 'BEG RCHE'
      label(4) = 'END RCHE'
      label(5) = 'CONTACT '
      label(6) = 'COELESCE'
      label(7) = 'COMENV  '
      label(8) = 'GNTAGE  '
      label(9) = 'NO REMNT'
      label(10) = 'MAX TIME'
      label(11) = 'DISRUPT '
      label(12) = 'BEG SYMB'
      label(13) = 'END SYMB'
      label(14) = 'BEG BSS'
*
*
* Open the input file - list of binary initial parameters. 
*
      OPEN(10,file='binaries.in',status='unknown')
      READ(10,*)nm1
*
* Open the output files. 
*
      OPEN(11,file='binaries.out',status='unknown')
      OPEN(12,file='search.out',status='unknown')
*
      do i = 1,nm1
*
* Read in parameters and set coefficients which depend on metallicity. 
*
         READ(10,*)m1,m2,tb,ecc,z,tmax
*
* Tanikawa modifies metallicity
*
         if(z.lt.0.0001) then
            zeta = nint(log10(z/0.02))
            call setMetallicity(zeta)
            z = 0.0001
         endif
         if(.true. .and. (0.00019 .lt. z .and. z .lt. 0.00021)) then
!         if(.false. .and. (0.00019 .lt. z .and. z .lt. 0.00021)) then
            zeta = nint(log10(z/0.02))
            call setMetallicity(zeta)
            write(0,*) 'Caution: This is Tanikawa FFs (log(Z/Zsun)=-2).'
         endif
*
         CALL zcnsts(z,zpars)
*
         ecc0 = ecc
         tb0 = tb/yeardy
         sep0 = aursun*(tb0*tb0*(mass(1) + mass(2)))**(1.d0/3.d0)
         tb0 = tb
*
* Initialize the binary. 
*
         kstar(1) = 1
         mass0(1) = m1
         mass(1) = m1
         massc(1) = 0.0
         ospin(1) = 0.0
         epoch(1) = 0.0
*
         kstar(2) = 1
         mass0(2) = m2
         mass(2) = m2
         massc(2) = 0.0
         ospin(2) = 0.0
         epoch(2) = 0.0
*
         tphys = 0.0
         tphysf = tmax
         dtp = 0.0
*
* Evolve the binary. 
*
         CALL evolv2(kstar,mass0,mass,rad,lum,massc,radc,
     &               menv,renv,ospin,epoch,tms,
     &               tphys,tphysf,dtp,z,zpars,tb,ecc)
*
* Search the BCM array for the formation of binaries of 
* interest (data on unit 12 if detected) and also output 
* the final state of the binary (unit 11). 
*
* In this example we will search for CVs. 
*
         jj = 0
         t1 = -1.0
         t2 = -1.0
 30      jj = jj + 1
         if(bcm(jj,1).lt.0.0) goto 40
         kw = INT(bcm(jj,2))
         kw2 = INT(bcm(jj,16))
*
         i1 = 15
         i2 = 29
         if(kw.gt.kw2)then
            kdum = kw2
            kw2 = kw
            kw = kdum
            i2 = 15
            i1 = 29
         endif 
*
         if(kw.le.1.and.bcm(jj,i1).ge.1.0)then
            if(kw2.ge.10.and.kw2.le.12)then
               if(t1.lt.0.0)then
                  t1 = bcm(jj,1)
                  kwx = kw
                  kwx2 = kw2
                  mx = bcm(jj,i1-11)
                  mx2 = bcm(jj,i2-11)
                  tbx = bcm(jj,30)
                  eccx = bcm(jj,32)
               endif
            endif
         endif
*
         if(t1.gt.0.0.and.(bcm(jj,i1).lt.1.0.or.
     &      kw.ne.kwx.or.kw2.ne.kwx2))then
            if(t2.lt.0.0)then
               t2 = bcm(jj,1)
               if(t2.gt.(t1+tol))then
                  WRITE(12,112)m1,m2,ecc0,tb0,t1,t2,kwx,kwx2,
     &                         mx,mx2,eccx,tbx
               endif
               t1 = -1.0
               t2 = -1.0
            endif
         endif
*
         goto 30
 40      continue
*
         if(t1.gt.0.0)then
            if(t2.lt.0.0) t2 = tmax
            WRITE(12,112)m1,m2,ecc0,tb0,t1,t2,kwx,kwx2,mx,mx2,eccx,tbx
         endif
*
         jj = jj - 1
         kw = INT(bcm(jj,2))
         kw2 = INT(bcm(jj,16))
         mx = bcm(jj,4)
         mx2 = bcm(jj,18)
         tbx = bcm(jj,30)*yeardy
         eccx = bcm(jj,32)
         WRITE(11,111)tmax,kw,kw2,mx,mx2,eccx,tbx
*
*
* Tanikawa adds here
*
!         if(askInUseOrNot()) then
            write(filename2,"('./output/binary',1i7.7'.txt')") i
            open(4,file='errorwrite.txt',status='unknown')
            open(24,file=filename2, status='unknown')
            j = 0
 50         j = j + 1
            if(bcm(j,1).lt.0.0)then
               bcm(j-1,1) = bcm(j,1)
               j = j - 1
            endif
            kw = INT(bcm(j,2))
            kw2 = INT(bcm(j,16))
            if(i.lt.10000)then
            endif
            if(bcm(j,1).ge.0.0) goto 50
            CLOSE(23)
 99         FORMAT(f10.4,2i3,10f10.4,5e12.4,f7.3)
 999        FORMAT(f10.4,2f10.4,1p,2e12.4)
            j = 0
 52         j = j + 1
            if(bpp(j,1).lt.0.0) goto 60
            kstar(1) = INT(bpp(j,4))
            kstar(2) = INT(bpp(j,5))
            kw = INT(bpp(j,10))
            WRITE(24,100)(bpp(j,k),k=1,3),kstar,(bpp(j,k),k=6,9),
     &           label(kw)
            goto 52
 60         continue
            if(bpp(j-1,1).lt.15000)then
               write(4,*)'error',i,mass0(1),mass0(2)
            endif
            
 100        FORMAT(f11.4,2f9.3,2i3,f13.3,f6.2,2f8.3,2x,a8)
!      endif
*     
      enddo
*
 111  FORMAT(f10.1,2i3,3f8.3,1p,e14.6)
 112  FORMAT(3f8.3,1p,e14.6,0p,2f10.2,2i3,3f8.3,1p,e14.6)
      CLOSE(10)
      CLOSE(11)
      CLOSE(12)
*
************************************************************************
*
      STOP
      END
***
