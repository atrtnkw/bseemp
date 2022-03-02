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
      integer i,j,k,jj,nm1,iqt,isp,bbunch
      parameter(bbunch=10000)
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
      real*8 tphysi,mass1i,mass2i,sepi,ecci
      integer kw1i,kw2i
      character*4 snlabel(6)
      real*8 tmrg(2),tcon(2),tkhq(2),rmrg(2),emrg(2)
*
* A. Tanikawa corrects reinitialization 21/08/27
      logical reinit
      parameter(reinit=.false.)
      integer ip,jp
      real*8,dimension(2) :: cinc
      real*8,dimension(3) :: nov1
      real*8,dimension(2) :: ffb
      real*8,dimension(2) :: typesn
      real*8,dimension(2) :: arem
      integer catcherr,catchsn,catchce,catchmg,catchnr,catchbd
      real*8 msec(2),dmsec(2),jsec,djsec,esec,desec,cdeli
      real*8,dimension(3) :: delv,delj,dele
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
!      neta = 0.0 !0.5
!      bwind = 0.0
!      hewind = 0.0 !1.0
!      alpha1 = 1.0
!      lambda = 1.0
!      ceflag = 0
!      tflag = 1
!      ifflag = 0 
!      wdflag = 1 
!      bhflag = 0 !1
!      nsflag = 1
!      mxns = 3.0
!      idum = 3234
!      pts1 = 0.05
!      pts2 = 0.01
!      pts3 = 0.02
!      sigma = 190.0 !265.0
!      beta = 0.125
!      xi = 1.0 
!      acc2 = 1.5
!      epsnov = 0.001
!      eddfac = 10.0
!      gamma = -1.0
      open(22,file='header.in',status='unknown')
* Tanikawa adds here 21/04/20
      read(22,*)z
*
      read(22,*)neta,bwind,hewind,alpha1,lambda,betaacc
      read(22,*)ceflag,tflag,ifflag,wdflag,bhflag,
     &     nsflag,psflag,mxns,idum
      read(22,*)NewStarModel,WindEnhanced,
     &     RadiusShrinkage,NewDynTide,NewMassTransfer
      read(22,*)pts1,pts2,pts3
      read(22,*)sigma,beta,xi,acc2,epsnov,eddfac,gamma
      close(22)
      write(0,*)'neta:',    neta
      write(0,*)'hewind:',  hewind
      write(0,*)'alpha1:',  alpha1
      write(0,*)'lambda:',  lambda
      write(0,*)'betaacc:', betaacc
      write(0,*)'bhflag: ', bhflag
      write(0,*)'nsflag: ', nsflag
      write(0,*)'psflag: ', psflag
      write(0,*)'NewStarModel: ',    NewStarModel
      write(0,*)'WindEnhanced: ',    WindEnhanced
      write(0,*)'RadiusShrinkage: ', RadiusShrinkage
      write(0,*)'NewDynTide: ',      NewDynTide
      write(0,*)'NewMassTransfer: ', NewMassTransfer
* Output run parameter
      open(4,file='./output/runparameter.dat',status='unknown')
      write(4,*) 'bbunch', bbunch
      close(4)
*
* Set the seed for the random number generator. 
*
      if(idum.gt.0) idum = -idum
*
* Set the collision matrix.
*
      CALL instar
*
* Tanikawa adds here
*
* Tanikawa change here 20/12/22
!      label(1) = 'INITIAL '
!      label(2) = 'KW CHNGE'
!      label(3) = 'BEG RCHE'
!      label(4) = 'END RCHE'
!      label(5) = 'CONTACT '
!      label(6) = 'COELESCE'
!      label(7) = 'COMENV  '
!      label(8) = 'GNTAGE  '
!      label(9) = 'NO REMNT'
!      label(10) = 'MAX TIME'
!      label(11) = 'DISRUPT '
!      label(12) = 'BEG SYMB'
!      label(13) = 'END SYMB'
!      label(14) = 'BEG BSS'
      label(1) = 'INITIAL '
      label(2) = 'KW_CHNGE'
      label(3) = 'BEG_RCHE'
      label(4) = 'END_RCHE'
      label(5) = 'CONTACT '
      label(6) = 'COELESCE'
      label(7) = 'COMENV  '
      label(8) = 'GNTAGE  '
      label(9) = 'NO_REMNT'
      label(10) = 'MAX_TIME'
      label(11) = 'DISRUPT '
      label(12) = 'BEG_SYMB'
      label(13) = 'END_SYMB'
      label(14) = 'BEG_BSS'
      snlabel(1) = 'NOSN'
      snlabel(2) = 'CCSN'
      snlabel(3) = 'DC  '
      snlabel(4) = 'PPI '
      snlabel(5) = 'PISN'
      snlabel(6) = 'AIC '
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
* Tanikawa change here 21/04/20
!         READ(10,*)m1,m2,tb,ecc,z,tmax,tphys
         READ(10,*)m1,m2,tb,ecc,tmax,tphys
*
         CALL zcnsts(z,zpars)
*
         ecc0 = ecc
         tb0 = tb/yeardy
!         sep0 = aursun*(tb0*tb0*(mass(1) + mass(2)))**(1.d0/3.d0)
         sep0 = aursun*(tb0*tb0*(m1 + m2))**(1.d0/3.d0)
         tb0 = tb
*
* Initialize the binary. 
*
         kstar(1) = 1
         mass0(1) = m1
         mass(1) = m1
         massc(1) = 0.0
         ospin(1) = 0.0
         epoch(1) = tphys ! 0.0
*
         kstar(2) = 1
         mass0(2) = m2
         mass(2) = m2
         massc(2) = 0.0
         ospin(2) = 0.0
         epoch(2) = tphys ! 0.0
*
!         tphys = 0.0
         tphysf = tmax
         dtp = 0.0
* A. Tanikawa adds here for merger inflation 21/10/14
         tmrg(1) = -1.
         tmrg(2) = -1.
         tcon(1) = -1.
         tcon(2) = -1.
         rmrg(1) = 0.
         rmrg(2) = 0.
         emrg(1) = 0.
         emrg(2) = 0.
*
*
* Evolve the binary. 
*
         tphysi = tphys
         mass1i = mass(1)
         mass2i = mass(2)
         kw1i   = kstar(1)
         kw2i   = kstar(2)
         sepi   = sep0 !(6.673d-8*((mass1i+mass2i)*1.989e33))**(1/3.)
         ecci   = ecc0
         !write(*,*) i
         flush(0)
         CALL evolv2(kstar,mass0,mass,rad,lum,massc,radc,
     &               menv,renv,ospin,epoch,tms,
     &               tphys,tphysf,dtp,z,zpars,tb,ecc,cdeli,
* A. Tanikawa adds here for merger inflation 21/10/14
     &               tmrg,tcon,tkhq,rmrg,emrg,
*
* A. Tanikawa corrects reinitialization 21/08/27
     &     reinit,ip,jp,ffb,typesn,arem,cinc,nov1,
     &     catcherr,catchsn,catchce,catchmg,catchnr,catchbd,
     &     msec,dmsec,jsec,djsec,esec,desec,delv,delj,dele)
*
         !write(*,*) "#"
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
         WRITE(11,*)tmax,kw,kw2,mx,mx2,eccx,tbx,m1,m2,ecc0,tb0
         isp = mod(i,bbunch)
         iqt = i/bbunch
         if(isp.eq.1) then
            if(iqt.ne.0) then
               close(24)
               close(4)
            endif
            write(filename2,"('./output/binary',1i7.7'.txt')") iqt+1
            open(4,file='errorwrite.txt',status='unknown',
     &           position='append')
            open(24,file=filename2, status='unknown')
         endif
         write(24,*) '###', i
         write(24,100)tphysi,mass1i,mass2i,kw1i,kw2i,
     &        sepi,ecci,0.,0.,label(1),0.,0.,1.0,1.0,
     &        snlabel(1),snlabel(1)
         j = 0
 50      j = j + 1
         if(bcm(j,1).lt.0.0)then
            bcm(j-1,1) = bcm(j,1)
            j = j - 1
         endif
         kw = INT(bcm(j,2))
         kw2 = INT(bcm(j,16))
         if(i.lt.10000)then
         endif
         if(bcm(j,1).ge.0.0) goto 50
 99      FORMAT(f10.4,2i3,10f10.4,5e12.4,f7.3)
 999     FORMAT(f10.4,2f10.4,1p,2e12.4)
         j = 0
 52      j = j + 1
         if(bpp(j,1).lt.0.0) goto 60
         kstar(1) = INT(bpp(j,4))
         kstar(2) = INT(bpp(j,5))
         kw = INT(bpp(j,10))
         WRITE(24,100)(bpp(j,k),k=1,3),kstar,(bpp(j,k),k=6,9),
     &        label(kw),bpp(j,11),bpp(j,12),bpp(j,13),bpp(j,14),
     &        snlabel(int(bpp(j,15))),snlabel(int(bpp(j,16))),
     &        bpp(j,17),bpp(j,18)
         goto 52
 60      continue
         write(24,*)
         flush(24)
         if(bpp(j-1,1).lt.15000 .and. int(bpp(j-1,10)).ne.9)then
            write(4,*)'error',i,mass0(1),mass0(2)
         endif
!         close(24)
!         close(4)
 100     FORMAT(f11.4,2f10.3,2i3,E13.4,f13.6,2f8.3,2x,a8,2E12.4,
     &        f10.3,f10.3,a5,a5,2f10.3)
      enddo
 111  FORMAT(f10.1,2i3,3f8.3,1p,e14.6)
 112  FORMAT(3f8.3,1p,e14.6,0p,2f10.2,2i3,3f8.3,1p,e14.6)
      CLOSE(10)
      CLOSE(11)
      CLOSE(12)
      close(24)
      close(4)
*
************************************************************************
*
      STOP
      END
***
