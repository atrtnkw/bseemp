***
      PROGRAM bse
***
*
* Evolves a binary by calling evolv2.f 
* (see header of subroutine for algorithm description). 
*
* Required input is described below. 
***
* See Tout et al., MNRAS, 1997, 291, 732 for a description of many of the
* processes in this code as well as the relevant references mentioned
* within the code.
* Updated reference is:
*           Hurley J.R., Tout C.A., & Pols O.R., 2002, MNRAS, 329, 897
* (please use this one).
***
* For single star evolution see Hurley, Pols & Tout, 2000, MNRAS, 315, 543.
* or Hurley, 2000, PhD Thesis, University of Cambridge (Chapter 2).
* The binary evolution algorithm is described in Chapter 3 of the thesis.
***
*
*           B I N A R Y
*           ***********
*
*       Roche lobe overflow.
*       --------------------
*
*       Developed by Jarrod Hurley, IOA, Cambridge.
*       .........................................................
*
*       Advice by Christopher Tout, Onno Pols & Sverre Aarseth.
*       ++++++++++++++++++++++++++++++++++++++++++++++++++
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
      integer kw,kw2,kstar(2),j,k,time
*
      real*8 mass0(2),mass(2),z,zpars(20)
      real*8 epoch(2),tms(2),tphys,tphysf,dtp,aj
      real*8 rad(2),lum(2),ospin(2)
      real*8 massc(2),radc(2),menv(2),renv(2)
      real*8 tb,ecc,yearsc
      PARAMETER(yearsc=3.1557d+07)
      CHARACTER*8 label(14)
*
* Tanikawa adds variables and header
*
      real*8 zeta
      include 'cppinterface.h'
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
* Input:
*
* mass is in solar units.
* tphysf is the maximum evolution time in Myr.
* tb is the orbital period in days.
* kstar is the stellar type: 0 or 1 on the ZAMS - unless in evolved state. 
* z is metallicity in the range 0.0001 -> 0.03 where 0.02 is Population I.
* eccentricity can be anywhere in the range 0.0 -> 1.0.
*
* neta is the Reimers mass-loss coefficent (neta*4x10^-13: 0.5 normally). 
* bwind is the binary enhanced mass loss parameter (inactive for single).
* hewind is a helium star mass loss factor (1.0 normally).
* alpha1 is the common-envelope efficiency parameter (1.0).  
* lambda is the binding energy factor for common envelope evolution (0.5).
*
* ceflag > 0 activates spin-energy correction in common-envelope (0). #defunct#
* ceflag = 3 activates de Kool common-envelope model (0). 
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
* If you enter a negative kstar then parameters for an evolved star are
* required in the order of:
* current age, initial mass and spin rate, 
* otherwise the star will start on the ZAMS.
*
      OPEN(22,file='binary.in', status='old')
      READ(22,*)mass0(1),mass0(2),tphysf,tb,kstar(1),kstar(2),z,ecc
      READ(22,*)neta,bwind,hewind,alpha1,lambda,betaacc
      READ(22,*)ceflag,tflag,ifflag,wdflag,bhflag,
     &     nsflag,psflag,mxns,idum
      read(22,*)NewStarModel,WindEnhanced,
     &     RadiusShrinkage,NewDynTide,NewMassTransfer
      READ(22,*)pts1,pts2,pts3
      READ(22,*)sigma,beta,xi,acc2,epsnov,eddfac,gamma
      if(kstar(1).lt.0.or.kstar(2).lt.0)then
         READ(22,*)tphys
         READ(22,*)aj,mass(1),ospin(1)
         epoch(1) = tphys - aj
         kstar(1) = ABS(kstar(1))
         READ(22,*)aj,mass(2),ospin(2)
         epoch(2) = tphys - aj
         kstar(2) = ABS(kstar(2))
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
      else
*
* Initialize the parameters.
* Set the initial spin of the stars. If ospin is zero (actually < 0.001)
* at time zero then evolv2 will set an appropriate ZAMS spin. If 
* ospin is greater than zero then it will start with that spin regardless
* of the time. If you want to start at time zero with negligible spin 
* then I suggest using a negligible value (but greater than 0.001).
* If ospin is negative then the stars will be in co-rotation with the orbit.
*
         tphys = 0.d0
         mass(1) = mass0(1)
         epoch(1) = 0.d0
         ospin(1) = 0.d0
         mass(2) = mass0(2)
         epoch(2) = 0.d0
         ospin(2) = 0.d0
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
      endif
      if(idum.gt.0) idum = -idum
      CLOSE(22)
      write(0,*)'NewStarModel: ',    NewStarModel
      write(0,*)'WindEnhanced: ',    WindEnhanced
      write(0,*)'RadiusShrinkage: ', RadiusShrinkage
      write(0,*)'NewDynTide: ',      NewDynTide
      write(0,*)'NewMassTransfer: ', NewMassTransfer
*
* Note that this routine can be used to evolve a single star if you 
* simply set mass(2) = 0.0 or tb = 0.0 (setting both is advised as  
* well as some dummy value for ecc). 
*
************************************************************************
*
* Set parameters which depend on the metallicity 
*
      CALL zcnsts(z,zpars)
*
* Set the collision matrix.
*
      CALL instar
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
* Set the data-save parameter. If dtp is zero then the parameters of the 
* star will be stored in the bcm array at each timestep otherwise they 
* will be stored at intervals of dtp. Setting dtp equal to tphysf will 
* store data only at the start and end while a value of dtp greater than 
* tphysf will mean that no data is stored.
*
      dtp = 0.d0
*
* Evolve the binary.
* 
      CALL evolv2(kstar,mass0,mass,rad,lum,massc,radc,
     &            menv,renv,ospin,epoch,tms,
     &            tphys,tphysf,dtp,z,zpars,tb,ecc,cdeli,
* A. Tanikawa adds here for merger inflation 21/10/14
     &            tmrg,tcon,tkhq,rmrg,emrg,
*
* A. Tanikawa corrects reinitialization 21/08/27
     &     reinit,ip,jp,ffb,typesn,arem,cinc,nov1,
     &     catcherr,catchsn,catchce,catchmg,catchnr,catchbd,
     &     msec,dmsec,jsec,djsec,esec,desec,delv,delj,dele)
*
*
************************************************************************
* Output:
* First check that bcm is not empty.
*
      if(bcm(1,1).lt.0.0) goto 50
*
* The bcm array stores the stellar and orbital parameters at the 
* specified output times. The parameters are (in order of storage):
*
*    Time, 
*    [stellar type, initial mass, current mass, log10(L), log10(r),
*    log10(Teff), core mass, core radius, mass of any convective 
*    envelope, radius of the envelope, epoch, spin, mass loss rate and 
*    ratio of radius to roche lobe radius (repeated for secondary)],
*    period, separation, eccentricity.
*
      OPEN(23,file='binary.dat', status='unknown')
      j = 0
 30   j = j + 1
      if(bcm(j,1).lt.0.0)then
         bcm(j-1,1) = bcm(j,1)
         j = j - 1
      endif
      kw = INT(bcm(j,2))
      kw2 = INT(bcm(j,16))
      WRITE(23,99)bcm(j,1),
     &     kw,kw2,
     &     bcm(j, 4),bcm(j,18), ! mt
     &     bcm(j, 8),bcm(j,22), ! mc
     &     bcm(j, 6),bcm(j,20), ! rad
     &     bcm(j,15),bcm(j,29), ! rad/rol
     &     bcm(j, 5),bcm(j,19), ! lumin
     &     bcm(j, 7),bcm(j,21), ! teff
     &     bcm(j,13),bcm(j,27), ! dimensionless spin ! ospin
     &     bcm(j,14),bcm(j,28), ! dmt-dmr
     &     bcm(j,33),bcm(j,34), ! k2str
     &     bcm(j,31),bcm(j,32)  ! sep,ecc
      if(bcm(j,1).ge.0.0) goto 30
      CLOSE(23)
 99   FORMAT(f10.4,2i3,12f10.4,7e12.4,f7.3)
 999  FORMAT(f10.4,2f10.4,1p,2e12.4)
*
* The bpp array acts as a log, storing parameters at each change
* of evolution stage.
*
 50   j = 0
      WRITE(*,*)'     TIME      M1       M2   K1 K2        SEP    ECC',  
     &          '  R1/ROL1 R2/ROL2  TYPE'
 52   j = j + 1
      if(bpp(j,1).lt.0.0) goto 60
      kstar(1) = INT(bpp(j,4))
      kstar(2) = INT(bpp(j,5))
      kw = INT(bpp(j,10))
* A. Tanikawa modifies here for spin 20/03/03
!      WRITE(*,100)(bpp(j,k),k=1,3),kstar,(bpp(j,k),k=6,9),label(kw)
      WRITE(*,100)(bpp(j,k),k=1,3),kstar,(bpp(j,k),k=6,9),label(kw),
     &     bpp(j,11),bpp(j,12),bpp(j,13),bpp(j,14),
     &     snlabel(int(bpp(j,15))),snlabel(int(bpp(j,16))),
     &     bpp(j,17),bpp(j,18)
*
      goto 52
 60   continue
* A. Tanikawa modifies here for spin 20/03/03
! 100  FORMAT(f11.4,2f9.3,2i3,f13.3,f6.2,2f8.3,2x,a8)
! 100  FORMAT(f11.4,2f9.3,2i3,f13.3,f6.2,2f8.3,2x,a8,f9.3,f9.3)
! 100  FORMAT(f11.4,2f9.3,2i3,f13.3,f12.6,2f8.3,2x,a8,f9.3,f9.3,
! 100  FORMAT(f11.4,2f9.3,2i3,E12.4,f12.6,2f8.3,2x,a8,2E12.4,
! 100  FORMAT(f14.4,2f9.3,2i3,f13.3,f12.6,2f8.3,2x,a8,f9.3,f9.3,
!     &     2f8.3,2a5)
 100  FORMAT(f11.4,2f10.3,2i3,E13.4,f13.6,2e12.3,2x,a8,2E12.4,
     &     f9.3,f9.3,a5,a5,2f10.3)
*
      WRITE(*,*)
*
************************************************************************
*
      STOP
      END
***
