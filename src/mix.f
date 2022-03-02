***
      SUBROUTINE MIX(M0,M,AJ,KS,ZPARS,krol,
* A. Tanikawa adds here 21/10/14
     &     tphys,tmrg,tcon,tkhq,rmrg,emrg)
*
*
*     Author : J. R. Hurley
*     Date :   7th July 1998
*
*       Evolution parameters for mixed star.
*       ------------------------------------
*
* Tanikawa adds this announcement.
*
      use iso_c_binding
*
      implicit none
*
      INTEGER KS(2),I1,I2,K1,K2,KW,ICASE,krol(2)
      INTEGER KTYPE(0:14,0:14)
      COMMON /TYPES/ KTYPE
      REAL*8 M0(2),M(2),AJ(2),ZPARS(20)
      REAL*8 TSCLS(20),LUMS(10),GB(10),TMS1,TMS2,TMS3,TN
      REAL*8 M01,M02,M03,M1,M2,M3,AGE1,AGE2,AGE3,MC3,MCH
      PARAMETER(MCH=1.44D0)
      REAL*8 NETA,BWIND,HEWIND,MXNS,betaacc
      COMMON /VALUE1/ NETA,BWIND,HEWIND,MXNS,betaacc
*
* A. Tanikawa adds here 21/10/14
*
      include 'cppinterface.h'
      real*8 tphys, q
      real*8 tmrg(2), tcon(2), tkhq(2), rmrg(2), emrg(2)
*
*
*
* A. Tanikawa adds here to prevent HGCE 21/02/22
      icase = ktype(krol(1),krol(2))
*
*       Define global indices with body #I1 being most evolved.
      IF(KS(1).GE.KS(2))THEN
          I1 = 1
          I2 = 2
      ELSE
          I1 = 2
          I2 = 1
      END IF
*
*       Specify case index for collision treatment.
      K1 = KS(I1)
      K2 = KS(I2)
* A. Tanikawa comments out here to prevent HGCE 21/02/22
!      ICASE = KTYPE(K1,K2)
*
*     if(icase.gt.100) WRITE(66,*)' MIX ERROR ICASE>100 ',icase,k1,k2
*
*       Determine evolution time scales for first star.
      M01 = M0(I1)
      M1 = M(I1)
      AGE1 = AJ(I1)
      CALL star(K1,M01,M1,TMS1,TN,TSCLS,LUMS,GB,ZPARS)
*
*       Obtain time scales for second star.
      M02 = M0(I2)
      M2 = M(I2)
      AGE2 = AJ(I2)
      CALL star(K2,M02,M2,TMS2,TN,TSCLS,LUMS,GB,ZPARS)
*
*       Check for planetary systems - defined as HeWDs and low-mass WDs!
      IF(K1.EQ.10.AND.M1.LT.0.05)THEN
         ICASE = K2
         IF(K2.LE.1)THEN
            ICASE = 1
            AGE1 = 0.D0
         ENDIF
      ELSEIF(K1.GE.11.AND.M1.LT.0.5.AND.ICASE.EQ.6)THEN
         ICASE = 9
      ENDIF
      IF(K2.EQ.10.AND.M2.LT.0.05)THEN
         ICASE = K1
         IF(K1.LE.1)THEN
            ICASE = 1
            AGE2 = 0.D0
         ENDIF
      ENDIF
*
*       Specify total mass.
      M3 = M1 + M2
      M03 = M01 + M02
      KW = ICASE
      AGE3 = 0.d0
*
*       Restrict merged stars to masses less than 100 Msun. 
      if(askInUseOrNot()) then
      else
         IF(M3.GE.100.D0)THEN
            M3 = 99.D0
            M03 = MIN(M03,M3)
         ENDIF
      endif
*
*       Evaluate apparent age and other parameters.
*
      IF(ICASE.EQ.1)THEN
*       Specify new age based on complete mixing.
         IF(K1.EQ.7) KW = 7
         CALL star(KW,M03,M3,TMS3,TN,TSCLS,LUMS,GB,ZPARS)
         !AGE3 = 0.1d0*TMS3*(AGE1*M1/TMS1 + AGE2*M2/TMS2)/M3
         AGE3 = TMS3*(AGE1*M1/TMS1 + AGE2*M2/TMS2)/M3
* A. Tanikawa adds here to prevent HGCE 21/02/22
      ELSEIF(ICASE.EQ.2)THEN
         CALL star(KW,M03,M3,TMS3,TN,TSCLS,LUMS,GB,ZPARS)
         AGE3 = 0.1d0*TMS3*(AGE1*M1/TMS1 + AGE2*M2/TMS2)/M3
*
      ELSEIF(ICASE.EQ.3.OR.ICASE.EQ.6.OR.ICASE.EQ.9)THEN
         MC3 = M1
         CALL gntage(MC3,M3,KW,ZPARS,M03,AGE3)
* A. Tanikawa changes here to prevent HGCE 21/02/22
!      ELSEIF(ICASE.EQ.4)THEN
      ELSEIF(ICASE.EQ.4.or.ICASE.EQ.5)THEN
*
         MC3 = M1
         AGE3 = AGE1/TMS1
         CALL gntage(MC3,M3,KW,ZPARS,M03,AGE3)
      ELSEIF(ICASE.EQ.7)THEN
         CALL star(KW,M03,M3,TMS3,TN,TSCLS,LUMS,GB,ZPARS)
         AGE3 = TMS3*(AGE2*M2/TMS2)/M3
      ELSEIF(ICASE.LE.12)THEN
*       Ensure that a new WD has the initial mass set correctly.
         M03 = M3
         IF(ICASE.LT.12.AND.M3.GE.MCH)THEN
            M3 = 0.D0
            KW = 15
         ENDIF
      ELSEIF(ICASE.EQ.13.OR.ICASE.EQ.14)THEN
*       Set unstable Thorne-Zytkow object with fast mass loss of envelope 
*       unless the less evolved star is a WD, NS or BH. 
         IF(K2.LT.10)THEN
            M03 = M1
            M3 = M1
         ENDIF
         IF(ICASE.EQ.13.AND.M3.GT.MXNS) KW = 14
      ELSEIF(ICASE.EQ.15)THEN
         M3 = 0.D0
      ELSEIF(ICASE.GT.100)THEN
*       Common envelope case which should only be used after COMENV.
         KW = K1
         AGE3 = AGE1
         M3 = M1
         M03 = M01
      ELSE
*       This should not be reached.
        KW = 1
        M03 = M3
      ENDIF
* A. Tanikawa adds this subroutine 21/10/14
      if(kw.le.1) then
         if(m(1).gt.m(2)) then
            q = m(2)/m(1)
         else
            q = m(1)/m(2)
         endif
         call mixquantity(tphys,kw,m03,age3,m3,q,
     &        tmrg,tcon,tkhq,rmrg,emrg)
      endif
*
*
* Put the result in *1.
*
      KS(1) = KW
      KS(2) = 15
      M(1) = M3
      M(2) = 0.D0
      M0(1) = M03
      AJ(1) = AGE3
*
      RETURN
      END
***

* A. Tanikawa adds this subroutine 21/10/14
      subroutine mixquantity0(tphys,kw,m0,aj,mt,q,tmrg,tcon,rmrg,emrg)
      use iso_c_binding
      implicit none
*
      real*8 TkhSunMyr
      parameter (TkhSunMyr=1.570d1)
      real*8 qc
      parameter (qc=3.0d-1)
      real*8 tphys, q
      integer kw
      real*8 m0,aj,mt,tm,tn,tscls(20),lums(10),gb(10),zpars(20)
      real*8 r,lum,mc,rc,menv,renv,k2,ffb,tsn
      real*8 jspin,arem
      real*8 ragbf
      real*8 rg
      real*8 tmrg(2), tcon(2), rmrg(2), emrg(2)
      real*8 tmrg3, rmrg3, tcon3, emrg3
      real*8 uagb, ubse, umrg
      include 'cppinterface.h'
*
      tmrg3 = -1.
      call star(kw,m0,mt,tm,tn,tscls,lums,gb,zpars)
      call hrdiag(m0,aj,mt,tm,tn,tscls,lums,gb,zpars,
     &     r,lum,kw,mc,rc,menv,renv,k2,ffb,tsn,
     &     jspin,arem,
     &     tphys,tmrg3,0.,0.)
      tmrg3 = tphys
      if(askInUseOrNot()
     &     .and. askInScopeOfApplication(m0)) then
         rg = 10**(getRadiusRedPhase(mt,lum))
      else
         rg = ragbf(mt,lum,zpars(2))
      endif
      ubse   = - mt**2/(2.*r)
      uagb   = - mt**2/(2.*rg)
      emrg3 = uagb - ubse
      if(q .lt. qc) then
         emrg3 = emrg3 * ((qc/q)*uagb-ubse)/(uagb-ubse)
      endif
      umrg = ubse + emrg3
      if(tmrg(1).ge.0.) then
         umrg  = umrg  + (1.-min(1.,(tphys-tmrg(1))/tcon(1)))*emrg(1)
         emrg3 = emrg3 + (1.-min(1.,(tphys-tmrg(1))/tcon(1)))*emrg(1)
      endif
      if(tmrg(2).ge.0.) then
         umrg  = umrg  + (1.-min(1.,(tphys-tmrg(2))/tcon(2)))*emrg(2)
         emrg3 = emrg3 + (1.-min(1.,(tphys-tmrg(2))/tcon(2)))*emrg(2)
      endif
      umrg  = max(ubse,min(umrg,uagb))
      rmrg3 = - mt**2/(2.*umrg)
      tcon3 = max(0.,(1.-(r/rmrg3))*mt**2/(r*lum)*TkhSunMyr)
      tmrg(1) = tmrg3
      tmrg(2) = -1.
      tcon(1) = tcon3
      tcon(2) = -1.
      rmrg(1) = rmrg3
      rmrg(2) = -1.
      emrg(1) = emrg3
      emrg(2) = -1.
*      
      return
      end
*

* A. Tanikawa adds this subroutine 21/10/14
      subroutine mixquantity(tphys,kw,m0,aj,mt,q,
     &     tmrg,tcon,tkhq,rmrg,emrg)
      use iso_c_binding
      implicit none
*
      real*8 TkhSunMyr
      parameter (TkhSunMyr=1.570d1)
      real*8 qc1,qc2,fqc2
      parameter (qc1=3.0d-1,fqc2=0.1)
      real*8 tphys, q
      integer kw,kwold
      real*8 m0,aj,mt,tm,tn,tscls(20),lums(10),gb(10),zpars(20)
      real*8 r,lum,mc,rc,menv,renv,k2,ffb,tsn
      real*8 jspin,arem
      real*8 ragbf
      real*8 rg
      real*8 tmrg(2), tcon(2), tkhq(2), rmrg(2), emrg(2)
      real*8 tmrg3, rmrg3, tcon3, emrg3
      real*8 uagb, ubse, umrg
      include 'cppinterface.h'
*
      tmrg3 = -1.
      kwold = kw
      call star(kw,m0,mt,tm,tn,tscls,lums,gb,zpars)
      call hrdiag(m0,aj,mt,tm,tn,tscls,lums,gb,zpars,
     &     r,lum,kw,mc,rc,menv,renv,k2,ffb,tsn,
     &     jspin,arem,
     &     tphys,tmrg3,0.,0.)
      tmrg3 = tphys
      if(askInUseOrNot()
     &     .and. askInScopeOfApplication(m0)) then
         rg = 10**(getRadiusRedPhase(mt,lum))
      else
         rg = ragbf(mt,lum,zpars(2))
      endif
      ubse  = - mt**2/(2.*r)
      uagb  = - mt**2/(2.*rg)
      qc2   = qc1 * uagb / ((fqc2 * (uagb-ubse) + ubse))
      emrg3 = uagb - ubse
      if(q .lt. qc2) then
         emrg3 = emrg3 * ((qc1/qc2)*uagb-ubse)/(uagb-ubse) * q/qc2
      elseif(q .lt. qc1) then
         emrg3 = emrg3 * ((qc1/q)*uagb-ubse)/(uagb-ubse)
      endif
      umrg = ubse + emrg3
      if(tmrg(1).ge.0. 
     &     .and. 1.-(tphys-tmrg(1))/(tkhq(1)) .ge. 0) then
         umrg  = umrg  + emrg(1)
         emrg3 = emrg3 + emrg(1)
      endif
      if(tmrg(2).ge.0.
     &     .and. 1.-(tphys-tmrg(2))/(tkhq(2)) .ge. 0) then
         umrg  = umrg  + emrg(2)
         emrg3 = emrg3 + emrg(2)
      endif
      umrg  = min(umrg,uagb)
      if(umrg .lt. ubse .and. kw.eq.kwold) then
         write(0,*) 'Error: umrg < ubse in mixquantity.'
         write(0,*) 'ubse umrg uagb', ubse, umrg, uagb
         stop
      endif
      rmrg3 = - mt**2/(2.*umrg)
      tcon3 = max(0.,(1.-(r/rmrg3))*mt**2/(r*lum)*TkhSunMyr)
      tmrg(1) = tmrg3
      tmrg(2) = -1.
      tcon(1) = tcon3
      tcon(2) = -1.
      tkhq(1) = mt**2/(r*lum)*TkhSunMyr*q
      tkhq(2) = -1.
      rmrg(1) = rmrg3
      rmrg(2) = -1.
      emrg(1) = emrg3
      emrg(2) = -1.
      if(kw.ne.kwold) then
         tmrg(1) = -1.
         tcon(1) = -1.
         tkhq(1) = -1.
         rmrg(1) = -1.
         emrg(1) = -1.
      endif
*      
      return
      end
*
