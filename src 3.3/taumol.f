C      path:      $Source$
C     author:    $Author: kcadyper $
C     revision:  $Revision: 22602 $
C     created:   $Date: 2013-11-11 11:21:35 -0500 (Mon, 11 Nov 2013) $
C
C  --------------------------------------------------------------------------
C |                                                                          |
C |  Copyright 2002, 2003, Atmospheric & Environmental Research, Inc. (AER). |
C |  This software may be used, copied, or redistributed as long as it is    |
C |  not sold and this copyright notice is reproduced on each copy made.     |
C |  This model is provided as is without any express or implied warranties. |
C |                       (http://www.rtweb.aer.com/)                        |
C |                                                                          |
C  --------------------------------------------------------------------------

*******************************************************************************
*                                                                             *
*                  Optical depths developed for the                           *
*                                                                             *
*                RAPID RADIATIVE TRANSFER MODEL (RRTM)                        *
*                                                                             *
*                                                                             *
*            ATMOSPHERIC AND ENVIRONMENTAL RESEARCH, INC.                     *
*                        131 HARTWELL AVENUE                                  *
*                        LEXINGTON, MA 02421                                  *
*                                                                             *
*                                                                             *
*                           ELI J. MLAWER                                     * 
*                         JENNIFER DELAMERE                                   * 
*                         STEVEN J. TAUBMAN                                   *
*                         SHEPARD A. CLOUGH                                   *
*                                                                             *
*                                                                             *
*                                                                             *
*                                                                             *
*                       email:  mlawer@aer.com                                *
*                       email:  jdelamer@aer.com                              *
*                                                                             *
*        The authors wish to acknowledge the contributions of the             *
*        following people:  Karen Cady-Pereira, Patrick D. Brown,             *  
*        Michael J. Iacono, Ronald E. Farren, Luke Chen, Robert Bergstrom.    *
*                                                                             *
*******************************************************************************
*     TAUMOL                                                                  *
*                                                                             *
*     This file contains the subroutines TAUGBn (where n goes from            *
*     1 to 16).  TAUGBn calculates the optical depths and Planck fractions    *
*     per g-value and layer for band n.                                       *
*                                                                             *
*  Output:  optical depths (unitless)                                         *
*           fractions needed to compute Planck functions at every layer       *
*               and g-value                                                   *
*                                                                             *
*     COMMON /TAUGCOM/  TAUG(MXLAY,MG)                                        *
*     COMMON /PLANKG/   FRACS(MXLAY,MG)                                       *
*                                                                             *
*  Input                                                                      *
*                                                                             *
*     COMMON /FEATURES/ NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)                  *
*     COMMON /PRECISE/  ONEMINUS                                              *
*     COMMON /PROFILE/  NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),                    *
*     &                 PZ(0:MXLAY),TZ(0:MXLAY)                               *
*     COMMON /PROFDATA/ LAYTROP,                                              *
*    &                  COLH2O(MXLAY),COLCO2(MXLAY),COLO3(MXLAY),             *
*    &                  COLN2O(MXLAY),COLCO(MXLAY),COLCH4(MXLAY),             *
*    &                  COLO2(MXLAY)
*     COMMON /INTFAC/   FAC00(MXLAY),FAC01(MXLAY),                            *
*    &                  FAC10(MXLAY),FAC11(MXLAY)                             *
*     COMMON /INTIND/   JP(MXLAY),JT(MXLAY),JT1(MXLAY)                        *
*     COMMON /SELF/     SELFFAC(MXLAY), SELFFRAC(MXLAY), INDSELF(MXLAY)       *
*                                                                             *
*     Description:                                                            *
*     NG(IBAND) - number of g-values in band IBAND                            *
*     NSPA(IBAND) - for the lower atmosphere, the number of reference         *
*                   atmospheres that are stored for band IBAND per            *
*                   pressure level and temperature.  Each of these            *
*                   atmospheres has different relative amounts of the         *
*                   key species for the band (i.e. different binary           *
*                   species parameters).                                      *
*     NSPB(IBAND) - same for upper atmosphere                                 *
*     ONEMINUS - since problems are caused in some cases by interpolation     *
*                parameters equal to or greater than 1, for these cases       *
*                these parameters are set to this value, slightly < 1.        *
*     PAVEL - layer pressures (mb)                                            *
*     TAVEL - layer temperatures (degrees K)                                  *
*     PZ - level pressures (mb)                                               *
*     TZ - level temperatures (degrees K)                                     *
*     LAYTROP - layer at which switch is made from one combination of         *
*               key species to another                                        *
*     COLH2O, COLCO2, COLO3, COLN2O, COLCH4 - column amounts of water         *
*               vapor,carbon dioxide, ozone, nitrous ozide, methane,          *
*               respectively (molecules/cm**2)                                *

*     FACij(LAY) - for layer LAY, these are factors that are needed to        *
*                  compute the interpolation factors that multiply the        *
*                  appropriate reference k-values.  A value of 0 (1) for      *
*                  i,j indicates that the corresponding factor multiplies     *
*                  reference k-value for the lower (higher) of the two        *
*                  appropriate temperatures, and altitudes, respectively.     *
*     JP - the index of the lower (in altitude) of the two appropriate        *
*          reference pressure levels needed for interpolation                 *
*     JT, JT1 - the indices of the lower of the two appropriate reference     *
*               temperatures needed for interpolation (for pressure           *
*               levels JP and JP+1, respectively)                             *
*     SELFFAC - scale factor needed for water vapor self-continuum, equals    *
*               (water vapor density)/(atmospheric density at 296K and        *
*               1013 mb)                                                      *
*     SELFFRAC - factor needed for temperature interpolation of reference     *
*                water vapor self-continuum data                              *
*     INDSELF - index of the lower of the two appropriate reference           *
*               temperatures needed for the self-continuum interpolation      *
*     FORFAC  - scale factor needed for water vapor foreign-continuum.        *
*     FORFRAC - factor needed for temperature interpolation of reference      *
*                water vapor foreign-continuum data                           *
*     INDFOR  - index of the lower of the two appropriate reference           *
*               temperatures needed for the foreign-continuum interpolation   *
*                                                                             *
*  Data input                                                                 *
*     COMMON /Kn/ KA(NSPA(n),5,13,MG), KB(NSPB(n),5,13:59,MG), SELFREF(10,MG),*
*                 FORREF(4,MG), KA_M'MGAS', KB_M'MGAS'                        *
*        (note:  n is the band number,'MGAS' is the species name of the minor *
*         gas)                                                                *
*                                                                             *
*     Description:                                                            *
*     KA - k-values for low reference atmospheres (key-species only)          *
*          (units: cm**2/molecule)                                            *
*     KB - k-values for high reference atmospheres (key-species only)         *
*          (units: cm**2/molecule)                                            *
*     KA_M'MGAS' - k-values for low reference atmosphere minor species        *
*          (units: cm**2/molecule)                                            *
*     KB_M'MGAS' - k-values for high reference atmosphere minor species       *
*          (units: cm**2/molecule)                                            *
*     SELFREF - k-values for water vapor self-continuum for reference         *
*               atmospheres (used below LAYTROP)                              *
*               (units: cm**2/molecule)                                       *
*     FORREF  - k-values for water vapor foreign-continuum for reference      *
*               atmospheres (used below/above LAYTROP)                        *
*               (units: cm**2/molecule)                                       *
*                                                                             *
*     DIMENSION ABSA(65*NSPA(n),MG), ABSB(235*NSPB(n),MG)                     *
*     EQUIVALENCE (KA,ABSA),(KB,ABSB)                                         *
*                                                                             *
*******************************************************************************

      SUBROUTINE TAUGB1A

      use variables !NJE

C     Written by Eli J. Mlawer, Atmospheric & Environmental Research.

C     BAND 1:  10-120 cm-1 (key - H2O, CO2; low minor - N2)

C     NOTE: Previous versions of RRTM BAND 1: 
C           10-250 cm-1 (low - H2O; high - H2O)


      PARAMETER (mg=16, mxlay=203, MXMOL = 38, NBANDS=16)
      PARAMETER (NTEMP=14,NETA=9,NREF=71,NCO2=2)

C  Output

      COMMON /TAUGCOM/  TAUG(MXLAY,MG)
      COMMON /PLANKG/   FRACS(MXLAY,MG)                                       

C  Input

      COMMON /FEATURES/ NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)
      COMMON /PRECISE/  ONEMINUS
      COMMON /PROFILE/  NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
     &                  PZ(0:MXLAY),TZ(0:MXLAY)
      COMMON /SPECIES/  COLDRY(MXLAY),WKL(MXMOL,MXLAY),WBROAD(MXLAY),
     &                  COLMOL(MXLAY),NMOL
      COMMON /PROFDATA/ LAYTROP, LAY_SWITCH_KABS,                                  
     &                  COLH2O(MXLAY),COLCO2(MXLAY),COLO3(MXLAY),  
     &                  COLN2O(MXLAY),COLCO(MXLAY),COLCH4(MXLAY),  
     &                  COLO2(MXLAY),COLBRD(MXLAY)
      COMMON /INTFAC/   FAC00(MXLAY),FAC01(MXLAY),                            
     &                  FAC10(MXLAY),FAC11(MXLAY)                             
      COMMON /CO2FAC/   co2_fac0(mxlay),co2_fac1(mxlay)
      COMMON /STD_REF/  PREF(nref),PREFLOG(nref),TREF(nref),
     &                  CHI_std(7,nref)
      COMMON /REFRAT_ETA/ RAT_H2OCO2(MXLAY),RAT_H2OCO2_1(MXLAY),
     &                  RAT_H2OO3(MXLAY),RAT_H2OO3_1(MXLAY),
     &                  RAT_H2ON2O(MXLAY),RAT_H2ON2O_1(MXLAY),
     &                  RAT_H2OCH4(MXLAY),RAT_H2OCH4_1(MXLAY),
     &                  RAT_N2OCO2(MXLAY),RAT_N2OCO2_1(MXLAY),
     &                  RAT_O3CO2(MXLAY),RAT_O3CO2_1(MXLAY),
     &                  RAT_CO2O3(MXLAY),RAT_CO2O3_1(MXLAY),
     &                  RAT_O2_DRY(MXLAY)

      COMMON /INTIND/   JP(MXLAY),JT(MXLAY),JT1(MXLAY)
      COMMON /SELF/     SELFFAC(MXLAY), SELFFRAC(MXLAY), INDSELF(MXLAY)
      COMMON /MINOR/    MINORFRAC(MXLAY), INDMINOR(MXLAY), 
     &                  SCALEMINOR(MXLAY),SCALEMINORN2(MXLAY)
      COMMON /FOREIGN/  FORFAC(MXLAY), FORFRAC(MXLAY), INDFOR(MXLAY)
      COMMON /CO2/      INDCO2(MXLAY)
      COMMON /K1A/       KA(NETA,NTEMP,NREF,MG),  FORREF(4,MG),
     &                  SELFREF(10,MG), K_CO2(NETA,NTEMP,NREF,MG),
     &                  KA_MN2(NETA,19,MG),KA_MN2_O2(NETA,19,MG)

      COMMON /CVRTAU/    HNAMTAU,HVRTAU

      CHARACTER*18       HNAMTAU,HVRTAU

      DIMENSION ABSA(8946,MG),ABSB(8946,MG),ABSC(8946,MG)
      DIMENSION FRACREFA(MG,9)

      integer js(nco2),js1(nco2),ind0(nco2),ind1(nco2)
      real fs(nco2),fs1(nco2)

      real fac000(nco2),fac100(nco2),fac200(nco2)
      real fac010(nco2),fac110(nco2),fac210(nco2)
      real fac001(nco2),fac101(nco2),fac201(nco2)
      real fac011(nco2),fac111(nco2),fac211(nco2)

      real specparm(nco2),specparm1(nco2)
      real speccomb(nco2),speccomb1(nco2)


      real tau_base(mg)
      real tau_major(nco2),tau_major1(nco2)

      real n2m1, n2m2, n2_o2_m1, n2_o2_m2



C Minor gas mapping levels:
C     LOWER - N2, P = 142.5490 mbar, T = 215.70 K

      EQUIVALENCE (KA,ABSB)
      EQUIVALENCE (K_CO2,ABSC)
      REAL KA, KB, K_CO2, KA_MN2,  KA_MN2_O2, MINORFRAC

C     Compute the optical depth by interpolating in ln(pressure) and 
C     temperature.  Below LAYTROP, the water vapor self-continuum and
C     foreign continuum is interpolated (in temperature) separately.



C Planck fraction mapping level: P = 212.7250 mbar, T = 223.06 K

      DATA (FRACREFA(IG, 1),IG=1,16) /
     &3.3477E-01,2.5211E-01,1.7614E-01,1.1300E-01,6.5709E-02,3.4091E-02,
     &1.5507E-02,6.0911E-03,2.1143E-03,1.5333E-04,1.1695E-04,8.6326E-05,
     &5.9999E-05,3.6776E-05,1.3580E-05,1.8984E-06/
      DATA (FRACREFA(IG, 2),IG=1,16) /
     &3.2715E-01,2.3040E-01,1.5374E-01,9.5000E-02,5.7570E-02,2.8449E-02,
     &1.5650E-02,2.3792E-02,4.6803E-02,5.9431E-03,4.7199E-03,4.2703E-03,
     &3.4534E-03,2.0563E-03,8.6882E-04,1.3891E-04/
      DATA (FRACREFA(IG, 3),IG=1,16) /
     &3.2283E-01,2.1888E-01,1.4098E-01,8.9186E-02,4.8561E-02,2.9305E-02,
     &1.4610E-02,6.4264E-02,4.9761E-02,6.0255E-03,4.7799E-03,4.2862E-03,
     &3.4627E-03,2.0624E-03,8.6896E-04,1.3891E-04/
      DATA (FRACREFA(IG, 4),IG=1,16) /
     &3.1868E-01,1.9795E-01,1.2850E-01,8.7230E-02,4.5005E-02,2.7377E-02,
     &5.1296E-02,7.1593E-02,5.0697E-02,6.0446E-03,4.8001E-03,4.2940E-03,
     &3.4646E-03,2.0647E-03,8.6902E-04,1.3891E-04/
      DATA (FRACREFA(IG, 5),IG=1,16) /
     &3.1081E-01,1.8930E-01,1.0491E-01,7.1181E-02,4.8246E-02,4.0642E-02,
     &8.7630E-02,7.4400E-02,5.1178E-02,6.0565E-03,4.8099E-03,4.2965E-03,
     &3.4664E-03,2.0656E-03,8.6902E-04,1.3891E-04/
      DATA (FRACREFA(IG, 6),IG=1,16) /
     &3.0589E-01,1.6094E-01,9.9906E-02,5.3576E-02,3.4738E-02,9.8089E-02,
     &9.7484E-02,7.6194E-02,5.1456E-02,6.0624E-03,4.8159E-03,4.2984E-03,
     &3.4672E-03,2.0663E-03,8.6902E-04,1.3891E-04/
      DATA (FRACREFA(IG, 7),IG=1,16) /
     &2.8064E-01,1.3610E-01,8.3887E-02,4.7695E-02,8.5042E-02,1.1333E-01,
     &1.0253E-01,7.7404E-02,5.1639E-02,6.0661E-03,4.8189E-03,4.3003E-03,
     &3.4672E-03,2.0678E-03,8.6896E-04,1.3891E-04/
      DATA (FRACREFA(IG, 8),IG=1,16) /
     &2.0584E-01,6.8266E-02,6.6115E-02,1.5524E-01,1.2614E-01,1.2129E-01,
     &1.0537E-01,7.8219E-02,5.1771E-02,6.0692E-03,4.8224E-03,4.3005E-03,
     &3.4683E-03,2.0680E-03,8.6896E-04,1.3891E-04/
      DATA (FRACREFA(IG, 9),IG=1,16) /
     &1.7722E-02,9.5719E-02,1.7511E-01,1.7510E-01,1.4952E-01,1.2667E-01,
     &1.0771E-01,7.8828E-02,5.1870E-02,6.0696E-03,4.8261E-03,4.3009E-03,
     &3.4689E-03,2.0674E-03,8.6896E-04,1.3891E-04/

      co2_lo_hi_rat = 0.04*3.85e-4/0.95

      HVRTAU = '$Revision: 22602 $'

      idump = 1



C     Calculate reference ratio to be used in calculation of Planck
C     fraction 

      base_rat = 15./9999.
      iplanck_lay = 21

      REFRAT_PLANCK_A = CHI_STD(1,iplanck_lay)/CHI_STD(2,iplanck_lay)
      if (pref(iplanck_lay).le. 120.) then
	 rnew_rat = base_rat
      else if (pref(iplanck_lay).le.800) then 
	 rnew_rat = base_rat*(pref(iplanck_lay)/120.)**(-3.7)
      else
	 rnew_rat= base_rat*(800./120.)**(-3.7)
      endif
      refrat_planck_a = refrat_planck_a*rnew_rat

C P = 142.549 mb
      REFRAT_M_A = CHI_STD(1,23)/CHI_STD(2,23)

C     Compute the optical depth by interpolating in ln(pressure), 
C     temperature, and appropriate species.  Below LAYTROP, the water
C     vapor self-continuum and foreign continuum is interpolated 
C     (in temperature) separately.  

      HVRTAU = '$Revision: 22602 $'
      base_rat = 15./9999.

      DO 2500 LAY = 1, NLAYERS
         if (pavel(lay).le. 120.) then
            rnew_rat = base_rat
         else if (pavel(lay).le.800) then 
            rnew_rat = base_rat*(pavel(lay)/120.)**(-3.7)
         else
            rnew_rat= base_rat*(800./120.)**(-3.7)
         endif

c Set up calculation for low and/or high CO2 values
         select case (indco2(lay))
         case (0)
             ico2_0 = 1
             ico2_1 = 1
         case(2)
             ico2_0 = 2
             ico2_1 = 2
         case(1)
             ico2_0 = 1
             ico2_1 = 2
          case default
             print *, 'indco2 in taugb1a has an invalid value'
             stop
         end select

        do ic=ico2_0,ico2_1
         if (ic.eq.1) then
            scale_eta=rnew_rat 
         else 
            scale_eta= co2_lo_hi_rat
         end if

         SPECCOMB(ic) = COLH2O(LAY) + scale_eta*
     &           RAT_H2OCO2(LAY)*COLCO2(LAY)
         SPECPARM(ic) = COLH2O(LAY)/SPECCOMB(ic)
         IF (SPECPARM(ic) .GE. ONEMINUS) SPECPARM(ic) = ONEMINUS
         SPECMULT = 8.*(SPECPARM(ic))
         JS(ic) = 1 + INT(SPECMULT)
         FS(ic) = AMOD(SPECMULT,1.0)

         SPECCOMB1(ic) = COLH2O(LAY)+scale_eta
     &        *RAT_H2OCO2_1(LAY)*COLCO2(LAY)
         SPECPARM1(ic) = COLH2O(LAY)/SPECCOMB1(ic)
         IF (SPECPARM1(ic) .GE. ONEMINUS) SPECPARM1(ic) = ONEMINUS
         SPECMULT1 = 8.*(SPECPARM1(ic))
         JS1(ic) = 1 + INT(SPECMULT1)
         FS1(ic) = AMOD(SPECMULT1,1.0)

         SPECCOMB_PLANCK = COLH2O(LAY)+
     &         REFRAT_PLANCK_A*COLCO2(LAY)
     &                *(pavel(lay)/pref(iplanck_lay))
         SPECPARM_PLANCK = COLH2O(LAY)/SPECCOMB_PLANCK
         IF (SPECPARM_PLANCK .GE. ONEMINUS) SPECPARM_PLANCK=ONEMINUS
         SPECMULT_PLANCK = 8.*SPECPARM_PLANCK
         JPL= 1 + INT(SPECMULT_PLANCK)
         FPL = AMOD(SPECMULT_PLANCK,1.0)

         SPECCOMB_MN2 = COLH2O(LAY) + REFRAT_M_A*COLCO2(LAY)
         SPECPARM_MN2 = COLH2O(LAY)/SPECCOMB_MN2
         IF (SPECPARM_MN2 .GE. ONEMINUS) SPECPARM_MN2 = ONEMINUS
         SPECMULT_MN2 = 8.*SPECPARM_MN2
         JMN2 = 1 + INT(SPECMULT_MN2)
         FMN2 = AMOD(SPECMULT_MN2,1.0)
         FMN2MF = MINORFRAC(LAY)*FMN2

         IND0(ic) = ((JP(LAY)-1)*NTEMP+(JT(LAY)-1))*NSPA(1) + JS(ic)
         IND1(ic) = (JP(LAY)*NTEMP+(JT1(LAY)-1))*NSPA(1) + JS1(ic)
         INDS = INDSELF(LAY)
         INDF = INDFOR(LAY)
         INDM = INDMINOR(LAY)

         IF (SPECPARM(ic) .LT. 0.125) THEN
            P = FS(ic) - 1
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC000(ic) = FK0*FAC00(LAY)
            FAC100(ic) = FK1*FAC00(LAY)
            FAC200(ic) = FK2*FAC00(LAY)
            FAC010(ic) = FK0*FAC10(LAY)
            FAC110(ic) = FK1*FAC10(LAY)
            FAC210(ic) = FK2*FAC10(LAY)
         ELSE IF (SPECPARM(ic) .GT. 0.875) THEN
            P = -FS(ic)
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC000(ic) = FK0*FAC00(LAY)
            FAC100(ic) = FK1*FAC00(LAY)
            FAC200(ic) = FK2*FAC00(LAY)
            FAC010(ic) = FK0*FAC10(LAY)
            FAC110(ic) = FK1*FAC10(LAY)
            FAC210(ic) = FK2*FAC10(LAY)
	     corradj = 1.
	     wv_log = alog10(wkl(1,lay)/coldry(lay))
	     if(wv_log.gt.-5.AND.wv_log.le.-3.0) then
	       corradj = 1.+0.05*(wv_log+3.0)
	     else if (wv_log.gt.-5.5.AND.wv_log.le.-5.) then
	       corradj = 0.9
	     else if (wv_log.gt.-6.5.AND.wv_log.le.-5.5) then
	       corradj = 1.-0.1*(wv_log+6.5)
	     endif
	     corradj = 1.
         ELSE
            FAC000(ic) = (1. - FS(ic)) * FAC00(LAY)
            FAC010(ic) = (1. - FS(ic)) * FAC10(LAY)
            FAC100(ic) = FS(ic) * FAC00(LAY)
            FAC110(ic) = FS(ic) * FAC10(LAY)
         END IF

       
         IF (SPECPARM1(ic) .LT. 0.125 ) THEN
            P = FS1(ic) - 1
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC001(ic) = FK0*FAC01(LAY)
            FAC101(ic) = FK1*FAC01(LAY)
            FAC201(ic) = FK2*FAC01(LAY)
            FAC011(ic) = FK0*FAC11(LAY)
            FAC111(ic) = FK1*FAC11(LAY)
            FAC211(ic) = FK2*FAC11(LAY)
         ELSE IF (SPECPARM1(ic) .GT. 0.875 ) THEN
            P = -FS1(ic)
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC001(ic) = FK0*FAC01(LAY)
            FAC101(ic) = FK1*FAC01(LAY)
            FAC201(ic) = FK2*FAC01(LAY)
            FAC011(ic) = FK0*FAC11(LAY)
            FAC111(ic) = FK1*FAC11(LAY)
            FAC211(ic) = FK2*FAC11(LAY)
              corradj1 = 1.
              wv_log = alog10(wkl(1,lay)/coldry(lay))
	     if(wv_log.gt.-5.AND.wv_log.le.-3.0) then
	       corradj1 = 1.+0.05*(wv_log+3.0)
	     else if (wv_log.gt.-5.5.AND.wv_log.le.-5.) then
	       corradj1 = 0.9
	     else if (wv_log.gt.-6.5.AND.wv_log.le.-5.5) then
	       corradj1 = 1.-0.1*(wv_log+6.5)
	     endif
              corradj1 = 1.
         ELSE
            FAC001(ic) = (1. - FS1(ic)) * FAC01(LAY)
            FAC011(ic) = (1. - FS1(ic)) * FAC11(LAY)
            FAC101(ic) = FS1(ic) * FAC01(LAY)
            FAC111(ic) = FS1(ic) * FAC11(LAY)
         END IF

        end do   ! end ic loop for general layer properties

    
         DO 2000 IG = 1, NG(1)
	   DO IC=ico2_0,ico2_1
	    if(ic.eq.1) then
	       absa = absb
	    else
		absa = absc
	    endif 

            IF (SPECPARM(ic).LT.0.125) THEN 
               TAU_MAJOR(ic) = SPECCOMB(ic) *
     &          (FAC000(ic)* ABSA(ind0(ic),IG) +
     &          FAC100(ic)* ABSA(ind0(ic)+1,IG) +
     &          FAC200(ic)* ABSA(ind0(ic)+2,IG) +
     &          FAC010(ic)* ABSA(ind0(ic)+9,IG) +
     &          FAC110(ic)* ABSA(ind0(ic)+10,IG) +
     &          FAC210(ic)* ABSA(ind0(ic)+11,IG))
           ELSE IF (SPECPARM(ic).GT.0.875) THEN 
                TAU_MAJOR(ic) = CORRADJ*speccomb(ic) * 
     &              (FAC200(ic)* ABSA(ind0(ic)-1,IG) +
     &              FAC100(ic)* ABSA(ind0(ic),IG) +
     &              FAC000(ic)* ABSA(ind0(ic)+1,IG) +
     &              FAC210(ic)* ABSA(ind0(ic)+8,IG) +
     &              FAC110(ic)* ABSA(ind0(ic)+9,IG) +
     &              FAC010(ic)* ABSA(ind0(ic)+10,IG))
           ELSE 
               TAU_MAJOR(ic) = speccomb(ic) * 
     &              (FAC000(ic) * ABSA(IND0(ic),IG) +
     &              FAC100(ic) * ABSA(IND0(ic)+1,IG) +
     &              FAC010(ic) * ABSA(IND0(ic)+9,IG) +
     &              FAC110(ic) * ABSA(IND0(ic)+10,IG))
           END IF

           IF (SPECPARM1(ic).LT.0.125) THEN 
                TAU_MAJOR1(ic) =  speccomb1(ic) *
     &              (FAC001(ic) * ABSA(IND1(ic),IG) + 
     &              FAC101(ic) * ABSA(IND1(ic)+1,IG) +
     &              FAC201(ic) * ABSA(IND1(ic)+2,IG) +
     &              FAC011(ic) * ABSA(IND1(ic)+9,IG) +
     &              FAC111(ic) * ABSA(IND1(ic)+10,IG) +
     &              FAC211(ic) * ABSA(IND1(ic)+11,IG))

           ELSE IF (SPECPARM1(ic).GT.0.875) THEN 
                TAU_MAJOR1(ic) = CORRADJ1*speccomb1(ic) * 
     &              (FAC201(ic) * ABSA(IND1(ic)-1,IG) +
     &              FAC101(ic) * ABSA(IND1(ic),IG) +
     &              FAC001(ic) * ABSA(IND1(ic)+1,IG) +
     &              FAC211(ic) * ABSA(IND1(ic)+8,IG) +
     &              FAC111(ic) * ABSA(IND1(ic)+9,IG) +
     &              FAC011(ic) * ABSA(IND1(ic)+10,IG))
           ELSE
               TAU_MAJOR1(ic) = speccomb1(ic) *
     &              (FAC001(ic) * ABSA(IND1(ic),IG) +
     &              FAC101(ic) * ABSA(IND1(ic)+1,IG) +
     &              FAC011(ic) * ABSA(IND1(ic)+9,IG) +
     &              FAC111(ic) * ABSA(IND1(ic)+10,IG))

           END IF
          end do    ! co2 loop

C  Do  interpolation for high CO2 cases if needed
          select case (indco2(lay) )
          case (0)
             TAUG(LAY,IG) = TAU_MAJOR(1)+TAU_MAJOR1(1)
          case (1)
             TAUG(LAY,IG) = co2_fac0(lay)*(TAU_MAJOR(1)+TAU_MAJOR1(1)) + 
     &                      co2_fac1(lay)*(TAU_MAJOR(2)+TAU_MAJOR1(2)) 
          case (2)
             TAUG(LAY,IG) = TAU_MAJOR(2)+TAU_MAJOR1(2)
          case default
             print *, 'indco2 in taugb1a has an invalid value'
             stop
          end select


           FRACS(LAY,IG) = FRACREFA(IG,JPL) + FPL *
     &              (FRACREFA(IG,JPL+1)-FRACREFA(IG,JPL))

 2000          CONTINUE

c Complete calculation 

            DO 2040 IG = 1, NG(1)
               TAUSELF = SELFFAC(LAY)* (SELFREF(INDS,IG) +
     &              SELFFRAC(LAY)  *
     &              (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
               TAUFOR = FORFAC(LAY) * (FORREF(INDF,IG) +
     &              FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &              FORREF(INDF,IG))) 

               N2M1 = KA_MN2(JMN2,INDM,IG) + FMN2*
     &              (KA_MN2(JMN2+1,INDM,IG) -
     &              KA_MN2(JMN2,INDM,IG))
               N2M2 = KA_MN2(JMN2,INDM+1,IG) + FMN2*
     &              (KA_MN2(JMN2+1,INDM+1,IG) -
     &              KA_MN2(JMN2,INDM+1,IG))
               ABSN2 = N2M1 + MINORFRAC(LAY) *
     &              (N2M2 - N2M1)

               N2_O2_M1 = KA_MN2_O2(JMN2,INDM,IG) + FMN2*
     &              (KA_MN2_O2(JMN2+1,INDM,IG) -
     &              KA_MN2_O2(JMN2,INDM,IG))
               N2_O2_M2 = KA_MN2_O2(JMN2,INDM+1,IG) + FMN2*
     &              (KA_MN2_O2(JMN2+1,INDM+1,IG) -
     &              KA_MN2_O2(JMN2,INDM+1,IG))
               ABSN2_O2 = N2_O2_M1 + MINORFRAC(LAY) *
     &              (N2_O2_M2 - N2_O2_M1)


               ratio = (chi_std(7,jp(lay))-rat_o2_dry(lay))/
     &                      chi_std(7,jp(lay))
               tauminor = ratio*absn2+(1.0-ratio)*absn2_o2
               tauminor = tauminor*scaleminor(lay)*colbrd(lay)

C                taug(lay,ig) = taug(lay,ig) 
C      &              + TAUSELF + TAUFOR + tauminor !NJE

               taug(lay,ig) = taug(lay,ig) + tauminor
               if( h2o_sb == 1 ) then 
                  taug(lay,ig) = taug(lay,ig) + TAUSELF
                endif
               if( h2o_for == 1 ) then 
                  taug(lay,ig) = taug(lay,ig) + TAUFOR
                endif
                 

 2040       continue
 2500    continue

      RETURN
      END


      SUBROUTINE TAUGB1
      use variables !NJE

C     Written by Eli J. Mlawer, Atmospheric & Environmental Research.

C     BAND 1:  120-350 cm-1 (key - H2O, CO2; low minor - N2)

C     NOTE: Previous versions of RRTM BAND 1: 
C           10-250 cm-1 (low - H2O; high - H2O)


      PARAMETER (mg=16, mxlay=203, MXMOL = 38, NBANDS=16)
      PARAMETER (NTEMP=14,NETA=9,NREF=71,NCO2=2)

C  Output

      COMMON /TAUGCOM/  TAUG(MXLAY,MG)
      COMMON /PLANKG/   FRACS(MXLAY,MG)                                       

C  Input

      COMMON /FEATURES/ NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)
      COMMON /PRECISE/  ONEMINUS
      COMMON /PROFILE/  NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
     &                  PZ(0:MXLAY),TZ(0:MXLAY)
      COMMON /SPECIES/  COLDRY(MXLAY),WKL(MXMOL,MXLAY),WBROAD(MXLAY),
     &                  COLMOL(MXLAY),NMOL
      COMMON /PROFDATA/ LAYTROP, LAY_SWITCH_KABS,                                  
     &                  COLH2O(MXLAY),COLCO2(MXLAY),COLO3(MXLAY),  
     &                  COLN2O(MXLAY),COLCO(MXLAY),COLCH4(MXLAY),  
     &                  COLO2(MXLAY),COLBRD(MXLAY)
      COMMON /INTFAC/   FAC00(MXLAY),FAC01(MXLAY),                            
     &                  FAC10(MXLAY),FAC11(MXLAY)                             
      COMMON /CO2FAC/   co2_fac0(mxlay),co2_fac1(mxlay)
      COMMON /STD_REF/  PREF(nref),PREFLOG(nref),TREF(nref),
     &                  CHI_std(7,nref)
      COMMON /REFRAT_ETA/ RAT_H2OCO2(MXLAY),RAT_H2OCO2_1(MXLAY),
     &                  RAT_H2OO3(MXLAY),RAT_H2OO3_1(MXLAY),
     &                  RAT_H2ON2O(MXLAY),RAT_H2ON2O_1(MXLAY),
     &                  RAT_H2OCH4(MXLAY),RAT_H2OCH4_1(MXLAY),
     &                  RAT_N2OCO2(MXLAY),RAT_N2OCO2_1(MXLAY),
     &                  RAT_O3CO2(MXLAY),RAT_O3CO2_1(MXLAY),
     &                  RAT_CO2O3(MXLAY),RAT_CO2O3_1(MXLAY),
     &                  RAT_O2_DRY(MXLAY)

      COMMON /INTIND/   JP(MXLAY),JT(MXLAY),JT1(MXLAY)
      COMMON /SELF/     SELFFAC(MXLAY), SELFFRAC(MXLAY), INDSELF(MXLAY)
      COMMON /MINOR/    MINORFRAC(MXLAY), INDMINOR(MXLAY), 
     &                  SCALEMINOR(MXLAY),SCALEMINORN2(MXLAY)
      COMMON /FOREIGN/  FORFAC(MXLAY), FORFRAC(MXLAY), INDFOR(MXLAY)
      COMMON /CO2/      INDCO2(MXLAY)
      COMMON /K1/       KA(NETA,NTEMP,NREF,MG),  FORREF(4,MG),
     &                  SELFREF(10,MG), K_CO2(NETA,NTEMP,NREF,MG),
     &                  KA_MN2(NETA,19,MG),KA_MN2_O2(NETA,19,MG)

      COMMON /CVRTAU/    HNAMTAU,HVRTAU

      CHARACTER*18       HNAMTAU,HVRTAU

      DIMENSION ABSA(8946,MG),ABSB(8946,MG),ABSC(8946,MG)
      DIMENSION FRACREFA(MG,9)

      integer js(nco2),js1(nco2),ind0(nco2),ind1(nco2)
      real fs(nco2),fs1(nco2)

      real fac000(nco2),fac100(nco2),fac200(nco2)
      real fac010(nco2),fac110(nco2),fac210(nco2)
      real fac001(nco2),fac101(nco2),fac201(nco2)
      real fac011(nco2),fac111(nco2),fac211(nco2)

      real specparm(nco2),specparm1(nco2)
      real speccomb(nco2),speccomb1(nco2)


      real tau_base(mg)
      real tau_major(nco2),tau_major1(nco2)

      real n2m1, n2m2, n2_o2_m1, n2_o2_m2


      EQUIVALENCE (KA,ABSB)
      EQUIVALENCE (K_CO2,ABSC)
      REAL KA, KB, K_CO2, KA_MN2,  KA_MN2_O2, MINORFRAC

C     Compute the optical depth by interpolating in ln(pressure) and 
C     temperature.  Below LAYTROP, the water vapor self-continuum and
C     foreign continuum is interpolated (in temperature) separately.



C Planck fraction mapping level: P = 1013.0 mbar, T = 290.0 K


      DATA (FRACREFA(IG, 1),IG=1,16) /
     &2.0896E-01,1.9095E-01,1.6568E-01,1.3613E-01,1.0573E-01,7.7513E-02,
     &5.3466E-02,3.4285E-02,2.0343E-02,2.0745E-03,1.6955E-03,1.3270E-03,
     &9.6610E-04,6.1207E-04,2.3005E-04,3.2342E-05/
      DATA (FRACREFA(IG, 2),IG=1,16) /
     &2.0300E-01,1.7912E-01,1.5195E-01,1.2377E-01,9.5101E-02,8.0125E-02,
     &6.7254E-02,5.1514E-02,3.5461E-02,3.8678E-03,3.1653E-03,2.4120E-03,
     &1.7230E-03,1.0427E-03,4.2869E-04,5.9007E-05/
      DATA (FRACREFA(IG, 3),IG=1,16) /
     &2.0238E-01,1.7919E-01,1.4763E-01,1.1820E-01,9.7590E-02,8.6027E-02,
     &6.8842E-02,5.1919E-02,3.5516E-02,3.8688E-03,3.1659E-03,2.4123E-03,
     &1.7230E-03,1.0428E-03,4.2869E-04,5.9007E-05/
      DATA (FRACREFA(IG, 4),IG=1,16) /
     &2.0178E-01,1.7966E-01,1.4463E-01,1.1714E-01,9.9764E-02,8.7302E-02,
     &6.9449E-02,5.2040E-02,3.5535E-02,3.8696E-03,3.1662E-03,2.4126E-03,
     &1.7228E-03,1.0430E-03,4.2875E-04,5.9007E-05/
      DATA (FRACREFA(IG, 5),IG=1,16) /
     &2.0155E-01,1.7983E-01,1.4357E-01,1.1652E-01,1.0057E-01,8.7907E-02,
     &6.9709E-02,5.2099E-02,3.5544E-02,3.8693E-03,3.1665E-03,2.4126E-03,
     &1.7229E-03,1.0427E-03,4.2886E-04,5.8969E-05/
      DATA (FRACREFA(IG, 6),IG=1,16) /
     &2.0141E-01,1.7993E-01,1.4310E-01,1.1597E-01,1.0111E-01,8.8242E-02,
     &6.9848E-02,5.2134E-02,3.5551E-02,3.8700E-03,3.1668E-03,2.4122E-03,
     &1.7231E-03,1.0428E-03,4.2883E-04,5.9007E-05/
      DATA (FRACREFA(IG, 7),IG=1,16) /
     &2.0133E-01,1.7999E-01,1.4269E-01,1.1569E-01,1.0152E-01,8.8438E-02,
     &6.9936E-02,5.2156E-02,3.5555E-02,3.8702E-03,3.1661E-03,2.4126E-03,
     &1.7228E-03,1.0427E-03,4.2893E-04,5.9007E-05/
      DATA (FRACREFA(IG, 8),IG=1,16) /
     &2.0128E-01,1.8002E-01,1.4133E-01,1.1661E-01,1.0176E-01,8.8558E-02,
     &7.0003E-02,5.2171E-02,3.5556E-02,3.8694E-03,3.1663E-03,2.4124E-03,
     &1.7231E-03,1.0425E-03,4.2884E-04,5.9081E-05/
      DATA (FRACREFA(IG, 9),IG=1,16) /
     &2.0155E-01,1.7983E-01,1.4357E-01,1.1651E-01,1.0058E-01,8.7909E-02,
     &6.9710E-02,5.2099E-02,3.5545E-02,3.8691E-03,3.1664E-03,2.4128E-03,
     &1.7228E-03,1.0429E-03,4.2885E-04,5.9007E-05/

      co2_lo_hi_rat = 0.04*3.85e-4/0.95

      HVRTAU = '$Revision: 22602 $'

      idump = 1


C     Calculate reference ratio to be used in calculation of Planck
C     fraction 


      base_rat = 15./9999.
      iplanck_lay = 21

      REFRAT_PLANCK_A = CHI_STD(1,iplanck_lay)/CHI_STD(2,iplanck_lay)
      if (pref(iplanck_lay).le. 120.) then
	 rnew_rat = base_rat
      else if (pref(iplanck_lay).le.800) then 
	 rnew_rat = base_rat*(pref(iplanck_lay)/120.)**(-3.7)
      else
	 rnew_rat= base_rat*(800./120.)**(-3.7)
      endif
      refrat_planck_a = refrat_planck_a*rnew_rat

C P = 142.549 mb
      REFRAT_M_A = CHI_STD(1,23)/CHI_STD(2,23)

C     Compute the optical depth by interpolating in ln(pressure), 
C     temperature, and appropriate species.  Below LAYTROP, the water
C     vapor self-continuum and foreign continuum is interpolated 
C     (in temperature) separately.  

      HVRTAU = '$Revision: 22602 $'
      base_rat = 15./9999.

      DO 2500 LAY = 1, NLAYERS
         if (pavel(lay).le. 120.) then
            rnew_rat = base_rat
         else if (pavel(lay).le.800) then 
            rnew_rat = base_rat*(pavel(lay)/120.)**(-3.7)
         else
            rnew_rat= base_rat*(800./120.)**(-3.7)
         endif

c Set up calculation for low and/or high CO2 values
         select case (indco2(lay))
         case (0)
             ico2_0 = 1
             ico2_1 = 1
         case(2)
             ico2_0 = 2
             ico2_1 = 2
         case(1)
             ico2_0 = 1
             ico2_1 = 2
          case default
             print *, 'indco2 in taugb1 has an invalid value'
             stop
         end select

        do ic=ico2_0,ico2_1
         if (ic.eq.1) then
            scale_eta=rnew_rat 
         else 
            scale_eta= co2_lo_hi_rat
         end if

         SPECCOMB(ic) = COLH2O(LAY) + scale_eta*
     &           RAT_H2OCO2(LAY)*COLCO2(LAY)
         SPECPARM(ic) = COLH2O(LAY)/SPECCOMB(ic)
         IF (SPECPARM(ic) .GE. ONEMINUS) SPECPARM(ic) = ONEMINUS
         SPECMULT = 8.*(SPECPARM(ic))
         JS(ic) = 1 + INT(SPECMULT)
         FS(ic) = AMOD(SPECMULT,1.0)

         SPECCOMB1(ic) = COLH2O(LAY)+scale_eta
     &        *RAT_H2OCO2_1(LAY)*COLCO2(LAY)
         SPECPARM1(ic) = COLH2O(LAY)/SPECCOMB1(ic)
         IF (SPECPARM1(ic) .GE. ONEMINUS) SPECPARM1(ic) = ONEMINUS
         SPECMULT1 = 8.*(SPECPARM1(ic))
         JS1(ic) = 1 + INT(SPECMULT1)
         FS1(ic) = AMOD(SPECMULT1,1.0)

         SPECCOMB_PLANCK = COLH2O(LAY)+
     &         REFRAT_PLANCK_A*COLCO2(LAY)
     &                *(pavel(lay)/pref(iplanck_lay))
         SPECPARM_PLANCK = COLH2O(LAY)/SPECCOMB_PLANCK
         IF (SPECPARM_PLANCK .GE. ONEMINUS) SPECPARM_PLANCK=ONEMINUS
         SPECMULT_PLANCK = 8.*SPECPARM_PLANCK
         JPL= 1 + INT(SPECMULT_PLANCK)
         FPL = AMOD(SPECMULT_PLANCK,1.0)

         SPECCOMB_MN2 = COLH2O(LAY) + REFRAT_M_A*COLCO2(LAY)
         SPECPARM_MN2 = COLH2O(LAY)/SPECCOMB_MN2
         IF (SPECPARM_MN2 .GE. ONEMINUS) SPECPARM_MN2 = ONEMINUS
         SPECMULT_MN2 = 8.*SPECPARM_MN2
         JMN2 = 1 + INT(SPECMULT_MN2)
         FMN2 = AMOD(SPECMULT_MN2,1.0)
         FMN2MF = MINORFRAC(LAY)*FMN2

         IND0(ic) = ((JP(LAY)-1)*NTEMP+(JT(LAY)-1))*NSPA(1) + JS(ic)
         IND1(ic) = (JP(LAY)*NTEMP+(JT1(LAY)-1))*NSPA(1) + JS1(ic)
         INDS = INDSELF(LAY)
         INDF = INDFOR(LAY)
         INDM = INDMINOR(LAY)

         IF (SPECPARM(ic) .LT. 0.125) THEN
            P = FS(ic) - 1
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC000(ic) = FK0*FAC00(LAY)
            FAC100(ic) = FK1*FAC00(LAY)
            FAC200(ic) = FK2*FAC00(LAY)
            FAC010(ic) = FK0*FAC10(LAY)
            FAC110(ic) = FK1*FAC10(LAY)
            FAC210(ic) = FK2*FAC10(LAY)
         ELSE IF (SPECPARM(ic) .GT. 0.875) THEN
            P = -FS(ic)
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC000(ic) = FK0*FAC00(LAY)
            FAC100(ic) = FK1*FAC00(LAY)
            FAC200(ic) = FK2*FAC00(LAY)
            FAC010(ic) = FK0*FAC10(LAY)
            FAC110(ic) = FK1*FAC10(LAY)
            FAC210(ic) = FK2*FAC10(LAY)
	     corradj = 1.
	     wv_log = alog10(wkl(1,lay)/coldry(lay))
	     if(wv_log.gt.-5.AND.wv_log.le.-3.0) then
	       corradj = 1.+0.05*(wv_log+3.0)
	     else if (wv_log.gt.-5.5.AND.wv_log.le.-5.) then
	       corradj = 0.9
	     else if (wv_log.gt.-6.5.AND.wv_log.le.-5.5) then
	       corradj = 1.-0.1*(wv_log+6.5)
	     endif
	     corradj = 1.
         ELSE
            FAC000(ic) = (1. - FS(ic)) * FAC00(LAY)
            FAC010(ic) = (1. - FS(ic)) * FAC10(LAY)
            FAC100(ic) = FS(ic) * FAC00(LAY)
            FAC110(ic) = FS(ic) * FAC10(LAY)
         END IF

       
         IF (SPECPARM1(ic) .LT. 0.125 ) THEN
            P = FS1(ic) - 1
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC001(ic) = FK0*FAC01(LAY)
            FAC101(ic) = FK1*FAC01(LAY)
            FAC201(ic) = FK2*FAC01(LAY)
            FAC011(ic) = FK0*FAC11(LAY)
            FAC111(ic) = FK1*FAC11(LAY)
            FAC211(ic) = FK2*FAC11(LAY)
         ELSE IF (SPECPARM1(ic) .GT. 0.875 ) THEN
            P = -FS1(ic)
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC001(ic) = FK0*FAC01(LAY)
            FAC101(ic) = FK1*FAC01(LAY)
            FAC201(ic) = FK2*FAC01(LAY)
            FAC011(ic) = FK0*FAC11(LAY)
            FAC111(ic) = FK1*FAC11(LAY)
            FAC211(ic) = FK2*FAC11(LAY)
              corradj1 = 1.
              wv_log = alog10(wkl(1,lay)/coldry(lay))
	     if(wv_log.gt.-5.AND.wv_log.le.-3.0) then
	       corradj1 = 1.+0.05*(wv_log+3.0)
	     else if (wv_log.gt.-5.5.AND.wv_log.le.-5.) then
	       corradj1 = 0.9
	     else if (wv_log.gt.-6.5.AND.wv_log.le.-5.5) then
	       corradj1 = 1.-0.1*(wv_log+6.5)
	     endif
              corradj1 = 1.
         ELSE
            FAC001(ic) = (1. - FS1(ic)) * FAC01(LAY)
            FAC011(ic) = (1. - FS1(ic)) * FAC11(LAY)
            FAC101(ic) = FS1(ic) * FAC01(LAY)
            FAC111(ic) = FS1(ic) * FAC11(LAY)
         END IF

        end do   ! end ic loop for general layer properties

    
         DO 2000 IG = 1, NG(1)
	   DO IC=ico2_0,ico2_1
	    if(ic.eq.1) then
	       absa = absb
	    else
		absa = absc
	    endif 

            IF (SPECPARM(ic).LT.0.125) THEN 
               TAU_MAJOR(ic) = SPECCOMB(ic) *
     &          (FAC000(ic)* ABSA(ind0(ic),IG) +
     &          FAC100(ic)* ABSA(ind0(ic)+1,IG) +
     &          FAC200(ic)* ABSA(ind0(ic)+2,IG) +
     &          FAC010(ic)* ABSA(ind0(ic)+9,IG) +
     &          FAC110(ic)* ABSA(ind0(ic)+10,IG) +
     &          FAC210(ic)* ABSA(ind0(ic)+11,IG))
           ELSE IF (SPECPARM(ic).GT.0.875) THEN 
                TAU_MAJOR(ic) = CORRADJ*speccomb(ic) * 
     &              (FAC200(ic)* ABSA(ind0(ic)-1,IG) +
     &              FAC100(ic)* ABSA(ind0(ic),IG) +
     &              FAC000(ic)* ABSA(ind0(ic)+1,IG) +
     &              FAC210(ic)* ABSA(ind0(ic)+8,IG) +
     &              FAC110(ic)* ABSA(ind0(ic)+9,IG) +
     &              FAC010(ic)* ABSA(ind0(ic)+10,IG))
           ELSE 
               TAU_MAJOR(ic) = speccomb(ic) * 
     &              (FAC000(ic) * ABSA(IND0(ic),IG) +
     &              FAC100(ic) * ABSA(IND0(ic)+1,IG) +
     &              FAC010(ic) * ABSA(IND0(ic)+9,IG) +
     &              FAC110(ic) * ABSA(IND0(ic)+10,IG))
           END IF

           IF (SPECPARM1(ic).LT.0.125) THEN 
                TAU_MAJOR1(ic) =  speccomb1(ic) *
     &              (FAC001(ic) * ABSA(IND1(ic),IG) + 
     &              FAC101(ic) * ABSA(IND1(ic)+1,IG) +
     &              FAC201(ic) * ABSA(IND1(ic)+2,IG) +
     &              FAC011(ic) * ABSA(IND1(ic)+9,IG) +
     &              FAC111(ic) * ABSA(IND1(ic)+10,IG) +
     &              FAC211(ic) * ABSA(IND1(ic)+11,IG))

           ELSE IF (SPECPARM1(ic).GT.0.875) THEN 
                TAU_MAJOR1(ic) = CORRADJ1*speccomb1(ic) * 
     &              (FAC201(ic) * ABSA(IND1(ic)-1,IG) +
     &              FAC101(ic) * ABSA(IND1(ic),IG) +
     &              FAC001(ic) * ABSA(IND1(ic)+1,IG) +
     &              FAC211(ic) * ABSA(IND1(ic)+8,IG) +
     &              FAC111(ic) * ABSA(IND1(ic)+9,IG) +
     &              FAC011(ic) * ABSA(IND1(ic)+10,IG))
           ELSE
               TAU_MAJOR1(ic) = speccomb1(ic) *
     &              (FAC001(ic) * ABSA(IND1(ic),IG) +
     &              FAC101(ic) * ABSA(IND1(ic)+1,IG) +
     &              FAC011(ic) * ABSA(IND1(ic)+9,IG) +
     &              FAC111(ic) * ABSA(IND1(ic)+10,IG))

           END IF
          end do    ! co2 loop

C  Do  interpolation for high CO2 cases if needed
          select case (indco2(lay) )
          case (0)
             TAUG(LAY,IG) = TAU_MAJOR(1)+TAU_MAJOR1(1)
          case (1)
             TAUG(LAY,IG) = co2_fac0(lay)*(TAU_MAJOR(1)+TAU_MAJOR1(1)) + 
     &                      co2_fac1(lay)*(TAU_MAJOR(2)+TAU_MAJOR1(2)) 
          case (2)
             TAUG(LAY,IG) = TAU_MAJOR(2)+TAU_MAJOR1(2)
          case default
             print *, 'indco2 in taugb1 has an invalid value'
             stop
          end select


           FRACS(LAY,IG) = FRACREFA(IG,JPL) + FPL *
     &              (FRACREFA(IG,JPL+1)-FRACREFA(IG,JPL))

 2000          CONTINUE

c Complete calculation 

            DO 2040 IG = 1, NG(2)
               TAUSELF = SELFFAC(LAY)* (SELFREF(INDS,IG) +
     &              SELFFRAC(LAY)  *
     &              (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
               TAUFOR = FORFAC(LAY) * (FORREF(INDF,IG) +
     &              FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &              FORREF(INDF,IG))) 

               N2M1 = KA_MN2(JMN2,INDM,IG) + FMN2*
     &              (KA_MN2(JMN2+1,INDM,IG) -
     &              KA_MN2(JMN2,INDM,IG))
               N2M2 = KA_MN2(JMN2,INDM+1,IG) + FMN2*
     &              (KA_MN2(JMN2+1,INDM+1,IG) -
     &              KA_MN2(JMN2,INDM+1,IG))
               ABSN2 = N2M1 + MINORFRAC(LAY) *
     &              (N2M2 - N2M1)

               N2_O2_M1 = KA_MN2_O2(JMN2,INDM,IG) + FMN2*
     &              (KA_MN2_O2(JMN2+1,INDM,IG) -
     &              KA_MN2_O2(JMN2,INDM,IG))
               N2_O2_M2 = KA_MN2_O2(JMN2,INDM+1,IG) + FMN2*
     &              (KA_MN2_O2(JMN2+1,INDM+1,IG) -
     &              KA_MN2_O2(JMN2,INDM+1,IG))
               ABSN2_O2 = N2_O2_M1 + MINORFRAC(LAY) *
     &              (N2_O2_M2 - N2_O2_M1)


               ratio = (chi_std(7,jp(lay))-rat_o2_dry(lay))/
     &                      chi_std(7,jp(lay))
               tauminor = ratio*absn2+(1.0-ratio)*absn2_o2
               tauminor = tauminor*scaleminor(lay)*colbrd(lay)

C                taug(lay,ig) = taug(lay,ig) 
C      &              + TAUSELF + TAUFOR + tauminor

                taug(lay,ig) = taug(lay,ig) + tauminor
               if( h2o_sb == 1 ) then 
                  taug(lay,ig) = taug(lay,ig) + TAUSELF
                endif
               if( h2o_for == 1 ) then 
                  taug(lay,ig) = taug(lay,ig) + TAUFOR
                endif
                 

 2040       continue
 2500    continue

      RETURN
      END

C----------------------------------------------------------------------------


C----------------------------------------------------------------------------

      SUBROUTINE TAUGB2
      use variables !NJE

C     BAND 2:  350-500 cm-1 (key - H2O, CO2)

C     NOTE: Previous version of RRTM BAND 2: 
C           250 - 500 cm-1 (low - H2O; high - H2O)

      PARAMETER (mg=16, mxlay=203, MXMOL = 38, NBANDS=16)
      PARAMETER (NTEMP=14,NETA=9,NREF=71,NCO2=2)

C  Output

      COMMON /TAUGCOM/  TAUG(MXLAY,MG)
      COMMON /PLANKG/   FRACS(MXLAY,MG)                                       

C  Input

      COMMON /FEATURES/ NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)
      COMMON /PRECISE/  ONEMINUS
      COMMON /PROFILE/  NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
     &                  PZ(0:MXLAY),TZ(0:MXLAY)
      COMMON /SPECIES/  COLDRY(MXLAY),WKL(MXMOL,MXLAY),WBROAD(MXLAY),
     &                  COLMOL(MXLAY),NMOL
      COMMON /PROFDATA/ LAYTROP, LAY_SWITCH_KABS,                                  
     &                  COLH2O(MXLAY),COLCO2(MXLAY),COLO3(MXLAY),  
     &                  COLN2O(MXLAY),COLCO(MXLAY),COLCH4(MXLAY),  
     &                  COLO2(MXLAY),COLBRD(MXLAY)
      COMMON /INTFAC/   FAC00(MXLAY),FAC01(MXLAY),                            
     &                  FAC10(MXLAY),FAC11(MXLAY)                             
      COMMON /CO2FAC/   co2_fac0(mxlay), co2_fac1(mxlay)
      COMMON /STD_REF/  PREF(nref),PREFLOG(nref),TREF(nref),
     &                  CHI_std(7,nref)
      COMMON /REFRAT_ETA/ RAT_H2OCO2(MXLAY),RAT_H2OCO2_1(MXLAY),
     &                  RAT_H2OO3(MXLAY),RAT_H2OO3_1(MXLAY),
     &                  RAT_H2ON2O(MXLAY),RAT_H2ON2O_1(MXLAY),
     &                  RAT_H2OCH4(MXLAY),RAT_H2OCH4_1(MXLAY),
     &                  RAT_N2OCO2(MXLAY),RAT_N2OCO2_1(MXLAY),
     &                  RAT_O3CO2(MXLAY),RAT_O3CO2_1(MXLAY),
     &                  RAT_CO2O3(MXLAY),RAT_CO2O3_1(MXLAY),
     &                  RAT_O2_DRY(MXLAY)

      COMMON /INTIND/   JP(MXLAY),JT(MXLAY),JT1(MXLAY)
      COMMON /SELF/     SELFFAC(MXLAY), SELFFRAC(MXLAY), INDSELF(MXLAY)
      COMMON /FOREIGN/  FORFAC(MXLAY), FORFRAC(MXLAY), INDFOR(MXLAY)
      COMMON /CO2/      INDCO2(MXLAY)
      COMMON /K2/       KA(NETA,NTEMP,NREF,MG),  FORREF(4,MG),
     &                  SELFREF(10,MG), K_CO2(NETA,NTEMP,NREF,MG)

      COMMON /CVRTAU/    HNAMTAU,HVRTAU

      CHARACTER*18       HNAMTAU,HVRTAU

      DIMENSION ABSA(8946,MG)
      DIMENSION FRACREFA(MG,9),ABSB(8946,MG),ABSC(8946,MG)

      real tau_hi_co2(2,2)

      real KA, K_CO2

      integer js(nco2),js1(nco2),ind0(nco2),ind1(nco2)
      real fs(nco2),fs1(nco2)

      real fac000(nco2),fac100(nco2),fac200(nco2)
      real fac010(nco2),fac110(nco2),fac210(nco2)
      real fac001(nco2),fac101(nco2),fac201(nco2)
      real fac011(nco2),fac111(nco2),fac211(nco2)

      real specparm(nco2),specparm1(nco2)
      real speccomb(nco2),speccomb1(nco2)


      real tau_base(mg)
      real tau_major(nco2),tau_major1(nco2)
      real corradj(nco2),corradj1(nco2)

      EQUIVALENCE (KA,ABSB)
      EQUIVALENCE (K_CO2,ABSC)

C Planck fraction mapping level: P = 1053.630 mbar, T = 294.2 K

      DATA (FRACREFA(IG, 1),IG=1,16) /
     &1.3986E-01,1.4172E-01,1.4023E-01,1.3428E-01,1.2291E-01,1.0809E-01,
     &8.9240E-02,6.7334E-02,4.2578E-02,4.0592E-03,3.3529E-03,2.6470E-03,
     &1.9400E-03,1.2348E-03,4.6529E-04,6.5464E-05/
      DATA (FRACREFA(IG, 2),IG=1,16) /
     &1.6255E-01,1.5264E-01,1.4167E-01,1.3005E-01,1.1659E-01,1.0052E-01,
     &8.0070E-02,6.0686E-02,4.0624E-02,4.3839E-03,3.5997E-03,2.7640E-03,
     &2.0180E-03,1.2898E-03,4.7761E-04,6.9319E-05/
      DATA (FRACREFA(IG, 3),IG=1,16) /
     &1.6326E-01,1.5222E-01,1.4197E-01,1.2963E-01,1.1646E-01,1.0051E-01,
     &8.0060E-02,6.0675E-02,4.0618E-02,4.3825E-03,3.5998E-03,2.7619E-03,
     &2.0169E-03,1.2904E-03,4.7764E-04,6.9319E-05/
      DATA (FRACREFA(IG, 4),IG=1,16) /
     &1.6348E-01,1.5228E-01,1.4234E-01,1.2923E-01,1.1629E-01,1.0047E-01,
     &8.0066E-02,6.0657E-02,4.0608E-02,4.3804E-03,3.6001E-03,2.7586E-03,
     &2.0153E-03,1.2911E-03,4.7772E-04,6.9319E-05/
      DATA (FRACREFA(IG, 5),IG=1,16) /
     &1.6364E-01,1.5256E-01,1.4310E-01,1.2839E-01,1.1600E-01,1.0042E-01,
     &8.0076E-02,6.0631E-02,4.0588E-02,4.3788E-03,3.6025E-03,2.7509E-03,
     &2.0132E-03,1.2918E-03,4.7778E-04,6.9319E-05/
      DATA (FRACREFA(IG, 6),IG=1,16) /
     &1.6381E-01,1.5279E-01,1.4384E-01,1.2784E-01,1.1558E-01,1.0031E-01,
     &8.0136E-02,6.0548E-02,4.0577E-02,4.3662E-03,3.5962E-03,2.7494E-03,
     &2.0110E-03,1.2914E-03,4.7794E-04,6.9319E-05/
      DATA (FRACREFA(IG, 7),IG=1,16) /
     &1.6398E-01,1.5358E-01,1.4359E-01,1.2793E-01,1.1510E-01,1.0011E-01,
     &8.0254E-02,6.0393E-02,4.0540E-02,4.3509E-03,3.5819E-03,2.7521E-03,
     &2.0056E-03,1.2905E-03,4.7823E-04,6.9319E-05/
      DATA (FRACREFA(IG, 8),IG=1,16) /
     &1.6419E-01,1.5551E-01,1.4335E-01,1.2761E-01,1.1427E-01,9.9810E-02,
     &8.0223E-02,6.0239E-02,4.0322E-02,4.3468E-03,3.5517E-03,2.7638E-03,
     &1.9823E-03,1.2880E-03,4.7903E-04,6.9319E-05/
      DATA (FRACREFA(IG, 9),IG=1,16) /
     &1.6365E-01,1.5256E-01,1.4310E-01,1.2839E-01,1.1600E-01,1.0042E-01,
     &8.0076E-02,6.0632E-02,4.0588E-02,4.3786E-03,3.6024E-03,2.7512E-03,
     &2.0132E-03,1.2918E-03,4.7778E-04,6.9319E-05/

      co2_lo_hi_rat = 0.04*3.85e-4/0.95


      idump = 1

C     Calculate reference ratio to be used in calculation of Planck
C     fraction 

C     P = 1053 mb
      REFRAT_PLANCK_A = CHI_STD(1,13)/CHI_STD(2,13)

C     Compute the optical depth by interpolating in ln(pressure), 
C     temperature, and appropriate species.  Below LAYTROP, the water
C     vapor self-continuum and foreign continuum is interpolated 
C     (in temperature) separately.  

      HVRTAU = '$Revision: 22602 $'

      DO 2500 LAY = 1, NLAYERS

C  Set corradj values to empirically correct for underestimation of OD in high co2 case
      corradj(1) = 1.0
      corradj1(1) = 1.0
      corradj(2) = 1.06
      corradj1(2) = 1.06
c     if (tavel(lay).ge.200) then
c         corradj(2) = (tavel(lay)-200.)*0.00625+0.95
c         corradj1(2) = (tavel(lay)-200.)*0.00625+0.95
c      endif

c Set up calculation for low and/or high CO2 values
      select case (indco2(lay))
      case (0)
	  ico2_0 = 1
	  ico2_1 = 1
      case(2)
	  ico2_0 = 2
	  ico2_1 = 2
      case(1)
	  ico2_0 = 1
	  ico2_1 = 2
       case default
	  print *, 'indco2 in taugb2 has an invalid value'
	  stop
      end select

      do ic=ico2_0,ico2_1
         if (ic.eq.1) then
            scale_eta=1.0
         else
            scale_eta= co2_lo_hi_rat
         end if

         SPECCOMB(ic) = COLH2O(LAY) + scale_eta*
     &           RAT_H2OCO2(LAY)*COLCO2(LAY)
         SPECPARM(ic) = COLH2O(LAY)/SPECCOMB(ic)
         IF (SPECPARM(ic) .GE. ONEMINUS) SPECPARM(ic) = ONEMINUS
         SPECMULT = 8.*(SPECPARM(ic))
         JS(ic) = 1 + INT(SPECMULT)
         FS(ic) = AMOD(SPECMULT,1.0)

         SPECCOMB1(ic) = COLH2O(LAY)+scale_eta
     &        *RAT_H2OCO2_1(LAY)*COLCO2(LAY)
         SPECPARM1(ic) = COLH2O(LAY)/SPECCOMB1(ic)
         IF (SPECPARM1(ic) .GE. ONEMINUS) SPECPARM1(ic) = ONEMINUS
         SPECMULT1 = 8.*(SPECPARM1(ic))
         JS1(ic) = 1 + INT(SPECMULT1)
         FS1(ic) = AMOD(SPECMULT1,1.0)

         SPECCOMB_PLANCK = COLH2O(LAY)+
     &         REFRAT_PLANCK_A*COLCO2(LAY)
         SPECPARM_PLANCK = COLH2O(LAY)/SPECCOMB_PLANCK
         IF (SPECPARM_PLANCK .GE. ONEMINUS) SPECPARM_PLANCK=ONEMINUS
         SPECMULT_PLANCK = 8.*SPECPARM_PLANCK
         JPL= 1 + INT(SPECMULT_PLANCK)
         FPL = AMOD(SPECMULT_PLANCK,1.0)

         IND0(ic) = ((JP(LAY)-1)*NTEMP+(JT(LAY)-1))*NSPA(2) + JS(ic)
         IND1(ic) = (JP(LAY)*NTEMP+(JT1(LAY)-1))*NSPA(2) + JS1(ic)
         INDS = INDSELF(LAY)
         INDF = INDFOR(LAY)

         IF (SPECPARM(ic) .LT. 0.125) THEN
            P = FS(ic) - 1
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC000(ic) = FK0*FAC00(LAY)
            FAC100(ic) = FK1*FAC00(LAY)
            FAC200(ic) = FK2*FAC00(LAY)
            FAC010(ic) = FK0*FAC10(LAY)
            FAC110(ic) = FK1*FAC10(LAY)
            FAC210(ic) = FK2*FAC10(LAY)
         ELSE IF (SPECPARM(ic) .GT. 0.875) THEN
            P = -FS(ic)
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC000(ic) = FK0*FAC00(LAY)
            FAC100(ic) = FK1*FAC00(LAY)
            FAC200(ic) = FK2*FAC00(LAY)
            FAC010(ic) = FK0*FAC10(LAY)
            FAC110(ic) = FK1*FAC10(LAY)
            FAC210(ic) = FK2*FAC10(LAY)
         ELSE
            FAC000(ic) = (1. - FS(ic)) * FAC00(LAY)
            FAC010(ic) = (1. - FS(ic)) * FAC10(LAY)
            FAC100(ic) = FS(ic) * FAC00(LAY)
            FAC110(ic) = FS(ic) * FAC10(LAY)
         END IF

       
         IF (SPECPARM1(ic) .LT. 0.125 ) THEN
            P = FS1(ic) - 1
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC001(ic) = FK0*FAC01(LAY)
            FAC101(ic) = FK1*FAC01(LAY)
            FAC201(ic) = FK2*FAC01(LAY)
            FAC011(ic) = FK0*FAC11(LAY)
            FAC111(ic) = FK1*FAC11(LAY)
            FAC211(ic) = FK2*FAC11(LAY)
         ELSE IF (SPECPARM1(ic) .GT. 0.875 ) THEN
            P = -FS1(ic)
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC001(ic) = FK0*FAC01(LAY)
            FAC101(ic) = FK1*FAC01(LAY)
            FAC201(ic) = FK2*FAC01(LAY)
            FAC011(ic) = FK0*FAC11(LAY)
            FAC111(ic) = FK1*FAC11(LAY)
            FAC211(ic) = FK2*FAC11(LAY)
        ELSE
            FAC001(ic) = (1. - FS1(ic)) * FAC01(LAY)
            FAC011(ic) = (1. - FS1(ic)) * FAC11(LAY)
            FAC101(ic) = FS1(ic) * FAC01(LAY)
            FAC111(ic) = FS1(ic) * FAC11(LAY)
         END IF

        end do   ! end ic loop for general layer properties

    
         DO 2000 IG = 1, NG(2)
	   DO IC=ico2_0,ico2_1
	    if(ic.eq.1) then
	       absa = absb
	    else
		absa = absc
	    endif 

            IF (SPECPARM(ic).LT.0.125) THEN 
               TAU_MAJOR(ic) = SPECCOMB(ic) *
     &          (FAC000(ic)* ABSA(ind0(ic),IG) +
     &          FAC100(ic)* ABSA(ind0(ic)+1,IG) +
     &          FAC200(ic)* ABSA(ind0(ic)+2,IG) +
     &          FAC010(ic)* ABSA(ind0(ic)+9,IG) +
     &          FAC110(ic)* ABSA(ind0(ic)+10,IG) +
     &          FAC210(ic)* ABSA(ind0(ic)+11,IG))
           ELSE IF (SPECPARM(ic).GT.0.875) THEN 
                TAU_MAJOR(ic) = speccomb(ic) * 
     &              (FAC200(ic)* ABSA(ind0(ic)-1,IG) +
     &              FAC100(ic)* ABSA(ind0(ic),IG) +
     &              FAC000(ic)* ABSA(ind0(ic)+1,IG) +
     &              FAC210(ic)* ABSA(ind0(ic)+8,IG) +
     &              FAC110(ic)* ABSA(ind0(ic)+9,IG) +
     &              FAC010(ic)* ABSA(ind0(ic)+10,IG))
           ELSE 
               TAU_MAJOR(ic) = speccomb(ic) * 
     &              (FAC000(ic) * ABSA(IND0(ic),IG) +
     &              FAC100(ic) * ABSA(IND0(ic)+1,IG) +
     &              FAC010(ic) * ABSA(IND0(ic)+9,IG) +
     &              FAC110(ic) * ABSA(IND0(ic)+10,IG))
           END IF

           IF (SPECPARM1(ic).LT.0.125) THEN 
                TAU_MAJOR1(ic) =  speccomb1(ic) *
     &              (FAC001(ic) * ABSA(IND1(ic),IG) + 
     &              FAC101(ic) * ABSA(IND1(ic)+1,IG) +
     &              FAC201(ic) * ABSA(IND1(ic)+2,IG) +
     &              FAC011(ic) * ABSA(IND1(ic)+9,IG) +
     &              FAC111(ic) * ABSA(IND1(ic)+10,IG) +
     &              FAC211(ic) * ABSA(IND1(ic)+11,IG))

           ELSE IF (SPECPARM1(ic).GT.0.875) THEN 
                TAU_MAJOR1(ic) = speccomb1(ic) * 
     &              (FAC201(ic) * ABSA(IND1(ic)-1,IG) +
     &              FAC101(ic) * ABSA(IND1(ic),IG) +
     &              FAC001(ic) * ABSA(IND1(ic)+1,IG) +
     &              FAC211(ic) * ABSA(IND1(ic)+8,IG) +
     &              FAC111(ic) * ABSA(IND1(ic)+9,IG) +
     &              FAC011(ic) * ABSA(IND1(ic)+10,IG))
           ELSE
               TAU_MAJOR1(ic) = speccomb1(ic) *
     &              (FAC001(ic) * ABSA(IND1(ic),IG) +
     &              FAC101(ic) * ABSA(IND1(ic)+1,IG) +
     &              FAC011(ic) * ABSA(IND1(ic)+9,IG) +
     &              FAC111(ic) * ABSA(IND1(ic)+10,IG))

           END IF

C  Adjust CO2 optical depth 
         tau_major(ic) = tau_major(ic)*corradj(ic)
         tau_major1(ic) = tau_major1(ic)*corradj1(ic)

          end do    ! co2 loop


C  Do  interpolation for high CO2 cases if needed
          
          select case (indco2(lay) )
          case (0)
             TAUG(LAY,IG) = TAU_MAJOR(1)+TAU_MAJOR1(1)
          case (1)
             TAUG(LAY,IG) = co2_fac0(lay)*(TAU_MAJOR(1)+TAU_MAJOR1(1)) + 
     &                      co2_fac1(lay)*(TAU_MAJOR(2)+TAU_MAJOR1(2)) 
          case (2)
             TAUG(LAY,IG) = TAU_MAJOR(2)+TAU_MAJOR1(2)
          case default
             print *, 'indco2 in taugb2 has an invalid value'
             stop
          end select


           FRACS(LAY,IG) = FRACREFA(IG,JPL) + FPL *
     &              (FRACREFA(IG,JPL+1)-FRACREFA(IG,JPL))

 2000          CONTINUE

c Complete calculation 

            DO 2040 IG = 1, NG(2)
               TAUSELF = SELFFAC(LAY)* (SELFREF(INDS,IG) +
     &              SELFFRAC(LAY)  *
     &              (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
               TAUFOR = (1.0*co2_fac0(lay)+1.8*co2_fac1(lay))*
     &              FORFAC(LAY) * (FORREF(INDF,IG) +
     &              FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &              FORREF(INDF,IG))) 

C                    taug(lay,ig) = taug(lay,ig) 
C      &              + TAUSELF + TAUFOR 


                  taug(lay,ig) = taug(lay,ig)
               if( h2o_sb == 1 ) then 
                  taug(lay,ig) = taug(lay,ig) + TAUSELF
                endif
               if( h2o_for == 1 ) then 
                  taug(lay,ig) = taug(lay,ig) + TAUFOR
                endif


 2040    continue

 2500 CONTINUE

      RETURN
      END

C----------------------------------------------------------------------------

      SUBROUTINE TAUGB3
      use variables !NJE

C     BAND 3:  500-630 cm-1 (low key - H2O,CO2; low minor - n2o)
C                           (high key - H2O,CO2; high minor - n2o)

      PARAMETER (mg=16, mxlay=203, MXMOL=38, NBANDS=16)
      PARAMETER (NTEMP=14,NETA=9,NREF=71,NCO2=2)

C  Output

      COMMON /TAUGCOM/  TAUG(MXLAY,MG)
      COMMON /PLANKG/   FRACS(MXLAY,MG)                                       

C  Input

      COMMON /FEATURES/ NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)
      COMMON /PRECISE/  ONEMINUS
      COMMON /PROFILE/  NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
     &                  PZ(0:MXLAY),TZ(0:MXLAY)
      COMMON /PROFDATA/ LAYTROP, LAY_SWITCH_KABS,                                  
     &                  COLH2O(MXLAY),COLCO2(MXLAY),COLO3(MXLAY),  
     &                  COLN2O(MXLAY),COLCO(MXLAY),COLCH4(MXLAY),  
     &                  COLO2(MXLAY),COLBRD(MXLAY)
      COMMON /SPECIES/  COLDRY(MXLAY),WKL(MXMOL,MXLAY),WBROAD(MXLAY),
     &                  COLMOL(MXLAY),NMOL
      COMMON /INTFAC/   FAC00(MXLAY),FAC01(MXLAY),                            
     &                  FAC10(MXLAY),FAC11(MXLAY)  
      COMMON /CO2FAC/   co2_fac0(mxlay), co2_fac1(mxlay)
      COMMON /STD_REF/  PREF(NREF),PREFLOG(NREF),TREF(NREF),
     &                  CHI_STD(7,NREF)            
      COMMON /INTIND/   JP(MXLAY),JT(MXLAY),JT1(MXLAY)
      COMMON /REFRAT_ETA/ RAT_H2OCO2(MXLAY),RAT_H2OCO2_1(MXLAY),
     &                  RAT_H2OO3(MXLAY),RAT_H2OO3_1(MXLAY),
     &                  RAT_H2ON2O(MXLAY),RAT_H2ON2O_1(MXLAY),
     &                  RAT_H2OCH4(MXLAY),RAT_H2OCH4_1(MXLAY),
     &                  RAT_N2OCO2(MXLAY),RAT_N2OCO2_1(MXLAY),
     &                  RAT_O3CO2(MXLAY),RAT_O3CO2_1(MXLAY),
     &                  RAT_CO2O3(MXLAY),RAT_CO2O3_1(MXLAY),
     &                  RAT_O2_DRY(MXLAY)

      COMMON /SELF/     SELFFAC(MXLAY), SELFFRAC(MXLAY), INDSELF(MXLAY)
      COMMON /FOREIGN/  FORFAC(MXLAY), FORFRAC(MXLAY), INDFOR(MXLAY)
      COMMON /MINOR/    MINORFRAC(MXLAY), INDMINOR(MXLAY), 
     &                  SCALEMINOR(MXLAY),SCALEMINORN2(MXLAY)
      COMMON /CO2/      INDCO2(MXLAY)
      COMMON /K3/       KA(NETA,NTEMP,NREF,MG), FORREF(4,MG),
     &                  SELFREF(10,MG), KA_MN2O(9,19,MG), 
     &                  K_CO2(NETA,NTEMP,NREF,MG)

      COMMON /CVRTAU/    HNAMTAU,HVRTAU

      CHARACTER*18       HNAMTAU,HVRTAU

      REAL KA,K_CO2
      REAL KA_MN2O, MINORFRAC
      REAL N2OM1,N2OM2

      real tau_hi_co2(2,2)

      DIMENSION ABSA(8946,MG),ABSB(8946,MG),ABSC(8946,MG)
      DIMENSION FRACREFA(MG,9)

      integer js(nco2),js1(nco2),ind0(nco2),ind1(nco2)
      real fs(nco2),fs1(nco2)

      real fac000(nco2),fac100(nco2),fac200(nco2)
      real fac010(nco2),fac110(nco2),fac210(nco2)
      real fac001(nco2),fac101(nco2),fac201(nco2)
      real fac011(nco2),fac111(nco2),fac211(nco2)

      real specparm(nco2),specparm1(nco2)
      real speccomb(nco2),speccomb1(nco2)


      real tau_base(mg)
      real tau_major(nco2),tau_major1(nco2)

      EQUIVALENCE (KA,ABSB)
      EQUIVALENCE (K_CO2,ABSC)

C Planck fraction mapping level: P=212.7250 mbar, T = 215.0 K

      DATA (FRACREFA(IG, 1),IG=1,16) /
     &1.6446E-01,1.5641E-01,1.4586E-01,1.3242E-01,1.1605E-01,9.6031E-02,
     &7.7646E-02,5.7616E-02,3.9096E-02,4.2471E-03,3.5126E-03,2.7743E-03,
     &2.0276E-03,1.2947E-03,4.8794E-04,6.8658E-05/
      DATA (FRACREFA(IG, 2),IG=1,16) /
     &1.6059E-01,1.5634E-01,1.4655E-01,1.3382E-01,1.1659E-01,9.6465E-02,
     &7.8143E-02,5.7879E-02,3.9209E-02,4.2472E-03,3.5125E-03,2.7743E-03,
     &2.0276E-03,1.2947E-03,4.8794E-04,6.8658E-05/
      DATA (FRACREFA(IG, 3),IG=1,16) /
     &1.5951E-01,1.5593E-01,1.4700E-01,1.3402E-01,1.1674E-01,9.6557E-02,
     &7.8340E-02,5.8089E-02,3.9392E-02,4.2472E-03,3.5126E-03,2.7743E-03,
     &2.0276E-03,1.2947E-03,4.8794E-04,6.8658E-05/
      DATA (FRACREFA(IG, 4),IG=1,16) /
     &1.5894E-01,1.5560E-01,1.4709E-01,1.3403E-01,1.1689E-01,9.6683E-02,
     &7.8499E-02,5.8239E-02,3.9591E-02,4.2597E-03,3.5126E-03,2.7743E-03,
     &2.0276E-03,1.2947E-03,4.8794E-04,6.8658E-05/
      DATA (FRACREFA(IG, 5),IG=1,16) /
     &1.5844E-01,1.5517E-01,1.4745E-01,1.3378E-01,1.1695E-01,9.6869E-02,
     &7.8685E-02,5.8409E-02,3.9741E-02,4.2901E-03,3.5555E-03,2.7743E-03,
     &2.0276E-03,1.2947E-03,4.8794E-04,6.8658E-05/
      DATA (FRACREFA(IG, 6),IG=1,16) /
     &1.5765E-01,1.5540E-01,1.4720E-01,1.3374E-01,1.1679E-01,9.7229E-02,
     &7.8852E-02,5.8630E-02,3.9946E-02,4.2847E-03,3.5535E-03,2.8251E-03,
     &2.0401E-03,1.2947E-03,4.8794E-04,6.8658E-05/
      DATA (FRACREFA(IG, 7),IG=1,16) /
     &1.5687E-01,1.5526E-01,1.4685E-01,1.3401E-01,1.1643E-01,9.7625E-02,
     &7.9163E-02,5.8915E-02,4.0168E-02,4.3458E-03,3.5684E-03,2.8108E-03,
     &2.1144E-03,1.2996E-03,4.8794E-04,6.8658E-05/
      DATA (FRACREFA(IG, 8),IG=1,16) /
     &1.5628E-01,1.5424E-01,1.4646E-01,1.3376E-01,1.1600E-01,9.8597E-02,
     &7.9630E-02,5.9562E-02,4.0480E-02,4.4003E-03,3.6225E-03,2.8956E-03,
     &2.1223E-03,1.3466E-03,5.1546E-04,7.9453E-05/
      DATA (FRACREFA(IG, 9),IG=1,16) /
     &1.4154E-01,1.4468E-01,1.4298E-01,1.3457E-01,1.2142E-01,1.0519E-01,
     &8.5661E-02,6.4416E-02,4.3398E-02,4.7210E-03,3.8999E-03,3.0833E-03,
     &2.2996E-03,1.4899E-03,5.7453E-04,7.9453E-05/

      co2_lo_hi_rat = 0.04*3.85e-4/0.95

C Minor gas mapping levels:
C     LOWER - N2O, P = 706.272 mbar, T = 275.0 K


C     P = 212.725 mb
      REFRAT_PLANCK_A = CHI_STD(1,21)/CHI_STD(2,21)

C     P = 706.270mb
      REFRAT_M_A = CHI_STD(1,15)/CHI_STD(2,15)

C     Compute the optical depth by interpolating in ln(pressure) and 
C     temperature, and appropriate species.  Below LAYTROP, the water vapor 
C     self-continuum and foreign continuum is interpolated (in temperature) 
C     separately.

      HVRTAU = '$Revision: 22602 $'

      DO 2500 LAY = 1, NLAYERS

c Set up calculation for low and/or high CO2 values
      select case (indco2(lay))
      case (0)
          ico2_0 = 1
          ico2_1 = 1
      case(2)
          ico2_0 = 2
          ico2_1 = 2
      case(1)
          ico2_0 = 1
          ico2_1 = 2
       case default
          print *, 'indco2 in taugb3 has an invalid value'
          stop
      end select

      do ic=ico2_0,ico2_1
         if (ic.eq.1) then
            scale_eta=1.0
         else
            scale_eta= co2_lo_hi_rat
         end if

         SPECCOMB(ic) = COLH2O(LAY) + scale_eta*
     &           RAT_H2OCO2(LAY)*COLCO2(LAY)
         SPECPARM(ic) = COLH2O(LAY)/SPECCOMB(ic)
         IF (SPECPARM(ic) .GE. ONEMINUS) SPECPARM(ic) = ONEMINUS
         SPECMULT = 8.*(SPECPARM(ic))
         JS(ic) = 1 + INT(SPECMULT)
         FS(ic) = AMOD(SPECMULT,1.0)

         SPECCOMB1(ic) = COLH2O(LAY)+scale_eta
     &        *RAT_H2OCO2_1(LAY)*COLCO2(LAY)
         SPECPARM1(ic) = COLH2O(LAY)/SPECCOMB1(ic)
         IF (SPECPARM1(ic) .GE. ONEMINUS) SPECPARM1(ic) = ONEMINUS
         SPECMULT1 = 8.*(SPECPARM1(ic))
         JS1(ic) = 1 + INT(SPECMULT1)
         FS1(ic) = AMOD(SPECMULT1,1.0)

         SPECCOMB_MN2O = COLH2O(LAY) + REFRAT_M_A*COLCO2(LAY)
         SPECPARM_MN2O = COLH2O(LAY)/SPECCOMB_MN2O
         IF (SPECPARM_MN2O .GE. ONEMINUS) SPECPARM_MN2O = ONEMINUS
         SPECMULT_MN2O = 8.*SPECPARM_MN2O
         JMN2O = 1 + INT(SPECMULT_MN2O)
         FMN2O = AMOD(SPECMULT_MN2O,1.0)
         FMN2OMF = MINORFRAC(LAY)*FMN2O

c     In atmospheres where the amount of N2O is too great to be considered
c     a minor species, adjust the column amount of N2O by an empirical factor 
c     to obtain the proper contribution.
         CHI_N2O = COLN2O(LAY)/COLDRY(LAY)
         RATN2O = 1.E20*CHI_N2O/CHI_STD(4,JP(LAY)+1)
         IF (RATN2O .GT. 1.5) THEN
            ADJFAC = 0.5+(RATN2O-0.5)**0.65
            ADJCOLN2O = ADJFAC*CHI_STD(4,JP(LAY)+1)*COLDRY(LAY)*1.E-20
         ELSE
            ADJCOLN2O = COLN2O(LAY)
         ENDIF

         SPECCOMB_PLANCK = COLH2O(LAY)+REFRAT_PLANCK_A*COLCO2(LAY)
         SPECPARM_PLANCK = COLH2O(LAY)/SPECCOMB_PLANCK
         IF (SPECPARM_PLANCK .GE. ONEMINUS) SPECPARM_PLANCK=ONEMINUS
         SPECMULT_PLANCK = 8.*SPECPARM_PLANCK
         JPL= 1 + INT(SPECMULT_PLANCK)
         FPL = AMOD(SPECMULT_PLANCK,1.0)

         IND0(ic) = ((JP(LAY)-1)*NTEMP+(JT(LAY)-1))*NSPA(3) + JS(ic)
         IND1(ic) = (JP(LAY)*NTEMP+(JT1(LAY)-1))*NSPA(3) + JS1(ic)
         INDS = INDSELF(LAY)
         INDF = INDFOR(LAY)
         INDM = INDMINOR(LAY)


         IF (SPECPARM(ic) .LT. 0.125) THEN
            P = FS(ic) - 1
            P4 = P**4
            P8 = P4**2
            P4h = (P8+50*FS(ic)*P4)/(1.+50*FS(ic))
            FK0 = P4h
            FK1 = 1 - P - P4h - P4
            FK2 = P + P4
            FAC000(ic) = FK0*FAC00(LAY)
            FAC100(ic) = FK1*FAC00(LAY)
            FAC200(ic) = FK2*FAC00(LAY)
            FAC010(ic) = FK0*FAC10(LAY)
            FAC110(ic) = FK1*FAC10(LAY)
            FAC210(ic) = FK2*FAC10(LAY)
         ELSE IF (SPECPARM(ic) .GT. 0.875) THEN
            P = -FS(ic)
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC000(ic) = FK0*FAC00(LAY)
            FAC100(ic) = FK1*FAC00(LAY)
            FAC200(ic) = FK2*FAC00(LAY)
            FAC010(ic) = FK0*FAC10(LAY)
            FAC110(ic) = FK1*FAC10(LAY)
            FAC210(ic) = FK2*FAC10(LAY)
         ELSE
            FAC000(ic) = (1. - FS(ic)) * FAC00(LAY)
            FAC010(ic) = (1. - FS(ic)) * FAC10(LAY)
            FAC100(ic) = FS(ic) * FAC00(LAY)
            FAC110(ic) = FS(ic) * FAC10(LAY)
         END IF

       
         IF (SPECPARM1(ic) .LT. 0.125 ) THEN
            P = FS1(ic) - 1
            P4 = P**4
            P8 = P4**2
            P4h = (P8+50*FS1(ic)*P4)/(1.+50*FS1(ic))
            FK0 = P4h
            FK1 = 1 - P - P4h - P4
            FK2 = P + P4
            FAC001(ic) = FK0*FAC01(LAY)
            FAC101(ic) = FK1*FAC01(LAY)
            FAC201(ic) = FK2*FAC01(LAY)
            FAC011(ic) = FK0*FAC11(LAY)
            FAC111(ic) = FK1*FAC11(LAY)
            FAC211(ic) = FK2*FAC11(LAY)
         ELSE IF (SPECPARM1(ic) .GT. 0.875 ) THEN
            P = -FS1(ic)
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC001(ic) = FK0*FAC01(LAY)
            FAC101(ic) = FK1*FAC01(LAY)
            FAC201(ic) = FK2*FAC01(LAY)
            FAC011(ic) = FK0*FAC11(LAY)
            FAC111(ic) = FK1*FAC11(LAY)
            FAC211(ic) = FK2*FAC11(LAY)
         ELSE
            FAC001(ic) = (1. - FS1(ic)) * FAC01(LAY)
            FAC011(ic) = (1. - FS1(ic)) * FAC11(LAY)
            FAC101(ic) = FS1(ic) * FAC01(LAY)
            FAC111(ic) = FS1(ic) * FAC11(LAY)
         END IF

        end do   ! end ic loop for general layer properties

    
         DO 2000 IG = 1, NG(3)
	   DO IC=ico2_0,ico2_1
	    if(ic.eq.1) then
	       absa = absb
	    else
		absa = absc
	    endif 

            IF (SPECPARM(ic).LT.0.125) THEN 
               TAU_MAJOR(ic) = SPECCOMB(ic) *
     &          (FAC000(ic)* ABSA(ind0(ic),IG) +
     &          FAC100(ic)* ABSA(ind0(ic)+1,IG) +
     &          FAC200(ic)* ABSA(ind0(ic)+2,IG) +
     &          FAC010(ic)* ABSA(ind0(ic)+9,IG) +
     &          FAC110(ic)* ABSA(ind0(ic)+10,IG) +
     &          FAC210(ic)* ABSA(ind0(ic)+11,IG))
           ELSE IF (SPECPARM(ic).GT.0.875) THEN 
                TAU_MAJOR(ic) = speccomb(ic) * 
     &              (FAC200(ic)* ABSA(ind0(ic)-1,IG) +
     &              FAC100(ic)* ABSA(ind0(ic),IG) +
     &              FAC000(ic)* ABSA(ind0(ic)+1,IG) +
     &              FAC210(ic)* ABSA(ind0(ic)+8,IG) +
     &              FAC110(ic)* ABSA(ind0(ic)+9,IG) +
     &              FAC010(ic)* ABSA(ind0(ic)+10,IG))
           ELSE 
               TAU_MAJOR(ic) = speccomb(ic) * 
     &              (FAC000(ic) * ABSA(IND0(ic),IG) +
     &              FAC100(ic) * ABSA(IND0(ic)+1,IG) +
     &              FAC010(ic) * ABSA(IND0(ic)+9,IG) +
     &              FAC110(ic) * ABSA(IND0(ic)+10,IG))
           END IF

           IF (SPECPARM1(ic).LT.0.125) THEN 
                TAU_MAJOR1(ic) =  speccomb1(ic) *
     &              (FAC001(ic) * ABSA(IND1(ic),IG) + 
     &              FAC101(ic) * ABSA(IND1(ic)+1,IG) +
     &              FAC201(ic) * ABSA(IND1(ic)+2,IG) +
     &              FAC011(ic) * ABSA(IND1(ic)+9,IG) +
     &              FAC111(ic) * ABSA(IND1(ic)+10,IG) +
     &              FAC211(ic) * ABSA(IND1(ic)+11,IG))

           ELSE IF (SPECPARM1(ic).GT.0.875) THEN 
                TAU_MAJOR1(ic) = speccomb1(ic) * 
     &              (FAC201(ic) * ABSA(IND1(ic)-1,IG) +
     &              FAC101(ic) * ABSA(IND1(ic),IG) +
     &              FAC001(ic) * ABSA(IND1(ic)+1,IG) +
     &              FAC211(ic) * ABSA(IND1(ic)+8,IG) +
     &              FAC111(ic) * ABSA(IND1(ic)+9,IG) +
     &              FAC011(ic) * ABSA(IND1(ic)+10,IG))
           ELSE
               TAU_MAJOR1(ic) = speccomb1(ic) *
     &              (FAC001(ic) * ABSA(IND1(ic),IG) +
     &              FAC101(ic) * ABSA(IND1(ic)+1,IG) +
     &              FAC011(ic) * ABSA(IND1(ic)+9,IG) +
     &              FAC111(ic) * ABSA(IND1(ic)+10,IG))

           END IF
          end do    ! co2 loop

C  Do  interpolation for high CO2 cases if needed
          select case (indco2(lay) )
          case (0)
             TAUG(LAY,IG) = TAU_MAJOR(1)+TAU_MAJOR1(1)
          case (1)
             TAUG(LAY,IG) = co2_fac0(lay)*(TAU_MAJOR(1)+TAU_MAJOR1(1)) + 
     &                      co2_fac1(lay)*(TAU_MAJOR(2)+TAU_MAJOR1(2)) 
          case (2)
             TAUG(LAY,IG) = TAU_MAJOR(2)+TAU_MAJOR1(2)
          case default
             print *, 'indco2 in taugb3 has an invalid value'
             stop
          end select


           FRACS(LAY,IG) = FRACREFA(IG,JPL) + FPL *
     &              (FRACREFA(IG,JPL+1)-FRACREFA(IG,JPL))

 2000          CONTINUE
c Complete calculation and do  interpolation for high CO2 cases if needed

        DO 2040 IG = 1, NG(3)
               TAUSELF = SELFFAC(LAY)* (SELFREF(INDS,IG) +
     &              SELFFRAC(LAY)  *
     &              (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
               TAUFOR = FORFAC(LAY) * (FORREF(INDF,IG) +
     &              FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &              FORREF(INDF,IG))) 
               N2OM1 = KA_MN2O(JMN2O,INDM,IG) + FMN2O*
     &              (KA_MN2O(JMN2O+1,INDM,IG) - 
     &              KA_MN2O(JMN2O,INDM,IG))
               N2OM2 = KA_MN2O(JMN2O,INDM+1,IG) + FMN2O*
     &              (KA_MN2O(JMN2O+1,INDM+1,IG) - 
     &              KA_MN2O(JMN2O,INDM+1,IG))
               ABSN2O = N2OM1 + MINORFRAC(LAY) *
     &              (N2OM2 - N2OM1)

C                    taug(lay,ig) = taug(lay,ig)
C      &              + TAUSELF + TAUFOR
C      &              + ADJCOLN2O*ABSN2O   

                  taug(lay,ig) = taug(lay,ig) + ADJCOLN2O*ABSN2O
               if( h2o_sb == 1 ) then 
                  taug(lay,ig) = taug(lay,ig) + TAUSELF
                endif
               if( h2o_for == 1 ) then 
                  taug(lay,ig) = taug(lay,ig) + TAUFOR
                endif

 2040     CONTINUE
 2500 CONTINUE


      RETURN
      END

C----------------------------------------------------------------------------

C----------------------------------------------------------------------------

      SUBROUTINE TAUGB4
      use variables !NJE

C     BAND 4:  630-700 cm-1 (low key - H2O,CO2; high - none)

      PARAMETER (mg=16, mxlay=203, NBANDS=16,NREF=71)
      PARAMETER (NTEMP=14,NETA=9)

C  Output

      COMMON /TAUGCOM/  TAUG(MXLAY,MG)
      COMMON /PLANKG/   FRACS(MXLAY,MG)

C  Input

      COMMON /FEATURES/ NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)
      COMMON /PRECISE/  ONEMINUS
      COMMON /PROFILE/  NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
     &                  PZ(0:MXLAY),TZ(0:MXLAY)
      COMMON /PROFDATA/ laytrop,lay_switch_kabs,                                   
     &                  COLH2O(MXLAY),COLCO2(MXLAY),COLO3(MXLAY),  
     &                  COLN2O(MXLAY),COLCO(MXLAY),COLCH4(MXLAY),  
     &                  COLO2(MXLAY),COLBRD(MXLAY)
      COMMON /STD_REF/  PREF(NREF),PREFLOG(NREF),TREF(NREF),
     &                  CHI_STD(7,NREF)
      COMMON /INTFAC/   FAC00(MXLAY),FAC01(MXLAY),                            
     &                  FAC10(MXLAY),FAC11(MXLAY)                             
      COMMON /INTIND/   JP(MXLAY),JT(MXLAY),JT1(MXLAY)
      COMMON /REFRAT_ETA/ RAT_H2OCO2(MXLAY),RAT_H2OCO2_1(MXLAY),
     &                  RAT_H2OO3(MXLAY),RAT_H2OO3_1(MXLAY),
     &                  RAT_H2ON2O(MXLAY),RAT_H2ON2O_1(MXLAY),
     &                  RAT_H2OCH4(MXLAY),RAT_H2OCH4_1(MXLAY),
     &                  RAT_N2OCO2(MXLAY),RAT_N2OCO2_1(MXLAY),
     &                  RAT_O3CO2(MXLAY),RAT_O3CO2_1(MXLAY),
     &                  RAT_CO2O3(MXLAY),RAT_CO2O3_1(MXLAY),
     &                  RAT_O2_DRY(MXLAY)

      COMMON /FOREIGN/  FORFAC(MXLAY), FORFRAC(MXLAY), INDFOR(MXLAY)
      COMMON /SELF/     SELFFAC(MXLAY), SELFFRAC(MXLAY), INDSELF(MXLAY)
      COMMON /K4/       KA(NETA,NTEMP,NREF,MG), FORREF(4,MG), 
     &                  SELFREF(10,MG)

      COMMON /CVRTAU/    HNAMTAU,HVRTAU

      CHARACTER*18       HNAMTAU,HVRTAU

      DIMENSION ABSA(8946,MG)
      DIMENSION FRACREFA(MG,9)

C Planck fraction mapping level : P = 707.27 mbar, T = 230.0 K

      DATA (FRACREFA(IG, 1),IG=1,16) /
     &1.5577E-01,1.4973E-01,1.4163E-01,1.3106E-01,1.1764E-01,1.0137E-01,
     &8.2819E-02,6.2234E-02,4.2101E-02,4.5934E-03,3.8094E-03,3.0124E-03,
     &2.2106E-03,1.4082E-03,5.3085E-04,7.4698E-05/
      DATA (FRACREFA(IG, 2),IG=1,16) /
     &1.5576E-01,1.4973E-01,1.4162E-01,1.3107E-01,1.1764E-01,1.0139E-01,
     &8.2820E-02,6.2235E-02,4.2101E-02,4.5934E-03,3.8094E-03,3.0124E-03,
     &2.2106E-03,1.4082E-03,5.3085E-04,7.4698E-05/
      DATA (FRACREFA(IG, 3),IG=1,16) /
     &1.5586E-01,1.4972E-01,1.4150E-01,1.3109E-01,1.1763E-01,1.0140E-01,
     &8.2820E-02,6.2236E-02,4.2101E-02,4.5934E-03,3.8094E-03,3.0124E-03,
     &2.2106E-03,1.4082E-03,5.3085E-04,7.4698E-05/
      DATA (FRACREFA(IG, 4),IG=1,16) /
     &1.5596E-01,1.4952E-01,1.4157E-01,1.3111E-01,1.1764E-01,1.0140E-01,
     &8.2832E-02,6.2237E-02,4.2101E-02,4.5934E-03,3.8094E-03,3.0124E-03,
     &2.2106E-03,1.4082E-03,5.3085E-04,7.4698E-05/
      DATA (FRACREFA(IG, 5),IG=1,16) /
     &1.5595E-01,1.4949E-01,1.4144E-01,1.3123E-01,1.1768E-01,1.0137E-01,
     &8.2872E-02,6.2239E-02,4.2101E-02,4.5934E-03,3.8094E-03,3.0124E-03,
     &2.2106E-03,1.4082E-03,5.3085E-04,7.4698E-05/
      DATA (FRACREFA(IG, 6),IG=1,16) /
     &1.5585E-01,1.4955E-01,1.4140E-01,1.3126E-01,1.1769E-01,1.0138E-01,
     &8.2900E-02,6.2243E-02,4.2101E-02,4.5934E-03,3.8094E-03,3.0124E-03,
     &2.2106E-03,1.4082E-03,5.3085E-04,7.4698E-05/
      DATA (FRACREFA(IG, 7),IG=1,16) /
     &1.5543E-01,1.4981E-01,1.4139E-01,1.3127E-01,1.1770E-01,1.0149E-01,
     &8.2877E-02,6.2297E-02,4.2101E-02,4.5934E-03,3.8094E-03,3.0124E-03,
     &2.2106E-03,1.4082E-03,5.3085E-04,7.4698E-05/
      DATA (FRACREFA(IG, 8),IG=1,16) /
     &1.5320E-01,1.5094E-01,1.4184E-01,1.3151E-01,1.1777E-01,1.0158E-01,
     &8.3009E-02,6.2370E-02,4.2128E-02,4.5935E-03,3.8094E-03,3.0124E-03,
     &2.2106E-03,1.4082E-03,5.3085E-04,7.4698E-05/
      DATA (FRACREFA(IG, 9),IG=1,16) /
     &1.4678E-01,1.4599E-01,1.4341E-01,1.3281E-01,1.2073E-01,1.0405E-01,
     &8.4456E-02,6.3165E-02,4.2631E-02,4.6667E-03,3.8945E-03,3.0677E-03,
     &2.2527E-03,1.4546E-03,5.6166E-04,7.9032E-05/

      EQUIVALENCE (KA,ABSA),(KB,ABSB)
      REAL KA,KB

C     P =   317.348
c     REFRAT_PLANCK_A = CHI_STD(1,11)/CHI_STD(2,11)

C     P = 707.27
      REFRAT_PLANCK_A = CHI_STD(1,15)/CHI_STD(2,15)

C     Compute the optical depth by interpolating in ln(pressure) and 
C     temperature, and appropriate species; the water 
C     vapor self-continuum and foreign continuum is interpolated (in temperature) 
C     separately.

      HVRTAU = '$Revision: 22602 $'

      DO 2500 LAY = 1, NLAYERS

         SPECCOMB = COLH2O(LAY) + RAT_H2OCO2(LAY)*COLCO2(LAY)
         SPECPARM = COLH2O(LAY)/SPECCOMB

         IF (SPECPARM .GE. ONEMINUS) SPECPARM = ONEMINUS
         SPECMULT = 8.*(SPECPARM)
         JS = 1 + INT(SPECMULT)
         FS = AMOD(SPECMULT,1.0)

         SPECCOMB1 = COLH2O(LAY) + RAT_H2OCO2_1(LAY)*COLCO2(LAY)
         SPECPARM1 = COLH2O(LAY)/SPECCOMB1
         IF (SPECPARM1 .GE. ONEMINUS) SPECPARM1 = ONEMINUS
         SPECMULT1 = 8.*(SPECPARM1)
         JS1 = 1 + INT(SPECMULT1)
         FS1 = AMOD(SPECMULT1,1.0)

         SPECCOMB_PLANCK = COLH2O(LAY)+REFRAT_PLANCK_A*COLCO2(LAY)
         SPECPARM_PLANCK = COLH2O(LAY)/SPECCOMB_PLANCK
         IF (SPECPARM_PLANCK .GE. ONEMINUS) SPECPARM_PLANCK=ONEMINUS
         SPECMULT_PLANCK = 8.*SPECPARM_PLANCK
         JPL= 1 + INT(SPECMULT_PLANCK)
         FPL = AMOD(SPECMULT_PLANCK,1.0)

         IND0 = ((JP(LAY)-1)*NTEMP+(JT(LAY)-1))*NSPA(4) + JS
         IND1 = (JP(LAY)*NTEMP+(JT1(LAY)-1))*NSPA(4) + JS1
         INDS = INDSELF(LAY)
         INDF = INDFOR(LAY)

         IF (SPECPARM .LT. 0.125 .AND. SPECPARM1 .LT. 0.125) THEN
            P = FS - 1
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC000 = FK0*FAC00(LAY)
            FAC100 = FK1*FAC00(LAY)
            FAC200 = FK2*FAC00(LAY)
            FAC010 = FK0*FAC10(LAY)
            FAC110 = FK1*FAC10(LAY)
            FAC210 = FK2*FAC10(LAY)

            P = FS1 - 1
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC001 = FK0*FAC01(LAY)
            FAC101 = FK1*FAC01(LAY)
            FAC201 = FK2*FAC01(LAY)
            FAC011 = FK0*FAC11(LAY)
            FAC111 = FK1*FAC11(LAY)
            FAC211 = FK2*FAC11(LAY)

            DO 2000 IG = 1, NG(4)
               TAUSELF = SELFFAC(LAY)* (SELFREF(INDS,IG) +
     &              SELFFRAC(LAY)  *
     &              (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
               TAUFOR =  FORFAC(LAY) * (FORREF(INDF,IG) +
     &              FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &              FORREF(INDF,IG))) 
               TAUG(LAY,IG) = SPECCOMB *
     &              (FAC000 * ABSA(IND0,IG) +
     &              FAC100 * ABSA(IND0+1,IG) +
     &              FAC200 * ABSA(IND0+2,IG) +
     &              FAC010 * ABSA(IND0+9,IG) +
     &              FAC110 * ABSA(IND0+10,IG) +
     &              FAC210 * ABSA(IND0+11,IG))
     &              + SPECCOMB1 *
     &              (FAC001 * ABSA(IND1,IG) + 
     &              FAC101 * ABSA(IND1+1,IG) +
     &              FAC201 * ABSA(IND1+2,IG) +
     &              FAC011 * ABSA(IND1+9,IG) +
     &              FAC111 * ABSA(IND1+10,IG) +
     &              FAC211 * ABSA(IND1+11,IG)) 
C      &              + TAUSELF + TAUFOR

              taug(lay,ig) = taug(lay,ig)
               if( h2o_sb == 1 ) then 
                  taug(lay,ig) = taug(lay,ig) + TAUSELF
                endif
               if( h2o_for == 1 ) then 
                  taug(lay,ig) = taug(lay,ig) + TAUFOR
                endif

               FRACS(LAY,IG) = FRACREFA(IG,JPL) + FPL *
     &          (FRACREFA(IG,JPL+1)-FRACREFA(IG,JPL))
 2000          CONTINUE
         ELSE IF (SPECPARM .GT. 0.875 .AND. SPECPARM1 .GT. 0.875) THEN
            P = -FS 
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC000 = FK0*FAC00(LAY)
            FAC100 = FK1*FAC00(LAY)
            FAC200 = FK2*FAC00(LAY)
            FAC010 = FK0*FAC10(LAY)
            FAC110 = FK1*FAC10(LAY)
            FAC210 = FK2*FAC10(LAY)

            P = -FS1 
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC001 = FK0*FAC01(LAY)
            FAC101 = FK1*FAC01(LAY)
            FAC201 = FK2*FAC01(LAY)
            FAC011 = FK0*FAC11(LAY)
            FAC111 = FK1*FAC11(LAY)
            FAC211 = FK2*FAC11(LAY)
            DO 2010 IG = 1, NG(4)
               TAUSELF = SELFFAC(LAY)* (SELFREF(INDS,IG) +
     &              SELFFRAC(LAY)  *
     &              (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
               TAUFOR =  FORFAC(LAY) * (FORREF(INDF,IG) +
     &              FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &              FORREF(INDF,IG))) 
               TAUG(LAY,IG) = SPECCOMB * 
     &              (FAC200 * ABSA(IND0-1,IG) +
     &              FAC100 * ABSA(IND0,IG) +
     &              FAC000 * ABSA(IND0+1,IG) +
     &              FAC210 * ABSA(IND0+8,IG) +
     &              FAC110 * ABSA(IND0+9,IG) +
     &              FAC010 * ABSA(IND0+10,IG))
     &              + SPECCOMB1 *
     &              (FAC201 * ABSA(IND1-1,IG) +
     &              FAC101 * ABSA(IND1,IG) +
     &              FAC001 * ABSA(IND1+1,IG) +
     &              FAC211 * ABSA(IND1+8,IG) +
     &              FAC111 * ABSA(IND1+9,IG) +
     &              FAC011 * ABSA(IND1+10,IG))
C      &              + TAUSELF + TAUFOR

              taug(lay,ig) = taug(lay,ig)
               if( h2o_sb == 1 ) then 
                  taug(lay,ig) = taug(lay,ig) + TAUSELF
                endif
               if( h2o_for == 1 ) then 
                  taug(lay,ig) = taug(lay,ig) + TAUFOR
                endif

               FRACS(LAY,IG) = FRACREFA(IG,JPL) + FPL *
     &          (FRACREFA(IG,JPL+1)-FRACREFA(IG,JPL))
 2010           CONTINUE
         ELSE
            FAC000 = (1. - FS) * FAC00(LAY)
            FAC010 = (1. - FS) * FAC10(LAY)
            FAC100 = FS * FAC00(LAY)
            FAC110 = FS * FAC10(LAY)

            FAC001 = (1. - FS1) * FAC01(LAY)
            FAC011 = (1. - FS1) * FAC11(LAY)
            FAC101 = FS1 * FAC01(LAY)
            FAC111 = FS1 * FAC11(LAY)

            DO 2020 IG = 1, NG(4)
               TAUSELF = SELFFAC(LAY)* (SELFREF(INDS,IG) +
     &              SELFFRAC(LAY)  *
     &              (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
               TAUFOR =  FORFAC(LAY) * (FORREF(INDF,IG) +
     &              FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &              FORREF(INDF,IG))) 
               TAUG(LAY,IG) = SPECCOMB * 
     &              (FAC000 * ABSA(IND0,IG) +
     &              FAC100 * ABSA(IND0+1,IG) +
     &              FAC010 * ABSA(IND0+9,IG) +
     &              FAC110 * ABSA(IND0+10,IG))
     &              + SPECCOMB1 *
     &              (FAC001 * ABSA(IND1,IG) + 
     &              FAC101 * ABSA(IND1+1,IG) +
     &              FAC011 * ABSA(IND1+9,IG) +
     &              FAC111 * ABSA(IND1+10,IG)) 
C      &              + TAUSELF + TAUFOR

                taug(lay,ig) = taug(lay,ig)
               if( h2o_sb == 1 ) then 
                  taug(lay,ig) = taug(lay,ig) + TAUSELF
                endif
               if( h2o_for == 1 ) then 
                  taug(lay,ig) = taug(lay,ig) + TAUFOR
                endif

               FRACS(LAY,IG) = FRACREFA(IG,JPL) + FPL *
     &          (FRACREFA(IG,JPL+1)-FRACREFA(IG,JPL))
 2020          CONTINUE
        ENDIF
 2500 CONTINUE

      RETURN
      END
C----------------------------------------------------------------------------

      SUBROUTINE TAUGB5
      use variables !NJE

C     BAND 5:  700-820 cm-1 (low key - H2O,CO2; low minor - O3)


      PARAMETER (mg=16, mxlay=203, MXMOL=38, MAXXSEC=4,NBANDS=16)
      PARAMETER (NTEMP=14,NETA=9,NREF=71,NCO2=2)

C  Output

      COMMON /TAUGCOM/  TAUG(MXLAY,MG)
      COMMON /PLANKG/   FRACS(MXLAY,MG)                                       

C  Input

      COMMON /FEATURES/ NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)
      COMMON /PRECISE/  ONEMINUS
      COMMON /PROFILE/  NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
     &                  PZ(0:MXLAY),TZ(0:MXLAY)
      COMMON /PROFDATA/ LAYTROP, LAY_SWITCH_KABS,                                  
     &                  COLH2O(MXLAY),COLCO2(MXLAY),COLO3(MXLAY),
     &                  COLN2O(MXLAY),COLCO(MXLAY),COLCH4(MXLAY),
     &                  COLO2(MXLAY),COLBRD(MXLAY)
      COMMON /SPECIES/  COLDRY(MXLAY),WKL(MXMOL,MXLAY),WBROAD(MXLAY),
     &                  COLMOL(MXLAY),NMOL
      COMMON /INTFAC/   FAC00(MXLAY),FAC01(MXLAY),                            
     &                  FAC10(MXLAY),FAC11(MXLAY)
      COMMON /CO2FAC/   co2_fac0(mxlay), co2_fac1(mxlay)
      COMMON /STD_REF/  PREF(NREF),PREFLOG(NREF),TREF(NREF),
     &                  CHI_STD(7,NREF)
      COMMON /XSEC/     WX(MAXXSEC,MXLAY)
      COMMON /INTIND/   JP(MXLAY),JT(MXLAY),JT1(MXLAY)
      COMMON /REFRAT_ETA/ RAT_H2OCO2(MXLAY),RAT_H2OCO2_1(MXLAY),
     &                  RAT_H2OO3(MXLAY),RAT_H2OO3_1(MXLAY),
     &                  RAT_H2ON2O(MXLAY),RAT_H2ON2O_1(MXLAY),
     &                  RAT_H2OCH4(MXLAY),RAT_H2OCH4_1(MXLAY),
     &                  RAT_N2OCO2(MXLAY),RAT_N2OCO2_1(MXLAY),
     &                  RAT_O3CO2(MXLAY),RAT_O3CO2_1(MXLAY),
     &                  RAT_CO2O3(MXLAY),RAT_CO2O3_1(MXLAY),
     &                  RAT_O2_DRY(MXLAY)

      COMMON /SELF/     SELFFAC(MXLAY), SELFFRAC(MXLAY), INDSELF(MXLAY)
      COMMON /FOREIGN/  FORFAC(MXLAY), FORFRAC(MXLAY), INDFOR(MXLAY)
      COMMON /MINOR/    MINORFRAC(MXLAY), INDMINOR(MXLAY),
     &                  SCALEMINOR(MXLAY),SCALEMINORN2(MXLAY)
      COMMON /CO2/      INDCO2(MXLAY)
      COMMON /K5/       KA(NETA,NTEMP,NREF,MG), FORREF(4,MG),
     &                  SELFREF(10,MG), KA_MO3(9,19,MG),
     &                  K_CO2(NETA,NTEMP,NREF,MG)

      COMMON /CVRTAU/    HNAMTAU,HVRTAU

      CHARACTER*18       HNAMTAU,HVRTAU

      REAL KA,K_CO2
      REAL KA_MO3, MINORFRAC
      REAL O2M1,O3M2

      DIMENSION ABSA(8946,MG),ABSB(8946,MG),ABSC(8946,MG)
      DIMENSION FRACREFA(MG,9),CCL4(MG)

      EQUIVALENCE (KA,ABSB)
      EQUIVALENCE (K_CO2,ABSC)
      integer js(nco2),js1(nco2),ind0(nco2),ind1(nco2)
      real fs(nco2),fs1(nco2)

      real fac000(nco2),fac100(nco2),fac200(nco2)
      real fac010(nco2),fac110(nco2),fac210(nco2)
      real fac001(nco2),fac101(nco2),fac201(nco2)
      real fac011(nco2),fac111(nco2),fac211(nco2)

      real specparm(nco2),specparm1(nco2)
      real speccomb(nco2),speccomb1(nco2)


      real tau_base(mg)
      real tau_major(nco2),tau_major1(nco2)


C Planck fraction mapping level : P = 473.42 mb, T = 260.0

      DATA (FRACREFA(IG, 1),IG=1,16) /
     &1.4123E-01,1.4220E-01,1.3797E-01,1.3093E-01,1.2241E-01,1.0697E-01,
     &8.8706E-02,6.7130E-02,4.5505E-02,4.9860E-03,4.1205E-03,3.2547E-03,
     &2.3787E-03,1.5464E-03,5.8411E-04,8.2259E-05/
      DATA (FRACREFA(IG, 2),IG=1,16) /
     &1.4159E-01,1.4265E-01,1.3779E-01,1.3068E-01,1.2221E-01,1.0680E-01,
     &8.8684E-02,6.7141E-02,4.5503E-02,4.9860E-03,4.1205E-03,3.2547E-03,
     &2.3787E-03,1.5464E-03,5.8411E-04,8.2259E-05/
      DATA (FRACREFA(IG, 3),IG=1,16) /
     &1.4164E-01,1.4291E-01,1.3780E-01,1.3080E-01,1.2201E-01,1.0668E-01,
     &8.8547E-02,6.7145E-02,4.5507E-02,4.9868E-03,4.1205E-03,3.2547E-03,
     &2.3787E-03,1.5464E-03,5.8411E-04,8.2259E-05/
      DATA (FRACREFA(IG, 4),IG=1,16) /
     &1.4163E-01,1.4329E-01,1.3776E-01,1.3104E-01,1.2182E-01,1.0652E-01,
     &8.8383E-02,6.7099E-02,4.5508E-02,4.9872E-03,4.1205E-03,3.2547E-03,
     &2.3787E-03,1.5464E-03,5.8411E-04,8.2259E-05/
      DATA (FRACREFA(IG, 5),IG=1,16) /
     &1.4160E-01,1.4361E-01,1.3775E-01,1.3127E-01,1.2174E-01,1.0630E-01,
     &8.8356E-02,6.6923E-02,4.5504E-02,4.9875E-03,4.1205E-03,3.2547E-03,
     &2.3787E-03,1.5464E-03,5.8411E-04,8.2259E-05/
      DATA (FRACREFA(IG, 6),IG=1,16) /
     &1.4151E-01,1.4396E-01,1.3774E-01,1.3154E-01,1.2173E-01,1.0611E-01,
     &8.8245E-02,6.6712E-02,4.5491E-02,4.9903E-03,4.1206E-03,3.2547E-03,
     &2.3787E-03,1.5464E-03,5.8411E-04,8.2259E-05/
      DATA (FRACREFA(IG, 7),IG=1,16) /
     &1.4126E-01,1.4429E-01,1.3790E-01,1.3200E-01,1.2162E-01,1.0606E-01,
     &8.7972E-02,6.6633E-02,4.5293E-02,4.9913E-03,4.1249E-03,3.2547E-03,
     &2.3787E-03,1.5464E-03,5.8411E-04,8.2259E-05/
      DATA (FRACREFA(IG, 8),IG=1,16) /
     &1.4081E-01,1.4413E-01,1.3875E-01,1.3271E-01,1.2154E-01,1.0615E-01,
     &8.7636E-02,6.6318E-02,4.4995E-02,4.9805E-03,4.1186E-03,3.2597E-03,
     &2.3811E-03,1.5464E-03,5.8411E-04,8.2259E-05/
      DATA (FRACREFA(IG, 9),IG=1,16) /
     &1.4162E-01,1.4564E-01,1.4355E-01,1.3565E-01,1.2181E-01,1.0418E-01,
     &8.4864E-02,6.3577E-02,4.2940E-02,4.7289E-03,3.9138E-03,3.0872E-03,
     &2.2645E-03,1.5187E-03,5.7499E-04,8.1659E-05/


C Minor gas mapping level :
C     LOWER - O3, P = 317.34 mbar, T = 245. K
C     LOWER - CCL4

      DATA CCL4/
     &     26.1407,  53.9776,  63.8085,  36.1701,
     &     15.4099, 10.23116,  4.82948,  5.03836,
     &     1.75558,0.,0.,0.,
     &     0.,0.,0.,0./



      co2_lo_hi_rat = 0.04*3.85e-4/0.95

C     Calculate reference ratio to be used in calculation of Planck
C     fraction in lower/upper atmosphere.

C     P = 473.420 mb
      REFRAT_PLANCK_A = CHI_STD(1,17)/CHI_STD(2,17)

C     P = 317.3480
      REFRAT_M_A = CHI_STD(1,19)/CHI_STD(2,19)

C     Compute the optical depth by interpolating in ln(pressure) and 
C     temperature, and appropriate species; the 
C     water vapor self-continuum and foreign continuum is 
C     interpolated (in temperature) separately.

      DO 2500 LAY = 1, NLAYERS

c Set up calculation for low and/or high CO2 values
      select case (indco2(lay))
      case (0)
          ico2_0 = 1
          ico2_1 = 1
      case(2)
          ico2_0 = 2
          ico2_1 = 2
      case(1)
          ico2_0 = 1
          ico2_1 = 2
       case default
          print *, 'indco2 in taugb3 has an invalid value'
          stop
      end select

      do ic=ico2_0,ico2_1
         if (ic.eq.1) then
            scale_eta=1.0
         else
            scale_eta= co2_lo_hi_rat
         end if

         SPECCOMB(ic) = COLH2O(LAY) + scale_eta*
     &           RAT_H2OCO2(LAY)*COLCO2(LAY)
         SPECPARM(ic) = COLH2O(LAY)/SPECCOMB(ic)
         IF (SPECPARM(ic) .GE. ONEMINUS) SPECPARM(ic) = ONEMINUS
         SPECMULT = 8.*(SPECPARM(ic))
         JS(ic) = 1 + INT(SPECMULT)
         FS(ic) = AMOD(SPECMULT,1.0)

         SPECCOMB1(ic) = COLH2O(LAY)+scale_eta
     &        *RAT_H2OCO2_1(LAY)*COLCO2(LAY)
         SPECPARM1(ic) = COLH2O(LAY)/SPECCOMB1(ic)
         IF (SPECPARM1(ic) .GE. ONEMINUS) SPECPARM1(ic) = ONEMINUS
         SPECMULT1 = 8.*(SPECPARM1(ic))
         JS1(ic) = 1 + INT(SPECMULT1)
         FS1(ic) = AMOD(SPECMULT1,1.0)

         SPECCOMB_MO3 = COLH2O(LAY) + REFRAT_M_A*COLCO2(LAY)
         SPECPARM_MO3 = COLH2O(LAY)/SPECCOMB_MO3
         IF (SPECPARM_MO3 .GE. ONEMINUS) SPECPARM_MO3 = ONEMINUS
         SPECMULT_MO3 = 8.*SPECPARM_MO3
         JMO3 = 1 + INT(SPECMULT_MO3)
         FMO3 = AMOD(SPECMULT_MO3,1.0)
         FMO3MF = MINORFRAC(LAY)*FMO3

         SPECCOMB_PLANCK = COLH2O(LAY)+REFRAT_PLANCK_A*COLCO2(LAY)
         SPECPARM_PLANCK = COLH2O(LAY)/SPECCOMB_PLANCK
         IF (SPECPARM_PLANCK .GE. ONEMINUS) SPECPARM_PLANCK=ONEMINUS
         SPECMULT_PLANCK = 8.*SPECPARM_PLANCK
         JPL= 1 + INT(SPECMULT_PLANCK)
         FPL = AMOD(SPECMULT_PLANCK,1.0)

         IND0(ic) = ((JP(LAY)-1)*NTEMP+(JT(LAY)-1))*NSPA(5) + JS(ic)
         IND1(ic) = (JP(LAY)*NTEMP+(JT1(LAY)-1))*NSPA(5) + JS1(ic)
         INDS = INDSELF(LAY)
         INDF = INDFOR(LAY)
         INDM = INDMINOR(LAY)


         IF (SPECPARM(ic) .LT. 0.125) THEN
            P = FS(ic) - 1
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P -  2.0*P4
            FK2 = P + P4
            FAC000(ic) = FK0*FAC00(LAY)
            FAC100(ic) = FK1*FAC00(LAY)
            FAC200(ic) = FK2*FAC00(LAY)
            FAC010(ic) = FK0*FAC10(LAY)
            FAC110(ic) = FK1*FAC10(LAY)
            FAC210(ic) = FK2*FAC10(LAY)
         ELSE IF (SPECPARM(ic) .GT. 0.875) THEN
            P = -FS(ic)
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC000(ic) = FK0*FAC00(LAY)
            FAC100(ic) = FK1*FAC00(LAY)
            FAC200(ic) = FK2*FAC00(LAY)
            FAC010(ic) = FK0*FAC10(LAY)
            FAC110(ic) = FK1*FAC10(LAY)
            FAC210(ic) = FK2*FAC10(LAY)
         ELSE
            FAC000(ic) = (1. - FS(ic)) * FAC00(LAY)
            FAC010(ic) = (1. - FS(ic)) * FAC10(LAY)
            FAC100(ic) = FS(ic) * FAC00(LAY)
            FAC110(ic) = FS(ic) * FAC10(LAY)
         END IF

       
         IF (SPECPARM1(ic) .LT. 0.125 ) THEN
            P = FS1(ic) - 1
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.*P4
            FK2 = P + P4
            FAC001(ic) = FK0*FAC01(LAY)
            FAC101(ic) = FK1*FAC01(LAY)
            FAC201(ic) = FK2*FAC01(LAY)
            FAC011(ic) = FK0*FAC11(LAY)
            FAC111(ic) = FK1*FAC11(LAY)
            FAC211(ic) = FK2*FAC11(LAY)
         ELSE IF (SPECPARM1(ic) .GT. 0.875 ) THEN
            P = -FS1(ic)
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC001(ic) = FK0*FAC01(LAY)
            FAC101(ic) = FK1*FAC01(LAY)
            FAC201(ic) = FK2*FAC01(LAY)
            FAC011(ic) = FK0*FAC11(LAY)
            FAC111(ic) = FK1*FAC11(LAY)
            FAC211(ic) = FK2*FAC11(LAY)
         ELSE
            FAC001(ic) = (1. - FS1(ic)) * FAC01(LAY)
            FAC011(ic) = (1. - FS1(ic)) * FAC11(LAY)
            FAC101(ic) = FS1(ic) * FAC01(LAY)
            FAC111(ic) = FS1(ic) * FAC11(LAY)
         END IF

        end do   ! end ic loop for general layer properties

    
         DO 2000 IG = 1, NG(5)
	   DO IC=ico2_0,ico2_1
	    if(ic.eq.1) then
	       absa = absb
	    else
		absa = absc
	    endif 

            IF (SPECPARM(ic).LT.0.125) THEN 
               TAU_MAJOR(ic) = SPECCOMB(ic) *
     &          (FAC000(ic)* ABSA(ind0(ic),IG) +
     &          FAC100(ic)* ABSA(ind0(ic)+1,IG) +
     &          FAC200(ic)* ABSA(ind0(ic)+2,IG) +
     &          FAC010(ic)* ABSA(ind0(ic)+9,IG) +
     &          FAC110(ic)* ABSA(ind0(ic)+10,IG) +
     &          FAC210(ic)* ABSA(ind0(ic)+11,IG))
           ELSE IF (SPECPARM(ic).GT.0.875) THEN 
                TAU_MAJOR(ic) = speccomb(ic) * 
     &              (FAC200(ic)* ABSA(ind0(ic)-1,IG) +
     &              FAC100(ic)* ABSA(ind0(ic),IG) +
     &              FAC000(ic)* ABSA(ind0(ic)+1,IG) +
     &              FAC210(ic)* ABSA(ind0(ic)+8,IG) +
     &              FAC110(ic)* ABSA(ind0(ic)+9,IG) +
     &              FAC010(ic)* ABSA(ind0(ic)+10,IG))
           ELSE 
               TAU_MAJOR(ic) = speccomb(ic) * 
     &              (FAC000(ic) * ABSA(IND0(ic),IG) +
     &              FAC100(ic) * ABSA(IND0(ic)+1,IG) +
     &              FAC010(ic) * ABSA(IND0(ic)+9,IG) +
     &              FAC110(ic) * ABSA(IND0(ic)+10,IG))
           END IF

           IF (SPECPARM1(ic).LT.0.125) THEN 
                TAU_MAJOR1(ic) =  speccomb1(ic) *
     &              (FAC001(ic) * ABSA(IND1(ic),IG) + 
     &              FAC101(ic) * ABSA(IND1(ic)+1,IG) +
     &              FAC201(ic) * ABSA(IND1(ic)+2,IG) +
     &              FAC011(ic) * ABSA(IND1(ic)+9,IG) +
     &              FAC111(ic) * ABSA(IND1(ic)+10,IG) +
     &              FAC211(ic) * ABSA(IND1(ic)+11,IG))

           ELSE IF (SPECPARM1(ic).GT.0.875) THEN 
                TAU_MAJOR1(ic) = speccomb1(ic) * 
     &              (FAC201(ic) * ABSA(IND1(ic)-1,IG) +
     &              FAC101(ic) * ABSA(IND1(ic),IG) +
     &              FAC001(ic) * ABSA(IND1(ic)+1,IG) +
     &              FAC211(ic) * ABSA(IND1(ic)+8,IG) +
     &              FAC111(ic) * ABSA(IND1(ic)+9,IG) +
     &              FAC011(ic) * ABSA(IND1(ic)+10,IG))
           ELSE
               TAU_MAJOR1(ic) = speccomb1(ic) *
     &              (FAC001(ic) * ABSA(IND1(ic),IG) +
     &              FAC101(ic) * ABSA(IND1(ic)+1,IG) +
     &              FAC011(ic) * ABSA(IND1(ic)+9,IG) +
     &              FAC111(ic) * ABSA(IND1(ic)+10,IG))

           END IF
          end do    ! co2 loop

C  Do  interpolation for high CO2 cases if needed
          select case (indco2(lay) )
          case (0)
             TAUG(LAY,IG) = TAU_MAJOR(1)+TAU_MAJOR1(1)
          case (1)
             TAUG(LAY,IG) = co2_fac0(lay)*(TAU_MAJOR(1)+TAU_MAJOR1(1)) + 
     &                      co2_fac1(lay)*(TAU_MAJOR(2)+TAU_MAJOR1(2)) 
          case (2)
             TAUG(LAY,IG) = TAU_MAJOR(2)+TAU_MAJOR1(2)
          case default
             print *, 'indco2 in taugb3 has an invalid value'
             stop
          end select


           FRACS(LAY,IG) = FRACREFA(IG,JPL) + FPL *
     &              (FRACREFA(IG,JPL+1)-FRACREFA(IG,JPL))

 2000          CONTINUE
c Complete calculation and do  interpolation for high CO2 cases if needed

        DO 2040 IG = 1, NG(5)
               TAUSELF = SELFFAC(LAY)* (SELFREF(INDS,IG) +
     &              SELFFRAC(LAY)  *
     &              (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
               TAUFOR = FORFAC(LAY) * (FORREF(INDF,IG) +
     &              FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &              FORREF(INDF,IG))) 
               O3M1 = KA_MO3(JMO3,INDM,IG) + FMO3*
     &              (KA_MO3(JMO3+1,INDM,IG) - 
     &              KA_MO3(JMO3,INDM,IG))
               O3M2 = KA_MO3(JMO3,INDM+1,IG) + FMO3*
     &              (KA_MO3(JMO3+1,INDM+1,IG) - 
     &              KA_MO3(JMO3,INDM+1,IG))
               ABSO3 = O3M1 + MINORFRAC(LAY) *
     &              (O3M2 - O3M1)

C                    taug(lay,ig) = taug(lay,ig)
C      &              + TAUSELF + TAUFOR
C      &              + ABSO3 * COLO3(lay)
C      &              + WX(1,LAY) * CCL4(IG)   

      taug(lay,ig)=taug(lay,ig)+ABSO3*COLO3(lay)+WX(1,LAY) * CCL4(IG)
               if( h2o_sb == 1 ) then 
                  taug(lay,ig) = taug(lay,ig) + TAUSELF
                endif
               if( h2o_for == 1 ) then 
                  taug(lay,ig) = taug(lay,ig) + TAUFOR
                endif


 2040     CONTINUE
 2500 CONTINUE


      RETURN
      END

C----------------------------------------------------------------------------

      SUBROUTINE TAUGB6
      use variables !NJE

C     BAND 6:  820-1000 cm-1 (low key - H2O, CO2)

      PARAMETER (mg=16, mxlay=203, MXMOL=38, MAXXSEC=4, NBANDS=16)
      PARAMETER (NTEMP=14,NETA=9,NREF=71)

C  Output

      COMMON /TAUGCOM/  TAUG(MXLAY,MG)
      COMMON /PLANKG/   FRACS(MXLAY,MG)                                       

C  Input

      COMMON /FEATURES/ NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)
      COMMON /PRECISE/  ONEMINUS
      COMMON /PROFILE/  NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
     &                  PZ(0:MXLAY),TZ(0:MXLAY)
      COMMON /PROFDATA/ laytrop,lay_switch_kabs,                                   
     &                  COLH2O(MXLAY),COLCO2(MXLAY),COLO3(MXLAY),  
     &                  COLN2O(MXLAY),COLCO(MXLAY),COLCH4(MXLAY),  
     &                  COLO2(MXLAY),COLBRD(MXLAY)
      COMMON /SPECIES/  COLDRY(MXLAY),WKL(MXMOL,MXLAY),WBROAD(MXLAY),
     &                  COLMOL(MXLAY),NMOL
      COMMON /XSEC/     WX(MAXXSEC,MXLAY)
      COMMON /INTFAC/   FAC00(MXLAY),FAC01(MXLAY),                            
     &                  FAC10(MXLAY),FAC11(MXLAY)                             
      COMMON /STD_REF/  PREF(NREF),PREFLOG(NREF),TREF(NREF),
     &                  CHI_STD(7,NREF)            
      COMMON /INTIND/   JP(MXLAY),JT(MXLAY),JT1(MXLAY)
      COMMON /REFRAT_ETA/ RAT_H2OCO2(MXLAY),RAT_H2OCO2_1(MXLAY),
     &                  RAT_H2OO3(MXLAY),RAT_H2OO3_1(MXLAY),
     &                  RAT_H2ON2O(MXLAY),RAT_H2ON2O_1(MXLAY),
     &                  RAT_H2OCH4(MXLAY),RAT_H2OCH4_1(MXLAY),
     &                  RAT_N2OCO2(MXLAY),RAT_N2OCO2_1(MXLAY),
     &                  RAT_O3CO2(MXLAY),RAT_O3CO2_1(MXLAY),
     &                  RAT_CO2O3(MXLAY),RAT_CO2O3_1(MXLAY),
     &                  RAT_O2_DRY(MXLAY)

      COMMON /SELF/     SELFFAC(MXLAY), SELFFRAC(MXLAY), INDSELF(MXLAY)
      COMMON /FOREIGN/  FORFAC(MXLAY), FORFRAC(MXLAY), INDFOR(MXLAY)
      COMMON /MINOR/    MINORFRAC(MXLAY), INDMINOR(MXLAY), 
     &                  SCALEMINOR(MXLAY),SCALEMINORN2(MXLAY)
      COMMON /K6/       KA(NETA,NTEMP,NREF,MG), FORREF(4,MG), 
     &                  SELFREF(10,MG), KA_MO3(9,19,MG)

      COMMON /CVRTAU/    HNAMTAU,HVRTAU

      CHARACTER*18       HNAMTAU,HVRTAU

      DIMENSION ABSA(8946,MG)
      DIMENSION FRACREFA(MG,9),CFC11ADJ(MG), CFC12(MG)
      
C Planck fraction mapping level : P = 473.4280 mb, T = 259.83 K
      DATA (FRACREFA(IG, 1),IG=1,16) /
     &1.6014E-01,1.6000E-01,1.4997E-01,1.3454E-01,1.1407E-01,9.3753E-02,
     &7.7356E-02,5.7670E-02,3.8568E-02,4.1091E-03,3.3944E-03,2.6833E-03,
     &1.9683E-03,1.2522E-03,4.6873E-04,6.4478E-05/
      DATA (FRACREFA(IG, 2),IG=1,16) /
     &1.5168E-01,1.5729E-01,1.4981E-01,1.3453E-01,1.1504E-01,9.7118E-02,
     &7.9436E-02,6.0106E-02,4.0389E-02,4.3069E-03,3.5822E-03,2.8217E-03,
     &2.0505E-03,1.2729E-03,4.8142E-04,8.4482E-05/
      DATA (FRACREFA(IG, 3),IG=1,16) /
     &1.4393E-01,1.5499E-01,1.4967E-01,1.3696E-01,1.1755E-01,9.9527E-02,
     &8.0457E-02,6.0822E-02,4.0712E-02,4.4423E-03,3.7223E-03,2.8904E-03,
     &2.1504E-03,1.5053E-03,5.8637E-04,8.4791E-05/
      DATA (FRACREFA(IG, 4),IG=1,16) /
     &1.3845E-01,1.5049E-01,1.4653E-01,1.3915E-01,1.2347E-01,1.0172E-01,
     &8.1735E-02,6.1408E-02,4.1157E-02,4.4477E-03,3.7149E-03,3.1734E-03,
     &2.3490E-03,1.5298E-03,5.8630E-04,8.4791E-05/
      DATA (FRACREFA(IG, 5),IG=1,16) /
     &1.3286E-01,1.4437E-01,1.4614E-01,1.3923E-01,1.2783E-01,1.0615E-01,
     &8.3419E-02,6.2112E-02,4.1550E-02,4.5069E-03,3.9996E-03,3.2349E-03,
     &2.4012E-03,1.5302E-03,5.8637E-04,8.4791E-05/
      DATA (FRACREFA(IG, 6),IG=1,16) /
     &1.3001E-01,1.4125E-01,1.4137E-01,1.3796E-01,1.2876E-01,1.1237E-01,
     &8.6543E-02,6.2983E-02,4.1916E-02,4.9001E-03,4.0793E-03,3.2556E-03,
     &2.4063E-03,1.5306E-03,5.8637E-04,8.4791E-05/
      DATA (FRACREFA(IG, 7),IG=1,16) /
     &1.2867E-01,1.3763E-01,1.3942E-01,1.3527E-01,1.2740E-01,1.1472E-01,
     &9.2455E-02,6.4544E-02,4.2975E-02,4.9481E-03,4.0622E-03,3.2872E-03,
     &2.4115E-03,1.5311E-03,5.8648E-04,8.4791E-05/
      DATA (FRACREFA(IG, 8),IG=1,16) /
     &1.2792E-01,1.3488E-01,1.3635E-01,1.3345E-01,1.2639E-01,1.1316E-01,
     &9.6040E-02,6.9781E-02,4.5070E-02,4.9608E-03,4.0828E-03,3.2936E-03,
     &2.4125E-03,1.5325E-03,5.8671E-04,8.4791E-05/
      DATA (FRACREFA(IG, 9),IG=1,16) /
     &1.2848E-01,1.3586E-01,1.3821E-01,1.3554E-01,1.2827E-01,1.1586E-01,
     &9.0858E-02,6.5663E-02,4.4340E-02,4.9394E-03,4.0739E-03,3.2919E-03,
     &2.4109E-03,1.5306E-03,5.8634E-04,8.4791E-05/


C Minor gas mapping level:
C     LOWER - O3,  P = 317.34 mbar, T = 245. K


      DATA CFC11ADJ/
     &     0.,  0., 36.7627,    150.757,    
     &     81.4109, 74.9112, 56.9325, 49.3226,  
     &     57.1074, 66.1202, 109.557, 89.0562,  
     &     149.865, 196.140, 258.393, 80.9923/   
      DATA CFC12/
     &     62.8368, 43.2626, 26.7549, 22.2487,
     &     23.5029, 34.8323, 26.2335, 23.2306,
     &     18.4062, 13.9534, 22.6268, 24.2604,
     &     30.0088, 26.3634, 15.8237, 57.5050/


      EQUIVALENCE (KA,ABSA),(KB,ABSB)
      REAL KA,KB
      REAL KA_MO3, MINORFRAC
      REAL O3M1, O3M2


C     Calculate reference ratio to be used in calculation of Planck
C     fraction in lower/upper atmosphere.

C     P = 473.420 mb
      REFRAT_PLANCK_A = CHI_STD(1,17)/CHI_STD(2,17)

C     P = 317.3480
      REFRAT_M_A = CHI_STD(1,19)/CHI_STD(2,19)

C     Compute the optical depth by interpolating in ln(pressure) and
C     temperature. The water vapor self-continuum and foreign continuum
C     is interpolated (in temperature) separately.  

      HVRTAU = '$Revision: 22602 $'
 
      DO 2500 LAY = 1, NLAYERS


         SPECCOMB = COLH2O(LAY) + RAT_H2OCO2(LAY)*COLCO2(LAY)
         SPECPARM = COLH2O(LAY)/SPECCOMB
         IF (SPECPARM .GE. ONEMINUS) SPECPARM = ONEMINUS
         SPECMULT = 8.*(SPECPARM)
         JS = 1 + INT(SPECMULT)
         FS = AMOD(SPECMULT,1.0)


         SPECCOMB1 = COLH2O(LAY) + RAT_H2OCO2_1(LAY)*COLCO2(LAY)
         SPECPARM1 = COLH2O(LAY)/SPECCOMB1
         IF (SPECPARM1 .GE. ONEMINUS) SPECPARM1 = ONEMINUS
         SPECMULT1 = 8.*(SPECPARM1)
         JS1 = 1 + INT(SPECMULT1)
         FS1 = AMOD(SPECMULT1,1.0)

         SPECCOMB_MO3 = COLH2O(LAY) + REFRAT_M_A*COLCO2(LAY)
         SPECPARM_MO3 = COLH2O(LAY)/SPECCOMB_MO3
         IF (SPECPARM_MO3 .GE. ONEMINUS) SPECPARM_MO3 = ONEMINUS
         SPECMULT_MO3 = 8.*SPECPARM_MO3
         JMO3 = 1 + INT(SPECMULT_MO3)
         FMO3 = AMOD(SPECMULT_MO3,1.0)

         SPECCOMB_PLANCK = COLH2O(LAY)+REFRAT_PLANCK_A*COLCO2(LAY)
         SPECPARM_PLANCK = COLH2O(LAY)/SPECCOMB_PLANCK
         IF (SPECPARM_PLANCK .GE. ONEMINUS) SPECPARM_PLANCK=ONEMINUS
         SPECMULT_PLANCK = 8.*SPECPARM_PLANCK
         JPL= 1 + INT(SPECMULT_PLANCK)
         FPL = AMOD(SPECMULT_PLANCK,1.0)

         IND0 = ((JP(LAY)-1)*NTEMP+(JT(LAY)-1))*NSPA(6) + 1
         IND1 = (JP(LAY)*NTEMP+(JT1(LAY)-1))*NSPA(6) + 1
         INDS = INDSELF(LAY)
         INDF = INDFOR(LAY)
         INDM = INDMINOR(LAY)


         IF (SPECPARM .LT. 0.125 .AND. SPECPARM1 .LT. 0.125) THEN
            P = FS - 1
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC000 = FK0*FAC00(LAY)
            FAC100 = FK1*FAC00(LAY)
            FAC200 = FK2*FAC00(LAY)
            FAC010 = FK0*FAC10(LAY)
            FAC110 = FK1*FAC10(LAY)
            FAC210 = FK2*FAC10(LAY)

            P = FS1 - 1
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC001 = FK0*FAC01(LAY)
            FAC101 = FK1*FAC01(LAY)
            FAC201 = FK2*FAC01(LAY)
            FAC011 = FK0*FAC11(LAY)
            FAC111 = FK1*FAC11(LAY)
            FAC211 = FK2*FAC11(LAY)


           
            DO 2000 IG = 1, NG(6)
               TAUSELF = SELFFAC(LAY)* (SELFREF(INDS,IG) +
     &              SELFFRAC(LAY)  *
     &              (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
               TAUFOR =  FORFAC(LAY) * (FORREF(INDF,IG) +
     &           FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &           FORREF(INDF,IG))) 
               O3M1 = KA_MO3(JMO3,INDM,IG) + FMO3*
     &              (KA_MO3(JMO3+1,INDM,IG)-KA_MO3(JMO3,INDM,IG))

               O3M2 = KA_MO3(JMO3,INDM+1,IG) + FMO3*
     &              (KA_MO3(JMO3+1,INDM+1,IG)-KA_MO3(JMO3,INDM+1,IG))
               ABSO3 = O3M1 + MINORFRAC(LAY)*(O3M2-O3M1)
               TAUG(LAY,IG) = SPECCOMB *
     &              (FAC000 * ABSA(IND0,IG) +
     &              FAC100 * ABSA(IND0+1,IG) +
     &              FAC200 * ABSA(IND0+2,IG) +
     &              FAC010 * ABSA(IND0+9,IG) +
     &              FAC110 * ABSA(IND0+10,IG) +
     &              FAC210 * ABSA(IND0+11,IG))
     &              + SPECCOMB1 *
     &              (FAC001 * ABSA(IND1,IG) + 
     &              FAC101 * ABSA(IND1+1,IG) +
     &              FAC201 * ABSA(IND1+2,IG) +
     &              FAC011 * ABSA(IND1+9,IG) +
     &              FAC111 * ABSA(IND1+10,IG) +
     &              FAC211 * ABSA(IND1+11,IG)) 
C      &              + TAUSELF + TAUFOR
     &              + ABSO3*COLO3(LAY)

              taug(lay,ig) = taug(lay,ig)
               if( h2o_sb == 1 ) then 
                  taug(lay,ig) = taug(lay,ig) + TAUSELF
                endif
               if( h2o_for == 1 ) then 
                  taug(lay,ig) = taug(lay,ig) + TAUFOR
                endif

               FRACS(LAY,IG) = FRACREFA(IG,JPL) + FPL *
     &              (FRACREFA(IG,JPL+1)-FRACREFA(IG,JPL))
 2000          CONTINUE
      ELSE IF (SPECPARM .GT. 0.875 .AND. SPECPARM1 .GT. 0.875) THEN
            P = -FS 
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC000 = FK0*FAC00(LAY)
            FAC100 = FK1*FAC00(LAY)
            FAC200 = FK2*FAC00(LAY)
            FAC010 = FK0*FAC10(LAY)
            FAC110 = FK1*FAC10(LAY)
            FAC210 = FK2*FAC10(LAY)

            P = -FS1 
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC001 = FK0*FAC01(LAY)
            FAC101 = FK1*FAC01(LAY)
            FAC201 = FK2*FAC01(LAY)
            FAC011 = FK0*FAC11(LAY)
            FAC111 = FK1*FAC11(LAY)
            FAC211 = FK2*FAC11(LAY)

            DO 2010 IG = 1, NG(6)
               TAUSELF = SELFFAC(LAY)* (SELFREF(INDS,IG) +
     &              SELFFRAC(LAY)  *
     &              (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
               TAUFOR = FORFAC(LAY) * (FORREF(INDF,IG) +
     &              FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &              FORREF(INDF,IG))) 
               O3M1 = KA_MO3(JMO3,INDM,IG) + FMO3*
     &              (KA_MO3(JMO3+1,INDM,IG)-KA_MO3(JMO3,INDM,IG))
               O3M2 = KA_MO3(JMO3,INDM+1,IG) + FMO3*
     &              (KA_MO3(JMO3+1,INDM+1,IG)-KA_MO3(JMO3,INDM+1,IG))
               ABSO3 = O3M1 + MINORFRAC(LAY)*(O3M2-O3M1)
               TAUG(LAY,IG) = SPECCOMB * 
     &              (FAC200 * ABSA(IND0-1,IG) +
     &              FAC100 * ABSA(IND0,IG) +
     &              FAC000 * ABSA(IND0+1,IG) +
     &              FAC210 * ABSA(IND0+8,IG) +
     &              FAC110 * ABSA(IND0+9,IG) +
     &              FAC010 * ABSA(IND0+10,IG))
     &              + SPECCOMB1 *
     &              (FAC201 * ABSA(IND1-1,IG) +
     &              FAC101 * ABSA(IND1,IG) +
     &              FAC001 * ABSA(IND1+1,IG) +
     &              FAC211 * ABSA(IND1+8,IG) +
     &              FAC111 * ABSA(IND1+9,IG) +
     &              FAC011 * ABSA(IND1+10,IG))
C      &              + TAUSELF+ TAUFOR
     &              + ABSO3*COLO3(LAY)

              taug(lay,ig) = taug(lay,ig)
               if( h2o_sb == 1 ) then 
                  taug(lay,ig) = taug(lay,ig) + TAUSELF
                endif
               if( h2o_for == 1 ) then 
                  taug(lay,ig) = taug(lay,ig) + TAUFOR
                endif

                FRACS(LAY,IG) = FRACREFA(IG,JPL) + FPL *
     &               (FRACREFA(IG,JPL+1)-FRACREFA(IG,JPL))
 2010           CONTINUE
       ELSE
            FAC000 = (1. - FS) * FAC00(LAY)
            FAC010 = (1. - FS) * FAC10(LAY)
            FAC100 = FS * FAC00(LAY)
            FAC110 = FS * FAC10(LAY)

            FAC001 = (1. - FS1) * FAC01(LAY)
            FAC011 = (1. - FS1) * FAC11(LAY)
            FAC101 = FS1 * FAC01(LAY)
            FAC111 = FS1 * FAC11(LAY)

            DO 2020 IG = 1, NG(6)
               TAUSELF = SELFFAC(LAY)* (SELFREF(INDS,IG) +
     &              SELFFRAC(LAY)  *
     &              (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
               TAUFOR =  FORFAC(LAY) * (FORREF(INDF,IG) +
     &              FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &           FORREF(INDF,IG))) 
               O3M1 = KA_MO3(JMO3,INDM,IG) + FMO3*
     &              (KA_MO3(JMO3+1,INDM,IG)-KA_MO3(JMO3,INDM,IG))
               O3M2 = KA_MO3(JMO3,INDM+1,IG) + FMO3*
     &              (KA_MO3(JMO3+1,INDM+1,IG)-KA_MO3(JMO3,INDM+1,IG))
               ABSO3 = O3M1 + MINORFRAC(LAY)*(O3M2-O3M1)
               TAUG(LAY,IG) = SPECCOMB * 
     &              (FAC000 * ABSA(IND0,IG) +
     &              FAC100 * ABSA(IND0+1,IG) +
     &              FAC010 * ABSA(IND0+9,IG) +
     &              FAC110 * ABSA(IND0+10,IG))
     &              + SPECCOMB1 *
     &              (FAC001 * ABSA(IND1,IG) + 
     &              FAC101 * ABSA(IND1+1,IG) +
     &              FAC011 * ABSA(IND1+9,IG) +
     &              FAC111 * ABSA(IND1+10,IG)) 
C      &              + TAUSELF + TAUFOR
     &              + ABSO3*COLO3(LAY)

              taug(lay,ig) = taug(lay,ig)
               if( h2o_sb == 1 ) then 
                  taug(lay,ig) = taug(lay,ig) + TAUSELF
                endif
               if( h2o_for == 1 ) then 
                  taug(lay,ig) = taug(lay,ig) + TAUFOR
                endif

            FRACS(LAY,IG) = FRACREFA(IG,JPL) + FPL *
     &          (FRACREFA(IG,JPL+1)-FRACREFA(IG,JPL))
 2020          CONTINUE
      ENDIF
 2500 CONTINUE

      RETURN
      END

C----------------------------------------------------------------------------

      SUBROUTINE TAUGB7
      use variables !NJE

C     BAND 7:  1000-1110 cm-1 (low key - H2O,CO2; low minor - O3)

      PARAMETER (mg=16, mxlay=203, MXMOL=38, NBANDS=16)
      PARAMETER (NTEMP=14,NETA=9,NREF=71)

C  Output

      COMMON /TAUGCOM/  TAUG(MXLAY,MG)
      COMMON /PLANKG/   FRACS(MXLAY,MG)                                       

C  Input

      COMMON /FEATURES/ NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)
      COMMON /PRECISE/  ONEMINUS
      COMMON /PROFILE/  NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
     &                  PZ(0:MXLAY),TZ(0:MXLAY)
      COMMON /PROFDATA/ laytrop,lay_switch_kabs,                                   
     &                  COLH2O(MXLAY),COLCO2(MXLAY),COLO3(MXLAY),  
     &                  COLN2O(MXLAY),COLCO(MXLAY),COLCH4(MXLAY),  
     &                  COLO2(MXLAY),COLBRD(MXLAY)
      COMMON /SPECIES/  COLDRY(MXLAY),WKL(MXMOL,MXLAY),WBROAD(MXLAY),
     &                  COLMOL(MXLAY),NMOL
      COMMON /INTFAC/   FAC00(MXLAY),FAC01(MXLAY),                            
     &                  FAC10(MXLAY),FAC11(MXLAY) 
      COMMON /STD_REF/  PREF(NREF),PREFLOG(NREF),TREF(NREF),
     &                  CHI_STD(7,NREF)
      COMMON /INTIND/   JP(MXLAY),JT(MXLAY),JT1(MXLAY)
      COMMON /REFRAT_ETA/ RAT_H2OCO2(MXLAY),RAT_H2OCO2_1(MXLAY),
     &                  RAT_H2OO3(MXLAY),RAT_H2OO3_1(MXLAY),
     &                  RAT_H2ON2O(MXLAY),RAT_H2ON2O_1(MXLAY),
     &                  RAT_H2OCH4(MXLAY),RAT_H2OCH4_1(MXLAY),
     &                  RAT_N2OCO2(MXLAY),RAT_N2OCO2_1(MXLAY),
     &                  RAT_O3CO2(MXLAY),RAT_O3CO2_1(MXLAY),
     &                  RAT_CO2O3(MXLAY),RAT_CO2O3_1(MXLAY),
     &                  RAT_O2_DRY(MXLAY)

      COMMON /SELF/     SELFFAC(MXLAY), SELFFRAC(MXLAY), INDSELF(MXLAY)
      COMMON /FOREIGN/  FORFAC(MXLAY), FORFRAC(MXLAY), INDFOR(MXLAY)
      COMMON /MINOR/    MINORFRAC(MXLAY), INDMINOR(MXLAY), 
     &                  SCALEMINOR(MXLAY),SCALEMINORN2(MXLAY)
      COMMON /K7A/       KA(NETA,NTEMP,NREF,MG),  FORREF(4,MG),
     &                  SELFREF(10,MG), KA_MO3(9,19,MG)

      COMMON /CVRTAU/    HNAMTAU,HVRTAU

      CHARACTER*18       HNAMTAU,HVRTAU

      DIMENSION ABSA(8946,MG)
      DIMENSION FRACREFA(MG,9)
      REAL KA_MO3, MINORFRAC

c Planck fraction mapping level : P = 473.92   260.

      DATA (FRACREFA(IG, 1),IG=1,16) /
     &1.5423E-01,1.5550E-01,1.4524E-01,1.3193E-01,1.1549E-01,9.8793E-02,
     &8.0830E-02,6.1175E-02,4.1503E-02,4.5203E-03,3.7376E-03,2.9497E-03,
     &2.1679E-03,1.3711E-03,5.0276E-04,7.0727E-05/
      DATA (FRACREFA(IG, 2),IG=1,16) /
     &1.6116E-01,1.5165E-01,1.4299E-01,1.3072E-01,1.1536E-01,9.9276E-02,
     &8.1112E-02,6.1231E-02,4.1270E-02,4.5021E-03,3.7210E-03,2.9407E-03,
     &2.1554E-03,1.3503E-03,4.9564E-04,7.2587E-05/
      DATA (FRACREFA(IG, 3),IG=1,16) /
     &1.6184E-01,1.5236E-01,1.4180E-01,1.2974E-01,1.1550E-01,9.9550E-02,
     &8.1419E-02,6.1245E-02,4.1456E-02,4.4696E-03,3.7074E-03,2.9383E-03,
     &2.1200E-03,1.3006E-03,4.8651E-04,7.2587E-05/
      DATA (FRACREFA(IG, 4),IG=1,16) /
     &1.6123E-01,1.5272E-01,1.4223E-01,1.2896E-01,1.1543E-01,9.9868E-02,
     &8.1569E-02,6.1412E-02,4.1598E-02,4.4763E-03,3.6974E-03,2.8916E-03,
     &2.0816E-03,1.2856E-03,4.8619E-04,7.2587E-05/
      DATA (FRACREFA(IG, 5),IG=1,16) /
     &1.5888E-01,1.5408E-01,1.4246E-01,1.2931E-01,1.1522E-01,9.9966E-02,
     &8.1942E-02,6.1560E-02,4.1629E-02,4.5020E-03,3.7020E-03,2.8168E-03,
     &2.0775E-03,1.2851E-03,4.8604E-04,7.2587E-05/
      DATA (FRACREFA(IG, 6),IG=1,16) /
     &1.5202E-01,1.5750E-01,1.4282E-01,1.3115E-01,1.1580E-01,1.0015E-01,
     &8.2260E-02,6.1709E-02,4.1575E-02,4.5385E-03,3.7410E-03,2.8106E-03,
     &2.0772E-03,1.2847E-03,4.8604E-04,7.2587E-05/
      DATA (FRACREFA(IG, 7),IG=1,16) /
     &1.4095E-01,1.5159E-01,1.5081E-01,1.3617E-01,1.1824E-01,1.0090E-01,
     &8.2696E-02,6.1848E-02,4.1740E-02,4.5632E-03,3.7729E-03,2.8086E-03,
     &2.0771E-03,1.2846E-03,4.8595E-04,7.2587E-05/
      DATA (FRACREFA(IG, 8),IG=1,16) /
     &1.3721E-01,1.4327E-01,1.4306E-01,1.3696E-01,1.2814E-01,1.0751E-01,
     &8.4263E-02,6.2518E-02,4.2002E-02,4.5657E-03,3.7729E-03,2.8033E-03,
     &2.0773E-03,1.2845E-03,4.8595E-04,7.2587E-05/
      DATA (FRACREFA(IG, 9),IG=1,16) /
     &1.4283E-01,1.4965E-01,1.4486E-01,1.3951E-01,1.2076E-01,1.0130E-01,
     &8.1593E-02,6.2420E-02,4.2013E-02,4.5600E-03,3.7722E-03,2.8134E-03,
     &2.0769E-03,1.2843E-03,4.8586E-04,7.2587E-05/


C Minor gas mapping level :
C     LOWER - O3, P = 317.43 mbar, T= 245. K

      EQUIVALENCE (KA,ABSA),(KB,ABSB)
      REAL KA,KB

C     Make sure this version of taugb7 is needed
      if(lay_switch_kabs.le.1) return

C     Calculate reference ratio to be used in calculation of Planck
C     fraction in lower atmosphere.

C     P = 473.92 mb
      REFRAT_PLANCK_A = CHI_STD(1,17)/CHI_STD(2,17)

C     P = 314.3480 mb
      REFRAT_M_A = CHI_STD(1,19)/CHI_STD(2,19)

C     Compute the optical depth by interpolating in ln(pressure), 
C     temperature, and appropriate species; the water
C     vapor self-continuum and foreign continuum is interpolated 
C     (in temperature) separately. 

      HVRTAU = '$Revision: 22602 $'

      DO 2500 LAY = 1, LAY_SWITCH_KABS

         SPECCOMB = COLH2O(LAY) + RAT_H2OCO2(LAY)*COLCO2(LAY)
         SPECPARM = COLH2O(LAY)/SPECCOMB
         IF (SPECPARM .GE. ONEMINUS) SPECPARM = ONEMINUS
         SPECMULT = 8.*(SPECPARM)
         JS = 1 + INT(SPECMULT)
         FS = AMOD(SPECMULT,1.0)

         SPECCOMB1 = COLH2O(LAY) + RAT_H2OCO2_1(LAY)*COLCO2(LAY)
         SPECPARM1 = COLH2O(LAY)/SPECCOMB1
         IF (SPECPARM1 .GE. ONEMINUS) SPECPARM1 = ONEMINUS
         SPECMULT1 = 8.*(SPECPARM1)
         JS1 = 1 + INT(SPECMULT1)
         FS1 = AMOD(SPECMULT1,1.0)

         SPECCOMB_MO3 = COLH2O(LAY) + REFRAT_M_A*COLCO2(LAY)
         SPECPARM_MO3 = COLH2O(LAY)/SPECCOMB_MO3
         IF (SPECPARM_MO3 .GE. ONEMINUS) SPECPARM_MO3 = ONEMINUS
         SPECMULT_MO3 = 8.*SPECPARM_MO3

         JMO3 = 1 + INT(SPECMULT_MO3)
         FMO3 = AMOD(SPECMULT_MO3,1.0)

         SPECCOMB_PLANCK = COLH2O(LAY)+REFRAT_PLANCK_A*COLCO2(LAY)
         SPECPARM_PLANCK = COLH2O(LAY)/SPECCOMB_PLANCK
         IF (SPECPARM_PLANCK .GE. ONEMINUS) SPECPARM_PLANCK=ONEMINUS
         SPECMULT_PLANCK = 8.*SPECPARM_PLANCK
         JPL= 1 + INT(SPECMULT_PLANCK)
         FPL = AMOD(SPECMULT_PLANCK,1.0)

         IND0 = ((JP(LAY)-1)*NTEMP+(JT(LAY)-1))*NSPA(7) + JS
         IND1 = (JP(LAY)*NTEMP+(JT1(LAY)-1))*NSPA(7) + JS1
         INDS = INDSELF(LAY)
         INDF = INDFOR(LAY)
         INDM = INDMINOR(LAY)


         IF (SPECPARM .LT. 0.125 .AND. SPECPARM1 .LT. 0.125) THEN
            P = FS - 1
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC000 = FK0*FAC00(LAY)
            FAC100 = FK1*FAC00(LAY)
            FAC200 = FK2*FAC00(LAY)
            FAC010 = FK0*FAC10(LAY)
            FAC110 = FK1*FAC10(LAY)
            FAC210 = FK2*FAC10(LAY)

            P = FS1 - 1
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC001 = FK0*FAC01(LAY)
            FAC101 = FK1*FAC01(LAY)
            FAC201 = FK2*FAC01(LAY)
            FAC011 = FK0*FAC11(LAY)
            FAC111 = FK1*FAC11(LAY)
            FAC211 = FK2*FAC11(LAY)


           
            DO 2000 IG = 1, NG(7)
               TAUSELF = SELFFAC(LAY)* (SELFREF(INDS,IG) +
     &              SELFFRAC(LAY)  *
     &              (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
               TAUFOR =  FORFAC(LAY) * (FORREF(INDF,IG) +
     &           FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &           FORREF(INDF,IG))) 
               O3M1 = KA_MO3(JMO3,INDM,IG) + FMO3*
     &              (KA_MO3(JMO3+1,INDM,IG)-KA_MO3(JMO3,INDM,IG))

               O3M2 = KA_MO3(JMO3,INDM+1,IG) + FMO3*
     &              (KA_MO3(JMO3+1,INDM+1,IG)-KA_MO3(JMO3,INDM+1,IG))
               ABSO3 = O3M1 + MINORFRAC(LAY)*(O3M2-O3M1)
               TAUG(LAY,IG) = SPECCOMB *
     &              (FAC000 * ABSA(IND0,IG) +
     &              FAC100 * ABSA(IND0+1,IG) +
     &              FAC200 * ABSA(IND0+2,IG) +
     &              FAC010 * ABSA(IND0+9,IG) +
     &              FAC110 * ABSA(IND0+10,IG) +
     &              FAC210 * ABSA(IND0+11,IG))
     &              + SPECCOMB1 *
     &              (FAC001 * ABSA(IND1,IG) + 
     &              FAC101 * ABSA(IND1+1,IG) +
     &              FAC201 * ABSA(IND1+2,IG) +
     &              FAC011 * ABSA(IND1+9,IG) +
     &              FAC111 * ABSA(IND1+10,IG) +
     &              FAC211 * ABSA(IND1+11,IG)) 
C      &              + TAUSELF + TAUFOR
     &              + ABSO3*COLO3(LAY)

              taug(lay,ig) = taug(lay,ig)
               if( h2o_sb == 1 ) then 
                  taug(lay,ig) = taug(lay,ig) + TAUSELF
                endif
               if( h2o_for == 1 ) then 
                  taug(lay,ig) = taug(lay,ig) + TAUFOR
                endif

               FRACS(LAY,IG) = FRACREFA(IG,JPL) + FPL *
     &              (FRACREFA(IG,JPL+1)-FRACREFA(IG,JPL))
 2000          CONTINUE
      ELSE IF (SPECPARM .GT. 0.875 .AND. SPECPARM1 .GT. 0.875) THEN
            P = -FS 
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC000 = FK0*FAC00(LAY)
            FAC100 = FK1*FAC00(LAY)
            FAC200 = FK2*FAC00(LAY)
            FAC010 = FK0*FAC10(LAY)
            FAC110 = FK1*FAC10(LAY)
            FAC210 = FK2*FAC10(LAY)

            P = -FS1 
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC001 = FK0*FAC01(LAY)
            FAC101 = FK1*FAC01(LAY)
            FAC201 = FK2*FAC01(LAY)
            FAC011 = FK0*FAC11(LAY)
            FAC111 = FK1*FAC11(LAY)
            FAC211 = FK2*FAC11(LAY)

            DO 2010 IG = 1, NG(7)
               TAUSELF = SELFFAC(LAY)* (SELFREF(INDS,IG) +
     &              SELFFRAC(LAY)  *
     &              (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
               TAUFOR = FORFAC(LAY) * (FORREF(INDF,IG) +
     &              FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &              FORREF(INDF,IG))) 
               O3M1 = KA_MO3(JMO3,INDM,IG) + FMO3*
     &              (KA_MO3(JMO3+1,INDM,IG)-KA_MO3(JMO3,INDM,IG))
               O3M2 = KA_MO3(JMO3,INDM+1,IG) + FMO3*
     &              (KA_MO3(JMO3+1,INDM+1,IG)-KA_MO3(JMO3,INDM+1,IG))
               ABSO3 = O3M1 + MINORFRAC(LAY)*(O3M2-O3M1)
               TAUG(LAY,IG) = SPECCOMB * 
     &              (FAC200 * ABSA(IND0-1,IG) +
     &              FAC100 * ABSA(IND0,IG) +
     &              FAC000 * ABSA(IND0+1,IG) +
     &              FAC210 * ABSA(IND0+8,IG) +
     &              FAC110 * ABSA(IND0+9,IG) +
     &              FAC010 * ABSA(IND0+10,IG))
     &              + SPECCOMB1 *
     &              (FAC201 * ABSA(IND1-1,IG) +
     &              FAC101 * ABSA(IND1,IG) +
     &              FAC001 * ABSA(IND1+1,IG) +
     &              FAC211 * ABSA(IND1+8,IG) +
     &              FAC111 * ABSA(IND1+9,IG) +
     &              FAC011 * ABSA(IND1+10,IG))
C      &              + TAUSELF+ TAUFOR
     &              + ABSO3*COLO3(LAY)

              taug(lay,ig) = taug(lay,ig)
               if( h2o_sb == 1 ) then 
                  taug(lay,ig) = taug(lay,ig) + TAUSELF
                endif
               if( h2o_for == 1 ) then 
                  taug(lay,ig) = taug(lay,ig) + TAUFOR
                endif

                FRACS(LAY,IG) = FRACREFA(IG,JPL) + FPL *
     &               (FRACREFA(IG,JPL+1)-FRACREFA(IG,JPL))
 2010           CONTINUE
       ELSE
            FAC000 = (1. - FS) * FAC00(LAY)
            FAC010 = (1. - FS) * FAC10(LAY)
            FAC100 = FS * FAC00(LAY)
            FAC110 = FS * FAC10(LAY)

            FAC001 = (1. - FS1) * FAC01(LAY)
            FAC011 = (1. - FS1) * FAC11(LAY)
            FAC101 = FS1 * FAC01(LAY)
            FAC111 = FS1 * FAC11(LAY)

            DO 2020 IG = 1, NG(7)
               TAUSELF = SELFFAC(LAY)* (SELFREF(INDS,IG) +
     &              SELFFRAC(LAY)  *
     &              (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
               TAUFOR =  FORFAC(LAY) * (FORREF(INDF,IG) +
     &              FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &           FORREF(INDF,IG))) 
               O3M1 = KA_MO3(JMO3,INDM,IG) + FMO3*
     &              (KA_MO3(JMO3+1,INDM,IG)-KA_MO3(JMO3,INDM,IG))
               O3M2 = KA_MO3(JMO3,INDM+1,IG) + FMO3*
     &              (KA_MO3(JMO3+1,INDM+1,IG)-KA_MO3(JMO3,INDM+1,IG))
               ABSO3 = O3M1 + MINORFRAC(LAY)*(O3M2-O3M1)
               TAUG(LAY,IG) = SPECCOMB * 
     &              (FAC000 * ABSA(IND0,IG) +
     &              FAC100 * ABSA(IND0+1,IG) +
     &              FAC010 * ABSA(IND0+9,IG) +
     &              FAC110 * ABSA(IND0+10,IG))
     &              + SPECCOMB1 *
     &              (FAC001 * ABSA(IND1,IG) + 
     &              FAC101 * ABSA(IND1+1,IG) +
     &              FAC011 * ABSA(IND1+9,IG) +
     &              FAC111 * ABSA(IND1+10,IG)) 
C      &              + TAUSELF + TAUFOR
     &              + ABSO3*COLO3(LAY)

              taug(lay,ig) = taug(lay,ig)
               if( h2o_sb == 1 ) then 
                  taug(lay,ig) = taug(lay,ig) + TAUSELF
                endif
               if( h2o_for == 1 ) then 
                  taug(lay,ig) = taug(lay,ig) + TAUFOR
                endif

            FRACS(LAY,IG) = FRACREFA(IG,JPL) + FPL *
     &          (FRACREFA(IG,JPL+1)-FRACREFA(IG,JPL))
 2020    CONTINUE
      ENDIF
 2500 CONTINUE


      RETURN
      END

C----------------------------------------------------------------------------

      SUBROUTINE TAUGB7B
      use variables !NJE

C     BAND 7:  1000-1110 cm-1 (low key - CO2,O3; low minor - H2O)

      PARAMETER (mg=16, mxlay=203, MXMOL=38, NBANDS=16)
      PARAMETER (NTEMP=14,NETA=9,NREF=71)

C  Output

      COMMON /TAUGCOM/  TAUG(MXLAY,MG)
      COMMON /PLANKG/   FRACS(MXLAY,MG)                                       

C  Input

      COMMON /FEATURES/ NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)
      COMMON /PRECISE/  ONEMINUS
      COMMON /PROFILE/  NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
     &                  PZ(0:MXLAY),TZ(0:MXLAY)
      COMMON /PROFDATA/ laytrop,lay_switch_kabs,                                   
     &                  COLH2O(MXLAY),COLCO2(MXLAY),COLO3(MXLAY),  
     &                  COLN2O(MXLAY),COLCO(MXLAY),COLCH4(MXLAY),  
     &                  COLO2(MXLAY),COLBRD(MXLAY)
      COMMON /SPECIES/  COLDRY(MXLAY),WKL(MXMOL,MXLAY),WBROAD(MXLAY),
     &                  COLMOL(MXLAY),NMOL
      COMMON /INTFAC/   FAC00(MXLAY),FAC01(MXLAY),                            
     &                  FAC10(MXLAY),FAC11(MXLAY) 
      COMMON /STD_REF/  PREF(NREF),PREFLOG(NREF),TREF(NREF),
     &                  CHI_STD(7,NREF)
      COMMON /INTIND/   JP(MXLAY),JT(MXLAY),JT1(MXLAY)
      COMMON /REFRAT_ETA/ RAT_H2OCO2(MXLAY),RAT_H2OCO2_1(MXLAY),
     &                  RAT_H2OO3(MXLAY),RAT_H2OO3_1(MXLAY),
     &                  RAT_H2ON2O(MXLAY),RAT_H2ON2O_1(MXLAY),
     &                  RAT_H2OCH4(MXLAY),RAT_H2OCH4_1(MXLAY),
     &                  RAT_N2OCO2(MXLAY),RAT_N2OCO2_1(MXLAY),
     &                  RAT_O3CO2(MXLAY),RAT_O3CO2_1(MXLAY),
     &                  RAT_CO2O3(MXLAY),RAT_CO2O3_1(MXLAY),
     &                  RAT_O2_DRY(MXLAY)

      COMMON /SELF/     SELFFAC(MXLAY), SELFFRAC(MXLAY), INDSELF(MXLAY)
      COMMON /FOREIGN/  FORFAC(MXLAY), FORFRAC(MXLAY), INDFOR(MXLAY)
      COMMON /MINOR/    MINORFRAC(MXLAY), INDMINOR(MXLAY), 
     &                  SCALEMINOR(MXLAY),SCALEMINORN2(MXLAY)
      COMMON /K7B/       KA(NETA,NTEMP,NREF,MG),  FORREF(4,MG),
     &                  SELFREF(10,MG), KA_MH2O(9,19,MG)

      COMMON /CVRTAU/    HNAMTAU,HVRTAU

      CHARACTER*18       HNAMTAU,HVRTAU

      DIMENSION ABSA(8946,MG)
      DIMENSION FRACREFA(MG,9)
      REAL KA_MH2O, MINORFRAC

c Planck fraction mapping level : P = 2.611   275.

      DATA (FRACREFA(IG, 1),IG=1,16) /
     &1.3989E-01,1.4465E-01,1.4424E-01,1.3510E-01,1.2173E-01,1.0489E-01,
     &8.5562E-02,6.4110E-02,4.3700E-02,4.7596E-03,3.9305E-03,3.1041E-03,
     &2.2716E-03,1.4443E-03,5.4213E-04,7.5337E-05/
      DATA (FRACREFA(IG, 2),IG=1,16) /
     &1.3993E-01,1.4481E-01,1.4422E-01,1.3504E-01,1.2169E-01,1.0486E-01,
     &8.5535E-02,6.4082E-02,4.3699E-02,4.7596E-03,3.9305E-03,3.1041E-03,
     &2.2716E-03,1.4443E-03,5.4213E-04,7.5337E-05/
      DATA (FRACREFA(IG, 3),IG=1,16) /
     &1.3994E-01,1.4489E-01,1.4422E-01,1.3504E-01,1.2164E-01,1.0485E-01,
     &8.5532E-02,6.4070E-02,4.3687E-02,4.7595E-03,3.9305E-03,3.1041E-03,
     &2.2716E-03,1.4443E-03,5.4213E-04,7.5337E-05/
      DATA (FRACREFA(IG, 4),IG=1,16) /
     &1.3995E-01,1.4498E-01,1.4420E-01,1.3503E-01,1.2164E-01,1.0481E-01,
     &8.5528E-02,6.4072E-02,4.3669E-02,4.7595E-03,3.9305E-03,3.1041E-03,
     &2.2716E-03,1.4443E-03,5.4213E-04,7.5337E-05/
      DATA (FRACREFA(IG, 5),IG=1,16) /
     &1.3997E-01,1.4509E-01,1.4417E-01,1.3500E-01,1.2164E-01,1.0479E-01,
     &8.5502E-02,6.4071E-02,4.3655E-02,4.7595E-03,3.9305E-03,3.1041E-03,
     &2.2716E-03,1.4443E-03,5.4213E-04,7.5337E-05/
      DATA (FRACREFA(IG, 6),IG=1,16) /
     &1.3999E-01,1.4521E-01,1.4415E-01,1.3495E-01,1.2162E-01,1.0478E-01,
     &8.5460E-02,6.4066E-02,4.3641E-02,4.7595E-03,3.9306E-03,3.1041E-03,
     &2.2716E-03,1.4443E-03,5.4213E-04,7.5338E-05/
      DATA (FRACREFA(IG, 7),IG=1,16) /
     &1.4003E-01,1.4538E-01,1.4414E-01,1.3489E-01,1.2157E-01,1.0477E-01,
     &8.5432E-02,6.4036E-02,4.3626E-02,4.7564E-03,3.9305E-03,3.1041E-03,
     &2.2716E-03,1.4443E-03,5.4213E-04,7.5338E-05/
      DATA (FRACREFA(IG, 8),IG=1,16) /
     &1.4012E-01,1.4571E-01,1.4415E-01,1.3478E-01,1.2146E-01,1.0469E-01,
     &8.5412E-02,6.3962E-02,4.3630E-02,4.7439E-03,3.9163E-03,3.1019E-03,
     &2.2716E-03,1.4443E-03,5.4213E-04,7.5338E-05/
      DATA (FRACREFA(IG, 9),IG=1,16) /
     &1.6537E-01,1.5336E-01,1.4139E-01,1.2813E-01,1.1463E-01,9.8793E-02,
     &8.0816E-02,6.0882E-02,4.1369E-02,4.4993E-03,3.7236E-03,2.9513E-03,
     &2.1291E-03,1.3720E-03,5.1608E-04,7.1965E-05/


C Minor gas mapping level :
C     LOWER - H2O, P = 1286.91 mbar, T= 275. K

      EQUIVALENCE (KA,ABSA),(KB,ABSB)
      REAL KA,KB

C     Make sure this version of taugb7 is needed
      if(lay_switch_kabs.ge.nlayers) return

C     Calculate reference ratio to be used in calculation of Planck
C     fraction 

C     P = 2.611
      REFRAT_PLANCK_A = CHI_STD(2,37)/CHI_STD(3,37)

C     P = 1286.91
      REFRAT_M_A = CHI_STD(2,6)/CHI_STD(3,6)

C     Compute the optical depth by interpolating in ln(pressure), 
C     temperature, and appropriate species; the water
C     vapor self-continuum and foreign continuum is interpolated 
C     (in temperature) separately. 

      HVRTAU = '$Revision: 22602 $'

      DO 2500 LAY = lay_switch_kabs+1, NLAYERS

         SPECCOMB = COLCO2(LAY) + RAT_CO2O3(LAY)*COLO3(LAY)
         SPECPARM = COLCO2(LAY)/SPECCOMB
         IF (SPECPARM .GE. ONEMINUS) SPECPARM = ONEMINUS
         SPECMULT = 8.*(SPECPARM)
         JS = 1 + INT(SPECMULT)
         FS = AMOD(SPECMULT,1.0)

         SPECCOMB1 = COLCO2(LAY) + RAT_CO2O3_1(LAY)*COLO3(LAY)
         SPECPARM1 = COLCO2(LAY)/SPECCOMB1
         IF (SPECPARM1 .GE. ONEMINUS) SPECPARM1 = ONEMINUS
         SPECMULT1 = 8.*(SPECPARM1)
         JS1 = 1 + INT(SPECMULT1)
         FS1 = AMOD(SPECMULT1,1.0)

         SPECCOMB_MH2O = COLCO2(LAY) + REFRAT_M_A*COLO3(LAY)
         SPECPARM_MH2O = COLCO2(LAY)/SPECCOMB_MH2O
         IF (SPECPARM_MH2O .GE. ONEMINUS) SPECPARM_MH2O = ONEMINUS
         SPECMULT_MH2O = 8.*SPECPARM_MH2O

         JMH2O = 1 + INT(SPECMULT_MH2O)
         FMH2O = AMOD(SPECMULT_MH2O,1.0)

         SPECCOMB_PLANCK = COLCO2(LAY)+REFRAT_PLANCK_A*COLO3(LAY)
         SPECPARM_PLANCK = COLCO2(LAY)/SPECCOMB_PLANCK
         IF (SPECPARM_PLANCK .GE. ONEMINUS) SPECPARM_PLANCK=ONEMINUS
         SPECMULT_PLANCK = 8.*SPECPARM_PLANCK
         JPL= 1 + INT(SPECMULT_PLANCK)
         FPL = AMOD(SPECMULT_PLANCK,1.0)

         IND0 = ((JP(LAY)-1)*NTEMP+(JT(LAY)-1))*NSPA(7) + JS
         IND1 = (JP(LAY)*NTEMP+(JT1(LAY)-1))*NSPA(7) + JS1
         INDS = INDSELF(LAY)
         INDF = INDFOR(LAY)
         INDM = INDMINOR(LAY)


         IF (SPECPARM .LT. 0.125 .AND. SPECPARM1 .LT. 0.125) THEN
            P = FS - 1
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC000 = FK0*FAC00(LAY)
            FAC100 = FK1*FAC00(LAY)
            FAC200 = FK2*FAC00(LAY)
            FAC010 = FK0*FAC10(LAY)
            FAC110 = FK1*FAC10(LAY)
            FAC210 = FK2*FAC10(LAY)

            P = FS1 - 1
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC001 = FK0*FAC01(LAY)
            FAC101 = FK1*FAC01(LAY)
            FAC201 = FK2*FAC01(LAY)
            FAC011 = FK0*FAC11(LAY)
            FAC111 = FK1*FAC11(LAY)
            FAC211 = FK2*FAC11(LAY)


           
            DO 2000 IG = 1, NG(7)
               TAUSELF = SELFFAC(LAY)* (SELFREF(INDS,IG) +
     &              SELFFRAC(LAY)  *
     &              (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
               TAUFOR =  FORFAC(LAY) * (FORREF(INDF,IG) +
     &           FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &           FORREF(INDF,IG))) 
               H2OM1 = ka_mh2o(Jmh2o,INDM,IG) + Fmh2o*
     &              (ka_mh2o(Jmh2o+1,INDM,IG)-ka_mh2o(Jmh2o,INDM,IG))

               H2OM2 = ka_mh2o(Jmh2o,INDM+1,IG) + Fmh2o*
     &            (ka_mh2o(Jmh2o+1,INDM+1,IG)-ka_mh2o(Jmh2o,INDM+1,IG))
               ABSH2O = H2OM1 + MINORFRAC(LAY)*(H2OM2-H2OM1)
               TAUG(LAY,IG) = SPECCOMB *
     &              (FAC000 * ABSA(IND0,IG) +
     &              FAC100 * ABSA(IND0+1,IG) +
     &              FAC200 * ABSA(IND0+2,IG) +
     &              FAC010 * ABSA(IND0+9,IG) +
     &              FAC110 * ABSA(IND0+10,IG) +
     &              FAC210 * ABSA(IND0+11,IG))
     &              + SPECCOMB1 *
     &              (FAC001 * ABSA(IND1,IG) + 
     &              FAC101 * ABSA(IND1+1,IG) +
     &              FAC201 * ABSA(IND1+2,IG) +
     &              FAC011 * ABSA(IND1+9,IG) +
     &              FAC111 * ABSA(IND1+10,IG) +
     &              FAC211 * ABSA(IND1+11,IG)) 
C      &              + TAUSELF + TAUFOR
     &              + ABSH2O*COLH2O(LAY)

              taug(lay,ig) = taug(lay,ig)
               if( h2o_sb == 1 ) then 
                  taug(lay,ig) = taug(lay,ig) + TAUSELF
                endif
               if( h2o_for == 1 ) then 
                  taug(lay,ig) = taug(lay,ig) + TAUFOR
                endif

               FRACS(LAY,IG) = FRACREFA(IG,JPL) + FPL *
     &              (FRACREFA(IG,JPL+1)-FRACREFA(IG,JPL))
 2000          CONTINUE
      ELSE IF (SPECPARM .GT. 0.875 .AND. SPECPARM1 .GT. 0.875) THEN
            P = -FS 
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC000 = FK0*FAC00(LAY)
            FAC100 = FK1*FAC00(LAY)
            FAC200 = FK2*FAC00(LAY)
            FAC010 = FK0*FAC10(LAY)
            FAC110 = FK1*FAC10(LAY)
            FAC210 = FK2*FAC10(LAY)

            P = -FS1 
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC001 = FK0*FAC01(LAY)
            FAC101 = FK1*FAC01(LAY)
            FAC201 = FK2*FAC01(LAY)
            FAC011 = FK0*FAC11(LAY)
            FAC111 = FK1*FAC11(LAY)
            FAC211 = FK2*FAC11(LAY)

            DO 2010 IG = 1, NG(7)
               TAUSELF = SELFFAC(LAY)* (SELFREF(INDS,IG) +
     &              SELFFRAC(LAY)  *
     &              (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
               TAUFOR = FORFAC(LAY) * (FORREF(INDF,IG) +
     &              FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &              FORREF(INDF,IG))) 
               H2OM1 = KA_MH2O(Jmh2o,INDM,IG) + Fmh2o*
     &              (KA_MH2O(Jmh2o+1,INDM,IG)-KA_MH2O(Jmh2o,INDM,IG))
               H2OM2 = KA_MH2O(Jmh2o,INDM+1,IG) + Fmh2o*
     &            (KA_MH2O(Jmh2o+1,INDM+1,IG)-KA_MH2O(Jmh2o,INDM+1,IG))
               ABSH2O = H2OM1 + MINORFRAC(LAY)*(H2OM2-H2OM1)
               TAUG(LAY,IG) = SPECCOMB * 
     &              (FAC200 * ABSA(IND0-1,IG) +
     &              FAC100 * ABSA(IND0,IG) +
     &              FAC000 * ABSA(IND0+1,IG) +
     &              FAC210 * ABSA(IND0+8,IG) +
     &              FAC110 * ABSA(IND0+9,IG) +
     &              FAC010 * ABSA(IND0+10,IG))
     &              + SPECCOMB1 *
     &              (FAC201 * ABSA(IND1-1,IG) +
     &              FAC101 * ABSA(IND1,IG) +
     &              FAC001 * ABSA(IND1+1,IG) +
     &              FAC211 * ABSA(IND1+8,IG) +
     &              FAC111 * ABSA(IND1+9,IG) +
     &              FAC011 * ABSA(IND1+10,IG))
C      &              + TAUSELF+ TAUFOR
     &              + ABSH2O*COLH2O(LAY)

              taug(lay,ig) = taug(lay,ig)
               if( h2o_sb == 1 ) then 
                  taug(lay,ig) = taug(lay,ig) + TAUSELF
                endif
               if( h2o_for == 1 ) then 
                  taug(lay,ig) = taug(lay,ig) + TAUFOR
                endif

                FRACS(LAY,IG) = FRACREFA(IG,JPL) + FPL * 
     &               (FRACREFA(IG,JPL+1)-FRACREFA(IG,JPL))
 2010           CONTINUE
       ELSE
            FAC000 = (1. - FS) * FAC00(LAY)
            FAC010 = (1. - FS) * FAC10(LAY)
            FAC100 = FS * FAC00(LAY)
            FAC110 = FS * FAC10(LAY)

            FAC001 = (1. - FS1) * FAC01(LAY)
            FAC011 = (1. - FS1) * FAC11(LAY)
            FAC101 = FS1 * FAC01(LAY)
            FAC111 = FS1 * FAC11(LAY)


            DO 2020 IG = 1, NG(7)
               TAUSELF = SELFFAC(LAY)* (SELFREF(INDS,IG) +
     &              SELFFRAC(LAY)  *
     &              (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
               TAUFOR =  FORFAC(LAY) * (FORREF(INDF,IG) +
     &              FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &           FORREF(INDF,IG))) 
               H2OM1 = KA_MH2O(Jmh2o,INDM,IG) + Fmh2o*
     &              (KA_MH2O(Jmh2o+1,INDM,IG)-KA_MH2O(Jmh2o,INDM,IG))
               H2OM2 = KA_MH2O(Jmh2o,INDM+1,IG) + Fmh2o*
     &            (KA_MH2O(Jmh2o+1,INDM+1,IG)-KA_MH2O(Jmh2o,INDM+1,IG))
               ABSH2O = H2OM1 + MINORFRAC(LAY)*(H2OM2-H2OM1)
               TAUG(LAY,IG) = SPECCOMB * 
     &              (FAC000 * ABSA(IND0,IG) +
     &              FAC100 * ABSA(IND0+1,IG) +
     &              FAC010 * ABSA(IND0+9,IG) +
     &              FAC110 * ABSA(IND0+10,IG))
     &              + SPECCOMB1 *
     &              (FAC001 * ABSA(IND1,IG) + 
     &              FAC101 * ABSA(IND1+1,IG) +
     &              FAC011 * ABSA(IND1+9,IG) +
     &              FAC111 * ABSA(IND1+10,IG)) 
C      &              + TAUSELF + TAUFOR
     &              + ABSH2O*COLH2O(LAY)

              taug(lay,ig) = taug(lay,ig)
               if( h2o_sb == 1 ) then 
                  taug(lay,ig) = taug(lay,ig) + TAUSELF
                endif
               if( h2o_for == 1 ) then 
                  taug(lay,ig) = taug(lay,ig) + TAUFOR
                endif

            FRACS(LAY,IG) = FRACREFA(IG,JPL) + FPL *
     &          (FRACREFA(IG,JPL+1)-FRACREFA(IG,JPL))
 2020    CONTINUE
      ENDIF
 2500 CONTINUE


      RETURN
      END

*******************************************************************************
*******************************************************************************

      SUBROUTINE TAUGB8
      use variables !NJE

C     BAND 8:  1110-1210 cm-1 (key - H2O,CH4; minor - O3,N2O)

      PARAMETER (mg=16, mxlay=203, MXMOL=38, MAXXSEC=4,NBANDS=16)
      PARAMETER (NTEMP=14,NETA=9,NREF=71)

C  Output

      COMMON /TAUGCOM/  TAUG(MXLAY,MG)
      COMMON /PLANKG/   FRACS(MXLAY,MG)                                       

C  Input

      COMMON /FEATURES/ NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)
      COMMON /PRECISE/  ONEMINUS
      COMMON /PROFILE/  NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
     &                  PZ(0:MXLAY),TZ(0:MXLAY)
      COMMON /PROFDATA/ laytrop,lay_switch_kabs,                                   
     &                  COLH2O(MXLAY),COLCO2(MXLAY),COLO3(MXLAY),
     &                  COLN2O(MXLAY),COLCO(MXLAY),COLCH4(MXLAY),
     &                  COLO2(MXLAY),COLBRD(MXLAY)
      COMMON /SPECIES/  COLDRY(MXLAY),WKL(MXMOL,MXLAY),WBROAD(MXLAY),
     &                  COLMOL(MXLAY),NMOL
      COMMON /XSEC/     WX(MAXXSEC,MXLAY)
      COMMON /INTFAC/   FAC00(MXLAY),FAC01(MXLAY),                            
     &                  FAC10(MXLAY),FAC11(MXLAY)
      COMMON /STD_REF/  PREF(NREF),PREFLOG(NREF),TREF(NREF),
     &                  CHI_STD(7,NREF)
      COMMON /INTIND/   JP(MXLAY),JT(MXLAY),JT1(MXLAY)
      COMMON /REFRAT_ETA/ RAT_H2OCO2(MXLAY),RAT_H2OCO2_1(MXLAY),
     &                  RAT_H2OO3(MXLAY),RAT_H2OO3_1(MXLAY),
     &                  RAT_H2ON2O(MXLAY),RAT_H2ON2O_1(MXLAY),
     &                  RAT_H2OCH4(MXLAY),RAT_H2OCH4_1(MXLAY),
     &                  RAT_N2OCO2(MXLAY),RAT_N2OCO2_1(MXLAY),
     &                  RAT_O3CO2(MXLAY),RAT_O3CO2_1(MXLAY),
     &                  RAT_CO2O3(MXLAY),RAT_CO2O3_1(MXLAY),
     &                  RAT_O2_DRY(MXLAY)

      COMMON /SELF/     SELFFAC(MXLAY), SELFFRAC(MXLAY), INDSELF(MXLAY)
      COMMON /FOREIGN/  FORFAC(MXLAY), FORFRAC(MXLAY), INDFOR(MXLAY)
      COMMON /MINOR/    MINORFRAC(MXLAY), INDMINOR(MXLAY),
     &                  SCALEMINOR(MXLAY),SCALEMINORN2(MXLAY)
      COMMON /K8/       KA(NETA,NTEMP,NREF,MG),  FORREF(4,MG),
     &                  SELFREF(10,MG), KA_MO3(9,19,MG), 
     &                  KA_MN2O(9,19,MG), KA_MCO2(9,19,MG)

      COMMON /CVRTAU/    HNAMTAU,HVRTAU

      CHARACTER*18       HNAMTAU,HVRTAU

      DIMENSION ABSA(8946,MG)
      DIMENSION FRACREFA(MG,9),CFC12(MG),CFC22ADJ(MG)
      REAL KA_MO3, KA_MN2O,MINORFRAC, KA, KA_MCO2
      real n2om1,n2om2

C Planck fraction mapping level : P=473.4280 mb, T = 260.0 K


      DATA (FRACREFA(IG, 1),IG=1,16) /
     &1.6843E-01,1.5594E-01,1.4392E-01,1.2952E-01,1.1542E-01,9.8378E-02,
     &7.8444E-02,5.7633E-02,3.8404E-02,4.1056E-03,3.3956E-03,2.6766E-03,
     &1.9590E-03,1.2488E-03,4.6554E-04,6.5899E-05/
      DATA (FRACREFA(IG, 2),IG=1,16) /
     &1.6601E-01,1.5726E-01,1.4261E-01,1.2846E-01,1.1327E-01,9.7139E-02,
     &7.8947E-02,5.9627E-02,4.1204E-02,4.6260E-03,3.8416E-03,2.9571E-03,
     &2.1244E-03,1.3332E-03,5.0700E-04,7.1636E-05/
      DATA (FRACREFA(IG, 3),IG=1,16) /
     &1.6546E-01,1.5779E-01,1.4228E-01,1.2739E-01,1.1250E-01,9.7975E-02,
     &7.9144E-02,6.0201E-02,4.1794E-02,4.6266E-03,3.8408E-03,2.9569E-03,
     &2.1244E-03,1.3332E-03,5.0700E-04,7.1636E-05/
      DATA (FRACREFA(IG, 4),IG=1,16) /
     &1.6503E-01,1.5811E-01,1.4153E-01,1.2693E-01,1.1288E-01,9.8103E-02,
     &7.9303E-02,6.0609E-02,4.2026E-02,4.6266E-03,3.8408E-03,2.9564E-03,
     &2.1245E-03,1.3331E-03,5.0704E-04,7.1636E-05/
      DATA (FRACREFA(IG, 5),IG=1,16) /
     &1.6462E-01,1.5806E-01,1.4114E-01,1.2665E-01,1.1344E-01,9.8031E-02,
     &7.9591E-02,6.0942E-02,4.2075E-02,4.6270E-03,3.8406E-03,2.9560E-03,
     &2.1244E-03,1.3330E-03,5.0709E-04,7.1636E-05/
      DATA (FRACREFA(IG, 6),IG=1,16) /
     &1.6375E-01,1.5845E-01,1.4068E-01,1.2673E-01,1.1374E-01,9.7974E-02,
     &7.9936E-02,6.1202E-02,4.2087E-02,4.6279E-03,3.8397E-03,2.9554E-03,
     &2.1245E-03,1.3328E-03,5.0712E-04,7.1636E-05/
      DATA (FRACREFA(IG, 7),IG=1,16) /
     &1.6275E-01,1.5901E-01,1.4052E-01,1.2656E-01,1.1412E-01,9.7810E-02,
     &8.0351E-02,6.1332E-02,4.2097E-02,4.6295E-03,3.8383E-03,2.9538E-03,
     &2.1248E-03,1.3325E-03,5.0725E-04,7.1636E-05/
      DATA (FRACREFA(IG, 8),IG=1,16) /
     &1.6102E-01,1.5999E-01,1.4042E-01,1.2667E-01,1.1426E-01,9.7926E-02,
     &8.0793E-02,6.1355E-02,4.2115E-02,4.6352E-03,3.8339E-03,2.9489E-03,
     &2.1254E-03,1.3312E-03,5.0758E-04,7.1636E-05/
      DATA (FRACREFA(IG, 9),IG=1,16) /
     &1.6252E-01,1.5926E-01,1.3835E-01,1.2702E-01,1.1481E-01,9.8114E-02,
     &8.1005E-02,6.1341E-02,4.2112E-02,4.6272E-03,3.8406E-03,2.9560E-03,
     &2.1244E-03,1.3330E-03,5.0709E-04,7.1636E-05/

C Minor gas mapping level:
C     LOWER - O3,  P = 317.348 mb, T = 245.00 K
C     LOWER - N2O, P = 706.2720 mb, should be T= 278.94 K but used t=245 for coding simplicity)
C     LOWER - CO2, P = 862.640 mb, T-245.00K
C     LOWER - CFC12,CFC11

      DATA CFC12/
     &     85.4027, 89.4696, 74.0959, 67.7480,
     &     61.2444, 59.9073, 60.8296, 63.0998,
     &     59.6110, 64.0735, 57.2622, 58.9721,
     &     43.5505, 26.1192, 32.7023, 32.8667/
C     Original CFC22 is multiplied by 1.485 to account for the 780-850 cm-1 
C     and 1290-1335 cm-1 bands.
      DATA CFC22ADJ/
     &     135.335, 89.6642, 76.2375, 65.9748,
     &     63.1164, 60.2935, 64.0299, 75.4264,
     &     51.3018, 7.07911, 5.86928, 0.398693,
     &     2.82885, 9.12751, 6.28271, 0./

      EQUIVALENCE (KA,ABSA),(KB,ABSB)

C     Calculate reference ratio to be used in calculation of Planck
C     fraction 
      REFRAT_PLANCK_A = CHI_STD(1,17)/CHI_STD(6,17)


C     Calculate reference ratio to be used in minor species absorption calculation 
      REFRAT_M_O3 = CHI_STD(1,19)/CHI_STD(6,19)
      REFRAT_M_N2O = CHI_STD(1,15)/CHI_STD(6,15)
      REFRAT_M_CO2 = CHI_STD(1,14)/CHI_STD(6,14)

C     Compute the optical depth by interpolating in ln(pressure) and 
C     temperature, and appropriate species.  Below LAYTROP, the water vapor 
C     self-continuum and foreign continuum is interpolated (in temperature) 
C     separately.

      HVRTAU = '$Revision: 22602 $'

      idump = 1

      DO 2500 LAY = 1, NLAYERS

         SPECCOMB = COLH2O(LAY) + RAT_H2OCH4(LAY)*COLCH4(LAY)
         SPECPARM = COLH2O(LAY)/SPECCOMB
         IF (SPECPARM .GE. ONEMINUS) SPECPARM = ONEMINUS
         SPECMULT = 8.*(SPECPARM)
         JS = 1 + INT(SPECMULT)
         FS = AMOD(SPECMULT,1.0)

         SPECCOMB1 = COLH2O(LAY) + RAT_H2OCH4_1(LAY)*COLCH4(LAY)
         SPECPARM1 = COLH2O(LAY)/SPECCOMB1
         IF (SPECPARM1 .GE. ONEMINUS) SPECPARM1 = ONEMINUS
         SPECMULT1 = 8.*(SPECPARM1)
         JS1 = 1 + INT(SPECMULT1)
         FS1 = AMOD(SPECMULT1,1.0)

         SPECCOMB_MO3 = COLH2O(LAY) + REFRAT_M_O3*COLCH4(LAY)
         SPECPARM_MO3 = COLH2O(LAY)/SPECCOMB_MO3
         IF (SPECPARM_MO3 .GE. ONEMINUS) SPECPARM_MO3 = ONEMINUS
         SPECMULT_MO3 = 8.*SPECPARM_MO3

         JMO3 = 1 + INT(SPECMULT_MO3)
         FMO3 = AMOD(SPECMULT_MO3,1.0)

         SPECCOMB_MN2O = COLH2O(LAY) + REFRAT_M_N2O*COLCH4(LAY)
         SPECPARM_MN2O = COLH2O(LAY)/SPECCOMB_MN2O
         IF (SPECPARM_MN2O .GE. ONEMINUS) SPECPARM_MN2O = ONEMINUS
         SPECMULT_MN2O = 8.*SPECPARM_MN2O

         JMN2O = 1 + INT(SPECMULT_MN2O)
         FMN2O = AMOD(SPECMULT_MN2O,1.0)

         SPECCOMB_MCO2 = COLH2O(LAY) + REFRAT_M_CO2*COLCH4(LAY)
         SPECPARM_MCO2 = COLH2O(LAY)/SPECCOMB_MCO2
         IF (SPECPARM_MCO2 .GE. ONEMINUS) SPECPARM_MCO2 = ONEMINUS
         SPECMULT_MCO2 = 8.*SPECPARM_MCO2

         JMCO2 = 1 + INT(SPECMULT_MCO2)
         FMCO2 = AMOD(SPECMULT_MCO2,1.0)

c  adjust CO2 column amount to reflect different P and T environment
         wk_co2 = 1.0e20*colco2(lay)/coldry(lay)
         ADJCOL2 = COLCO2(LAY)*
     &           (PAVEL(LAY)*(1.0+0.3*wk_co2)/862.64)*
     &            (245./TAVEL(LAY))**1.7

         ADJCOL2 = adjcol2/(1.0+0.16+1.0e-3*wk_co2*pavel(lay))

         SPECCOMB_PLANCK = COLH2O(LAY)+REFRAT_PLANCK_A*COLCH4(LAY)
         SPECPARM_PLANCK = COLH2O(LAY)/SPECCOMB_PLANCK
         IF (SPECPARM_PLANCK .GE. ONEMINUS) SPECPARM_PLANCK=ONEMINUS
         SPECMULT_PLANCK = 8.*SPECPARM_PLANCK
         JPL= 1 + INT(SPECMULT_PLANCK)
         FPL = AMOD(SPECMULT_PLANCK,1.0)

         IND0 = ((JP(LAY)-1)*NTEMP+(JT(LAY)-1))*NSPA(8) + JS
         IND1 = (JP(LAY)*NTEMP+(JT1(LAY)-1))*NSPA(8) + JS1
         INDS = INDSELF(LAY)
         INDF = INDFOR(LAY)
         INDM = INDMINOR(LAY)

 

         IF (SPECPARM .LT. 0.125 .AND. SPECPARM1 .LT. 0.125) THEN
            P = FS - 1
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC000 = FK0*FAC00(LAY)
            FAC100 = FK1*FAC00(LAY)
            FAC200 = FK2*FAC00(LAY)
            FAC010 = FK0*FAC10(LAY)
            FAC110 = FK1*FAC10(LAY)
            FAC210 = FK2*FAC10(LAY)

            P = FS1 - 1
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC001 = FK0*FAC01(LAY)
            FAC101 = FK1*FAC01(LAY)
            FAC201 = FK2*FAC01(LAY)
            FAC011 = FK0*FAC11(LAY)
            FAC111 = FK1*FAC11(LAY)
            FAC211 = FK2*FAC11(LAY)


           
            DO 2000 IG = 1, NG(8)
               TAUSELF = SELFFAC(LAY)* (SELFREF(INDS,IG) +
     &              SELFFRAC(LAY)  *
     &              (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
               TAUFOR =  FORFAC(LAY) * (FORREF(INDF,IG) +
     &           FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &           FORREF(INDF,IG))) 

               O3M1 = KA_MO3(JMO3,INDM,IG) + FMO3*
     &              (KA_MO3(JMO3+1,INDM,IG)-KA_MO3(JMO3,INDM,IG))

               O3M2 = KA_MO3(JMO3,INDM+1,IG) + FMO3*
     &              (KA_MO3(JMO3+1,INDM+1,IG)-KA_MO3(JMO3,INDM+1,IG))
               ABSO3 = O3M1 + MINORFRAC(LAY)*(O3M2-O3M1)

               N2OM1 = KA_MN2O(JMN2O,INDM,IG) + FMN2O*
     &              (KA_MN2O(JMN2O+1,INDM,IG)-KA_MN2O(JMN2O,INDM,IG))
               N2OM2 = KA_MN2O(JMN2O,INDM+1,IG) + FMN2O*
     &             (KA_MN2O(JMN2O+1,INDM+1,IG)-KA_MN2O(JMN2O,INDM+1,IG))
               ABSN2O = N2OM1 + MINORFRAC(LAY)*(N2OM2-N2OM1)

               CO2M1 = KA_MCO2(JMCO2,INDM,IG) + FMCO2*
     &              (KA_MCO2(JMCO2+1,INDM,IG)-KA_MCO2(JMCO2,INDM,IG))
               CO2M2 = KA_MCO2(JMCO2,INDM+1,IG) + FMCO2*
     &             (KA_MCO2(JMCO2+1,INDM+1,IG)-KA_MCO2(JMCO2,INDM+1,IG))
               ABSCO2 = CO2M1 + MINORFRAC(LAY)*(CO2M2-CO2M1)

               CO2M1 = KA_MCO2(JMCO2,9,IG) + FMCO2*
     &              (KA_MCO2(JMCO2+1,9,IG)-KA_MCO2(JMCO2,9,IG))
               CO2M2 = KA_MCO2(JMCO2,INDM+1,IG) + FMCO2*
     &             (KA_MCO2(JMCO2+1,INDM+1,IG)-KA_MCO2(JMCO2,INDM+1,IG))
               ABSCO2 = CO2M1 

               TAUG(LAY,IG) = SPECCOMB *
     &              (FAC000 * ABSA(IND0,IG) +
     &              FAC100 * ABSA(IND0+1,IG) +
     &              FAC200 * ABSA(IND0+2,IG) +
     &              FAC010 * ABSA(IND0+9,IG) +
     &              FAC110 * ABSA(IND0+10,IG) +
     &              FAC210 * ABSA(IND0+11,IG))
     &              + SPECCOMB1 *
     &              (FAC001 * ABSA(IND1,IG) + 
     &              FAC101 * ABSA(IND1+1,IG) +
     &              FAC201 * ABSA(IND1+2,IG) +
     &              FAC011 * ABSA(IND1+9,IG) +
     &              FAC111 * ABSA(IND1+10,IG) +
     &              FAC211 * ABSA(IND1+11,IG)) 
C      &              + TAUSELF + TAUFOR
     &              + ABSO3*COLO3(LAY) 
     &              + ABSN2O*COLN2O(LAY)
     &              + ABSCO2*adjcol2
     &           + WX(3,LAY) * CFC12(IG)
     &           + WX(4,LAY) * CFC22ADJ(IG)

              taug(lay,ig) = taug(lay,ig)
               if( h2o_sb == 1 ) then 
                  taug(lay,ig) = taug(lay,ig) + TAUSELF
                endif
               if( h2o_for == 1 ) then 
                  taug(lay,ig) = taug(lay,ig) + TAUFOR
                endif


               FRACS(LAY,IG) = FRACREFA(IG,JPL) + FPL *
     &              (FRACREFA(IG,JPL+1)-FRACREFA(IG,JPL))
 2000          CONTINUE
      ELSE IF (SPECPARM .GT. 0.875 .AND. SPECPARM1 .GT. 0.875) THEN
            P = -FS 
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC000 = FK0*FAC00(LAY)
            FAC100 = FK1*FAC00(LAY)
            FAC200 = FK2*FAC00(LAY)
            FAC010 = FK0*FAC10(LAY)
            FAC110 = FK1*FAC10(LAY)
            FAC210 = FK2*FAC10(LAY)

            P = -FS1 
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC001 = FK0*FAC01(LAY)
            FAC101 = FK1*FAC01(LAY)
            FAC201 = FK2*FAC01(LAY)
            FAC011 = FK0*FAC11(LAY)
            FAC111 = FK1*FAC11(LAY)
            FAC211 = FK2*FAC11(LAY)

            DO 2010 IG = 1, NG(8)
               TAUSELF = SELFFAC(LAY)* (SELFREF(INDS,IG) +
     &              SELFFRAC(LAY)  *
     &              (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
               TAUFOR = FORFAC(LAY) * (FORREF(INDF,IG) +
     &              FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &              FORREF(INDF,IG))) 

               O3M1 = KA_MO3(JMO3,INDM,IG) + FMO3*
     &              (KA_MO3(JMO3+1,INDM,IG)-KA_MO3(JMO3,INDM,IG))
               O3M2 = KA_MO3(JMO3,INDM+1,IG) + FMO3*
     &              (KA_MO3(JMO3+1,INDM+1,IG)-KA_MO3(JMO3,INDM+1,IG))
               ABSO3 = O3M1 + MINORFRAC(LAY)*(O3M2-O3M1)

               N2OM1 = KA_MN2O(JMN2O,INDM,IG) + FMN2O*
     &              (KA_MN2O(JMN2O+1,INDM,IG)-KA_MN2O(JMN2O,INDM,IG))
               N2OM2 = KA_MN2O(JMN2O,INDM+1,IG) + FMN2O*
     &             (KA_MN2O(JMN2O+1,INDM+1,IG)-KA_MN2O(JMN2O,INDM+1,IG))
               ABSN2O = N2OM1 + MINORFRAC(LAY)*(N2OM2-N2OM1)

               CO2M1 = KA_MCO2(JMCO2,9,IG) + FMCO2*
     &              (KA_MCO2(JMCO2+1,9,IG)-KA_MCO2(JMCO2,9,IG))
               CO2M2 = KA_MCO2(JMCO2,INDM+1,IG) + FMCO2*
     &             (KA_MCO2(JMCO2+1,INDM+1,IG)-KA_MCO2(JMCO2,INDM+1,IG))
               ABSCO2 = CO2M1 

               TAUG(LAY,IG) = SPECCOMB * 
     &              (FAC200 * ABSA(IND0-1,IG) +
     &              FAC100 * ABSA(IND0,IG) +
     &              FAC000 * ABSA(IND0+1,IG) +
     &              FAC210 * ABSA(IND0+8,IG) +
     &              FAC110 * ABSA(IND0+9,IG) +
     &              FAC010 * ABSA(IND0+10,IG))
     &              + SPECCOMB1 *
     &              (FAC201 * ABSA(IND1-1,IG) +
     &              FAC101 * ABSA(IND1,IG) +
     &              FAC001 * ABSA(IND1+1,IG) +
     &              FAC211 * ABSA(IND1+8,IG) +
     &              FAC111 * ABSA(IND1+9,IG) +
     &              FAC011 * ABSA(IND1+10,IG))
C      &              + TAUSELF+ TAUFOR
     &              + ABSO3*COLO3(LAY)
     &              + ABSN2O*COLN2O(LAY)
     &              + ABSCO2*adjcol2
     &           + WX(3,LAY) * CFC12(IG)
     &           + WX(4,LAY) * CFC22ADJ(IG)

              taug(lay,ig) = taug(lay,ig)
               if( h2o_sb == 1 ) then 
                  taug(lay,ig) = taug(lay,ig) + TAUSELF
                endif
               if( h2o_for == 1 ) then 
                  taug(lay,ig) = taug(lay,ig) + TAUFOR
                endif

                FRACS(LAY,IG) = FRACREFA(IG,JPL) + FPL *
     &               (FRACREFA(IG,JPL+1)-FRACREFA(IG,JPL))
 2010           CONTINUE
       ELSE
            FAC000 = (1. - FS) * FAC00(LAY)
            FAC010 = (1. - FS) * FAC10(LAY)
            FAC100 = FS * FAC00(LAY)
            FAC110 = FS * FAC10(LAY)

            FAC001 = (1. - FS1) * FAC01(LAY)
            FAC011 = (1. - FS1) * FAC11(LAY)
            FAC101 = FS1 * FAC01(LAY)
            FAC111 = FS1 * FAC11(LAY)

            DO 2020 IG = 1, NG(8)
               TAUSELF = SELFFAC(LAY)* (SELFREF(INDS,IG) +
     &              SELFFRAC(LAY)  *
     &              (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
               TAUFOR =  FORFAC(LAY) * (FORREF(INDF,IG) +
     &              FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &           FORREF(INDF,IG))) 

               O3M1 = KA_MO3(JMO3,INDM,IG) + FMO3*
     &              (KA_MO3(JMO3+1,INDM,IG)-KA_MO3(JMO3,INDM,IG))
               O3M2 = KA_MO3(JMO3,INDM+1,IG) + FMO3*
     &              (KA_MO3(JMO3+1,INDM+1,IG)-KA_MO3(JMO3,INDM+1,IG))
               ABSO3 = O3M1 + MINORFRAC(LAY)*(O3M2-O3M1)

               N2OM1 = KA_MN2O(JMN2O,INDM,IG) + FMN2O*
     &              (KA_MN2O(JMN2O+1,INDM,IG)-KA_MN2O(JMN2O,INDM,IG))
               N2OM2 = KA_MN2O(JMN2O,INDM+1,IG) + FMN2O*
     &             (KA_MN2O(JMN2O+1,INDM+1,IG)-KA_MN2O(JMN2O,INDM+1,IG))
               ABSN2O = N2OM1 + MINORFRAC(LAY)*(N2OM2-N2OM1)

               CO2M1 = KA_MCO2(JMCO2,INDM,IG) + FMCO2*
     &              (KA_MCO2(JMCO2+1,INDM,IG)-KA_MCO2(JMCO2,INDM,IG))
               CO2M2 = KA_MCO2(JMCO2,INDM+1,IG) + FMCO2*
     &             (KA_MCO2(JMCO2+1,INDM+1,IG)-KA_MCO2(JMCO2,INDM+1,IG))
               ABSCO2 = CO2M1 + MINORFRAC(LAY)*(CO2M2-CO2M1)

               CO2M1 = KA_MCO2(JMCO2,9,IG) + FMCO2*
     &              (KA_MCO2(JMCO2+1,9,IG)-KA_MCO2(JMCO2,9,IG))
               CO2M2 = KA_MCO2(JMCO2,INDM+1,IG) + FMCO2*
     &             (KA_MCO2(JMCO2+1,INDM+1,IG)-KA_MCO2(JMCO2,INDM+1,IG))
               ABSCO2 = CO2M1 

               TAUG(LAY,IG) = SPECCOMB * 
     &              (FAC000 * ABSA(IND0,IG) +
     &              FAC100 * ABSA(IND0+1,IG) +
     &              FAC010 * ABSA(IND0+9,IG) +
     &              FAC110 * ABSA(IND0+10,IG))
     &              + SPECCOMB1 *
     &              (FAC001 * ABSA(IND1,IG) + 
     &              FAC101 * ABSA(IND1+1,IG) +
     &              FAC011 * ABSA(IND1+9,IG) +
     &              FAC111 * ABSA(IND1+10,IG)) 
C      &              + TAUSELF + TAUFOR
     &              + ABSO3*COLO3(LAY)
     &              + ABSN2O*COLN2O(LAY)
     &              + ABSCO2*adjcol2
     &           + WX(3,LAY) * CFC12(IG)
     &           + WX(4,LAY) * CFC22ADJ(IG)

              taug(lay,ig) = taug(lay,ig)
               if( h2o_sb == 1 ) then 
                  taug(lay,ig) = taug(lay,ig) + TAUSELF
                endif
               if( h2o_for == 1 ) then 
                  taug(lay,ig) = taug(lay,ig) + TAUFOR
                endif

            FRACS(LAY,IG) = FRACREFA(IG,JPL) + FPL *
     &          (FRACREFA(IG,JPL+1)-FRACREFA(IG,JPL))
 2020    CONTINUE
      ENDIF
 2500 continue

      RETURN
      END

C----------------------------------------------------------------------------

      SUBROUTINE TAUGB9
      use variables !NJE

C     BAND 9:  1210-1430 cm-1 (key - H2O,CO2; minor - N2O,CH4)

      PARAMETER (mg=16, mxlay=203, MXMOL=38, NBANDS=16)
      PARAMETER (NTEMP=14,NETA=9,NREF=71)

C  Output

      COMMON /TAUGCOM/  TAUG(MXLAY,MG)
      COMMON /PLANKG/   FRACS(MXLAY,MG)                                       

C  Input

      COMMON /FEATURES/ NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)
      COMMON /PRECISE/  ONEMINUS
      COMMON /PROFILE/  NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
     &                  PZ(0:MXLAY),TZ(0:MXLAY)
      COMMON /PROFDATA/ laytrop,lay_switch_kabs,                                   
     &                  COLH2O(MXLAY),COLCO2(MXLAY),COLO3(MXLAY),  
     &                  COLN2O(MXLAY),COLCO(MXLAY),COLCH4(MXLAY),  
     &                  COLO2(MXLAY),COLBRD(MXLAY)
      COMMON /SPECIES/  COLDRY(MXLAY),WKL(MXMOL,MXLAY),WBROAD(MXLAY),
     &                  COLMOL(MXLAY),NMOL
      COMMON /STD_REF/  PREF(nref),PREFLOG(nref),TREF(nref),
     &                  CHI_std(7,nref)
      COMMON /INTFAC/   FAC00(MXLAY),FAC01(MXLAY),                            
     &                  FAC10(MXLAY),FAC11(MXLAY)                             
      COMMON /INTIND/   JP(MXLAY),JT(MXLAY),JT1(MXLAY)
      COMMON /REFRAT_ETA/ RAT_H2OCO2(MXLAY),RAT_H2OCO2_1(MXLAY),
     &                  RAT_H2OO3(MXLAY),RAT_H2OO3_1(MXLAY),
     &                  RAT_H2ON2O(MXLAY),RAT_H2ON2O_1(MXLAY),
     &                  RAT_H2OCH4(MXLAY),RAT_H2OCH4_1(MXLAY),
     &                  RAT_N2OCO2(MXLAY),RAT_N2OCO2_1(MXLAY),
     &                  RAT_O3CO2(MXLAY),RAT_O3CO2_1(MXLAY),
     &                  RAT_CO2O3(MXLAY),RAT_CO2O3_1(MXLAY),
     &                  RAT_O2_DRY(MXLAY)

      COMMON /SELF/     SELFFAC(MXLAY), SELFFRAC(MXLAY), INDSELF(MXLAY)
      COMMON /FOREIGN/  FORFAC(MXLAY), FORFRAC(MXLAY), INDFOR(MXLAY)
      COMMON /MINOR/    MINORFRAC(MXLAY), INDMINOR(MXLAY), 
     &                  SCALEMINOR(MXLAY),SCALEMINORN2(MXLAY)
      COMMON /K9/       KA(NETA,NTEMP,NREF,MG),  FORREF(4,MG),
     &                  SELFREF(10,MG), KA_MCH4(9,19,MG), 
     &                  KA_MN2O(9,19,MG)

      COMMON /CVRTAU/    HNAMTAU,HVRTAU

      CHARACTER*18       HNAMTAU,HVRTAU

      DIMENSION ABSA(8946,MG)
      DIMENSION FRACREFA(MG,9)
      REAL KA_MCH4, KA_MN2O,MINORFRAC, KA
      real n2om1,n2om2


C Planck fractions mapping level : P=862.64mb, T = 275.0 K

      DATA (FRACREFA(IG, 1),IG=1,16) /
     &1.7060E-01,1.4112E-01,1.4419E-01,1.3136E-01,1.1688E-01,9.5923E-02,
     &8.0494E-02,6.1948E-02,4.2527E-02,4.6517E-03,3.7686E-03,2.8586E-03,
     &1.9595E-03,1.2147E-03,4.4023E-04,6.1939E-05/
      DATA (FRACREFA(IG, 2),IG=1,16) /
     &1.8801E-01,1.7137E-01,1.4961E-01,1.2905E-01,1.0515E-01,8.6685E-02,
     &7.0370E-02,5.2706E-02,3.4833E-02,3.6606E-03,2.9639E-03,2.3686E-03,
     &1.6884E-03,1.0752E-03,3.9673E-04,5.4616E-05/
      DATA (FRACREFA(IG, 3),IG=1,16) /
     &1.8844E-01,1.7117E-01,1.4944E-01,1.2905E-01,1.0516E-01,8.6669E-02,
     &7.0341E-02,5.2704E-02,3.4824E-02,3.6555E-03,2.9658E-03,2.3684E-03,
     &1.6885E-03,1.0753E-03,3.9645E-04,5.4451E-05/
      DATA (FRACREFA(IG, 4),IG=1,16) /
     &1.8893E-01,1.7077E-01,1.4943E-01,1.2906E-01,1.0516E-01,8.6636E-02,
     &7.0314E-02,5.2692E-02,3.4815E-02,3.6511E-03,2.9650E-03,2.3682E-03,
     &1.6887E-03,1.0754E-03,3.9612E-04,5.4616E-05/
      DATA (FRACREFA(IG, 5),IG=1,16) /
     &1.8953E-01,1.7032E-01,1.4941E-01,1.2903E-01,1.0518E-01,8.6589E-02,
     &7.0271E-02,5.2682E-02,3.4799E-02,3.6483E-03,2.9623E-03,2.3678E-03,
     &1.6891E-03,1.0753E-03,3.9564E-04,5.4617E-05/
      DATA (FRACREFA(IG, 6),IG=1,16) /
     &1.9052E-01,1.6963E-01,1.4936E-01,1.2894E-01,1.0527E-01,8.6501E-02,
     &7.0162E-02,5.2671E-02,3.4767E-02,3.6430E-03,2.9579E-03,2.3673E-03,
     &1.6898E-03,1.0754E-03,3.9466E-04,5.4617E-05/
      DATA (FRACREFA(IG, 7),IG=1,16) /
     &1.9149E-01,1.6934E-01,1.4902E-01,1.2896E-01,1.0525E-01,8.6457E-02,
     &6.9985E-02,5.2619E-02,3.4724E-02,3.6138E-03,2.9606E-03,2.3655E-03,
     &1.6915E-03,1.0752E-03,3.9271E-04,5.4617E-05/
      DATA (FRACREFA(IG, 8),IG=1,16) /
     &1.9257E-01,1.6861E-01,1.4970E-01,1.2903E-01,1.0513E-01,8.6275E-02,
     &6.9588E-02,5.2354E-02,3.4635E-02,3.5751E-03,2.9611E-03,2.3666E-03,
     &1.6925E-03,1.0758E-03,3.8815E-04,5.4617E-05/
      DATA (FRACREFA(IG, 9),IG=1,16) /
     &1.8953E-01,1.7033E-01,1.4940E-01,1.2903E-01,1.0518E-01,8.6589E-02,
     &7.0272E-02,5.2681E-02,3.4799E-02,3.6483E-03,2.9623E-03,2.3679E-03,
     &1.6890E-03,1.0753E-03,3.9564E-04,5.4617E-05/


C Minor gas mapping level :
C     LOWER - N2O,CH4, P = 706.272 mbar, T = 275. K

      EQUIVALENCE (KA,ABSA)

      idump = 1

C     Calculate reference ratio to be used in calculation of Planck
C     fraction 

C     P = 212 mb
      REFRAT_PLANCK_A = CHI_STD(1,14)/CHI_STD(2,14)

C     Calculate reference ratio to be used in calculation of minor species absorption
C     P = 706.272 mb 
      REFRAT_M_A = CHI_STD(1,15)/CHI_STD(2,15)

C     Compute the optical depth by interpolating in ln(pressure), 
C     temperature, and appropriate species.  Below LAYTROP, the water
C     vapor self-continuum and foreign continuum is interpolated 
C     (in temperature) separately.  

      HVRTAU = '$Revision: 22602 $'

      DO 2500 LAY = 1, NLAYERS

         SPECCOMB = COLH2O(LAY) + RAT_H2OCO2(LAY)*COLCO2(LAY)
         SPECPARM = COLH2O(LAY)/SPECCOMB
         IF (SPECPARM .GE. ONEMINUS) SPECPARM = ONEMINUS
         SPECMULT = 8.*(SPECPARM)
         JS = 1 + INT(SPECMULT)
         FS = AMOD(SPECMULT,1.0)

         SPECCOMB1 = COLH2O(LAY) + RAT_H2OCO2_1(LAY)*COLCO2(LAY)
         SPECPARM1 = COLH2O(LAY)/SPECCOMB1
         IF (SPECPARM1 .GE. ONEMINUS) SPECPARM1 = ONEMINUS
         SPECMULT1 = 8.*(SPECPARM1)
         JS1 = 1 + INT(SPECMULT1)
         FS1 = AMOD(SPECMULT1,1.0)

         SPECCOMB_MN2O = COLH2O(LAY) + REFRAT_M_A*COLCO2(LAY)
         SPECPARM_MN2O = COLH2O(LAY)/SPECCOMB_MN2O
         IF (SPECPARM_MN2O .GE. ONEMINUS) SPECPARM_MN2O = ONEMINUS
         SPECMULT_MN2O = 8.*SPECPARM_MN2O
         JMN2O = 1 + INT(SPECMULT_MN2O)
         FMN2O = AMOD(SPECMULT_MN2O,1.0)

         SPECCOMB_MCH4 = COLH2O(LAY) + REFRAT_M_A*COLCO2(LAY)
         SPECPARM_MCH4 = COLH2O(LAY)/SPECCOMB_MCH4
         IF (SPECPARM_MCH4 .GE. ONEMINUS) SPECPARM_MCH4 = ONEMINUS
         SPECMULT_MCH4 = 8.*SPECPARM_MCH4
         JMCH4 = 1 + INT(SPECMULT_MCH4)
         FMCH4 = AMOD(SPECMULT_MCH4,1.0)

c     In atmospheres where the amount of N2O is too great to be considered
c     a minor species, adjust the column amount of N2O by an empirical factor 
c     to obtain the proper contribution.
         CHI_N2O = COLN2O(LAY)/(COLDRY(LAY))
         RATN2O = 1.E20*CHI_N2O/CHI_STD(4,JP(LAY)+1)
         iF (RATN2O .GT. 0.224) THEN
            ADJFAC = 0.55*(RATN2O)**0.6
            ADJCOLN2O = ADJFAC*CHI_STD(4,JP(LAY)+1)*COLDRY(LAY)*1.E-20
         ELSE
            ADJCOLN2O = COLN2O(LAY)
         ENDIF

         CHI_CH4 = COLCH4(LAY)/(COLDRY(LAY))
         RATCH4 = 1.E20*CHI_CH4/CHI_STD(6,JP(LAY)+1)
         iF (RATCH4 .GT. 0.284) THEN
            ADJFAC = 0.42*(RATCH4)**0.56
            ADJCOLCH4 = ADJFAC*CHI_STD(6,JP(LAY)+1)*COLDRY(LAY)*1.E-20
         ELSE
            ADJCOLCH4 = COLCH4(LAY)
         ENDIF

         SPECCOMB_PLANCK = COLH2O(LAY)+REFRAT_PLANCK_A*COLCO2(LAY)
         SPECPARM_PLANCK = COLH2O(LAY)/SPECCOMB_PLANCK
         IF (SPECPARM_PLANCK .GE. ONEMINUS) SPECPARM_PLANCK=ONEMINUS
         SPECMULT_PLANCK = 8.*SPECPARM_PLANCK
         JPL= 1 + INT(SPECMULT_PLANCK)
         FPL = AMOD(SPECMULT_PLANCK,1.0)

         IND0 = ((JP(LAY)-1)*NTEMP+(JT(LAY)-1))*NSPA(9) + JS
         IND1 = (JP(LAY)*NTEMP+(JT1(LAY)-1))*NSPA(9) + JS1
         INDS = INDSELF(LAY)
         INDF = INDFOR(LAY)
         INDM = INDMINOR(LAY)

         IF (SPECPARM .LT. 0.125 .AND. SPECPARM1 .LT. 0.125) THEN
            P = FS - 1
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC000 = FK0*FAC00(LAY)
            FAC100 = FK1*FAC00(LAY)
            FAC200 = FK2*FAC00(LAY)
            FAC010 = FK0*FAC10(LAY)
            FAC110 = FK1*FAC10(LAY)
            FAC210 = FK2*FAC10(LAY)

            P = FS1 - 1
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC001 = FK0*FAC01(LAY)
            FAC101 = FK1*FAC01(LAY)
            FAC201 = FK2*FAC01(LAY)
            FAC011 = FK0*FAC11(LAY)
            FAC111 = FK1*FAC11(LAY)
            FAC211 = FK2*FAC11(LAY)

            DO 2000 IG = 1, NG(9)
               TAUSELF = SELFFAC(LAY)* (SELFREF(INDS,IG) +
     &              SELFFRAC(LAY)  *
     &              (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
               TAUFOR = FORFAC(LAY) * (FORREF(INDF,IG) +
     &              FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &              FORREF(INDF,IG))) 

               N2OM1 = KA_MN2O(JMN2O,INDM,IG) + FMN2O*
     &              (KA_MN2O(JMN2O+1,INDM,IG)-
     &              KA_MN2O(JMN2O,INDM,IG))
               N2OM2 = KA_MN2O(JMN2O,INDM+1,IG) + FMN2O*
     &              (KA_MN2O(JMN2O+1,INDM+1,IG)-
     &              KA_MN2O(JMN2O,INDM+1,IG))
               ABSN2O = N2OM1 + MINORFRAC(LAY)
     &              * (N2OM2 - N2OM1)
               CH4M1 = KA_MCH4(JMCH4,INDM,IG) + FMCH4*
     &              (KA_MCH4(JMCH4+1,INDM,IG)-
     &              KA_MCH4(JMCH4,INDM,IG))
               CH4M2 = KA_MCH4(JMCH4,INDM+1,IG) + FMCH4*
     &              (KA_MCH4(JMCH4+1,INDM+1,IG)-
     &              KA_MCH4(JMCH4,INDM+1,IG))
               ABSCH4 = CH4M1 + MINORFRAC(LAY)
     &              * (CH4M2 - CH4M1)

               TAUG(LAY,IG) = SPECCOMB *
     &          (FAC000 * ABSA(IND0,IG) +
     &          FAC100 * ABSA(IND0+1,IG) +
     &          FAC200 * ABSA(IND0+2,IG) +
     &          FAC010 * ABSA(IND0+9,IG) +
     &          FAC110 * ABSA(IND0+10,IG) +
     &          FAC210 * ABSA(IND0+11,IG))
     &          + SPECCOMB1 *
     &          (FAC001 * ABSA(IND1,IG) + 
     &          FAC101 * ABSA(IND1+1,IG) +
     &          FAC201 * ABSA(IND1+2,IG) +
     &          FAC011 * ABSA(IND1+9,IG) +
     &          FAC111 * ABSA(IND1+10,IG) +
     &          FAC211 * ABSA(IND1+11,IG)) 
C      &          + TAUSELF + TAUFOR
     &          + ADJCOLN2O*ABSN2O            
     &          + ADJCOLCH4*ABSCH4  

                taug(lay,ig) = taug(lay,ig)
               if( h2o_sb == 1 ) then 
                  taug(lay,ig) = taug(lay,ig) + TAUSELF
                endif
               if( h2o_for == 1 ) then 
                  taug(lay,ig) = taug(lay,ig) + TAUFOR
                endif          

               FRACS(LAY,IG) = FRACREFA(IG,JPL) + FPL *
     &              (FRACREFA(IG,JPL+1)-FRACREFA(IG,JPL))
 2000          CONTINUE
      ELSE IF (SPECPARM .GT. 0.875 .AND. SPECPARM1 .GT. 0.875) THEN
            P = -FS 
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC000 = FK0*FAC00(LAY)
            FAC100 = FK1*FAC00(LAY)
            FAC200 = FK2*FAC00(LAY)
            FAC010 = FK0*FAC10(LAY)
            FAC110 = FK1*FAC10(LAY)
            FAC210 = FK2*FAC10(LAY)

            P = -FS1 
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC001 = FK0*FAC01(LAY)
            FAC101 = FK1*FAC01(LAY)
            FAC201 = FK2*FAC01(LAY)
            FAC011 = FK0*FAC11(LAY)
            FAC111 = FK1*FAC11(LAY)
            FAC211 = FK2*FAC11(LAY)

            DO 2010 IG = 1, NG(9)
               TAUSELF = SELFFAC(LAY)* (SELFREF(INDS,IG) +
     &              SELFFRAC(LAY)  *
     &              (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
               TAUFOR = FORFAC(LAY) * (FORREF(INDF,IG) +
     &              FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &              FORREF(INDF,IG))) 
               N2OM1 = KA_MN2O(JMN2O,INDM,IG) + FMN2O*
     &              (KA_MN2O(JMN2O+1,INDM,IG)-
     &              KA_MN2O(JMN2O,INDM,IG))
               N2OM2 = KA_MN2O(JMN2O,INDM+1,IG) + FMN2O*
     &              (KA_MN2O(JMN2O+1,INDM+1,IG)-
     &              KA_MN2O(JMN2O,INDM+1,IG))
               ABSN2O = N2OM1 + MINORFRAC(LAY)
     &              * (N2OM2 - N2OM1)
               CH4M1 = KA_MCH4(JMCH4,INDM,IG) + FMCH4*
     &              (KA_MCH4(JMCH4+1,INDM,IG)-
     &              KA_MCH4(JMCH4,INDM,IG))
               CH4M2 = KA_MCH4(JMCH4,INDM+1,IG) + FMCH4*
     &              (KA_MCH4(JMCH4+1,INDM+1,IG)-
     &              KA_MCH4(JMCH4,INDM+1,IG))
               ABSCH4 = CH4M1 + MINORFRAC(LAY)
     &              * (CH4M2 - CH4M1)
               TAUG(LAY,IG) = SPECCOMB * 
     &              (FAC200 * ABSA(IND0-1,IG) +
     &              FAC100 * ABSA(IND0,IG) +
     &              FAC000 * ABSA(IND0+1,IG) +
     &              FAC210 * ABSA(IND0+8,IG) +
     &              FAC110 * ABSA(IND0+9,IG) +
     &              FAC010 * ABSA(IND0+10,IG))
     &              + SPECCOMB1 *
     &              (FAC201 * ABSA(IND1-1,IG) +
     &              FAC101 * ABSA(IND1,IG) +
     &              FAC001 * ABSA(IND1+1,IG) +
     &              FAC211 * ABSA(IND1+8,IG) +
     &              FAC111 * ABSA(IND1+9,IG) +
     &              FAC011 * ABSA(IND1+10,IG))
C      &              + TAUSELF + TAUFOR
     &              + ADJCOLN2O*ABSN2O            
     &              + ADJCOLCH4*ABSCH4  


                taug(lay,ig) = taug(lay,ig)
               if( h2o_sb == 1 ) then 
                  taug(lay,ig) = taug(lay,ig) + TAUSELF
                endif
               if( h2o_for == 1 ) then 
                  taug(lay,ig) = taug(lay,ig) + TAUFOR
                endif          


               FRACS(LAY,IG) = FRACREFA(IG,JPL) + FPL *
     &              (FRACREFA(IG,JPL+1)-FRACREFA(IG,JPL))
 2010           CONTINUE
       ELSE
            FAC000 = (1. - FS) * FAC00(LAY)
            FAC010 = (1. - FS) * FAC10(LAY)
            FAC100 = FS * FAC00(LAY)
            FAC110 = FS * FAC10(LAY)

            FAC001 = (1. - FS1) * FAC01(LAY)
            FAC011 = (1. - FS1) * FAC11(LAY)
            FAC101 = FS1 * FAC01(LAY)
            FAC111 = FS1 * FAC11(LAY)

            DO 2020 IG = 1, NG(9)
               TAUSELF = SELFFAC(LAY)* (SELFREF(INDS,IG) +
     &              SELFFRAC(LAY)  *
     &              (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
               TAUFOR = FORFAC(LAY) * (FORREF(INDF,IG) +
     &              FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &              FORREF(INDF,IG))) 
               N2OM1 = KA_MN2O(JMN2O,INDM,IG) + FMN2O*
     &              (KA_MN2O(JMN2O+1,INDM,IG)-
     &              KA_MN2O(JMN2O,INDM,IG))
               N2OM2 = KA_MN2O(JMN2O,INDM+1,IG) + FMN2O*
     &              (KA_MN2O(JMN2O+1,INDM+1,IG)-
     &              KA_MN2O(JMN2O,INDM+1,IG))
               ABSN2O = N2OM1 + MINORFRAC(LAY)
     &              * (N2OM2 - N2OM1)
               CH4M1 = KA_MCH4(JMCH4,INDM,IG) + FMCH4*
     &              (KA_MCH4(JMCH4+1,INDM,IG)-
     &              KA_MCH4(JMCH4,INDM,IG))
               CH4M2 = KA_MCH4(JMCH4,INDM+1,IG) + FMCH4*
     &              (KA_MCH4(JMCH4+1,INDM+1,IG)-
     &              KA_MCH4(JMCH4,INDM+1,IG))
               ABSCH4 = CH4M1 + MINORFRAC(LAY)
     &              * (CH4M2 - CH4M1)
               TAUG(LAY,IG) = SPECCOMB * 
     &              (FAC000 * ABSA(IND0,IG) +
     &              FAC100 * ABSA(IND0+1,IG) +
     &              FAC010 * ABSA(IND0+9,IG) +
     &              FAC110 * ABSA(IND0+10,IG))
     &              + SPECCOMB1 *
     &              (FAC001 * ABSA(IND1,IG) + 
     &              FAC101 * ABSA(IND1+1,IG) +
     &              FAC011 * ABSA(IND1+9,IG) +
     &              FAC111 * ABSA(IND1+10,IG)) 
C      &              + TAUSELF + TAUFOR
     &              + ADJCOLN2O*ABSN2O            
     &              + ADJCOLCH4*ABSCH4  


                taug(lay,ig) = taug(lay,ig)
               if( h2o_sb == 1 ) then 
                  taug(lay,ig) = taug(lay,ig) + TAUSELF
                endif
               if( h2o_for == 1 ) then 
                  taug(lay,ig) = taug(lay,ig) + TAUFOR
                endif          


               FRACS(LAY,IG) = FRACREFA(IG,JPL) + FPL *
     &              (FRACREFA(IG,JPL+1)-FRACREFA(IG,JPL))
 2020          CONTINUE
      ENDIF
 2500 CONTINUE

      RETURN
      END

*******************************************************************************


      SUBROUTINE TAUGB10
      use variables !NJE

C     BAND 10:  1390-1480 cm-1 (low key - H2O, CO2; minor - O2)

      PARAMETER (mg=16, mxlay=203, MXMOL=38, NBANDS=16)
      PARAMETER (NTEMP=14,NETA=9,NREF=71)

C  Output

      COMMON /TAUGCOM/  TAUG(MXLAY,MG)
      COMMON /PLANKG/   FRACS(MXLAY,MG)                                       

C  Input

      COMMON /FEATURES/ NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)
      COMMON /PRECISE/  ONEMINUS
      COMMON /PROFILE/  NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
     &                  PZ(0:MXLAY),TZ(0:MXLAY)
      COMMON /PROFDATA/ laytrop,lay_switch_kabs,                                   
     &                  COLH2O(MXLAY),COLCO2(MXLAY),COLO3(MXLAY),  
     &                  COLN2O(MXLAY),COLCO(MXLAY),COLCH4(MXLAY),  
     &                  COLO2(MXLAY),COLBRD(MXLAY)
      COMMON /SPECIES/  COLDRY(MXLAY),WKL(MXMOL,MXLAY),WBROAD(MXLAY),
     &                  COLMOL(MXLAY),NMOL
      COMMON /STD_REF/  PREF(nref),PREFLOG(nref),TREF(nref),
     &                  CHI_std(7,nref)
      COMMON /INTFAC/   FAC00(MXLAY),FAC01(MXLAY),                            
     &                  FAC10(MXLAY),FAC11(MXLAY)                             
      COMMON /INTIND/   JP(MXLAY),JT(MXLAY),JT1(MXLAY)
      COMMON /REFRAT_ETA/ RAT_H2OCO2(MXLAY),RAT_H2OCO2_1(MXLAY),
     &                  RAT_H2OO3(MXLAY),RAT_H2OO3_1(MXLAY),
     &                  RAT_H2ON2O(MXLAY),RAT_H2ON2O_1(MXLAY),
     &                  RAT_H2OCH4(MXLAY),RAT_H2OCH4_1(MXLAY),
     &                  RAT_N2OCO2(MXLAY),RAT_N2OCO2_1(MXLAY),
     &                  RAT_O3CO2(MXLAY),RAT_O3CO2_1(MXLAY),
     &                  RAT_CO2O3(MXLAY),RAT_CO2O3_1(MXLAY),
     &                  RAT_O2_DRY(MXLAY)

      COMMON /SELF/     SELFFAC(MXLAY), SELFFRAC(MXLAY), INDSELF(MXLAY)
      COMMON /FOREIGN/  FORFAC(MXLAY), FORFRAC(MXLAY), INDFOR(MXLAY)
      COMMON /MINOR/    MINORFRAC(MXLAY), INDMINOR(MXLAY), 
     &                  SCALEMINOR(MXLAY),SCALEMINORN2(MXLAY)
      COMMON /K10/      KA(NETA,NTEMP,NREF,MG),  FORREF(4,MG),
     &                  SELFREF(10,MG), KA_MO2(9,19,MG)

      COMMON /CVRTAU/    HNAMTAU,HVRTAU

      CHARACTER*18       HNAMTAU,HVRTAU

      DIMENSION ABSA(8946,MG)
      DIMENSION FRACREFA(MG,9)
      REAL KA_MO2, MINORFRAC, KA


C Planck fractions mapping level : P=1053.64mb, T = 290.0 K


      DATA (FRACREFA(IG, 1),IG=1,16) /
     &1.0786E-01,1.0050E-01,1.0045E-01,1.2434E-01,1.3090E-01,1.2888E-01,
     &1.1690E-01,9.5425E-02,6.8555E-02,7.6594E-03,6.3661E-03,5.0512E-03,
     &3.7164E-03,2.3718E-03,8.9499E-04,1.2598E-04/
      DATA (FRACREFA(IG, 2),IG=1,16) /
     &1.5150E-01,1.4872E-01,1.3894E-01,1.3067E-01,1.1684E-01,1.0302E-01,
     &8.6335E-02,6.4241E-02,4.3690E-02,4.8419E-03,3.9535E-03,3.0970E-03,
     &2.3706E-03,1.3405E-03,3.9988E-04,5.6349E-05/
      DATA (FRACREFA(IG, 3),IG=1,16) /
     &1.5163E-01,1.4862E-01,1.3896E-01,1.3061E-01,1.1682E-01,1.0303E-01,
     &8.6328E-02,6.4274E-02,4.3682E-02,4.8530E-03,3.9469E-03,3.0924E-03,
     &2.3740E-03,1.3397E-03,4.0000E-04,5.6377E-05/
      DATA (FRACREFA(IG, 4),IG=1,16) /
     &1.5183E-01,1.4845E-01,1.3906E-01,1.3048E-01,1.1677E-01,1.0304E-01,
     &8.6288E-02,6.4319E-02,4.3707E-02,4.8769E-03,3.9158E-03,3.0895E-03,
     &2.3754E-03,1.3389E-03,4.0024E-04,5.6424E-05/
      DATA (FRACREFA(IG, 5),IG=1,16) /
     &1.5218E-01,1.4814E-01,1.3916E-01,1.3032E-01,1.1677E-01,1.0299E-01,
     &8.6307E-02,6.4373E-02,4.3721E-02,4.8840E-03,3.8980E-03,3.0941E-03,
     &2.3701E-03,1.3378E-03,4.0056E-04,5.6375E-05/
      DATA (FRACREFA(IG, 6),IG=1,16) /
     &1.5254E-01,1.4795E-01,1.3910E-01,1.3016E-01,1.1679E-01,1.0301E-01,
     &8.6286E-02,6.4383E-02,4.3738E-02,4.9170E-03,3.8715E-03,3.1161E-03,
     &2.3463E-03,1.3370E-03,4.0120E-04,5.6555E-05/
      DATA (FRACREFA(IG, 7),IG=1,16) /
     &1.5316E-01,1.4791E-01,1.3869E-01,1.2984E-01,1.1686E-01,1.0295E-01,
     &8.6340E-02,6.4546E-02,4.3650E-02,4.9716E-03,3.8079E-03,3.1799E-03,
     &2.2865E-03,1.3426E-03,4.0257E-04,5.6803E-05/
      DATA (FRACREFA(IG, 8),IG=1,16) /
     &1.5369E-01,1.4809E-01,1.3786E-01,1.2989E-01,1.1661E-01,1.0323E-01,
     &8.6198E-02,6.4991E-02,4.3641E-02,4.7036E-03,3.8506E-03,3.2522E-03,
     &2.2311E-03,1.2977E-03,4.0860E-04,5.7625E-05/
      DATA (FRACREFA(IG, 9),IG=1,16) /
     &1.5218E-01,1.4814E-01,1.3916E-01,1.3032E-01,1.1677E-01,1.0298E-01,
     &8.6311E-02,6.4370E-02,4.3722E-02,4.8839E-03,3.8982E-03,3.0941E-03,
     &2.3700E-03,1.3378E-03,4.0058E-04,5.6473E-05/

C Minor gas mapping level :
C     O2 - 706.272 mbar - 275K

      EQUIVALENCE (KA,ABSA)

      idump = 1

C     Calculate reference ratio to be used in calculation of Planck
C     fraction 

C     P = 1053.63 mb
      REFRAT_PLANCK_A = CHI_STD(1,13)/CHI_STD(2,13)

C     Calculate reference ratio to be used in calculation of minor species absorption
C     P = 706.272 mb 
      REFRAT_M_A = CHI_STD(1,15)/CHI_STD(2,15)

C     Compute the optical depth by interpolating in ln(pressure), 
C     temperature, and appropriate species.  Below LAYTROP, the water
C     vapor self-continuum and foreign continuum is interpolated 
C     (in temperature) separately.  

      HVRTAU = '$Revision: 22602 $'

      DO 2500 LAY = 1, NLAYERS

         SPECCOMB = COLH2O(LAY) + RAT_H2OCO2(LAY)*COLCO2(LAY)
         SPECPARM = COLH2O(LAY)/SPECCOMB
         IF (SPECPARM .GE. ONEMINUS) SPECPARM = ONEMINUS
         SPECMULT = 8.*(SPECPARM)
         JS = 1 + INT(SPECMULT)
         FS = AMOD(SPECMULT,1.0)

         SPECCOMB1 = COLH2O(LAY) + RAT_H2OCO2_1(LAY)*COLCO2(LAY)
         SPECPARM1 = COLH2O(LAY)/SPECCOMB1
         IF (SPECPARM1 .GE. ONEMINUS) SPECPARM1 = ONEMINUS
         SPECMULT1 = 8.*(SPECPARM1)
         JS1 = 1 + INT(SPECMULT1)
         FS1 = AMOD(SPECMULT1,1.0)

         SPECCOMB_MO2 = COLH2O(LAY) + REFRAT_M_A*COLCO2(LAY)
         SPECPARM_MO2 = COLH2O(LAY)/SPECCOMB_MO2
         IF (SPECPARM_MO2 .GE. ONEMINUS) SPECPARM_MO2 = ONEMINUS
         SPECMULT_MO2 = 8.*SPECPARM_MO2
         JMO2 = 1 + INT(SPECMULT_MO2)
         FMO2 = AMOD(SPECMULT_MO2,1.0)

c     In atmospheres where the amount of N2O is too great to be considered
c     a minor species, adjust the column amount of N2O by an empirical factor 
c     to obtain the proper contribution.
c        CHI_N2O = COLN2O(LAY)/(COLDRY(LAY))
c        RATN2O = 1.E20*CHI_N2O/CHI_STD(4,JP(LAY)+1)
c        iF (RATN2O .GT. 0.224) THEN
c           ADJFAC = 0.55*(RATN2O)**0.6
c           ADJCOLN2O = ADJFAC*CHI_STD(4,JP(LAY)+1)*COLDRY(LAY)*1.E-20
c        ELSE
            ADJCOLO2 = COLO2(LAY)
c        ENDIF


         SPECCOMB_PLANCK = COLH2O(LAY)+REFRAT_PLANCK_A*COLCO2(LAY)
         SPECPARM_PLANCK = COLH2O(LAY)/SPECCOMB_PLANCK
         IF (SPECPARM_PLANCK .GE. ONEMINUS) SPECPARM_PLANCK=ONEMINUS
         SPECMULT_PLANCK = 8.*SPECPARM_PLANCK
         JPL= 1 + INT(SPECMULT_PLANCK)
         FPL = AMOD(SPECMULT_PLANCK,1.0)

         IND0 = ((JP(LAY)-1)*NTEMP+(JT(LAY)-1))*NSPA(10) + JS
         IND1 = (JP(LAY)*NTEMP+(JT1(LAY)-1))*NSPA(10) + JS1
         INDS = INDSELF(LAY)
         INDF = INDFOR(LAY)
         INDM = INDMINOR(LAY)

         IF (SPECPARM .LT. 0.125 .AND. SPECPARM1 .LT. 0.125) THEN
            P = FS - 1
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC000 = FK0*FAC00(LAY)
            FAC100 = FK1*FAC00(LAY)
            FAC200 = FK2*FAC00(LAY)
            FAC010 = FK0*FAC10(LAY)
            FAC110 = FK1*FAC10(LAY)
            FAC210 = FK2*FAC10(LAY)

            P = FS1 - 1
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC001 = FK0*FAC01(LAY)
            FAC101 = FK1*FAC01(LAY)
            FAC201 = FK2*FAC01(LAY)
            FAC011 = FK0*FAC11(LAY)
            FAC111 = FK1*FAC11(LAY)
            FAC211 = FK2*FAC11(LAY)

            DO 2000 IG = 1, NG(10)
               TAUSELF = SELFFAC(LAY)* (SELFREF(INDS,IG) +
     &              SELFFRAC(LAY)  *
     &              (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
               TAUFOR = FORFAC(LAY) * (FORREF(INDF,IG) +
     &              FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &              FORREF(INDF,IG))) 

               O2M1 = KA_MO2(JMO2,INDM,IG) + FMO2*
     &              (KA_MO2(JMO2+1,INDM,IG)-
     &              KA_MO2(JMO2,INDM,IG))
               O2M2 = KA_MO2(JMO2,INDM+1,IG) + FMO2*
     &              (KA_MO2(JMO2+1,INDM+1,IG)-
     &              KA_MO2(JMO2,INDM+1,IG))
               ABSO2 = O2M1 + MINORFRAC(LAY)
     &              * (O2M2 - O2M1)

               TAUG(LAY,IG) = SPECCOMB *
     &          (FAC000 * ABSA(IND0,IG) +
     &          FAC100 * ABSA(IND0+1,IG) +
     &          FAC200 * ABSA(IND0+2,IG) +
     &          FAC010 * ABSA(IND0+9,IG) +
     &          FAC110 * ABSA(IND0+10,IG) +
     &          FAC210 * ABSA(IND0+11,IG))
     &          + SPECCOMB1 *
     &          (FAC001 * ABSA(IND1,IG) + 
     &          FAC101 * ABSA(IND1+1,IG) +
     &          FAC201 * ABSA(IND1+2,IG) +
     &          FAC011 * ABSA(IND1+9,IG) +
     &          FAC111 * ABSA(IND1+10,IG) +
     &          FAC211 * ABSA(IND1+11,IG)) 
C      &          + TAUSELF + TAUFOR
     &          + ADJCOLO2*ABSO2    


                taug(lay,ig) = taug(lay,ig)
               if( h2o_sb == 1 ) then 
                  taug(lay,ig) = taug(lay,ig) + TAUSELF
                endif
               if( h2o_for == 1 ) then 
                  taug(lay,ig) = taug(lay,ig) + TAUFOR
                endif        

               FRACS(LAY,IG) = FRACREFA(IG,JPL) + FPL *
     &              (FRACREFA(IG,JPL+1)-FRACREFA(IG,JPL))
 2000          CONTINUE
      ELSE IF (SPECPARM .GT. 0.875 .AND. SPECPARM1 .GT. 0.875) THEN
            P = -FS 
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC000 = FK0*FAC00(LAY)
            FAC100 = FK1*FAC00(LAY)
            FAC200 = FK2*FAC00(LAY)
            FAC010 = FK0*FAC10(LAY)
            FAC110 = FK1*FAC10(LAY)
            FAC210 = FK2*FAC10(LAY)

            P = -FS1 
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC001 = FK0*FAC01(LAY)
            FAC101 = FK1*FAC01(LAY)
            FAC201 = FK2*FAC01(LAY)
            FAC011 = FK0*FAC11(LAY)
            FAC111 = FK1*FAC11(LAY)
            FAC211 = FK2*FAC11(LAY)

            DO 2010 IG = 1, NG(10)
               TAUSELF = SELFFAC(LAY)* (SELFREF(INDS,IG) +
     &              SELFFRAC(LAY)  *
     &              (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
               TAUFOR = FORFAC(LAY) * (FORREF(INDF,IG) +
     &              FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &              FORREF(INDF,IG))) 

               O2M1 = KA_MO2(JMO2,INDM,IG) + FMO2*
     &              (KA_MO2(JMO2+1,INDM,IG)-
     &              KA_MO2(JMO2,INDM,IG))
               O2M2 = KA_MO2(JMO2,INDM+1,IG) + FMO2*
     &              (KA_MO2(JMO2+1,INDM+1,IG)-
     &              KA_MO2(JMO2,INDM+1,IG))
               ABSO2 = O2M1 + MINORFRAC(LAY)
     &              * (O2M2 - O2M1)

               TAUG(LAY,IG) = SPECCOMB * 
     &              (FAC200 * ABSA(IND0-1,IG) +
     &              FAC100 * ABSA(IND0,IG) +
     &              FAC000 * ABSA(IND0+1,IG) +
     &              FAC210 * ABSA(IND0+8,IG) +
     &              FAC110 * ABSA(IND0+9,IG) +
     &              FAC010 * ABSA(IND0+10,IG))
     &              + SPECCOMB1 *
     &              (FAC201 * ABSA(IND1-1,IG) +
     &              FAC101 * ABSA(IND1,IG) +
     &              FAC001 * ABSA(IND1+1,IG) +
     &              FAC211 * ABSA(IND1+8,IG) +
     &              FAC111 * ABSA(IND1+9,IG) +
     &              FAC011 * ABSA(IND1+10,IG))
C      &              + TAUSELF + TAUFOR
     &              + ADJCOLO2*ABSO2            


                taug(lay,ig) = taug(lay,ig)
               if( h2o_sb == 1 ) then 
                  taug(lay,ig) = taug(lay,ig) + TAUSELF
                endif
               if( h2o_for == 1 ) then 
                  taug(lay,ig) = taug(lay,ig) + TAUFOR
                endif

c


               FRACS(LAY,IG) = FRACREFA(IG,JPL) + FPL *
     &              (FRACREFA(IG,JPL+1)-FRACREFA(IG,JPL))
 2010           CONTINUE
       ELSE
            FAC000 = (1. - FS) * FAC00(LAY)
            FAC010 = (1. - FS) * FAC10(LAY)
            FAC100 = FS * FAC00(LAY)
            FAC110 = FS * FAC10(LAY)

            FAC001 = (1. - FS1) * FAC01(LAY)
            FAC011 = (1. - FS1) * FAC11(LAY)
            FAC101 = FS1 * FAC01(LAY)
            FAC111 = FS1 * FAC11(LAY)

            DO 2020 IG = 1, NG(10)
               TAUSELF = SELFFAC(LAY)* (SELFREF(INDS,IG) +
     &              SELFFRAC(LAY)  *
     &              (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
               TAUFOR = FORFAC(LAY) * (FORREF(INDF,IG) +
     &              FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &              FORREF(INDF,IG))) 

               O2M1 = KA_MO2(JMO2,INDM,IG) + FMO2*
     &              (KA_MO2(JMO2+1,INDM,IG)-
     &              KA_MO2(JMO2,INDM,IG))
               O2M2 = KA_MO2(JMO2,INDM+1,IG) + FMO2*
     &              (KA_MO2(JMO2+1,INDM+1,IG)-
     &              KA_MO2(JMO2,INDM+1,IG))
               ABSO2 = O2M1 + MINORFRAC(LAY)
     &              * (O2M2 - O2M1)

               TAUG(LAY,IG) = SPECCOMB * 
     &              (FAC000 * ABSA(IND0,IG) +
     &              FAC100 * ABSA(IND0+1,IG) +
     &              FAC010 * ABSA(IND0+9,IG) +
     &              FAC110 * ABSA(IND0+10,IG))
     &              + SPECCOMB1 *
     &              (FAC001 * ABSA(IND1,IG) + 
     &              FAC101 * ABSA(IND1+1,IG) +
     &              FAC011 * ABSA(IND1+9,IG) +
     &              FAC111 * ABSA(IND1+10,IG)) 
C      &              + TAUSELF + TAUFOR
     &              + ADJCOLO2*ABSO2    

                taug(lay,ig) = taug(lay,ig)
               if( h2o_sb == 1 ) then 
                  taug(lay,ig) = taug(lay,ig) + TAUSELF
                endif
               if( h2o_for == 1 ) then 
                  taug(lay,ig) = taug(lay,ig) + TAUFOR
                endif        


               FRACS(LAY,IG) = FRACREFA(IG,JPL) + FPL *
     &              (FRACREFA(IG,JPL+1)-FRACREFA(IG,JPL))
 2020          CONTINUE
      ENDIF
 2500 CONTINUE

      RETURN
      END

C----------------------------------------------------------------------------

      SUBROUTINE TAUGB11
      use variables !NJE

      return
      end
 
C----------------------------------------------------------------------------

      SUBROUTINE TAUGB12
      use variables !NJE

C     BAND 12:  1800-2080 cm-1 (low - H2O,CO2)

      PARAMETER (mg=16, mxlay=203, NBANDS=16)
      PARAMETER (NTEMP=14,NETA=9,NREF=71)

C  Output

      COMMON /TAUGCOM/  TAUG(MXLAY,MG)
      COMMON /PLANKG/   FRACS(MXLAY,MG)                                       

C  Input

      COMMON /FEATURES/ NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)
      COMMON /PRECISE/  ONEMINUS
      COMMON /PROFILE/  NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
     &                  PZ(0:MXLAY),TZ(0:MXLAY)
      COMMON /PROFDATA/ laytrop,lay_switch_kabs,                                   
     &                  COLH2O(MXLAY),COLCO2(MXLAY),COLO3(MXLAY),  
     &                  COLN2O(MXLAY),COLCO(MXLAY),COLCH4(MXLAY),  
     &                  COLO2(MXLAY),COLBRD(MXLAY)
      COMMON /INTFAC/   FAC00(MXLAY),FAC01(MXLAY),                            
     &                  FAC10(MXLAY),FAC11(MXLAY)  
      COMMON /STD_REF/  PREF(NREF),PREFLOG(NREF),TREF(NREF),
     &                   CHI_STD(7,NREF)            
      COMMON /INTIND/   JP(MXLAY),JT(MXLAY),JT1(MXLAY)
      COMMON /REFRAT_ETA/ RAT_H2OCO2(MXLAY),RAT_H2OCO2_1(MXLAY),
     &                  RAT_H2OO3(MXLAY),RAT_H2OO3_1(MXLAY),
     &                  RAT_H2ON2O(MXLAY),RAT_H2ON2O_1(MXLAY),
     &                  RAT_H2OCH4(MXLAY),RAT_H2OCH4_1(MXLAY),
     &                  RAT_N2OCO2(MXLAY),RAT_N2OCO2_1(MXLAY),
     &                  RAT_O3CO2(MXLAY),RAT_O3CO2_1(MXLAY),
     &                  RAT_CO2O3(MXLAY),RAT_CO2O3_1(MXLAY),
     &                  RAT_O2_DRY(MXLAY)

      COMMON /SELF/     SELFFAC(MXLAY), SELFFRAC(MXLAY), INDSELF(MXLAY)
      COMMON /FOREIGN/  FORFAC(MXLAY), FORFRAC(MXLAY), INDFOR(MXLAY)
      COMMON /K12/      KA(NETA,NTEMP,NREF,MG),FORREF(4,MG),
     &                  SELFREF(10,MG)

      COMMON /CVRTAU/    HNAMTAU,HVRTAU

      CHARACTER*18       HNAMTAU,HVRTAU

      DIMENSION ABSA(8946,MG)
      DIMENSION FRACREFA(MG,9)

      REAL KA

C P = 174.1640 mbar, T= 215. K

      DATA (FRACREFA(IG, 1),IG=1,16) /
     &2.5190E-01,1.5140E-01,1.3760E-01,1.2021E-01,1.1057E-01,9.2104E-02,
     &5.7709E-02,4.0961E-02,2.8274E-02,3.2992E-03,2.4129E-03,1.4916E-03,
     &1.0859E-03,6.8566E-04,2.5569E-04,3.4972E-05/
      DATA (FRACREFA(IG, 2),IG=1,16) /
     &1.3468E-01,1.7766E-01,1.5817E-01,1.4588E-01,1.2452E-01,9.7226E-02,
     &6.4780E-02,5.0278E-02,3.5444E-02,3.7503E-03,3.0111E-03,2.1535E-03,
     &1.3752E-03,7.4498E-04,2.9461E-04,3.4972E-05/
      DATA (FRACREFA(IG, 3),IG=1,16) /
     &1.2155E-01,1.6994E-01,1.6330E-01,1.4806E-01,1.2945E-01,9.3241E-02,
     &6.9891E-02,5.4079E-02,3.7785E-02,4.0741E-03,3.0804E-03,2.4195E-03,
     &1.6697E-03,8.9880E-04,4.4761E-04,1.1626E-04/
      DATA (FRACREFA(IG, 4),IG=1,16) /
     &1.1715E-01,1.5390E-01,1.7087E-01,1.4914E-01,1.3273E-01,9.1537E-02,
     &7.3955E-02,5.7319E-02,3.9513E-02,4.2990E-03,3.3657E-03,2.4614E-03,
     &1.8346E-03,1.0056E-03,8.0215E-04,1.1626E-04/
      DATA (FRACREFA(IG, 5),IG=1,16) /
     &1.1481E-01,1.4646E-01,1.6841E-01,1.5422E-01,1.3045E-01,9.2200E-02,
     &7.7785E-02,5.9576E-02,4.1182E-02,4.4005E-03,3.5787E-03,2.6275E-03,
     &1.7782E-03,1.6034E-03,8.0218E-04,1.1626E-04/
      DATA (FRACREFA(IG, 6),IG=1,16) /
     &1.1369E-01,1.4009E-01,1.6723E-01,1.5441E-01,1.2533E-01,9.6626E-02,
     &8.1513E-02,6.2364E-02,4.2790E-02,4.4853E-03,3.7922E-03,2.8102E-03,
     &1.8585E-03,2.0998E-03,8.0215E-04,1.1626E-04/
      DATA (FRACREFA(IG, 7),IG=1,16) /
     &1.1214E-01,1.3189E-01,1.6339E-01,1.5531E-01,1.1866E-01,1.0563E-01,
     &8.5671E-02,6.5592E-02,4.4528E-02,4.7181E-03,3.8600E-03,2.9093E-03,
     &2.6891E-03,2.0998E-03,8.0213E-04,1.1626E-04/
      DATA (FRACREFA(IG, 8),IG=1,16) /
     &1.1015E-01,1.2704E-01,1.5158E-01,1.4482E-01,1.2773E-01,1.1273E-01,
     &9.0416E-02,6.9718E-02,4.6916E-02,4.9311E-03,3.8491E-03,4.0150E-03,
     &3.0925E-03,2.0999E-03,8.0207E-04,1.1626E-04/
      DATA (FRACREFA(IG, 9),IG=1,16) /
     &8.3975E-02,1.0375E-01,1.2161E-01,1.4135E-01,1.5020E-01,1.3197E-01,
     &1.0916E-01,8.0808E-02,5.5595E-02,6.1208E-03,5.1188E-03,4.2231E-03,
     &3.0937E-03,2.0998E-03,8.0227E-04,1.1626E-04/

      EQUIVALENCE (KA,ABSA)


C     Calculate reference ratio to be used in calculation of Planck
C     fraction in lower/upper atmosphere.

C     P =   174.164 mb 
      REFRAT_PLANCK_A = CHI_STD(1,22)/CHI_STD(2,22)

C     Compute the optical depth by interpolating in ln(pressure), 
C     temperature, and appropriate species.  Below LAYTROP, the water
C     vapor self-continuum adn foreign continuum is interpolated 
C     (in temperature) separately.  

      HVRTAU = '$Revision: 22602 $'

      DO 2500 LAY = 1, NLAYERS

         SPECCOMB = COLH2O(LAY) + RAT_H2OCO2(LAY)*COLCO2(LAY)
         SPECPARM = COLH2O(LAY)/SPECCOMB
         IF (SPECPARM .GE. ONEMINUS) SPECPARM = ONEMINUS
         SPECMULT = 8.*(SPECPARM)
         JS = 1 + INT(SPECMULT)
         FS = AMOD(SPECMULT,1.0)

         SPECCOMB1 = COLH2O(LAY) + RAT_H2OCO2_1(LAY)*COLCO2(LAY)
         SPECPARM1 = COLH2O(LAY)/SPECCOMB1
         IF (SPECPARM1 .GE. ONEMINUS) SPECPARM1 = ONEMINUS
         SPECMULT1 = 8.*(SPECPARM1)
         JS1 = 1 + INT(SPECMULT1)
         FS1 = AMOD(SPECMULT1,1.0)

         SPECCOMB_PLANCK = COLH2O(LAY)+REFRAT_PLANCK_A*COLCO2(LAY)
         SPECPARM_PLANCK = COLH2O(LAY)/SPECCOMB_PLANCK
         IF (SPECPARM_PLANCK .GE. ONEMINUS) SPECPARM_PLANCK=ONEMINUS
         SPECMULT_PLANCK = 8.*SPECPARM_PLANCK
         JPL= 1 + INT(SPECMULT_PLANCK)
         FPL = AMOD(SPECMULT_PLANCK,1.0)

         IND0 = ((JP(LAY)-1)*NTEMP+(JT(LAY)-1))*NSPA(12) + JS
         IND1 = (JP(LAY)*NTEMP+(JT1(LAY)-1))*NSPA(12) + JS1
         INDS = INDSELF(LAY)
         INDF = INDFOR(LAY)

         IF (SPECPARM .LT. 0.125 .AND. SPECPARM1 .LT. 0.125) THEN
            P = FS - 1
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC000 = FK0*FAC00(LAY)
            FAC100 = FK1*FAC00(LAY)
            FAC200 = FK2*FAC00(LAY)
            FAC010 = FK0*FAC10(LAY)
            FAC110 = FK1*FAC10(LAY)
            FAC210 = FK2*FAC10(LAY)

            P = FS1 - 1
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC001 = FK0*FAC01(LAY)
            FAC101 = FK1*FAC01(LAY)
            FAC201 = FK2*FAC01(LAY)
            FAC011 = FK0*FAC11(LAY)
            FAC111 = FK1*FAC11(LAY)
            FAC211 = FK2*FAC11(LAY)


            DO 2000 IG = 1, NG(12)
               TAUSELF = SELFFAC(LAY)* (SELFREF(INDS,IG) +
     &              SELFFRAC(LAY)  *
     &              (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
               TAUFOR = FORFAC(LAY) * (FORREF(INDF,IG) +
     &              FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &              FORREF(INDF,IG))) 
               TAUG(LAY,IG) = SPECCOMB *
     &              (FAC000 * ABSA(IND0,IG) +
     &              FAC100 * ABSA(IND0+1,IG) +
     &              FAC200 * ABSA(IND0+2,IG) +
     &              FAC010 * ABSA(IND0+9,IG) +
     &              FAC110 * ABSA(IND0+10,IG) +
     &              FAC210 * ABSA(IND0+11,IG))
     &              + SPECCOMB1 *
     &              (FAC001 * ABSA(IND1,IG) + 
     &              FAC101 * ABSA(IND1+1,IG) +
     &              FAC201 * ABSA(IND1+2,IG) +
     &              FAC011 * ABSA(IND1+9,IG) +
     &              FAC111 * ABSA(IND1+10,IG) +
     &              FAC211 * ABSA(IND1+11,IG)) 
C      &              + TAUSELF + TAUFOR

            taug(lay,ig) = taug(lay,ig)
               if( h2o_sb == 1 ) then 
                  taug(lay,ig) = taug(lay,ig) + TAUSELF
                endif
               if( h2o_for == 1 ) then 
                  taug(lay,ig) = taug(lay,ig) + TAUFOR
                endif

       
               FRACS(LAY,IG) = FRACREFA(IG,JPL) + FPL *
     &              (FRACREFA(IG,JPL+1)-FRACREFA(IG,JPL))
 2000          CONTINUE
      ELSE IF (SPECPARM .GT. 0.875 .AND. SPECPARM1 .GT. 0.875) THEN
            P = -FS 
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC000 = FK0*FAC00(LAY)
            FAC100 = FK1*FAC00(LAY)
            FAC200 = FK2*FAC00(LAY)
            FAC010 = FK0*FAC10(LAY)
            FAC110 = FK1*FAC10(LAY)
            FAC210 = FK2*FAC10(LAY)

            P = -FS1 
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC001 = FK0*FAC01(LAY)
            FAC101 = FK1*FAC01(LAY)
            FAC201 = FK2*FAC01(LAY)
            FAC011 = FK0*FAC11(LAY)
            FAC111 = FK1*FAC11(LAY)
            FAC211 = FK2*FAC11(LAY)
            DO 2010 IG = 1, NG(12)
               TAUSELF = SELFFAC(LAY)* (SELFREF(INDS,IG) +
     &              SELFFRAC(LAY)  *
     &           (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
               TAUFOR =  FORFAC(LAY) * (FORREF(INDF,IG) +
     &           FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &           FORREF(INDF,IG))) 
               TAUG(LAY,IG) = SPECCOMB * 
     &              (FAC200 * ABSA(IND0-1,IG) +
     &              FAC100 * ABSA(IND0,IG) +
     &              FAC000 * ABSA(IND0+1,IG) +
     &              FAC210 * ABSA(IND0+8,IG) +
     &              FAC110 * ABSA(IND0+9,IG) +
     &              FAC010 * ABSA(IND0+10,IG))
     &              + SPECCOMB1 *
     &              (FAC201 * ABSA(IND1-1,IG) +
     &              FAC101 * ABSA(IND1,IG) +
     &              FAC001 * ABSA(IND1+1,IG) +
     &              FAC211 * ABSA(IND1+8,IG) +
     &              FAC111 * ABSA(IND1+9,IG) +
     &              FAC011 * ABSA(IND1+10,IG))
C      &              + TAUSELF + TAUFOR

            taug(lay,ig) = taug(lay,ig)
               if( h2o_sb == 1 ) then 
                  taug(lay,ig) = taug(lay,ig) + TAUSELF
                endif
               if( h2o_for == 1 ) then 
                  taug(lay,ig) = taug(lay,ig) + TAUFOR
                endif

               FRACS(LAY,IG) = FRACREFA(IG,JPL) + FPL *
     &              (FRACREFA(IG,JPL+1)-FRACREFA(IG,JPL))
 2010           CONTINUE
       ELSE
            FAC000 = (1. - FS) * FAC00(LAY)
            FAC010 = (1. - FS) * FAC10(LAY)
            FAC100 = FS * FAC00(LAY)
            FAC110 = FS * FAC10(LAY)

            FAC001 = (1. - FS1) * FAC01(LAY)
            FAC011 = (1. - FS1) * FAC11(LAY)
            FAC101 = FS1 * FAC01(LAY)
            FAC111 = FS1 * FAC11(LAY)

            DO 2020 IG = 1, NG(12)
               TAUSELF = SELFFAC(LAY)* (SELFREF(INDS,IG) +
     &              SELFFRAC(LAY)  *
     &              (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
               TAUFOR =  FORFAC(LAY) * (FORREF(INDF,IG) +
     &              FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &              FORREF(INDF,IG))) 
               TAUG(LAY,IG) = SPECCOMB * 
     &              (FAC000 * ABSA(IND0,IG) +
     &              FAC100 * ABSA(IND0+1,IG) +
     &              FAC010 * ABSA(IND0+9,IG) +
     &              FAC110 * ABSA(IND0+10,IG))
     &              + SPECCOMB1 *
     &              (FAC001 * ABSA(IND1,IG) + 
     &              FAC101 * ABSA(IND1+1,IG) +
     &              FAC011 * ABSA(IND1+9,IG) +
     &              FAC111 * ABSA(IND1+10,IG)) 
C      &              + TAUSELF + TAUFOR
  
              taug(lay,ig) = taug(lay,ig)
               if( h2o_sb == 1 ) then 
                  taug(lay,ig) = taug(lay,ig) + TAUSELF
                endif
               if( h2o_for == 1 ) then 
                  taug(lay,ig) = taug(lay,ig) + TAUFOR
                endif

               FRACS(LAY,IG) = FRACREFA(IG,JPL) + FPL *
     &              (FRACREFA(IG,JPL+1)-FRACREFA(IG,JPL))
 2020          CONTINUE
      ENDIF
 2500 CONTINUE


      RETURN
      END

C----------------------------------------------------------------------------

      SUBROUTINE TAUGB13
      use variables !NJE

C     BAND 13:  2080-2250 cm-1 (low key - H2O,CO2; low minor - N2O )

      PARAMETER (mg=16, mxlay=203, MXMOL=38,NBANDS=16)
      PARAMETER (NTEMP=14,NETA=9,NREF=71)

C  Output

      COMMON /TAUGCOM/  TAUG(MXLAY,MG)
      COMMON /PLANKG/   FRACS(MXLAY,MG)                                       

C  Input

      COMMON /FEATURES/ NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)
      COMMON /PRECISE/  ONEMINUS
      COMMON /PROFILE/  NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
     &                  PZ(0:MXLAY),TZ(0:MXLAY)
      COMMON /PROFDATA/ laytrop,lay_switch_kabs,                                   
     &                  COLH2O(MXLAY),COLCO2(MXLAY),COLO3(MXLAY),  
     &                  COLN2O(MXLAY),COLCO(MXLAY),COLCH4(MXLAY),  
     &                  COLO2(MXLAY),COLBRD(MXLAY)
      COMMON /SPECIES/  COLDRY(MXLAY),WKL(MXMOL,MXLAY),WBROAD(MXLAY),
     &                  COLMOL(MXLAY),NMOL
      COMMON /STD_REF/  PREF(NREF),PREFLOG(NREF),TREF(NREF),
     &                  CHI_STD(7,NREF)
      COMMON /INTFAC/   FAC00(MXLAY),FAC01(MXLAY),                            
     &                  FAC10(MXLAY),FAC11(MXLAY)                             
      COMMON /INTIND/   JP(MXLAY),JT(MXLAY),JT1(MXLAY)
      COMMON /REFRAT_ETA/ RAT_H2OCO2(MXLAY),RAT_H2OCO2_1(MXLAY),
     &                  RAT_H2OO3(MXLAY),RAT_H2OO3_1(MXLAY),
     &                  RAT_H2ON2O(MXLAY),RAT_H2ON2O_1(MXLAY),
     &                  RAT_H2OCH4(MXLAY),RAT_H2OCH4_1(MXLAY),
     &                  RAT_N2OCO2(MXLAY),RAT_N2OCO2_1(MXLAY),
     &                  RAT_O3CO2(MXLAY),RAT_O3CO2_1(MXLAY),
     &                  RAT_CO2O3(MXLAY),RAT_CO2O3_1(MXLAY),
     &                  RAT_O2_DRY(MXLAY)

      COMMON /SELF/     SELFFAC(MXLAY), SELFFRAC(MXLAY), INDSELF(MXLAY)
      COMMON /FOREIGN/  FORFAC(MXLAY), FORFRAC(MXLAY), INDFOR(MXLAY)
      COMMON /MINOR/    MINORFRAC(MXLAY), INDMINOR(MXLAY), 
     &                  SCALEMINOR(MXLAY),SCALEMINORN2(MXLAY)
      COMMON /K13/      KA(NETA,NTEMP,NREF,MG),FORREF(4,MG),
     &                  SELFREF(10,MG),  KA_MN2O(9,19,MG)

      COMMON /CVRTAU/    HNAMTAU,HVRTAU

      CHARACTER*18       HNAMTAU,HVRTAU

      DIMENSION ABSA(8946,MG)
      DIMENSION FRACREFA(MG,9)
      REAL KA,KA_MN2O,MINORFRAC

      REAL N2OM1, N2OM2

C Planck fraction mapping level : P=473.4280 mb, T = 260. K      

      DATA (FRACREFA(IG, 1),IG=1,16) /
     &1.8325E-01,1.6959E-01,1.5171E-01,1.3094E-01,1.1022E-01,9.1281E-02,
     &7.1660E-02,5.0014E-02,3.0249E-02,3.2768E-03,2.7053E-03,2.1258E-03,
     &1.5590E-03,9.8653E-04,3.6841E-04,5.1629E-05/
      DATA (FRACREFA(IG, 2),IG=1,16) /
     &1.8031E-01,1.6787E-01,1.4814E-01,1.3078E-01,1.1123E-01,9.3148E-02,
     &7.4331E-02,5.1934E-02,3.1032E-02,3.3342E-03,2.7974E-03,2.1258E-03,
     &1.5590E-03,9.8654E-04,3.6841E-04,5.1629E-05/
      DATA (FRACREFA(IG, 3),IG=1,16) /
     &1.7878E-01,1.6382E-01,1.4833E-01,1.3122E-01,1.1235E-01,9.4057E-02,
     &7.5174E-02,5.3081E-02,3.1740E-02,3.4109E-03,2.8156E-03,2.2073E-03,
     &1.6178E-03,9.8653E-04,3.6842E-04,5.1629E-05/
      DATA (FRACREFA(IG, 4),IG=1,16) /
     &1.7744E-01,1.6120E-01,1.4810E-01,1.2952E-01,1.1474E-01,9.5644E-02,
     &7.4855E-02,5.3864E-02,3.2905E-02,3.4827E-03,2.8927E-03,2.2478E-03,
     &1.6796E-03,1.0104E-03,3.6842E-04,5.1629E-05/
      DATA (FRACREFA(IG, 5),IG=1,16) /
     &1.7467E-01,1.5867E-01,1.4785E-01,1.3074E-01,1.1520E-01,9.6766E-02,
     &7.5963E-02,5.3864E-02,3.4074E-02,3.6818E-03,2.9942E-03,2.3194E-03,
     &1.6965E-03,1.1012E-03,3.6842E-04,5.1629E-05/
      DATA (FRACREFA(IG, 6),IG=1,16) /
     &1.7127E-01,1.5801E-01,1.4380E-01,1.3249E-01,1.1597E-01,9.9468E-02,
     &7.6588E-02,5.4240E-02,3.5512E-02,3.7076E-03,3.1748E-03,2.4112E-03,
     &1.7642E-03,1.1124E-03,4.4147E-04,5.1629E-05/
      DATA (FRACREFA(IG, 7),IG=1,16) /
     &1.6421E-01,1.5660E-01,1.4256E-01,1.3432E-01,1.1797E-01,1.0109E-01,
     &7.8029E-02,5.4346E-02,3.7304E-02,4.1167E-03,3.2043E-03,2.6446E-03,
     &1.8864E-03,1.1487E-03,4.6845E-04,1.0013E-04/
      DATA (FRACREFA(IG, 8),IG=1,16) /
     &1.5203E-01,1.4845E-01,1.4663E-01,1.3737E-01,1.2053E-01,1.0365E-01,
     &8.0781E-02,5.5946E-02,3.9497E-02,4.4959E-03,3.7080E-03,2.7979E-03,
     &2.1134E-03,1.2953E-03,6.1320E-04,1.0013E-04/
      DATA (FRACREFA(IG, 9),IG=1,16) /
     &1.8676E-01,1.6692E-01,1.4609E-01,1.2753E-01,1.0662E-01,8.4446E-02,
     &6.4705E-02,5.2567E-02,4.5691E-02,5.3163E-03,4.5644E-03,3.6256E-03,
     &2.6427E-03,1.7549E-03,6.7040E-04,1.0013E-04/

C Minor gas mapping levels :
C     LOWER - N2O, P = 473.428, T = 260.0 K

      EQUIVALENCE (KA,ABSA)


C     Calculate reference ratio to be used in calculation of Planck
C     fraction in lower/upper atmosphere.

C     P = 473.420 mb (Level 17)
      REFRAT_PLANCK_A = CHI_STD(1,17)/CHI_STD(2,17)

C     P = 473.428  (Level 17)
      REFRAT_M_A = CHI_STD(1,17)/CHI_STD(2,17)


C     Compute the optical depth by interpolating in ln(pressure), 
C     temperature, and appropriate species.  Below LAYTROP, the water
C     vapor self-continuum and foreign continuum is interpolated 
C     (in temperature) separately.  

      HVRTAU = '$Revision: 22602 $'

      DO 2500 LAY = 1, NLAYERS

         SPECCOMB = COLH2O(LAY) + RAT_H2OCO2(LAY)*COLCO2(LAY)
         SPECPARM = COLH2O(LAY)/SPECCOMB
         IF (SPECPARM .GE. ONEMINUS) SPECPARM = ONEMINUS
         SPECMULT = 8.*(SPECPARM)
         JS = 1 + INT(SPECMULT)
         FS = AMOD(SPECMULT,1.0)


         SPECCOMB1 = COLH2O(LAY) + RAT_H2OCO2_1(LAY)*COLCO2(LAY)
         SPECPARM1 = COLH2O(LAY)/SPECCOMB1
         IF (SPECPARM1 .GE. ONEMINUS) SPECPARM1 = ONEMINUS
         SPECMULT1 = 8.*(SPECPARM1)
         JS1 = 1 + INT(SPECMULT1)
         FS1 = AMOD(SPECMULT1,1.0)

         SPECCOMB_MN2O = COLH2O(LAY) + REFRAT_M_A*COLCO2(LAY)
         SPECPARM_MN2O = COLH2O(LAY)/SPECCOMB_MN2O
         IF (SPECPARM_MN2O .GE. ONEMINUS) SPECPARM_MN2O = ONEMINUS
         SPECMULT_MN2O = 8.*SPECPARM_MN2O
         JMN2O = 1 + INT(SPECMULT_MN2O)
         FMN2O = AMOD(SPECMULT_MN2O,1.0)


         SPECCOMB_PLANCK = COLH2O(LAY)+REFRAT_PLANCK_A*COLCO2(LAY)
         SPECPARM_PLANCK = COLH2O(LAY)/SPECCOMB_PLANCK
         IF (SPECPARM_PLANCK .GE. ONEMINUS) SPECPARM_PLANCK=ONEMINUS
         SPECMULT_PLANCK = 8.*SPECPARM_PLANCK
         JPL= 1 + INT(SPECMULT_PLANCK)
         FPL = AMOD(SPECMULT_PLANCK,1.0)

         IND0 = ((JP(LAY)-1)*NTEMP+(JT(LAY)-1))*NSPA(13) + JS
         IND1 = (JP(LAY)*NTEMP+(JT1(LAY)-1))*NSPA(13) + JS1
         INDS = INDSELF(LAY)
         INDF = INDFOR(LAY)
         INDM = INDMINOR(LAY)

         IF (SPECPARM .LT. 0.125 .AND. SPECPARM1 .LT. 0.125) THEN
            P = FS - 1
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC000 = FK0*FAC00(LAY)
            FAC100 = FK1*FAC00(LAY)
            FAC200 = FK2*FAC00(LAY)
            FAC010 = FK0*FAC10(LAY)
            FAC110 = FK1*FAC10(LAY)
            FAC210 = FK2*FAC10(LAY)

            P = FS1 - 1
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC001 = FK0*FAC01(LAY)
            FAC101 = FK1*FAC01(LAY)
            FAC201 = FK2*FAC01(LAY)
            FAC011 = FK0*FAC11(LAY)
            FAC111 = FK1*FAC11(LAY)
            FAC211 = FK2*FAC11(LAY)

            DO 2000 IG = 1, NG(13)
               TAUSELF = SELFFAC(LAY)* (SELFREF(INDS,IG) +
     &              SELFFRAC(LAY)  *
     &              (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
               TAUFOR = FORFAC(LAY) * (FORREF(INDF,IG) +
     &              FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &              FORREF(INDF,IG))) 
               N2OM1 = KA_MN2O(JMN2O,INDM,IG) + FMN2O*
     &              (KA_MN2O(JMN2O+1,INDM,IG)-
     &              KA_MN2O(JMN2O,INDM,IG))
               N2OM2 = KA_MN2O(JMN2O,INDM+1,IG) + FMN2O*
     &              (KA_MN2O(JMN2O+1,INDM+1,IG)-
     &              KA_MN2O(JMN2O,INDM+1,IG))
               ABSN2O = N2OM1 + 
     &              MINORFRAC(LAY) * (N2OM2 - N2OM1)

               TAUG(LAY,IG) = SPECCOMB *
     &              (FAC000 * ABSA(IND0,IG) +
     &              FAC100 * ABSA(IND0+1,IG) +
     &              FAC200 * ABSA(IND0+2,IG) +
     &              FAC010 * ABSA(IND0+9,IG) +
     &              FAC110 * ABSA(IND0+10,IG) +
     &              FAC210 * ABSA(IND0+11,IG))+
     &              SPECCOMB1 *
     &              (FAC001 * ABSA(IND1,IG) + 
     &              FAC101 * ABSA(IND1+1,IG) +
     &              FAC201 * ABSA(IND1+2,IG) +
     &              FAC011 * ABSA(IND1+9,IG) +
     &              FAC111 * ABSA(IND1+10,IG) +
     &              FAC211 * ABSA(IND1+11,IG)) 
C      &              + TAUSELF + TAUFOR
     &              + COLN2O(LAY)*ABSN2O


                taug(lay,ig) = taug(lay,ig)
               if( h2o_sb == 1 ) then 
                  taug(lay,ig) = taug(lay,ig) + TAUSELF
                endif
               if( h2o_for == 1 ) then 
                  taug(lay,ig) = taug(lay,ig) + TAUFOR
                endif

               FRACS(LAY,IG) = FRACREFA(IG,JPL) + FPL *
     &              (FRACREFA(IG,JPL+1)-FRACREFA(IG,JPL))
 2000          CONTINUE
      ELSE IF (SPECPARM .GT. 0.875 .AND. SPECPARM1 .GT. 0.875) THEN
            P = -FS 
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC000 = FK0*FAC00(LAY)
            FAC100 = FK1*FAC00(LAY)
            FAC200 = FK2*FAC00(LAY)
            FAC010 = FK0*FAC10(LAY)
            FAC110 = FK1*FAC10(LAY)
            FAC210 = FK2*FAC10(LAY)

            P = -FS1 
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC001 = FK0*FAC01(LAY)
            FAC101 = FK1*FAC01(LAY)
            FAC201 = FK2*FAC01(LAY)
            FAC011 = FK0*FAC11(LAY)
            FAC111 = FK1*FAC11(LAY)
            FAC211 = FK2*FAC11(LAY)

            DO 2010 IG = 1, NG(13)
               TAUSELF = SELFFAC(LAY)* (SELFREF(INDS,IG) +
     &              SELFFRAC(LAY)  *
     &              (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
               TAUFOR = FORFAC(LAY) * (FORREF(INDF,IG) +
     &              FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &              FORREF(INDF,IG))) 
               N2OM1 = KA_MN2O(JMN2O,INDM,IG) + FMN2O*
     &              (KA_MN2O(JMN2O+1,INDM,IG)-
     &              KA_MN2O(JMN2O,INDM,IG))
               N2OM2 = KA_MN2O(JMN2O,INDM+1,IG) + FMN2O*
     &              (KA_MN2O(JMN2O+1,INDM+1,IG)-
     &              KA_MN2O(JMN2O,INDM+1,IG))
               ABSN2O = N2OM1 + 
     &              MINORFRAC(LAY) * (N2OM2 - N2OM1)

               TAUG(LAY,IG) = SPECCOMB * 
     &              (FAC200 * ABSA(IND0-1,IG) +
     &              FAC100 * ABSA(IND0,IG) +
     &              FAC000 * ABSA(IND0+1,IG) +
     &              FAC210 * ABSA(IND0+8,IG) +
     &              FAC110 * ABSA(IND0+9,IG) +
     &              FAC010 * ABSA(IND0+10,IG))
     &              + SPECCOMB1 *
     &              (FAC201 * ABSA(IND1-1,IG) +
     &              FAC101 * ABSA(IND1,IG) +
     &              FAC001 * ABSA(IND1+1,IG) +
     &              FAC211 * ABSA(IND1+8,IG) +
     &              FAC111 * ABSA(IND1+9,IG) +
     &              FAC011 * ABSA(IND1+10,IG))
C      &              + TAUSELF + TAUFOR
     &              + COLN2O(LAY)*ABSN2O

              taug(lay,ig) = taug(lay,ig)
               if( h2o_sb == 1 ) then 
                  taug(lay,ig) = taug(lay,ig) + TAUSELF
                endif
               if( h2o_for == 1 ) then 
                  taug(lay,ig) = taug(lay,ig) + TAUFOR
                endif

                FRACS(LAY,IG) = FRACREFA(IG,JPL) + FPL *
     &               (FRACREFA(IG,JPL+1)-FRACREFA(IG,JPL))
 2010           CONTINUE
       ELSE
            FAC000 = (1. - FS) * FAC00(LAY)
            FAC010 = (1. - FS) * FAC10(LAY)
            FAC100 = FS * FAC00(LAY)
            FAC110 = FS * FAC10(LAY)

            FAC001 = (1. - FS1) * FAC01(LAY)
            FAC011 = (1. - FS1) * FAC11(LAY)
            FAC101 = FS1 * FAC01(LAY)
            FAC111 = FS1 * FAC11(LAY)

            DO 2020 IG = 1, NG(13)
               TAUSELF = SELFFAC(LAY)* (SELFREF(INDS,IG) +
     &              SELFFRAC(LAY)  *
     &              (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
               TAUFOR = FORFAC(LAY) * (FORREF(INDF,IG) +
     &              FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &              FORREF(INDF,IG))) 
               N2OM1 = KA_MN2O(JMN2O,INDM,IG) + FMN2O*
     &              (KA_MN2O(JMN2O+1,INDM,IG)-
     &              KA_MN2O(JMN2O,INDM,IG))
               N2OM2 = KA_MN2O(JMN2O,INDM+1,IG) + FMN2O*
     &              (KA_MN2O(JMN2O+1,INDM+1,IG)-
     &              KA_MN2O(JMN2O,INDM+1,IG))
               ABSN2O = N2OM1 + 
     &              MINORFRAC(LAY) * (N2OM2 - N2OM1)

               TAUG(LAY,IG) = SPECCOMB * 
     &              (FAC000 * ABSA(IND0,IG) +
     &              FAC100 * ABSA(IND0+1,IG) +
     &              FAC010 * ABSA(IND0+9,IG) +
     &              FAC110 * ABSA(IND0+10,IG))
     &              + SPECCOMB1 *
     &              (FAC001 * ABSA(IND1,IG) + 
     &              FAC101 * ABSA(IND1+1,IG) +
     &              FAC011 * ABSA(IND1+9,IG) +
     &              FAC111 * ABSA(IND1+10,IG)) 
C      &              + TAUSELF + TAUFOR
     &              + COLN2O(LAY)*ABSN2O


                taug(lay,ig) = taug(lay,ig)
               if( h2o_sb == 1 ) then 
                  taug(lay,ig) = taug(lay,ig) + TAUSELF
                endif
               if( h2o_for == 1 ) then 
                  taug(lay,ig) = taug(lay,ig) + TAUFOR
                endif

               FRACS(LAY,IG) = FRACREFA(IG,JPL) + FPL *
     &              (FRACREFA(IG,JPL+1)-FRACREFA(IG,JPL))

 2020          CONTINUE
      ENDIF
 2500 CONTINUE

      RETURN
      END

*******************************************************************************

      SUBROUTINE TAUGB14
      use variables !NJE

C     BAND 14:  2250-2700 cm-1 (low - CO2 - minor N2)

      PARAMETER (mg=16, mxlay=203, NBANDS=16)
      PARAMETER (NTEMP=14,NETA=9,NREF=71)

C  Output

      COMMON /TAUGCOM/  TAUG(MXLAY,MG)
      COMMON /PLANKG/   FRACS(MXLAY,MG)                                       

C  Input

      COMMON /FEATURES/ NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)
      COMMON /PROFILE/  NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
     &                  PZ(0:MXLAY),TZ(0:MXLAY)
      COMMON /PROFDATA/ laytrop,lay_switch_kabs,                                   
     &                  COLH2O(MXLAY),COLCO2(MXLAY),COLO3(MXLAY),  
     &                  COLN2O(MXLAY),COLCO(MXLAY),COLCH4(MXLAY),  
     &                  COLO2(MXLAY),COLBRD(MXLAY)
      COMMON /INTFAC/   FAC00(MXLAY),FAC01(MXLAY),                            
     &                  FAC10(MXLAY),FAC11(MXLAY)                             
      COMMON /INTIND/   JP(MXLAY),JT(MXLAY),JT1(MXLAY)
      COMMON /STD_REF/  PREF(nref),PREFLOG(nref),TREF(nref),
     &                  CHI_std(7,nref)
      COMMON /REFRAT_ETA/ RAT_H2OCO2(MXLAY),RAT_H2OCO2_1(MXLAY),
     &                  RAT_H2OO3(MXLAY),RAT_H2OO3_1(MXLAY),
     &                  RAT_H2ON2O(MXLAY),RAT_H2ON2O_1(MXLAY),
     &                  RAT_H2OCH4(MXLAY),RAT_H2OCH4_1(MXLAY),
     &                  RAT_N2OCO2(MXLAY),RAT_N2OCO2_1(MXLAY),
     &                  RAT_O3CO2(MXLAY),RAT_O3CO2_1(MXLAY),
     &                  RAT_CO2O3(MXLAY),RAT_CO2O3_1(MXLAY),
     &                  RAT_O2_DRY(MXLAY)

      COMMON /SELF/     SELFFAC(MXLAY), SELFFRAC(MXLAY), INDSELF(MXLAY)
      COMMON /FOREIGN/  FORFAC(MXLAY), FORFRAC(MXLAY), INDFOR(MXLAY)
      COMMON /MINOR/    MINORFRAC(MXLAY), INDMINOR(MXLAY),
     &                  SCALEMINOR(MXLAY),SCALEMINORN2(MXLAY)
      COMMON /K14/      KA(NTEMP,NREF,MG), FORREF(4,MG), SELFREF(10,MG),
     &                  KA_MN2(19,MG),KA_MN2_O2(19,MG)

      COMMON /CVRTAU/    HNAMTAU,HVRTAU

      CHARACTER*18       HNAMTAU,HVRTAU

      DIMENSION ABSA(994,MG)
      DIMENSION FRACREFA(MG)

C Planck fraction mapping level : P = 142.5940 mb, T = 215 K
       DATA FRACREFA /
     &8.1830E-02,8.0214E-02,4.5215E-02,1.5383E-01,2.4637E-01,1.6529E-01,
     &7.5689E-02,8.3903E-02,4.9078E-02,5.4375E-03,4.4760E-03,3.6267E-03,
     &2.6343E-03,1.6825E-03,6.3361E-04,8.9164E-05/



      Equivalence (KA,ABSA),(KB,ABSB)
      REAL KA, KA_MN2, KA_MN2_O2, MINORFRAC

C     Compute the optical depth by interpolating in ln(pressure) and 
C     temperature.  Below LAYTROP, the water vapor self-continuum 
C     and foreign continuum is interpolated (in temperature) separately.  

      HVRTAU = '$Revision: 22602 $'

      DO 2500 LAY = 1, NLAYERS
         IND0 = ((JP(LAY)-1)*NTEMP+(JT(LAY)-1))*NSPA(14) + 1
         IND1 = (JP(LAY)*NTEMP+(JT1(LAY)-1))*NSPA(14) + 1
         INDS = INDSELF(LAY)
         INDF = INDFOR(LAY)
         INDM = INDMINOR(LAY)

         DO 2000 IG = 1, NG(14)
            TAUSELF = SELFFAC(LAY) * (SELFREF(INDS,IG) + 
     &           SELFFRAC(LAY) *
     &           (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
            TAUFOR =  FORFAC(LAY) * (FORREF(INDF,IG) +
     &           FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &           FORREF(INDF,IG))) 

               ABSN2 = KA_MN2(INDM,IG) + MINORFRAC(LAY) *
     &                KA_MN2(INDM+1,IG)             

               ABSN2_O2 = KA_MN2_O2(INDM,IG) + MINORFRAC(LAY) * 
     &                    KA_MN2_O2(INDM+1,IG)

               ratio = (chi_std(7,jp(lay))-rat_o2_dry(lay))/
     &                 chi_std(7,jp(lay))
               tauminor = ratio*absn2+(1.0-ratio)*absn2_o2
               tauminor = tauminor*scaleminor(lay)*colbrd(lay)

            TAUG(LAY,IG) = COLCO2(LAY) *
     &          (FAC00(LAY) * ABSA(IND0,IG) +
     &           FAC10(LAY) * ABSA(IND0+1,IG) +
     &           FAC01(LAY) * ABSA(IND1,IG) + 
     &           FAC11(LAY) * ABSA(IND1+1,IG))
C      &           + TAUSELF + TAUFOR
     &           + tauminor

            taug(lay,ig) = taug(lay,ig)
               if( h2o_sb == 1 ) then 
                  taug(lay,ig) = taug(lay,ig) + TAUSELF
                endif
               if( h2o_for == 1 ) then 
                  taug(lay,ig) = taug(lay,ig) + TAUFOR
                endif

            FRACS(LAY,IG) = FRACREFA(IG)
 2000    CONTINUE
 2500 CONTINUE

      RETURN
      END

C----------------------------------------------------------------------------


C----------------------------------------------------------------------------


      SUBROUTINE TAUGB15
      use variables !NJE

C     BAND 15:  2700-3250 cm-1 (low key- H2O,CO2 - minor CH4)
c     (used to be band 16, more or less)

      PARAMETER (mg=16, mxlay=203, NBANDS=16)
      PARAMETER (NTEMP=14,NETA=9,NREF=71)

C  Output

      COMMON /TAUGCOM/  TAUG(MXLAY,MG)
      COMMON /PLANKG/   FRACS(MXLAY,MG)                                       

C  Input

      COMMON /FEATURES/ NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)
      COMMON /PRECISE/  ONEMINUS
      COMMON /PROFILE/  NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
     &                  PZ(0:MXLAY),TZ(0:MXLAY)
      COMMON /PROFDATA/ laytrop,lay_switch_kabs,                                   
     &                  COLH2O(MXLAY),COLCO2(MXLAY),COLO3(MXLAY),  
     &                  COLN2O(MXLAY),COLCO(MXLAY),COLCH4(MXLAY),  
     &                  COLO2(MXLAY),COLBRD(MXLAY)
      COMMON /STD_REF/  PREF(NREF),PREFLOG(NREF),TREF(NREF),
     &                  CHI_STD(7,NREF)
      COMMON /INTFAC/   FAC00(MXLAY),FAC01(MXLAY),                            
     &                  FAC10(MXLAY),FAC11(MXLAY)                             
      COMMON /INTIND/   JP(MXLAY),JT(MXLAY),JT1(MXLAY)
      COMMON /MINOR/    MINORFRAC(MXLAY), INDMINOR(MXLAY), 
     &                  SCALEMINOR(MXLAY),SCALEMINORN2(MXLAY)
      COMMON /REFRAT_ETA/ RAT_H2OCO2(MXLAY),RAT_H2OCO2_1(MXLAY),
     &                  RAT_H2OO3(MXLAY),RAT_H2OO3_1(MXLAY),
     &                  RAT_H2ON2O(MXLAY),RAT_H2ON2O_1(MXLAY),
     &                  RAT_H2OCH4(MXLAY),RAT_H2OCH4_1(MXLAY),
     &                  RAT_N2OCO2(MXLAY),RAT_N2OCO2_1(MXLAY),
     &                  RAT_O3CO2(MXLAY),RAT_O3CO2_1(MXLAY),
     &                  RAT_CO2O3(MXLAY),RAT_CO2O3_1(MXLAY),
     &                  RAT_O2_DRY(MXLAY)

      COMMON /FOREIGN/  FORFAC(MXLAY), FORFRAC(MXLAY), INDFOR(MXLAY)
      COMMON /SELF/     SELFFAC(MXLAY), SELFFRAC(MXLAY), INDSELF(MXLAY)
      COMMON /K15/      KA(NETA,NTEMP,NREF,MG), FORREF(4,MG),
     &                  SELFREF(10,MG), KA_MCH4(9,19,MG)


      COMMON /CVRTAU/    HNAMTAU,HVRTAU

      CHARACTER*18       HNAMTAU,HVRTAU

      DIMENSION ABSA(8946,MG)
      DIMENSION FRACREFA(MG,9), FRACREFB(MG)

C Planck fraction mapping level: P = 387.6100 mbar, T = 245.0 K


      EQUIVALENCE (KA,ABSA)
      REAL KA,KA_MCH4, MINORFRAC

      DATA (FRACREFA(IG, 1),IG=1,16) /
     &1.1050E-01,1.0526E-01,1.1114E-01,1.4218E-01,8.0046E-02,2.0454E-01,
     &1.0385E-01,7.6428E-02,5.2340E-02,6.0266E-03,5.0414E-03,1.7703E-03,
     &4.6542E-04,2.9102E-04,1.0077E-04,1.4191E-05/
      DATA (FRACREFA(IG, 2),IG=1,16) /
     &2.2118E-01,2.5703E-01,1.9177E-01,1.1935E-01,7.3870E-02,5.4036E-02,
     &3.9797E-02,2.4102E-02,1.4024E-02,1.4854E-03,1.1901E-03,9.7400E-04,
     &6.6151E-04,3.9649E-04,1.2659E-04,1.4837E-05/
      DATA (FRACREFA(IG, 3),IG=1,16) /
     &2.2766E-01,2.5458E-01,1.9080E-01,1.1564E-01,7.4330E-02,5.4154E-02,
     &3.9887E-02,2.4086E-02,1.4026E-02,1.4856E-03,1.1902E-03,9.7388E-04,
     &6.6166E-04,3.9632E-04,1.2657E-04,1.4837E-05/
      DATA (FRACREFA(IG, 4),IG=1,16) /
     &2.3064E-01,2.5372E-01,1.8889E-01,1.1513E-01,7.4592E-02,5.4215E-02,
     &3.9880E-02,2.4072E-02,1.4027E-02,1.4855E-03,1.1904E-03,9.7366E-04,
     &6.6181E-04,3.9619E-04,1.2653E-04,1.4837E-05/
      DATA (FRACREFA(IG, 5),IG=1,16) /
     &2.3265E-01,2.5254E-01,1.8827E-01,1.1482E-01,7.4677E-02,5.4242E-02,
     &3.9870E-02,2.4059E-02,1.4026E-02,1.4858E-03,1.1903E-03,9.7362E-04,
     &6.6184E-04,3.9615E-04,1.2645E-04,1.4837E-05/
      DATA (FRACREFA(IG, 6),IG=1,16) /
     &2.3408E-01,2.5163E-01,1.8811E-01,1.1448E-01,7.4732E-02,5.4205E-02,
     &3.9852E-02,2.4035E-02,1.4024E-02,1.4860E-03,1.1904E-03,9.7330E-04,
     &6.6228E-04,3.9572E-04,1.2635E-04,1.4840E-05/
      DATA (FRACREFA(IG, 7),IG=1,16) /
     &2.3560E-01,2.5047E-01,1.8910E-01,1.1352E-01,7.4513E-02,5.4127E-02,
     &3.9818E-02,2.3985E-02,1.4024E-02,1.4869E-03,1.1905E-03,9.7264E-04,
     &6.6302E-04,3.9483E-04,1.2625E-04,1.4840E-05/
      DATA (FRACREFA(IG, 8),IG=1,16) /
     &2.3684E-01,2.4957E-01,1.9233E-01,1.1127E-01,7.3712E-02,5.3890E-02,
     &3.9713E-02,2.3815E-02,1.4022E-02,1.4886E-03,1.1902E-03,9.7122E-04,
     &6.6251E-04,3.9511E-04,1.2565E-04,1.4851E-05/
      DATA (FRACREFA(IG, 9),IG=1,16) /
     &2.3656E-01,2.4978E-01,1.8663E-01,1.1488E-01,7.5115E-02,5.4237E-02,
     &3.9867E-02,2.4057E-02,1.4028E-02,1.4857E-03,1.1903E-03,9.7374E-04,
     &6.6182E-04,3.9611E-04,1.2648E-04,1.4837E-05/

C     Calculate reference ratio to be used in calculation of Planck
C     fraction in lower atmosphere.

C     P = 387. mb , 245K
      REFRAT_PLANCK_A = chi_std(1,18)/chi_std(2,18)


C     Calculate reference ratio to be used in calculation of minor species absorption
C     P =  387 mbar, 245K 
      REFRAT_M_A = CHI_STD(1,18)/CHI_STD(2,18)

C     Compute the optical depth by interpolating in ln(pressure), 
C     temperature,and appropriate species.  Below LAYTROP, the water
C     vapor self-continuum and foreign continuum is interpolated 
C     (in temperature) separately.  

      HVRTAU = '$Revision: 22602 $'

      DO 2500 LAY = 1, NLAYERS


         SPECCOMB = COLH2O(LAY) + RAT_H2OCO2(LAY)*COLCO2(LAY)
         SPECPARM = COLH2O(LAY)/SPECCOMB
         IF (SPECPARM .GE. ONEMINUS) SPECPARM = ONEMINUS
         SPECMULT = 8.*(SPECPARM)
         JS = 1 + INT(SPECMULT)
         FS = AMOD(SPECMULT,1.0)

         SPECCOMB1 = COLH2O(LAY) + RAT_H2OCO2_1(LAY)*COLCO2(LAY)
         SPECPARM1 = COLH2O(LAY)/SPECCOMB1
         IF (SPECPARM1 .GE. ONEMINUS) SPECPARM1 = ONEMINUS
         SPECMULT1 = 8.*(SPECPARM1)
         JS1 = 1 + INT(SPECMULT1)
         FS1 = AMOD(SPECMULT1,1.0)

         SPECCOMB_PLANCK = COLH2O(LAY)+REFRAT_PLANCK_A*COLCO2(LAY)
         SPECPARM_PLANCK = COLH2O(LAY)/SPECCOMB_PLANCK
         IF (SPECPARM_PLANCK .GE. ONEMINUS) SPECPARM_PLANCK=ONEMINUS
         SPECMULT_PLANCK = 8.*SPECPARM_PLANCK
         JPL= 1 + INT(SPECMULT_PLANCK)
         FPL = AMOD(SPECMULT_PLANCK,1.0)

         SPECCOMB_MCH4 = COLH2O(LAY) + REFRAT_M_A*COLCO2(LAY)
         SPECPARM_MCH4 = COLH2O(LAY)/SPECCOMB_MCH4
         IF (SPECPARM_MCH4 .GE. ONEMINUS) SPECPARM_MCH4 = ONEMINUS
         SPECMULT_MCH4 = 8.*SPECPARM_MCH4
         JMCH4 = 1 + INT(SPECMULT_MCH4)
         FMCH4 = AMOD(SPECMULT_MCH4,1.0)

         IND0 = ((JP(LAY)-1)*NTEMP+(JT(LAY)-1))*NSPA(15) + JS
         IND1 = (JP(LAY)*NTEMP+(JT1(LAY)-1))*NSPA(15) + JS1
         INDS = INDSELF(LAY)
         INDF = INDFOR(LAY)
         INDM = INDMINOR(LAY)
         IF (SPECPARM .LT. 0.125 .AND. SPECPARM1 .LT. 0.125) THEN
            P = FS - 1
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC000 = FK0*FAC00(LAY)
            FAC100 = FK1*FAC00(LAY)
            FAC200 = FK2*FAC00(LAY)
            FAC010 = FK0*FAC10(LAY)
            FAC110 = FK1*FAC10(LAY)
            FAC210 = FK2*FAC10(LAY)

            P = FS1 - 1
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC001 = FK0*FAC01(LAY)
            FAC101 = FK1*FAC01(LAY)
            FAC201 = FK2*FAC01(LAY)
            FAC011 = FK0*FAC11(LAY)
            FAC111 = FK1*FAC11(LAY)
            FAC211 = FK2*FAC11(LAY)

            DO 2000 IG = 1, NG(15)
               TAUSELF = SELFFAC(LAY)* (SELFREF(INDS,IG) +
     &              SELFFRAC(LAY)  *
     &              (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
               TAUFOR =  FORFAC(LAY) * (FORREF(INDF,IG) +
     &              FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &              FORREF(INDF,IG))) 

               CH4M1 = KA_MCH4(JMCH4,INDM,IG) + FMCH4*
     &              (KA_MCH4(JMCH4+1,INDM,IG)-
     &              KA_MCH4(JMCH4,INDM,IG))
               CH4M2 = KA_MCH4(JMCH4,INDM+1,IG) + FMCH4*
     &              (KA_MCH4(JMCH4+1,INDM+1,IG)-
     &              KA_MCH4(JMCH4,INDM+1,IG))
               ABSCH4 = CH4M1 + MINORFRAC(LAY)
     &              * (CH4M2 - CH4M1)
               TAUCH4 = COLCH4(LAY)*ABSCH4

               TAUG(LAY,IG) = SPECCOMB *
     &              (FAC000 * ABSA(IND0,IG) +
     &              FAC100 * ABSA(IND0+1,IG) +
     &              FAC200 * ABSA(IND0+2,IG) +
     &              FAC010 * ABSA(IND0+9,IG) +
     &              FAC110 * ABSA(IND0+10,IG) +
     &              FAC210 * ABSA(IND0+11,IG))
     &              + SPECCOMB1 *
     &              (FAC001 * ABSA(IND1,IG) + 
     &              FAC101 * ABSA(IND1+1,IG) +
     &              FAC201 * ABSA(IND1+2,IG) +
     &              FAC011 * ABSA(IND1+9,IG) +
     &              FAC111 * ABSA(IND1+10,IG) +
     &              FAC211 * ABSA(IND1+11,IG)) 
C      &              + TAUSELF + TAUFOR + TAUCH4 !NJE
     &              +TAUCH4

                taug(lay,ig) = taug(lay,ig)
               if( h2o_sb == 1 ) then 
                  taug(lay,ig) = taug(lay,ig) + TAUSELF
                endif
               if( h2o_for == 1 ) then 
                  taug(lay,ig) = taug(lay,ig) + TAUFOR
                endif

               FRACS(LAY,IG) = FRACREFA(IG,JPL) + FPL *
     &              (FRACREFA(IG,JPL+1)-FRACREFA(IG,JPL))
 2000          CONTINUE
      ELSE IF (SPECPARM .GT. 0.875 .AND. SPECPARM1 .GT. 0.875) THEN
            P = -FS 
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC000 = FK0*FAC00(LAY)
            FAC100 = FK1*FAC00(LAY)
            FAC200 = FK2*FAC00(LAY)
            FAC010 = FK0*FAC10(LAY)
            FAC110 = FK1*FAC10(LAY)
            FAC210 = FK2*FAC10(LAY)

            P = -FS1 
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC001 = FK0*FAC01(LAY)
            FAC101 = FK1*FAC01(LAY)
            FAC201 = FK2*FAC01(LAY)
            FAC011 = FK0*FAC11(LAY)
            FAC111 = FK1*FAC11(LAY)
            FAC211 = FK2*FAC11(LAY)
            DO 2010 IG = 1, NG(15)
               TAUSELF = SELFFAC(LAY)* (SELFREF(INDS,IG) +
     &              SELFFRAC(LAY)  *
     &              (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
               TAUFOR = FORFAC(LAY) * (FORREF(INDF,IG) +
     &              FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &              FORREF(INDF,IG))) 
               CH4M1 = KA_MCH4(JMCH4,INDM,IG) + FMCH4*
     &              (KA_MCH4(JMCH4+1,INDM,IG)-
     &              KA_MCH4(JMCH4,INDM,IG))
               CH4M2 = KA_MCH4(JMCH4,INDM+1,IG) + FMCH4*
     &              (KA_MCH4(JMCH4+1,INDM+1,IG)-
     &              KA_MCH4(JMCH4,INDM+1,IG))
               ABSCH4 = CH4M1 + MINORFRAC(LAY)
     &              * (CH4M2 - CH4M1)
               TAUCH4 = COLCH4(LAY)*ABSCH4
               TAUG(LAY,IG) = SPECCOMB * 
     &              (FAC200 * ABSA(IND0-1,IG) +
     &              FAC100 * ABSA(IND0,IG) +
     &              FAC000 * ABSA(IND0+1,IG) +
     &              FAC210 * ABSA(IND0+8,IG) +
     &              FAC110 * ABSA(IND0+9,IG) +
     &              FAC010 * ABSA(IND0+10,IG))
     &              + SPECCOMB1 *
     &              (FAC201 * ABSA(IND1-1,IG) +
     &              FAC101 * ABSA(IND1,IG) +
     &              FAC001 * ABSA(IND1+1,IG) +
     &              FAC211 * ABSA(IND1+8,IG) +
     &              FAC111 * ABSA(IND1+9,IG) +
     &              FAC011 * ABSA(IND1+10,IG))
C      &              + TAUSELF + TAUFOR + TAUCH4 !NJE
     &              +TAUCH4 

              taug(lay,ig) = taug(lay,ig)
               if( h2o_sb == 1 ) then 
                  taug(lay,ig) = taug(lay,ig) + TAUSELF
                endif
               if( h2o_for == 1 ) then 
                  taug(lay,ig) = taug(lay,ig) + TAUFOR
                endif

               FRACS(LAY,IG) = FRACREFA(IG,JPL) + FPL *
     &              (FRACREFA(IG,JPL+1)-FRACREFA(IG,JPL))
 2010           CONTINUE
       ELSE
            FAC000 = (1. - FS) * FAC00(LAY)
            FAC010 = (1. - FS) * FAC10(LAY)
            FAC100 = FS * FAC00(LAY)
            FAC110 = FS * FAC10(LAY)

            FAC001 = (1. - FS1) * FAC01(LAY)
            FAC011 = (1. - FS1) * FAC11(LAY)
            FAC101 = FS1 * FAC01(LAY)
            FAC111 = FS1 * FAC11(LAY)
            DO 2020 IG = 1, NG(15)
               TAUSELF = SELFFAC(LAY)* (SELFREF(INDS,IG) +
     &              SELFFRAC(LAY)  *
     &              (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
               TAUFOR = FORFAC(LAY) * (FORREF(INDF,IG) +
     &              FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &              FORREF(INDF,IG))) 
               CH4M1 = KA_MCH4(JMCH4,INDM,IG) + FMCH4*
     &              (KA_MCH4(JMCH4+1,INDM,IG)-
     &              KA_MCH4(JMCH4,INDM,IG))
               CH4M2 = KA_MCH4(JMCH4,INDM+1,IG) + FMCH4*
     &              (KA_MCH4(JMCH4+1,INDM+1,IG)-
     &              KA_MCH4(JMCH4,INDM+1,IG))
               ABSCH4 = CH4M1 + MINORFRAC(LAY)
     &              * (CH4M2 - CH4M1)
               TAUCH4 = COLCH4(LAY)*ABSCH4

               TAUG(LAY,IG) = SPECCOMB * 
     &              (FAC000 * ABSA(IND0,IG) +
     &              FAC100 * ABSA(IND0+1,IG) +
     &              FAC010 * ABSA(IND0+9,IG) +
     &              FAC110 * ABSA(IND0+10,IG))
     &              + SPECCOMB1 *
     &              (FAC001 * ABSA(IND1,IG) + 
     &              FAC101 * ABSA(IND1+1,IG) +
     &              FAC011 * ABSA(IND1+9,IG) +
     &              FAC111 * ABSA(IND1+10,IG)) 
C      &              + TAUSELF + TAUFOR + TAUCH4
     &              +TAUCH4

              taug(lay,ig) = taug(lay,ig)
               if( h2o_sb == 1 ) then 
                  taug(lay,ig) = taug(lay,ig) + TAUSELF
                endif
               if( h2o_for == 1 ) then 
                  taug(lay,ig) = taug(lay,ig) + TAUFOR
                endif              

               FRACS(LAY,IG) = FRACREFA(IG,JPL) + FPL *
     &              (FRACREFA(IG,JPL+1)-FRACREFA(IG,JPL))

 2020          CONTINUE
      ENDIF

 2500 CONTINUE

      RETURN
      END
C----------------------------------------------------------------------------

      SUBROUTINE TAUGB16
      use variables !NJE

      return
      end
 
C----------------------------------------------------------------------------
