C     path:      $Source$
C     author:    $Author: cadyp $
C     revision:  $Revision: 11240 $
C     created:   $Date: 2010-12-28 15:23:15 -0500 (Tue, 28 Dec 2010) $
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

      Subroutine CLDPROP

C     Purpose:  Compute the cloud optical depth(s) for each cloudy
C               layer.

      PARAMETER (MXLAY=203)
      PARAMETER (NBANDS = 16)
      PARAMETER (N_ICE_M = 300)

      COMMON /CONTROL/  NUMANGS, ISCAT, NSTR, 
     &                  IOUT, ISTART, IEND, ICLD
      COMMON /IFIL/     IRD,IPR,IPU,IDUM(15)
      COMMON /BANDS/     WAVENUM1(NBANDS),WAVENUM2(NBANDS),
     &                   DELWAVE(NBANDS)
      COMMON /PROFILE/  NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
     &                  PZ(0:MXLAY),TZ(0:MXLAY)
      COMMON /CLOUDIN/   ICD,ICLDATM,INFLAG,
     &     CLDDAT1(MXLAY),CLDDAT2(MXLAY),
     &     ICEFLAG,LIQFLAG,CLDDAT3(MXLAY),CLDDAT4(MXLAY)
      COMMON /CLOUDDAT/ NCBANDS,CLDFRAC(MXLAY),
     &     TAUCLOUD(MXLAY,NBANDS), 
     &     SSACLOUD(MXLAY,NBANDS),
     &     XMOM(0:16,MXLAY,NBANDS),
     &     TAUTOT(NBANDS)

      COMMON / CLDOPTPROPS /ABSCLD1, ABSLIQ0,ABSICE0(2), 
     &     d_eff_m(n_ice_m),ABSICE4( N_ICE_M,15),
     &     ABSLIQ1(58,15),
     &     EXTICE2(43,16),SSAICE2(43,16),ASYICE2(43,16),
     &     EXTICE3(46,16),SSAICE3(46,16),ASYICE3(46,16)

      COMMON /CVRCLD/    HNAMCLD,HVRCLD

      CHARACTER*18       HNAMCLD,HVRCLD

      DIMENSION ABSCOICE(NBANDS), ABSCOLIQ(NBANDS)
      DIMENSION EXTCOLIQ(NBANDS),SSACOLIQ(NBANDS),GLIQ(NBANDS)
      DIMENSION EXTCOICE(NBANDS),SSACOICE(NBANDS),GICE(NBANDS)
      DIMENSION FORWLIQ(NBANDS),FORWICE(NBANDS)

      DIMENSION IPAT(16,0:2)

      DATA EPS /1.E-6/

      DATA IPAT /1,1,1,1,1,1,1,1,1, 1, 1, 1, 1, 1, 1, 1,
     &           1,2,3,3,3,4,4,4,5, 5, 5, 5, 5, 5, 5, 5,
     &           1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16/

C     Explanation of the method for each value of INFLAG.  Values of
C     0  for INFLAG do not distingish being liquid and ice clouds.
C     INFLAG = 2 does distinguish between liquid and ice clouds, and
C     requires further user input to specify the method to be used to 
C     compute the aborption due to each
C     INFLAG = 0:  For each cloudy layer, the cloud fraction and (gray)
C                  cloud optical depth (and if ISCAT=2) single-scattering albedo,
C                  and phase-function moments are input.	 
C     INFLAG = 2:  For each cloudy layer, the cloud fraction, cloud 
C                  water path (g/m2), and cloud ice fraction are input.
C     ICEFLAG = 4:  The ice effective diameter (microns) is input and the
C                     optical properties due to ice clouds are computed from
C                     values provided by Mitchell
C                 :reference_1 = "Mitchell, D.L., A. Macke and Y. Liu, 1996: Modeling cirrus clouds. Part II:
C Treatment of radiative properties. J. Atmos. Sci., 53, 2967-2988" ;
C                 :reference_2 = "Mitchell, D.L., 2000: Parameterization of the Mie extinction and absorption
C coefficients for water clouds. J. Atmos. Sci., 57, 1311-1326" ;
C                 :reference_3 = "Mitchell, D.L., 2002: Effective diameter in radiation transfer: General
C definition, applications, and limitations. J. Atmos. Sci., 59, 2330-2346" ;
C                 :reference_4 = "Mitchell et al. 2006: Testing and comparing the modified anomalous diffraction
C approximation. J. Atmos. Sci., 63, 2948-2962" ;

C       LIQFLAG = 1:  The water droplet effective radius (microns) is input 
C                     and the optical depths due to water clouds are computed 
C                     as in Hu and Stamnes, J., Clim., 6, 728-742, (1993).
C                     The values for absorption coefficients appropriate for
C                     the spectral bands in RRTM have been obtained for a 
C                     range of effective radii by an averaging procedure 
C                     based on the work of J. Pinto (private communication).
C                     Linear interpolation is used to get the absorption 
C                     coefficients for the input effective radius.
C     INFLAG = 10:  For each cloudy layer, the cloud fraction and 
C                   cloud optical depth (and if ISCAT=2) single-scattering albedo,
C                   and phase-function moments are input for each band.	 

C Note: ISCAT=2 is currently not implemented for PRRTM


      HVRCLD = '$Revision: 11240 $'

      ICLDATM = 0
      NCBANDS = 1

      IF (ISCAT .EQ. 0 .OR. ISCAT .EQ. 1) THEN
         
         IF (INFLAG .EQ. 1) THEN 
            print *, 'INFLAG=1 not a valid option for PRRTM'
         ELSEIF (INFLAG .EQ. 2) THEN
            WRITE(ICD,8900) '  LAY','  BAND  ',
     &           '  WVN1  ','  WVN2  ',
     &           '   ICE OD    ','   LIQ OD    ',
     &           '  TOT ABS OD ','  TOT ABS OD '

            WRITE(ICD,8901) '  BY BLOCK  '
         ENDIF

         NPRELAY = 0
         DO 3000 LAY = 1, NLAYERS
            IF (CLDFRAC(LAY) .GE. EPS) THEN
               IF (CLDFRAC(LAY) .NE. 1.0 .AND. ISCAT .GE. 1) 
     &              STOP 'CLDFRAC MUST BE 1 WHEN USING DISORT'
               ICLDATM = 1

C           Ice clouds and water clouds combined.
               IF(INFLAG .EQ. 1) THEN
                  CWP = CLDDAT1(LAY)
		  NCBANDS = 15
                  DO 1000 IB = 1, NCBANDS                  
                     TAUCLOUD(LAY,IB) = ABSCLD1 * CWP
                     TAUTOT(IB) = TAUTOT(IB) + 
     &                    TAUCLOUD(LAY,IB)                     
                     WRITE(ICD,9002) LAY,IB,WAVENUM1(IB),WAVENUM2(IB),
     &                    TAUCLOUD(LAY,IB),
     &                    TAUTOT(IB)
 1000             CONTINUE
                  
C           Separate treatement of ice clouds and water clouds.
               ELSEIF(INFLAG .EQ. 2) THEN
                  CWP = CLDDAT1(LAY)
                  FICE = CLDDAT2(LAY)
                  RADICE = CLDDAT3(LAY)
                  
C              Calculation of absorption coefficients due to ice clouds.
                  IF (FICE .EQ. 0.0) THEN
                     ABSCOICE(1) = 0.0
                     ICEPAT = 0
                  ELSEIF (ICEFLAG .EQ. 4) THEN
                     diam_ice = RADICE
                     IF (diam_ice .LT. 0.62 .OR. diam_ice .GT. 777.58) 
     &               STOP     'ICE RADIUS OUT OF BOUNDS'
                     NCBANDS = 15
                     if (diam_ice.le.d_eff_m(1) ) then
                        k_d_eff = 2
                        wd = 1.0
                        onemwd = 0.
                     elseif (diam_ice.ge.d_eff_m(n_ice_m)) then 
                        k_d_eff = n_ice_m
                        wd =  0.
                        onemwd = 1.0
                     else
			do ir=2,n_ice_m
                          k_d_eff = ir
                          if (d_eff_m(ir) > diam_ice) exit
                        enddo
                     wd = (d_eff_m(k_d_eff)-diam_ice)/(d_eff_m(k_d_eff)-
     &                      d_eff_m(k_d_eff-1))
                     onemwd = 1.0-wd
                     endif
C Convert absorption from m2/kg to m2/g
                     DO 2200 IB = 1, NCBANDS
                        ABSCOICE(IB) = FICE * 1.0e-3 * 
     &                       (wd*ABSICE4(k_d_eff-1,IB) +
     &                        onemwd*ABSICE4(k_d_eff,IB) ) 
 2200                CONTINUE
                     ICEPAT = 2
                  ENDIF
                  
C     Calculation of absorption coefficients due to water clouds.
                  FLIQ = 1. - FICE
                  IF (FLIQ .EQ. 0.0) THEN
                     ABSCOLIQ(1) = 0.0
                     LIQPAT = 0
                     IF (ICEPAT .EQ. 1) ICEPAT = 2
                  ELSEIF (LIQFLAG .EQ. 1) THEN
                     RADLIQ = CLDDAT4(LAY)
                     IF (RADLIQ .LT. 2.5 .OR. RADLIQ .GT. 60.) STOP
     &                    'LIQUID EFFECTIVE RADIUS OUT OF BOUNDS'
                     INDEX = RADLIQ - 1.5
                     IF (INDEX .EQ. 58) INDEX = 57
                     IF (INDEX .EQ. 0) INDEX = 1
                     FINT = RADLIQ - 1.5 - INDEX
                     NCBANDS = 15
                     DO 2400 IB = 1, NCBANDS
                        ABSCOLIQ(IB) = FLIQ * 
     &                       (ABSLIQ1(INDEX,IB) + FINT *
     &                       (ABSLIQ1(INDEX+1,IB) - 
     &                       (ABSLIQ1(INDEX,IB))))
 2400                CONTINUE
                     LIQPAT = 2
                  ENDIF
                  
C KEEP RUNNING TOTAL OF OPTICAL DEPTH OF EACH CLOUD.  WHEN
C A CLEAR LAYER SEPARATES CLOUDY LAYERS, RESET COUNTER

                  IF (LAY .NE. NPRELAY+1) TAUTOT(:) = 0.0

                  DO 2800 IB = 1, NCBANDS
                     TAULIQ = CWP * ABSCOLIQ(IPAT(IB,LIQPAT))
                     TAUICE = CWP * ABSCOICE(IPAT(IB,ICEPAT))
                     TAUCLOUD(LAY,IB) = TAULIQ +
     &                    TAUICE
                     TAUTOT(IB) = TAUTOT(IB) + 
     &                    TAUCLOUD(LAY,IB)
                     IF (NCBANDS .EQ. 1) THEN
                         DO 2600 IBPR =1, 16                             
                             WRITE(ICD,9000) LAY,IBPR,
     &                           WAVENUM1(IBPR),WAVENUM2(IBPR),
     &                           TAUICE,TAULIQ,
     &                           TAUCLOUD(LAY,IB),
     &                           TAUTOT(IB)
 2600                    CONTINUE
                     ELSEIF (NCBANDS .EQ. 5) THEN
                         DO 2700 IBPR = 1, 16
                             IF (IPAT(IBPR,1) .EQ. IB) THEN
                                 WRITE(ICD,9000) LAY,IBPR,
     &                               WAVENUM1(IBPR),WAVENUM2(IBPR),
     &                               TAUICE,TAULIQ,
     &                               TAUCLOUD(LAY,IB),
     &                               TAUTOT(IB)
                             ENDIF
 2700                    CONTINUE
                     ELSE
                         WRITE(ICD,9000) LAY,IB,
     &                       WAVENUM1(IB),WAVENUM2(IB),
     &                       TAUICE,TAULIQ,
     &                       TAUCLOUD(LAY,IB),
     &                       TAUTOT(IB)
                     ENDIF
 2800             CONTINUE
               ENDIF
               NPRELAY = LAY
            ENDIF
 3000    CONTINUE

      ELSEIF (ISCAT .EQ. 2) THEN
         print *, 'ISCAT=2 not implemented for PRRTM'
         STOP
         WRITE(ICD,8902) '  LAY','  BAND  ',
     &        '  WVN1  ','  WVN2  ',
     &        '   ICE OD    ','   LIQ OD   ',
     &        '  TOT EXT OD ',' TOT EXT OD ',
     &        '  ICE SSA    ',' LIQ  SSA   ',
     &        '     ICE     ','    LIQ     '
         WRITE(ICD,8903) '   BY BLOCK  ','  ASYM  FAC  ',
     &        '  ASYM  FAC  '

         NPRELAY = 0
         DO 6000 LAY = 1, NLAYERS

            IF (CLDFRAC(LAY) .GE. EPS) THEN
               IF (CLDFRAC(LAY) .NE. 1.0 ) 
     &              STOP 'CLDFRAC MUST BE 1 WHEN ISCAT=1,2'
               ICLDATM = 1
               
               IF (INFLAG .EQ. 1) THEN
                  PRINT*,'INFLAG OPTION 1 NOT ALLOWED WITH ISCAT=2'
                  STOP
               ENDIF
               
C     Separate treatment of ice clouds and water clouds.
               IF (INFLAG .EQ. 2) THEN
                  CWP = CLDDAT1(LAY)
                  FICE = CLDDAT2(LAY)
                  RADICE = CLDDAT3(LAY)
                  
                  IF (FICE .NE. 1.0) THEN
                     PRINT*,'ICE FRACTION MUST BE SET TO 1.0'
                     STOP
                  ENDIF
                  
                  IF (ICEFLAG .NE. 2 .AND. ICEFLAG .NE. 3) THEN
                     PRINT*,'ICEFLAG MUST BE SET TO 2 OR 3 FOR
     &                    SCATTERING CALCULATIONS',iceflag
                     STOP
                  ENDIF

                  IF (ICEFLAG .EQ. 2) THEN
                     IF (RADICE .LT. 5.0 .OR. RADICE .GT. 131.)
     &                    STOP 'ICE RADIUS OUT OF BOUNDS'
                     NCBANDS = 15
                     FACTOR = (RADICE - 2.)/3.
                     INDEX = INT(FACTOR)
                     IF (INDEX .EQ. 43) INDEX = 42
                     FINT = FACTOR - FLOAT(INDEX)
                     DO 5200 IB = 1, NCBANDS
                        EXTCOICE(IB) = FICE * 
     &                       (EXTICE2(INDEX,IB) + FINT *
     &                       (EXTICE2(INDEX+1,IB) - 
     &                       EXTICE2(INDEX,IB)))
                        SSACOICE(IB) = SSAICE2(INDEX,IB) + FINT *
     &                       (SSAICE2(INDEX+1,IB) - 
     &                       SSAICE2(INDEX,IB))
                        GICE(IB) = ASYICE2(INDEX,IB) + FINT *
     &                       (ASYICE2(INDEX+1,IB) - 
     &                       ASYICE2(INDEX,IB))
 5200                CONTINUE
                     ICEPAT = 2
                  ENDIF

                  IF (ICEFLAG .EQ. 3) THEN
                     IF (RADICE .LT. 5.0 .OR. RADICE .GT. 140.)
     &                    STOP 'ICE RADIUS OUT OF BOUNDS'
                     NCBANDS = 15
                     FACTOR = (RADICE-2)/3.
                     INDEX = INT(FACTOR)
                     IF (INDEX .EQ. 46) INDEX = 45
                     FINT = FACTOR - FLOAT(INDEX)
                     DO 5300 IB = 1, NCBANDS
                        EXTCOICE(IB) = FICE * 
     &                       (EXTICE3(INDEX,IB) + FINT *
     &                       (EXTICE3(INDEX+1,IB) - 
     &                       EXTICE3(INDEX,IB)))
                        SSACOICE(IB) = SSAICE3(INDEX,IB) + FINT *
     &                       (SSAICE3(INDEX+1,IB) - 
     &                       SSAICE3(INDEX,IB))
                        GICE(IB) = ASYICE3(INDEX,IB) + FINT *
     &                       (ASYICE3(INDEX+1,IB) - 
     &                       ASYICE3(INDEX,IB))
 5300                CONTINUE
                     ICEPAT = 2
                  ENDIF
                  
C     Calculation of water clouds
                FLIQ = 1. - FICE
                IF (FLIQ .NE. 0.0) STOP
     &               'LIQUID PARTICLES NOT PERMITTED, NO SCATTERING 
     &               PROPS AVAILABLE'

C KEEP RUNNING TOTAL OF OPTICAL DEPTH OF EACH CLOUD.  WHEN
C A CLEAR LAYER SEPARATES CLOUDY LAYERS, RESET COUNTER

                  IF (LAY .NE. NPRELAY+1) TAUTOT(:) = 0.0

c     Calculation of optical properties of ice-only cloud.

                DO 5800 IB = 1, NCBANDS
                   TAUCLOUD(LAY,IB) = CWP * EXTCOICE(IB)
                   SSACLOUD(LAY,IB) = SSACOICE(IB)
                   XMOM(0,LAY,IB) = 1.0
                   DO 5750 ISTR = 1, NSTR
                      XMOM(ISTR,LAY,IB) = GICE(IB)**ISTR
 5750              CONTINUE
                   TAUTOT(IB) = TAUTOT(IB) + 
     &                  TAUCLOUD(LAY,IB)
                   WRITE(ICD,9001) LAY,IB,WAVENUM1(IB),WAVENUM2(IB),
     &                  TAUCLOUD(LAY,IB),0.0,TAUCLOUD(LAY,IB),
     &                  TAUTOT(IB),SSACLOUD(LAY,IB),
     &                  0.0,GICE(IB),0.0
 5800           CONTINUE
             ENDIF   
             NPRELAY = LAY
          ENDIF
 6000  CONTINUE
      ENDIF

 7998 FORMAT(A5,A6,A8,A8,2(A13))
 7999 FORMAT(2X,3X,1X,4X,1X,7X,1X,7X,1X,
     &     1(12X,1X),A13)
 8900 FORMAT(A5,A6,A8,A8,4(A13))
 8901 FORMAT(2X,3X,1X,4X,1X,7X,1X,7X,1X,
     &     12X,1X,12X,1X,12X,1X,A13)
 8902 FORMAT(A5,A6,A8,A8,8(A13))
 8903 FORMAT((2X,3X),(1X,4X,1X),2(7X,1X),
     &     3(12X,1X),A13,2(12X,1X),A13,A13)

 9000 FORMAT((2X,I3),(1X,I4,1X),2(F7.1,1X),4(E12.5,1X))
 9001 FORMAT(2X,I3,1X,I4,1X,2(F7.1,1X),4(E12.5,1X),
     &     4(E12.5,1X))
 9002 FORMAT((2X,I3),(1X,I4,1X),2(F7.1,1X),2(E12.5,1X))

      RETURN
      END

      BLOCK DATA CLDPARAMS
      PARAMETER (n_ice_m=300)

      COMMON / CLDOPTPROPS /ABSCLD1, ABSLIQ0,ABSICE0(2), 
     &     d_eff_m(n_ice_m),ABSICE4( N_ICE_M,15),
     &     ABSLIQ1(58,15),
     &     EXTICE2(43,16),SSAICE2(43,16),ASYICE2(43,16),
     &     EXTICE3(46,16),SSAICE3(46,16),ASYICE3(46,16)


C     ABSCLDn is the liquid water absorption coefficient (m2/g). 
C     For INFLAG = 1.
      DATA ABSCLD1 /0.0602410/
C  
C     Everything below is for INFLAG = 2.

C     ABSICEn(J,IB) are the parameters needed to compute the liquid water 
C     absorption coefficient in spectral region IB for ICEFLAG=n.  The 
C     units of ABSICEn(1,IB) are m2/g and ABSICEn(2,IB) has units 
C     (microns (m2/g)).
C     For ICEFLAG = 0.
      DATA ABSICE0 /0.005,  1.0/

C     For ICEFLAG = 4.  In each band, the absorption
C     coefficients are listed for a range of effective radii from  0.62
C     to 777.58.0 microns 
C     Mitchell, D.L., A. Macke and Y. Liu, 1996: Modeling cirrus clouds. Part II: Treatment of radiative properties. J. Atmos. Sci., 53, 2967-2988
C     ABSORPTION UNITS (ABS COEF): [(m^2)/(kg )]



C     For LIQFLAG = 0.
      DATA ABSLIQ0 /0.0903614/

C     For LIQFLAG = 1.  In each band, the absorption
C     coefficients are listed for a range of effective radii from 2.5
C     to 59.5 microns in increments of 1.0 micron.

      DATA (ABSLIQ1(I, 1),I=1,58) /
c     BAND  1
     & 2.99061E-02, 3.21129E-02, 3.36937E-02, 3.53329E-02, 3.69950E-02,
     & 3.86469E-02, 4.02559E-02, 4.17892E-02, 4.32128E-02, 4.44916E-02,
     & 4.68165E-02, 4.75607E-02, 4.76033E-02, 4.73015E-02, 4.68306E-02,
     & 4.62777E-02, 4.56860E-02, 4.50753E-02, 4.44532E-02, 4.38205E-02,
     & 4.31741E-02, 4.25088E-02, 4.18180E-02, 4.10941E-02, 4.03291E-02,
     & 3.95142E-02, 3.86405E-02, 3.76982E-02, 3.64729E-02, 3.54831E-02,
     & 3.45437E-02, 3.36482E-02, 3.27919E-02, 3.19705E-02, 3.11805E-02,
     & 3.04191E-02, 2.96838E-02, 2.89724E-02, 2.82832E-02, 2.76144E-02,
     & 2.69647E-02, 2.63327E-02, 2.57174E-02, 2.51176E-02, 2.45325E-02,
     & 2.39613E-02, 2.34031E-02, 2.28573E-02, 2.23231E-02, 2.18001E-02,
     & 2.12876E-02, 2.07852E-02, 2.02924E-02, 1.98087E-02, 1.93337E-02,
     & 1.88671E-02, 1.84084E-02, 1.79574E-02/
      DATA (ABSLIQ1(I, 2),I=1,58) /
c     BAND  2
     & 2.89979E-02, 8.54713E-02, 9.65387E-02, 9.66821E-02, 9.47110E-02,
     & 9.23208E-02, 8.99218E-02, 8.76115E-02, 8.53955E-02, 8.32473E-02,
     & 8.01803E-02, 7.51209E-02, 7.07937E-02, 6.70051E-02, 6.36276E-02,
     & 6.05728E-02, 5.77775E-02, 5.51949E-02, 5.27894E-02, 5.05336E-02,
     & 4.84056E-02, 4.63880E-02, 4.44665E-02, 4.26295E-02, 4.08670E-02,
     & 3.91710E-02, 3.75343E-02, 3.59510E-02, 3.48701E-02, 3.36285E-02,
     & 3.24607E-02, 3.13604E-02, 3.03216E-02, 2.93395E-02, 2.84094E-02,
     & 2.75273E-02, 2.66895E-02, 2.58928E-02, 2.51342E-02, 2.44110E-02,
     & 2.37207E-02, 2.30613E-02, 2.24305E-02, 2.18266E-02, 2.12480E-02,
     & 2.06929E-02, 2.01601E-02, 1.96482E-02, 1.91560E-02, 1.86824E-02,
     & 1.82263E-02, 1.77868E-02, 1.73630E-02, 1.69540E-02, 1.65592E-02,
     & 1.61778E-02, 1.58090E-02, 1.54524E-02/
      DATA (ABSLIQ1(I, 3),I=1,58) /
c     BAND  3
     & 2.19486E-01, 1.80687E-01, 1.59150E-01, 1.44731E-01, 1.33703E-01,
     & 1.24355E-01, 1.15756E-01, 1.07318E-01, 9.86119E-02, 8.92739E-02,
     & 8.34911E-02, 7.70773E-02, 7.15240E-02, 6.66615E-02, 6.23641E-02,
     & 5.85359E-02, 5.51020E-02, 5.20032E-02, 4.91916E-02, 4.66283E-02,
     & 4.42813E-02, 4.21236E-02, 4.01330E-02, 3.82905E-02, 3.65797E-02,
     & 3.49869E-02, 3.35002E-02, 3.21090E-02, 3.08957E-02, 2.97601E-02,
     & 2.86966E-02, 2.76984E-02, 2.67599E-02, 2.58758E-02, 2.50416E-02,
     & 2.42532E-02, 2.35070E-02, 2.27997E-02, 2.21284E-02, 2.14904E-02,
     & 2.08834E-02, 2.03051E-02, 1.97536E-02, 1.92271E-02, 1.87239E-02,
     & 1.82425E-02, 1.77816E-02, 1.73399E-02, 1.69162E-02, 1.65094E-02,
     & 1.61187E-02, 1.57430E-02, 1.53815E-02, 1.50334E-02, 1.46981E-02,
     & 1.43748E-02, 1.40628E-02, 1.37617E-02/
      DATA (ABSLIQ1(I, 4),I=1,58) /
c     BAND  4
     & 2.95174E-01, 2.34765E-01, 1.98038E-01, 1.72114E-01, 1.52083E-01,
     & 1.35654E-01, 1.21613E-01, 1.09252E-01, 9.81263E-02, 8.79448E-02,
     & 8.12566E-02, 7.44563E-02, 6.86374E-02, 6.36042E-02, 5.92094E-02,
     & 5.53402E-02, 5.19087E-02, 4.88455E-02, 4.60951E-02, 4.36124E-02,
     & 4.13607E-02, 3.93096E-02, 3.74338E-02, 3.57119E-02, 3.41261E-02,
     & 3.26610E-02, 3.13036E-02, 3.00425E-02, 2.88497E-02, 2.78077E-02,
     & 2.68317E-02, 2.59158E-02, 2.50545E-02, 2.42430E-02, 2.34772E-02,
     & 2.27533E-02, 2.20679E-02, 2.14181E-02, 2.08011E-02, 2.02145E-02,
     & 1.96561E-02, 1.91239E-02, 1.86161E-02, 1.81311E-02, 1.76673E-02,
     & 1.72234E-02, 1.67981E-02, 1.63903E-02, 1.59989E-02, 1.56230E-02,
     & 1.52615E-02, 1.49138E-02, 1.45791E-02, 1.42565E-02, 1.39455E-02,
     & 1.36455E-02, 1.33559E-02, 1.30761E-02/
      DATA (ABSLIQ1(I, 5),I=1,58) /
c     BAND  5
     & 3.00925E-01, 2.36949E-01, 1.96947E-01, 1.68692E-01, 1.47190E-01,
     & 1.29986E-01, 1.15719E-01, 1.03568E-01, 9.30028E-02, 8.36658E-02,
     & 7.71075E-02, 7.07002E-02, 6.52284E-02, 6.05024E-02, 5.63801E-02,
     & 5.27534E-02, 4.95384E-02, 4.66690E-02, 4.40925E-02, 4.17664E-02,
     & 3.96559E-02, 3.77326E-02, 3.59727E-02, 3.43561E-02, 3.28662E-02,
     & 3.14885E-02, 3.02110E-02, 2.90231E-02, 2.78948E-02, 2.69109E-02,
     & 2.59884E-02, 2.51217E-02, 2.43058E-02, 2.35364E-02, 2.28096E-02,
     & 2.21218E-02, 2.14700E-02, 2.08515E-02, 2.02636E-02, 1.97041E-02,
     & 1.91711E-02, 1.86625E-02, 1.81769E-02, 1.77126E-02, 1.72683E-02,
     & 1.68426E-02, 1.64344E-02, 1.60427E-02, 1.56664E-02, 1.53046E-02,
     & 1.49565E-02, 1.46214E-02, 1.42985E-02, 1.39871E-02, 1.36866E-02,
     & 1.33965E-02, 1.31162E-02, 1.28453E-02/
      DATA (ABSLIQ1(I, 6),I=1,58) /
c     BAND  6
     & 2.64691E-01, 2.12018E-01, 1.78009E-01, 1.53539E-01, 1.34721E-01,
     & 1.19580E-01, 1.06996E-01, 9.62772E-02, 8.69710E-02, 7.87670E-02,
     & 7.29272E-02, 6.70920E-02, 6.20977E-02, 5.77732E-02, 5.39910E-02,
     & 5.06538E-02, 4.76866E-02, 4.50301E-02, 4.26374E-02, 4.04704E-02,
     & 3.84981E-02, 3.66948E-02, 3.50394E-02, 3.35141E-02, 3.21038E-02,
     & 3.07957E-02, 2.95788E-02, 2.84438E-02, 2.73790E-02, 2.64390E-02,
     & 2.55565E-02, 2.47263E-02, 2.39437E-02, 2.32047E-02, 2.25056E-02,
     & 2.18433E-02, 2.12149E-02, 2.06177E-02, 2.00495E-02, 1.95081E-02,
     & 1.89917E-02, 1.84984E-02, 1.80269E-02, 1.75755E-02, 1.71431E-02,
     & 1.67283E-02, 1.63303E-02, 1.59478E-02, 1.55801E-02, 1.52262E-02,
     & 1.48853E-02, 1.45568E-02, 1.42400E-02, 1.39342E-02, 1.36388E-02,
     & 1.33533E-02, 1.30773E-02, 1.28102E-02/
      DATA (ABSLIQ1(I, 7),I=1,58) /
c     BAND  7
     & 9.08506E-02, 1.01239E-01, 9.40300E-02, 8.65253E-02, 8.05098E-02,
     & 7.55990E-02, 7.14174E-02, 6.77038E-02, 6.42726E-02, 6.09809E-02,
     & 5.71685E-02, 5.39597E-02, 5.11468E-02, 4.86383E-02, 4.63722E-02,
     & 4.43046E-02, 4.24029E-02, 4.06421E-02, 3.90028E-02, 3.74693E-02,
     & 3.60290E-02, 3.46712E-02, 3.33871E-02, 3.21692E-02, 3.10111E-02,
     & 2.99073E-02, 2.88530E-02, 2.78439E-02, 2.69610E-02, 2.61906E-02,
     & 2.54452E-02, 2.47263E-02, 2.40345E-02, 2.33699E-02, 2.27320E-02,
     & 2.21201E-02, 2.15333E-02, 2.09707E-02, 2.04313E-02, 1.99138E-02,
     & 1.94173E-02, 1.89408E-02, 1.84832E-02, 1.80435E-02, 1.76209E-02,
     & 1.72144E-02, 1.68233E-02, 1.64467E-02, 1.60838E-02, 1.57341E-02,
     & 1.53967E-02, 1.50712E-02, 1.47568E-02, 1.44531E-02, 1.41596E-02,
     & 1.38757E-02, 1.36010E-02, 1.33350E-02/
      DATA (ABSLIQ1(I, 8),I=1,58) /
c     BAND  8
     & 1.13356E-01, 7.75511E-02, 6.95359E-02, 6.63696E-02, 6.44460E-02,
     & 6.28438E-02, 6.11750E-02, 5.91589E-02, 5.64787E-02, 5.26915E-02,
     & 4.95621E-02, 4.84945E-02, 4.70016E-02, 4.53210E-02, 4.35851E-02,
     & 4.18675E-02, 4.02078E-02, 3.86263E-02, 3.71317E-02, 3.57262E-02,
     & 3.44081E-02, 3.31737E-02, 3.20183E-02, 3.09367E-02, 2.99236E-02,
     & 2.89738E-02, 2.80825E-02, 2.72450E-02, 2.62009E-02, 2.54798E-02,
     & 2.47768E-02, 2.40947E-02, 2.34352E-02, 2.27992E-02, 2.21869E-02,
     & 2.15981E-02, 2.10323E-02, 2.04889E-02, 1.99670E-02, 1.94658E-02,
     & 1.89844E-02, 1.85220E-02, 1.80776E-02, 1.76503E-02, 1.72394E-02,
     & 1.68441E-02, 1.64635E-02, 1.60969E-02, 1.57436E-02, 1.54030E-02,
     & 1.50744E-02, 1.47573E-02, 1.44510E-02, 1.41552E-02, 1.38692E-02,
     & 1.35926E-02, 1.33249E-02, 1.30658E-02/
      DATA (ABSLIQ1(I, 9),I=1,58) /
c     BAND  9
     & 1.03555E-01, 6.71391E-02, 6.25829E-02, 6.13137E-02, 6.04141E-02,
     & 5.93269E-02, 5.78448E-02, 5.58020E-02, 5.29702E-02, 4.89536E-02,
     & 4.60045E-02, 4.55945E-02, 4.45514E-02, 4.31883E-02, 4.16855E-02,
     & 4.01466E-02, 3.86309E-02, 3.71712E-02, 3.57846E-02, 3.44786E-02,
     & 3.32550E-02, 3.21123E-02, 3.10473E-02, 3.00555E-02, 2.91321E-02,
     & 2.82723E-02, 2.74713E-02, 2.67246E-02, 2.56779E-02, 2.49769E-02,
     & 2.42975E-02, 2.36410E-02, 2.30078E-02, 2.23977E-02, 2.18106E-02,
     & 2.12458E-02, 2.07028E-02, 2.01806E-02, 1.96786E-02, 1.91958E-02,
     & 1.87315E-02, 1.82848E-02, 1.78548E-02, 1.74409E-02, 1.70421E-02,
     & 1.66578E-02, 1.62874E-02, 1.59300E-02, 1.55852E-02, 1.52523E-02,
     & 1.49307E-02, 1.46199E-02, 1.43194E-02, 1.40287E-02, 1.37474E-02,
     & 1.34750E-02, 1.32111E-02, 1.29554E-02/
      DATA (ABSLIQ1(I,10),I=1,58) /
c     BAND 10
     & 7.10016E-02, 6.91408E-02, 6.73308E-02, 6.55993E-02, 6.39260E-02,
     & 6.22378E-02, 6.03977E-02, 5.81895E-02, 5.52966E-02, 5.12752E-02,
     & 4.84390E-02, 4.73228E-02, 4.57790E-02, 4.40534E-02, 4.22825E-02,
     & 4.05412E-02, 3.88696E-02, 3.72872E-02, 3.58018E-02, 3.44147E-02,
     & 3.31230E-02, 3.19220E-02, 3.08059E-02, 2.97688E-02, 2.88044E-02,
     & 2.79070E-02, 2.70711E-02, 2.62916E-02, 2.52462E-02, 2.45359E-02,
     & 2.38524E-02, 2.31955E-02, 2.25647E-02, 2.19592E-02, 2.13782E-02,
     & 2.08207E-02, 2.02857E-02, 1.97721E-02, 1.92791E-02, 1.88055E-02,
     & 1.83504E-02, 1.79130E-02, 1.74922E-02, 1.70873E-02, 1.66975E-02,
     & 1.63219E-02, 1.59600E-02, 1.56109E-02, 1.52742E-02, 1.49490E-02,
     & 1.46350E-02, 1.43316E-02, 1.40382E-02, 1.37544E-02, 1.34797E-02,
     & 1.32138E-02, 1.29561E-02, 1.27064E-02/
      DATA (ABSLIQ1(I,11),I=1,58) /
c     BAND 11
     & 1.40942E-01, 1.26465E-01, 1.15557E-01, 1.06654E-01, 9.89432E-02,
     & 9.19288E-02, 8.52715E-02, 7.87196E-02, 7.20748E-02, 6.51725E-02,
     & 6.21602E-02, 5.83267E-02, 5.47462E-02, 5.14681E-02, 4.84919E-02,
     & 4.57977E-02, 4.33586E-02, 4.11471E-02, 3.91375E-02, 3.73061E-02,
     & 3.56324E-02, 3.40980E-02, 3.26872E-02, 3.13861E-02, 3.01826E-02,
     & 2.90663E-02, 2.80279E-02, 2.70593E-02, 2.59818E-02, 2.51507E-02,
     & 2.43641E-02, 2.36191E-02, 2.29128E-02, 2.22423E-02, 2.16053E-02,
     & 2.09994E-02, 2.04224E-02, 1.98724E-02, 1.93476E-02, 1.88463E-02,
     & 1.83670E-02, 1.79084E-02, 1.74690E-02, 1.70478E-02, 1.66436E-02,
     & 1.62554E-02, 1.58823E-02, 1.55235E-02, 1.51780E-02, 1.48452E-02,
     & 1.45245E-02, 1.42150E-02, 1.39163E-02, 1.36278E-02, 1.33490E-02,
     & 1.30794E-02, 1.28185E-02, 1.25660E-02/
      DATA (ABSLIQ1(I,12),I=1,58) /
c     BAND 12
     & 3.67457E-02, 3.85757E-02, 3.97071E-02, 4.02679E-02, 4.02976E-02,
     & 3.98096E-02, 3.88080E-02, 3.72933E-02, 3.52648E-02, 3.27210E-02,
     & 3.31136E-02, 3.21018E-02, 3.10994E-02, 3.01244E-02, 2.91861E-02,
     & 2.82892E-02, 2.74356E-02, 2.66249E-02, 2.58562E-02, 2.51277E-02,
     & 2.44374E-02, 2.37831E-02, 2.31627E-02, 2.25740E-02, 2.20150E-02,
     & 2.14838E-02, 2.09784E-02, 2.04974E-02, 1.99116E-02, 1.95101E-02,
     & 1.91157E-02, 1.87292E-02, 1.83511E-02, 1.79818E-02, 1.76215E-02,
     & 1.72703E-02, 1.69281E-02, 1.65950E-02, 1.62707E-02, 1.59552E-02,
     & 1.56482E-02, 1.53496E-02, 1.50591E-02, 1.47765E-02, 1.45016E-02,
     & 1.42342E-02, 1.39739E-02, 1.37207E-02, 1.34742E-02, 1.32342E-02,
     & 1.30005E-02, 1.27729E-02, 1.25511E-02, 1.23351E-02, 1.21246E-02,
     & 1.19193E-02, 1.17192E-02, 1.15240E-02/
      DATA (ABSLIQ1(I,13),I=1,58) /
c     BAND 13
     & 3.11868E-02, 4.48357E-02, 4.90224E-02, 4.96406E-02, 4.86806E-02,
     & 4.69610E-02, 4.48630E-02, 4.25795E-02, 4.02138E-02, 3.78236E-02,
     & 3.74266E-02, 3.60384E-02, 3.47074E-02, 3.34434E-02, 3.22499E-02,
     & 3.11264E-02, 3.00704E-02, 2.90784E-02, 2.81463E-02, 2.72702E-02,
     & 2.64460E-02, 2.56698E-02, 2.49381E-02, 2.42475E-02, 2.35948E-02,
     & 2.29774E-02, 2.23925E-02, 2.18379E-02, 2.11793E-02, 2.07076E-02,
     & 2.02470E-02, 1.97981E-02, 1.93613E-02, 1.89367E-02, 1.85243E-02,
     & 1.81240E-02, 1.77356E-02, 1.73588E-02, 1.69935E-02, 1.66392E-02,
     & 1.62956E-02, 1.59624E-02, 1.56393E-02, 1.53259E-02, 1.50219E-02,
     & 1.47268E-02, 1.44404E-02, 1.41624E-02, 1.38925E-02, 1.36302E-02,
     & 1.33755E-02, 1.31278E-02, 1.28871E-02, 1.26530E-02, 1.24253E-02,
     & 1.22038E-02, 1.19881E-02, 1.17782E-02/
      DATA (ABSLIQ1(I,14),I=1,58) /
c     BAND 14
     & 7.85078E-03, 2.49775E-02, 2.89588E-02, 2.93774E-02, 2.86988E-02,
     & 2.77048E-02, 2.66706E-02, 2.56950E-02, 2.48100E-02, 2.40212E-02,
     & 2.34529E-02, 2.27775E-02, 2.21332E-02, 2.15219E-02, 2.09440E-02,
     & 2.03983E-02, 1.98834E-02, 1.93976E-02, 1.89390E-02, 1.85056E-02,
     & 1.80958E-02, 1.77079E-02, 1.73402E-02, 1.69913E-02, 1.66599E-02,
     & 1.63448E-02, 1.60448E-02, 1.57589E-02, 1.53705E-02, 1.51354E-02,
     & 1.49010E-02, 1.46681E-02, 1.44374E-02, 1.42094E-02, 1.39845E-02,
     & 1.37631E-02, 1.35452E-02, 1.33311E-02, 1.31209E-02, 1.29147E-02,
     & 1.27124E-02, 1.25141E-02, 1.23198E-02, 1.21295E-02, 1.19431E-02,
     & 1.17605E-02, 1.15817E-02, 1.14067E-02, 1.12353E-02, 1.10674E-02,
     & 1.09031E-02, 1.07421E-02, 1.05845E-02, 1.04301E-02, 1.02789E-02,
     & 1.01308E-02, 9.98566E-03, 9.84346E-03/
      DATA (ABSLIQ1(I,15),I=1,58) /
c     BAND 15
     & 7.78660E-02, 1.19597E-01, 1.05497E-01, 9.06054E-02, 7.88828E-02,
     & 6.99360E-02, 6.30540E-02, 5.76763E-02, 5.34044E-02, 4.99583E-02,
     & 4.57177E-02, 4.28364E-02, 4.03049E-02, 3.80642E-02, 3.60676E-02,
     & 3.42775E-02, 3.26637E-02, 3.12016E-02, 2.98708E-02, 2.86544E-02,
     & 2.75383E-02, 2.65105E-02, 2.55610E-02, 2.46810E-02, 2.38633E-02,
     & 2.31013E-02, 2.23895E-02, 2.17231E-02, 2.09880E-02, 2.04280E-02,
     & 1.98931E-02, 1.93817E-02, 1.88925E-02, 1.84241E-02, 1.79754E-02,
     & 1.75451E-02, 1.71322E-02, 1.67357E-02, 1.63547E-02, 1.59883E-02,
     & 1.56356E-02, 1.52959E-02, 1.49686E-02, 1.46529E-02, 1.43483E-02,
     & 1.40541E-02, 1.37698E-02, 1.34950E-02, 1.32292E-02, 1.29719E-02,
     & 1.27227E-02, 1.24813E-02, 1.22472E-02, 1.20201E-02, 1.17998E-02,
     & 1.15859E-02, 1.13782E-02, 1.11763E-02/

C     Spherical Ice Particle Parameterization
C     EXTINCTION UNITS (EXT COEF/IWC): [(m^-1)/(g m^-3)]
      DATA (EXTICE2(I,1),I=1,43) /
C    BAND 1
     &1.024095e-01,8.941945e-02,8.047835e-02,7.378238e-02,6.846938e-02,
     &6.408248e-02,6.035569e-02,5.712182e-02,5.426928e-02,5.172011e-02,
     &4.941777e-02,4.731999e-02,4.539440e-02,4.361567e-02,4.196359e-02,
     &4.042182e-02,3.897697e-02,3.761791e-02,3.633530e-02,3.512124e-02,
     &3.396896e-02,3.287266e-02,3.182731e-02,3.082849e-02,2.987237e-02,
     &2.895553e-02,2.807497e-02,2.722800e-02,2.641222e-02,2.562548e-02,
     &2.486583e-02,2.413152e-02,2.342097e-02,2.273271e-02,2.206544e-02,
     &2.141794e-02,2.078911e-02,2.017793e-02,1.958346e-02,1.900483e-02,
     &1.844125e-02,1.789198e-02,1.735631e-02/
      DATA (EXTICE2(I,2),I=1,43) /
C    BAND 2
     &1.565614e-01,1.317435e-01,1.158132e-01,1.041811e-01,9.506821e-02,
     &8.760376e-02,8.129839e-02,7.585057e-02,7.106177e-02,6.679452e-02,
     &6.294987e-02,5.945427e-02,5.625161e-02,5.329815e-02,5.055916e-02,
     &4.800661e-02,4.561755e-02,4.337298e-02,4.125699e-02,3.925613e-02,
     &3.735892e-02,3.555550e-02,3.383735e-02,3.219702e-02,3.062800e-02,
     &2.912455e-02,2.768159e-02,2.629459e-02,2.495953e-02,2.367277e-02,
     &2.243106e-02,2.123145e-02,2.007127e-02,1.894810e-02,1.785972e-02,
     &1.680411e-02,1.577943e-02,1.478396e-02,1.381614e-02,1.287453e-02,
     &1.195779e-02,1.106467e-02,1.019404e-02/
      DATA (EXTICE2(I,3),I=1,43) /
C    BAND 3
     &2.968902e-01,2.116436e-01,1.670425e-01,1.388669e-01,1.191456e-01,
     &1.044167e-01,9.291190e-02,8.362511e-02,7.593776e-02,6.944668e-02,
     &6.387677e-02,5.903330e-02,5.477418e-02,5.099304e-02,4.760855e-02,
     &4.455733e-02,4.178919e-02,3.926383e-02,3.694847e-02,3.481615e-02,
     &3.284447e-02,3.101465e-02,2.931082e-02,2.771947e-02,2.622900e-02,
     &2.482941e-02,2.351202e-02,2.226926e-02,2.109450e-02,1.998187e-02,
     &1.892621e-02,1.792292e-02,1.696789e-02,1.605746e-02,1.518833e-02,
     &1.435755e-02,1.356242e-02,1.280053e-02,1.206967e-02,1.136783e-02,
     &1.069318e-02,1.004405e-02,9.418899e-03/
      DATA (EXTICE2(I,4),I=1,43) /
C    BAND 4
     &3.846587e-01,2.479981e-01,1.835502e-01,1.457513e-01,1.207923e-01,
     &1.030276e-01,8.971024e-02,7.933926e-02,7.102391e-02,6.420128e-02,
     &5.849786e-02,5.365580e-02,4.949125e-02,4.586950e-02,4.268954e-02,
     &3.987410e-02,3.736306e-02,3.510892e-02,3.307361e-02,3.122631e-02,
     &2.954173e-02,2.799898e-02,2.658062e-02,2.527196e-02,2.406056e-02,
     &2.293580e-02,2.188858e-02,2.091103e-02,1.999630e-02,1.913844e-02,
     &1.833223e-02,1.757307e-02,1.685689e-02,1.618009e-02,1.553944e-02,
     &1.493210e-02,1.435548e-02,1.380728e-02,1.328541e-02,1.278799e-02,
     &1.231331e-02,1.185983e-02,1.142613e-02/
      DATA (EXTICE2(I,5),I=1,43) /
C    BAND 5
     &3.836455e-01,2.424484e-01,1.772772e-01,1.396330e-01,1.150705e-01,
     &9.775863e-02,8.488845e-02,7.493832e-02,6.701161e-02,6.054541e-02,
     &5.516830e-02,5.062520e-02,4.673511e-02,4.336599e-02,4.041923e-02,
     &3.781967e-02,3.550906e-02,3.344150e-02,3.158034e-02,2.989599e-02,
     &2.836426e-02,2.696518e-02,2.568214e-02,2.450121e-02,2.341058e-02,
     &2.240023e-02,2.146155e-02,2.058714e-02,1.977057e-02,1.900626e-02,
     &1.828931e-02,1.761544e-02,1.698084e-02,1.638217e-02,1.581643e-02,
     &1.528097e-02,1.477340e-02,1.429159e-02,1.383361e-02,1.339773e-02,
     &1.298237e-02,1.258611e-02,1.220766e-02/
      DATA (EXTICE2(I,6),I=1,43) /
C    BAND 6
     &2.463346e-01,1.731218e-01,1.359502e-01,1.129018e-01,9.697616e-02,
     &8.519425e-02,7.605765e-02,6.872406e-02,6.268077e-02,5.759639e-02,
     &5.324640e-02,4.947285e-02,4.616112e-02,4.322579e-02,4.060185e-02,
     &3.823881e-02,3.609686e-02,3.414409e-02,3.235464e-02,3.070730e-02,
     &2.918447e-02,2.777146e-02,2.645583e-02,2.522706e-02,2.407611e-02,
     &2.299519e-02,2.197757e-02,2.101737e-02,2.010945e-02,1.924928e-02,
     &1.843286e-02,1.765663e-02,1.691744e-02,1.621245e-02,1.553914e-02,
     &1.489521e-02,1.427861e-02,1.368747e-02,1.312010e-02,1.257496e-02,
     &1.205064e-02,1.154585e-02,1.105943e-02/
      DATA (EXTICE2(I,7),I=1,43) /
C    BAND 7
     &2.332844e-01,1.753972e-01,1.432347e-01,1.220420e-01,1.067171e-01,
     &9.496334e-02,8.557395e-02,7.784593e-02,7.133827e-02,6.575843e-02,
     &6.090362e-02,5.662826e-02,5.282472e-02,4.941149e-02,4.632555e-02,
     &4.351730e-02,4.094709e-02,3.858278e-02,3.639802e-02,3.437096e-02,
     &3.248331e-02,3.071964e-02,2.906680e-02,2.751353e-02,2.605011e-02,
     &2.466811e-02,2.336015e-02,2.211980e-02,2.094133e-02,1.981973e-02,
     &1.875050e-02,1.772963e-02,1.675356e-02,1.581905e-02,1.492320e-02,
     &1.406338e-02,1.323721e-02,1.244252e-02,1.167733e-02,1.093985e-02,
     &1.022841e-02,9.541504e-03,8.877727e-03/
      DATA (EXTICE2(I,8),I=1,43) /
C    BAND 8
     &3.163919e-01,2.164187e-01,1.662972e-01,1.356062e-01,1.146515e-01,
     &9.932297e-02,8.756212e-02,7.821700e-02,7.058960e-02,6.423083e-02,
     &5.883780e-02,5.419837e-02,5.015924e-02,4.660673e-02,4.345460e-02,
     &4.063621e-02,3.809917e-02,3.580171e-02,3.371008e-02,3.179670e-02,
     &3.003880e-02,2.841738e-02,2.691650e-02,2.552264e-02,2.422427e-02,
     &2.301149e-02,2.187574e-02,2.080962e-02,1.980662e-02,1.886106e-02,
     &1.796794e-02,1.712282e-02,1.632176e-02,1.556126e-02,1.483817e-02,
     &1.414969e-02,1.349327e-02,1.286663e-02,1.226769e-02,1.169458e-02,
     &1.114559e-02,1.061916e-02,1.011386e-02/
      DATA (EXTICE2(I,9),I=1,43) /
C    BAND 9
     &3.936120e-01,2.470191e-01,1.797882e-01,1.411182e-01,1.159652e-01,
     &9.828077e-02,8.516020e-02,7.503375e-02,6.697841e-02,6.041569e-02,
     &5.496451e-02,5.036352e-02,4.642748e-02,4.302142e-02,4.004462e-02,
     &3.742042e-02,3.508943e-02,3.300489e-02,3.112952e-02,2.943321e-02,
     &2.789136e-02,2.648370e-02,2.519337e-02,2.400622e-02,2.291028e-02,
     &2.189540e-02,2.095285e-02,2.007513e-02,1.925574e-02,1.848903e-02,
     &1.777005e-02,1.709446e-02,1.645842e-02,1.585854e-02,1.529181e-02,
     &1.475555e-02,1.424733e-02,1.376502e-02,1.330667e-02,1.287053e-02,
     &1.245501e-02,1.205867e-02,1.168022e-02/
      DATA (EXTICE2(I,10),I=1,43) /
C    BAND 10
     &4.254883e-01,2.572374e-01,1.828685e-01,1.411963e-01,1.146376e-01,
     &9.627516e-02,8.284435e-02,7.260635e-02,6.455138e-02,5.805351e-02,
     &5.270441e-02,4.822656e-02,4.442483e-02,4.115808e-02,3.832175e-02,
     &3.583677e-02,3.364222e-02,3.169044e-02,2.994361e-02,2.837136e-02,
     &2.694900e-02,2.565626e-02,2.447637e-02,2.339529e-02,2.240124e-02,
     &2.148420e-02,2.063566e-02,1.984827e-02,1.911574e-02,1.843257e-02,
     &1.779398e-02,1.719579e-02,1.663432e-02,1.610633e-02,1.560893e-02,
     &1.513957e-02,1.469597e-02,1.427608e-02,1.387807e-02,1.350029e-02,
     &1.314125e-02,1.279961e-02,1.247414e-02/
      DATA (EXTICE2(I,11),I=1,43) /
C    BAND 11
     &4.395814e-01,2.605045e-01,1.829231e-01,1.400554e-01,1.130281e-01,
     &9.450563e-02,8.105812e-02,7.087317e-02,6.290529e-02,5.651022e-02,
     &5.126987e-02,4.690139e-02,4.320677e-02,4.004337e-02,3.730588e-02,
     &3.491492e-02,3.280956e-02,3.094223e-02,2.927534e-02,2.777872e-02,
     &2.642795e-02,2.520302e-02,2.408740e-02,2.306730e-02,2.213115e-02,
     &2.126916e-02,2.047299e-02,1.973550e-02,1.905053e-02,1.841276e-02,
     &1.781754e-02,1.726083e-02,1.673906e-02,1.624911e-02,1.578818e-02,
     &1.535382e-02,1.494384e-02,1.455627e-02,1.418935e-02,1.384150e-02,
     &1.351130e-02,1.319746e-02,1.289882e-02/
      DATA (EXTICE2(I,12),I=1,43) /
C    BAND 12
     &5.011192e-01,2.791706e-01,1.886085e-01,1.406223e-01,1.113326e-01,
     &9.178395e-02,7.790561e-02,6.759635e-02,5.966820e-02,5.340185e-02,
     &4.833776e-02,4.416940e-02,4.068494e-02,3.773356e-02,3.520513e-02,
     &3.301747e-02,3.110810e-02,2.942870e-02,2.794136e-02,2.661592e-02,
     &2.542815e-02,2.435834e-02,2.339031e-02,2.251066e-02,2.170821e-02,
     &2.097355e-02,2.029873e-02,1.967696e-02,1.910244e-02,1.857014e-02,
     &1.807575e-02,1.761548e-02,1.718604e-02,1.678454e-02,1.640844e-02,
     &1.605547e-02,1.572364e-02,1.541118e-02,1.511649e-02,1.483815e-02,
     &1.457488e-02,1.432555e-02,1.408910e-02/
      DATA (EXTICE2(I,13),I=1,43) /
C    BAND 13
     &4.923918e-01,2.742336e-01,1.852807e-01,1.381734e-01,1.094335e-01,
     &9.025950e-02,7.665201e-02,6.654720e-02,5.877857e-02,5.263995e-02,
     &4.768033e-02,4.359892e-02,4.018790e-02,3.729932e-02,3.482517e-02,
     &3.268489e-02,3.081720e-02,2.917473e-02,2.772033e-02,2.642446e-02,
     &2.526337e-02,2.421773e-02,2.327170e-02,2.241216e-02,2.162816e-02,
     &2.091049e-02,2.025135e-02,1.964410e-02,1.908305e-02,1.856331e-02,
     &1.808062e-02,1.763130e-02,1.721212e-02,1.682026e-02,1.645321e-02,
     &1.610878e-02,1.578501e-02,1.548016e-02,1.519267e-02,1.492116e-02,
     &1.466438e-02,1.442121e-02,1.419062e-02/
      DATA (EXTICE2(I,14),I=1,43) /
C    BAND 14
     &4.869234e-01,2.718080e-01,1.839165e-01,1.373055e-01,1.088377e-01,
     &8.982848e-02,7.632799e-02,6.629626e-02,5.857948e-02,5.247880e-02,
     &4.754761e-02,4.348792e-02,4.009378e-02,3.721849e-02,3.475494e-02,
     &3.262317e-02,3.076238e-02,2.912556e-02,2.767578e-02,2.638373e-02,
     &2.522579e-02,2.418276e-02,2.323891e-02,2.238118e-02,2.159868e-02,
     &2.088225e-02,2.022413e-02,1.961773e-02,1.905737e-02,1.853819e-02,
     &1.805595e-02,1.760698e-02,1.718807e-02,1.679640e-02,1.642949e-02,
     &1.608514e-02,1.576140e-02,1.545655e-02,1.516903e-02,1.489746e-02,
     &1.464059e-02,1.439730e-02,1.416658e-02/
      DATA (EXTICE2(I,15),I=1,43) /
C    BAND 15
     &4.712597e-01,2.658087e-01,1.810589e-01,1.358108e-01,1.080308e-01,
     &8.940151e-02,7.612280e-02,6.622477e-02,5.858962e-02,5.253839e-02,
     &4.763608e-02,4.359183e-02,4.020412e-02,3.732921e-02,3.486192e-02,
     &3.272361e-02,3.085441e-02,2.920791e-02,2.774769e-02,2.644471e-02,
     &2.527560e-02,2.422135e-02,2.326630e-02,2.239751e-02,2.160413e-02,
     &2.087706e-02,2.020855e-02,1.959204e-02,1.902185e-02,1.849312e-02,
     &1.800163e-02,1.754369e-02,1.711609e-02,1.671601e-02,1.634095e-02,
     &1.598871e-02,1.565735e-02,1.534510e-02,1.505042e-02,1.477192e-02,
     &1.450833e-02,1.425854e-02,1.402151e-02/
      DATA (EXTICE2(I,16),I=1,43) /
C    BAND 16
     &4.101824e-01,2.435514e-01,1.713697e-01,1.314865e-01,1.063406e-01,
     &8.910701e-02,7.659480e-02,6.711784e-02,5.970353e-02,5.375249e-02,
     &4.887577e-02,4.481025e-02,4.137171e-02,3.842744e-02,3.587948e-02,
     &3.365396e-02,3.169419e-02,2.995593e-02,2.840419e-02,2.701091e-02,
     &2.575336e-02,2.461293e-02,2.357423e-02,2.262443e-02,2.175276e-02,
     &2.095012e-02,2.020875e-02,1.952199e-02,1.888412e-02,1.829018e-02,
     &1.773586e-02,1.721738e-02,1.673144e-02,1.627510e-02,1.584579e-02,
     &1.544122e-02,1.505934e-02,1.469833e-02,1.435654e-02,1.403251e-02,
     &1.372492e-02,1.343255e-02,1.315433e-02/
C     SINGLE-SCATTERING ALBEDO: Unitless
      DATA (SSAICE2(I,1),I=1,43) /
C    BAND 1
     &2.384497e-01,2.909285e-01,3.267788e-01,3.540130e-01,3.759747e-01,
     &3.943837e-01,4.102388e-01,4.241709e-01,4.366037e-01,4.478358e-01,
     &4.580852e-01,4.675163e-01,4.762559e-01,4.844041e-01,4.920411e-01,
     &4.992324e-01,5.060323e-01,5.124859e-01,5.186316e-01,5.245023e-01,
     &5.301260e-01,5.355275e-01,5.407281e-01,5.457469e-01,5.506008e-01,
     &5.553048e-01,5.598726e-01,5.643165e-01,5.686477e-01,5.728765e-01,
     &5.770124e-01,5.810643e-01,5.850404e-01,5.889484e-01,5.927958e-01,
     &5.965895e-01,6.003362e-01,6.040423e-01,6.077142e-01,6.113580e-01,
     &6.149798e-01,6.185855e-01,6.221813e-01/
      DATA (SSAICE2(I,2),I=1,43) /
C    BAND 2
     &8.221222e-01,7.943076e-01,7.738458e-01,7.572275e-01,7.430030e-01,
     &7.304254e-01,7.190571e-01,7.086179e-01,6.989172e-01,6.898189e-01,
     &6.812222e-01,6.730501e-01,6.652427e-01,6.577517e-01,6.505382e-01,
     &6.435701e-01,6.368202e-01,6.302659e-01,6.238876e-01,6.176684e-01,
     &6.115934e-01,6.056497e-01,5.998256e-01,5.941108e-01,5.884958e-01,
     &5.829720e-01,5.775313e-01,5.721662e-01,5.668698e-01,5.616352e-01,
     &5.564557e-01,5.513249e-01,5.462362e-01,5.411828e-01,5.361576e-01,
     &5.311530e-01,5.261609e-01,5.211719e-01,5.161756e-01,5.111597e-01,
     &5.061093e-01,5.010065e-01,4.958288e-01/
      DATA (SSAICE2(I,3),I=1,43) /
C    BAND 3
     &6.411479e-01,6.217354e-01,6.080982e-01,5.975189e-01,5.888491e-01,
     &5.814910e-01,5.750920e-01,5.694265e-01,5.643407e-01,5.597253e-01,
     &5.554996e-01,5.516022e-01,5.479855e-01,5.446115e-01,5.414499e-01,
     &5.384755e-01,5.356676e-01,5.330090e-01,5.304849e-01,5.280827e-01,
     &5.257918e-01,5.236027e-01,5.215074e-01,5.194988e-01,5.175707e-01,
     &5.157175e-01,5.139344e-01,5.122171e-01,5.105619e-01,5.089654e-01,
     &5.074245e-01,5.059368e-01,5.044998e-01,5.031117e-01,5.017706e-01,
     &5.004753e-01,4.992246e-01,4.980175e-01,4.968537e-01,4.957327e-01,
     &4.946547e-01,4.936202e-01,4.926301e-01/
      DATA (SSAICE2(I,4),I=1,43) /
C    BAND 4
     &5.308658e-01,5.286308e-01,5.270808e-01,5.259007e-01,5.249552e-01,
     &5.241732e-01,5.235120e-01,5.229445e-01,5.224518e-01,5.220204e-01,
     &5.216405e-01,5.213044e-01,5.210062e-01,5.207411e-01,5.205054e-01,
     &5.202958e-01,5.201099e-01,5.199454e-01,5.198004e-01,5.196735e-01,
     &5.195632e-01,5.194684e-01,5.193881e-01,5.193214e-01,5.192675e-01,
     &5.192259e-01,5.191958e-01,5.191767e-01,5.191683e-01,5.191702e-01,
     &5.191819e-01,5.192032e-01,5.192338e-01,5.192736e-01,5.193222e-01,
     &5.193797e-01,5.194458e-01,5.195204e-01,5.196036e-01,5.196951e-01,
     &5.197952e-01,5.199036e-01,5.200205e-01/
      DATA (SSAICE2(I,5),I=1,43) /
C    BAND 5
     &4.443291e-01,4.591129e-01,4.693524e-01,4.772410e-01,4.836841e-01,
     &4.891456e-01,4.938956e-01,4.981055e-01,5.018908e-01,5.053336e-01,
     &5.084939e-01,5.114170e-01,5.141381e-01,5.166852e-01,5.190805e-01,
     &5.213424e-01,5.234860e-01,5.255239e-01,5.274670e-01,5.293243e-01,
     &5.311038e-01,5.328121e-01,5.344554e-01,5.360388e-01,5.375669e-01,
     &5.390438e-01,5.404732e-01,5.418583e-01,5.432020e-01,5.445071e-01,
     &5.457758e-01,5.470105e-01,5.482129e-01,5.493851e-01,5.505286e-01,
     &5.516449e-01,5.527355e-01,5.538017e-01,5.548446e-01,5.558653e-01,
     &5.568650e-01,5.578445e-01,5.588047e-01/
      DATA (SSAICE2(I,6),I=1,43) /
C    BAND 6
     &3.723265e-01,3.996998e-01,4.181439e-01,4.320349e-01,4.431625e-01,
     &4.524352e-01,4.603781e-01,4.673218e-01,4.734879e-01,4.790323e-01,
     &4.840687e-01,4.886825e-01,4.929397e-01,4.968922e-01,5.005815e-01,
     &5.040415e-01,5.073000e-01,5.103805e-01,5.133026e-01,5.160831e-01,
     &5.187363e-01,5.212749e-01,5.237097e-01,5.260504e-01,5.283055e-01,
     &5.304825e-01,5.325884e-01,5.346293e-01,5.366107e-01,5.385378e-01,
     &5.404153e-01,5.422477e-01,5.440388e-01,5.457926e-01,5.475127e-01,
     &5.492024e-01,5.508652e-01,5.525042e-01,5.541224e-01,5.557230e-01,
     &5.573090e-01,5.588835e-01,5.604495e-01/
      DATA (SSAICE2(I,7),I=1,43) /
C    BAND 7
     &6.749288e-01,6.475680e-01,6.291381e-01,6.152112e-01,6.040010e-01,
     &5.946084e-01,5.865167e-01,5.794023e-01,5.730486e-01,5.673037e-01,
     &5.620572e-01,5.572259e-01,5.527461e-01,5.485674e-01,5.446498e-01,
     &5.409604e-01,5.374723e-01,5.341632e-01,5.310141e-01,5.280089e-01,
     &5.251338e-01,5.223770e-01,5.197280e-01,5.171778e-01,5.147184e-01,
     &5.123427e-01,5.100445e-01,5.078180e-01,5.056582e-01,5.035605e-01,
     &5.015208e-01,4.995353e-01,4.976005e-01,4.957133e-01,4.938707e-01,
     &4.920701e-01,4.903088e-01,4.885844e-01,4.868948e-01,4.852377e-01,
     &4.836110e-01,4.820127e-01,4.804408e-01/
      DATA (SSAICE2(I,8),I=1,43) /
C    BAND 8
     &7.148414e-01,6.801247e-01,6.565983e-01,6.387794e-01,6.244318e-01,
     &6.124207e-01,6.020902e-01,5.930271e-01,5.849538e-01,5.776752e-01,
     &5.710486e-01,5.649666e-01,5.593463e-01,5.541225e-01,5.492428e-01,
     &5.446646e-01,5.403528e-01,5.362779e-01,5.324151e-01,5.287436e-01,
     &5.252450e-01,5.219039e-01,5.187066e-01,5.156411e-01,5.126970e-01,
     &5.098650e-01,5.071367e-01,5.045048e-01,5.019627e-01,4.995044e-01,
     &4.971244e-01,4.948179e-01,4.925804e-01,4.904078e-01,4.882964e-01,
     &4.862427e-01,4.842437e-01,4.822964e-01,4.803980e-01,4.785462e-01,
     &4.767386e-01,4.749730e-01,4.732474e-01/
      DATA (SSAICE2(I,9),I=1,43) /
C    BAND 9
     &6.712279e-01,6.442293e-01,6.257659e-01,6.116928e-01,6.003067e-01,
     &5.907386e-01,5.824839e-01,5.752233e-01,5.687419e-01,5.628879e-01,
     &5.575502e-01,5.526449e-01,5.481072e-01,5.438859e-01,5.399399e-01,
     &5.362355e-01,5.327451e-01,5.294455e-01,5.263172e-01,5.233434e-01,
     &5.205098e-01,5.178041e-01,5.152154e-01,5.127343e-01,5.103524e-01,
     &5.080623e-01,5.058575e-01,5.037321e-01,5.016807e-01,4.996987e-01,
     &4.977816e-01,4.959258e-01,4.941275e-01,4.923836e-01,4.906910e-01,
     &4.890472e-01,4.874496e-01,4.858958e-01,4.843839e-01,4.829118e-01,
     &4.814778e-01,4.800801e-01,4.787172e-01/
      DATA (SSAICE2(I,10),I=1,43) /
C    BAND 10
     &6.254590e-01,6.055970e-01,5.921137e-01,5.818892e-01,5.736492e-01,
     &5.667460e-01,5.608054e-01,5.555910e-01,5.509442e-01,5.467533e-01,
     &5.429368e-01,5.394330e-01,5.361945e-01,5.331840e-01,5.303713e-01,
     &5.277321e-01,5.252462e-01,5.228967e-01,5.206694e-01,5.185522e-01,
     &5.165348e-01,5.146082e-01,5.127644e-01,5.109968e-01,5.092993e-01,
     &5.076665e-01,5.060936e-01,5.045765e-01,5.031113e-01,5.016946e-01,
     &5.003232e-01,4.989945e-01,4.977057e-01,4.964546e-01,4.952390e-01,
     &4.940570e-01,4.929068e-01,4.917867e-01,4.906952e-01,4.896308e-01,
     &4.885923e-01,4.875785e-01,4.865881e-01/
      DATA (SSAICE2(I,11),I=1,43) /
C    BAND 11
     &6.232263e-01,6.037961e-01,5.906828e-01,5.807779e-01,5.728182e-01,
     &5.661645e-01,5.604483e-01,5.554375e-01,5.509768e-01,5.469570e-01,
     &5.432984e-01,5.399411e-01,5.368390e-01,5.339556e-01,5.312619e-01,
     &5.287341e-01,5.263529e-01,5.241018e-01,5.219671e-01,5.199373e-01,
     &5.180022e-01,5.161534e-01,5.143831e-01,5.126849e-01,5.110531e-01,
     &5.094823e-01,5.079682e-01,5.065066e-01,5.050939e-01,5.037269e-01,
     &5.024025e-01,5.011181e-01,4.998713e-01,4.986598e-01,4.974815e-01,
     &4.963348e-01,4.952177e-01,4.941288e-01,4.930666e-01,4.920297e-01,
     &4.910170e-01,4.900272e-01,4.890593e-01/
      DATA (SSAICE2(I,12),I=1,43) /
C    BAND 12
     &8.165189e-01,7.690707e-01,7.369135e-01,7.125590e-01,6.929511e-01,
     &6.765386e-01,6.624248e-01,6.500445e-01,6.390180e-01,6.290783e-01,
     &6.200302e-01,6.117269e-01,6.040549e-01,5.969249e-01,5.902653e-01,
     &5.840177e-01,5.781341e-01,5.725743e-01,5.673045e-01,5.622957e-01,
     &5.575233e-01,5.529660e-01,5.486050e-01,5.444242e-01,5.404091e-01,
     &5.365471e-01,5.328269e-01,5.292383e-01,5.257724e-01,5.224210e-01,
     &5.191768e-01,5.160329e-01,5.129835e-01,5.100229e-01,5.071461e-01,
     &5.043485e-01,5.016258e-01,4.989740e-01,4.963895e-01,4.938691e-01,
     &4.914094e-01,4.890078e-01,4.866614e-01/
      DATA (SSAICE2(I,13),I=1,43) /
C    BAND 13
     &7.081550e-01,6.764607e-01,6.549873e-01,6.387269e-01,6.256366e-01,
     &6.146798e-01,6.052573e-01,5.969915e-01,5.896290e-01,5.829913e-01,
     &5.769484e-01,5.714020e-01,5.662765e-01,5.615123e-01,5.570617e-01,
     &5.528856e-01,5.489522e-01,5.452345e-01,5.417100e-01,5.383595e-01,
     &5.351665e-01,5.321168e-01,5.291979e-01,5.263991e-01,5.237107e-01,
     &5.211244e-01,5.186325e-01,5.162284e-01,5.139061e-01,5.116601e-01,
     &5.094855e-01,5.073778e-01,5.053332e-01,5.033477e-01,5.014182e-01,
     &4.995414e-01,4.977147e-01,4.959352e-01,4.942007e-01,4.925089e-01,
     &4.908577e-01,4.892452e-01,4.876697e-01/
      DATA (SSAICE2(I,14),I=1,43) /
C    BAND 14
     &7.353033e-01,6.997772e-01,6.756819e-01,6.574216e-01,6.427122e-01,
     &6.303940e-01,6.197965e-01,6.104970e-01,6.022115e-01,5.947404e-01,
     &5.879375e-01,5.816930e-01,5.759219e-01,5.705575e-01,5.655461e-01,
     &5.608441e-01,5.564153e-01,5.522298e-01,5.482621e-01,5.444907e-01,
     &5.408970e-01,5.374650e-01,5.341808e-01,5.310321e-01,5.280082e-01,
     &5.250995e-01,5.222976e-01,5.195949e-01,5.169845e-01,5.144605e-01,
     &5.120172e-01,5.096496e-01,5.073532e-01,5.051238e-01,5.029576e-01,
     &5.008511e-01,4.988011e-01,4.968047e-01,4.948591e-01,4.929617e-01,
     &4.911103e-01,4.893027e-01,4.875368e-01/
      DATA (SSAICE2(I,15),I=1,43) /
C    BAND 15
     &8.248475e-01,7.814674e-01,7.518947e-01,7.294008e-01,7.112284e-01,
     &6.959732e-01,6.828218e-01,6.712600e-01,6.609423e-01,6.516250e-01,
     &6.431297e-01,6.353220e-01,6.280981e-01,6.213760e-01,6.150901e-01,
     &6.091866e-01,6.036215e-01,5.983576e-01,5.933637e-01,5.886134e-01,
     &5.840837e-01,5.797549e-01,5.756098e-01,5.716333e-01,5.678121e-01,
     &5.641345e-01,5.605899e-01,5.571691e-01,5.538635e-01,5.506656e-01,
     &5.475687e-01,5.445663e-01,5.416530e-01,5.388234e-01,5.360730e-01,
     &5.333974e-01,5.307925e-01,5.282548e-01,5.257807e-01,5.233673e-01,
     &5.210114e-01,5.187106e-01,5.164621e-01/
      DATA (SSAICE2(I,16),I=1,43) /
C    BAND 16
     &6.630615e-01,6.451169e-01,6.333696e-01,6.246927e-01,6.178420e-01,
     &6.121976e-01,6.074069e-01,6.032505e-01,5.995830e-01,5.963030e-01,
     &5.933372e-01,5.906311e-01,5.881427e-01,5.858395e-01,5.836955e-01,
     &5.816896e-01,5.798046e-01,5.780264e-01,5.763429e-01,5.747441e-01,
     &5.732213e-01,5.717672e-01,5.703754e-01,5.690403e-01,5.677571e-01,
     &5.665215e-01,5.653297e-01,5.641782e-01,5.630643e-01,5.619850e-01,
     &5.609381e-01,5.599214e-01,5.589328e-01,5.579707e-01,5.570333e-01,
     &5.561193e-01,5.552272e-01,5.543558e-01,5.535041e-01,5.526708e-01,
     &5.518551e-01,5.510561e-01,5.502729e-01/
C     ASYMMETRY FACTOR: Unitless
      DATA (ASYICE2(I,1),I=1,43) /
C    BAND 1
     &2.255639e-01,4.645627e-01,5.756219e-01,6.411863e-01,6.850450e-01,
     &7.167515e-01,7.409211e-01,7.600705e-01,7.756950e-01,7.887423e-01,
     &7.998436e-01,8.094368e-01,8.178357e-01,8.252715e-01,8.319187e-01,
     &8.379115e-01,8.433551e-01,8.483330e-01,8.529127e-01,8.571493e-01,
     &8.610881e-01,8.647670e-01,8.682179e-01,8.714677e-01,8.745396e-01,
     &8.774535e-01,8.802267e-01,8.828741e-01,8.854091e-01,8.878434e-01,
     &8.901874e-01,8.924504e-01,8.946407e-01,8.967659e-01,8.988330e-01,
     &9.008482e-01,9.028173e-01,9.047457e-01,9.066384e-01,9.085002e-01,
     &9.103353e-01,9.121481e-01,9.139426e-01/
      DATA (ASYICE2(I,2),I=1,43) /
C    BAND 2
     &5.393286e-01,6.558766e-01,7.164199e-01,7.545486e-01,7.811948e-01,
     &8.010778e-01,8.165998e-01,8.291243e-01,8.394884e-01,8.482370e-01,
     &8.557418e-01,8.622656e-01,8.680002e-01,8.730891e-01,8.776420e-01,
     &8.817442e-01,8.854634e-01,8.888538e-01,8.919594e-01,8.948166e-01,
     &8.974552e-01,8.999005e-01,9.021735e-01,9.042922e-01,9.062719e-01,
     &9.081256e-01,9.098647e-01,9.114988e-01,9.130362e-01,9.144842e-01,
     &9.158490e-01,9.171357e-01,9.183488e-01,9.194918e-01,9.205677e-01,
     &9.215785e-01,9.225256e-01,9.234093e-01,9.242292e-01,9.249837e-01,
     &9.256698e-01,9.262828e-01,9.268157e-01/
      DATA (ASYICE2(I,3),I=1,43) /
C    BAND 3
     &6.402550e-01,7.366100e-01,7.861283e-01,8.170823e-01,8.385946e-01,
     &8.545773e-01,8.670110e-01,8.770146e-01,8.852724e-01,8.922285e-01,
     &8.981848e-01,9.033542e-01,9.078918e-01,9.119133e-01,9.155069e-01,
     &9.187415e-01,9.216713e-01,9.243398e-01,9.267824e-01,9.290280e-01,
     &9.311009e-01,9.330210e-01,9.348055e-01,9.364687e-01,9.380231e-01,
     &9.394792e-01,9.408463e-01,9.421324e-01,9.433444e-01,9.444886e-01,
     &9.455702e-01,9.465940e-01,9.475642e-01,9.484844e-01,9.493581e-01,
     &9.501880e-01,9.509766e-01,9.517263e-01,9.524388e-01,9.531158e-01,
     &9.537587e-01,9.543686e-01,9.549462e-01/
      DATA (ASYICE2(I,4),I=1,43) /
C    BAND 4
     &6.868425e-01,7.885874e-01,8.342997e-01,8.602518e-01,8.769749e-01,
     &8.886475e-01,8.972569e-01,9.038686e-01,9.091054e-01,9.133553e-01,
     &9.168731e-01,9.198324e-01,9.223563e-01,9.245338e-01,9.264314e-01,
     &9.280994e-01,9.295769e-01,9.308944e-01,9.320762e-01,9.331421e-01,
     &9.341079e-01,9.349870e-01,9.357901e-01,9.365266e-01,9.372040e-01,
     &9.378290e-01,9.384073e-01,9.389435e-01,9.394420e-01,9.399062e-01,
     &9.403395e-01,9.407445e-01,9.411237e-01,9.414794e-01,9.418132e-01,
     &9.421271e-01,9.424225e-01,9.427007e-01,9.429630e-01,9.432104e-01,
     &9.434440e-01,9.436646e-01,9.438730e-01/
      DATA (ASYICE2(I,5),I=1,43) /
C    BAND 5
     &7.273309e-01,8.266992e-01,8.655715e-01,8.855658e-01,8.974914e-01,
     &9.053017e-01,9.107583e-01,9.147555e-01,9.177919e-01,9.201655e-01,
     &9.220647e-01,9.236137e-01,9.248978e-01,9.259770e-01,9.268948e-01,
     &9.276835e-01,9.283676e-01,9.289656e-01,9.294923e-01,9.299591e-01,
     &9.303752e-01,9.307482e-01,9.310841e-01,9.313880e-01,9.316640e-01,
     &9.319156e-01,9.321458e-01,9.323571e-01,9.325516e-01,9.327311e-01,
     &9.328973e-01,9.330515e-01,9.331948e-01,9.333284e-01,9.334532e-01,
     &9.335698e-01,9.336792e-01,9.337819e-01,9.338784e-01,9.339694e-01,
     &9.340551e-01,9.341361e-01,9.342127e-01/
      DATA (ASYICE2(I,6),I=1,43) /
C    BAND 6
     &7.887503e-01,8.770704e-01,9.104538e-01,9.272855e-01,9.371982e-01,
     &9.436352e-01,9.481052e-01,9.513645e-01,9.538305e-01,9.557506e-01,
     &9.572803e-01,9.585216e-01,9.595441e-01,9.603968e-01,9.611150e-01,
     &9.617249e-01,9.622461e-01,9.626938e-01,9.630796e-01,9.634129e-01,
     &9.637009e-01,9.639497e-01,9.641639e-01,9.643477e-01,9.645041e-01,
     &9.646359e-01,9.647453e-01,9.648341e-01,9.649039e-01,9.649559e-01,
     &9.649912e-01,9.650105e-01,9.650147e-01,9.650041e-01,9.649792e-01,
     &9.649403e-01,9.648876e-01,9.648210e-01,9.647406e-01,9.646461e-01,
     &9.645374e-01,9.644141e-01,9.642758e-01/
      DATA (ASYICE2(I,7),I=1,43) /
C    BAND 7
     &8.378017e-01,8.852928e-01,9.090067e-01,9.234678e-01,9.333026e-01,
     &9.404700e-01,9.459500e-01,9.502900e-01,9.538212e-01,9.567564e-01,
     &9.592387e-01,9.613686e-01,9.632181e-01,9.648409e-01,9.662775e-01,
     &9.675592e-01,9.687106e-01,9.697512e-01,9.706969e-01,9.715604e-01,
     &9.723525e-01,9.730820e-01,9.737564e-01,9.743818e-01,9.749637e-01,
     &9.755068e-01,9.760148e-01,9.764914e-01,9.769395e-01,9.773618e-01,
     &9.777606e-01,9.781380e-01,9.784958e-01,9.788357e-01,9.791591e-01,
     &9.794674e-01,9.797619e-01,9.800437e-01,9.803137e-01,9.805730e-01,
     &9.808224e-01,9.810629e-01,9.812952e-01/
      DATA (ASYICE2(I,8),I=1,43) /
C    BAND 8
     &8.410085e-01,8.742948e-01,8.935566e-01,9.065956e-01,9.162142e-01,
     &9.237079e-01,9.297715e-01,9.348164e-01,9.391043e-01,9.428109e-01,
     &9.460592e-01,9.489383e-01,9.515147e-01,9.538391e-01,9.559510e-01,
     &9.578816e-01,9.596561e-01,9.612949e-01,9.628148e-01,9.642301e-01,
     &9.655524e-01,9.667918e-01,9.679568e-01,9.690549e-01,9.700923e-01,
     &9.710747e-01,9.720070e-01,9.728933e-01,9.737376e-01,9.745431e-01,
     &9.753129e-01,9.760497e-01,9.767558e-01,9.774336e-01,9.780849e-01,
     &9.787115e-01,9.793151e-01,9.798972e-01,9.804592e-01,9.810022e-01,
     &9.815276e-01,9.820363e-01,9.825294e-01/
      DATA (ASYICE2(I,9),I=1,43) /
C    BAND 9
     &8.447005e-01,8.780132e-01,8.968896e-01,9.095104e-01,9.187442e-01,
     &9.258958e-01,9.316569e-01,9.364336e-01,9.404822e-01,9.439738e-01,
     &9.470275e-01,9.497295e-01,9.521436e-01,9.543184e-01,9.562917e-01,
     &9.580932e-01,9.597468e-01,9.612720e-01,9.626849e-01,9.639986e-01,
     &9.652243e-01,9.663716e-01,9.674484e-01,9.684618e-01,9.694176e-01,
     &9.703212e-01,9.711770e-01,9.719892e-01,9.727611e-01,9.734961e-01,
     &9.741968e-01,9.748658e-01,9.755053e-01,9.761174e-01,9.767038e-01,
     &9.772662e-01,9.778062e-01,9.783250e-01,9.788239e-01,9.793041e-01,
     &9.797666e-01,9.802124e-01,9.806423e-01/
      DATA (ASYICE2(I,10),I=1,43) /
C    BAND 10
     &8.564947e-01,8.915851e-01,9.102676e-01,9.221993e-01,9.306152e-01,
     &9.369369e-01,9.418971e-01,9.459156e-01,9.492519e-01,9.520761e-01,
     &9.545046e-01,9.566202e-01,9.584835e-01,9.601400e-01,9.616244e-01,
     &9.629641e-01,9.641805e-01,9.652912e-01,9.663102e-01,9.672493e-01,
     &9.681181e-01,9.689247e-01,9.696762e-01,9.703783e-01,9.710361e-01,
     &9.716540e-01,9.722357e-01,9.727846e-01,9.733035e-01,9.737950e-01,
     &9.742615e-01,9.747048e-01,9.751268e-01,9.755291e-01,9.759132e-01,
     &9.762803e-01,9.766317e-01,9.769684e-01,9.772913e-01,9.776014e-01,
     &9.778994e-01,9.781862e-01,9.784623e-01/
      DATA (ASYICE2(I,11),I=1,43) /
C    BAND 11
     &8.697900e-01,9.013810e-01,9.181518e-01,9.288524e-01,9.363985e-01,
     &9.420679e-01,9.465182e-01,9.501255e-01,9.531224e-01,9.556611e-01,
     &9.578459e-01,9.597507e-01,9.614299e-01,9.629240e-01,9.642643e-01,
     &9.654751e-01,9.665758e-01,9.675818e-01,9.685059e-01,9.693585e-01,
     &9.701484e-01,9.708827e-01,9.715676e-01,9.722085e-01,9.728097e-01,
     &9.733753e-01,9.739087e-01,9.744127e-01,9.748899e-01,9.753428e-01,
     &9.757732e-01,9.761830e-01,9.765738e-01,9.769470e-01,9.773040e-01,
     &9.776459e-01,9.779737e-01,9.782883e-01,9.785908e-01,9.788818e-01,
     &9.791621e-01,9.794323e-01,9.796931e-01/
      DATA (ASYICE2(I,12),I=1,43) /
C    BAND 12
     &8.283310e-01,8.543591e-01,8.712236e-01,8.835970e-01,8.933168e-01,
     &9.012909e-01,9.080327e-01,9.138601e-01,9.189835e-01,9.235486e-01,
     &9.276609e-01,9.313988e-01,9.348222e-01,9.379780e-01,9.409034e-01,
     &9.436284e-01,9.461776e-01,9.485715e-01,9.508272e-01,9.529591e-01,
     &9.549796e-01,9.568992e-01,9.587273e-01,9.604716e-01,9.621394e-01,
     &9.637367e-01,9.652691e-01,9.667413e-01,9.681578e-01,9.695225e-01,
     &9.708388e-01,9.721100e-01,9.733388e-01,9.745280e-01,9.756799e-01,
     &9.767967e-01,9.778803e-01,9.789326e-01,9.799553e-01,9.809499e-01,
     &9.819179e-01,9.828606e-01,9.837793e-01/
      DATA (ASYICE2(I,13),I=1,43) /
C    BAND 13
     &8.306222e-01,8.696458e-01,8.911483e-01,9.052375e-01,9.153823e-01,
     &9.231361e-01,9.293121e-01,9.343824e-01,9.386425e-01,9.422880e-01,
     &9.454540e-01,9.482376e-01,9.507103e-01,9.529263e-01,9.549272e-01,
     &9.567459e-01,9.584087e-01,9.599368e-01,9.613475e-01,9.626553e-01,
     &9.638721e-01,9.650082e-01,9.660721e-01,9.670713e-01,9.680121e-01,
     &9.689001e-01,9.697400e-01,9.705361e-01,9.712922e-01,9.720114e-01,
     &9.726968e-01,9.733509e-01,9.739760e-01,9.745744e-01,9.751477e-01,
     &9.756979e-01,9.762263e-01,9.767345e-01,9.772236e-01,9.776949e-01,
     &9.781495e-01,9.785882e-01,9.790121e-01/
      DATA (ASYICE2(I,14),I=1,43) /
C    BAND 14
     &8.240566e-01,8.644102e-01,8.868070e-01,9.015382e-01,9.121690e-01,
     &9.203056e-01,9.267921e-01,9.321201e-01,9.365980e-01,9.404302e-01,
     &9.437583e-01,9.466840e-01,9.492823e-01,9.516101e-01,9.537112e-01,
     &9.556203e-01,9.573649e-01,9.589674e-01,9.604460e-01,9.618159e-01,
     &9.630899e-01,9.642786e-01,9.653911e-01,9.664352e-01,9.674178e-01,
     &9.683445e-01,9.692205e-01,9.700502e-01,9.708376e-01,9.715862e-01,
     &9.722990e-01,9.729789e-01,9.736282e-01,9.742491e-01,9.748438e-01,
     &9.754140e-01,9.759613e-01,9.764872e-01,9.769931e-01,9.774802e-01,
     &9.779496e-01,9.784025e-01,9.788396e-01/
      DATA (ASYICE2(I,15),I=1,43) /
C    BAND 15
     &7.868651e-01,8.322305e-01,8.577581e-01,8.747072e-01,8.870264e-01,
     &8.965099e-01,9.041069e-01,9.103734e-01,9.156595e-01,9.201985e-01,
     &9.241524e-01,9.276378e-01,9.307412e-01,9.335281e-01,9.360494e-01,
     &9.383451e-01,9.404472e-01,9.423817e-01,9.441700e-01,9.458298e-01,
     &9.473759e-01,9.488208e-01,9.501753e-01,9.514485e-01,9.526482e-01,
     &9.537815e-01,9.548541e-01,9.558715e-01,9.568383e-01,9.577585e-01,
     &9.586359e-01,9.594736e-01,9.602747e-01,9.610417e-01,9.617771e-01,
     &9.624829e-01,9.631611e-01,9.638136e-01,9.644418e-01,9.650473e-01,
     &9.656314e-01,9.661955e-01,9.667405e-01/
      DATA (ASYICE2(I,16),I=1,43) /
C    BAND 16
     &7.946655e-01,8.547685e-01,8.806016e-01,8.949880e-01,9.041676e-01,
     &9.105399e-01,9.152249e-01,9.188160e-01,9.216573e-01,9.239620e-01,
     &9.258695e-01,9.274745e-01,9.288441e-01,9.300267e-01,9.310584e-01,
     &9.319665e-01,9.327721e-01,9.334918e-01,9.341387e-01,9.347236e-01,
     &9.352551e-01,9.357402e-01,9.361850e-01,9.365942e-01,9.369722e-01,
     &9.373225e-01,9.376481e-01,9.379516e-01,9.382352e-01,9.385010e-01,
     &9.387505e-01,9.389854e-01,9.392070e-01,9.394163e-01,9.396145e-01,
     &9.398024e-01,9.399809e-01,9.401508e-01,9.403126e-01,9.404670e-01,
     &9.406144e-01,9.407555e-01,9.408906e-01/

C     Hexagonal Ice Particle Parameterization
C     EXTINCTION UNITS (EXT COEF/IWC): [(m^-1)/(g m^-3)]
      DATA (EXTICE3(I,1),I=1,46) /
C    BAND 1
     &3.408303e-03,5.644582e-02,8.802256e-02,9.139862e-02,8.751853e-02,
     &8.196726e-02,7.636347e-02,7.118487e-02,6.654102e-02,6.241830e-02,
     &5.876490e-02,5.552147e-02,5.263196e-02,5.004698e-02,4.772428e-02,
     &4.562809e-02,4.372831e-02,4.199961e-02,4.042058e-02,3.897310e-02,
     &3.764178e-02,3.641344e-02,3.527679e-02,3.422209e-02,3.324092e-02,
     &3.232593e-02,3.147072e-02,3.066968e-02,2.991787e-02,2.921093e-02,
     &2.854498e-02,2.791659e-02,2.732268e-02,2.676052e-02,2.622764e-02,
     &2.572182e-02,2.524106e-02,2.478356e-02,2.434767e-02,2.393190e-02,
     &2.353490e-02,2.315542e-02,2.279235e-02,2.244463e-02,2.211133e-02,
     &2.179156e-02/
      DATA (EXTICE3(I,2),I=1,46) /
C    BAND 2
     &1.672430e-03,9.734149e-02,1.408746e-01,1.412504e-01,1.319985e-01,
     &1.210834e-01,1.106662e-01,1.012940e-01,9.302056e-02,8.575111e-02,
     &7.935612e-02,7.370958e-02,6.870035e-02,6.423409e-02,6.023195e-02,
     &5.662833e-02,5.336865e-02,5.040738e-02,4.770635e-02,4.523343e-02,
     &4.296144e-02,4.086725e-02,3.893106e-02,3.713588e-02,3.546701e-02,
     &3.391173e-02,3.245891e-02,3.109885e-02,2.982301e-02,2.862385e-02,
     &2.749470e-02,2.642966e-02,2.542343e-02,2.447131e-02,2.356906e-02,
     &2.271289e-02,2.189937e-02,2.112540e-02,2.038818e-02,1.968516e-02,
     &1.901402e-02,1.837265e-02,1.775911e-02,1.717164e-02,1.660862e-02,
     &1.606856e-02/
      DATA (EXTICE3(I,3),I=1,46) /
C    BAND 3
     &1.972490e-01,2.471923e-01,2.190426e-01,1.886876e-01,1.635161e-01,
     &1.433107e-01,1.270080e-01,1.136759e-01,1.026128e-01,9.330520e-02,
     &8.537673e-02,7.854795e-02,7.260841e-02,6.739721e-02,6.278951e-02,
     &5.868713e-02,5.501189e-02,5.170080e-02,4.870261e-02,4.597520e-02,
     &4.348362e-02,4.119867e-02,3.909578e-02,3.715408e-02,3.535578e-02,
     &3.368560e-02,3.213036e-02,3.067862e-02,2.932039e-02,2.804694e-02,
     &2.685059e-02,2.572455e-02,2.466282e-02,2.366004e-02,2.271144e-02,
     &2.181276e-02,2.096016e-02,2.015018e-02,1.937972e-02,1.864596e-02,
     &1.794633e-02,1.727851e-02,1.664038e-02,1.603000e-02,1.544561e-02,
     &1.488558e-02/
      DATA (EXTICE3(I,4),I=1,46) /
C    BAND 4
     &3.811222e-01,2.987800e-01,2.363621e-01,1.936824e-01,1.633485e-01,
     &1.408443e-01,1.235383e-01,1.098370e-01,9.873002e-02,8.954892e-02,
     &8.183532e-02,7.526487e-02,6.960186e-02,6.467097e-02,6.033921e-02,
     &5.650387e-02,5.308439e-02,5.001673e-02,4.724933e-02,4.474022e-02,
     &4.245490e-02,4.036473e-02,3.844576e-02,3.667782e-02,3.504375e-02,
     &3.352892e-02,3.212075e-02,3.080835e-02,2.958230e-02,2.843434e-02,
     &2.735724e-02,2.634465e-02,2.539094e-02,2.449112e-02,2.364075e-02,
     &2.283587e-02,2.207291e-02,2.134870e-02,2.066036e-02,2.000527e-02,
     &1.938110e-02,1.878570e-02,1.821714e-02,1.767362e-02,1.715355e-02,
     &1.665542e-02/
      DATA (EXTICE3(I,5),I=1,46) /
C    BAND 5
     &3.995665e-01,2.934192e-01,2.273221e-01,1.845439e-01,1.549229e-01,
     &1.332808e-01,1.168041e-01,1.038520e-01,9.340782e-02,8.481006e-02,
     &7.761020e-02,7.149375e-02,6.623375e-02,6.166239e-02,5.765292e-02,
     &5.410788e-02,5.095109e-02,4.812212e-02,4.557248e-02,4.326279e-02,
     &4.116070e-02,3.923946e-02,3.747670e-02,3.585359e-02,3.435419e-02,
     &3.296488e-02,3.167396e-02,3.047135e-02,2.934829e-02,2.829713e-02,
     &2.731120e-02,2.638460e-02,2.551213e-02,2.468919e-02,2.391168e-02,
     &2.317594e-02,2.247869e-02,2.181699e-02,2.118818e-02,2.058988e-02,
     &2.001992e-02,1.947633e-02,1.895732e-02,1.846126e-02,1.798667e-02,
     &1.753218e-02/
      DATA (EXTICE3(I,6),I=1,46) /
C    BAND 6
     &1.813874e-01,1.978092e-01,1.715767e-01,1.471311e-01,1.276487e-01,
     &1.122980e-01,1.000450e-01,9.009510e-02,8.187953e-02,7.499315e-02,
     &6.914382e-02,6.411723e-02,5.975329e-02,5.593036e-02,5.255454e-02,
     &4.955227e-02,4.686517e-02,4.444635e-02,4.225773e-02,4.026807e-02,
     &3.845151e-02,3.678649e-02,3.525486e-02,3.384124e-02,3.253254e-02,
     &3.131752e-02,3.018649e-02,2.913104e-02,2.814387e-02,2.721856e-02,
     &2.634948e-02,2.553167e-02,2.476072e-02,2.403273e-02,2.334421e-02,
     &2.269203e-02,2.207339e-02,2.148578e-02,2.092691e-02,2.039474e-02,
     &1.988739e-02,1.940317e-02,1.894054e-02,1.849807e-02,1.807450e-02,
     &1.766862e-02/
      DATA (EXTICE3(I,7),I=1,46) /
C    BAND 7
     &4.608531e-02,1.749754e-01,1.750374e-01,1.581885e-01,1.407428e-01,
     &1.254718e-01,1.125654e-01,1.017008e-01,9.250502e-02,8.465627e-02,
     &7.789684e-02,7.202462e-02,6.688150e-02,6.234315e-02,5.831110e-02,
     &5.470656e-02,5.146596e-02,4.853753e-02,4.587872e-02,4.345428e-02,
     &4.123479e-02,3.919551e-02,3.731550e-02,3.557690e-02,3.396443e-02,
     &3.246491e-02,3.106693e-02,2.976055e-02,2.853708e-02,2.738892e-02,
     &2.630933e-02,2.529237e-02,2.433275e-02,2.342578e-02,2.256724e-02,
     &2.175337e-02,2.098078e-02,2.024641e-02,1.954749e-02,1.888153e-02,
     &1.824626e-02,1.763959e-02,1.705966e-02,1.650472e-02,1.597320e-02,
     &1.546365e-02/
      DATA (EXTICE3(I,8),I=1,46) /
C    BAND 8
     &2.757113e-01,2.671436e-01,2.242236e-01,1.887468e-01,1.615472e-01,
     &1.405234e-01,1.239339e-01,1.105653e-01,9.958709e-02,9.042264e-02,
     &8.266308e-02,7.601193e-02,7.024963e-02,6.521047e-02,6.076722e-02,
     &5.682062e-02,5.329219e-02,5.011908e-02,4.725037e-02,4.464441e-02,
     &4.226678e-02,4.008881e-02,3.808641e-02,3.623922e-02,3.452991e-02,
     &3.294361e-02,3.146754e-02,3.009060e-02,2.880315e-02,2.759675e-02,
     &2.646398e-02,2.539832e-02,2.439397e-02,2.344581e-02,2.254924e-02,
     &2.170017e-02,2.089493e-02,2.013021e-02,1.940303e-02,1.871069e-02,
     &1.805076e-02,1.742100e-02,1.681940e-02,1.624410e-02,1.569343e-02,
     &1.516584e-02/
      DATA (EXTICE3(I,9),I=1,46) /
C    BAND 9
     &4.502510e-01,3.193893e-01,2.440352e-01,1.965220e-01,1.640522e-01,
     &1.405164e-01,1.226933e-01,1.087361e-01,9.751395e-02,8.829643e-02,
     &8.059148e-02,7.405560e-02,6.844183e-02,6.356811e-02,5.929729e-02,
     &5.552411e-02,5.216644e-02,4.915928e-02,4.645050e-02,4.399781e-02,
     &4.176655e-02,3.972805e-02,3.785837e-02,3.613737e-02,3.454802e-02,
     &3.307576e-02,3.170811e-02,3.043432e-02,2.924505e-02,2.813215e-02,
     &2.708850e-02,2.610783e-02,2.518461e-02,2.431393e-02,2.349144e-02,
     &2.271324e-02,2.197585e-02,2.127615e-02,2.061131e-02,1.997879e-02,
     &1.937629e-02,1.880173e-02,1.825320e-02,1.772898e-02,1.722749e-02,
     &1.674727e-02/
      DATA (EXTICE3(I,10),I=1,46) /
C    BAND 10
     &5.006417e-01,3.306610e-01,2.458712e-01,1.953080e-01,1.617721e-01,
     &1.379144e-01,1.200780e-01,1.062405e-01,9.519348e-02,8.617061e-02,
     &7.866248e-02,7.231733e-02,6.688443e-02,6.218032e-02,5.806758e-02,
     &5.444133e-02,5.122008e-02,4.833958e-02,4.574849e-02,4.340528e-02,
     &4.127601e-02,3.933266e-02,3.755191e-02,3.591416e-02,3.440286e-02,
     &3.300390e-02,3.170522e-02,3.049640e-02,2.936844e-02,2.831347e-02,
     &2.732464e-02,2.639591e-02,2.552198e-02,2.469812e-02,2.392016e-02,
     &2.318436e-02,2.248739e-02,2.182625e-02,2.119825e-02,2.060095e-02,
     &2.003217e-02,1.948990e-02,1.897234e-02,1.847783e-02,1.800487e-02,
     &1.755207e-02/
      DATA (EXTICE3(I,11),I=1,46) /
C    BAND 11
     &5.218467e-01,3.358426e-01,2.470669e-01,1.951513e-01,1.610989e-01,
     &1.370472e-01,1.191559e-01,1.053276e-01,9.431944e-02,8.534878e-02,
     &7.789788e-02,7.161072e-02,6.623443e-02,6.158443e-02,5.752286e-02,
     &5.394468e-02,5.076843e-02,4.792999e-02,4.537818e-02,4.307168e-02,
     &4.097672e-02,3.906547e-02,3.731479e-02,3.570525e-02,3.422043e-02,
     &3.284639e-02,3.157116e-02,3.038446e-02,2.927736e-02,2.824213e-02,
     &2.727197e-02,2.636094e-02,2.550379e-02,2.469588e-02,2.393307e-02,
     &2.321170e-02,2.252846e-02,2.188042e-02,2.126493e-02,2.067958e-02,
     &2.012221e-02,1.959087e-02,1.908376e-02,1.859927e-02,1.813592e-02,
     &1.769235e-02/
      DATA (EXTICE3(I,12),I=1,46) /
C    BAND 12
     &5.974497e-01,3.649665e-01,2.622591e-01,2.044015e-01,1.672868e-01,
     &1.414571e-01,1.224453e-01,1.078669e-01,9.633309e-02,8.698039e-02,
     &7.924354e-02,7.273708e-02,6.718913e-02,6.240238e-02,5.823024e-02,
     &5.456145e-02,5.131011e-02,4.840878e-02,4.580383e-02,4.345204e-02,
     &4.131820e-02,3.937336e-02,3.759346e-02,3.595836e-02,3.445109e-02,
     &3.305720e-02,3.176437e-02,3.056198e-02,2.944085e-02,2.839302e-02,
     &2.741152e-02,2.649024e-02,2.562379e-02,2.480743e-02,2.403692e-02,
     &2.330851e-02,2.261883e-02,2.196487e-02,2.134392e-02,2.075355e-02,
     &2.019154e-02,1.965589e-02,1.914480e-02,1.865660e-02,1.818979e-02,
     &1.774300e-02/
      DATA (EXTICE3(I,13),I=1,46) /
C    BAND 13
     &6.043932e-01,3.625783e-01,2.584703e-01,2.006099e-01,1.638026e-01,
     &1.383330e-01,1.196641e-01,1.053939e-01,9.413223e-02,8.501864e-02,
     &7.749223e-02,7.117168e-02,6.578876e-02,6.114930e-02,5.710925e-02,
     &5.355953e-02,5.041600e-02,4.761274e-02,4.509735e-02,4.282767e-02,
     &4.076939e-02,3.889430e-02,3.717899e-02,3.560388e-02,3.415247e-02,
     &3.281075e-02,3.156674e-02,3.041015e-02,2.933208e-02,2.832481e-02,
     &2.738159e-02,2.649650e-02,2.566434e-02,2.488050e-02,2.414089e-02,
     &2.344188e-02,2.278022e-02,2.215300e-02,2.155760e-02,2.099166e-02,
     &2.045305e-02,1.993985e-02,1.945029e-02,1.898280e-02,1.853590e-02,
     &1.810828e-02/
      DATA (EXTICE3(I,14),I=1,46) /
C    BAND 14
     &6.072115e-01,3.603507e-01,2.556553e-01,1.979378e-01,1.614082e-01,
     &1.362193e-01,1.178036e-01,1.037544e-01,9.268423e-02,8.373670e-02,
     &7.635495e-02,7.016114e-02,6.488994e-02,6.034953e-02,5.639781e-02,
     &5.292726e-02,4.985507e-02,4.711636e-02,4.465963e-02,4.244346e-02,
     &4.043417e-02,3.860408e-02,3.693023e-02,3.539343e-02,3.397752e-02,
     &3.266875e-02,3.145541e-02,3.032744e-02,2.927611e-02,2.829387e-02,
     &2.737414e-02,2.651111e-02,2.569970e-02,2.493542e-02,2.421427e-02,
     &2.353269e-02,2.288752e-02,2.227591e-02,2.169531e-02,2.114341e-02,
     &2.061814e-02,2.011762e-02,1.964012e-02,1.918411e-02,1.874815e-02,
     &1.833095e-02/
      DATA (EXTICE3(I,15),I=1,46) /
C    BAND 15
     &6.177203e-01,3.627414e-01,2.560843e-01,1.977307e-01,1.609776e-01,
     &1.357202e-01,1.173003e-01,1.032749e-01,9.224026e-02,8.333247e-02,
     &7.599102e-02,6.983637e-02,6.460236e-02,6.009687e-02,5.617773e-02,
     &5.273750e-02,4.969348e-02,4.698095e-02,4.454857e-02,4.235508e-02,
     &4.036693e-02,3.855657e-02,3.690119e-02,3.538169e-02,3.398201e-02,
     &3.268851e-02,3.148955e-02,3.037512e-02,2.933659e-02,2.836647e-02,
     &2.745820e-02,2.660605e-02,2.580498e-02,2.505053e-02,2.433873e-02,
     &2.366608e-02,2.302942e-02,2.242594e-02,2.185312e-02,2.130867e-02,
     &2.079054e-02,2.029686e-02,1.982594e-02,1.937624e-02,1.894635e-02,
     &1.853500e-02/
      DATA (EXTICE3(I,16),I=1,46) /
C    BAND 16
     &5.194013e-01,3.215089e-01,2.327917e-01,1.824424e-01,1.499977e-01,
     &1.273492e-01,1.106421e-01,9.780982e-02,8.764435e-02,7.939266e-02,
     &7.256081e-02,6.681137e-02,6.190600e-02,5.767154e-02,5.397915e-02,
     &5.073102e-02,4.785151e-02,4.528125e-02,4.297296e-02,4.088853e-02,
     &3.899690e-02,3.727251e-02,3.569411e-02,3.424393e-02,3.290694e-02,
     &3.167040e-02,3.052340e-02,2.945654e-02,2.846172e-02,2.753188e-02,
     &2.666085e-02,2.584322e-02,2.507423e-02,2.434967e-02,2.366579e-02,
     &2.301926e-02,2.240711e-02,2.182666e-02,2.127551e-02,2.075150e-02,
     &2.025267e-02,1.977725e-02,1.932364e-02,1.889035e-02,1.847607e-02,
     &1.807956e-02/
C     SINGLE-SCATTERING ALBEDO: Unitless
      DATA (SSAICE3(I,1),I=1,46) /
C    BAND 1
     &8.733202e-02,1.733042e-01,2.494597e-01,2.853637e-01,3.129915e-01,
     &3.366261e-01,3.574925e-01,3.761110e-01,3.927678e-01,4.076523e-01,
     &4.209083e-01,4.326557e-01,4.430008e-01,4.520419e-01,4.598721e-01,
     &4.665813e-01,4.722567e-01,4.769845e-01,4.808494e-01,4.839354e-01,
     &4.863257e-01,4.881034e-01,4.893508e-01,4.901503e-01,4.905838e-01,
     &4.907333e-01,4.906804e-01,4.905066e-01,4.902935e-01,4.901225e-01,
     &4.900749e-01,4.902321e-01,4.906751e-01,4.914853e-01,4.927439e-01,
     &4.945318e-01,4.969304e-01,5.000205e-01,5.038834e-01,5.086000e-01,
     &5.142513e-01,5.209185e-01,5.286824e-01,5.376241e-01,5.478246e-01,
     &5.593649e-01/
      DATA (SSAICE3(I,2),I=1,46) /
C    BAND 2
     &7.617260e-01,8.272990e-01,8.134738e-01,8.040737e-01,7.953976e-01,
     &7.869914e-01,7.787330e-01,7.705796e-01,7.625153e-01,7.545359e-01,
     &7.466425e-01,7.388394e-01,7.311327e-01,7.235295e-01,7.160381e-01,
     &7.086672e-01,7.014261e-01,6.943246e-01,6.873731e-01,6.805821e-01,
     &6.739626e-01,6.675263e-01,6.612848e-01,6.552505e-01,6.494361e-01,
     &6.438547e-01,6.385199e-01,6.334458e-01,6.286469e-01,6.241384e-01,
     &6.199357e-01,6.160553e-01,6.125138e-01,6.093286e-01,6.065180e-01,
     &6.041006e-01,6.020961e-01,6.005248e-01,5.994078e-01,5.987672e-01,
     &5.986259e-01,5.990079e-01,5.999382e-01,6.014428e-01,6.035491e-01,
     &6.062855e-01/
      DATA (SSAICE3(I,3),I=1,46) /
C    BAND 3
     &6.485070e-01,6.545008e-01,6.483874e-01,6.411445e-01,6.338616e-01,
     &6.268166e-01,6.201055e-01,6.137647e-01,6.078060e-01,6.022297e-01,
     &5.970301e-01,5.921982e-01,5.877229e-01,5.835920e-01,5.797922e-01,
     &5.763099e-01,5.731308e-01,5.702407e-01,5.676249e-01,5.652685e-01,
     &5.631568e-01,5.612747e-01,5.596072e-01,5.581391e-01,5.568553e-01,
     &5.557405e-01,5.547796e-01,5.539571e-01,5.532579e-01,5.526665e-01,
     &5.521676e-01,5.517458e-01,5.513856e-01,5.510717e-01,5.507886e-01,
     &5.505208e-01,5.502528e-01,5.499691e-01,5.496542e-01,5.492923e-01,
     &5.488681e-01,5.483658e-01,5.477697e-01,5.470643e-01,5.462337e-01,
     &5.452623e-01/
      DATA (SSAICE3(I,4),I=1,46) /
C    BAND 4
     &5.367012e-01,5.372179e-01,5.366731e-01,5.359685e-01,5.352768e-01,
     &5.346515e-01,5.341121e-01,5.336653e-01,5.333122e-01,5.330513e-01,
     &5.328793e-01,5.327924e-01,5.327859e-01,5.328549e-01,5.329943e-01,
     &5.331989e-01,5.334632e-01,5.337819e-01,5.341494e-01,5.345602e-01,
     &5.350088e-01,5.354897e-01,5.359972e-01,5.365259e-01,5.370701e-01,
     &5.376244e-01,5.381831e-01,5.387408e-01,5.392920e-01,5.398311e-01,
     &5.403526e-01,5.408512e-01,5.413213e-01,5.417575e-01,5.421545e-01,
     &5.425067e-01,5.428090e-01,5.430558e-01,5.432420e-01,5.433623e-01,
     &5.434113e-01,5.433838e-01,5.432748e-01,5.430789e-01,5.427911e-01,
     &5.424063e-01/
      DATA (SSAICE3(I,5),I=1,46) /
C    BAND 5
     &4.144473e-01,4.233085e-01,4.317951e-01,4.397881e-01,4.472768e-01,
     &4.542704e-01,4.607843e-01,4.668359e-01,4.724439e-01,4.776272e-01,
     &4.824055e-01,4.867983e-01,4.908252e-01,4.945063e-01,4.978612e-01,
     &5.009100e-01,5.036726e-01,5.061689e-01,5.084192e-01,5.104433e-01,
     &5.122613e-01,5.138934e-01,5.153598e-01,5.166804e-01,5.178756e-01,
     &5.189655e-01,5.199702e-01,5.209102e-01,5.218055e-01,5.226765e-01,
     &5.235435e-01,5.244268e-01,5.253468e-01,5.263238e-01,5.273783e-01,
     &5.285306e-01,5.298012e-01,5.312106e-01,5.327792e-01,5.345277e-01,
     &5.364765e-01,5.386462e-01,5.410574e-01,5.437308e-01,5.466871e-01,
     &5.499469e-01/
      DATA (SSAICE3(I,6),I=1,46) /
C    BAND 6
     &3.685509e-01,4.062125e-01,4.219575e-01,4.336348e-01,4.434522e-01,
     &4.520497e-01,4.596957e-01,4.665324e-01,4.726498e-01,4.781130e-01,
     &4.829739e-01,4.872771e-01,4.910627e-01,4.943681e-01,4.972286e-01,
     &4.996785e-01,5.017514e-01,5.034799e-01,5.048966e-01,5.060334e-01,
     &5.069224e-01,5.075951e-01,5.080833e-01,5.084185e-01,5.086320e-01,
     &5.087554e-01,5.088200e-01,5.088572e-01,5.088985e-01,5.089752e-01,
     &5.091187e-01,5.093605e-01,5.097319e-01,5.102646e-01,5.109899e-01,
     &5.119394e-01,5.131447e-01,5.146373e-01,5.164489e-01,5.186112e-01,
     &5.211558e-01,5.241145e-01,5.275191e-01,5.314014e-01,5.357933e-01,
     &5.407268e-01/
      DATA (SSAICE3(I,7),I=1,46) /
C    BAND 7
     &7.347648e-01,6.945725e-01,6.844409e-01,6.757818e-01,6.676970e-01,
     &6.599913e-01,6.525964e-01,6.454817e-01,6.386312e-01,6.320353e-01,
     &6.256876e-01,6.195836e-01,6.137200e-01,6.080940e-01,6.027033e-01,
     &5.975460e-01,5.926203e-01,5.879246e-01,5.834574e-01,5.792174e-01,
     &5.752032e-01,5.714137e-01,5.678476e-01,5.645040e-01,5.613817e-01,
     &5.584797e-01,5.557971e-01,5.533330e-01,5.510865e-01,5.490569e-01,
     &5.472434e-01,5.456452e-01,5.442618e-01,5.430926e-01,5.421369e-01,
     &5.413943e-01,5.408645e-01,5.405470e-01,5.404415e-01,5.405479e-01,
     &5.408659e-01,5.413956e-01,5.421369e-01,5.430900e-01,5.442550e-01,
     &5.456321e-01/
      DATA (SSAICE3(I,8),I=1,46) /
C    BAND 8
     &7.565911e-01,7.410307e-01,7.267244e-01,7.132696e-01,7.005889e-01,
     &6.886431e-01,6.774024e-01,6.668406e-01,6.569327e-01,6.476543e-01,
     &6.389815e-01,6.308907e-01,6.233581e-01,6.163602e-01,6.098736e-01,
     &6.038748e-01,5.983403e-01,5.932468e-01,5.885708e-01,5.842888e-01,
     &5.803776e-01,5.768135e-01,5.735734e-01,5.706336e-01,5.679708e-01,
     &5.655616e-01,5.633825e-01,5.614101e-01,5.596210e-01,5.579915e-01,
     &5.564984e-01,5.551181e-01,5.538271e-01,5.526020e-01,5.514191e-01,
     &5.502550e-01,5.490862e-01,5.478890e-01,5.466400e-01,5.453154e-01,
     &5.438918e-01,5.423455e-01,5.406528e-01,5.387902e-01,5.367338e-01,
     &5.344600e-01/
      DATA (SSAICE3(I,9),I=1,46) /
C    BAND 9
     &7.253131e-01,7.112494e-01,6.974209e-01,6.843066e-01,6.719929e-01,
     &6.604896e-01,6.497830e-01,6.398502e-01,6.306641e-01,6.221956e-01,
     &6.144141e-01,6.072887e-01,6.007877e-01,5.948792e-01,5.895312e-01,
     &5.847116e-01,5.803881e-01,5.765282e-01,5.730996e-01,5.700698e-01,
     &5.674064e-01,5.650768e-01,5.630484e-01,5.612888e-01,5.597653e-01,
     &5.584452e-01,5.572961e-01,5.562853e-01,5.553800e-01,5.545478e-01,
     &5.537558e-01,5.529714e-01,5.521620e-01,5.512948e-01,5.503371e-01,
     &5.492561e-01,5.480192e-01,5.465935e-01,5.449463e-01,5.430449e-01,
     &5.408564e-01,5.383480e-01,5.354869e-01,5.322403e-01,5.285753e-01,
     &5.244591e-01/
      DATA (SSAICE3(I,10),I=1,46) /
C    BAND 10
     &6.692312e-01,6.569887e-01,6.455728e-01,6.349744e-01,6.251679e-01,
     &6.161241e-01,6.078124e-01,6.002017e-01,5.932608e-01,5.869582e-01,
     &5.812625e-01,5.761422e-01,5.715658e-01,5.675017e-01,5.639183e-01,
     &5.607841e-01,5.580676e-01,5.557370e-01,5.537609e-01,5.521078e-01,
     &5.507459e-01,5.496437e-01,5.487696e-01,5.480921e-01,5.475796e-01,
     &5.472004e-01,5.469230e-01,5.467159e-01,5.465474e-01,5.463859e-01,
     &5.461999e-01,5.459577e-01,5.456279e-01,5.451788e-01,5.445788e-01,
     &5.437964e-01,5.427999e-01,5.415578e-01,5.400386e-01,5.382105e-01,
     &5.360422e-01,5.335019e-01,5.305581e-01,5.271792e-01,5.233336e-01,
     &5.189898e-01/
      DATA (SSAICE3(I,11),I=1,46) /
C    BAND 11
     &6.597440e-01,6.486312e-01,6.385754e-01,6.293389e-01,6.208334e-01,
     &6.130062e-01,6.058168e-01,5.992306e-01,5.932156e-01,5.877415e-01,
     &5.827791e-01,5.782999e-01,5.742760e-01,5.706800e-01,5.674845e-01,
     &5.646626e-01,5.621876e-01,5.600328e-01,5.581718e-01,5.565782e-01,
     &5.552258e-01,5.540883e-01,5.531398e-01,5.523542e-01,5.517054e-01,
     &5.511675e-01,5.507147e-01,5.503210e-01,5.499606e-01,5.496077e-01,
     &5.492364e-01,5.488211e-01,5.483359e-01,5.477550e-01,5.470529e-01,
     &5.462036e-01,5.451816e-01,5.439610e-01,5.425163e-01,5.408217e-01,
     &5.388515e-01,5.365800e-01,5.339817e-01,5.310307e-01,5.277016e-01,
     &5.239685e-01/
      DATA (SSAICE3(I,12),I=1,46) /
C    BAND 12
     &8.415691e-01,8.237634e-01,8.070239e-01,7.912304e-01,7.763216e-01,
     &7.622533e-01,7.489874e-01,7.364895e-01,7.247271e-01,7.136691e-01,
     &7.032852e-01,6.935461e-01,6.844229e-01,6.758871e-01,6.679108e-01,
     &6.604662e-01,6.535258e-01,6.470623e-01,6.410487e-01,6.354581e-01,
     &6.302637e-01,6.254387e-01,6.209567e-01,6.167911e-01,6.129154e-01,
     &6.093033e-01,6.059286e-01,6.027648e-01,5.997858e-01,5.969653e-01,
     &5.942772e-01,5.916953e-01,5.891935e-01,5.867455e-01,5.843253e-01,
     &5.819067e-01,5.794637e-01,5.769700e-01,5.743997e-01,5.717265e-01,
     &5.689243e-01,5.659670e-01,5.628284e-01,5.594825e-01,5.559030e-01,
     &5.520638e-01/
      DATA (SSAICE3(I,13),I=1,46) /
C    BAND 13
     &7.418384e-01,7.252704e-01,7.097171e-01,6.951785e-01,6.816320e-01,
     &6.690487e-01,6.573967e-01,6.466425e-01,6.367511e-01,6.276869e-01,
     &6.194133e-01,6.118931e-01,6.050886e-01,5.989616e-01,5.934736e-01,
     &5.885856e-01,5.842585e-01,5.804528e-01,5.771289e-01,5.742470e-01,
     &5.717672e-01,5.696493e-01,5.678531e-01,5.663385e-01,5.650650e-01,
     &5.639922e-01,5.630797e-01,5.622870e-01,5.615735e-01,5.608986e-01,
     &5.602218e-01,5.595025e-01,5.586999e-01,5.577736e-01,5.566826e-01,
     &5.553865e-01,5.538444e-01,5.520156e-01,5.498592e-01,5.473345e-01,
     &5.444005e-01,5.410162e-01,5.371407e-01,5.327328e-01,5.277513e-01,
     &5.221549e-01/
      DATA (SSAICE3(I,14),I=1,46) /
C    BAND 14
     &7.431625e-01,7.253592e-01,7.089350e-01,6.937690e-01,6.797727e-01,
     &6.668739e-01,6.550099e-01,6.441244e-01,6.341649e-01,6.250822e-01,
     &6.168287e-01,6.093588e-01,6.026277e-01,5.965913e-01,5.912065e-01,
     &5.864305e-01,5.822207e-01,5.785351e-01,5.753318e-01,5.725690e-01,
     &5.702053e-01,5.681992e-01,5.665094e-01,5.650948e-01,5.639141e-01,
     &5.629263e-01,5.620905e-01,5.613656e-01,5.607109e-01,5.600853e-01,
     &5.594482e-01,5.587586e-01,5.579758e-01,5.570591e-01,5.559676e-01,
     &5.546606e-01,5.530972e-01,5.512367e-01,5.490383e-01,5.464611e-01,
     &5.434641e-01,5.400063e-01,5.360468e-01,5.315444e-01,5.264578e-01,
     &5.207459e-01/
      DATA (SSAICE3(I,15),I=1,46) /
C    BAND 15
     &8.214523e-01,8.021560e-01,7.840431e-01,7.670541e-01,7.511382e-01,
     &7.362494e-01,7.223443e-01,7.093816e-01,6.973208e-01,6.861227e-01,
     &6.757482e-01,6.661588e-01,6.573164e-01,6.491829e-01,6.417204e-01,
     &6.348911e-01,6.286574e-01,6.229816e-01,6.178260e-01,6.131531e-01,
     &6.089254e-01,6.051053e-01,6.016552e-01,5.985377e-01,5.957152e-01,
     &5.931503e-01,5.908054e-01,5.886430e-01,5.866256e-01,5.847156e-01,
     &5.828757e-01,5.810683e-01,5.792558e-01,5.774008e-01,5.754657e-01,
     &5.734130e-01,5.712052e-01,5.688048e-01,5.661741e-01,5.632757e-01,
     &5.600721e-01,5.565255e-01,5.525985e-01,5.482534e-01,5.434526e-01,
     &5.381586e-01/
      DATA (SSAICE3(I,16),I=1,46) /
C    BAND 16
     &6.749442e-01,6.649947e-01,6.565828e-01,6.489928e-01,6.420046e-01,
     &6.355231e-01,6.294964e-01,6.238901e-01,6.186783e-01,6.138395e-01,
     &6.093543e-01,6.052049e-01,6.013742e-01,5.978457e-01,5.946030e-01,
     &5.916302e-01,5.889115e-01,5.864310e-01,5.841731e-01,5.821221e-01,
     &5.802624e-01,5.785785e-01,5.770549e-01,5.756759e-01,5.744262e-01,
     &5.732901e-01,5.722524e-01,5.712974e-01,5.704097e-01,5.695739e-01,
     &5.687747e-01,5.679964e-01,5.672238e-01,5.664415e-01,5.656340e-01,
     &5.647860e-01,5.638821e-01,5.629070e-01,5.618452e-01,5.606815e-01,
     &5.594006e-01,5.579870e-01,5.564255e-01,5.547008e-01,5.527976e-01,
     &5.507005e-01/
C     ASYMMETRY FACTOR: Unitless
      DATA (ASYICE3(I,1),I=1,46) /
C    BAND 1
     &6.596879e-01,6.405129e-01,6.397929e-01,6.577098e-01,6.788112e-01,
     &6.996961e-01,7.195365e-01,7.380555e-01,7.551580e-01,7.708248e-01,
     &7.850743e-01,7.979456e-01,8.094902e-01,8.197678e-01,8.288438e-01,
     &8.367874e-01,8.436711e-01,8.495695e-01,8.545592e-01,8.587185e-01,
     &8.621270e-01,8.648654e-01,8.670153e-01,8.686596e-01,8.698817e-01,
     &8.707659e-01,8.713971e-01,8.718612e-01,8.722444e-01,8.726337e-01,
     &8.731167e-01,8.737814e-01,8.747166e-01,8.760115e-01,8.777560e-01,
     &8.800403e-01,8.829553e-01,8.865923e-01,8.910430e-01,8.963998e-01,
     &9.027552e-01,9.102023e-01,9.188345e-01,9.287454e-01,9.400292e-01,
     &9.523461e-01/
      DATA (ASYICE3(I,2),I=1,46) /
C    BAND 2
     &7.389769e-01,7.411960e-01,7.490019e-01,7.564984e-01,7.637359e-01,
     &7.707210e-01,7.774570e-01,7.839465e-01,7.901930e-01,7.962001e-01,
     &8.019724e-01,8.075150e-01,8.128335e-01,8.179342e-01,8.228237e-01,
     &8.275093e-01,8.319986e-01,8.362997e-01,8.404208e-01,8.443707e-01,
     &8.481583e-01,8.517929e-01,8.552840e-01,8.586412e-01,8.618744e-01,
     &8.649936e-01,8.680092e-01,8.709315e-01,8.737713e-01,8.765394e-01,
     &8.792470e-01,8.819057e-01,8.845275e-01,8.871247e-01,8.897104e-01,
     &8.922982e-01,8.949025e-01,8.975389e-01,9.002236e-01,9.029745e-01,
     &9.058105e-01,9.087524e-01,9.118223e-01,9.150447e-01,9.184456e-01,
     &9.220537e-01/
      DATA (ASYICE3(I,3),I=1,46) /
C    BAND 3
     &7.533538e-01,7.631828e-01,7.742293e-01,7.848919e-01,7.950315e-01,
     &8.046253e-01,8.136762e-01,8.221958e-01,8.301991e-01,8.377028e-01,
     &8.447247e-01,8.512829e-01,8.573959e-01,8.630826e-01,8.683620e-01,
     &8.732532e-01,8.777756e-01,8.819487e-01,8.857920e-01,8.893253e-01,
     &8.925684e-01,8.955411e-01,8.982635e-01,9.007556e-01,9.030377e-01,
     &9.051300e-01,9.070528e-01,9.088266e-01,9.104718e-01,9.120090e-01,
     &9.134588e-01,9.148418e-01,9.161788e-01,9.174904e-01,9.187976e-01,
     &9.201211e-01,9.214818e-01,9.229007e-01,9.243986e-01,9.259967e-01,
     &9.277160e-01,9.295774e-01,9.316023e-01,9.338116e-01,9.362266e-01,
     &9.388685e-01/
      DATA (ASYICE3(I,4),I=1,46) /
C    BAND 4
     &7.840105e-01,7.959983e-01,8.074838e-01,8.182465e-01,8.282675e-01,
     &8.375608e-01,8.461497e-01,8.540609e-01,8.613228e-01,8.679643e-01,
     &8.740148e-01,8.795041e-01,8.844620e-01,8.889182e-01,8.929028e-01,
     &8.964456e-01,8.995766e-01,9.023257e-01,9.047230e-01,9.067984e-01,
     &9.085817e-01,9.101031e-01,9.113922e-01,9.124791e-01,9.133937e-01,
     &9.141659e-01,9.148254e-01,9.154023e-01,9.159263e-01,9.164274e-01,
     &9.169353e-01,9.174799e-01,9.180910e-01,9.187986e-01,9.196324e-01,
     &9.206223e-01,9.217981e-01,9.231897e-01,9.248270e-01,9.267399e-01,
     &9.289582e-01,9.315119e-01,9.344308e-01,9.377451e-01,9.414845e-01,
     &9.456792e-01/
      DATA (ASYICE3(I,5),I=1,46) /
C    BAND 5
     &8.304283e-01,8.399040e-01,8.485749e-01,8.565538e-01,8.638881e-01,
     &8.706126e-01,8.767578e-01,8.823528e-01,8.874253e-01,8.920026e-01,
     &8.961113e-01,8.997780e-01,9.030287e-01,9.058894e-01,9.083861e-01,
     &9.105442e-01,9.123895e-01,9.139473e-01,9.152430e-01,9.163020e-01,
     &9.171497e-01,9.178111e-01,9.183116e-01,9.186762e-01,9.189303e-01,
     &9.190989e-01,9.192071e-01,9.192802e-01,9.193432e-01,9.194211e-01,
     &9.195392e-01,9.197225e-01,9.199961e-01,9.203850e-01,9.209143e-01,
     &9.216090e-01,9.224942e-01,9.235947e-01,9.249355e-01,9.265416e-01,
     &9.284376e-01,9.306484e-01,9.331987e-01,9.361131e-01,9.394162e-01,
     &9.431324e-01/
      DATA (ASYICE3(I,6),I=1,46) /
C    BAND 6
     &8.792712e-01,8.938493e-01,9.015383e-01,9.078755e-01,9.134852e-01,
     &9.185471e-01,9.231387e-01,9.273046e-01,9.310755e-01,9.344765e-01,
     &9.375293e-01,9.402543e-01,9.426709e-01,9.447982e-01,9.466549e-01,
     &9.482596e-01,9.496311e-01,9.507879e-01,9.517486e-01,9.525319e-01,
     &9.531565e-01,9.536411e-01,9.540046e-01,9.542656e-01,9.544431e-01,
     &9.545559e-01,9.546229e-01,9.546631e-01,9.546954e-01,9.547387e-01,
     &9.548121e-01,9.549345e-01,9.551251e-01,9.554028e-01,9.557866e-01,
     &9.562958e-01,9.569494e-01,9.577665e-01,9.587663e-01,9.599680e-01,
     &9.613908e-01,9.630539e-01,9.649767e-01,9.671785e-01,9.696786e-01,
     &9.724966e-01/
      DATA (ASYICE3(I,7),I=1,46) /
C    BAND 7
     &8.992856e-01,9.093587e-01,9.140642e-01,9.183077e-01,9.222541e-01,
     &9.259455e-01,9.294018e-01,9.326363e-01,9.356600e-01,9.384824e-01,
     &9.411129e-01,9.435602e-01,9.458332e-01,9.479406e-01,9.498907e-01,
     &9.516924e-01,9.533540e-01,9.548842e-01,9.562914e-01,9.575842e-01,
     &9.587710e-01,9.598605e-01,9.608610e-01,9.617811e-01,9.626294e-01,
     &9.634142e-01,9.641441e-01,9.648277e-01,9.654732e-01,9.660894e-01,
     &9.666845e-01,9.672673e-01,9.678460e-01,9.684293e-01,9.690257e-01,
     &9.696436e-01,9.702917e-01,9.709786e-01,9.717127e-01,9.725029e-01,
     &9.733577e-01,9.742858e-01,9.752962e-01,9.763975e-01,9.775987e-01,
     &9.789088e-01/
      DATA (ASYICE3(I,8),I=1,46) /
C    BAND 8
     &8.768357e-01,8.832796e-01,8.887092e-01,8.937467e-01,8.984873e-01,
     &9.029645e-01,9.071961e-01,9.111946e-01,9.149701e-01,9.185316e-01,
     &9.218877e-01,9.250464e-01,9.280157e-01,9.308033e-01,9.334170e-01,
     &9.358643e-01,9.381529e-01,9.402904e-01,9.422842e-01,9.441418e-01,
     &9.458708e-01,9.474786e-01,9.489727e-01,9.503604e-01,9.516494e-01,
     &9.528469e-01,9.539604e-01,9.549973e-01,9.559651e-01,9.568711e-01,
     &9.577229e-01,9.585276e-01,9.592929e-01,9.600261e-01,9.607346e-01,
     &9.614258e-01,9.621071e-01,9.627859e-01,9.634698e-01,9.641660e-01,
     &9.648820e-01,9.656253e-01,9.664032e-01,9.672233e-01,9.680930e-01,
     &9.690198e-01/
      DATA (ASYICE3(I,9),I=1,46) /
C    BAND 9
     &8.609711e-01,8.676429e-01,8.740194e-01,8.800768e-01,8.858172e-01,
     &8.912476e-01,8.963767e-01,9.012137e-01,9.057682e-01,9.100498e-01,
     &9.140683e-01,9.178335e-01,9.213552e-01,9.246434e-01,9.277079e-01,
     &9.305587e-01,9.332058e-01,9.356590e-01,9.379284e-01,9.400240e-01,
     &9.419558e-01,9.437340e-01,9.453685e-01,9.468694e-01,9.482469e-01,
     &9.495111e-01,9.506721e-01,9.517401e-01,9.527252e-01,9.536376e-01,
     &9.544876e-01,9.552853e-01,9.560409e-01,9.567647e-01,9.574670e-01,
     &9.581578e-01,9.588476e-01,9.595466e-01,9.602650e-01,9.610130e-01,
     &9.618010e-01,9.626392e-01,9.635379e-01,9.645074e-01,9.655578e-01,
     &9.666994e-01/
      DATA (ASYICE3(I,10),I=1,46) /
C    BAND 10
     &8.805232e-01,8.872100e-01,8.934945e-01,8.993817e-01,9.048833e-01,
     &9.100128e-01,9.147839e-01,9.192106e-01,9.233070e-01,9.270873e-01,
     &9.305656e-01,9.337561e-01,9.366732e-01,9.393308e-01,9.417434e-01,
     &9.439252e-01,9.458903e-01,9.476530e-01,9.492276e-01,9.506283e-01,
     &9.518694e-01,9.529650e-01,9.539295e-01,9.547771e-01,9.555221e-01,
     &9.561786e-01,9.567610e-01,9.572836e-01,9.577605e-01,9.582059e-01,
     &9.586343e-01,9.590598e-01,9.594966e-01,9.599591e-01,9.604615e-01,
     &9.610179e-01,9.616428e-01,9.623503e-01,9.631546e-01,9.640701e-01,
     &9.651110e-01,9.662916e-01,9.676260e-01,9.691286e-01,9.708136e-01,
     &9.726952e-01/
      DATA (ASYICE3(I,11),I=1,46) /
C    BAND 11
     &8.860103e-01,8.920722e-01,8.976788e-01,9.029125e-01,9.078046e-01,
     &9.123744e-01,9.166373e-01,9.206068e-01,9.242956e-01,9.277160e-01,
     &9.308801e-01,9.337998e-01,9.364867e-01,9.389526e-01,9.412092e-01,
     &9.432681e-01,9.451410e-01,9.468395e-01,9.483754e-01,9.497603e-01,
     &9.510059e-01,9.521239e-01,9.531262e-01,9.540244e-01,9.548305e-01,
     &9.555563e-01,9.562136e-01,9.568142e-01,9.573703e-01,9.578935e-01,
     &9.583959e-01,9.588895e-01,9.593861e-01,9.598976e-01,9.604362e-01,
     &9.610135e-01,9.616417e-01,9.623325e-01,9.630979e-01,9.639496e-01,
     &9.648996e-01,9.659594e-01,9.671409e-01,9.684557e-01,9.699153e-01,
     &9.715312e-01/
      DATA (ASYICE3(I,12),I=1,46) /
C    BAND 12
     &8.097123e-01,8.172929e-01,8.245083e-01,8.314251e-01,8.380674e-01,
     &8.444481e-01,8.505763e-01,8.564596e-01,8.621044e-01,8.675169e-01,
     &8.727030e-01,8.776684e-01,8.824186e-01,8.869590e-01,8.912951e-01,
     &8.954321e-01,8.993753e-01,9.031299e-01,9.067009e-01,9.100936e-01,
     &9.133130e-01,9.163640e-01,9.192519e-01,9.219815e-01,9.245578e-01,
     &9.269860e-01,9.292709e-01,9.314175e-01,9.334309e-01,9.353161e-01,
     &9.370781e-01,9.387218e-01,9.402523e-01,9.416746e-01,9.429938e-01,
     &9.442149e-01,9.453428e-01,9.463827e-01,9.473394e-01,9.482180e-01,
     &9.490235e-01,9.497606e-01,9.504344e-01,9.510495e-01,9.516106e-01,
     &9.521225e-01/
      DATA (ASYICE3(I,13),I=1,46) /
C    BAND 13
     &8.377231e-01,8.469574e-01,8.556844e-01,8.638984e-01,8.716085e-01,
     &8.788282e-01,8.855727e-01,8.918582e-01,8.977012e-01,9.031189e-01,
     &9.081287e-01,9.127481e-01,9.169950e-01,9.208873e-01,9.244432e-01,
     &9.276807e-01,9.306183e-01,9.332744e-01,9.356674e-01,9.378159e-01,
     &9.397385e-01,9.414539e-01,9.429808e-01,9.443380e-01,9.455443e-01,
     &9.466185e-01,9.475794e-01,9.484460e-01,9.492372e-01,9.499718e-01,
     &9.506689e-01,9.513473e-01,9.520260e-01,9.527240e-01,9.534602e-01,
     &9.542535e-01,9.551230e-01,9.560875e-01,9.571659e-01,9.583773e-01,
     &9.597404e-01,9.612742e-01,9.629975e-01,9.649292e-01,9.670880e-01,
     &9.694928e-01/
      DATA (ASYICE3(I,14),I=1,46) /
C    BAND 14
     &8.344994e-01,8.443682e-01,8.536201e-01,8.622890e-01,8.703985e-01,
     &8.779698e-01,8.850230e-01,8.915780e-01,8.976545e-01,9.032725e-01,
     &9.084518e-01,9.132124e-01,9.175742e-01,9.215575e-01,9.251823e-01,
     &9.284689e-01,9.314374e-01,9.341082e-01,9.365017e-01,9.386382e-01,
     &9.405382e-01,9.422220e-01,9.437101e-01,9.450231e-01,9.461814e-01,
     &9.472056e-01,9.481162e-01,9.489338e-01,9.496790e-01,9.503724e-01,
     &9.510345e-01,9.516860e-01,9.523475e-01,9.530396e-01,9.537830e-01,
     &9.545982e-01,9.555059e-01,9.565267e-01,9.576813e-01,9.589902e-01,
     &9.604741e-01,9.621535e-01,9.640490e-01,9.661813e-01,9.685709e-01,
     &9.712382e-01/
      DATA (ASYICE3(I,15),I=1,46) /
C    BAND 15
     &7.853957e-01,7.968767e-01,8.077401e-01,8.180018e-01,8.276797e-01,
     &8.367923e-01,8.453586e-01,8.533977e-01,8.609288e-01,8.679711e-01,
     &8.745440e-01,8.806667e-01,8.863587e-01,8.916394e-01,8.965282e-01,
     &9.010445e-01,9.052079e-01,9.090379e-01,9.125540e-01,9.157758e-01,
     &9.187228e-01,9.214147e-01,9.238711e-01,9.261117e-01,9.281562e-01,
     &9.300242e-01,9.317354e-01,9.333098e-01,9.347669e-01,9.361267e-01,
     &9.374088e-01,9.386331e-01,9.398195e-01,9.409877e-01,9.421576e-01,
     &9.433489e-01,9.445817e-01,9.458755e-01,9.472504e-01,9.487260e-01,
     &9.503221e-01,9.520585e-01,9.539550e-01,9.560313e-01,9.583069e-01,
     &9.608016e-01/
      DATA (ASYICE3(I,16),I=1,46) /
C    BAND 16
     &8.340752e-01,8.435170e-01,8.517487e-01,8.592064e-01,8.660387e-01,
     &8.723204e-01,8.780997e-01,8.834137e-01,8.882934e-01,8.927662e-01,
     &8.968577e-01,9.005914e-01,9.039899e-01,9.070745e-01,9.098659e-01,
     &9.123836e-01,9.146466e-01,9.166734e-01,9.184817e-01,9.200886e-01,
     &9.215109e-01,9.227648e-01,9.238661e-01,9.248304e-01,9.256727e-01,
     &9.264078e-01,9.270505e-01,9.276150e-01,9.281156e-01,9.285662e-01,
     &9.289806e-01,9.293726e-01,9.297557e-01,9.301435e-01,9.305491e-01,
     &9.309859e-01,9.314671e-01,9.320055e-01,9.326140e-01,9.333053e-01,
     &9.340919e-01,9.349861e-01,9.360000e-01,9.371451e-01,9.384329e-01,
     &9.398744e-01/

      DATA (D_EFF_M(I),I=1,300) /
     & 6.20000e-01, 3.74000e+00, 7.14000e+00, 1.06700e+01, 1.43000e+01,
     & 1.79800e+01, 2.16900e+01, 2.53000e+01, 2.86600e+01, 3.17500e+01,
     & 3.44900e+01, 3.68800e+01, 3.89300e+01, 4.07000e+01, 4.22300e+01,
     & 4.35800e+01, 4.48000e+01, 4.59400e+01, 4.70200e+01, 4.80700e+01,
     & 4.91200e+01, 5.01600e+01, 5.12200e+01, 5.23000e+01, 5.34000e+01,
     & 5.45300e+01, 5.56800e+01, 5.68600e+01, 5.80700e+01, 5.93000e+01,
     & 6.05500e+01, 6.18300e+01, 6.31300e+01, 6.44400e+01, 6.57800e+01,
     & 6.71300e+01, 6.84900e+01, 6.98800e+01, 7.12700e+01, 7.26800e+01,
     & 7.41000e+01, 7.55300e+01, 7.69700e+01, 7.84200e+01, 7.98800e+01,
     & 8.13400e+01, 8.28200e+01, 8.43000e+01, 8.58000e+01, 8.73000e+01,
     & 8.88100e+01, 9.03200e+01, 9.18500e+01, 9.33800e+01, 9.49200e+01,
     & 9.64700e+01, 9.80300e+01, 9.96000e+01, 1.01180e+02, 1.02760e+02,
     & 1.04360e+02, 1.05970e+02, 1.07590e+02, 1.09220e+02, 1.10860e+02,
     & 1.12520e+02, 1.14190e+02, 1.15870e+02, 1.17580e+02, 1.19290e+02,
     & 1.21030e+02, 1.22780e+02, 1.24560e+02, 1.26350e+02, 1.28170e+02,
     & 1.30010e+02, 1.31880e+02, 1.33770e+02, 1.35690e+02, 1.37640e+02,
     & 1.39630e+02, 1.41640e+02, 1.43690e+02, 1.45780e+02, 1.47910e+02,
     & 1.50080e+02, 1.52290e+02, 1.54550e+02, 1.56860e+02, 1.59220e+02,
     & 1.61640e+02, 1.64110e+02, 1.66650e+02, 1.69250e+02, 1.71910e+02,
     & 1.74650e+02, 1.77450e+02, 1.80340e+02, 1.83310e+02, 1.86360e+02,
     & 1.89500e+02, 1.92730e+02, 1.96050e+02, 1.99470e+02, 2.03000e+02,
     & 2.06630e+02, 2.10360e+02, 2.14210e+02, 2.18170e+02, 2.22240e+02,
     & 2.26430e+02, 2.30740e+02, 2.35160e+02, 2.39700e+02, 2.44360e+02,
     & 2.49130e+02, 2.54000e+02, 2.58990e+02, 2.64070e+02, 2.69250e+02,
     & 2.74510e+02, 2.79850e+02, 2.85260e+02, 2.90730e+02, 2.96230e+02,
     & 3.01770e+02, 3.07330e+02, 3.12890e+02, 3.18440e+02, 3.23970e+02,
     & 3.29450e+02, 3.34880e+02, 3.40240e+02, 3.45520e+02, 3.50700e+02,
     & 3.55790e+02, 3.60750e+02, 3.65600e+02, 3.70320e+02, 3.74920e+02,
     & 3.79370e+02, 3.83690e+02, 3.87880e+02, 3.91940e+02, 3.95860e+02,
     & 3.99660e+02, 4.03340e+02, 4.06900e+02, 4.10360e+02, 4.13710e+02,
     & 4.16970e+02, 4.20150e+02, 4.23240e+02, 4.26260e+02, 4.29210e+02,
     & 4.32110e+02, 4.34950e+02, 4.37750e+02, 4.40500e+02, 4.43220e+02,
     & 4.45910e+02, 4.48560e+02, 4.51200e+02, 4.53810e+02, 4.56410e+02,
     & 4.58990e+02, 4.61550e+02, 4.64110e+02, 4.66650e+02, 4.69190e+02,
     & 4.71720e+02, 4.74240e+02, 4.76760e+02, 4.79270e+02, 4.81780e+02,
     & 4.84290e+02, 4.86790e+02, 4.89290e+02, 4.91780e+02, 4.94270e+02,
     & 4.96760e+02, 4.99250e+02, 5.01730e+02, 5.04220e+02, 5.06690e+02,
     & 5.09170e+02, 5.11650e+02, 5.14120e+02, 5.16590e+02, 5.19050e+02,
     & 5.21520e+02, 5.23980e+02, 5.26440e+02, 5.28900e+02, 5.31360e+02,
     & 5.33810e+02, 5.36260e+02, 5.38710e+02, 5.41160e+02, 5.43600e+02,
     & 5.46050e+02, 5.48490e+02, 5.50920e+02, 5.53360e+02, 5.55790e+02,
     & 5.58220e+02, 5.60650e+02, 5.63080e+02, 5.65510e+02, 5.67930e+02,
     & 5.70350e+02, 5.72770e+02, 5.75180e+02, 5.77600e+02, 5.80010e+02,
     & 5.82420e+02, 5.84830e+02, 5.87230e+02, 5.89630e+02, 5.92030e+02,
     & 5.94430e+02, 5.96830e+02, 5.99220e+02, 6.01620e+02, 6.04010e+02,
     & 6.06400e+02, 6.08780e+02, 6.11170e+02, 6.13550e+02, 6.15930e+02,
     & 6.18310e+02, 6.20680e+02, 6.23050e+02, 6.25430e+02, 6.27800e+02,
     & 6.30160e+02, 6.32530e+02, 6.34890e+02, 6.37250e+02, 6.39610e+02,
     & 6.41970e+02, 6.44320e+02, 6.46680e+02, 6.49030e+02, 6.51380e+02,
     & 6.53720e+02, 6.56070e+02, 6.58410e+02, 6.60750e+02, 6.63090e+02,
     & 6.65430e+02, 6.67760e+02, 6.70100e+02, 6.72430e+02, 6.74760e+02,
     & 6.77080e+02, 6.79410e+02, 6.81730e+02, 6.84050e+02, 6.86370e+02,
     & 6.88690e+02, 6.91000e+02, 6.93320e+02, 6.95630e+02, 6.97940e+02,
     & 7.00250e+02, 7.02550e+02, 7.04850e+02, 7.07160e+02, 7.09460e+02,
     & 7.11750e+02, 7.14050e+02, 7.16340e+02, 7.18640e+02, 7.20930e+02,
     & 7.23220e+02, 7.25500e+02, 7.27790e+02, 7.30070e+02, 7.32350e+02,
     & 7.34630e+02, 7.36910e+02, 7.39180e+02, 7.41450e+02, 7.43730e+02,
     & 7.46000e+02, 7.48260e+02, 7.50530e+02, 7.52790e+02, 7.55060e+02,
     & 7.57320e+02, 7.59580e+02, 7.61830e+02, 7.64090e+02, 7.66340e+02,
     & 7.68590e+02, 7.70840e+02, 7.73090e+02, 7.75330e+02, 7.77580e+02/
      DATA (ABSICE4(I, 1),I=1,300) /
c     BAND  1
     & 5.68640e+00, 5.69170e+00, 5.80320e+00, 6.40460e+00, 7.69140e+00,
     & 8.77600e+00, 9.34800e+00, 9.76880e+00, 1.01240e+01, 1.04020e+01,
     & 1.04960e+01, 1.02370e+01, 1.00340e+01, 9.87950e+00, 9.76300e+00,
     & 9.67270e+00, 9.59880e+00, 9.53380e+00, 9.47220e+00, 9.40990e+00,
     & 9.34460e+00, 9.27460e+00, 9.19920e+00, 9.11830e+00, 9.03200e+00,
     & 8.94070e+00, 8.84510e+00, 8.74560e+00, 8.64310e+00, 8.53800e+00,
     & 8.44460e+00, 8.38840e+00, 8.32800e+00, 8.26390e+00, 8.19670e+00,
     & 8.12680e+00, 8.05460e+00, 7.98040e+00, 7.90470e+00, 7.82760e+00,
     & 7.74950e+00, 7.67050e+00, 7.59080e+00, 7.51070e+00, 7.43030e+00,
     & 7.34970e+00, 7.26900e+00, 7.18840e+00, 7.10790e+00, 7.02760e+00,
     & 6.94760e+00, 6.86790e+00, 6.78870e+00, 6.70990e+00, 6.63160e+00,
     & 6.55370e+00, 6.47650e+00, 6.39980e+00, 6.32370e+00, 6.24820e+00,
     & 6.17850e+00, 6.12870e+00, 6.07920e+00, 6.02980e+00, 5.98070e+00,
     & 5.93170e+00, 5.88300e+00, 5.83440e+00, 5.78600e+00, 5.73780e+00,
     & 5.69280e+00, 5.64800e+00, 5.60330e+00, 5.55860e+00, 5.51390e+00,
     & 5.46910e+00, 5.42440e+00, 5.37960e+00, 5.33480e+00, 5.28990e+00,
     & 5.24490e+00, 5.19980e+00, 5.15460e+00, 5.10920e+00, 5.06360e+00,
     & 5.01790e+00, 4.97200e+00, 4.92580e+00, 4.87940e+00, 4.83280e+00,
     & 4.78590e+00, 4.73870e+00, 4.69120e+00, 4.64340e+00, 4.59520e+00,
     & 4.54680e+00, 4.49800e+00, 4.44880e+00, 4.39930e+00, 4.34950e+00,
     & 4.29930e+00, 4.24880e+00, 4.19790e+00, 4.14680e+00, 4.09530e+00,
     & 4.04360e+00, 3.99170e+00, 3.93950e+00, 3.88720e+00, 3.83470e+00,
     & 3.78220e+00, 3.72960e+00, 3.67710e+00, 3.62470e+00, 3.57240e+00,
     & 3.52040e+00, 3.46870e+00, 3.41740e+00, 3.36660e+00, 3.31640e+00,
     & 3.26680e+00, 3.21800e+00, 3.17000e+00, 3.12290e+00, 3.07680e+00,
     & 3.03180e+00, 2.98790e+00, 2.94530e+00, 2.90390e+00, 2.86380e+00,
     & 2.82510e+00, 2.78780e+00, 2.75190e+00, 2.71750e+00, 2.68450e+00,
     & 2.65290e+00, 2.62270e+00, 2.59390e+00, 2.56640e+00, 2.54030e+00,
     & 2.51540e+00, 2.49170e+00, 2.46920e+00, 2.44780e+00, 2.42750e+00,
     & 2.40810e+00, 2.38960e+00, 2.37200e+00, 2.35510e+00, 2.33900e+00,
     & 2.32350e+00, 2.30870e+00, 2.29440e+00, 2.28060e+00, 2.26730e+00,
     & 2.25440e+00, 2.24180e+00, 2.22960e+00, 2.21780e+00, 2.20610e+00,
     & 2.19480e+00, 2.18370e+00, 2.17280e+00, 2.16200e+00, 2.15150e+00,
     & 2.14110e+00, 2.13080e+00, 2.12070e+00, 2.11080e+00, 2.10090e+00,
     & 2.09120e+00, 2.08150e+00, 2.07200e+00, 2.06260e+00, 2.05330e+00,
     & 2.04400e+00, 2.03490e+00, 2.02590e+00, 2.01690e+00, 2.00800e+00,
     & 1.99920e+00, 1.99050e+00, 1.98190e+00, 1.97330e+00, 1.96480e+00,
     & 1.95640e+00, 1.94810e+00, 1.93990e+00, 1.93170e+00, 1.92360e+00,
     & 1.91550e+00, 1.90760e+00, 1.89970e+00, 1.89180e+00, 1.88410e+00,
     & 1.87640e+00, 1.86880e+00, 1.86120e+00, 1.85370e+00, 1.84630e+00,
     & 1.83890e+00, 1.83160e+00, 1.82430e+00, 1.81720e+00, 1.81000e+00,
     & 1.80300e+00, 1.79590e+00, 1.78900e+00, 1.78210e+00, 1.77520e+00,
     & 1.76840e+00, 1.76170e+00, 1.75500e+00, 1.74840e+00, 1.74180e+00,
     & 1.73530e+00, 1.72880e+00, 1.72240e+00, 1.71600e+00, 1.70970e+00,
     & 1.70340e+00, 1.69720e+00, 1.69100e+00, 1.68490e+00, 1.67880e+00,
     & 1.67280e+00, 1.66680e+00, 1.66080e+00, 1.65490e+00, 1.64910e+00,
     & 1.64330e+00, 1.63750e+00, 1.63180e+00, 1.62610e+00, 1.62040e+00,
     & 1.61480e+00, 1.60930e+00, 1.60370e+00, 1.59830e+00, 1.59280e+00,
     & 1.58740e+00, 1.58200e+00, 1.57670e+00, 1.57140e+00, 1.56620e+00,
     & 1.56090e+00, 1.55580e+00, 1.55060e+00, 1.54550e+00, 1.54040e+00,
     & 1.53540e+00, 1.53040e+00, 1.52540e+00, 1.52050e+00, 1.51560e+00,
     & 1.51070e+00, 1.50590e+00, 1.50110e+00, 1.49630e+00, 1.49150e+00,
     & 1.48680e+00, 1.48220e+00, 1.47750e+00, 1.47290e+00, 1.46830e+00,
     & 1.46380e+00, 1.45920e+00, 1.45470e+00, 1.45030e+00, 1.44580e+00,
     & 1.44140e+00, 1.43700e+00, 1.43270e+00, 1.42840e+00, 1.42410e+00,
     & 1.41980e+00, 1.41550e+00, 1.41130e+00, 1.40710e+00, 1.40300e+00,
     & 1.39880e+00, 1.39470e+00, 1.39060e+00, 1.38660e+00, 1.38250e+00,
     & 1.37850e+00, 1.37450e+00, 1.37050e+00, 1.36660e+00, 1.36270e+00,
     & 1.35880e+00, 1.35490e+00, 1.35110e+00, 1.34730e+00, 1.34350e+00,
     & 1.33970e+00, 1.33590e+00, 1.33220e+00, 1.32850e+00, 1.32480e+00/
      DATA (ABSICE4(I, 2),I=1,300) /
c     BAND  2
     & 6.46430e+01, 6.54170e+01, 6.78060e+01, 6.71830e+01, 6.40400e+01,
     & 6.01070e+01, 5.58560e+01, 5.19450e+01, 4.86460e+01, 4.59370e+01,
     & 4.35810e+01, 4.12560e+01, 3.94710e+01, 3.80700e+01, 3.69390e+01,
     & 3.59940e+01, 3.51760e+01, 3.44400e+01, 3.37580e+01, 3.31090e+01,
     & 3.24790e+01, 3.18600e+01, 3.12470e+01, 3.06370e+01, 3.00310e+01,
     & 2.94280e+01, 2.88290e+01, 2.82350e+01, 2.76480e+01, 2.70680e+01,
     & 2.65070e+01, 2.59780e+01, 2.54560e+01, 2.49420e+01, 2.44380e+01,
     & 2.39430e+01, 2.34590e+01, 2.29860e+01, 2.25240e+01, 2.20730e+01,
     & 2.16350e+01, 2.12070e+01, 2.07910e+01, 2.03860e+01, 1.99930e+01,
     & 1.96100e+01, 1.92390e+01, 1.88770e+01, 1.85260e+01, 1.81850e+01,
     & 1.78530e+01, 1.75300e+01, 1.72170e+01, 1.69120e+01, 1.66150e+01,
     & 1.63270e+01, 1.60460e+01, 1.57730e+01, 1.55070e+01, 1.52470e+01,
     & 1.50010e+01, 1.47830e+01, 1.45700e+01, 1.43610e+01, 1.41570e+01,
     & 1.39560e+01, 1.37590e+01, 1.35660e+01, 1.33760e+01, 1.31890e+01,
     & 1.30090e+01, 1.28310e+01, 1.26570e+01, 1.24840e+01, 1.23140e+01,
     & 1.21470e+01, 1.19810e+01, 1.18180e+01, 1.16560e+01, 1.14960e+01,
     & 1.13380e+01, 1.11820e+01, 1.10270e+01, 1.08730e+01, 1.07210e+01,
     & 1.05690e+01, 1.04190e+01, 1.02700e+01, 1.01220e+01, 9.97500e+00,
     & 9.82870e+00, 9.68320e+00, 9.53830e+00, 9.39410e+00, 9.25060e+00,
     & 9.10770e+00, 8.96530e+00, 8.82350e+00, 8.68220e+00, 8.54150e+00,
     & 8.40130e+00, 8.26160e+00, 8.12250e+00, 7.98410e+00, 7.84630e+00,
     & 7.70910e+00, 7.57280e+00, 7.43730e+00, 7.30270e+00, 7.16910e+00,
     & 7.03670e+00, 6.90550e+00, 6.77570e+00, 6.64730e+00, 6.52060e+00,
     & 6.39570e+00, 6.27260e+00, 6.15170e+00, 6.03290e+00, 5.91660e+00,
     & 5.80270e+00, 5.69150e+00, 5.58320e+00, 5.47770e+00, 5.37530e+00,
     & 5.27610e+00, 5.18010e+00, 5.08750e+00, 4.99820e+00, 4.91240e+00,
     & 4.83010e+00, 4.75120e+00, 4.67570e+00, 4.60370e+00, 4.53510e+00,
     & 4.46970e+00, 4.40760e+00, 4.34870e+00, 4.29270e+00, 4.23960e+00,
     & 4.18930e+00, 4.14170e+00, 4.09650e+00, 4.05370e+00, 4.01310e+00,
     & 3.97460e+00, 3.93790e+00, 3.90310e+00, 3.86990e+00, 3.83810e+00,
     & 3.80780e+00, 3.77870e+00, 3.75080e+00, 3.72390e+00, 3.69800e+00,
     & 3.67290e+00, 3.64860e+00, 3.62510e+00, 3.60210e+00, 3.57980e+00,
     & 3.55790e+00, 3.53660e+00, 3.51570e+00, 3.49520e+00, 3.47510e+00,
     & 3.45530e+00, 3.43580e+00, 3.41670e+00, 3.39780e+00, 3.37920e+00,
     & 3.36080e+00, 3.34270e+00, 3.32480e+00, 3.30710e+00, 3.28960e+00,
     & 3.27240e+00, 3.25530e+00, 3.23850e+00, 3.22180e+00, 3.20530e+00,
     & 3.18900e+00, 3.17290e+00, 3.15700e+00, 3.14120e+00, 3.12560e+00,
     & 3.11020e+00, 3.09490e+00, 3.07980e+00, 3.06490e+00, 3.05010e+00,
     & 3.03540e+00, 3.02090e+00, 3.00660e+00, 2.99240e+00, 2.97830e+00,
     & 2.96440e+00, 2.95070e+00, 2.93700e+00, 2.92350e+00, 2.91020e+00,
     & 2.89690e+00, 2.88380e+00, 2.87080e+00, 2.85800e+00, 2.84530e+00,
     & 2.83270e+00, 2.82020e+00, 2.80780e+00, 2.79560e+00, 2.78340e+00,
     & 2.77140e+00, 2.75950e+00, 2.74770e+00, 2.73600e+00, 2.72440e+00,
     & 2.71300e+00, 2.70160e+00, 2.69030e+00, 2.67920e+00, 2.66810e+00,
     & 2.65720e+00, 2.64630e+00, 2.63550e+00, 2.62480e+00, 2.61430e+00,
     & 2.60380e+00, 2.59340e+00, 2.58310e+00, 2.57290e+00, 2.56270e+00,
     & 2.55270e+00, 2.54270e+00, 2.53280e+00, 2.52310e+00, 2.51340e+00,
     & 2.50370e+00, 2.49420e+00, 2.48470e+00, 2.47530e+00, 2.46600e+00,
     & 2.45680e+00, 2.44760e+00, 2.43850e+00, 2.42950e+00, 2.42060e+00,
     & 2.41170e+00, 2.40290e+00, 2.39420e+00, 2.38560e+00, 2.37700e+00,
     & 2.36850e+00, 2.36000e+00, 2.35160e+00, 2.34330e+00, 2.33500e+00,
     & 2.32690e+00, 2.31870e+00, 2.31070e+00, 2.30270e+00, 2.29470e+00,
     & 2.28680e+00, 2.27900e+00, 2.27120e+00, 2.26350e+00, 2.25590e+00,
     & 2.24830e+00, 2.24080e+00, 2.23330e+00, 2.22580e+00, 2.21850e+00,
     & 2.21120e+00, 2.20390e+00, 2.19670e+00, 2.18950e+00, 2.18240e+00,
     & 2.17540e+00, 2.16840e+00, 2.16140e+00, 2.15450e+00, 2.14760e+00,
     & 2.14080e+00, 2.13410e+00, 2.12730e+00, 2.12070e+00, 2.11400e+00,
     & 2.10750e+00, 2.10090e+00, 2.09450e+00, 2.08800e+00, 2.08160e+00,
     & 2.07530e+00, 2.06890e+00, 2.06270e+00, 2.05640e+00, 2.05030e+00,
     & 2.04410e+00, 2.03800e+00, 2.03200e+00, 2.02590e+00, 2.01990e+00/
      DATA (ABSICE4(I, 3),I=1,300) /
c     BAND  3
     & 1.55080e+01, 2.00760e+01, 2.39750e+01, 2.40710e+01, 2.41230e+01,
     & 2.41430e+01, 2.40460e+01, 2.38140e+01, 2.35220e+01, 2.32550e+01,
     & 2.28350e+01, 2.19160e+01, 2.12320e+01, 2.07200e+01, 2.03280e+01,
     & 2.00190e+01, 1.97620e+01, 1.95370e+01, 1.93290e+01, 1.91260e+01,
     & 1.89230e+01, 1.87160e+01, 1.85020e+01, 1.82800e+01, 1.80510e+01,
     & 1.78160e+01, 1.75740e+01, 1.73280e+01, 1.70790e+01, 1.68280e+01,
     & 1.65970e+01, 1.64270e+01, 1.62520e+01, 1.60730e+01, 1.58910e+01,
     & 1.57070e+01, 1.55210e+01, 1.53340e+01, 1.51480e+01, 1.49610e+01,
     & 1.47760e+01, 1.45910e+01, 1.44080e+01, 1.42260e+01, 1.40470e+01,
     & 1.38690e+01, 1.36940e+01, 1.35210e+01, 1.33500e+01, 1.31810e+01,
     & 1.30150e+01, 1.28520e+01, 1.26910e+01, 1.25320e+01, 1.23760e+01,
     & 1.22230e+01, 1.20720e+01, 1.19230e+01, 1.17760e+01, 1.16320e+01,
     & 1.14950e+01, 1.13780e+01, 1.12620e+01, 1.11480e+01, 1.10340e+01,
     & 1.09220e+01, 1.08110e+01, 1.07000e+01, 1.05900e+01, 1.04820e+01,
     & 1.03760e+01, 1.02710e+01, 1.01670e+01, 1.00630e+01, 9.96020e+00,
     & 9.85740e+00, 9.75490e+00, 9.65270e+00, 9.55080e+00, 9.44910e+00,
     & 9.34760e+00, 9.24610e+00, 9.14480e+00, 9.04340e+00, 8.94210e+00,
     & 8.84070e+00, 8.73920e+00, 8.63760e+00, 8.53570e+00, 8.43370e+00,
     & 8.33140e+00, 8.22890e+00, 8.12600e+00, 8.02290e+00, 7.91930e+00,
     & 7.81540e+00, 7.71120e+00, 7.60660e+00, 7.50160e+00, 7.39620e+00,
     & 7.29050e+00, 7.18450e+00, 7.07820e+00, 6.97160e+00, 6.86480e+00,
     & 6.75790e+00, 6.65080e+00, 6.54380e+00, 6.43680e+00, 6.33000e+00,
     & 6.22340e+00, 6.11720e+00, 6.01160e+00, 5.90650e+00, 5.80230e+00,
     & 5.69890e+00, 5.59660e+00, 5.49560e+00, 5.39590e+00, 5.29780e+00,
     & 5.20130e+00, 5.10680e+00, 5.01420e+00, 4.92380e+00, 4.83560e+00,
     & 4.74990e+00, 4.66670e+00, 4.58610e+00, 4.50830e+00, 4.43320e+00,
     & 4.36090e+00, 4.29150e+00, 4.22500e+00, 4.16140e+00, 4.10060e+00,
     & 4.04260e+00, 3.98730e+00, 3.93480e+00, 3.88490e+00, 3.83740e+00,
     & 3.79240e+00, 3.74970e+00, 3.70920e+00, 3.67080e+00, 3.63430e+00,
     & 3.59960e+00, 3.56660e+00, 3.53520e+00, 3.50520e+00, 3.47650e+00,
     & 3.44910e+00, 3.42280e+00, 3.39760e+00, 3.37320e+00, 3.34970e+00,
     & 3.32700e+00, 3.30500e+00, 3.28360e+00, 3.26270e+00, 3.24240e+00,
     & 3.22260e+00, 3.20310e+00, 3.18410e+00, 3.16540e+00, 3.14710e+00,
     & 3.12900e+00, 3.11130e+00, 3.09380e+00, 3.07650e+00, 3.05950e+00,
     & 3.04270e+00, 3.02610e+00, 3.00970e+00, 2.99350e+00, 2.97740e+00,
     & 2.96160e+00, 2.94590e+00, 2.93040e+00, 2.91510e+00, 2.90000e+00,
     & 2.88500e+00, 2.87010e+00, 2.85540e+00, 2.84090e+00, 2.82650e+00,
     & 2.81230e+00, 2.79820e+00, 2.78420e+00, 2.77040e+00, 2.75670e+00,
     & 2.74320e+00, 2.72980e+00, 2.71650e+00, 2.70330e+00, 2.69030e+00,
     & 2.67740e+00, 2.66470e+00, 2.65200e+00, 2.63950e+00, 2.62710e+00,
     & 2.61480e+00, 2.60260e+00, 2.59060e+00, 2.57860e+00, 2.56680e+00,
     & 2.55510e+00, 2.54340e+00, 2.53190e+00, 2.52050e+00, 2.50920e+00,
     & 2.49800e+00, 2.48690e+00, 2.47590e+00, 2.46500e+00, 2.45420e+00,
     & 2.44350e+00, 2.43290e+00, 2.42240e+00, 2.41190e+00, 2.40160e+00,
     & 2.39140e+00, 2.38120e+00, 2.37110e+00, 2.36120e+00, 2.35130e+00,
     & 2.34140e+00, 2.33170e+00, 2.32210e+00, 2.31250e+00, 2.30300e+00,
     & 2.29360e+00, 2.28430e+00, 2.27500e+00, 2.26580e+00, 2.25670e+00,
     & 2.24770e+00, 2.23880e+00, 2.22990e+00, 2.22110e+00, 2.21230e+00,
     & 2.20370e+00, 2.19510e+00, 2.18650e+00, 2.17810e+00, 2.16970e+00,
     & 2.16140e+00, 2.15310e+00, 2.14490e+00, 2.13680e+00, 2.12870e+00,
     & 2.12070e+00, 2.11270e+00, 2.10480e+00, 2.09700e+00, 2.08920e+00,
     & 2.08150e+00, 2.07390e+00, 2.06630e+00, 2.05880e+00, 2.05130e+00,
     & 2.04390e+00, 2.03650e+00, 2.02920e+00, 2.02190e+00, 2.01470e+00,
     & 2.00760e+00, 2.00050e+00, 1.99340e+00, 1.98640e+00, 1.97950e+00,
     & 1.97260e+00, 1.96570e+00, 1.95890e+00, 1.95220e+00, 1.94550e+00,
     & 1.93880e+00, 1.93220e+00, 1.92560e+00, 1.91910e+00, 1.91260e+00,
     & 1.90620e+00, 1.89980e+00, 1.89350e+00, 1.88720e+00, 1.88090e+00,
     & 1.87470e+00, 1.86860e+00, 1.86240e+00, 1.85640e+00, 1.85030e+00,
     & 1.84430e+00, 1.83840e+00, 1.83240e+00, 1.82650e+00, 1.82070e+00,
     & 1.81490e+00, 1.80910e+00, 1.80340e+00, 1.79770e+00, 1.79210e+00/
      DATA (ABSICE4(I, 4),I=1,300) /
c     BAND  4
     & 5.12570e+01, 7.51580e+01, 7.83780e+01, 7.38680e+01, 6.92030e+01,
     & 6.42050e+01, 5.90680e+01, 5.43470e+01, 5.04010e+01, 4.72020e+01,
     & 4.44960e+01, 4.19710e+01, 4.00380e+01, 3.85250e+01, 3.73070e+01,
     & 3.62940e+01, 3.54190e+01, 3.46370e+01, 3.39150e+01, 3.32320e+01,
     & 3.25730e+01, 3.19290e+01, 3.12940e+01, 3.06670e+01, 3.00450e+01,
     & 2.94300e+01, 2.88210e+01, 2.82190e+01, 2.76260e+01, 2.70430e+01,
     & 2.64780e+01, 2.59420e+01, 2.54150e+01, 2.48990e+01, 2.43930e+01,
     & 2.39000e+01, 2.34180e+01, 2.29480e+01, 2.24900e+01, 2.20450e+01,
     & 2.16130e+01, 2.11920e+01, 2.07840e+01, 2.03880e+01, 2.00030e+01,
     & 1.96300e+01, 1.92670e+01, 1.89160e+01, 1.85740e+01, 1.82430e+01,
     & 1.79210e+01, 1.76090e+01, 1.73050e+01, 1.70110e+01, 1.67240e+01,
     & 1.64450e+01, 1.61740e+01, 1.59110e+01, 1.56540e+01, 1.54040e+01,
     & 1.51640e+01, 1.49440e+01, 1.47290e+01, 1.45180e+01, 1.43110e+01,
     & 1.41090e+01, 1.39100e+01, 1.37140e+01, 1.35230e+01, 1.33340e+01,
     & 1.31510e+01, 1.29700e+01, 1.27920e+01, 1.26170e+01, 1.24440e+01,
     & 1.22730e+01, 1.21050e+01, 1.19380e+01, 1.17740e+01, 1.16110e+01,
     & 1.14510e+01, 1.12920e+01, 1.11340e+01, 1.09780e+01, 1.08230e+01,
     & 1.06690e+01, 1.05170e+01, 1.03660e+01, 1.02150e+01, 1.00660e+01,
     & 9.91770e+00, 9.77000e+00, 9.62320e+00, 9.47700e+00, 9.33160e+00,
     & 9.18680e+00, 9.04260e+00, 8.89890e+00, 8.75590e+00, 8.61350e+00,
     & 8.47160e+00, 8.33030e+00, 8.18970e+00, 8.04970e+00, 7.91030e+00,
     & 7.77180e+00, 7.63400e+00, 7.49710e+00, 7.36120e+00, 7.22640e+00,
     & 7.09270e+00, 6.96030e+00, 6.82930e+00, 6.69980e+00, 6.57200e+00,
     & 6.44600e+00, 6.32200e+00, 6.20010e+00, 6.08040e+00, 5.96320e+00,
     & 5.84850e+00, 5.73650e+00, 5.62730e+00, 5.52110e+00, 5.41800e+00,
     & 5.31810e+00, 5.22150e+00, 5.12830e+00, 5.03850e+00, 4.95210e+00,
     & 4.86920e+00, 4.78990e+00, 4.71400e+00, 4.64150e+00, 4.57250e+00,
     & 4.50680e+00, 4.44430e+00, 4.38500e+00, 4.32870e+00, 4.27540e+00,
     & 4.22480e+00, 4.17690e+00, 4.13150e+00, 4.08840e+00, 4.04760e+00,
     & 4.00880e+00, 3.97200e+00, 3.93700e+00, 3.90360e+00, 3.87170e+00,
     & 3.84120e+00, 3.81190e+00, 3.78390e+00, 3.75690e+00, 3.73080e+00,
     & 3.70560e+00, 3.68120e+00, 3.65750e+00, 3.63440e+00, 3.61190e+00,
     & 3.59000e+00, 3.56860e+00, 3.54750e+00, 3.52690e+00, 3.50670e+00,
     & 3.48680e+00, 3.46730e+00, 3.44800e+00, 3.42900e+00, 3.41030e+00,
     & 3.39190e+00, 3.37370e+00, 3.35570e+00, 3.33790e+00, 3.32040e+00,
     & 3.30300e+00, 3.28590e+00, 3.26900e+00, 3.25220e+00, 3.23570e+00,
     & 3.21930e+00, 3.20310e+00, 3.18710e+00, 3.17120e+00, 3.15550e+00,
     & 3.14000e+00, 3.12470e+00, 3.10950e+00, 3.09450e+00, 3.07960e+00,
     & 3.06490e+00, 3.05040e+00, 3.03590e+00, 3.02170e+00, 3.00760e+00,
     & 2.99360e+00, 2.97980e+00, 2.96610e+00, 2.95250e+00, 2.93910e+00,
     & 2.92580e+00, 2.91260e+00, 2.89960e+00, 2.88670e+00, 2.87390e+00,
     & 2.86120e+00, 2.84870e+00, 2.83630e+00, 2.82390e+00, 2.81180e+00,
     & 2.79970e+00, 2.78770e+00, 2.77590e+00, 2.76410e+00, 2.75250e+00,
     & 2.74100e+00, 2.72960e+00, 2.71820e+00, 2.70700e+00, 2.69590e+00,
     & 2.68490e+00, 2.67400e+00, 2.66320e+00, 2.65240e+00, 2.64180e+00,
     & 2.63130e+00, 2.62080e+00, 2.61050e+00, 2.60020e+00, 2.59000e+00,
     & 2.57990e+00, 2.56990e+00, 2.56000e+00, 2.55020e+00, 2.54040e+00,
     & 2.53080e+00, 2.52120e+00, 2.51170e+00, 2.50220e+00, 2.49290e+00,
     & 2.48360e+00, 2.47440e+00, 2.46530e+00, 2.45620e+00, 2.44730e+00,
     & 2.43830e+00, 2.42950e+00, 2.42070e+00, 2.41210e+00, 2.40340e+00,
     & 2.39490e+00, 2.38640e+00, 2.37800e+00, 2.36960e+00, 2.36130e+00,
     & 2.35310e+00, 2.34490e+00, 2.33680e+00, 2.32880e+00, 2.32080e+00,
     & 2.31290e+00, 2.30500e+00, 2.29720e+00, 2.28950e+00, 2.28180e+00,
     & 2.27420e+00, 2.26660e+00, 2.25910e+00, 2.25160e+00, 2.24420e+00,
     & 2.23680e+00, 2.22950e+00, 2.22230e+00, 2.21510e+00, 2.20800e+00,
     & 2.20090e+00, 2.19380e+00, 2.18690e+00, 2.17990e+00, 2.17300e+00,
     & 2.16620e+00, 2.15940e+00, 2.15260e+00, 2.14590e+00, 2.13930e+00,
     & 2.13270e+00, 2.12610e+00, 2.11960e+00, 2.11310e+00, 2.10670e+00,
     & 2.10030e+00, 2.09400e+00, 2.08770e+00, 2.08140e+00, 2.07520e+00,
     & 2.06900e+00, 2.06290e+00, 2.05680e+00, 2.05070e+00, 2.04470e+00/
      DATA (ABSICE4(I, 5),I=1,300) /
c     BAND  5
     & 1.05870e+02, 1.54530e+02, 1.43170e+02, 1.23790e+02, 1.06640e+02,
     & 9.16950e+01, 7.92810e+01, 6.95260e+01, 6.20770e+01, 5.63170e+01,
     & 5.18680e+01, 4.84540e+01, 4.58250e+01, 4.37530e+01, 4.20740e+01,
     & 4.06720e+01, 3.94640e+01, 3.83900e+01, 3.74100e+01, 3.64960e+01,
     & 3.56280e+01, 3.47950e+01, 3.39890e+01, 3.32050e+01, 3.24390e+01,
     & 3.16920e+01, 3.09620e+01, 3.02490e+01, 2.95540e+01, 2.88770e+01,
     & 2.82150e+01, 2.75650e+01, 2.69350e+01, 2.63250e+01, 2.57350e+01,
     & 2.51650e+01, 2.46140e+01, 2.40810e+01, 2.35680e+01, 2.30720e+01,
     & 2.25930e+01, 2.21320e+01, 2.16860e+01, 2.12560e+01, 2.08400e+01,
     & 2.04390e+01, 2.00520e+01, 1.96770e+01, 1.93150e+01, 1.89650e+01,
     & 1.86260e+01, 1.82980e+01, 1.79800e+01, 1.76720e+01, 1.73730e+01,
     & 1.70840e+01, 1.68020e+01, 1.65290e+01, 1.62630e+01, 1.60050e+01,
     & 1.57550e+01, 1.55180e+01, 1.52860e+01, 1.50590e+01, 1.48380e+01,
     & 1.46210e+01, 1.44090e+01, 1.42010e+01, 1.39970e+01, 1.37970e+01,
     & 1.36020e+01, 1.34100e+01, 1.32210e+01, 1.30350e+01, 1.28520e+01,
     & 1.26720e+01, 1.24950e+01, 1.23200e+01, 1.21480e+01, 1.19770e+01,
     & 1.18090e+01, 1.16430e+01, 1.14790e+01, 1.13160e+01, 1.11550e+01,
     & 1.09950e+01, 1.08370e+01, 1.06800e+01, 1.05250e+01, 1.03700e+01,
     & 1.02170e+01, 1.00640e+01, 9.91310e+00, 9.76250e+00, 9.61280e+00,
     & 9.46380e+00, 9.31560e+00, 9.16810e+00, 9.02140e+00, 8.87530e+00,
     & 8.72990e+00, 8.58520e+00, 8.44120e+00, 8.29800e+00, 8.15550e+00,
     & 8.01390e+00, 7.87310e+00, 7.73340e+00, 7.59460e+00, 7.45700e+00,
     & 7.32060e+00, 7.18560e+00, 7.05200e+00, 6.92000e+00, 6.78980e+00,
     & 6.66140e+00, 6.53510e+00, 6.41090e+00, 6.28900e+00, 6.16960e+00,
     & 6.05280e+00, 5.93880e+00, 5.82770e+00, 5.71960e+00, 5.61470e+00,
     & 5.51300e+00, 5.41470e+00, 5.31980e+00, 5.22840e+00, 5.14050e+00,
     & 5.05620e+00, 4.97540e+00, 4.89810e+00, 4.82440e+00, 4.75410e+00,
     & 4.68730e+00, 4.62370e+00, 4.56330e+00, 4.50600e+00, 4.45170e+00,
     & 4.40020e+00, 4.35140e+00, 4.30520e+00, 4.26140e+00, 4.21980e+00,
     & 4.18040e+00, 4.14290e+00, 4.10720e+00, 4.07320e+00, 4.04070e+00,
     & 4.00970e+00, 3.97990e+00, 3.95130e+00, 3.92380e+00, 3.89720e+00,
     & 3.87160e+00, 3.84670e+00, 3.82260e+00, 3.79910e+00, 3.77620e+00,
     & 3.75390e+00, 3.73200e+00, 3.71060e+00, 3.68960e+00, 3.66900e+00,
     & 3.64880e+00, 3.62880e+00, 3.60920e+00, 3.58990e+00, 3.57080e+00,
     & 3.55200e+00, 3.53350e+00, 3.51520e+00, 3.49710e+00, 3.47920e+00,
     & 3.46150e+00, 3.44410e+00, 3.42680e+00, 3.40980e+00, 3.39290e+00,
     & 3.37620e+00, 3.35970e+00, 3.34340e+00, 3.32730e+00, 3.31130e+00,
     & 3.29550e+00, 3.27990e+00, 3.26440e+00, 3.24910e+00, 3.23390e+00,
     & 3.21900e+00, 3.20410e+00, 3.18940e+00, 3.17490e+00, 3.16050e+00,
     & 3.14630e+00, 3.13220e+00, 3.11820e+00, 3.10440e+00, 3.09070e+00,
     & 3.07710e+00, 3.06370e+00, 3.05040e+00, 3.03730e+00, 3.02420e+00,
     & 3.01130e+00, 2.99850e+00, 2.98590e+00, 2.97330e+00, 2.96090e+00,
     & 2.94860e+00, 2.93640e+00, 2.92430e+00, 2.91240e+00, 2.90050e+00,
     & 2.88870e+00, 2.87710e+00, 2.86550e+00, 2.85410e+00, 2.84280e+00,
     & 2.83150e+00, 2.82040e+00, 2.80940e+00, 2.79840e+00, 2.78760e+00,
     & 2.77680e+00, 2.76620e+00, 2.75560e+00, 2.74510e+00, 2.73480e+00,
     & 2.72450e+00, 2.71430e+00, 2.70410e+00, 2.69410e+00, 2.68420e+00,
     & 2.67430e+00, 2.66450e+00, 2.65480e+00, 2.64520e+00, 2.63560e+00,
     & 2.62620e+00, 2.61680e+00, 2.60750e+00, 2.59820e+00, 2.58910e+00,
     & 2.58000e+00, 2.57090e+00, 2.56200e+00, 2.55310e+00, 2.54430e+00,
     & 2.53560e+00, 2.52690e+00, 2.51830e+00, 2.50980e+00, 2.50130e+00,
     & 2.49290e+00, 2.48450e+00, 2.47630e+00, 2.46810e+00, 2.45990e+00,
     & 2.45180e+00, 2.44380e+00, 2.43580e+00, 2.42790e+00, 2.42000e+00,
     & 2.41230e+00, 2.40450e+00, 2.39680e+00, 2.38920e+00, 2.38160e+00,
     & 2.37410e+00, 2.36670e+00, 2.35930e+00, 2.35190e+00, 2.34460e+00,
     & 2.33740e+00, 2.33020e+00, 2.32300e+00, 2.31590e+00, 2.30890e+00,
     & 2.30190e+00, 2.29490e+00, 2.28800e+00, 2.28120e+00, 2.27440e+00,
     & 2.26760e+00, 2.26090e+00, 2.25430e+00, 2.24760e+00, 2.24110e+00,
     & 2.23450e+00, 2.22800e+00, 2.22160e+00, 2.21520e+00, 2.20880e+00,
     & 2.20250e+00, 2.19620e+00, 2.19000e+00, 2.18380e+00, 2.17770e+00/
      DATA (ABSICE4(I, 6),I=1,300) /
c     BAND  6
     & 2.65420e+02, 2.95170e+02, 2.20410e+02, 1.67070e+02, 1.30490e+02,
     & 1.04820e+02, 8.65120e+01, 7.34310e+01, 6.40540e+01, 5.71290e+01,
     & 5.20490e+01, 4.84190e+01, 4.56370e+01, 4.34540e+01, 4.16950e+01,
     & 4.02360e+01, 3.89870e+01, 3.78870e+01, 3.68910e+01, 3.59690e+01,
     & 3.51020e+01, 3.42740e+01, 3.34780e+01, 3.27080e+01, 3.19590e+01,
     & 3.12310e+01, 3.05210e+01, 2.98310e+01, 2.91590e+01, 2.85050e+01,
     & 2.78670e+01, 2.72390e+01, 2.66320e+01, 2.60450e+01, 2.54770e+01,
     & 2.49300e+01, 2.44010e+01, 2.38900e+01, 2.33970e+01, 2.29210e+01,
     & 2.24620e+01, 2.20190e+01, 2.15900e+01, 2.11770e+01, 2.07770e+01,
     & 2.03900e+01, 2.00170e+01, 1.96550e+01, 1.93050e+01, 1.89660e+01,
     & 1.86380e+01, 1.83190e+01, 1.80110e+01, 1.77110e+01, 1.74200e+01,
     & 1.71380e+01, 1.68630e+01, 1.65950e+01, 1.63350e+01, 1.60820e+01,
     & 1.58360e+01, 1.55970e+01, 1.53650e+01, 1.51370e+01, 1.49150e+01,
     & 1.46980e+01, 1.44850e+01, 1.42760e+01, 1.40710e+01, 1.38710e+01,
     & 1.36740e+01, 1.34810e+01, 1.32910e+01, 1.31050e+01, 1.29210e+01,
     & 1.27400e+01, 1.25620e+01, 1.23860e+01, 1.22130e+01, 1.20420e+01,
     & 1.18730e+01, 1.17060e+01, 1.15410e+01, 1.13780e+01, 1.12170e+01,
     & 1.10570e+01, 1.08980e+01, 1.07410e+01, 1.05840e+01, 1.04300e+01,
     & 1.02760e+01, 1.01230e+01, 9.97100e+00, 9.82000e+00, 9.66990e+00,
     & 9.52050e+00, 9.37190e+00, 9.22410e+00, 9.07690e+00, 8.93040e+00,
     & 8.78460e+00, 8.63960e+00, 8.49520e+00, 8.35160e+00, 8.20880e+00,
     & 8.06680e+00, 7.92560e+00, 7.78550e+00, 7.64640e+00, 7.50840e+00,
     & 7.37160e+00, 7.23630e+00, 7.10230e+00, 6.97000e+00, 6.83940e+00,
     & 6.71060e+00, 6.58390e+00, 6.45940e+00, 6.33720e+00, 6.21740e+00,
     & 6.10030e+00, 5.98590e+00, 5.87450e+00, 5.76610e+00, 5.66080e+00,
     & 5.55880e+00, 5.46020e+00, 5.36500e+00, 5.27330e+00, 5.18510e+00,
     & 5.10050e+00, 5.01950e+00, 4.94200e+00, 4.86800e+00, 4.79750e+00,
     & 4.73030e+00, 4.66650e+00, 4.60600e+00, 4.54850e+00, 4.49400e+00,
     & 4.44230e+00, 4.39330e+00, 4.34690e+00, 4.30300e+00, 4.26120e+00,
     & 4.22160e+00, 4.18400e+00, 4.14820e+00, 4.11410e+00, 4.08150e+00,
     & 4.05030e+00, 4.02040e+00, 3.99170e+00, 3.96410e+00, 3.93740e+00,
     & 3.91170e+00, 3.88670e+00, 3.86250e+00, 3.83890e+00, 3.81590e+00,
     & 3.79350e+00, 3.77160e+00, 3.75010e+00, 3.72900e+00, 3.70830e+00,
     & 3.68800e+00, 3.66800e+00, 3.64830e+00, 3.62890e+00, 3.60980e+00,
     & 3.59090e+00, 3.57230e+00, 3.55390e+00, 3.53570e+00, 3.51780e+00,
     & 3.50000e+00, 3.48250e+00, 3.46520e+00, 3.44810e+00, 3.43110e+00,
     & 3.41440e+00, 3.39780e+00, 3.38140e+00, 3.36520e+00, 3.34920e+00,
     & 3.33330e+00, 3.31760e+00, 3.30210e+00, 3.28670e+00, 3.27150e+00,
     & 3.25640e+00, 3.24150e+00, 3.22680e+00, 3.21220e+00, 3.19770e+00,
     & 3.18340e+00, 3.16930e+00, 3.15530e+00, 3.14140e+00, 3.12760e+00,
     & 3.11400e+00, 3.10050e+00, 3.08720e+00, 3.07400e+00, 3.06090e+00,
     & 3.04790e+00, 3.03510e+00, 3.02240e+00, 3.00980e+00, 2.99730e+00,
     & 2.98490e+00, 2.97270e+00, 2.96050e+00, 2.94850e+00, 2.93660e+00,
     & 2.92480e+00, 2.91310e+00, 2.90150e+00, 2.89000e+00, 2.87860e+00,
     & 2.86730e+00, 2.85610e+00, 2.84510e+00, 2.83410e+00, 2.82320e+00,
     & 2.81240e+00, 2.80170e+00, 2.79110e+00, 2.78050e+00, 2.77010e+00,
     & 2.75980e+00, 2.74950e+00, 2.73930e+00, 2.72930e+00, 2.71930e+00,
     & 2.70940e+00, 2.69950e+00, 2.68980e+00, 2.68010e+00, 2.67050e+00,
     & 2.66100e+00, 2.65160e+00, 2.64220e+00, 2.63290e+00, 2.62370e+00,
     & 2.61460e+00, 2.60550e+00, 2.59650e+00, 2.58760e+00, 2.57870e+00,
     & 2.57000e+00, 2.56130e+00, 2.55260e+00, 2.54400e+00, 2.53550e+00,
     & 2.52710e+00, 2.51870e+00, 2.51040e+00, 2.50210e+00, 2.49390e+00,
     & 2.48580e+00, 2.47770e+00, 2.46970e+00, 2.46180e+00, 2.45390e+00,
     & 2.44600e+00, 2.43830e+00, 2.43050e+00, 2.42290e+00, 2.41530e+00,
     & 2.40770e+00, 2.40020e+00, 2.39280e+00, 2.38540e+00, 2.37800e+00,
     & 2.37080e+00, 2.36350e+00, 2.35630e+00, 2.34920e+00, 2.34210e+00,
     & 2.33510e+00, 2.32810e+00, 2.32120e+00, 2.31430e+00, 2.30750e+00,
     & 2.30070e+00, 2.29390e+00, 2.28720e+00, 2.28060e+00, 2.27390e+00,
     & 2.26740e+00, 2.26090e+00, 2.25440e+00, 2.24790e+00, 2.24160e+00,
     & 2.23520e+00, 2.22890e+00, 2.22260e+00, 2.21640e+00, 2.21020e+00/
      DATA (ABSICE4(I, 7),I=1,300) /
c     BAND  7
     & 2.45570e+02, 2.16460e+02, 1.63790e+02, 1.29070e+02, 1.04830e+02,
     & 8.70320e+01, 7.38120e+01, 6.40900e+01, 5.69860e+01, 5.16660e+01,
     & 4.76520e+01, 4.45840e+01, 4.22400e+01, 4.04080e+01, 3.89380e+01,
     & 3.77210e+01, 3.66800e+01, 3.57610e+01, 3.49250e+01, 3.41470e+01,
     & 3.34090e+01, 3.26980e+01, 3.20090e+01, 3.13360e+01, 3.06770e+01,
     & 3.00310e+01, 2.93970e+01, 2.87760e+01, 2.81670e+01, 2.75730e+01,
     & 2.69930e+01, 2.64280e+01, 2.58770e+01, 2.53420e+01, 2.48220e+01,
     & 2.43160e+01, 2.38260e+01, 2.33500e+01, 2.28890e+01, 2.24420e+01,
     & 2.20090e+01, 2.15890e+01, 2.11820e+01, 2.07880e+01, 2.04060e+01,
     & 2.00350e+01, 1.96760e+01, 1.93280e+01, 1.89910e+01, 1.86630e+01,
     & 1.83450e+01, 1.80360e+01, 1.77360e+01, 1.74440e+01, 1.71600e+01,
     & 1.68840e+01, 1.66160e+01, 1.63540e+01, 1.60990e+01, 1.58500e+01,
     & 1.56080e+01, 1.53740e+01, 1.51440e+01, 1.49200e+01, 1.47010e+01,
     & 1.44860e+01, 1.42750e+01, 1.40690e+01, 1.38670e+01, 1.36690e+01,
     & 1.34740e+01, 1.32830e+01, 1.30950e+01, 1.29100e+01, 1.27280e+01,
     & 1.25490e+01, 1.23730e+01, 1.21990e+01, 1.20270e+01, 1.18570e+01,
     & 1.16900e+01, 1.15250e+01, 1.13610e+01, 1.11990e+01, 1.10390e+01,
     & 1.08800e+01, 1.07220e+01, 1.05660e+01, 1.04110e+01, 1.02570e+01,
     & 1.01050e+01, 9.95290e+00, 9.80210e+00, 9.65220e+00, 9.50310e+00,
     & 9.35480e+00, 9.20730e+00, 9.06050e+00, 8.91430e+00, 8.76890e+00,
     & 8.62410e+00, 8.48010e+00, 8.33670e+00, 8.19410e+00, 8.05230e+00,
     & 7.91130e+00, 7.77120e+00, 7.63210e+00, 7.49400e+00, 7.35700e+00,
     & 7.22130e+00, 7.08690e+00, 6.95400e+00, 6.82270e+00, 6.69310e+00,
     & 6.56540e+00, 6.43970e+00, 6.31610e+00, 6.19490e+00, 6.07610e+00,
     & 5.96000e+00, 5.84660e+00, 5.73610e+00, 5.62860e+00, 5.52420e+00,
     & 5.42310e+00, 5.32540e+00, 5.23100e+00, 5.14020e+00, 5.05280e+00,
     & 4.96900e+00, 4.88870e+00, 4.81190e+00, 4.73860e+00, 4.66880e+00,
     & 4.60230e+00, 4.53920e+00, 4.47920e+00, 4.42230e+00, 4.36830e+00,
     & 4.31720e+00, 4.26870e+00, 4.22280e+00, 4.17920e+00, 4.13800e+00,
     & 4.09880e+00, 4.06150e+00, 4.02610e+00, 3.99230e+00, 3.96010e+00,
     & 3.92920e+00, 3.89970e+00, 3.87130e+00, 3.84390e+00, 3.81760e+00,
     & 3.79210e+00, 3.76740e+00, 3.74340e+00, 3.72010e+00, 3.69740e+00,
     & 3.67520e+00, 3.65350e+00, 3.63230e+00, 3.61140e+00, 3.59100e+00,
     & 3.57090e+00, 3.55110e+00, 3.53160e+00, 3.51240e+00, 3.49350e+00,
     & 3.47490e+00, 3.45650e+00, 3.43830e+00, 3.42030e+00, 3.40260e+00,
     & 3.38510e+00, 3.36770e+00, 3.35060e+00, 3.33370e+00, 3.31690e+00,
     & 3.30040e+00, 3.28400e+00, 3.26780e+00, 3.25180e+00, 3.23600e+00,
     & 3.22030e+00, 3.20480e+00, 3.18940e+00, 3.17420e+00, 3.15920e+00,
     & 3.14430e+00, 3.12960e+00, 3.11510e+00, 3.10060e+00, 3.08640e+00,
     & 3.07220e+00, 3.05830e+00, 3.04440e+00, 3.03070e+00, 3.01710e+00,
     & 3.00370e+00, 2.99040e+00, 2.97720e+00, 2.96420e+00, 2.95120e+00,
     & 2.93840e+00, 2.92580e+00, 2.91320e+00, 2.90080e+00, 2.88840e+00,
     & 2.87620e+00, 2.86410e+00, 2.85220e+00, 2.84030e+00, 2.82850e+00,
     & 2.81690e+00, 2.80530e+00, 2.79390e+00, 2.78260e+00, 2.77130e+00,
     & 2.76020e+00, 2.74920e+00, 2.73820e+00, 2.72740e+00, 2.71660e+00,
     & 2.70600e+00, 2.69540e+00, 2.68490e+00, 2.67460e+00, 2.66430e+00,
     & 2.65410e+00, 2.64400e+00, 2.63390e+00, 2.62400e+00, 2.61410e+00,
     & 2.60440e+00, 2.59470e+00, 2.58500e+00, 2.57550e+00, 2.56610e+00,
     & 2.55670e+00, 2.54740e+00, 2.53810e+00, 2.52900e+00, 2.51990e+00,
     & 2.51090e+00, 2.50200e+00, 2.49310e+00, 2.48430e+00, 2.47560e+00,
     & 2.46690e+00, 2.45840e+00, 2.44980e+00, 2.44140e+00, 2.43300e+00,
     & 2.42470e+00, 2.41640e+00, 2.40820e+00, 2.40010e+00, 2.39200e+00,
     & 2.38400e+00, 2.37600e+00, 2.36820e+00, 2.36030e+00, 2.35260e+00,
     & 2.34480e+00, 2.33720e+00, 2.32960e+00, 2.32200e+00, 2.31450e+00,
     & 2.30710e+00, 2.29970e+00, 2.29240e+00, 2.28510e+00, 2.27790e+00,
     & 2.27070e+00, 2.26360e+00, 2.25650e+00, 2.24950e+00, 2.24250e+00,
     & 2.23560e+00, 2.22870e+00, 2.22190e+00, 2.21510e+00, 2.20840e+00,
     & 2.20170e+00, 2.19510e+00, 2.18850e+00, 2.18190e+00, 2.17540e+00,
     & 2.16890e+00, 2.16250e+00, 2.15620e+00, 2.14980e+00, 2.14350e+00,
     & 2.13730e+00, 2.13110e+00, 2.12490e+00, 2.11880e+00, 2.11270e+00/
      DATA (ABSICE4(I, 8),I=1,300) /
c     BAND  8
     & 4.86200e+01, 6.25690e+01, 6.17720e+01, 6.06430e+01, 5.75190e+01,
     & 5.32040e+01, 4.90490e+01, 4.55030e+01, 4.26290e+01, 4.03290e+01,
     & 3.83350e+01, 3.63180e+01, 3.47900e+01, 3.36090e+01, 3.26710e+01,
     & 3.19010e+01, 3.12430e+01, 3.06580e+01, 3.01200e+01, 2.96090e+01,
     & 2.91140e+01, 2.86260e+01, 2.81420e+01, 2.76580e+01, 2.71760e+01,
     & 2.66930e+01, 2.62120e+01, 2.57340e+01, 2.52600e+01, 2.47900e+01,
     & 2.43370e+01, 2.39210e+01, 2.35080e+01, 2.31000e+01, 2.26980e+01,
     & 2.23020e+01, 2.19130e+01, 2.15310e+01, 2.11570e+01, 2.07910e+01,
     & 2.04340e+01, 2.00840e+01, 1.97430e+01, 1.94100e+01, 1.90850e+01,
     & 1.87680e+01, 1.84590e+01, 1.81580e+01, 1.78640e+01, 1.75770e+01,
     & 1.72980e+01, 1.70250e+01, 1.67600e+01, 1.65000e+01, 1.62470e+01,
     & 1.60000e+01, 1.57580e+01, 1.55230e+01, 1.52920e+01, 1.50670e+01,
     & 1.48480e+01, 1.46380e+01, 1.44330e+01, 1.42320e+01, 1.40340e+01,
     & 1.38390e+01, 1.36480e+01, 1.34600e+01, 1.32760e+01, 1.30940e+01,
     & 1.29160e+01, 1.27400e+01, 1.25670e+01, 1.23960e+01, 1.22270e+01,
     & 1.20600e+01, 1.18960e+01, 1.17330e+01, 1.15730e+01, 1.14130e+01,
     & 1.12560e+01, 1.11000e+01, 1.09450e+01, 1.07920e+01, 1.06400e+01,
     & 1.04890e+01, 1.03390e+01, 1.01910e+01, 1.00430e+01, 9.89580e+00,
     & 9.74960e+00, 9.60410e+00, 9.45940e+00, 9.31520e+00, 9.17180e+00,
     & 9.02880e+00, 8.88650e+00, 8.74470e+00, 8.60340e+00, 8.46270e+00,
     & 8.32250e+00, 8.18280e+00, 8.04380e+00, 7.90530e+00, 7.76750e+00,
     & 7.63050e+00, 7.49420e+00, 7.35870e+00, 7.22420e+00, 7.09070e+00,
     & 6.95840e+00, 6.82730e+00, 6.69760e+00, 6.56940e+00, 6.44280e+00,
     & 6.31810e+00, 6.19520e+00, 6.07450e+00, 5.95600e+00, 5.83980e+00,
     & 5.72620e+00, 5.61530e+00, 5.50720e+00, 5.40200e+00, 5.29990e+00,
     & 5.20100e+00, 5.10530e+00, 5.01290e+00, 4.92400e+00, 4.83840e+00,
     & 4.75640e+00, 4.67780e+00, 4.60260e+00, 4.53090e+00, 4.46250e+00,
     & 4.39740e+00, 4.33560e+00, 4.27690e+00, 4.22110e+00, 4.16830e+00,
     & 4.11830e+00, 4.07080e+00, 4.02590e+00, 3.98330e+00, 3.94290e+00,
     & 3.90450e+00, 3.86800e+00, 3.83340e+00, 3.80030e+00, 3.76870e+00,
     & 3.73860e+00, 3.70960e+00, 3.68190e+00, 3.65510e+00, 3.62930e+00,
     & 3.60440e+00, 3.58020e+00, 3.55680e+00, 3.53400e+00, 3.51170e+00,
     & 3.49000e+00, 3.46880e+00, 3.44800e+00, 3.42760e+00, 3.40760e+00,
     & 3.38790e+00, 3.36860e+00, 3.34950e+00, 3.33070e+00, 3.31220e+00,
     & 3.29400e+00, 3.27600e+00, 3.25820e+00, 3.24060e+00, 3.22320e+00,
     & 3.20610e+00, 3.18910e+00, 3.17240e+00, 3.15580e+00, 3.13940e+00,
     & 3.12320e+00, 3.10720e+00, 3.09140e+00, 3.07570e+00, 3.06020e+00,
     & 3.04490e+00, 3.02970e+00, 3.01470e+00, 2.99980e+00, 2.98510e+00,
     & 2.97060e+00, 2.95620e+00, 2.94190e+00, 2.92780e+00, 2.91390e+00,
     & 2.90010e+00, 2.88640e+00, 2.87280e+00, 2.85940e+00, 2.84610e+00,
     & 2.83300e+00, 2.82000e+00, 2.80710e+00, 2.79430e+00, 2.78170e+00,
     & 2.76920e+00, 2.75680e+00, 2.74450e+00, 2.73230e+00, 2.72030e+00,
     & 2.70840e+00, 2.69650e+00, 2.68480e+00, 2.67320e+00, 2.66170e+00,
     & 2.65030e+00, 2.63910e+00, 2.62790e+00, 2.61680e+00, 2.60580e+00,
     & 2.59490e+00, 2.58410e+00, 2.57340e+00, 2.56280e+00, 2.55230e+00,
     & 2.54190e+00, 2.53160e+00, 2.52140e+00, 2.51120e+00, 2.50120e+00,
     & 2.49120e+00, 2.48140e+00, 2.47160e+00, 2.46180e+00, 2.45220e+00,
     & 2.44270e+00, 2.43320e+00, 2.42380e+00, 2.41450e+00, 2.40520e+00,
     & 2.39610e+00, 2.38700e+00, 2.37800e+00, 2.36910e+00, 2.36020e+00,
     & 2.35140e+00, 2.34270e+00, 2.33400e+00, 2.32540e+00, 2.31690e+00,
     & 2.30850e+00, 2.30010e+00, 2.29180e+00, 2.28350e+00, 2.27530e+00,
     & 2.26720e+00, 2.25920e+00, 2.25120e+00, 2.24320e+00, 2.23530e+00,
     & 2.22750e+00, 2.21980e+00, 2.21210e+00, 2.20440e+00, 2.19680e+00,
     & 2.18930e+00, 2.18180e+00, 2.17440e+00, 2.16710e+00, 2.15980e+00,
     & 2.15250e+00, 2.14530e+00, 2.13820e+00, 2.13110e+00, 2.12400e+00,
     & 2.11700e+00, 2.11010e+00, 2.10320e+00, 2.09630e+00, 2.08950e+00,
     & 2.08280e+00, 2.07610e+00, 2.06940e+00, 2.06280e+00, 2.05630e+00,
     & 2.04970e+00, 2.04330e+00, 2.03680e+00, 2.03050e+00, 2.02410e+00,
     & 2.01780e+00, 2.01160e+00, 2.00540e+00, 1.99920e+00, 1.99310e+00,
     & 1.98700e+00, 1.98090e+00, 1.97490e+00, 1.96900e+00, 1.96300e+00/
      DATA (ABSICE4(I, 9),I=1,300) /
c     BAND  9
     & 4.67570e+01, 6.31640e+01, 6.32440e+01, 6.20880e+01, 5.81630e+01,
     & 5.35130e+01, 4.92800e+01, 4.56820e+01, 4.27650e+01, 4.04290e+01,
     & 3.84100e+01, 3.63780e+01, 3.48380e+01, 3.36470e+01, 3.27020e+01,
     & 3.19240e+01, 3.12600e+01, 3.06700e+01, 3.01260e+01, 2.96110e+01,
     & 2.91120e+01, 2.86200e+01, 2.81330e+01, 2.76470e+01, 2.71610e+01,
     & 2.66770e+01, 2.61950e+01, 2.57150e+01, 2.52390e+01, 2.47680e+01,
     & 2.43150e+01, 2.38980e+01, 2.34840e+01, 2.30760e+01, 2.26730e+01,
     & 2.22770e+01, 2.18880e+01, 2.15070e+01, 2.11330e+01, 2.07680e+01,
     & 2.04110e+01, 2.00620e+01, 1.97210e+01, 1.93890e+01, 1.90650e+01,
     & 1.87490e+01, 1.84410e+01, 1.81410e+01, 1.78480e+01, 1.75620e+01,
     & 1.72840e+01, 1.70120e+01, 1.67470e+01, 1.64890e+01, 1.62360e+01,
     & 1.59900e+01, 1.57500e+01, 1.55150e+01, 1.52850e+01, 1.50600e+01,
     & 1.48420e+01, 1.46330e+01, 1.44280e+01, 1.42270e+01, 1.40300e+01,
     & 1.38360e+01, 1.36450e+01, 1.34580e+01, 1.32740e+01, 1.30920e+01,
     & 1.29140e+01, 1.27390e+01, 1.25660e+01, 1.23950e+01, 1.22270e+01,
     & 1.20610e+01, 1.18960e+01, 1.17340e+01, 1.15730e+01, 1.14140e+01,
     & 1.12570e+01, 1.11010e+01, 1.09470e+01, 1.07940e+01, 1.06420e+01,
     & 1.04910e+01, 1.03420e+01, 1.01930e+01, 1.00450e+01, 9.89820e+00,
     & 9.75210e+00, 9.60680e+00, 9.46210e+00, 9.31810e+00, 9.17470e+00,
     & 9.03190e+00, 8.88960e+00, 8.74790e+00, 8.60670e+00, 8.46600e+00,
     & 8.32590e+00, 8.18630e+00, 8.04720e+00, 7.90880e+00, 7.77110e+00,
     & 7.63410e+00, 7.49780e+00, 7.36240e+00, 7.22790e+00, 7.09450e+00,
     & 6.96210e+00, 6.83110e+00, 6.70140e+00, 6.57320e+00, 6.44670e+00,
     & 6.32190e+00, 6.19910e+00, 6.07830e+00, 5.95980e+00, 5.84370e+00,
     & 5.73010e+00, 5.61920e+00, 5.51100e+00, 5.40590e+00, 5.30380e+00,
     & 5.20480e+00, 5.10910e+00, 5.01680e+00, 4.92780e+00, 4.84230e+00,
     & 4.76020e+00, 4.68160e+00, 4.60640e+00, 4.53470e+00, 4.46630e+00,
     & 4.40120e+00, 4.33940e+00, 4.28060e+00, 4.22490e+00, 4.17210e+00,
     & 4.12200e+00, 4.07450e+00, 4.02960e+00, 3.98700e+00, 3.94650e+00,
     & 3.90820e+00, 3.87170e+00, 3.83700e+00, 3.80400e+00, 3.77240e+00,
     & 3.74220e+00, 3.71330e+00, 3.68550e+00, 3.65880e+00, 3.63300e+00,
     & 3.60800e+00, 3.58380e+00, 3.56040e+00, 3.53760e+00, 3.51530e+00,
     & 3.49360e+00, 3.47240e+00, 3.45160e+00, 3.43120e+00, 3.41120e+00,
     & 3.39150e+00, 3.37210e+00, 3.35310e+00, 3.33430e+00, 3.31580e+00,
     & 3.29750e+00, 3.27950e+00, 3.26170e+00, 3.24410e+00, 3.22680e+00,
     & 3.20960e+00, 3.19270e+00, 3.17590e+00, 3.15930e+00, 3.14300e+00,
     & 3.12680e+00, 3.11070e+00, 3.09490e+00, 3.07920e+00, 3.06370e+00,
     & 3.04840e+00, 3.03320e+00, 3.01820e+00, 3.00330e+00, 2.98860e+00,
     & 2.97400e+00, 2.95960e+00, 2.94540e+00, 2.93130e+00, 2.91730e+00,
     & 2.90350e+00, 2.88980e+00, 2.87630e+00, 2.86290e+00, 2.84960e+00,
     & 2.83640e+00, 2.82340e+00, 2.81050e+00, 2.79770e+00, 2.78510e+00,
     & 2.77260e+00, 2.76020e+00, 2.74790e+00, 2.73570e+00, 2.72370e+00,
     & 2.71170e+00, 2.69990e+00, 2.68820e+00, 2.67660e+00, 2.66510e+00,
     & 2.65370e+00, 2.64240e+00, 2.63120e+00, 2.62010e+00, 2.60910e+00,
     & 2.59830e+00, 2.58750e+00, 2.57680e+00, 2.56620e+00, 2.55570e+00,
     & 2.54530e+00, 2.53490e+00, 2.52470e+00, 2.51460e+00, 2.50450e+00,
     & 2.49450e+00, 2.48460e+00, 2.47480e+00, 2.46510e+00, 2.45550e+00,
     & 2.44590e+00, 2.43650e+00, 2.42710e+00, 2.41780e+00, 2.40850e+00,
     & 2.39930e+00, 2.39030e+00, 2.38120e+00, 2.37230e+00, 2.36340e+00,
     & 2.35460e+00, 2.34590e+00, 2.33730e+00, 2.32870e+00, 2.32010e+00,
     & 2.31170e+00, 2.30330e+00, 2.29500e+00, 2.28670e+00, 2.27850e+00,
     & 2.27040e+00, 2.26240e+00, 2.25440e+00, 2.24640e+00, 2.23850e+00,
     & 2.23070e+00, 2.22300e+00, 2.21520e+00, 2.20760e+00, 2.20000e+00,
     & 2.19250e+00, 2.18500e+00, 2.17760e+00, 2.17020e+00, 2.16290e+00,
     & 2.15570e+00, 2.14850e+00, 2.14130e+00, 2.13420e+00, 2.12720e+00,
     & 2.12020e+00, 2.11320e+00, 2.10630e+00, 2.09950e+00, 2.09270e+00,
     & 2.08590e+00, 2.07920e+00, 2.07250e+00, 2.06590e+00, 2.05940e+00,
     & 2.05290e+00, 2.04640e+00, 2.03990e+00, 2.03360e+00, 2.02720e+00,
     & 2.02090e+00, 2.01470e+00, 2.00840e+00, 2.00230e+00, 1.99610e+00,
     & 1.99010e+00, 1.98400e+00, 1.97800e+00, 1.97200e+00, 1.96610e+00/
      DATA (ABSICE4(I,10),I=1,300) /
c     BAND 10
     & 6.21430e+01, 8.45450e+01, 8.31970e+01, 7.82550e+01, 7.05800e+01,
     & 6.34110e+01, 5.72600e+01, 5.21860e+01, 4.81470e+01, 4.49400e+01,
     & 4.22870e+01, 3.98910e+01, 3.80630e+01, 3.66400e+01, 3.55010e+01,
     & 3.45590e+01, 3.37520e+01, 3.30340e+01, 3.23750e+01, 3.17530e+01,
     & 3.11560e+01, 3.05730e+01, 2.99990e+01, 2.94320e+01, 2.88700e+01,
     & 2.83130e+01, 2.77620e+01, 2.72170e+01, 2.66800e+01, 2.61510e+01,
     & 2.56380e+01, 2.51520e+01, 2.46740e+01, 2.42050e+01, 2.37460e+01,
     & 2.32970e+01, 2.28590e+01, 2.24310e+01, 2.20130e+01, 2.16070e+01,
     & 2.12110e+01, 2.08270e+01, 2.04520e+01, 2.00880e+01, 1.97340e+01,
     & 1.93900e+01, 1.90550e+01, 1.87300e+01, 1.84140e+01, 1.81060e+01,
     & 1.78070e+01, 1.75150e+01, 1.72320e+01, 1.69560e+01, 1.66870e+01,
     & 1.64250e+01, 1.61690e+01, 1.59200e+01, 1.56770e+01, 1.54390e+01,
     & 1.52080e+01, 1.49870e+01, 1.47690e+01, 1.45560e+01, 1.43480e+01,
     & 1.41430e+01, 1.39420e+01, 1.37450e+01, 1.35520e+01, 1.33610e+01,
     & 1.31750e+01, 1.29910e+01, 1.28100e+01, 1.26320e+01, 1.24570e+01,
     & 1.22840e+01, 1.21130e+01, 1.19450e+01, 1.17780e+01, 1.16140e+01,
     & 1.14510e+01, 1.12900e+01, 1.11310e+01, 1.09730e+01, 1.08160e+01,
     & 1.06610e+01, 1.05070e+01, 1.03550e+01, 1.02030e+01, 1.00520e+01,
     & 9.90280e+00, 9.75400e+00, 9.60610e+00, 9.45890e+00, 9.31240e+00,
     & 9.16670e+00, 9.02160e+00, 8.87710e+00, 8.73330e+00, 8.59010e+00,
     & 8.44760e+00, 8.30560e+00, 8.16430e+00, 8.02380e+00, 7.88390e+00,
     & 7.74480e+00, 7.60660e+00, 7.46930e+00, 7.33300e+00, 7.19770e+00,
     & 7.06370e+00, 6.93100e+00, 6.79970e+00, 6.67000e+00, 6.54190e+00,
     & 6.41570e+00, 6.29150e+00, 6.16940e+00, 6.04960e+00, 5.93230e+00,
     & 5.81750e+00, 5.70540e+00, 5.59620e+00, 5.49000e+00, 5.38680e+00,
     & 5.28690e+00, 5.19030e+00, 5.09700e+00, 5.00720e+00, 4.92090e+00,
     & 4.83800e+00, 4.75870e+00, 4.68280e+00, 4.61040e+00, 4.54140e+00,
     & 4.47580e+00, 4.41330e+00, 4.35410e+00, 4.29780e+00, 4.24450e+00,
     & 4.19400e+00, 4.14610e+00, 4.10080e+00, 4.05780e+00, 4.01700e+00,
     & 3.97830e+00, 3.94150e+00, 3.90650e+00, 3.87310e+00, 3.84130e+00,
     & 3.81080e+00, 3.78160e+00, 3.75360e+00, 3.72660e+00, 3.70060e+00,
     & 3.67540e+00, 3.65110e+00, 3.62740e+00, 3.60440e+00, 3.58190e+00,
     & 3.56000e+00, 3.53860e+00, 3.51760e+00, 3.49710e+00, 3.47690e+00,
     & 3.45700e+00, 3.43750e+00, 3.41830e+00, 3.39930e+00, 3.38060e+00,
     & 3.36220e+00, 3.34400e+00, 3.32610e+00, 3.30840e+00, 3.29090e+00,
     & 3.27360e+00, 3.25650e+00, 3.23960e+00, 3.22280e+00, 3.20630e+00,
     & 3.19000e+00, 3.17380e+00, 3.15780e+00, 3.14200e+00, 3.12640e+00,
     & 3.11090e+00, 3.09560e+00, 3.08050e+00, 3.06550e+00, 3.05060e+00,
     & 3.03600e+00, 3.02140e+00, 3.00710e+00, 2.99280e+00, 2.97870e+00,
     & 2.96480e+00, 2.95100e+00, 2.93730e+00, 2.92380e+00, 2.91040e+00,
     & 2.89720e+00, 2.88400e+00, 2.87100e+00, 2.85810e+00, 2.84540e+00,
     & 2.83280e+00, 2.82030e+00, 2.80790e+00, 2.79560e+00, 2.78340e+00,
     & 2.77140e+00, 2.75950e+00, 2.74770e+00, 2.73600e+00, 2.72440e+00,
     & 2.71290e+00, 2.70150e+00, 2.69020e+00, 2.67900e+00, 2.66790e+00,
     & 2.65690e+00, 2.64610e+00, 2.63530e+00, 2.62460e+00, 2.61400e+00,
     & 2.60350e+00, 2.59310e+00, 2.58270e+00, 2.57250e+00, 2.56240e+00,
     & 2.55230e+00, 2.54230e+00, 2.53250e+00, 2.52270e+00, 2.51290e+00,
     & 2.50330e+00, 2.49370e+00, 2.48430e+00, 2.47490e+00, 2.46550e+00,
     & 2.45630e+00, 2.44710e+00, 2.43800e+00, 2.42900e+00, 2.42010e+00,
     & 2.41120e+00, 2.40240e+00, 2.39360e+00, 2.38500e+00, 2.37640e+00,
     & 2.36790e+00, 2.35940e+00, 2.35100e+00, 2.34270e+00, 2.33440e+00,
     & 2.32620e+00, 2.31810e+00, 2.31000e+00, 2.30200e+00, 2.29400e+00,
     & 2.28610e+00, 2.27830e+00, 2.27050e+00, 2.26280e+00, 2.25520e+00,
     & 2.24760e+00, 2.24000e+00, 2.23250e+00, 2.22510e+00, 2.21770e+00,
     & 2.21040e+00, 2.20310e+00, 2.19590e+00, 2.18880e+00, 2.18160e+00,
     & 2.17460e+00, 2.16760e+00, 2.16060e+00, 2.15370e+00, 2.14680e+00,
     & 2.14000e+00, 2.13320e+00, 2.12650e+00, 2.11980e+00, 2.11320e+00,
     & 2.10660e+00, 2.10010e+00, 2.09360e+00, 2.08720e+00, 2.08080e+00,
     & 2.07440e+00, 2.06810e+00, 2.06180e+00, 2.05560e+00, 2.04940e+00,
     & 2.04320e+00, 2.03710e+00, 2.03110e+00, 2.02500e+00, 2.01910e+00/
      DATA (ABSICE4(I,11),I=1,300) /
c     BAND 11
     & 9.13640e+01, 1.17980e+02, 1.11220e+02, 9.78830e+01, 8.46530e+01,
     & 7.37940e+01, 6.49920e+01, 5.80300e+01, 5.26520e+01, 4.84680e+01,
     & 4.51560e+01, 4.24270e+01, 4.03390e+01, 3.87080e+01, 3.73990e+01,
     & 3.63140e+01, 3.53840e+01, 3.45590e+01, 3.38050e+01, 3.31000e+01,
     & 3.24260e+01, 3.17730e+01, 3.11360e+01, 3.05100e+01, 2.98940e+01,
     & 2.92880e+01, 2.86910e+01, 2.81030e+01, 2.75260e+01, 2.69610e+01,
     & 2.64090e+01, 2.58780e+01, 2.53590e+01, 2.48520e+01, 2.43580e+01,
     & 2.38770e+01, 2.34090e+01, 2.29540e+01, 2.25120e+01, 2.20820e+01,
     & 2.16650e+01, 2.12600e+01, 2.08680e+01, 2.04870e+01, 2.01170e+01,
     & 1.97580e+01, 1.94100e+01, 1.90720e+01, 1.87440e+01, 1.84250e+01,
     & 1.81150e+01, 1.78140e+01, 1.75210e+01, 1.72360e+01, 1.69590e+01,
     & 1.66890e+01, 1.64260e+01, 1.61700e+01, 1.59200e+01, 1.56770e+01,
     & 1.54390e+01, 1.52100e+01, 1.49850e+01, 1.47650e+01, 1.45500e+01,
     & 1.43390e+01, 1.41320e+01, 1.39290e+01, 1.37300e+01, 1.35350e+01,
     & 1.33430e+01, 1.31550e+01, 1.29700e+01, 1.27870e+01, 1.26080e+01,
     & 1.24310e+01, 1.22570e+01, 1.20850e+01, 1.19150e+01, 1.17470e+01,
     & 1.15810e+01, 1.14180e+01, 1.12560e+01, 1.10950e+01, 1.09360e+01,
     & 1.07790e+01, 1.06230e+01, 1.04680e+01, 1.03140e+01, 1.01620e+01,
     & 1.00100e+01, 9.85950e+00, 9.70990e+00, 9.56100e+00, 9.41290e+00,
     & 9.26560e+00, 9.11910e+00, 8.97320e+00, 8.82800e+00, 8.68340e+00,
     & 8.53950e+00, 8.39630e+00, 8.25380e+00, 8.11210e+00, 7.97100e+00,
     & 7.83090e+00, 7.69150e+00, 7.55320e+00, 7.41580e+00, 7.27960e+00,
     & 7.14460e+00, 7.01090e+00, 6.87870e+00, 6.74810e+00, 6.61910e+00,
     & 6.49210e+00, 6.36700e+00, 6.24410e+00, 6.12350e+00, 6.00540e+00,
     & 5.88980e+00, 5.77700e+00, 5.66710e+00, 5.56020e+00, 5.45640e+00,
     & 5.35580e+00, 5.25860e+00, 5.16480e+00, 5.07440e+00, 4.98750e+00,
     & 4.90410e+00, 4.82430e+00, 4.74790e+00, 4.67510e+00, 4.60560e+00,
     & 4.53950e+00, 4.47670e+00, 4.41700e+00, 4.36040e+00, 4.30680e+00,
     & 4.25590e+00, 4.20780e+00, 4.16210e+00, 4.11880e+00, 4.07780e+00,
     & 4.03880e+00, 4.00180e+00, 3.96660e+00, 3.93300e+00, 3.90090e+00,
     & 3.87030e+00, 3.84090e+00, 3.81270e+00, 3.78550e+00, 3.75930e+00,
     & 3.73400e+00, 3.70940e+00, 3.68560e+00, 3.66240e+00, 3.63980e+00,
     & 3.61780e+00, 3.59620e+00, 3.57510e+00, 3.55440e+00, 3.53410e+00,
     & 3.51410e+00, 3.49440e+00, 3.47510e+00, 3.45600e+00, 3.43720e+00,
     & 3.41870e+00, 3.40040e+00, 3.38230e+00, 3.36450e+00, 3.34680e+00,
     & 3.32940e+00, 3.31220e+00, 3.29520e+00, 3.27830e+00, 3.26170e+00,
     & 3.24530e+00, 3.22900e+00, 3.21290e+00, 3.19700e+00, 3.18120e+00,
     & 3.16570e+00, 3.15030e+00, 3.13500e+00, 3.11990e+00, 3.10500e+00,
     & 3.09020e+00, 3.07560e+00, 3.06110e+00, 3.04680e+00, 3.03260e+00,
     & 3.01860e+00, 3.00470e+00, 2.99090e+00, 2.97730e+00, 2.96380e+00,
     & 2.95050e+00, 2.93720e+00, 2.92410e+00, 2.91120e+00, 2.89830e+00,
     & 2.88560e+00, 2.87300e+00, 2.86060e+00, 2.84820e+00, 2.83600e+00,
     & 2.82380e+00, 2.81180e+00, 2.79990e+00, 2.78810e+00, 2.77650e+00,
     & 2.76490e+00, 2.75340e+00, 2.74200e+00, 2.73080e+00, 2.71960e+00,
     & 2.70860e+00, 2.69760e+00, 2.68670e+00, 2.67600e+00, 2.66530e+00,
     & 2.65470e+00, 2.64420e+00, 2.63380e+00, 2.62350e+00, 2.61330e+00,
     & 2.60320e+00, 2.59310e+00, 2.58320e+00, 2.57330e+00, 2.56350e+00,
     & 2.55380e+00, 2.54420e+00, 2.53460e+00, 2.52520e+00, 2.51580e+00,
     & 2.50650e+00, 2.49720e+00, 2.48810e+00, 2.47900e+00, 2.47000e+00,
     & 2.46100e+00, 2.45220e+00, 2.44340e+00, 2.43460e+00, 2.42600e+00,
     & 2.41740e+00, 2.40890e+00, 2.40040e+00, 2.39200e+00, 2.38370e+00,
     & 2.37540e+00, 2.36720e+00, 2.35910e+00, 2.35100e+00, 2.34300e+00,
     & 2.33500e+00, 2.32720e+00, 2.31930e+00, 2.31150e+00, 2.30380e+00,
     & 2.29620e+00, 2.28860e+00, 2.28100e+00, 2.27350e+00, 2.26610e+00,
     & 2.25870e+00, 2.25140e+00, 2.24410e+00, 2.23690e+00, 2.22970e+00,
     & 2.22260e+00, 2.21550e+00, 2.20850e+00, 2.20160e+00, 2.19460e+00,
     & 2.18780e+00, 2.18090e+00, 2.17420e+00, 2.16740e+00, 2.16080e+00,
     & 2.15410e+00, 2.14750e+00, 2.14100e+00, 2.13450e+00, 2.12800e+00,
     & 2.12160e+00, 2.11530e+00, 2.10890e+00, 2.10270e+00, 2.09640e+00,
     & 2.09020e+00, 2.08410e+00, 2.07790e+00, 2.07190e+00, 2.06580e+00/
      DATA (ABSICE4(I,12),I=1,300) /
c     BAND 12
     & 3.23110e+01, 4.58690e+01, 4.84580e+01, 4.59140e+01, 4.28200e+01,
     & 3.99020e+01, 3.72630e+01, 3.50130e+01, 3.32060e+01, 3.17950e+01,
     & 3.05270e+01, 2.90510e+01, 2.79490e+01, 2.71140e+01, 2.64640e+01,
     & 2.59410e+01, 2.55020e+01, 2.51180e+01, 2.47650e+01, 2.44300e+01,
     & 2.41040e+01, 2.37790e+01, 2.34530e+01, 2.31240e+01, 2.27910e+01,
     & 2.24550e+01, 2.21160e+01, 2.17760e+01, 2.14360e+01, 2.10970e+01,
     & 2.07760e+01, 2.05050e+01, 2.02330e+01, 1.99600e+01, 1.96880e+01,
     & 1.94170e+01, 1.91490e+01, 1.88820e+01, 1.86190e+01, 1.83600e+01,
     & 1.81040e+01, 1.78510e+01, 1.76030e+01, 1.73590e+01, 1.71190e+01,
     & 1.68840e+01, 1.66530e+01, 1.64250e+01, 1.62030e+01, 1.59840e+01,
     & 1.57690e+01, 1.55580e+01, 1.53520e+01, 1.51490e+01, 1.49490e+01,
     & 1.47540e+01, 1.45620e+01, 1.43730e+01, 1.41870e+01, 1.40050e+01,
     & 1.38260e+01, 1.36520e+01, 1.34800e+01, 1.33110e+01, 1.31440e+01,
     & 1.29790e+01, 1.28170e+01, 1.26570e+01, 1.24990e+01, 1.23430e+01,
     & 1.21890e+01, 1.20360e+01, 1.18850e+01, 1.17360e+01, 1.15880e+01,
     & 1.14420e+01, 1.12970e+01, 1.11530e+01, 1.10100e+01, 1.08690e+01,
     & 1.07280e+01, 1.05890e+01, 1.04500e+01, 1.03120e+01, 1.01750e+01,
     & 1.00380e+01, 9.90190e+00, 9.76640e+00, 9.63150e+00, 9.49690e+00,
     & 9.36280e+00, 9.22900e+00, 9.09550e+00, 8.96230e+00, 8.82940e+00,
     & 8.69670e+00, 8.56420e+00, 8.43190e+00, 8.29990e+00, 8.16810e+00,
     & 8.03650e+00, 7.90520e+00, 7.77420e+00, 7.64350e+00, 7.51320e+00,
     & 7.38340e+00, 7.25410e+00, 7.12530e+00, 6.99730e+00, 6.87000e+00,
     & 6.74370e+00, 6.61840e+00, 6.49420e+00, 6.37130e+00, 6.24980e+00,
     & 6.12990e+00, 6.01170e+00, 5.89540e+00, 5.78110e+00, 5.66900e+00,
     & 5.55930e+00, 5.45200e+00, 5.34740e+00, 5.24550e+00, 5.14650e+00,
     & 5.05060e+00, 4.95770e+00, 4.86800e+00, 4.78150e+00, 4.69830e+00,
     & 4.61850e+00, 4.54200e+00, 4.46880e+00, 4.39890e+00, 4.33230e+00,
     & 4.26880e+00, 4.20850e+00, 4.15120e+00, 4.09690e+00, 4.04530e+00,
     & 3.99640e+00, 3.95010e+00, 3.90620e+00, 3.86460e+00, 3.82510e+00,
     & 3.78760e+00, 3.75200e+00, 3.71810e+00, 3.68580e+00, 3.65490e+00,
     & 3.62540e+00, 3.59710e+00, 3.56990e+00, 3.54380e+00, 3.51850e+00,
     & 3.49410e+00, 3.47050e+00, 3.44760e+00, 3.42520e+00, 3.40350e+00,
     & 3.38220e+00, 3.36140e+00, 3.34110e+00, 3.32110e+00, 3.30150e+00,
     & 3.28230e+00, 3.26330e+00, 3.24460e+00, 3.22620e+00, 3.20810e+00,
     & 3.19020e+00, 3.17260e+00, 3.15510e+00, 3.13790e+00, 3.12090e+00,
     & 3.10410e+00, 3.08740e+00, 3.07100e+00, 3.05480e+00, 3.03870e+00,
     & 3.02280e+00, 3.00710e+00, 2.99160e+00, 2.97620e+00, 2.96100e+00,
     & 2.94590e+00, 2.93100e+00, 2.91630e+00, 2.90170e+00, 2.88730e+00,
     & 2.87300e+00, 2.85890e+00, 2.84490e+00, 2.83100e+00, 2.81730e+00,
     & 2.80370e+00, 2.79030e+00, 2.77700e+00, 2.76380e+00, 2.75080e+00,
     & 2.73790e+00, 2.72510e+00, 2.71240e+00, 2.69990e+00, 2.68750e+00,
     & 2.67520e+00, 2.66300e+00, 2.65090e+00, 2.63900e+00, 2.62720e+00,
     & 2.61540e+00, 2.60380e+00, 2.59230e+00, 2.58090e+00, 2.56960e+00,
     & 2.55840e+00, 2.54730e+00, 2.53630e+00, 2.52540e+00, 2.51460e+00,
     & 2.50390e+00, 2.49330e+00, 2.48280e+00, 2.47240e+00, 2.46210e+00,
     & 2.45180e+00, 2.44170e+00, 2.43160e+00, 2.42160e+00, 2.41180e+00,
     & 2.40200e+00, 2.39220e+00, 2.38260e+00, 2.37310e+00, 2.36360e+00,
     & 2.35420e+00, 2.34490e+00, 2.33560e+00, 2.32650e+00, 2.31740e+00,
     & 2.30840e+00, 2.29950e+00, 2.29060e+00, 2.28180e+00, 2.27310e+00,
     & 2.26440e+00, 2.25580e+00, 2.24730e+00, 2.23890e+00, 2.23050e+00,
     & 2.22220e+00, 2.21400e+00, 2.20580e+00, 2.19770e+00, 2.18960e+00,
     & 2.18160e+00, 2.17370e+00, 2.16580e+00, 2.15800e+00, 2.15030e+00,
     & 2.14260e+00, 2.13490e+00, 2.12740e+00, 2.11980e+00, 2.11240e+00,
     & 2.10500e+00, 2.09760e+00, 2.09030e+00, 2.08310e+00, 2.07590e+00,
     & 2.06880e+00, 2.06170e+00, 2.05470e+00, 2.04770e+00, 2.04070e+00,
     & 2.03390e+00, 2.02700e+00, 2.02020e+00, 2.01350e+00, 2.00680e+00,
     & 2.00020e+00, 1.99360e+00, 1.98700e+00, 1.98050e+00, 1.97410e+00,
     & 1.96770e+00, 1.96130e+00, 1.95500e+00, 1.94870e+00, 1.94250e+00,
     & 1.93630e+00, 1.93010e+00, 1.92400e+00, 1.91800e+00, 1.91190e+00,
     & 1.90590e+00, 1.90000e+00, 1.89410e+00, 1.88820e+00, 1.88240e+00/
      DATA (ABSICE4(I,13),I=1,300) /
c     BAND 13
     & 6.51160e+01, 8.71150e+01, 8.51680e+01, 7.59760e+01, 6.76210e+01,
     & 6.04480e+01, 5.44080e+01, 4.94970e+01, 4.56400e+01, 4.26150e+01,
     & 4.01370e+01, 3.78900e+01, 3.61840e+01, 3.48610e+01, 3.38090e+01,
     & 3.29430e+01, 3.22040e+01, 3.15510e+01, 3.09530e+01, 3.03900e+01,
     & 2.98480e+01, 2.93190e+01, 2.87990e+01, 2.82830e+01, 2.77710e+01,
     & 2.72630e+01, 2.67590e+01, 2.62590e+01, 2.57660e+01, 2.52790e+01,
     & 2.48090e+01, 2.43680e+01, 2.39340e+01, 2.35070e+01, 2.30870e+01,
     & 2.26760e+01, 2.22730e+01, 2.18790e+01, 2.14940e+01, 2.11180e+01,
     & 2.07510e+01, 2.03930e+01, 2.00450e+01, 1.97050e+01, 1.93740e+01,
     & 1.90510e+01, 1.87360e+01, 1.84300e+01, 1.81310e+01, 1.78400e+01,
     & 1.75570e+01, 1.72800e+01, 1.70110e+01, 1.67480e+01, 1.64910e+01,
     & 1.62400e+01, 1.59950e+01, 1.57560e+01, 1.55230e+01, 1.52940e+01,
     & 1.50710e+01, 1.48540e+01, 1.46410e+01, 1.44320e+01, 1.42280e+01,
     & 1.40270e+01, 1.38300e+01, 1.36370e+01, 1.34460e+01, 1.32590e+01,
     & 1.30760e+01, 1.28950e+01, 1.27170e+01, 1.25410e+01, 1.23680e+01,
     & 1.21970e+01, 1.20290e+01, 1.18620e+01, 1.16980e+01, 1.15350e+01,
     & 1.13740e+01, 1.12150e+01, 1.10570e+01, 1.09010e+01, 1.07460e+01,
     & 1.05930e+01, 1.04400e+01, 1.02890e+01, 1.01390e+01, 9.98970e+00,
     & 9.84140e+00, 9.69390e+00, 9.54710e+00, 9.40110e+00, 9.25580e+00,
     & 9.11110e+00, 8.96710e+00, 8.82360e+00, 8.68080e+00, 8.53860e+00,
     & 8.39690e+00, 8.25590e+00, 8.11550e+00, 7.97570e+00, 7.83670e+00,
     & 7.69840e+00, 7.56090e+00, 7.42430e+00, 7.28870e+00, 7.15420e+00,
     & 7.02080e+00, 6.88880e+00, 6.75810e+00, 6.62900e+00, 6.50160e+00,
     & 6.37590e+00, 6.25230e+00, 6.13070e+00, 6.01140e+00, 5.89450e+00,
     & 5.78020e+00, 5.66860e+00, 5.55980e+00, 5.45400e+00, 5.35130e+00,
     & 5.25180e+00, 5.15550e+00, 5.06260e+00, 4.97320e+00, 4.88720e+00,
     & 4.80460e+00, 4.72560e+00, 4.65000e+00, 4.57790e+00, 4.50910e+00,
     & 4.44370e+00, 4.38150e+00, 4.32240e+00, 4.26640e+00, 4.21330e+00,
     & 4.16300e+00, 4.11530e+00, 4.07010e+00, 4.02720e+00, 3.98660e+00,
     & 3.94800e+00, 3.91140e+00, 3.87650e+00, 3.84330e+00, 3.81150e+00,
     & 3.78120e+00, 3.75210e+00, 3.72420e+00, 3.69730e+00, 3.67130e+00,
     & 3.64630e+00, 3.62200e+00, 3.59840e+00, 3.57550e+00, 3.55310e+00,
     & 3.53130e+00, 3.50990e+00, 3.48900e+00, 3.46850e+00, 3.44840e+00,
     & 3.42860e+00, 3.40920e+00, 3.39000e+00, 3.37120e+00, 3.35250e+00,
     & 3.33420e+00, 3.31610e+00, 3.29820e+00, 3.28050e+00, 3.26310e+00,
     & 3.24580e+00, 3.22880e+00, 3.21200e+00, 3.19530e+00, 3.17880e+00,
     & 3.16260e+00, 3.14650e+00, 3.13050e+00, 3.11480e+00, 3.09920e+00,
     & 3.08380e+00, 3.06850e+00, 3.05340e+00, 3.03850e+00, 3.02370e+00,
     & 3.00910e+00, 2.99460e+00, 2.98030e+00, 2.96610e+00, 2.95210e+00,
     & 2.93820e+00, 2.92440e+00, 2.91080e+00, 2.89730e+00, 2.88400e+00,
     & 2.87080e+00, 2.85770e+00, 2.84470e+00, 2.83190e+00, 2.81920e+00,
     & 2.80660e+00, 2.79420e+00, 2.78180e+00, 2.76960e+00, 2.75750e+00,
     & 2.74550e+00, 2.73360e+00, 2.72180e+00, 2.71020e+00, 2.69860e+00,
     & 2.68710e+00, 2.67580e+00, 2.66460e+00, 2.65340e+00, 2.64240e+00,
     & 2.63140e+00, 2.62060e+00, 2.60980e+00, 2.59920e+00, 2.58860e+00,
     & 2.57820e+00, 2.56780e+00, 2.55750e+00, 2.54730e+00, 2.53720e+00,
     & 2.52720e+00, 2.51730e+00, 2.50740e+00, 2.49760e+00, 2.48800e+00,
     & 2.47840e+00, 2.46880e+00, 2.45940e+00, 2.45000e+00, 2.44070e+00,
     & 2.43150e+00, 2.42240e+00, 2.41330e+00, 2.40430e+00, 2.39540e+00,
     & 2.38660e+00, 2.37780e+00, 2.36910e+00, 2.36050e+00, 2.35190e+00,
     & 2.34340e+00, 2.33500e+00, 2.32660e+00, 2.31830e+00, 2.31010e+00,
     & 2.30190e+00, 2.29380e+00, 2.28580e+00, 2.27780e+00, 2.26990e+00,
     & 2.26200e+00, 2.25420e+00, 2.24650e+00, 2.23880e+00, 2.23120e+00,
     & 2.22360e+00, 2.21610e+00, 2.20860e+00, 2.20120e+00, 2.19390e+00,
     & 2.18660e+00, 2.17930e+00, 2.17220e+00, 2.16500e+00, 2.15790e+00,
     & 2.15090e+00, 2.14390e+00, 2.13700e+00, 2.13010e+00, 2.12330e+00,
     & 2.11650e+00, 2.10970e+00, 2.10300e+00, 2.09640e+00, 2.08980e+00,
     & 2.08320e+00, 2.07670e+00, 2.07030e+00, 2.06380e+00, 2.05750e+00,
     & 2.05110e+00, 2.04480e+00, 2.03860e+00, 2.03240e+00, 2.02620e+00,
     & 2.02010e+00, 2.01400e+00, 2.00800e+00, 2.00200e+00, 1.99600e+00/
      DATA (ABSICE4(I,14),I=1,300) /
c     BAND 14
     & 4.38370e+01, 5.87350e+01, 5.83400e+01, 5.33720e+01, 4.86930e+01,
     & 4.45320e+01, 4.09360e+01, 3.79670e+01, 3.56280e+01, 3.38110e+01,
     & 3.22560e+01, 3.06200e+01, 2.93930e+01, 2.84570e+01, 2.77250e+01,
     & 2.71330e+01, 2.66340e+01, 2.61960e+01, 2.57960e+01, 2.54180e+01,
     & 2.50510e+01, 2.46890e+01, 2.43270e+01, 2.39640e+01, 2.36000e+01,
     & 2.32340e+01, 2.28670e+01, 2.25000e+01, 2.21340e+01, 2.17700e+01,
     & 2.14240e+01, 2.11220e+01, 2.08210e+01, 2.05210e+01, 2.02230e+01,
     & 1.99280e+01, 1.96360e+01, 1.93480e+01, 1.90650e+01, 1.87850e+01,
     & 1.85110e+01, 1.82410e+01, 1.79760e+01, 1.77160e+01, 1.74610e+01,
     & 1.72110e+01, 1.69660e+01, 1.67260e+01, 1.64900e+01, 1.62600e+01,
     & 1.60340e+01, 1.58120e+01, 1.55950e+01, 1.53830e+01, 1.51740e+01,
     & 1.49700e+01, 1.47690e+01, 1.45720e+01, 1.43790e+01, 1.41890e+01,
     & 1.40030e+01, 1.38210e+01, 1.36420e+01, 1.34660e+01, 1.32930e+01,
     & 1.31220e+01, 1.29540e+01, 1.27880e+01, 1.26240e+01, 1.24630e+01,
     & 1.23040e+01, 1.21470e+01, 1.19920e+01, 1.18380e+01, 1.16860e+01,
     & 1.15360e+01, 1.13880e+01, 1.12410e+01, 1.10950e+01, 1.09500e+01,
     & 1.08070e+01, 1.06640e+01, 1.05230e+01, 1.03820e+01, 1.02430e+01,
     & 1.01040e+01, 9.96570e+00, 9.82820e+00, 9.69140e+00, 9.55500e+00,
     & 9.41920e+00, 9.28380e+00, 9.14890e+00, 9.01430e+00, 8.88000e+00,
     & 8.74600e+00, 8.61240e+00, 8.47900e+00, 8.34590e+00, 8.21310e+00,
     & 8.08060e+00, 7.94840e+00, 7.81660e+00, 7.68510e+00, 7.55400e+00,
     & 7.42350e+00, 7.29350e+00, 7.16420e+00, 7.03550e+00, 6.90770e+00,
     & 6.78080e+00, 6.65500e+00, 6.53040e+00, 6.40700e+00, 6.28510e+00,
     & 6.16480e+00, 6.04620e+00, 5.92950e+00, 5.81490e+00, 5.70240e+00,
     & 5.59240e+00, 5.48480e+00, 5.37990e+00, 5.27770e+00, 5.17850e+00,
     & 5.08220e+00, 4.98910e+00, 4.89910e+00, 4.81240e+00, 4.72910e+00,
     & 4.64900e+00, 4.57230e+00, 4.49890e+00, 4.42880e+00, 4.36200e+00,
     & 4.29840e+00, 4.23790e+00, 4.18050e+00, 4.12600e+00, 4.07430e+00,
     & 4.02520e+00, 3.97880e+00, 3.93480e+00, 3.89300e+00, 3.85350e+00,
     & 3.81590e+00, 3.78010e+00, 3.74610e+00, 3.71370e+00, 3.68280e+00,
     & 3.65320e+00, 3.62480e+00, 3.59760e+00, 3.57130e+00, 3.54600e+00,
     & 3.52160e+00, 3.49790e+00, 3.47490e+00, 3.45250e+00, 3.43060e+00,
     & 3.40930e+00, 3.38850e+00, 3.36810e+00, 3.34810e+00, 3.32840e+00,
     & 3.30910e+00, 3.29010e+00, 3.27140e+00, 3.25290e+00, 3.23470e+00,
     & 3.21680e+00, 3.19910e+00, 3.18160e+00, 3.16430e+00, 3.14720e+00,
     & 3.13040e+00, 3.11370e+00, 3.09720e+00, 3.08090e+00, 3.06480e+00,
     & 3.04890e+00, 3.03310e+00, 3.01750e+00, 3.00210e+00, 2.98690e+00,
     & 2.97180e+00, 2.95680e+00, 2.94210e+00, 2.92740e+00, 2.91300e+00,
     & 2.89860e+00, 2.88450e+00, 2.87040e+00, 2.85650e+00, 2.84280e+00,
     & 2.82920e+00, 2.81570e+00, 2.80240e+00, 2.78920e+00, 2.77610e+00,
     & 2.76310e+00, 2.75030e+00, 2.73760e+00, 2.72500e+00, 2.71260e+00,
     & 2.70020e+00, 2.68800e+00, 2.67590e+00, 2.66390e+00, 2.65210e+00,
     & 2.64030e+00, 2.62860e+00, 2.61710e+00, 2.60570e+00, 2.59430e+00,
     & 2.58310e+00, 2.57200e+00, 2.56090e+00, 2.55000e+00, 2.53920e+00,
     & 2.52840e+00, 2.51780e+00, 2.50730e+00, 2.49680e+00, 2.48650e+00,
     & 2.47620e+00, 2.46600e+00, 2.45590e+00, 2.44590e+00, 2.43600e+00,
     & 2.42620e+00, 2.41640e+00, 2.40680e+00, 2.39720e+00, 2.38770e+00,
     & 2.37830e+00, 2.36890e+00, 2.35960e+00, 2.35050e+00, 2.34130e+00,
     & 2.33230e+00, 2.32330e+00, 2.31440e+00, 2.30560e+00, 2.29690e+00,
     & 2.28820e+00, 2.27960e+00, 2.27100e+00, 2.26260e+00, 2.25420e+00,
     & 2.24580e+00, 2.23760e+00, 2.22930e+00, 2.22120e+00, 2.21310e+00,
     & 2.20510e+00, 2.19710e+00, 2.18930e+00, 2.18140e+00, 2.17360e+00,
     & 2.16590e+00, 2.15830e+00, 2.15070e+00, 2.14310e+00, 2.13560e+00,
     & 2.12820e+00, 2.12080e+00, 2.11350e+00, 2.10620e+00, 2.09900e+00,
     & 2.09190e+00, 2.08480e+00, 2.07770e+00, 2.07070e+00, 2.06370e+00,
     & 2.05680e+00, 2.05000e+00, 2.04320e+00, 2.03640e+00, 2.02970e+00,
     & 2.02300e+00, 2.01640e+00, 2.00980e+00, 2.00330e+00, 1.99680e+00,
     & 1.99040e+00, 1.98400e+00, 1.97770e+00, 1.97140e+00, 1.96510e+00,
     & 1.95890e+00, 1.95270e+00, 1.94660e+00, 1.94050e+00, 1.93440e+00,
     & 1.92840e+00, 1.92250e+00, 1.91650e+00, 1.91060e+00, 1.90480e+00/
      DATA (ABSICE4(I,15),I=1,300) /
c     BAND 15
     & 5.19930e+02, 2.73810e+02, 1.67330e+02, 1.19150e+02, 9.23010e+01,
     & 7.53010e+01, 6.36670e+01, 5.54570e+01, 4.95940e+01, 4.52800e+01,
     & 4.20290e+01, 3.94440e+01, 3.74820e+01, 3.59610e+01, 3.47510e+01,
     & 3.37570e+01, 3.29130e+01, 3.21700e+01, 3.14960e+01, 3.08670e+01,
     & 3.02700e+01, 2.96920e+01, 2.91280e+01, 2.85750e+01, 2.80300e+01,
     & 2.74930e+01, 2.69630e+01, 2.64410e+01, 2.59280e+01, 2.54230e+01,
     & 2.49340e+01, 2.44700e+01, 2.40160e+01, 2.35710e+01, 2.31360e+01,
     & 2.27110e+01, 2.22970e+01, 2.18940e+01, 2.15000e+01, 2.11180e+01,
     & 2.07450e+01, 2.03820e+01, 2.00300e+01, 1.96860e+01, 1.93530e+01,
     & 1.90280e+01, 1.87120e+01, 1.84050e+01, 1.81060e+01, 1.78140e+01,
     & 1.75310e+01, 1.72550e+01, 1.69860e+01, 1.67230e+01, 1.64670e+01,
     & 1.62180e+01, 1.59740e+01, 1.57370e+01, 1.55040e+01, 1.52770e+01,
     & 1.50550e+01, 1.48390e+01, 1.46260e+01, 1.44190e+01, 1.42150e+01,
     & 1.40150e+01, 1.38190e+01, 1.36270e+01, 1.34380e+01, 1.32530e+01,
     & 1.30700e+01, 1.28910e+01, 1.27140e+01, 1.25400e+01, 1.23680e+01,
     & 1.21990e+01, 1.20320e+01, 1.18670e+01, 1.17050e+01, 1.15440e+01,
     & 1.13840e+01, 1.12270e+01, 1.10710e+01, 1.09160e+01, 1.07630e+01,
     & 1.06110e+01, 1.04600e+01, 1.03110e+01, 1.01620e+01, 1.00140e+01,
     & 9.86750e+00, 9.72150e+00, 9.57620e+00, 9.43160e+00, 9.28770e+00,
     & 9.14440e+00, 9.00160e+00, 8.85950e+00, 8.71790e+00, 8.57690e+00,
     & 8.43640e+00, 8.29640e+00, 8.15710e+00, 8.01840e+00, 7.88030e+00,
     & 7.74290e+00, 7.60640e+00, 7.47060e+00, 7.33580e+00, 7.20200e+00,
     & 7.06940e+00, 6.93800e+00, 6.80790e+00, 6.67940e+00, 6.55240e+00,
     & 6.42730e+00, 6.30400e+00, 6.18290e+00, 6.06390e+00, 5.94730e+00,
     & 5.83330e+00, 5.72190e+00, 5.61330e+00, 5.50770e+00, 5.40510e+00,
     & 5.30570e+00, 5.20950e+00, 5.11670e+00, 5.02730e+00, 4.94130e+00,
     & 4.85870e+00, 4.77970e+00, 4.70410e+00, 4.63190e+00, 4.56310e+00,
     & 4.49760e+00, 4.43540e+00, 4.37630e+00, 4.32020e+00, 4.26700e+00,
     & 4.21660e+00, 4.16880e+00, 4.12360e+00, 4.08060e+00, 4.03990e+00,
     & 4.00130e+00, 3.96460e+00, 3.92960e+00, 3.89630e+00, 3.86450e+00,
     & 3.83410e+00, 3.80500e+00, 3.77700e+00, 3.75000e+00, 3.72400e+00,
     & 3.69890e+00, 3.67450e+00, 3.65090e+00, 3.62790e+00, 3.60550e+00,
     & 3.58360e+00, 3.56220e+00, 3.54120e+00, 3.52070e+00, 3.50050e+00,
     & 3.48070e+00, 3.46110e+00, 3.44190e+00, 3.42300e+00, 3.40430e+00,
     & 3.38590e+00, 3.36770e+00, 3.34980e+00, 3.33210e+00, 3.31460e+00,
     & 3.29730e+00, 3.28020e+00, 3.26330e+00, 3.24660e+00, 3.23000e+00,
     & 3.21370e+00, 3.19750e+00, 3.18160e+00, 3.16570e+00, 3.15010e+00,
     & 3.13460e+00, 3.11930e+00, 3.10420e+00, 3.08920e+00, 3.07430e+00,
     & 3.05960e+00, 3.04510e+00, 3.03070e+00, 3.01650e+00, 3.00240e+00,
     & 2.98850e+00, 2.97460e+00, 2.96100e+00, 2.94740e+00, 2.93400e+00,
     & 2.92080e+00, 2.90760e+00, 2.89460e+00, 2.88170e+00, 2.86900e+00,
     & 2.85630e+00, 2.84380e+00, 2.83140e+00, 2.81910e+00, 2.80700e+00,
     & 2.79490e+00, 2.78300e+00, 2.77110e+00, 2.75940e+00, 2.74780e+00,
     & 2.73630e+00, 2.72490e+00, 2.71360e+00, 2.70240e+00, 2.69130e+00,
     & 2.68030e+00, 2.66940e+00, 2.65860e+00, 2.64790e+00, 2.63730e+00,
     & 2.62680e+00, 2.61630e+00, 2.60600e+00, 2.59570e+00, 2.58560e+00,
     & 2.57550e+00, 2.56550e+00, 2.55560e+00, 2.54580e+00, 2.53610e+00,
     & 2.52640e+00, 2.51680e+00, 2.50730e+00, 2.49790e+00, 2.48860e+00,
     & 2.47930e+00, 2.47010e+00, 2.46100e+00, 2.45200e+00, 2.44300e+00,
     & 2.43410e+00, 2.42530e+00, 2.41650e+00, 2.40790e+00, 2.39930e+00,
     & 2.39070e+00, 2.38220e+00, 2.37380e+00, 2.36550e+00, 2.35720e+00,
     & 2.34900e+00, 2.34080e+00, 2.33270e+00, 2.32470e+00, 2.31670e+00,
     & 2.30880e+00, 2.30100e+00, 2.29320e+00, 2.28540e+00, 2.27780e+00,
     & 2.27010e+00, 2.26260e+00, 2.25510e+00, 2.24760e+00, 2.24020e+00,
     & 2.23290e+00, 2.22560e+00, 2.21840e+00, 2.21120e+00, 2.20400e+00,
     & 2.19700e+00, 2.18990e+00, 2.18300e+00, 2.17600e+00, 2.16910e+00,
     & 2.16230e+00, 2.15550e+00, 2.14880e+00, 2.14210e+00, 2.13540e+00,
     & 2.12880e+00, 2.12230e+00, 2.11580e+00, 2.10930e+00, 2.10290e+00,
     & 2.09650e+00, 2.09020e+00, 2.08390e+00, 2.07760e+00, 2.07140e+00,
     & 2.06530e+00, 2.05910e+00, 2.05310e+00, 2.04700e+00, 2.04100e+00/
      END
