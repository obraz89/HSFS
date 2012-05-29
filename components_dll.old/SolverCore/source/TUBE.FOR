      SUBROUTINE TUBE(Z,GB,GK,GRHO,GC)
C
C This subroutine calculates the functions GB(Z) and GK(Z), which are defined
C in JASA, 92, No. 2, 1995 for Pr=0.708, Gamma=1.4
C     
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 Z,GB,GK,GRHO,GC
      EPS=1.E-8
      IF(CDABS(Z).LT.EPS)THEN
      GB=0.
      GK=1.
      GC=0.
      GRHO=1.
      RETURN
      ENDIF
C      
      GAM=1.4
      PR=0.708
C      
C Coefficints
C
      A1=0.
      A2=0.487816
      A3=0.176809
      A4=0.034778
C
      B1=0.487816
      B2=0.272785
      B3=0.081597
      B4=0.014549
C
      W1=DSQRT(GAM/8.)
C
      GB=W1*Z*(1.+A2*Z+A3*Z**2+A4*Z**3)
     */(1.+B1*Z+B2*Z**2+B3*Z**3+B4*Z**4)
C
      A1=0.571125
      A2=0.281620
      A3=0.080889
      A4=0.012106
C
      B1=0.
      B2=0.571125
      B3=0.210930
      B4=0.040516
C
      W1=DSQRT(GAM*8.)
C
      GK=W1*(1.+A1*Z+A2*Z**2+A3*Z**3+A4*Z**4)
     */(1.+B2*Z+B3*Z**2+B4*Z**3)/Z
C
      GRHO=GK/GB
      GC=GK*GB
C      
      RETURN
      END                                               
