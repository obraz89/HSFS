
      SUBROUTINE POISK2(LP)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /HADY2/ Z(4)
     *       /HADY3/ SI,ZZZ0
     *       /ABOCT/ MAB,MAB1,MAB2
     *       /BASA6/ A6(10),BVB,DL1
      COMMON/OUT  /AAK,AASIG,AAG,AAM,AAXI,AAPB
	COMMON/JP/JPOISK
      COMPLEX*16 SI,Z0,F,F1,Z00,ZB,Z,ZZZ0
      DATA EPS /1.D-8/,D/1.D-7/
C------------ Global search for Mack second mode -------
	IF(JPOISK.EQ.1)THEN
	Z(2)=(0.,0.)	! wavenumber beta
	WLOW=0.1	! initial freuency
	WUP=0.6		! final frequncy
	AI=-5.E-3	! growth rate
	CFUP=1.		! maximal phase speed
	CFLOW=1.-1./AAM	! minimal phase speed
	NW=101
	NAR=101
	DW=(WUP-WLOW)/(NW-1)
	FMIN=1.E+9
	DO JW=1,NW
	W=WLOW+DW*(JW-1)
	ARLOW=W/CFUP
	ARUP=W/CFLOW
	DAR=(ARUP-ARLOW)/(NAR-1)
	WRITE(2,1002)
	DO JAR=2,NAR-1
	AR=ARLOW+DAR*(JAR-1)
	Z(3)=W
	Z(1)=DCMPLX(AR,AI)
      CALL HADSI
C	SI - complex, SI (минимум модуля - глобально нашли моду)
	FSI=CDABS(SI)
	IF(FSI.LT.FMIN)THEN
	FMIN=FSI
	ARMIN=AR
	AIMIN=AI
	WMIN=W
	ENDIF
	ENDDO
	ENDDO
	WRITE(*,1000)ARMIN,AIMIN,WMIN,WMIN/ARMIN,FMIN
	Z(1)=DCMPLX(ARMIN,AIMIN)
	Z(3)=WMIN
	ENDIF
C-------------------------------------------------------end of global search
 1000 FORMAT(1P,20E14.6)
 1002 FORMAT(/)
C
C Local search of eigenvalues
C SI=0 <->  собственные значения
C
      K=1
      CALL HADSI
      MAB2=1
      IF(MAB1.EQ.1) WRITE(*,1)Z(LP),SI,K
1     FORMAT(' Z=',2F15.11,' SI=',1P,2D10.3,' K=',I3)
      IF(MAB.LT.0) RETURN
      IF(CDABS(SI).LT.EPS) RETURN
      Z00=Z(LP)
      Z0=Z00
      F=SI
      Z(LP)=Z(LP)+D
      DO 2 K=2,24
      MAB2=K
      CALL HADSI
      IF(MAB1.EQ.1) WRITE(*,1)Z(LP),SI,K
      IF(MAB.LT.0) RETURN
      IF(CDABS(SI).LT.EPS) RETURN
      IF(CDABS(Z(LP)-Z0).LT.1.D-9) RETURN
      F1=(SI-F)/(Z(LP)-Z0)
      F1=F1/F
      IF(CDABS(F1).LT.1.D-3) THEN
         Z(LP)=Z00
         WRITE(6,'('' :POISK2 *** TOO GREAT VARIATION ***'')')
         MAB=-1
         RETURN
      ENDIF
      ZB=Z0
      Z0=Z(LP)
      F=SI
      Z(LP)=ZB-1.D0/F1
2     CONTINUE
      WRITE(6,'('' :POISK ***24 ITERATION ARE USED***'')')
      MAB=-1
      RETURN
      END
