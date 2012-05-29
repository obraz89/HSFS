      SUBROUTINE LHADY
C
C This is key solver of 3-D stability equations in 3D compressible boundary layer
C
C           MAB=0 -ASYMPTOTIC IS NOT ESTABLISHED; DIRECT PROBLEM;
C           MAB=1 -HOMOGENEOUS ASYMPTOTICS; DIRECT PROBLEM;
C           MAB=2   "-------------------"   CONJUGATE PROBLEM;
C           MAB=4 -ASYMPTOTICS ARE NOT ESTABLISHED "-----------"
C           MAB=5 -COMPUTATIONS TO UP DIRECTION
C           MAB=6 -ONLY ASYMPTOTICS ARE COMPUTED FI,LAM
C     COMPUTATION TO UP DIRECTION: MAB=3
C
C Developed by Fedrorov A.V., 12.11.1985
C Do not distribute without Fedorov permission
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 FU(8,201),FI(8,8),Q(8,4,201),CF(4),AS(4,4),SP,Z,E,A,B
     *   ,W,R,WA,XI,LAM(4),BI(4,4),S1,S2,VT(8,4),H(8,8),SI,AZ(4,8,4)
     *   ,F0(8),F1(8),F2(8,4),VMS(8,8),FIQ(8,4),ZZZ0
      EQUIVALENCE (FI(1,1),FIQ(1,1))
      REAL*8 M,MU,MU1,MU2,MG,MY,MYY,K,M2,MM,M1,KN
      DIMENSION IQ(201)
      LOGICAL LOG
C
      COMMON /HADY1/ Y1,D,NSB,L1
     *       /HADY2/ A,B,W,R
     *       /HADY3/ SI,ZZZ0
     *       /MACKY1/ Y1O,DDDO,C4,YCC,NO  !DDDO
     *       /HADY4/ FI,LAM
     *       /HADY5/ FU,CF
     *       /HADY6/ Q,BI,VT,AS
     *       /ABOCT/ MAB,MAB1,MMAB2
     *       /OSNOB/ Y,U,U1,U2,T,T1,T2,MU,MU1,MU2,WS,WS1,WS2
     *       /OUT  / K,PR,G,M,XXXI,PPPB
     
      COMMON /CF_PROFILE_DATA/ MF_I, MF_J, MF_X, MF_WE, MF_ANG
      integer MF_I, MF_J
      real*8 MF_X, MF_WE, MF_ANG
C
      DATA N,MF,L/201,8,4/
      N=NO
C      
      NS=NSB
      IF(NS.LT.0) NS=-NS
      ICON=0
      L10=L1
      IF(MAB.EQ.2.OR.MAB.EQ.4.OR.MAB.EQ.6) ICON=1
      IF(L1.LE.-1) WRITE(6,10)A,B,W,R,Y1,NS,MAB,K,PR,G,M,LAM,FIQ
10    FORMAT(' :LHADY :'/' A=',2(1PD10.3),' B=',
     *       2(1PD10.3),'W=',2(1PD10.3)/
     *     ' R=',2(1PD10.3),' Y1=',(1PD10.3),' NS=',I3,
     *     ' MAB=',I3,' K=',(1PD10.3)/,' PR='
     *       ,(1PD10.3),' G=',(1PD10.3),' M=',1PD10.3/
     *       ' LAM='/4(2X,2D9.2)/' FI=',8(/4(2X,2D9.2)))
C
C     Verification of parameters
C
      LOG=CDABS(A).LT.40..AND.CDABS(B).LT.40..AND.CDABS(W).LT.40.
     *    .AND.CDABS(A-W).GT.1.D-5.AND.Y1.GT.0..AND.NS.GE.1
     *    .AND.MAB.GE.0.AND.MAB.LE.6.AND.K.GE.0..AND.K.LT.1.D6
     *    .AND.PR.GT.0..AND.PR.LT.40..AND.G.GE.1..AND.G.LT.40.
     *    .AND.M.GE.0..AND.M.LT.100..AND.CDABS(R).GT.1.D-5
     *    .AND.CDABS(R).LT.1.E8
      IF(LOG) GO TO 1
      PRINT 20,LOG
20    FORMAT(' :LHADY **ERROR IN PARAMETERS** LOG=',L1)
97    WRITE(6,10)A,B,W,R,Y1,NS,MAB,K,PR,G,M,LAM,FIQ
      !to stop mess
      STOP
      MAB=-1
      RETURN
1     CONTINUE
      IF(MAB.EQ.3) GO TO 301
C
C     Additional constants
C
      D=Y1/(N-1)
      D1=-D/NS
      Y=Y1
      IF(MAB.NE.5) GO TO 46
      Y=0.
      D1=-D1
46    NFI=2*MF*L
      E=(0.,1.)
      DO 2 I=1,MF
      DO 2 J=1,MF
2     H(I,J)=0.
      MM=2./3.*(K-1.)
      F=MM+1.
      RM=2./3.*(K+2.)
      M2=G*M*M
      MG=(G-1.)*M*M
      ASSIGN 3 TO METC
      IF(MAB.EQ.1.OR.MAB.EQ.2.OR.MAB.EQ.6) GO TO 201
3     CONTINUE
C
C     Computation from the upper boundary down to the wall
C
      L1=0
      DO 4 I=1,L
      S=AMAXMU(FI(1,I),16)
      IF(S.GT.1.D-4.AND.S.LT.1.E9) GO TO 4
      PRINT 30,I,S
30    FORMAT(' :LHADY **ERROR IN ASYMPTOTICS ',I1,' S=',(1PD10.3))
      GO TO 97
4     CONTINUE
      ASSIGN 5 TO MET1
      GO TO 101
5     CONTINUE
      IND=0
6     IND=IND+1
      LNS=0
8     LNS=LNS+1
      IF(LNS.NE.1.AND.NSB.GT.0) GO TO 7
      IF(LNS.EQ.1) CALL UCOPY(FI,Q(1,1,N-IND+1),NFI)
      S=AMAXMU(FI,NFI)
      IF(S.LT.1000.)GO TO 7
      L1=L1+1
      IF(NSB.GT.0) IQ(L1)=N-IND+1
      IF(L10.LE.-2) WRITE(6,51)FIQ
51    FORMAT(' *** ORTOGONALIZATION *** FI=',4(/2(1P,8D10.3/)))
      CALL ORT(FI,MF,L,1.D-6)
      IF(MAB.GE.0) GO TO 7
      WRITE(6,40)IND,Y,FIQ
40    FORMAT(' :LHADY **ERROR AT ORTOG.** IND=',I3,
     *  ' Y=',F8.2,' FI=',4(/2(1P,8D10.3/)))
      GO TO 97
7     CONTINUE
      IF(L10.LE.-2) WRITE(6,52)IND,Y,U,U1,T,MU,FIQ
52    FORMAT(' LHADY IND=',I2,' Y=',F5.2,' U,U1,T,MU=',4F5.2,
     *       ' VECTORS P/K:',4(/2(1P,8D10.3/)))
      IF(IND.EQ.N) GO TO 28
C
C Start Runge-Kutta integration
C
      LN=0
      ASSIGN 11 TO MET2
9     LN=LN+1
      CALL UCOPY(FI(1,LN),F0,16)
      GO TO 103
11    CONTINUE
      DO 12 I=1,MF
      AZ(1,I,LN)=D1*F1(I)
12    F2(I,LN)=F0(I)+0.5*AZ(1,I,LN)
      IF(LN-4) 9,13,13
13    CONTINUE
      Y=Y+D1*0.5
      ASSIGN 14 TO MET1
      GO TO 101
14    CONTINUE
      LN=0
15    LN=LN+1
      ASSIGN 16 TO MET2
      CALL UCOPY(F2(1,LN),F0,16)
      GO TO 103
16    CONTINUE
      DO 17 I=1,MF
      AZ(2,I,LN)=D1*F1(I)
17    F2(I,LN)=FI(I,LN)+AZ(2,I,LN)*0.5
      CALL UCOPY(F2(1,LN),F0,16)
      ASSIGN 18 TO MET2
      GO TO 103
18    CONTINUE
      DO 19 I=1,MF
      AZ(3,I,LN)=D1*F1(I)
19    F2(I,LN)=FI(I,LN)+AZ(3,I,LN)
      IF(LN-4) 15,21,21
21    CONTINUE
      Y=Y+D1*0.5
      ASSIGN 22 TO MET1
      GO TO 101
22    CONTINUE
       LN=0
       ASSIGN 24 TO MET2
23     LN=LN+1
       CALL UCOPY(F2(1,LN),F0,16)
       GO TO 103
24     CONTINUE
       DO 25 I=1,MF
       AZ(4,I,LN)=D1*F1(I)
       FI(I,LN)=FI(I,LN)+(AZ(1,I,LN)+2.*(AZ(2,I,LN)+AZ(3,I,LN))
     *        +AZ(4,I,LN))/6.
25    CONTINUE
       IF(LN-4) 23,26,26
C    Runge-Kutta is done
26     CONTINUE
       IF(LNS-NS) 8,27,27
27     CONTINUE
      GO TO 6
28    CONTINUE
C
C     Calculation of functional
C
      IF(MAB.EQ.5) RETURN
      IF(ICON.EQ.0) GO TO 41
      DO 42 I=1,4
      F0(I)=0.
      AS(I,1)=FI(3,I)
      AS(I,2)=FI(2,I)
      AS(I,3)=FI(6,I)
      AS(I,4)=FI(8,I)
42    CONTINUE
      F0(1)=1.
      GO TO 43
41    DO 29 I=1,4
      F0(I)=0.
      AS(I,1)=FI(1,I)
      AS(I,2)=FI(2,I)
      AS(I,3)=FI(3,I)
29    AS(I,4)=FI(7,I)
      JB=2
      IF(MAB.EQ.1) JB=2
      F0(JB)=1.
43    CALL SYS(AS,F0,CF,4)
      IF(MAB.GE.0) GO TO 31
      PRINT 50,AS
50    FORMAT(' :LHADY **AS=',/,4(4(2X,1PD10.3,1PD10.3),/))
      GO TO 97
31    IF(ICON.EQ.0) GO TO 44
C     Conjugate problem
      SI=FI(4,1)*CF(1)+FI(4,2)*CF(2)+FI(4,3)*CF(3)+FI(4,4)*CF(4)
      ART=CDABS(SI)
      IF(ART.LT.1.E 4) GO TO 301
      PRINT 90,ART
90    FORMAT(' LHADY CONJUGATE PROBL. ***NEVJASKA=',(1PD10.3))
      MAB=-7
      GO TO 97
44    SI=FI(5,1)*CF(1)+FI(5,2)*CF(2)+FI(5,3)*CF(3)+FI(5,4)*CF(4)
      ART=CDABS(SI)
C     IF(MAB1.GE.1)  WRITE(6,80)ART,W,L1,MAB
80    FORMAT(' :SI=',1PE9.2,' W=',2D15.7,' L1=',I4,' MAB=',I4)
      RETURN
C
C     Matrix of stability equations
C
101   CALL OSN(N)
102   WA=A*U+B*WS-W
      M1=1./MU
      MY=MU1*T1
      MYY=MU2*T1*T1+MU1*T2
      XI=1./(R*M1+E*RM*M2*WA)
      TD=1./T
      UM=MU2*T1*U1+MU1*U2
C     1 -ST LINE:
      H(2,1)=1.
C     2-ND LINE:
      H(1,2)=E*R*TD*M1*WA+A*A+B*B
      H(2,2)=-M1*MY
      H(3,2)=R*TD*M1*U1-E*A*(M1*MY+F*TD*T1)
      H(4,2)=E*A*R*M1-F*M2*A*WA
      H(5,2)=F*A*TD*WA-M1*UM
      H(6,2)=-M1*MU1*U1
C     3-D LINE:
      H(1,3)=-E*A
      H(3,3)=TD*T1
      H(4,3)=-E*M2*WA
      H(5,3)=E*TD*WA
      H(7,3)=-E*B
C     4 -TH LINE:
      H(1,4)=-E*XI*A*(2.*M1*MY+RM*TD*T1)
      H(2,4)=-E*XI*A
      H(3,4)=XI*(-A*A-B*B+RM*TD*M1*MY*T1+RM*TD*T2-E*R*TD*M1*WA)
      H(4,4)=-E*XI*RM*M2*((M1*MY+TD*T1)*WA+A*U1+B*WS1)
      H(5,4)=E*XI*((M1*MU1+RM*TD)*(A*U1+B*WS1)+RM*TD*M1*MY*WA)
      H(6,4)=E*XI*RM*TD*WA
      H(7,4)=-E*XI*B*(2.*M1*MY+RM*TD*T1)
      H(8,4)=-E*XI*B
C     5 -TH LINE:
      H(6,5)=1.
C     6 -TH LINE:
      H(2,6)=-2.*PR*MG*U1
      H(3,6)=R*PR*TD*M1*T1-2.*E*PR*MG*(A*U1+B*WS1)
      H(4,6)=-E*R*PR*M1*MG*WA
      H(5,6)=E*R*PR*TD*M1*WA+A*A+B*B-MG*PR*M1*MU1*(U1*U1+WS1**2)-M1*MYY
      H(6,6)=-2.*M1*MY
      H(8,6)=-2.*PR*MG*WS1
C     7 -TH LINE:
      H(8,7)=1.
C     8 -TH LINE:
      H(3,8)=-E*B*(M1*MY+F*TD*T1)+R*M1*TD*WS1
      H(4,8)=E*R*B*M1-F*B*M2*WA
      H(5,8)=F*B*TD*WA-M1*(MU2*T1*WS1+MU1*WS2)
      H(6,8)=-M1*MU1*WS1
      H(7,8)=E*R*TD*M1*WA+A*A+B*B
      H(8,8)=-M1*MY
      GO TO MET1,(202,5,14,22)
103   CONTINUE
C.................................
C     MULTIPLICATION: F1(*)=H*F0
C................................
      IF(ICON.EQ.0) GO TO 106
C     CONJUGATE PROBLEM:
      DO 107 LT1=1,8
      F1(LT1)=0.
      DO 107 LT2=1,8
107   F1(LT1)=F1(LT1)-H(LT1,LT2)*F0(LT2)
      GO TO 105
106   CONTINUE
      DO 104 LT1=1,8
      F1(LT1)=0.
      DO 104 LT2=1,8
      F1(LT1)=F1(LT1)+H(LT2,LT1)*F0(LT2)
104   CONTINUE
105   GO TO MET2,(11,16,18,24)
201   CONTINUE
C
C     ASYMPTOTICS OUT OF THE LAYER
C
      U=1.
      U1=0.
      U2=0.
      !WS=0.
      WS = MF_WE
      WS1=0.
      WS2=0.
      T=1.
      T1=0.
      T2=0.
      MU=1.
      MU1=0.
      MU2=0.
*      Y=Y1
*      CALL OSN
      ASSIGN 202 TO MET1
      GO TO 102
202   CONTINUE
      BI(1,1)=H(1,2)
      BI(2,1)=H(4,2)
      BI(3,1)=H(5,2)
      BI(2,2)=H(4,2)*H(2,4)+H(4,3)*H(3,4)+H(4,6)*H(6,4)+H(4,8)*H(8,4)
      BI(3,2)=H(5,2)*H(2,4)+H(5,3)*H(3,4)+H(6,4)*H(5,6)+H(8,4)*H(5,8)
      BI(2,3)=H(4,6)
      BI(3,3)=H(5,6)
      BI(2,4)=H(4,8)
      BI(3,4)=H(5,8)
      BI(4,4)=H(1,2)
      LAM(1)=-CDSQRT(BI(1,1))
      S1=(BI(2,2)+BI(3,3))*0.5
      S2=CDSQRT((BI(2,2)-BI(3,3))**2*0.25+BI(3,2)*BI(2,3))
      LAM(2)=-CDSQRT(S1+S2)
      LAM(3)=-CDSQRT(S1-S2)
      LAM(4)=LAM(1)
      IF(ICON.EQ.0) GO TO 208
C     ASYMPTOTICS FOR CONJUGATE PROBLEM.
      ASSIGN 209 TO METC
      GO TO 208
209   CONTINUE
      DO 211 I=1,4
211   LAM(I)=-LAM(I)
      CALL UCOPY(FI,VMS,64)
      ASSIGN 213 TO METC
      GO TO 208
213   CONTINUE
      CALL UCOPY(FI,VMS(1,5),64)
      IF(MAB.NE.6) GO TO 217
      DO 218 I=1,4
218   LAM(I)=-LAM(I)
      CALL UCOPY(VMS,FI,128)
      RETURN
217   CALL UCOPY(VMS,AZ,128)
      DO 214 I=1,4
      DO 215 J=1,8
215   F0(J)=0.
      F0(I+4)=1.
      CALL SYS(VMS,F0,FI(1,I),8)
      CALL UCOPY(AZ,VMS,128)
214   CONTINUE
      DO 216 I=1,4
216   LAM(I)=-LAM(I)
      ASSIGN 3 TO METC
      DO 212 J=1,4
      DO 212 I=1,8
      VT(I,J)=-LAM(J)*FI(I,J)
      DO 212 KI=1,8
      VT(I,J)=VT(I,J)-H(I,KI)*FI(KI,J)
212   CONTINUE
C      GO TO 207
      GO TO 206
208   BI(1,1)=1.
      BI(2,1)=0.
      BI(3,1)=0.
      BI(4,1)=0.
      DO 203 J=2,3
      BI(1,J)=((LAM(J)**2-H(5,6))*H(4,2)+H(5,2)*H(4,6))/(H(1,2)-
     -LAM(J)**2)
      BI(2,J)=H(5,6)-LAM(J)**2
      BI(3,J)=-H(4,6)
      BI(4,J)=(H(4,6)*H(5,8)+(LAM(J)**2-H(5,6))*H(4,8))/(H(1,2)-
     -LAM(J)**2)
203   CONTINUE
      BI(1,4)=0.
      BI(2,4)=0.
      BI(3,4)=0.
      BI(4,4)=1.
      DO 204 J=1,4
      FI(1,J)=BI(1,J)
      FI(2,J)=LAM(J)*BI(1,J)
      FI(3,J)=(H(1,3)*BI(1,J)+H(4,3)*BI(2,J)+H(5,3)*BI(3,J)+H(7,3)*
     *BI(4,J))/LAM(J)
      FI(4,J)=BI(2,J)
      FI(5,J)=BI(3,J)
      FI(6,J)=LAM(J)*BI(3,J)
      FI(7,J)=BI(4,J)
      FI(8,J)=(H(4,8)*BI(2,J)+H(5,8)*BI(3,J)+H(7,8)*BI(4,J))/LAM(J)
204   CONTINUE
C
C  VERIFICATION OF ASYMPTOTICS. VT=H*FI(J)-LAM(J)*FI(J)
C                         J=1,.....,4.
C
      DO 205 J=1,4
      DO 205 I=1,8
      VT(I,J)=-LAM(J)*FI(I,J)
      DO 205 KI=1,8
      VT(I,J)=VT(I,J)+H(KI,I)*FI(KI,J)
205   CONTINUE
207   VTR=AMAXMU(VT,64)
      IF(VTR.LT.1.D-3) GO TO 206
      WRITE(6,210)VTR
210    FORMAT(' :LHADY ERROR IN ASYMPTOTICS; NEVJASKA',(1PD10.3))
      PRINT 220,VT,H
220   FORMAT(' VT=',/,8(4(2X,2E11.3),/),
     *       ' H=',/,16(4(2X,2E11.3),/))
      GO TO 97
206    CONTINUE
      GO TO METC,(3,209,213)
301   CONTINUE
C
C     RECONSTRUCTION OF SOLUTION:
C
      DO 399 IND=1,N
      IF(L1.EQ.0)GOTO 304
      IF(IND.NE.IQ(L1)) GO TO 304   ! fix incomplete drawing 
      CALL UCOPY(Q(1,1,IND),FI,NFI)
      CALL ORT(FI,MF,L,1.D-6)
C     PRINT 376
376   FORMAT(' **ORTOGONALISATION** ')
      IF(MAB.GE.0) GO TO 302
      PRINT 310,IND,IQ,((Q(I,J,IND),I=1,MF),J=1,L)
310    FORMAT('  :BBEPX*** ERROR AT ORT.',/,6(20I4,/)
     *        ,6(12(1PD10.3),/))
      GO TO 97
302    CONTINUE
      L1=L1-1
      DO 303 I=1,L
      FU(I,IND)=CF(I)
      DO 303 J=1,L
      AS(J,I)=SP(FI(1,I),Q(1,J,IND),MF)
303     CONTINUE
      CALL SYS(AS,FU(1,IND),CF,L)
304     CONTINUE
      DO 306 J=1,MF
      Z=0.
      DO 305 I=1,L
305   Z=Z+CF(I)*Q(J,I,IND)
306   FU(J,IND)=Z
C     PRINT 375,IND,(FU(I,IND),I=1,8)

375   FORMAT(' IND=',I4,' FU=',2(/,8D10.2))
399   CONTINUE
      CALL UCOPY(Q(1,1,201),FI,NFI)
      IND=1
C     PRINT 375,IND,(FU(I,IND),I=1,8)
      IF(L1.EQ.0) RETURN
C      PRINT 320,L1,FI,IQ
320   FORMAT(' ** UP: CALL DOWN, L1,IQ=',
     */,20(I5))
      GO TO 97
      END

      FUNCTION LAG(X0,X,F,N)                                            
C********************************************************************** 
C     POLINOM LAGRANG
C********************************************************************** 
      IMPLICIT REAL*8 (A-H,O-Z)                                         
      COMMON/ABOCT/ MAB,MAB1,MAB2
      REAL*8 LAG,X0,X(N),F(N),P,S                                   
      P=0.                                                              
      DO 2 I=1,N                                                        
      S=1.                                                              
      DO 1 J=1,N                                                        
      IF(J.EQ.I) GO TO 1                                                
      IF(DABS(X(I)-X(J)).GT.1.D-10) GO TO 3                            
      PRINT 4,X                                                         
4     FORMAT(' :LAG **POINTS COINSIDE *** X=',88(/1P5D14.7))  
      MAB=-1                                                            
      RETURN                                                            
3     S=S*(X0-X(J))/(X(I)-X(J))                                         
1     CONTINUE                                                          
2     P=P+F(I)*S                                                        
      LAG=P                                                             
      RETURN                                                            
      END                                                               

      SUBROUTINE HADSI
C**********************************************************************
C     Service of LHADY to calculate functional
C     NF=number of line taking out, JF=number of vector with C=1.
C**********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /HADY3/ SI,Z0
     *       /HADY2/ A,B,W,R
     *       /HADY4/ FI(8,8),LAM(4)
     *       /ABOCT/ MAB,MAB1,MAB2
     *       /SHOCK1/ARB(7),LS
      DIMENSION RS(8)
      COMPLEX*16 SI,FI,A,B,W,R,Z0,LAM
      L1=0  
      CALL LHADY
      RETURN
      END

      SUBROUTINE VGR
C
C   COMPUTATION OF GROOP VELOCITY (VA,VB)=(DW/DA,DW/DB) 
C                   AND RATIO VR=VB/VA
C        
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /HADY2/ A,B,W,R
     *       /HADY3/ SI,ZZZ0
     *       /ABOCT/ MAB,MAB1,MAB2
     *       /VGRC/ VA,VB,VR
      COMPLEX*16 SI,A,B,W,R,ZZZ0,A0,B0,W0,W1,VA,VB,VR
      DATA D /1.D-6/
      A0=A
      B0=B 
      W0=W
      A=A0+D
      CALL POISK2(3)
      W1=W
      A=A0-D
      W=W0
      CALL POISK2(3)
      VA=(W1-W)*0.5D0/D
      A=A0
      B=B0+D
      W=W0
      CALL POISK2(3)
      W1=W
      B=B0-D
      W=W0
      CALL POISK2(3)
      VB=(W1-W)*0.5D0/D
      A=A0
      B=B0
      W=W0
      VR=VB/VA
      EPS=DIMAG(VR)/DREAL(VR)
      EPS=DABS(EPS)
      IF(EPS.GT.1.D-3)THEN
!      WRITE(*,10)
!   10 FORMAT('VGR : VB/VA is complex!')
!      STOP
      ENDIF 
      RETURN
      END

      SUBROUTINE AMAX(EPS)
C
C Calculattion of maximum growth rate -IMAG(A) versus frequency W
C    
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /HADY2/ A,B,W,R
     *       /HADY3/ SI,ZZZ0
     *       /ABOCT/ MAB,MAB1,MAB2
      COMPLEX*16 SI,A,B,W,R,ZZZ0,A00
      MAB1=0
      W0=DREAL(W)
      CALL POISK2(1)
    1 SS0=-DIMAG(A)
      A00=A
      W0=DREAL(W)
      DW=W0*0.001
      W=W0+DW
      A=A00*(W/W0)
      CALL POISK2(1)
      SS1=-DIMAG(A)
      W=W0-DW
      A=A00*(W/W0)
      CALL POISK2(1)
      SS2=-DIMAG(A)
      DELW=(SS1-SS2)*0.5*DW/(SS1-SS0*2.+SS2)
      W=W0-DELW
      A=A00*(W/W0)
      CALL POISK2(1)
      SSS=-DIMAG(A)
      SSS=SSS-SS0
	WRITE(*,1000)W,A,SSS
      IF(SSS.GT.0.D0.AND.SSS.LT.EPS)GO TO 2
      GO TO 1
    2 CONTINUE
 1000 FORMAT(1P,100E13.5)
      RETURN
      END

!
      SUBROUTINE GR_MAX_LSP(point, b_range)
C
C Calculattion of maximum growth rate -IMAG(A)+gv_dir*IMAG(B) versus frequency W
C    
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /HADY2/ A,B,W,R
     *       /HADY3/ SI,ZZZ0
     *       /ABOCT/ MAB,MAB1,MAB2
     *       /LSP_OUT/ GR_DIR
      COMPLEX*16 SI,A,B,W,R,ZZZ0,A00
      real*8 step
      complex*16 point,b_range, GR_DIR
      print*, "-------------------------------------"
      point = B
      b_range = DCMPLX(DABS(DREAL(B)), DABS(DIMAG(B)))
      MAB1=0
      W0=DREAL(W)
      CALL LSP_MAX_REC(point,b_range)
    1 SS0=-GR_DIR
      A00=A
      W0=DREAL(W)
      DW=W0*0.001
      W=W0+DW
      A=A00*(W/W0)
      point = B
      b_range = DCMPLX(DABS(DREAL(B)), DABS(DIMAG(B)))
      CALL LSP_MAX_REC(point,b_range)
      SS1=-GR_DIR
      W=W0-DW
      A=A00*(W/W0)
      point = B
      b_range = DCMPLX(DABS(DREAL(B)), DABS(DIMAG(B)))
      CALL LSP_MAX_REC(point,b_range)
      SS2=-GR_DIR
      DELW=(SS1-SS2)*0.5*DW/(SS1-SS0*2.+SS2)
      W=W0-DELW
      A=A00*(W/W0)
      point = B
      b_range = DCMPLX(DABS(DREAL(B)), DABS(DIMAG(B)))
      CALL LSP_MAX_REC(point,b_range)
      SSS=-GR_DIR
      SSS=SSS-SS0
	WRITE(*,1000)W,A,SSS
      IF(SSS.GT.0.D0.AND.SSS.LT.1.0d-5)GO TO 2
      GO TO 1
    2 CONTINUE
 1000 FORMAT(1P,100E13.5)
      RETURN
      END
