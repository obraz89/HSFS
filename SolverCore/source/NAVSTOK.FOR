      SUBROUTINE NAVSTOK(NNS,XST)
C
C This subroutine adapt Navier-Stokes profiles, which are in
C common-block /NS/, for stability calculations
C
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NY=541)
      COMMON/NS/YYY(NY),U(NY),U1(NY),U2(NY),T(NY),T1(NY),T2(NY),
     *W(NY),W1(NY),W2(NY),ENT(NY),FI(NY),DFI(NY)
	COMMON/RHONS/ RHO(NY)
     *      /OSNOB/Y,BU,BU1,BU2,BT,BT1,BT2,BMU,BMU1,BMU2,BW,BW1,BW2
     *      /OSNV/BBV,BBV1,VBB,VBB1
     *      /IUPT/IUPT
	COMMON/DNS/AMINF,REINF,TINF,XLL,REE,REE1,UE
      COMMON/BASA3/PRA,T0,A3,A4,A5,M,BETA,CQ,XI,AMK,TETB,Y1,TW,G
     *       /BASA6/ A6(10),BV,DL1,TWWW
      COMMON/OUT/K,SIG,GG,MM,XXI,PB
     *      /MACK02/TELOC,TWALL,TFFF
     *      /MACKY1/ Y1O,DDDD,C,YC,NO
     *      /MACKNJ/NJOB,IPRINT,LOLD,IV
     *      /MACKXY/XFIRST,YFIRST,ZFIRST
     *      /HADY1/Y1H,DOD,NNNS,LLL1
     *      /MOMENT/DL2
      REAL*8 M,K,MM
      DIMENSION F(NY),V(NY),V1(NY)
c
c Determine the shock locus
c
	YS=YYY(NNS)
	JS=NNS
	J=0
  200 J=J+1
	anev=dabs(1.d0-T(j))
	IF(ANEV.LT.1.e-4)GOTO 201
	IF(J.LT.NNS)GOTO 200
	YS=YYY(J)
	JS=J
  201 CONTINUE
c
C Search of upper boundary Ye
C
C Specify ratio of U'(Y_e)/U'(0)
	EPS=1.E-2
C
	J=0
  100 J=J+1
	DERIF=DABS(U1(J)/U1(1))
	IF(DERIF.LT.EPS)GOTO 101
	IF(J.LT.NNS)GOTO 100
  101 YE=YYY(J)*2.D0
	IF(YE.GE.YS)THEN
	J=JS-4
	YE=YYY(J)
	GOTO 102
	ENDIF
	J=0
  103 J=J+1
	IF(YYY(J).LT.YE)GOTO 103
  102 CONTINUE
C
C Parameters at the upper boundary Ye
	N=J
	UE=U(N)
	TTE=T(N)
	RHOE=RHO(N)
	AME=AMINF*UE/DSQRT(TTE)
	MM=AME
	M=AME
      T0=TINF*(1.+AMINF**2*(GG-1.)*0.5)
C Viscosity
      IF(IUPT.NE.0)GO TO 224
      EMU=TTE**BV
      GO TO 223
224   CONTINUE
      X=TTE
      E=(GG-1.)*MM**2
      TE=E*0.5
      XX=110.4D0/T0
      XX=XX*(1.D0+TE)
      EMU=X**1.5*(1.D0+XX)/(X+XX)
223   CONTINUE
C
	REE=REINF*UE*RHOE*XST/EMU ! Local Reynolds number at station XST
	REE1=REINF*UE*RHOE/EMU/XLL ! Local unit Reynolds number
	REEL=REINF*UE*RHOE/EMU ! Local Reynolds number based on length XLL
      M=AME ! Local Mach number
	MM=AME ! Local Mach number
	TELOC=TTE*TINF ! temperature at the upper b.l. edge
	TWALL=T(1)/TTE ! Wall temperature Twall/Te
      BETA=PB
      IK=0
      MAB=0
      CQ=0.
      T0=TELOC*(1.+MM*MM*0.5*(GG-1.)) ! Stagnation temperature
      E=(GG-1.)*MM**2
      TE=E*0.5
      IP=0
      IW=0
      IF(NJOB.EQ.2) IW=1
      TETB=1.000
      XI=XXI
      AMK=M
      PRA=SIG
      TW=TWALL
      TWWW=T(1)
C
c Normalize y-coordinate to SQRT(NU_e*X/U_e)
c Normalize flow parameters to their values at Ye=YE
c   		
      AAA=DSQRT(REINF/REEL)
c
      DO 37 I=1,N
      YYY(I)=YYY(I)/AAA
	U(I)=U(I)/U(N)
      U1(I)=U1(I)*AAA/U(N)
      U2(I)=U2(I)*AAA**2/U(N)
	T(I)=T(I)/T(N)
      T1(I)=T1(I)*AAA/T(N)
      T2(I)=T2(I)*AAA**2/T(N)
      W1(I)=W1(I)*AAA
37    W2(I)=W2(I)*AAA**2
C
C Calculate the boundary-layer displacement thicknesses
      DO 1 I=1,N
1     F(I)=1.D0-U(I)/T(I)
      DL1=SIMPR(F,YYY,N)
      DO 2 I=1,N
2     F(I)=U(I)*(1.D0-U(I))/T(I)
      DL2=SIMPR(F,YYY,N)
      DO 3 I=1,N
3     F(I)=U(I)*(T(I)-1.D0)/T(I)
      DLH=SIMPR(F,YYY,N)
C
C Interpolate to uniform grid
      Y1=YYY(N)
	Y1O=Y1
	NN=NO
	Y1H=Y1
      D=Y1/(NN-1.D0)
      DO 4 I=1,NN
      Y=D*(I-1.D0)
4     F(I)=PARAB2(Y,YYY,U,N)
      DO 5 I=1,NN
5     U(I)=F(I)
      DO 6 I=1,NN
      Y=D*(I-1.D0)
6     F(I)=PARAB2(Y,YYY,U1,N)
      DO 7 I=1,NN
7     U1(I)=F(I)
      DO 8 I=1,NN
      Y=D*(I-1.D0)
8     F(I)=PARAB2(Y,YYY,U2,N)
      DO 9 I=1,NN
9     U2(I)=F(I)
      DO 10 I=1,NN
      Y=D*(I-1.D0)
10    F(I)=PARAB2(Y,YYY,T,N)
      DO 11 I=1,NN
11    T(I)=F(I)
      DO 12 I=1,NN
      Y=D*(I-1.D0)
12    F(I)=PARAB2(Y,YYY,T1,N)
      DO 13 I=1,N
13    T1(I)=F(I)
      DO 14 I=1,NN
      Y=D*(I-1.D0)
14    F(I)=PARAB2(Y,YYY,T2,N)
      DO 15 I=1,NN
15    T2(I)=F(I)
      DO 16 I=1,NN
      Y=D*(I-1.D0)
16    F(I)=PARAB2(Y,YYY,W,N)
      DO 17 I=1,NN
17    W(I)=F(I)
      DO 18 I=1,NN
      Y=D*(I-1.D0)
18    F(I)=PARAB2(Y,YYY,W1,N)
      DO 19 I=1,NN
19    W1(I)=F(I)
      DO 20 I=1,NN
      Y=D*(I-1.D0)
20    F(I)=PARAB2(Y,YYY,W2,N)
      DO 21 I=1,NN
21    W2(I)=F(I)
      DO 23 I=1,NN
      IF(IUPT.NE.0)GO TO 24
      ENT(I)=T(I)**BV
      FI(I)=BV*T(I)**(BV-1.D0)
      DFI(I)=BV*(BV-1.D0)*T(I)**(BV-2.D0)
      GO TO 23
24    CONTINUE
      X=T(I)
      XX=110.4/T0
      XX=XX*(1.+TE)
      ENT(I)=X**1.5*(1.+XX)/(X+XX)
      FI(I)=1.5*ENT(I)/X-ENT(I)/(X+XX)
      DFI(I)=1.5*(FI(I)-ENT(I)/X)/X-
     *(FI(I)-ENT(I)/(X+XX))/(X+XX)
23    CONTINUE
      F(1)=0.
      DO 336 I=1,NN-1
336   F(I+1)=F(I)+0.5*D*(U(I)/T(I)+U(I+1)/T(I+1))
      DO 337 I=1,NN
      Y=D*(I-1.D0)
      V(I)=0.5*(Y*U(I)-T(I)*F(I))
337   V1(I)=0.5*(Y*U1(I)-T1(I)*F(I))
	N=NN
      RETURN
c
C Entry to calculate mean flow parameters at Y
c
      ENTRY OSN
      BU=PARAB1(Y,DOD,U,NO)
      BU1=PARAB1(Y,DOD,U1,NO)
      BU2=PARAB1(Y,DOD,U2,NO)
      BW=PARAB1(Y,DOD,W,NO)
      BW1=PARAB1(Y,DOD,W1,NO)
      BW2=PARAB1(Y,DOD,W2,NO)
      BT=PARAB1(Y,DOD,T,NO)
      BT1=PARAB1(Y,DOD,T1,NO)
      BT2=PARAB1(Y,DOD,T2,NO)
      BMU=PARAB1(Y,DOD,ENT,NO)
      BMU1=PARAB1(Y,DOD,FI,NO)
      BMU2=PARAB1(Y,DOD,DFI,NO)
      BBV=PARAB1(Y,DOD,V,NO)
      BBV1=PARAB1(Y,DOD,V1,NO)
      IF(Y.LT.Y1) GO TO 28
      BU = 1.0
      BU1=0.
      BU2=0.
      BW = 0.0 !     -0.02761 !!!!!!!!!!!!!!!!!!!!!!!!!!!1
      BW1 = 0.0
      BW2 = 0.0
      BT=1.
      BT1=0.
      BT2=0.
      BMU=1.
c      VBB=V(N)
c      VBB1=V1(N)
28    CONTINUE
      RETURN
      END
