      SUBROUTINE ORT(FI,M,L,EPS)
C
C   Orthonormalization of complex vector FI of dimension M
C   EPS=minimal norma of the vector
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /ABOCT/ MAB
      COMPLEX*16 FI(M,L),SP,P
      CALL UTEST(FI,M*L*2)
      DO 2 I=1,L
      S=SN(FI(1,I),M)
      IF(S.GT.EPS) GO TO 1
      MAB=-1
      PRINT 10,I,M,L
10    FORMAT(' **ORT BEKTOP HOMEP',I3,' M,L=',2I3)
      RETURN
1     CONTINUE
      DO 2 J=1,M
2     FI(J,I)=FI(J,I)/S
      L1=L-1
      IF(L.EQ.1) RETURN
      DO 5 K=1,L1
      K1=K+1
      DO 5 I=K1,L
      P=-SP(FI(1,K),FI(1,I),M)
      DO 3 J=1,M
3     FI(J,I)=P*FI(J,K)+FI(J,I)
      S=SN(FI(1,I),M)
      IF(S.GT.EPS) GO TO 4
      MAB=-2
      PRINT 20,K,M,L
20    FORMAT(' ***ORT KOMPLANARN. K,M,L=',3I3)
      RETURN
4     CONTINUE
      DO 5 J=1,M
5     FI(J,I)=FI(J,I)/S
      RETURN
      END
