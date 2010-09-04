
      SUBROUTINE SYSR(A,B,X,N)
C
C     Solves the system of linear equations with real coefficients
C
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 A(N,N),B(N),X(N),T
      COMMON /ABOCT/MAB
      DATA EPS/1.D-10/
      N1=N-1
      DO 1 I=1,N
      C=0.
      K=0
      DO 2 J=I,N
      IF(DABS(A(I,J)).LE.C) GO TO 2
      K=J
      C=DABS(A(I,J))
2     CONTINUE
      IF(C.GT.EPS) GO TO 11
      WRITE(5,10)N
10     FORMAT('  ***SYS*** ERROR.-BAD MATRIX N=',I3)
      MAB=-1
      RETURN
11    IF(K.EQ.I)GO TO 4
      DO 3 J=I,N
      T=A(J,I)
      A(J,I)=A(J,K)
3     A(J,K)=T
      T=B(I)
      B(I)=B(K)
      B(K)=T
4     CONTINUE
      T=A(I,I)
      DO 5 J=I,N
5     A(J,I)=A(J,I)/T
      B(I)=B(I)/T
      IF(I.EQ.N) GO TO 1
      I1=I+1
      DO 6 J=I1,N
      T=A(I,J)
      DO 7 K=I,N
7     A(K,J)=A(K,J)-T*A(K,I)
6     B(J)=B(J)-B(I)*T
1     CONTINUE
      X(N)=B(N)
      DO 8 I=1,N1
      X(N-I)=B(N-I)
      NI=N-I+1
      DO 8 J=NI,N
8     X(N-I)=X(N-I)-A(J,N-I)*X(J)
      RETURN
      END
