
      FUNCTION PARAB1(X,D,F,N)
C**********************************************************************
C  parabolic interpolation with constant step
C    output parameters:
C    X   -coordinat
C    D   -step of integration
C    F   -massiv of function
C    N   -dimension of massiv.
C**********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 F(1)
      K=X/D+1
      IF(K.LT.1) K=1
      IF(K.LT.N-1) GO TO 1
      PARAB1=F(N)
      GO TO 2
1     X1=(K-1)*D
      X2=X1+D
      X3=X2+D
      PARAB1=(F(K)*(X2-X)*(X3-X)-2.*F(K+1)*(X1-X)*(X3-X)+F(K+2)*
     *(X2-X)*(X1-X))/(2.*D*D)
2     CONTINUE
      RETURN
      END

      FUNCTION PARAB2(X,Y,F,N)
C**********************************************************************
C  parabolic interpolation with varios step.
C      output parameters:
C  X   -point of interpolation
C  Y   -massiv of coordinates
C  F   -massiv of function;  N   -dimension of massiv
C**********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8   Y(N),F(N)
      IF(X.LE.Y(1)) THEN
        PARAB2=F(1)
        GO TO 3
      ENDIF
      DO 1 I=1,N-2
      K=I
      IF(Y(I).LE.X.AND.Y(I+1).GE.X) GO TO 2
1     CONTINUE
      PARAB2=F(N)
      GO TO 3
2     CONTINUE
      X1=Y(K)
      X2=Y(K+1)
      X3=Y(K+2)
      PARAB2=F(K)*(X-X2)*(X-X3)/((X1-X2)*(X1-X3))+F(K+1)*
     *(X-X1)*(X-X3)/((X2-X1)*(X2-X3))+F(K+2)*(X-X1)*(X-X2)/
     /((X3-X1)*(X3-X2))
3     CONTINUE
      RETURN
      END


      FUNCTION PARAB22(X,Y,F,N)
C**********************************************************************
C  parabolic interpolation with varios step.
C      output parameters:
C  X   -point of interpolation
C  Y   -massiv of coordinates
C  F   -massiv of function;  N   -dimension of massiv
C**********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8   Y(N),F(N)
      IF(X.GE.Y(1)) THEN
        PARAB2=F(1)
        GO TO 3
      ENDIF
      DO 1 I=1,N-2
      K=I
      IF(Y(I).GE.X.AND.Y(I+1).LE.X) GO TO 2
1     CONTINUE
      PARAB2=F(N)
      GO TO 3
2     CONTINUE
      X1=Y(K)
      X2=Y(K+1)
      X3=Y(K+2)
      PARAB22=F(K)*(X-X2)*(X-X3)/((X1-X2)*(X1-X3))+F(K+1)*
     *(X-X1)*(X-X3)/((X2-X1)*(X2-X3))+F(K+2)*(X-X1)*(X-X2)/
     /((X3-X1)*(X3-X2))
3     CONTINUE
      RETURN
      END
