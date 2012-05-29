
      SUBROUTINE SRKR(F,Y,H,M,OP)
      IMPLICIT REAL*8 (A-H,O-Z)
C***********************************************************
C   RUNGE-KUTTA METHOD OF 4 ORDER.
C   IT MAKES ONE STEP H IN INTEGRATING THE SYSTEM:
C     DF/DY=OP(F,Y);  F=VECTOR OF DIMENSION M.
C  OP(F,F1,Y) - SUBROUTINE WHICH FINDS RIGHT HAND SIDE
C  F1=DF/DY.
C  SUB. OP MUST BE DESCRIBED AS EXTERNAL.
C  OUTPUT : F(Y+H).
C  M MUST BE LESS THEN 9.
C***********************************************************
      REAL*8 F(M),A(4,9),F1(9),F2(9)
      CALL OP(F,F1,Y)
      DO 1 I=1,M
      A(1,I)=H*F1(I)
1     F2(I)=F(I)+A(1,I)/2.
      CALL OP(F2,F1,Y+H/2.)
      DO 2 I=1,M
      A(2,I)=H*F1(I)
2     F2(I)=F(I)+A(2,I)/2.
      CALL OP(F2,F1,Y+H/2.)
      DO 3 I=1,M
      A(3,I)=H*F1(I)
3     F2(I)=F(I)+A(3,I)
      CALL OP(F2,F1,Y+H)
      DO 4 I=1,M
      A(4,I)=H*F1(I)
4     F(I)=F(I)+(A(1,I)+2.*(A(2,I)+A(3,I))+A(4,I))/6.
      RETURN
      END

