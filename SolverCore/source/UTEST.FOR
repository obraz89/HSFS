
      SUBROUTINE UTEST(A,N)                                           
      IMPLICIT REAL*8 (A-H,O-Z)                                                      
      COMPLEX*16 A(N)
      DATA C/0./                                                  
      DO 1 I=1,N                                                      
      P=CDABS(A(I))                                                    
      IF(P.LT.1.D19) GO TO 2                                          
      K=N                                                             
      IF(K.GT.80) K=80                                                
      PRINT 10,N,I,A(I),(A(J),J=1,K)                                  
10    FORMAT(' PEREPOLNENIE-UTEST N=',I5,' A(',I4,')=',1P2D9.2,' A='/  
     *11(4(2X,2D9.2)/))                                               
      P=1./C                                                          
      STOP                                                            
2     IF(P.LT.1.D-19) A(I)=0.                                         
1     CONTINUE                                                        
      RETURN                                                          
      END                                                             

      FUNCTION AMAXMU(A,N)
C*******************************************************
C   COMPUTATION OF MAXIMUM (BY MODUL) ELEMENT OF MASSIV
C     AMAXMU= MAX (DABS(A(I)))
C             I=1,N
C*******************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(N)
      R=0.D0
      DO 1 I=1,N
      IF(R.GT.DABS(A(I))) GO TO 1
      R=DABS(A(I))
1     CONTINUE
      AMAXMU=R
      RETURN
      END

      FUNCTION SN(A,N)
C *************************************************
C COMPUTATION OF NORM OF VECTOR A WITH DIMENSION N.
C *************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 A(N)
      CALL UTEST(A,2*N)
      SN=0.
      DO 1 I=1,N
1     SN=SN+DCONJG(A(I))*A(I)
      SN=DSQRT(SN)
      RETURN
      END

      FUNCTION SP(A,B,M)
C ********************************
C
C  COMPUTATION OF SCALAR PRODUCT 
C  OF TWO VECTORS: (A,B)
C ********************************
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16  SP,A(M),B(M)
      SP=0.
      DO 1 I=1,M
1     SP=SP+DCONJG(A(I))*B(I)
      RETURN
      END
