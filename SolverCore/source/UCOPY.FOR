
      SUBROUTINE UCOPY(A,B,N)                                      
C
C   COPY OF N ELEMENTS OF VECTOR A INTO VECTOR  B 
C
      IMPLICIT REAL*8 (A-H,O-Z)                                    
      DIMENSION A(N),B(N)                                          
      CALL UTEST(A,N)                                              
      DO 1 I=1,N                                                   
      B(I)=A(I)                                                    
1     CONTINUE                                                     
      RETURN                                                       
      END                                                          

