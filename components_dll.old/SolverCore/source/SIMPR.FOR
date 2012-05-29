      FUNCTION SIMPR(F,Y,N)        
C ***********************************************
C
C  N MUST BE ODD AND > 4
C  SIM=INTEGRAL OF SYMPSON FROM "F"
C************************************************
      IMPLICIT REAL*8 (A-H,O-Z)                  
      REAL*8 F(N),Y(N)                           
      A=0.                                       
      N1=N-1                                     
      DO 1 I=1,N1                                
1     A=A+(F(I)+F(I+1))*(Y(I+1)-Y(I))*5.D-1      
      SIMPR=A                                    
      RETURN                                     
      END                                        

