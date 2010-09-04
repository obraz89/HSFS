      FUNCTION SIM(F,N,D)                                             
C **********************************************************                
C  
C  N MUST BE ODD AND > 4
C  SIM=INTEGRAL OF SIMPSON FROM OT "F", WITH STEP "D".    
C***********************************************************
      IMPLICIT REAL*8 (A-H,O-Z)                             
      COMPLEX*16 F(N),SIM,A                                 
      SIM=D/3.*(F(1)+F(N))                                  
      A=0.                                                  
      N1=N-1                                                
      N2=N-2                                                
      DO 1 I=2,N1,2                                         
1     A=A+F(I)                                              
      SIM=SIM+4.*D/3.*A                                     
      A=0.                                                  
      DO 2 I=3,N2,2                                         
2     A=A+F(I)                                              
      SIM=SIM+2.*D/3.*A                                     
      RETURN                                                
      END                                                   

