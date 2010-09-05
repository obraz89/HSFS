subroutine SEARCH_MAX_INSTAB_TIME
implicit real*8 (A-H, O-Z)
integer, parameter:: NY = 541
COMMON/NS/YNS(NY),UNS(NY),UNS1(NY),UNS2(NY),    &
                  TNS(NY),TNS1(NY),TNS2(NY),    &
                  WNS(NY),WNS1(NY),WNS2(NY),    &
                  ENTNS(NY),FINS(NY),DFINS(NY),XST(NY)
COMMON/RHONS/ RHO(NY)                           
COMMON/OUT  /K,SIG,G,M,XI,PB                    &
            /HADY1/Y1,D,NS,L1                   &
            /MACK02/TELOC,TWALL,TF              &
            /MACKY1/ Y1O,DO,C,YC,NO             &
            /MACKNJ/NJOB,IPRINT,LOLD,IV         &
            /MACKXY/XFIRST,YFIRST,ZFIRST        &
            /ABOCT/ MAB,MAB1,MAB2               &
            /BASA6/ A6(10),BVB,DL1              &
            /BASAP/UPP,DLP                      &
            /BASA3/A7(9),AMK,TETB,A8(4),III(2)  &   
            /HADY2/ A,B,W,R                     &
            /HADY5/ FU,CF                       &
            /HADY3/ SI,Z0                       &
            /AK/ AK                             &
            /OSNOB/YYY,U,U1,U2,T,T1,T2,         &
                   VMU,VMU1,VMU2,WS,WS1,WS2     &
            /VGRC/ VA,VB,VR                     &
            /osncq/cq                           &
            /IUPT/IUPT                          &
	        /DNS/AMINF,REINF,TINF,              &
	             XLL,REE,REE1,UEE               &
	        /JP/JPOISK
	!temp common
COMMON/CF_PROFILE_DATA/ MF_I, MF_J, MF_X, MF_WE, MF_ANG
COMPLEX*16 A,B,W,R,FUR,FU(8,201),           &
           CF(4),SI,W0,V,Z0,A0,             &
           AEXTR,VA,VB,VR,FREQ
           
REAL*8 K,M,MF_X, MF_WE, MF_ANG,PI              
integer MF_I, MF_J

!locals          
real*4, dimension(2) :: elapsed_time
complex*16 a_spat, b_spat, w_spat, a_time, b_time, w_time

PI=3.1415926535897932384626433832795D0
! Specify gas parameters
! IUPT=0 - power law for viscosity mu~T^BVB
! IUPT=1 - Satherland law for viscosity
! Second (bulk) viscosity coefficient (Mack recommends K=1.2)
G=1.4D0   ! specific heat ratio
SIG=0.72   ! Prandtl number
RVAP=287.1    ! gas constant per unit mass [J/kg/K]
IUPT=1
BVB=0.75
K=0.

end subroutine SEARCH_MAX_INSTAB_TIME