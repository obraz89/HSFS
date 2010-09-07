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
COMMON/ADDITIONAL_NS/ X_DIM, Y_DIM
COMMON/CONTROL/ REQ_TS_GLOB     ! this is bad 
COMMON/SOLVER_OUTPUT/ SIGMA_SPAT
COMPLEX*16 A,B,W,R,FUR,FU(8,201),           &
           CF(4),SI,W0,V,Z0,A0,             &
           AEXTR,VA,VB,VR,FREQ
           
REAL*8 K,M,MF_X, MF_WE, MF_ANG,PI, SIGMA_SPAT               
integer MF_I, MF_J, Y_DIM, X_DIM, REQ_TS_GLOB

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


! Additional mean-flow parameters
!
IBAS=0 ! If IBAS=1 then write mean-flow profiles in first X-station and stops
XI=0.D0
PB=0.D0
CQ=0.D0
XIR=XI*PI/180.
PB0=PB
IHH=0
JHH=0
IA=0
Z0=0.
MAB0=1
MAB=MAB0
MAB1=1
NJOB=2
IPRINT=0

! Specify parameters for stability solver
!
NO=150  !250  ! number of grid points
NS=1    ! coefficient increasing number of grid points 
IGR=100	  ! index of poisition where eigenfun is to be plotted
      
NYY=Y_DIM      
XXX=XST(1)
CALL NAVSTOK(NYY,XXX)

RE1=REE/(XXX*XLL)		! Local unit Reynolds number
RRR=DSQRT(REE)  ! Reynolds number for stability calculations
R=RRR
LP=1 ! if LP=1 then search of alpha at fixed omega !obsolete, unused

CW='V'
A0=A
W0=W
MAB0=1
MAB1=0
BBVB=BVB
MAB=MAB0
MAB2=MAB

DELS=RRR/RE1	! local length-scale
! if it is first station ....
IF (REQ_TS_GLOB==1) THEN
call TS_GLOBAL_TIME()
print *, "GLOBAL SEARCH:"
WRITE(*,'(1P,/," RRR=",E13.5,/," W=",2E13.5, &
             /,"  B=",2E13.5,/," A=",2E13.5,  &
             /," PHASE SPEED=",E13.5)')RRR,W,B,A,DREAL(W)/DREAL(A)

REQ_TS_GLOB=0
endif
CALL TS_WMAX_TIME_ALL()
call TS_TRANS_TO_SPATIAL(a_spat, b_spat, w_spat)
sigma_spat = -dimag(a_spat)/dels       
print *,'time elapsed:', etime(elapsed_time)

end subroutine SEARCH_MAX_INSTAB_TIME