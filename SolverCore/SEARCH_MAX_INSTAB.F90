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
            /MACKY1/ Y1O,DO,C,YC,NO             &
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
            /IUPT/IUPT                          &
	        /DNS/AMINF,REINF,TINF,              &
	             XLL,REE,REE1,UEE               &
	        /JP/JPOISK
	!temp common
COMMON/CF_PROFILE_DATA/ MF_I, MF_J, MF_X, MF_WE, MF_ANG
COMMON/ADDITIONAL_NS/ X_DIM, Y_DIM
COMMON/SOLVER_OUTPUT/ SIGMA_SPAT, A_SPAT, B_SPAT, W_SPAT 
COMPLEX*16 A,B,W,R,FUR,FU(8,201),           &
           CF(4),SI,V,Z0,                   &
           AEXTR,VA,VB,VR,FREQ
           
REAL*8 K,M,MF_X, MF_WE, MF_ANG               
integer MF_I, MF_J, Y_DIM, X_DIM, REQ_TS_GLOB, REQ_CF_GLOB

!locals          
real*4, dimension(2) :: elapsed_time
complex*16 SIGMA_SPAT, A_SPAT, B_SPAT, W_SPAT, a_time, b_time, w_time


!Z0=0. ! notused part of hady3
!MAB=1

! Specify parameters for stability solver
!
!NO=150  !250  ! number of grid points
!NS=1    ! coefficient increasing number of grid points             

!XXX=XST(1)
!CALL NAVSTOK(Y_DIM,XXX)

RE1=REE/(XST(1)*XLL)		! Local unit Reynolds number
RRR=DSQRT(REE)  ! Reynolds number for stability calculations
!R=RRR

DELS=RRR/RE1	! local length-scale

CALL TS_WMAX_TIME_ALL() ! the same for CF
!call TS_WMAX_TIME_STAT()
call TS_TRANS_TO_SPATIAL(A_SPAT, B_SPAT, W_SPAT)
sigma_spat = -dimag(a_spat)/dels       
print *,'time elapsed:', etime(elapsed_time)

end subroutine SEARCH_MAX_INSTAB_TIME