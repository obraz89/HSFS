! Global search of TS wave if JPOISK=1;
! Local search of most unstable waves
! SI=0 <->  eigenvalues found
! POISK_TS - LSP ! saptial approach
subroutine POISK_TS(LP)

IMPLICIT REAL*8 (A-H,O-Z)

common /HADY2/ Z(4)                         &
            /HADY3/ SI,ZZZ0                 &
            /ABOCT/ MAB,MAB1,MAB2           &
            /AK/ AK                         &
            /BASA6/ A6(10),BVB,DL1  
common/OUT  /AAK,AASIG,AAG,AAM,AAXI,AAPB
common/JP/JPOISK
complex*16 SI,Z0,F,F1,Z00,ZB,Z,ZZZ0

real*8, parameter :: pi = 3.141592653
real*8, parameter :: theta = pi/4.0
complex*16 point, b_range

data EPS /1.D-8/,D/1.D-7/

!------------ Global search of TS mode -------      

IF(JPOISK.EQ.1)THEN
WLOW=0.05	                    ! initial freuency
WUP=0.9		                    ! final frequncy
AI=-1.E-3	                    ! growth rates
BI=-1.E-4
CFUP=1.		                    ! maximum phase speed
CFLOW=1.-1./(AAM*cos(theta))	! minimum phase speed
IF (CFLOW<0.0) CFLOW = 0.1
NW=101
NAR=101
DW=(WUP-WLOW)/(NW-1)
FMIN=1.E+9
DO JW=1,NW
  W=WLOW+DW*(JW-1)
  ARLOW=W/CFUP
  ARUP=W/CFLOW
  DAR=(ARUP-ARLOW)/(NAR-1)
  WRITE(2,'(/)')
  DO JAR=2,NAR-1
    AR=ARLOW+DAR*(JAR-1)
    BR=AR*tan(theta)
    Z(3)=W
    Z(1)=DCMPLX(AR,AI)
    Z(2)=DCMPLX(BR,BI)
    CALL HADSI
    FSI=CDABS(SI)
    IF(FSI.LT.FMIN)THEN
    print *, "FSI=", FSI
      FMIN=FSI
      ARMIN=AR
      BRMIN=BR
      AIMIN=AI
      BIMIN=BI
      WMIN=W
    ENDIF
  ENDDO
ENDDO
!WRITE(*,'(1P,20E14.6)')ARMIN,AIMIN,WMIN,WMIN/ARMIN,FMIN
Z(1)=DCMPLX(ARMIN,AIMIN)
Z(2)=DCMPLX(BRMIN,BIMIN)
Z(3)=WMIN
JPOISK=0
PAUSE
ENDIF

!-------------------------------------------------------end of global search
! local search of TS mode
b_range = DCMPLX(0.1*DREAL(Z(2)), DABS(DIMAG(Z(2))))    ! range is changed
point= Z(2)
CALL LSP_MAX_REC(point, b_range)
pause
RETURN
END

!just wrapping for TIME_WMAX_STAT - search max wi vs (a,b)
!keeping wr const (monochrome maximization)
subroutine TS_WMAX_TIME_STAT()
use TIME_SEARCH
call TIME_WMAX_STAT()
end subroutine TS_WMAX_TIME_STAT

subroutine TS_WMAX_TIME_ALL()
use TIME_SEARCH
call TIME_WMAX_ALL()
end subroutine TS_WMAX_TIME_ALL

subroutine TS_TRANS_TO_SPATIAL(a_spat, b_spat, w_spat)
use TIME_SEARCH
complex*16, intent(inout) :: a_spat, b_spat, w_spat
call TRANS_TO_SPATIAL(a_spat, b_spat, w_spat)
end subroutine TS_TRANS_TO_SPATIAL

subroutine TS_SEARCH_DIFF_MODES(a_time, b_time, w_time_init)
use TIME_SEARCH
complex*16, intent(in) :: a_time, b_time, w_time_init
call SEARCH_DIFF_MODES(a_time, b_time, w_time_init)
end subroutine TS_SEARCH_DIFF_MODES



