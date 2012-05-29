subroutine SEARCH_INITIAL_INSTAB_TIME
implicit none
common /HADY2/A,B,W,R
complex*16 A,B,W,R
call TS_GLOBAL_TIME()
print *, "GLOBAL SEARCH:"
WRITE(*,'(1P,/," W=",2E13.5, &
             /," B=",2E13.5,/," A=",2E13.5,  &
             /," PHASE SPEED=",E13.5)') W,B,A,DREAL(W)/DREAL(A)

! some logic  to choose...

!call CF_GLOBAL_SEARCH()
!print *, "GLOBAL SEARCH:"
!WRITE(*,'(1P,/," RRR=",E13.5,/," W=",2E13.5, &
!             /,"  B=",2E13.5,/," A=",2E13.5,  &
!             /," PHASE SPEED=",E13.5)')RRR,W,B,A,DREAL(W)/DREAL(A)
!
end subroutine SEARCH_INITIAL_INSTAB_TIME