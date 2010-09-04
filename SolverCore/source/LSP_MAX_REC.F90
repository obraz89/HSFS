! SEE LSP_MAX.F90
! recursive modification (performance is critical)

recursive subroutine LSP_MAX_REC(point_init, b_range_init)
implicit none

complex*16, intent(in) :: point_init, b_range_init 

interface 
  function find_min(fun)
    real(8) fun(:)
    integer find_min
  end
end interface

COMMON /HADY2/ ALPHA, BETA, W, R              &    
       /HADY3/ SI,Z0                          &    
       /ABOCT/ MAB,MAB1,MAB2                  &  
       /JP/ JPOISK                            &  
       /DNS/ AMINF,REINF,TINF,XLL,REE,REE1,UEE
!temp COMMON
COMMON/LSP_OUT/GR_DIR
!
!-------------------------------------globals (COMMON)
complex*16 ALPHA, BETA, W, R
complex*16 SI,Z0
integer MAB, MAB1, MAB2, JPOISK       
real(8) AMINF,REINF,TINF,XLL,REE,REE1,UEE, GR_DIR
!-------------------------------------globals

complex*16 b_start, b_end,                    &
           b_prv, a_prv, point, b_range
complex*16, allocatable:: a(:,:), b(:,:)
real(8) gv_dir, db_r, db_i, step_r, step_i
real(8), allocatable :: nar(:), sl_a(:), sl_b(:)
integer, allocatable :: bp_ind(:)
integer n, i, j, ind

!"template" parameter
! greater n - more precise calculation, but drastic slow
n=5
step_r = DREAL(b_range_init)/real(n)
step_i = DIMAG(b_range_init)/real(n)
point = point_init
b_start=point_init - 0.5*b_range_init
b_end=point_init + 0.5*b_range_init
allocate (a(n,n), b(n, n),                & 
          sl_a(n), sl_b(n),               &
          bp_ind(n), nar(n))                 

db_r = dreal(b_end-b_start)/real(n-1)
db_i = dimag(b_end-b_start)/real(n-1)

do i=1,n
  do j=1, n
    b(i,j)=b_start + dcmplx(real(i-1)*db_r, real(j-1)*db_i)      
    BETA=b(i,j)
    CALL POISK2(1)  
    a(i,j)=ALPHA 
  end do
end do

do j=1,n
  do i=1, n
    sl_a(i)=DIMAG(a(i,j))
    sl_b(i)=DREAL(b(i,j))
  end do  
bp_ind(j) = find_min(sl_a)
ind = bp_ind(j)
gv_dir = DREAL(a(ind+1,j)-a(ind-1,j))/ &
         DREAL(b(ind+1,j)-b(ind-1,j))
nar(j) = DIMAG(a(ind,j))-gv_dir*DIMAG(b(ind,j))
end do

j=find_min(nar)
ALPHA=a(bp_ind(j),j)
BETA =b(bp_ind(j),j)
!temp
gv_dir = DREAL(a(bp_ind(j)+1,j)-a(bp_ind(j)-1,j))/ &
         DREAL(b(bp_ind(j)+1,j)-b(bp_ind(j)-1,j))
print *, "gv_dir:", gv_dir
GR_DIR = DIMAG(a(ind,j))-gv_dir*DIMAG(b(ind,j))
!

if (step_r<0.0001.AND.step_i<0.00002) then
  write (*,*) "LSP---A,B,W:"
  write (*,'(2E13.5,/,2E13.5,/,2E13.5,/,"---")')   &
  ALPHA, BETA, W
else
  !print *, 'iteration used, BETA=', BETA
  b_range = dcmplx(step_r,step_i)   !!!2.0
  point=BETA
  CALL LSP_MAX_REC(point, b_range)
endif
deallocate (a, b, sl_a, sl_b, bp_ind, nar)
end subroutine LSP_MAX_REC

function find_min(fun)
real(8) fun(:)
integer i, n, find_min
n = SIZE(fun)
find_min=1
do i=2, n
 if (fun(i)<=fun(find_min)) then
  find_min = i 
 end if
 if (find_min==1) find_min=2
 if(find_min==n) find_min=n-1
end do

end function find_min

