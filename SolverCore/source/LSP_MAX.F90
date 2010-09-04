!------------------------------------------------------------
! LSP_MAX - "Local Saddle Point" Maximization
!
! В пространстенной постановке решается задача нахождения
! инкремента нарастания наиболее неустойчивой 
! собственной моды и направления 
! ее нарастания (групповой скорости) монохроматического
! возмущения. Для устранения 
! неоднозначности в дисперсионном соотношении используются 
! соотношения:
! da_i/db_r=0;
! -(a_i + da_r/db_r * b_i) -> max
!------------------------------------------------------------
subroutine LSP_MAX
implicit none
interface 
  function find_min(fun)
    real(8) fun(:)
    integer find_min
  end
end interface

COMMON /HADY2/ ALPHA, BETA, W, R              &    
       /HADY3/ SI,Z0                          &    
       /ABOCT/ MAB,MAB1,MAB2                  &  ! ?
       /JP/ JPOISK                            &  ! ? 
       /DNS/ AMINF,REINF,TINF,XLL,REE,REE1,UEE
!-------------------------------------globals (COMMON)
complex*16 ALPHA, BETA, W, R
complex*16 SI,Z0
integer MAB, MAB1, MAB2, JPOISK       
real(8) AMINF,REINF,TINF,XLL,REE,REE1,UEE
!-------------------------------------globals

complex*16 b_start, b_end, b_prv,             &
            a_prv
complex*16, allocatable:: a(:,:), b(:,:)
real(8) gv_dir,db_r, db_i
real(8), allocatable :: nar(:),sl_a(:), sl_b(:)
integer, allocatable :: bp_ind(:)
integer nr, ni, i, j, ind
b_start=dcmplx(-0.2, -0.1)
b_end=dcmplx(0.2, 0.1)
nr=20
ni=20
allocate (a(nr,ni), b(nr, ni),                & 
          sl_a(nr), sl_b(nr),                 &
          bp_ind(ni), nar(ni))                 

db_r = dreal(b_end-b_start)/real(nr)
db_i = dimag(b_end-b_start)/real(ni)

do i=1,nr
  do j=1, ni
    b(i,j)=b_start + dcmplx(real(i-1)*db_r, real(j-1)*db_i)      
    BETA=b(i,j)
    CALL POISK2(1)
    a(i,j)=ALPHA 
  end do
end do




do j=1,ni
  do i=1, nr
    sl_a(i)=DIMAG(a(i,j))
    sl_b(i)=DREAL(b(i,j))
  end do  
bp_ind(j) = find_min(sl_a)
ind = bp_ind(j)
gv_dir = DREAL(a(ind+1,j)-a(ind-1,j))/ &
         DREAL(b(ind+1,j)-b(ind-1,j))
nar(j) = DIMAG(a(ind,j))-gv_dir*DIMAG(b(ind,j))
end do

!open(UNIT=123,FILE="LSP_example.dat")
!do i = 2, nr-1
!  do j = 2, ni-1
!  write(123, *) DREAL(B(i,j)), DIMAG(B(i,j)), -DIMAG(a(i,j)) +         &
!                                             DREAL(a(i+1,j)-a(i-1,j))/ &
!                                             DREAL(b(i+1,j)-b(i-1,j))* &
!                                             DIMAG(b(i,j))
!  end do
!end do
!close(123)

j=find_min(nar)
ALPHA=a(bp_ind(j),j)
BETA =b(bp_ind(j),j)
write (*,"(2E13.5,/,2E13.5,/,2E13.5)")   &
 ALPHA, BETA, DREAL(W)/DSQRT(REE),DREAL(W/ALPHA)

open(UNIT=1122, file="LSP.dat", ACCESS='APPEND')
write(1122,'(E15.3, E15.3)') ALPHA
close(1122)
deallocate (a, b, sl_a, sl_b, bp_ind, nar)
end subroutine LSP_MAX