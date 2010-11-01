module TIME_SEARCH
implicit none
real(8) TGT_WR                                             
contains
! moving along integration path
! values from previous spatial
!point should be adjusted for current point
subroutine ADJUST_INITIAL_VALUES()
implicit none
common /HADY2/ A, B, W, R           &
       /HADY3/ SI,ZZZ0                                     
complex*16 A, B, W, R, SI, ZZZ0
complex*16 a_def, b_def, w_def
real fsi_cur, fsi_prv, da, db, dwr, dwi
! to enter loop
fsi_prv = 10
fsi_cur = 1
do while (fsi_cur<fsi_prv)
call HADSI
fsi_prv = cdabs(SI)
a_def = A
b_def = B
w_def = W

da = 0.001*A
A = a_def + da
call HADSI
if (cdabs(SI)>fsi_prv) da = - da
A = a_def

db = 0.001*B
B = B + db
call HADSI
if (cdabs(SI)>fsi_prv) db =  - db
B = b_def

dwr = 0.001*dreal(W)
W = w_def + dwr
call HADSI
if (cdabs(SI)>fsi_prv) dwr =  - dwr 
W = w_def

dwi = 0.001*dimag(W)
W = w_def + dcmplx(0, dwi)
call HADSI
if (cdabs(SI)>fsi_prv) dwi =  - dwi
W = w_def
A = A + da
B = B + db
W = W + dcmplx(dwr, dwi)
call HADSI
fsi_cur = cdabs(SI)
print*, "internal adjust: FSI=", fsi_cur
end do
end subroutine
! search for maximum instability
! in time approach with keeping Wr constant
! criterion - Function GR_VEL_FUN = VBi*VAr - VAi*VBr
subroutine TIME_WMAX_STAT()
!use sys_functs
implicit none
common /HADY2/ A, B, W, R           &
       /VGRC/ VA,VB,VR             
complex*16 A, B, W, R,              &
           VA, VB, VR
complex*16 a_def, b_def,            &
           a_prv, b_prv, w_prv,     &
           a_cur, b_cur, w_cur
real(8) db_def, da, db,             &
        gr_fun_prv, gr_fun_cur,     & 
        wi_prv, wi_cur, wi_def
call POISK2(3)
a_def = A
b_def = B
wi_def = DIMAG(W)
db = 0.01*B             ! constant <>?
gr_fun_prv = GR_VEL_FUN()
wi_prv = DIMAG(W)
da = -DREAL(VB)/DREAL(VA)*db
B = b_def + db
A = a_def + da
gr_fun_cur = GR_VEL_FUN()
wi_cur = DIMAG(W)
!if (DABS(gr_fun_cur)>DABS(gr_fun_prv)) then
if (wi_cur<wi_prv) then
  db = -db
  da = -DREAL(VB)/DREAL(VA)*db
  B = b_def + db
  A = a_def + da
  gr_fun_cur = GR_VEL_FUN()
  wi_cur = DIMAG(W)
end if
!do while (DABS(gr_fun_cur)<DABS(gr_fun_prv)) 
do while (wi_prv<wi_cur) 
!print*, "WI, WR:", wi_cur, DREAL(W)
da = -DREAL(VB)/DREAL(VA)*db
B = B + db
A = A + da
gr_fun_prv = gr_fun_cur
wi_prv = wi_cur
gr_fun_cur = GR_VEL_FUN()
wi_cur = DIMAG(W)
if (DREAL(W)<0.0d0) W=DCMPLX(0.0d0, DIMAG(W)) 
end do 
contains
  function GR_VEL_FUN
  implicit none
  common /HADY2/ A, B, W, R               &
       /VGRC/ VA, VB, VR
  complex*16 A,B,W,R,                     &
           VA, VB, VR
  real(8) GR_VEL_FUN
  call POISK2(3)
  call VGR
  GR_VEL_FUN = DIMAG(VB)*DREAL(VA) - DIMAG(VA)*DREAL(VB)
  end function GR_VEL_FUN
end subroutine TIME_WMAX_STAT
!#######################################################
!maximize wi versus all frequencies wr
subroutine TIME_WMAX_ALL
common /HADY2/ A, B, W, R
complex*16 A, B, W, R
complex*16 a_prv, b_prv, w_prv,     &
           a_cur, b_cur, w_cur
real(8) dar
!call POISK2(3)
call ADJUST_INITIAL_VALUES()
call TIME_WMAX_STAT()
a_prv = A
b_prv = B
w_prv = W
dar = 0.005*DREAL(a_prv)
call TIME_SHIFT_FREQ(dar)
! if we fall in a wrong place
!if (dimag(W)<-3.0d-4) call TS_GLOBAL_TIME()
call TIME_WMAX_STAT()
a_cur = A
b_cur = B
w_cur = W
if (DIMAG(w_prv)>DIMAG(w_cur)) then
  A = a_prv
  B = b_prv
  W = w_prv
  dar=-dar
  call TIME_SHIFT_FREQ(dar)
  call TIME_WMAX_STAT()
  a_cur = A
  b_cur = B
  w_cur = W
end if
do while (DIMAG(w_cur)>DIMAG(w_prv))
  print *, "TIME_WMAX_ALL iter:", W
  A = a_cur
  B = b_cur
  W = w_cur
  a_prv = a_cur
  b_prv = b_cur
  w_prv = w_cur
  call TIME_SHIFT_FREQ(dar)
  call TIME_WMAX_STAT()
  a_cur = A
  b_cur = B
  w_cur = W
end do
end subroutine TIME_WMAX_ALL
!#######################################################
! if we know  TimEigenset (w,a,b)
! initial approach may be (w+dw/da*da, a+da, b)
subroutine TIME_SHIFT_FREQ(dar)
common /HADY2/ A, B, W, R               &
       /VGRC/ VA, VB, VR
complex*16 A,B,W,R,                     &
           VA, VB, VR
real(8), intent(in) :: dar
call POISK2(3)
call VGR()
W = W + VA*dar
A = A + dar
call POISK2(3)
end subroutine TIME_SHIFT_FREQ
!#######################################################
subroutine TRANS_TO_SPATIAL(a_spat, b_spat, w_spat)
implicit none
common /HADY2/ A, B, W, R           &
       /VGRC/ VA,VB,VR 
complex*16, intent(inout) :: a_spat, b_spat, w_spat             
complex*16 A, B, W, R,              &
           VA, VB, VR,              &
           a_save, b_save, w_save
complex*16 cos_fi, sin_fi,     &
           c, dk, dw, da, db
real(8) psi
call POISK2(3)
call VGR
print *, "group vel: VA, VB:", VA, VB
! set rough spat eigen values
!c = CDSQRT(VA**2.0 + VB**2.0)
!cos_fi = VA/c
!sin_fi = VB/c
!psi = atan(DREAL(VB)/DREAL(VA))
dw = DCMPLX(0.0d0, -DIMAG(W))
!dk = dw*1.0d0/(c*(cos_fi*dcos(psi) + sin_fi*dsin(psi)))
!da = dk*dcos(psi)
!db = dk*(dsin(psi)) ! ????
da = -1.0*DCMPLX(0.0d0, 1.0d0)*DIMAG(W)/VA
db = da*DREAL(VB)/DREAL(VA)
a_save = A
b_save = B
w_save = W
A = A + da
B = B + db
W = W + dw
! precise solution
call POISK2(1)
a_spat = A
b_spat = B
w_spat = W
A = a_save
B = b_save
W = w_save
end subroutine TRANS_TO_SPATIAL
!######################################
subroutine TRANS_TO_TIME(a_time, b_time, w_time)
implicit none
common /HADY2/ A, B, W, R           &
       /VGRC/ VA,VB,VR 
complex*16, intent(inout) :: a_time, b_time, w_time             
complex*16 A, B, W, R,              &
           VA, VB, VR,              &
           a_save, b_save, w_save
complex*16 cos_fi, sin_fi,     &
           c, dk, dw, da, db
real(8) psi
call POISK2(1)
call VGR
print *, "group vel: VA, VB:", VA, VB
! set rough spat eigen values
!c = CDSQRT(VA**2.0 + VB**2.0)
!cos_fi = VA/c
!sin_fi = VB/c
!psi = atan(DREAL(VB)/DREAL(VA))
da = DCMPLX(0.0d0, -DIMAG(A))
db = DCMPLX(0.0d0, -DIMAG(B))
dw = VA*da + VB*db
!dk = dw*1.0d0/(c*(cos_fi*dcos(psi) + sin_fi*dsin(psi)))
!da = dk*dcos(psi)
!db = dk*(dsin(psi)) ! ????
a_save = A
b_save = B
w_save = W
A = A + da
B = B + db
W = W + dw
! precise solution
call POISK2(3)
a_time = A
b_time = B
w_time = W
A = a_save
B = b_save
W = w_save
end subroutine TRANS_TO_TIME

subroutine SEARCH_DIFF_MODES(alpha, beta, w_init)
implicit none
common /HADY2/ A, B, W, R           &
       /VGRC/ VA,VB,VR              &
       /HADY3/ SI,Z0                &
       /JP/JPOISK 
complex*16, intent(in) :: alpha, beta, w_init
complex*16 A, B, W, R,              &
           VA, VB, VR,              &
           SI, Z0,                  &
           a_save, b_save, w_save,  &
           w_start, w_end
integer JPOISK, ir, ii,i, j

ir=1000
ii=1000
w_start = 2.0*dreal(w_init)
w_end = w_start + 2.0*w_init
a_save = A
b_save = B
w_save = W
A = alpha
B = beta
W = w_init
JPOISK=0
OPEN(file="output/search_diff_modes.dat", unit=4321)
WRITE (4321, '("DIFF_MODES: ","A=",2E13.5,/,"B=",2E13.5,/,"W_INIT=", 2E13.5,/,"Wr   Wi  FSI",/)') alpha, beta, w_init
do i=1, ir
write (6,fmt='(1E12.5)'), real(i)/real(ir)
  do j=1, ii 
    W = dcmplx(dreal(w_start)+real(i-1)/real(ir-1)*dreal(w_end - w_start),  &
               dimag(w_start)+real(j-1)/real(ii-1)*dimag(w_end - w_start))
    call HADSI
    WRITE (4321,'(3E13.5)')  W, CDABS(SI)  
  end do
end do
A = a_save
B = b_save
W = w_save
! debug
close(4321)
stop
end subroutine SEARCH_DIFF_MODES

end module TIME_SEARCH 