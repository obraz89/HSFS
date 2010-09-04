subroutine smooth_3D_profiles
include 'link_f90_static.h'

use csscv_int
use cs1gd_int

implicit none
real(8), parameter :: Re=1.424E+07, Ma=3.5, gamma=1.4
real(8), parameter :: L=0.381, Tinf=90.318
real(8), parameter :: pi = 3.1415926535897932384626433832795d0
real(8) :: alfa=5.d0*pi/180.d0

integer, parameter :: ny=161, n_bt=80, n_up=81
character*(*), parameter :: aoa_case="al=2", k_case="k=5"

! from NS solver
real(8), dimension(ny) :: y,u,w,p,t,r
! for bottom splines                          
real(8), dimension(n_bt) :: y_bt, u_bt, w_bt, t_bt, &
                            u1_bt, u2_bt,           &
                            w1_bt, w2_bt,           &
                            t1_bt, t2_bt,           &
                            ub_bt,wb_bt,tb_bt
real(8), dimension(4,n_bt) :: ucoef_bt,wcoef_bt,tcoef_bt
! for upper splines
real(8), dimension(n_up) :: y_up, u_up, w_up, t_up,  &
                            u1_up, u2_up,            &
                            w1_up, w2_up,            &
                            t1_up, t2_up,            &
                            ub_up,wb_up,tb_up
real(8), dimension(4,n_up) :: ucoef_up,wcoef_up,tcoef_up
integer size, vsize
real(8) ue(201), we(201), angle(201),x(201)
integer i,j, iequal
!------------------------------------------------------
OPEN(123, FILE='input/new_slices/KINGssd_'//k_case//'_sz=75_'//aoa_case//'_nr_new.dat',  &
     FORM='BINARY',STATUS='OLD')
OPEN (1,FILE='output/new/CF_'//k_case//'_'//aoa_case//'_LONG(75)_NR_new.dat')
read(123) size,vsize
print *, "mfdata file : nx,ny:", size, vsize
pause
do i = 1, size
  read(123) x(i),ue(i), we(i), angle(i)
  do j=1,ny
    READ(123) y(j)
    READ(123) u(j)
    READ(123) w(j)
    READ(123) p(j)
    READ(123) t(j)
  end do
  r(1:ny)=gamma*Ma*Ma*p(1:ny)/t(1:ny)
  y(1:ny) = dsqrt(Re)*y(1:ny) ! нормировка на погран слой
  y_bt(1:n_bt) = y(1:n_bt)
  y_up(1:n_up) = y(n_bt+1:ny)
  u_bt(1:n_bt) = u(1:n_bt)
  u_up(1:n_up) = u(n_bt+1:ny)
  w_bt(1:n_bt) = w(1:n_bt)
  w_up(1:n_up) = w(n_bt+1:ny)
  t_bt(1:n_bt) = t(1:n_bt)
  t_up(1:n_up) = t(n_bt+1:ny)
iequal = 0

!bottom splines
  call d_csscv(y_bt, u_bt, iequal, ub_bt, ucoef_bt) 
  call d_csscv(y_bt, w_bt, iequal, wb_bt, wcoef_bt) 
  call d_csscv(y_bt, t_bt, iequal, tb_bt, tcoef_bt) 

  call d_cs1gd(1, y_bt, ub_bt, ucoef_bt, u1_bt) 
  call d_cs1gd(1, y_bt, wb_bt, wcoef_bt, w1_bt) 
  call d_cs1gd(1, y_bt, tb_bt, tcoef_bt, t1_bt)

  call d_csscv(y_bt, u1_bt, iequal, ub_bt, ucoef_bt)
  call d_csscv(y_bt, w1_bt, iequal, wb_bt, wcoef_bt)
  call d_csscv(y_bt, t1_bt, iequal, tb_bt, tcoef_bt)
  
  call d_cs1gd(0, y_bt, ub_bt, ucoef_bt, u1_bt) 
  call d_cs1gd(0, y_bt, wb_bt, wcoef_bt, w1_bt)
  call d_cs1gd(0, y_bt, tb_bt, tcoef_bt, t1_bt)

  call d_cs1gd(1, y_bt, ub_bt, ucoef_bt, u2_bt) 
  call d_cs1gd(1, y_bt, wb_bt, wcoef_bt, w2_bt)
  call d_cs1gd(1, y_bt, tb_bt, tcoef_bt, t2_bt)
  
!upper slice
  call d_csscv(y_up, u_up, iequal, ub_up, ucoef_up) 
  call d_csscv(y_up, w_up, iequal, wb_up, wcoef_up) 
  call d_csscv(y_up, t_up, iequal, tb_up, tcoef_up) 

  call d_cs1gd(1, y_up, ub_up, ucoef_up, u1_up) 
  call d_cs1gd(1, y_up, wb_up, wcoef_up, w1_up) 
  call d_cs1gd(1, y_up, tb_up, tcoef_up, t1_up)

  call d_csscv(y_up, u1_up, iequal, ub_up, ucoef_up)
  call d_csscv(y_up, w1_up, iequal, wb_up, wcoef_up)
  call d_csscv(y_up, t1_up, iequal, tb_up, tcoef_up)
  
  call d_cs1gd(0, y_up, ub_up, ucoef_up, u1_up) 
  call d_cs1gd(0, y_up, wb_up, wcoef_up, w1_up)
  call d_cs1gd(0, y_up, tb_up, tcoef_up, t1_up)

  call d_cs1gd(1, y_up, ub_up, ucoef_up, u2_up) 
  call d_cs1gd(1, y_up, wb_up, wcoef_up, w2_up)
  call d_cs1gd(1, y_up, tb_up, tcoef_up, t2_up)
!---
! writing Ue, We ang angle 
  write(1, fmt='(2i10, 3E30.15)') i, 0, ue(i), we(i), angle(i)
! writing bottom part
  do j = 1,n_bt
    
	  write (1,fmt='(2i10, 12E30.15)') i,j,x(i),y_bt(j),  &
	  u_bt(j),u1_bt(j),u2_bt(j),                          &
	  w_bt(j),w1_bt(j),w2_bt(j),                          &
	  T_bt(j),T1_bt(j),T2_bt(j),                          &
	  r(j)        
  enddo
!writing upper part  
  do j = 1,n_up
    
	  write (1,fmt='(2i10, 12E30.15)') i,n_bt+j,x(i),y_up(j),  &
	  u_up(j),u1_up(j),u2_up(j),                               &
	  w_up(j),w1_up(j),w2_up(j),                               &
	  T_up(j),T1_up(j),T2_up(j),                               &
	  r(n_bt+j)        
  enddo
end do
close(123)
close(1)

end subroutine smooth_3D_profiles