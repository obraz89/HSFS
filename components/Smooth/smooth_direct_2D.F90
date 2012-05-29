! need to add input/ to all input pathes
subroutine process_direct

include 'link_f90_static.h'

use csscv_int
use cs1gd_int

implicit none
real(8), parameter :: Re=1.424E+07, Ma=3.5, gamma=1.4
real(8), parameter :: L=0.381, Tinf=90.318
real(8), parameter :: pi = 3.1415926535897932384626433832795d0
real(8) :: alfa=5.d0*pi/180.d0

integer, parameter :: nx=201, ny=201
character file_name*100, nx_str*(3), ny_str*(3)

real(8) xx(nx,ny),y(nx,ny),u(nx,ny),v(nx,ny),r(nx,ny),T(nx,ny), p(nx,ny)
real(8) x
real(8) ut(ny),ut1(ny),ut2(ny), T1(ny),T2(ny)
real(8) yd(ny), b1(ny),b2(ny),utc(4,ny),Tc(4,ny)


character(3) strn, xxx
character(4) unam,vnam,rnam,Tnam
integer i,j,k, iequal, i_spher
!------------------------------------------------------
write(nx_str, '(i3)') nx
write(ny_str, '(i3)') ny

OPEN(123, FILE='input/KING_2D.dat', FORM='BINARY',STATUS='OLD')
do j=1,ny
    do i=1,nx
        READ(123) xx(i,j),y(i,j)  ! сетка
    end do
end do

do j=1,ny
    do i=1,nx
        READ(123) u(i,j)
        READ(123) v(i,j)
        READ(123) p(i,j)
        READ(123) t(i,j)
    end do
end do
close(123)
!конец считывания

r=gamma*Ma*Ma*p(1:nx,1:ny)/t(1:nx,1:ny)


!NS mean flow profiles in the reverse order required by the 
!Stability Solver

file_name = 'output/KING_SSD_2D_'//nx_str//'x'//ny_str//'.dat'

open (1,file=file_name)
!write (1,fmt='(2a10, 9a30)') 'i','j','x','y',       &
!                             'u','u1','u2',         &
!							 'T','T1','T2','rho'
do i = nx-10,2,-1
	ut = u(i,1:ny)*dcos(alfa)+v(i,1:ny)*dsin(alfa) !профиль (скор вдоль конуса)
	
	x = xx(i,1)
	yd = (y(i,1:ny)-y(i,1))/cos(alfa)*dsqrt(Re) !нормировка на погранслой

  iequal = 0

  call d_csscv (yd, ut, iequal, b1, utc) !коэффициенты функции
  call d_csscv (yd, T(i,1:ny) , iequal, b2,  Tc)

  call d_cs1gd (1, yd, b1, utc, ut1) !первая производная
  call d_cs1gd (1, yd, b2, Tc ,  T1)

  call d_csscv (yd, ut1, iequal, b1, utc) !коэффициенты первой производной
  call d_csscv (yd, T1 , iequal, b2, Tc)
  
  call d_cs1gd (0, yd, b1, utc, ut1) !сглаженная первая производная
  call d_cs1gd (0, yd, b2, Tc ,  T1)

  call d_cs1gd (1, yd, b1, utc, ut2) !вторая производная
  call d_cs1gd (1, yd, b2, Tc , T2)
  write (1,fmt='(2i5, 3E30.15)') 0 ,0 ,0., 0., 0.
  do j = 1,ny
    
	  write (1,fmt='(2i10, 12E30.15)') nx-10-i+1,j,x,yd(j),  &
	                                  ut(j),ut1(j),ut2(j),  &
	                                  0.0, 0.0, 0.0,        &
	                                  T(i,j),T1(j),T2(j),   &
	                                  r(i,j)
!	  write (1,fmt='(4E30.15)') yd(j),ut(j),ut1(j),ut2(j)

  enddo

enddo
close(1) 
 
end subroutine process_direct