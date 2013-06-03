subroutine smooth_3D_profiles(y, u, ny, u1, u2)
include 'link_f90_static.h'

use csscv_int
use cs1gd_int

implicit none

! from NS solver
integer, intent(in) :: ny
real(8), dimension(ny) :: y, u, u1, u2
integer n_bt, n_up, j, iequal
! for bottom splines                          
real(8), allocatable :: y_bt(:), u_bt(:), u1_bt(:), u2_bt(:), ub_bt(:), ucoef_bt(:,:)
! for upper splines
real(8), allocatable :: y_up(:), u_up(:), u1_up(:), u2_up(:), ub_up(:), ucoef_up(:,:)
!------------------------------------------------------
n_bt = ny
!n_up = ny - n_bt
allocate (y_bt(n_bt), u_bt(n_bt), u1_bt(n_bt), u2_bt(n_bt), ub_bt(n_bt))
allocate (ucoef_bt(4,n_bt))

!allocate (y_up(n_up), u_up(n_up), u1_up(n_up), u2_up(n_up), ub_up(n_up))
!allocate (ucoef_up(4,n_up))

  u_bt(1:n_bt) = u(1:n_bt)
  y_bt(1:n_bt) = y(1:n_bt)
  
!  u_up(1:n_up) = u(n_bt+1:ny)
!  y_up(1:n_up) = y(n_bt+1:ny)

iequal = 0

!bottom splines
  call d_csscv(y_bt, u_bt, iequal, ub_bt, ucoef_bt) 
  call d_cs1gd(1, y_bt, ub_bt, ucoef_bt, u1_bt) 
  call d_csscv(y_bt, u1_bt, iequal, ub_bt, ucoef_bt)
  call d_cs1gd(0, y_bt, ub_bt, ucoef_bt, u1_bt) 
  call d_cs1gd(1, y_bt, ub_bt, ucoef_bt, u2_bt) 
  
!upper slice
!  call d_csscv(y_up, u_up, iequal, ub_up, ucoef_up)  
!  call d_cs1gd(1, y_up, ub_up, ucoef_up, u1_up) 
!  call d_csscv(y_up, u1_up, iequal, ub_up, ucoef_up)
!  call d_cs1gd(0, y_up, ub_up, ucoef_up, u1_up) 
!  call d_cs1gd(1, y_up, ub_up, ucoef_up, u2_up) 
!---
u1(1:n_bt) = u1_bt
!u1(n_bt+1:ny) = u1_up
u2(1:n_bt) = u2_bt
!u2(n_bt+1:ny) = u2_up

deallocate (u_bt, u1_bt, u2_bt, ub_bt, ucoef_bt) 
!deallocate (u_up, u1_up, u2_up, ub_up, ucoef_up)
end subroutine smooth_3D_profiles