! helper for gradient methods
! FOR TIME SEARCH ONLY
module HELPER
contains 
function DFSI_DX(ind)
implicit none
integer, intent(in):: ind
real*8 DFSI_DX
common /HADY2/ Z(4)                               &
       /HADY3/ SI,ZZZ0                       
complex*16 Z, SI, ZZZ0, z_def, si_def
real*8 dx,  fsi_rgt, fsi_lft
z_def = Z(ind)
si_def = SI
!dx = 0.0001*dabs(dreal(Z(ind)))
!if (dx<1.0d-8) dx = 1.0d-8
dx=1.0d-8
Z(ind)=z_def+dx
call HADSI
fsi_rgt=cdabs(SI)

Z(ind)=z_def-dx
call HADSI
fsi_lft=cdabs(SI)

DFSI_DX = (fsi_rgt - fsi_lft)/(2.0*dx)
Z(ind) = z_def
SI = si_def
end function DFSI_DX

subroutine NEWTON_STEP()
implicit none
common /HADY2/ Z(4)                               &
       /HADY3/ SI,ZZZ0                       
complex*16 Z, SI, ZZZ0, a_def, b_def
real*8 df_da, df_db, fsi_def, limiter
call HADSI
fsi_def = cdabs(SI)
df_da = DFSI_DX(1)
df_db = DFSI_DX(2)
limiter = 1.0
a_def = Z(1)
b_def = Z(2)
!do while (cdabs(SI)>=fsi_def)
Z(1) = a_def - limiter*cdabs(SI)/df_da
!Z(1) = a_def - limiter*cdabs(SI)*df_da/(df_da**2.0+df_db**2.0)
!Z(2) = b_def - limiter*cdabs(SI)*df_db/(df_da**2.0+df_db**2.0)
call HADSI
print*, "A:",dreal(Z(1))
!limiter = limiter/2.0
!end do
end subroutine NEWTON_STEP
end module HELPER