subroutine MACK_GLOBAL_SPATIAL
IMPLICIT none
COMMON /HADY2/ Z(4)                     &
       /HADY3/ SI,ZZZ0                  & 
       /OUT  /AAK,AASIG,AAG,AAM,AAXI,AAPB
COMPLEX*16 SI,Z,ZZZ0
real(8) AAK,AASIG,AAG,AAM,AAXI,AAPB
!</globals>
real(8) w_low, w_up, w, dw,             &
        ar_low, ar_up, ar, d_ar,        &
        cf_low, cf_up,                  &
        ai, ai_min, ar_min, w_min,      &
        fsi, fsi_min
integer n_w, n_ar,                      & 
        j_w, j_ar         
!------------ Global search for Mack second mode -------
Z(2)=(0.,0.)	        ! wavenumber beta
w_low=0.1	            ! initial freuency
w_up =0.6               ! final frequncy
ai=-5.E-3	            ! growth rate
cf_up=1.		        ! maximal phase speed
cf_low=1.-1./AAM	    ! minimal phase speed
n_w=101
n_ar=101
dw=(w_up-w_low)/(dble(n_w)-1.0)
fsi_min=1.E+9
DO j_w=1,n_w
  w=w_low+dw*(dble(j_w)-1.0)
  ar_low=w/cf_up
  ar_up=w/cf_low
  d_ar=(ar_up-ar_low)/(dble(n_ar)-1.0)
  WRITE(2,'(/)')
  DO j_ar=2,n_ar-1
    ar=ar_low+d_ar*(dble(j_ar)-1.0)
	Z(3)=w
    Z(1)=DCMPLX(ar,ai)
    CALL HADSI
    fsi=CDABS(SI)
    IF(fsi.LT.fsi_min)THEN
      fsi_min=fsi
      ar_min=ar
      ai_min=ai
      w_min=w
    ENDIF
  ENDDO
ENDDO
Z(1)=DCMPLX(ar_min,ai_min)
Z(3)=w_min
write(*,'("Mack Glob Srch Res:",/,      &
          "A: ",2E13.5,/,               &
          "B: ",2E13.5,/,               &
          "W: ",2E13.5 )')              &
      Z(1), Z(2), Z(3) 
end subroutine MACK_GLOBAL_SPATIAL
!#############################################
subroutine MACK_GLOBAL_TIME
IMPLICIT none
COMMON /HADY2/ Z(4)                     &
       /HADY3/ SI,ZZZ0                  & 
       /OUT  /AAK,AASIG,AAG,AAM,AAXI,AAPB
COMPLEX*16 SI,Z,ZZZ0, w, a_min, w_min
real(8) AAK,AASIG,AAG,AAM,AAXI,AAPB

real(8) wr_low, wr_up, dw_r,            &
        a_low, a_up, a, da, wi,         &
        cf_low, cf_up,                  &
        fsi, fsi_min
integer nw, na,                         & 
        jw, ja         
!------------ Global search for Mack second mode -------
Z(2)=(0.,0.)	        ! wavenumber beta
wr_low=0.1	            ! initial freuency
wr_up =1.0               ! final frequncy
wi=1.E-3	            ! growth rate
cf_up=1.		        ! maximal phase speed
cf_low=1.-1./AAM	    ! minimal phase speed
nw=101
na=101
dw_r=(wr_up-wr_low)/(dble(nw)-1.0)
fsi_min=1.E+9
DO jw=1,nw
  w=DCMPLX(wr_low+dw_r*(dble(jw)-1.0), wi)
  a_low=DREAL(w)/cf_up
  a_up= DREAL(w)/cf_low
  da=(a_up-a_low)/(dble(na)-1.0)
  DO ja=2,na-1
    a=a_low+da*(dble(ja)-1.0)
	Z(3)=w
    Z(1)=DCMPLX(a,0)
    CALL HADSI
    fsi=CDABS(SI)
    IF(fsi.LT.fsi_min)THEN
      print *,"GlobSrch FSI:", fsi
      fsi_min=fsi
      a_min=Z(1)
      w_min=Z(3)
    ENDIF
  ENDDO
ENDDO
Z(1)=a_min
Z(3)=w_min
write(*,'("Mack Time Glob Srch Res:",/, &
          "A: ",2E13.5,/,               &
          "B: ",2E13.5,/,               &
          "W: ",2E13.5 )')              &
      Z(1), Z(2), Z(3) 
end subroutine MACK_GLOBAL_TIME
!###########################################################
subroutine TS_GLOBAL_TIME()
implicit none
common /HADY2/ Z(4)                         &
       /HADY3/ SI,ZZZ0                      &
       /OUT  / AAK,AASIG,AAG,AAM,AAXI,AAPB  
complex*16 SI,Z,ZZZ0 
real(8) AAK,AASIG,AAG,AAM,AAXI,AAPB
!</globals>
real(8), parameter :: pi = 3.141592653
real(8), parameter :: theta = pi/3.0         
real(8), parameter :: wi = 2.0d-5 
real(8), parameter :: wr_start = 0.005
real(8), parameter :: wr_end   = 0.1
integer, parameter :: nw = 100, na = 100
real(8)  wr_cur, cfup, cflow,                      &
         a_start, a_end,                           &
         fsi_cur, fsi_min,                         &
         a_res, b_res, wr_res
integer i,j
CFUP = 1.0                        ! maximum phase speed
CFLOW = 1.-1./(AAM*cos(theta))	! minimum phase speed        
IF (CFLOW<0.0) CFLOW = 0.1
fsi_min = 1.0d+09
do i=1, nw
  wr_cur = wr_start + (wr_end - wr_start)/dble(nw)*dble(i)
  a_start = wr_cur/cfup
  a_end = wr_cur/cflow
  do j = 1, na
    Z(1) = a_start + (a_end - a_start)/dble(na)*dble(j)
    Z(2) = Z(1)*tan(theta)
    Z(3) = dcmplx(wr_cur, wi)
    call LHADY
    fsi_cur = cdabs(SI)
    if (fsi_cur<fsi_min) then
      print *,"GlobSrch FSI:", fsi_cur
      fsi_min = fsi_cur
      a_res = DREAL(Z(1))
      b_res = DREAL(Z(2))
      wr_res = DREAL(Z(3))
    end if
  end do
end do
Z(1) = dcmplx(a_res,0)
Z(2) = dcmplx(b_res,0)
Z(3) = dcmplx(wr_res, wi)
call POISK2(3)
end subroutine TS_GLOBAL_TIME
!###########################################################
subroutine CF_GLOBAL_SEARCH
implicit none
common /HADY2/ A, B, W, R                    &
       /HADY3/ SI,ZZZ0                       &
       /CF_PROFILE_DATA/ MF_I, MF_J,         &
                         MF_X, MF_WE, MF_ANG 
complex*16 A, B, W, R,                &
           SI, ZZZ0,                  &
           VA,VB,VR                  
integer MF_I, MF_J, i, j, k
real(8) MF_X, MF_WE, MF_ANG
integer, parameter :: na = 100,       &
                      nb = 100,       &
                      nw = 100
                      
real(8) psi
real(8) alpha, beta,                             &
        alpha_start, alpha_end, alpha_min,       &
        beta_start, beta_end, beta_min,          &
        da, db, fsi_min ,psi1        
complex*16 omega, dw,                            & 
           omega_start, omega_end, omega_min
! define coord frame rotation angle 
!(read from input NS profiles datafile)
psi = -MF_ANG   ! must be generalized
beta_start = 1.0d-3
beta_end = 5.0d-1
omega_start = DCMPLX(0.00d0, 0.0)
omega_end = DCMPLX(0.00d0, 0.02)
dw = (omega_end - omega_start)/dble(nw)
db = (beta_end - beta_start)/dble(nb)
fsi_min = 1.0d+9
do j = 1, nb
  do k = 1, nw
    B = DCMPLX(beta_start + dble(j)*db, 0.0d0)
    A = tan(psi)*B
    W = omega_start + dble(k)*dw
    call HADSI
    if (CDABS(SI)<fsi_min) then
      print *, "current fsi_min:", fsi_min
      write(*,'("A=",2E13.5,/,"B=",2E13.5,/,"W=",2E13.5,/)') A, B, W
      fsi_min = CDABS(SI)
      alpha_min = DREAL(A)
      beta_min = DREAL(B)
      omega_min = W
    end if
  end do
end do
A = alpha_min
B = beta_min
W = omega_min
call POISK2(3)
print '("CF Time Global Search Result:",/,                &
        "A=",2E13.5,/,"B=",2E13.5,/,"W=",2E13.5,/)',      &
        A, B, W  
end subroutine CF_GLOBAL_SEARCH