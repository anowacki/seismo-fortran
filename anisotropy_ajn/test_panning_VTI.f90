program test
! Test the CIJ_pnning_VTI function

use anisotropy_ajn

implicit none

real(8) :: C(6,6),vp,vs,rho,xi,phi

! At 2750 in AK135
vp = 13650.
vs = 7249.
rho = 5423.
xi = 1.3
phi = 0.8

C = CIJ_panning_VTI(vp,vs,rho,xi,phi)

write(*,*) C

end program test
