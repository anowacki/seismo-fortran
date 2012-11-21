program test
!  Test the CIJ_brow_chev subroutine

use anisotropy_ajn

real(8),dimension(6,6) :: C,CI,CX,CT,CO,CM,CR
real(8) :: pI,pX,pT,pO,pM,pR
real(8) :: rho
character(len=250) :: fname,fmt1,fmt2

if (iargc() /= 1) then
   write(0,'(a)') 'Usage: test_CIJ_brow_chev [.ecs file]'
   stop
endif

call getarg(1,fname)

!  Load constants
call CIJ_load(fname,C,rho)

!  Convert into GPa
C = C / 1.d9

!  Express constants according to B&C 2004
call CIJ_brow_chev(C,CI,CX,CT,CO,CM,CR) !,pI,pX,pT,pO,pM,pR)

!  Output format
fmt1 = '(a)' !,f6.2," %")'
fmt2 = '(6e12.3)'
fmt2 = '(6f7.1)'

write(*,'(a)') 'Input tensor:'
write(*,fmt2) C
write(*,*)

write(*,fmt1) 'Isotropic part:    ' !,pI*100._8
write(*,fmt2) CI
write(*,*)

write(*,fmt1) 'Hexagonal part:    ' !,pX*100._8
write(*,fmt2) CX
write(*,*)

write(*,fmt1) 'Tetragonal part:   ' !,pT*100._8
write(*,fmt2) CT
write(*,*)

write(*,fmt1) 'Orthorhombic part: ' !,pO*100._8
write(*,fmt2) CO
write(*,*)

write(*,fmt1) 'Monoclinic part:   ' !,pM*100._8
write(*,fmt2) CM
write(*,*)

write(*,fmt1) 'Triclinic part:    ' !,pR*100._8
write(*,fmt2) CR
write(*,*)

write(*,fmt1) 'Sum of all parts:  ' !,(pI + pX + pT + pO + pM + pR)*100._8
write(*,fmt2) CI + CX + CT + CO + CM + CR
write(*,*)


end program test