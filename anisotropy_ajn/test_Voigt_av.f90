program test_Voigt_av
!===============================================================================
!  Test the CIJ_Voigt_av subroutine

use anisotropy_ajn

implicit none

real(8) :: C(10,6,6),rh(10),VF(10)
real(8) :: CVoigt(6,6),rhVoigt
integer :: i

!  Set density to be the same for all values
rh = 5500.

!  Set equal weighting for all values--this will be normalised in the subroutine.
VF = 1.

!  Create 10 arrays filled with 1,2,3...10.  The Voigt average is an array filled
!  with the mean = (1+2+3+...+9+10)/10 = 5.5
do i=1,10
   C(i,:,:) = real(i)
enddo

call CIJ_Voigt_av(VF,C,rh,CVoigt,rhVoigt)

write(*,'(a)') '====== Testing CIJ_Voigt_av ======'
write(*,'(a)') 'Expected output:'
write(*,'(6(f4.1))') (5.5,i=1,36)
write(*,'(a)') ''
write(*,'(a)') 'Voigt average:'
write(*,'(6(f4.1))') CVoigt
write(*,'(a)') '=================================='

end program