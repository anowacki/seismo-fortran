program test_disp

use anisotropy_ajn, only: CIJ_disp

real(8) :: C(6,6)
integer :: i,j

do i=1,6
   do j=1,6
      C(i,j) = real(10*i + j)*10.**i
   enddo
enddo

write(*,*) 'Auto:'
call CIJ_disp(C)

write(*,*)
write(*,*) '1 dp, 10^6:'
call CIJ_disp(C, power=6, ndp=1)

write(*,*)
write(*,*) '5 dp, sci, stdout'
call CIJ_disp(C, ndp=5, expo=.true., unit=0)


end program test_disp