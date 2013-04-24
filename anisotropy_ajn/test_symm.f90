program test_CIJ_symm

use anisotropy_ajn, only: CIJ_symm

implicit none

real(8) :: C(6,6) = reshape((/11, 12, 13, 14, 15, 16, &
                               0, 22, 23, 24, 25, 26, &
                               0,  0, 33, 34, 35, 36, &
                               0,  0,  0, 44, 45, 46, &
                               0,  0,  0,  0, 55, 56, &
                               0,  0,  0,  0,  0, 66/), (/6,6/))

write(*,'(a)') 'Matrix before symmetrisation:'
write(*,'(6(f4.0,1x))') C

call CIJ_symm(C)

write(*,'(/a)') 'Matrix after symmetrisation:'
write(*,'(6(f4.0,1x))') C

end program test_CIJ_symm
