program test_read_string

use anisotropy_ajn, only: CIJ_disp, CIJ_read_string

implicit none

double precision :: C(6,6), rho
logical :: found_rho

write(*,*) "=== 36 ECs and rho ==="
call CIJ_read_string("11 12 13 14 15 16 " // &
                     "21 22 23 24 25 26 " // &
                     "31 32 33 34 35 36 " // &
                     "41 42 43 44 45 46 " // &
                     "51 52 53 54 55 56 " // &
                     "61 62 63 64 65 66 " // &
                     "1000", C, rho)
call CIJ_disp(C, power=0)
write(*,*) rho

write(*,*)
write(*,*) "=== 36 ECs (no rho) ==="
call CIJ_read_string("11 12 13 14 15 16 " // &
                     "21 22 23 24 25 26 " // &
                     "31 32 33 34 35 36 " // &
                     "41 42 43 44 45 46 " // &
                     "51 52 53 54 55 56 " // &
                     "61 62 63 64 65 66 ", C, rho, found_rho=found_rho)
call CIJ_disp(C, power=0)
write(*,*) rho, found_rho

write(*,*)
write(*,*) "=== 21 ECs and rho ==="
call CIJ_read_string("11 12 13 14 15 16 " // &
                     "   22 23 24 25 26 " // &
                     "      33 34 35 36 " // &
                     "         44 45 46 " // &
                     "            55 56 " // &
                     "               66 " // &
                     "1000", C, rho, found_rho=found_rho)
call CIJ_disp(C, power=0)
write(*,*) rho, found_rho

write(*,*)
write(*,*) "=== 21 ECs and rho ==="
call CIJ_read_string("11 12 13 14 15 16 " // &
                     "   22 23 24 25 26 " // &
                     "      33 34 35 36 " // &
                     "         44 45 46 " // &
                     "            55 56 " // &
                     "               66 ", C, rho)
call CIJ_disp(C, power=0)
write(*,*) rho

end program test_read_string
