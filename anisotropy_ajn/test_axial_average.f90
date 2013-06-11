!===============================================================================
program test
!===============================================================================
   
   use anisotropy_ajn
   
   implicit none
   
   integer, parameter :: rs = 8
   real(rs) :: C(6,6), Cave(6,6), rho
   
   ! Constants of San Carlos olivine, Abramson et al, 1997
   C(1,1) = 320.50e9
   C(1,2) = 68.10e9
   C(1,3) = 71.60e9
   C(1,4) = 0.00e9
   C(1,5) = 0.00e9
   C(1,6) = 0.00e9
   C(2,1) = 68.10e9
   C(2,2) = 196.50e9
   C(2,3) = 76.80e9
   C(2,4) = 0.00e9
   C(2,5) = 0.00e9
   C(2,6) = 0.00e9
   C(3,1) = 71.60e9
   C(3,2) = 76.80e9
   C(3,3) = 233.50e9
   C(3,4) = 0.00e9
   C(3,5) = 0.00e9
   C(3,6) = 0.00e9
   C(4,1) = 0.00e9
   C(4,2) = 0.00e9
   C(4,3) = 0.00e9
   C(4,4) = 64.00e9
   C(4,5) = 0.00e9
   C(4,6) = 0.00e9
   C(5,1) = 0.00e9
   C(5,2) = 0.00e9
   C(5,3) = 0.00e9
   C(5,4) = 0.00e9
   C(5,5) = 77.00e9
   C(5,6) = 0.00e9
   C(6,1) = 0.00e9
   C(6,2) = 0.00e9
   C(6,3) = 0.00e9
   C(6,4) = 0.00e9
   C(6,5) = 0.00e9
   C(6,6) = 78.70e9
   rho = 3355.e0
   
   ! Rotate a little bit to remove the perfect match when using few rotations
   call CIJ_rot3(C,13._rs,194._rs,38._rs,Cave)
   C = Cave
   
   ! Rotate about the 3-axis
   call CIJ_axial_average(C,3,Cave,nrot=10,ave_type='VRH')
   
   ! Write out constants density normalised
   write(*,*) Cave/rho

end program test
!-------------------------------------------------------------------------------