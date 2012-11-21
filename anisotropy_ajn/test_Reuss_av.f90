program test_Reuss_av
!===============================================================================
!  Test the CIJ_Reuss_av subroutine

use anisotropy_ajn

implicit none

real(8) :: C(5,6,6),rh(5),VF(5)
real(8) :: CReuss(6,6),rhReuss
real(8) :: CTrue(6,6)
integer :: n,ii,jj

!  Set density to be the same for all values
rh = 5500.

!  Set equal weighting for all values--this will be normalised in the subroutine.
VF = 1.

!  Create 5 arrays filled with:
!     1     1     1     1     1     1
!     1     2     2     2     2     2
!     1     2     3     3     3     3   x   {1..5}
!     1     2     3     4     4     4
!     1     2     3     4     5     5
!     1     2     3     4     5     6
!
!  The Reuss average is an array filled with the harmonic mean of each element:
CTrue = reshape( &
    [ 2.1898,   2.1898,   2.1898,   2.1898,   2.1898,   2.1898, &
      2.1898,   4.3796,   4.3796,   4.3796,   4.3796,   4.3796, &
      2.1898,   4.3796,   6.5693,   6.5693,   6.5693,   6.5693, &
      2.1898,   4.3796,   6.5693,   8.7591,   8.7591,   8.7591, &
      2.1898,   4.3796,   6.5693,   8.7591,  10.9489,  10.9489, &
      2.1898,   4.3796,   6.5693,   8.7591,  10.9489,  13.1387  ], [6,6])
!  If no F2008-compatible compiler around, change the [ ] for (/ /).

do n=1,5
   do ii=1,6
      do jj=ii,6
         C(n,ii,jj) = real(n*ii)
         C(n,jj,ii) = C(n,ii,jj)
      enddo
   enddo
enddo

call CIJ_Reuss_av(VF,C,rh,CReuss,rhReuss)

write(*,'(a)') '====== Testing CIJ_Reuss_av ======'
write(*,'(a)') 'Expected output:'
write(*,'(6(f5.1))') CTrue
write(*,'(a)') ''
write(*,'(a)') 'Reuss average:'
write(*,'(6(f5.1))') CReuss
write(*,'(a)') '=================================='

end program