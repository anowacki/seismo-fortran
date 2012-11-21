program test

use spherical_splines

implicit none
real(8) :: h,d,h_ref
real(8),parameter :: tol = 0.001_8
integer :: iostatus

!  Read values from file with precomputed half-widths (radians) and values of h
open(10,file='half-width.text')

iostatus = 0
do while (iostatus == 0)
   read(10,*,iostat=iostatus) h_ref,d
   if (iostatus < 0) exit  ! EOF
   if (iostatus > 0) then
      write(0,'(a)') 'Problem reading file'
      stop
   endif
   h = sph_splines_AP_width2h(d)
   if (abs(h-h_ref) > tol) write(*,'(a,f0.5)') 'h and h_ref differ by more than ',tol
enddo

close(10)

end program test
