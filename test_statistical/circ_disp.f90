program cd
!  Test the circ_disp function in statistical

use statistical

integer,parameter :: n=10
real(8) :: a(n)
integer :: i

!  Generate random angles in range 0 - 360
do i=1,n
   a(i) = rand()*360.
enddo

write(*,*) 'd of random angles:',circ_dispers(a)

!  Generate angles normally distributed about 0
do i=1,n
   a(i) = exp(-real(i-n)**2)*360.
enddo

write(*,*) 'd of ''normal'' angles:',circ_dispers(a)


end program
