program f
! Return the value of the f-function
use statistical
implicit none
real(8) :: x
integer :: m,n
character(len=80) :: arg

if (iargc() /= 3) then
   write(0,'(a)') 'Usage: f [x] [n] [m]'
   stop
endif

call getarg(1,arg) ; read(arg,*) x
call getarg(2,arg) ; read(arg,*) n
call getarg(3,arg) ; read(arg,*) m

write(*,*) f_dist(n,m,x)

end program f
