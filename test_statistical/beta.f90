program beta
use statistical
implicit none

integer,parameter :: rs=8
real(rs) :: p,q,dp,dq

dp = 0.1_rs
dq = 0.1_rs

! Set up p and q so they are never integers
p = -3._rs + dp/2._rs
do while (p < 3._rs)
   q = -3._rs + dq/2._rs
   do while (q < 3._rs)
      write(*,*) p,q,beta_func(p,q)
      q = q + dq
   enddo
   p = p + dp
   write(*,*)  ! For gnuplot
enddo

end program beta
