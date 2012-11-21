use statistical

implicit none

integer,parameter :: n=8
integer :: i
real(8) :: a(n),mu,mu1,mu2

do i=1,n
   a(i) = i
enddo

call circ_mean_bootstrap_smallN(a,mu,mu1,mu2)

write(*,*) mu,mu1,mu2

end