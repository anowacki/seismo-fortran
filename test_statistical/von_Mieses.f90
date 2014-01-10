program von_Mieses
!  Print out the von Mieses distribution over a circle for a given kappa
!  If not given on command line, kappa assumed=1

use statistical

real(8) :: theta,dtheta,kappa,mu
character(len=80) :: arg

mu = 0.
kappa = 1.

if (command_argument_count() == 1) then
   call get_command_argument(1,arg)
   read(arg,*) kappa
   if (kappa < 0.) then
      write(0,'(a)') 'von_Mieses: kappa must be >= 0'
      stop
   endif
endif

theta = 0.
dtheta = 1
do while (theta <= 360.)
   write(*,*) theta,circ_von_mieses(theta,mu,kappa)
   theta = theta + dtheta
enddo

end program