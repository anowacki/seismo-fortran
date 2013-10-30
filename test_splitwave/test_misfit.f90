! Test the misfit function
program test
   use splitwave
   real :: phi1 = 10., phi2 = 30., dt1 = 1., dt2 = 1.2
   real :: spol = 0.
   real :: freq = 0.1
   real :: noise = 0.
   character(len=1) :: wavetype = 'g'
   
   call get_command_line
   
   write(*,*) sw_misfit(phi1,dt1, phi2,dt2, spol, freq=freq, noise=noise, wavetype=wavetype)

contains
   subroutine get_command_line()
      implicit none
      character(len=250) :: arg
      integer :: iostat
      if (command_argument_count() /= 0 .and. command_argument_count() < 4) &
         call usage
      if (command_argument_count() == 4) then
         call get_command_argument(1,arg)
         read(arg,*,iostat=iostat) phi1
         if (iostat /= 0) call usage
         call get_command_argument(2,arg)
         read(arg,*,iostat=iostat) dt1
         if (iostat /= 0) call usage
         call get_command_argument(3,arg)
         read(arg,*,iostat=iostat) phi2
         if (iostat /= 0) call usage
         call get_command_argument(4,arg)
         read(arg,*,iostat=iostat) dt2
         if (iostat /= 0) call usage
      endif
   end subroutine get_command_line
   
   subroutine usage()
      implicit none
      write(0,'(a)') 'Usage: test_misfit ([phi1] [dt1] [phi2] [dt2])'
      stop
   end subroutine usage
end program test
