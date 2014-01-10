program f_dist_prog
!===============================================================================
!  Print out the f-distribution, for plotting
   use statistical
   implicit none
   integer,parameter :: rs=8
   integer :: m,n
   real(rs) :: x,dx,x1
   character(len=25) :: arg
   
   if (command_argument_count() /= 3) then
      write(*,'(a)') 'Usage: f_dist [n] [m] [upper value]'
      stop
   endif
   call get_command_argument(1,arg) ;  read(arg,*) n
   call get_command_argument(2,arg) ;  read(arg,*) m
   call get_command_argument(3,arg) ;  read(arg,*) x1
   
   x  = 0._rs
   dx = 0.01_rs
   
   do while (x <= x1)
      x = x + dx
      write(*,*) x, f_dist(n,m,x)
   enddo
end program f_dist_prog
