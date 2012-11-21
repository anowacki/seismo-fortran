program F_cum
!===============================================================================
!  Outputs values of the cumulative F distribution F(x;m,n)
   use statistical
   implicit none
   integer,parameter :: rs = 8
   real(rs) :: x
   integer  :: m,n
   character(len=80) :: arg
   
   if (iargc() /= 3) then
      write(*,'(a)') 'Usage: F_cum [x] [n] [m]'
      stop
   endif
   
   call getarg(1,arg) ;  read(arg,*) x
   call getarg(2,arg) ;  read(arg,*) n
   call getarg(3,arg) ;  read(arg,*) m
   
   write(*,*) f_dist_cum(n,m,x)
   
end program F_cum