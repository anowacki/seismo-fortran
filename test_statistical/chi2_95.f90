!===============================================================================
program chi2_95
!===============================================================================
!  Compute percentage points of the chi-2 probability distribution
   use statistical
   implicit none
   
   integer,parameter :: rs=8
   real(rs) :: p,v,g
   character(20) :: arg
   
   if (command_argument_count() /= 1) then
      write(*,'(a)') 'Usage: chi2_95 [degrees of freedom]'
      stop
   endif
   
   call get_command_argument(1,arg) ; read(arg,*) v
   p = 0.05
   g = log_gamma(v/2._rs)
   
   write(*,*) ppchi2(p,v,g)

end program chi2_95
