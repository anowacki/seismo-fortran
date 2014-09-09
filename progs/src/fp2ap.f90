!===============================================================================
! One of Andy Nowacki's Fortran utility modules for dealing with seismic
! anisotropy and other problems.
!
! Andy Nowacki <andy.nowacki@bristol.ac.uk>
!
! See the file LICENCE for licence details.
!===============================================================================
program fp2ap
!===============================================================================
!  Calculate the strike, dip and rake of the auxiliary plane, given the fualt plane
!  for a double-couple event.

   use moment_tensor, only: mt_faultplane2auxplane

   implicit none

   integer,parameter :: rs = 8
   character(80) :: arg
   real(rs) :: s1,d1,k1  ! Strike, dip, rake of first plane
   real(rs) :: s2,d2,k2  !                      second plane

   if (command_argument_count() /= 3) then
      write(0,'(a)') 'Usage: fp2ap [strike] [dip] [rake]'
      stop
   endif

!  Get parameters of first plane
   call get_command_argument(1, arg) ; read(arg,*) s1
   call get_command_argument(2, arg) ; read(arg,*) d1
   call get_command_argument(3, arg) ; read(arg,*) k1

   call mt_faultplane2auxplane(s1, d1, k1, s2, d2, k2)

   write(*,'(f6.1,f5.1,f7.1)') s2, d2, k2

end program fp2ap
!-------------------------------------------------------------------------------
