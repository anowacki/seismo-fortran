!===============================================================================
! One of Andy Nowacki's Fortran utility modules for dealing with seismic
! anisotropy and other problems.
!
! Andy Nowacki <andy.nowacki@bristol.ac.uk>
!
! See the file LICENCE for licence details.
!===============================================================================
program circ_test_orient
!===============================================================================
! Reading a list of orientational data, check whether they are statistically
! significantly aligned with a particular direction.  Give the chosen orientation
! to test normally; it is converted to orientational data (doubled) internally.
   use statistical

   implicit none

   integer,  parameter :: nmax = 10000
   real(8) :: angles(nmax), test_angle, alpha = 0.05
   integer :: iostat, n
   character(len=250) :: arg

   ! Get command line parameters
   if (command_argument_count() < 1) call usage
   call get_command_argument(1, arg)
   if (arg(1:2) == "-h") call usage
   read(arg,*,iostat=iostat) test_angle
   if (iostat /= 0) then
      write(0,'(a)') 'circ_test_orient: Error: Cannot get value of alpha from ' &
         //'argument "'//trim(arg)//'"'
      stop
   endif
   if (command_argument_count() == 2) then
      call get_command_argument(2, arg)
      read(arg,*,iostat=iostat) alpha
      if (iostat /= 0) then
         write(0,'(a)') 'circ_test_orient: Error: Cannot get value of alpha from ' &
            //'argument "'//trim(arg)//'"'
         stop
      endif
   endif

   n = 1
   do
      read(*,*,iostat=iostat) angles(n)
      if (iostat < 0) then
         n = n - 1
         exit
      endif
      if (iostat > 0) then
         write(0,'(a)') 'circ_test_orient: Error: Problem reading angle from stdin.'
         stop
      endif
      n = n + 1
   enddo

   if (circ_test_random_orient(angles(1:n),test_angle,alpha=alpha,degrees=.true.)) then
      write(*,'(a,f0.3,a)') 'Angles differ significantly from random towards ', &
         test_angle,' degrees: Reject H0 (random)'
   else
      write(*,'(a,f0.3,a)') 'Angles do not differ significantly from random towards ', &
         test_angle,' degrees: Accept H0 (random)'
   endif

contains
   !============================================================================
   subroutine usage
   !============================================================================
   write(0,'(a)') &
      'Usage: circ_test_orient [test angle] (alpha) < [orientations in degrees]', &
      'Tests whether a set of angles (on stdin) significantly differ from random', &
      'towards the test angle, in the 5% limit.', &
      'Optionally set the significance level by calling with alpha = 0.1|0.05|0.01|0.005|0.001|0.0001'
   stop
   end subroutine usage
   !----------------------------------------------------------------------------
end program circ_test_orient
!-------------------------------------------------------------------------------
