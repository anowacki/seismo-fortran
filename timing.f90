!===============================================================================
module timing
!===============================================================================
!  Subroutines for timing things inside a Fortran code.

   implicit none
   
   integer, parameter :: lu_stderr = 0, lu_stdout = 6
   real, private, save :: time1 = 0., time2 = 0.

contains

!===============================================================================
subroutine tick()
!===============================================================================
! Start the timer.
   call cpu_time(time1)
end subroutine tick
!-------------------------------------------------------------------------------

!===============================================================================
subroutine tock(out, message, fmt, lu, time)
!===============================================================================
! Stop the timer and print out the time unless out=.false.
! Optionally supply a format string to format the time with, and a message to
! precede the time.
! Optionally supply the logical unit number to write the output to.
! Optionally supply the elapsed time.
   logical, intent(in), optional :: out
   character(len=*), intent(in), optional :: message, fmt
   integer, intent(in), optional :: lu
   real, intent(out), optional :: time
   logical :: out_in
   character(len=250) :: fmt_in
   integer :: lu_in
   
   call cpu_time(time2)
   if (present(time)) time = time2 - time1
   out_in = .true.
   if (present(out)) out_in = out
   if (out_in) then
      fmt_in = '("Time elapased: ",f0.3)'
      if (present(fmt)) fmt_in = fmt
      lu_in = lu_stdout
      if (present(lu)) lu_in = lu
      if (present(message)) write(lu_in,'(a)') message
      write(lu_in,fmt_in) time2 - time1
   endif
end subroutine tock
!-------------------------------------------------------------------------------

end module timing
!-------------------------------------------------------------------------------
