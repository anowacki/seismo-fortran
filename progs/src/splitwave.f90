!===============================================================================
program splitwave
!===============================================================================

   use f90sac
   use splitwave

   implicit none
   
   type(SACtrace) :: t1, t2, t3
   real :: theta, dt
   character(len=250) :: f, f1, f2, f3
   integer :: iarg
   
   ! Check arguments
   if (command_argument_count() /= 3 .and. command_argument_count() /= 2) then
      write(0,'(a)') 'Usage:  splitwave [filename] [phi] [dt]'
      write(0,'(a)') '    or  splitwave [phi] [dt] (operates on ''wave.BH?'')'
      stop
   endif
   
   ! Get filename base
   if (command_argument_count() == 3) then
      call get_command_argument(1,f)
      iarg = 1
   else
      f = 'wave'
      iarg = 0
   endif
   
   ! Filenames
   f1 = trim(f) // '.BHN'
   f2 = trim(f) // '.BHE'
   f3 = trim(f) // '.BHZ'
   
   ! Get splitting parameters
   call get_command_argument(iarg+1,f)
   read(f,*) theta
   call f90sac_unwind(theta)
   
   call get_command_argument(iarg+2,f)
   read(f,*) dt
   
   ! Load traces
   call f90sac_readtrace(f1,t1)
   call f90sac_readtrace(f2,t2)
   call f90sac_readtrace(f3,t3)

   ! Apply splitting parameter
   call sw_splitN(t1, t2, t3, 1, (/theta/), (/dt/), quiet=.true.)

   ! Write out altered traces
   call f90sac_writetrace(f1,t1)
   call f90sac_writetrace(f2,t2)
   call f90sac_writetrace(f3,t3)

end program splitwave
!-------------------------------------------------------------------------------