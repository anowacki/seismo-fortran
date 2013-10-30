!==============================================================================
program splitwaveN
!==============================================================================
!  'Split' an artificial shear wave 'wave.BH?' using J. Wookey's f90sac routines.
!  Reads several sets of splitting operators from stdin.
!  Doesn't actually need the vertical component, which can be ignored.
!  
!  Usage:
!        splitwave (filename)
!        Reads splitting parameters from stdin
   
   use f90sac
   use splitwave
   
   implicit none
   
   type(SACtrace) :: t1, t2, t3
   real :: theta, dt
   character(len=250) :: f, f1, f2, f3
   integer :: iarg, iostatus
   
   
   if (command_argument_count() /= 0 .and. command_argument_count() /= 1) then
      write(0,'(a)') 'Usage:  splitwave [filename] < [splitting parameters file]'
      write(0,'(a)') '    or  splitwave (operates on ''wave.BH?'') < [splits]'
      write(0,'(a)') 'Splitting parameters file should contain [phi] [dt]'
      stop
   endif
   
!  Get filename base
   if (command_argument_count() == 1) then
      call get_command_argument(1,f)
      iarg = 1
   else
      f = 'wave'
      iarg = 0
   endif
   
!  Filenames
   f1 = trim(f) // '.BHN'
   f2 = trim(f) // '.BHE'
   f3 = trim(f) // '.BHZ'
   
!  Load traces
   call f90sac_readtrace(f1,t1)
   call f90sac_readtrace(f2,t2)
   call f90sac_readtrace(f3,t3)
   
   iostatus = 0
!  Loop over the splitting parameters
   do while (iostatus == 0)
      read(*,*,iostat=iostatus) theta,dt
      if (iostatus < 0) exit
      if (iostatus > 0) then
         write(0,'(a)') &
            'splitwaveN: Error: Some problem reading in splitting parameters.'
         stop
      endif
      call sw_splitN(t1,t2,t3, 1, (/theta/), (/dt/), quiet=.true.)
   enddo
   
!  Write out altered traces
   call f90sac_writetrace(f1,t1)
   call f90sac_writetrace(f2,t2)
   call f90sac_writetrace(f3,t3)
      
end program splitwaveN
!------------------------------------------------------------------------------
