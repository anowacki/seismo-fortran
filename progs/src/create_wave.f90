!===============================================================================
program create_wave
!===============================================================================
!  Creates a synthetic waveform polarised to taste, with the length and
!  delta requested.
!
!  Options are:
!      -type g(aussian1) r(icker): default gaussian first derivative
!      -f [dominant freq]        : default 0.1
!      -delta [sampling rate]    : default 0.05
!      -noise [amplitude]        : default 0.0
!      -spol  [source poln]      : default 0.0 (N)

   use f90sac
   use splitwave

   implicit none

   type(SACtrace) :: E, N, Z
   real :: delta, freq, noise, spol
   integer :: nerr, iarg, iskip
   character(len=4) :: wavetype
   character(len=250) :: arg, file, fname

   ! Set defaults
   fname = 'wave'
   freq = 0.1
   delta = 0.05
   noise = 0.0
   wavetype = 'g'
   spol = 0.0
   iskip = 0

   if (modulo(command_argument_count(),2) /= 0) then
      write(*,'(a)') 'create_wave'
      write(*,'(a)') '  Creates a synthetic waveform to split with splitwave.'
      write(*,'(a)') 'Options are:'
      write(*,'(a)') '    -name    [file name]         : default ''wave.BH?'' '
      write(*,'(a)') '    -t(ype)  g(aussian1)|r(icker): default gaussian first derivative'
      write(*,'(a)') '    -f(req)  [dominant freq]     : default 0.1'
      write(*,'(a)') '    -d(elta) [sampling rate]     : default 0.05'
      write(*,'(a)') '    -n(oise) [amplitude]         : default 0.0'
      write(*,'(a)') '    -s(pol)  [source poln.]      : default 0.0 (N)'
      stop
   endif

   ! Read options from command line
   do iarg=1,command_argument_count()
      if (iarg>iskip) then
         call get_command_argument(iarg,arg)
         if (arg(1:1) /= '-') then
            write(*,'(a,a)') 'Unrecognised option ',arg
            stop
         else if (arg(1:5) == '-name') then
            call get_command_argument(iarg+1,fname)
            iskip = iarg+1
         else if (arg(1:2) == '-f') then       ! dominant frequency
            call get_command_argument(iarg+1,arg)
            read(arg,*,err=901) freq
            iskip = iarg+1
         else if (arg(1:2) == '-n') then       ! noise
            call get_command_argument(iarg+1,arg)
            read(arg,*,err=901) noise
            iskip = iarg+1
            if (noise < 0.0) noise = abs(noise)
            if (noise > 1.0) then
               write(*,'(a,a)') 'Value for noise too large: ',noise
               write(*,'(a)')   'Please give a value between 0 and 1.'
               stop
            endif
         else if (arg(1:2) == '-s') then      ! spol
            call get_command_argument(iarg+1,arg)
            read(arg,*,err=901) spol
            iskip = iarg+1
         else if (arg(1:2) == '-d') then      ! delta
            call get_command_argument(iarg+1,arg)
            read(arg,*,err=901) delta
            iskip = iarg+1
         else if (arg(1:2) == '-t') then      ! type
            call get_command_argument(iarg+1,arg)
            iskip = iarg+1
            if (arg(1:1) == 'g' .or. arg(1:1) == 'G') wavetype = 'g'
            if (arg(1:1) == 'r' .or. arg(1:1) == 'R') wavetype = 'r'
            if (arg(1:1) /= 'g' .and. arg(1:1) /= 'r' .and. &
                arg(1:1) /= 'G' .and. arg(1:1) /= 'R') then
               write(*,'(a,a)') 'Unrecognised wave type ',arg
               write(*,'(a)')   'Specify g(aussian1) or r(icker).'
               stop
            endif
         endif
      endif
   enddo
   
   call sw_create_wave(E,N,Z, freq=freq, delta=delta, noise=noise, spol=spol)

   ! Write the traces out
   file = trim(fname) // '.BHN' ; call f90sac_writetrace(file,N)
   file = trim(fname) // '.BHE' ; call f90sac_writetrace(file,E)
   file = trim(fname) // '.BHZ' ; call f90sac_writetrace(file,Z)

   stop
901   write(*,'(a)') 'Bad value supplied.'
end program create_wave
!-------------------------------------------------------------------------------
