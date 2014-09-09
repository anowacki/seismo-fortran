!===============================================================================
! One of Andy Nowacki's Fortran utility modules for dealing with seismic
! anisotropy and other problems.
!
! Andy Nowacki <andy.nowacki@bristol.ac.uk>
!
! See the file LICENCE for licence details.
!===============================================================================
program split_misfit
!===============================================================================
!  split_misfit takes (sets) of two splitting parameters and gives a misfit
!  between the two.  It does so using a method developed by James Wookey and
!  Andrew Walker <j.wookey@bristol.ac.uk>, <andrew.walker@bristol.ac.uk>.  The
!  details are outlined in the code (splitwave.f90: sw_misfit), but are also
!  described at:
!     https://github.com/andreww/MSAT/
!  and more specifically in:
!     https://github.com/andreww/MSAT/blob/master/msat/MS_splitting_misfit.m

   use splitwave
   
   implicit none
   
   real :: phi1, dt1, phi2, dt2
   real :: freq = 0.1, delta = 0.1, noise = 0.05, spol = 0.
   character(len=250) :: arg
   character(len=1) :: wavetype = 'g'
   integer :: iostat = 0, iarg, narg
   logical :: read_stdin = .false., spol_set = .false., got_cl_splits = .false.
   
   call get_options
   
   iostat = 0
   do while (iostat == 0)
      if (read_stdin) then
         if (spol_set) then
            read(*,*,iostat=iostat) phi1, dt1, phi2, dt2
            call check_read('phi1, dt1, phi2, dt2')
         else
            read(*,*,iostat=iostat) phi1, dt1, phi2, dt2, spol
            call check_read('phi1, dt1, phi2, dt2, spol')
         endif
         if (iostat < 0) exit
      endif
      write(*,*) phi1, dt1, phi2, dt2, spol, &
         sw_misfit(phi1, dt1, phi2, dt2, spol=spol, freq=freq, &
         delta=delta, noise=noise, wavetype=wavetype)
      if (.not.read_stdin) exit
   enddo
   
contains

   !============================================================================
   subroutine get_options()
   !============================================================================
      implicit none
      
      narg = command_argument_count()
      ! Incorrect invocation
      if (narg == 0) call usage
      ! Check whether we're reading from stdin
      if (narg > 0) then
         call get_command_argument(narg,arg)
         if (arg == '-') then
            read_stdin = .true.
            write(0,'(a)') &
               'split_misfit: Reading values of phi1,dt1,phi2,dt2(,spol) from stdin...'
         endif
      endif
      ! Get other options
      iarg = 1
      arg_loop: do while (iarg <= narg)
         call get_command_argument(iarg,arg)
         select case(arg(1:2))
            case('-d')
               call check_nargs(1)
               call get_command_argument(iarg+1,arg)
               read(arg,*,iostat=iostat) delta
               call check_read('delta')
               iarg = iarg + 2
            case('-f')
               call check_nargs(1)
               call get_command_argument(iarg+1,arg)
               read(arg,*,iostat=iostat) freq
               call check_read('freq')
               iarg = iarg + 2
            case('-n')
               call check_nargs(1)
               call get_command_argument(iarg+1,arg)
               read(arg,*,iostat=iostat) noise
               call check_read('noise')
               iarg = iarg + 2
            case('-s')
               spol_set = .true.
               call get_command_argument(iarg+1,arg)
               read(arg,*,iostat=iostat) spol
               call check_read('spol')
               iarg = iarg + 2
            case('-t')
               call check_nargs(1)
               call get_command_argument(iarg+1,arg)
               wavetype = arg(1:1)
               call check_read('wavetype')
               iarg = iarg + 2
            case('-h')
               call usage
            case default
               if (.not.read_stdin) then
                  call get_cl_splits
                  exit arg_loop
               else
                  if (iarg == narg) exit arg_loop
                  ! Unknown option
                  write(0,'(a)') 'split_misfit: Unknown option "'//trim(arg)//'"'
                  call usage
               endif
         end select
      enddo arg_loop
      
      if (.not.read_stdin .and. .not.got_cl_splits) then
         write(0,'(a)') 'split_misfit: no command-line splitting parameters given'
         call usage
      endif
      
   end subroutine get_options
   !----------------------------------------------------------------------------
   
   !============================================================================
   subroutine get_cl_splits()
   !============================================================================
   !  Read the two splitting operators from the command line, starting at iarg
      implicit none
      if (iarg + 3 > narg) call usage
      
      read(arg,*,iostat=iostat) phi1
      call check_read('phi1')
      call get_command_argument(iarg+1,arg)
      read(arg,*,iostat=iostat) dt1
      call check_read('dt1')
      call get_command_argument(iarg+2,arg)
      read(arg,*,iostat=iostat) phi2
      call check_read('phi2')
      call get_command_argument(iarg+3,arg)
      read(arg,*,iostat=iostat) dt2
      call check_read('dt2')
      got_cl_splits = .true.
      if (iarg + 4 == narg) then
         call get_command_argument(iarg+4,arg)
         read(arg,*,iostat=iostat) spol
         call check_read('fixed spol')
      endif
   end subroutine get_cl_splits
   !----------------------------------------------------------------------------

   !============================================================================
   subroutine check_read(str)
   !============================================================================
      implicit none
      character(len=*), intent(in) :: str
      if (iostat > 0) then
         write(0,'(a)') &
            'split_misfit: Cannot get value of '//trim(str)//' from argument "' &
            //trim(arg)//'"'
            stop
!       else
!          write(0,'(a)') 'Got value of '//trim(str)//' as "'//trim(arg)//'"'
      endif
   end subroutine check_read
   !----------------------------------------------------------------------------
   
   !============================================================================
   subroutine check_nargs(n)
   !============================================================================
      implicit none
      integer, intent(in) :: n
      if (iarg + n > narg) then
         write(0,'(a)') 'split_misfit: Not enough arguments supplied for option ' &
         //arg(1:2)
         call usage
      endif
   end subroutine check_nargs
   !----------------------------------------------------------------------------
   
   !============================================================================
   subroutine usage()
   !============================================================================
      implicit none
      integer :: iout
      ! Set whether we write to stdout or stderr based on whether we asked to
      ! call this routine
      if (arg(1:2) == '-h') then
         iout = 6
      else
         iout = 0
      endif
      write(iout,'(a)') &
         'split_misfit calculates a misfit between two splitting operators', &
         '',&
         'Usage: split_misfit (options) [phi1] [dt1] [phi2] [dt2] (spol)', &
         '   or', &
         '       split_misfit (options) - < (phi1,dt1,phi2,dt2,(spol) on stdin)', &
         'Options:', &
         '   -d [delta]   : Set sampling interval of wave / s [0.1]', &
         '   -f [freq]    : Set dominant frequency of wave / Hz [0.1]', &
         '   -n [noise]   : Set fraction of random noise [0.05]', &
         '   -s [spol]    : Set spol for all pairs of splits [read from stdin or command line]', &
         '   -t ["g"|"r"] : Set wave type: (g)aussian first derivative or (r)icker [g]', &
         '   -h(elp)      : Prints this message', &
         '   -            : Read values from stdin'
      stop
   end subroutine usage
   !----------------------------------------------------------------------------

end program split_misfit
!-------------------------------------------------------------------------------