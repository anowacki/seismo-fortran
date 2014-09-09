!===============================================================================
! One of Andy Nowacki's Fortran utility modules for dealing with seismic
! anisotropy and other problems.
!
! Andy Nowacki <andy.nowacki@bristol.ac.uk>
!
! See the file LICENCE for licence details.
!===============================================================================
program split_misfit_ecs
!===============================================================================
! split_misfit_ecs compares a set of azimuths, inclinations and splitting
! parameters to those predicted by an elasticity tensor and computes the misfit
! as calculated by Wookey and Walker (MSAT).

   use anisotropy_ajn, only: CIJ_load
   use splitwave, only: sw_misfit_ecs
   
   implicit none
   
   integer, parameter :: nmax = 10000
   real(8), dimension(nmax) :: az, inc, phi, dt, spol, misfit
   real(8) :: freq = 0.1, delta = 0.05, noise = 0.05
   real(8) :: t ! Layer thickness
   real(8) :: C(6,6), rho
   character(len=3000) :: arg
   character(len=1) :: wavetype = 'g'
   integer :: iostat = 0, iarg, narg, n, i, j
   logical :: read_stdin = .false., spol_set = .false., got_cl_splits = .false., &
              thickness_set = .false.

   call get_options
   ! Get ECs from first line of stdin if necessary
   if (read_stdin) call read_ecs_stdin
   ! Fill arrays with angles and splits
   call read_splits_stdin
   ! Calculate splits
   if (thickness_set) then
      call sw_misfit_ecs(C, az(1:n), inc(1:n), phi(1:n), dt(1:n), spol(1:n), &
         misfit(1:n), t=t, freq=freq, delta=delta, noise=noise, wavetype=wavetype)
   else
      call sw_misfit_ecs(C, az(1:n), inc(1:n), phi(1:n), dt(1:n), spol(1:n), &
         misfit(1:n), freq=freq, delta=delta, noise=noise, wavetype=wavetype)
   endif
   ! Write out results
   do i=1,n
      write(*,*) az(i), inc(i), phi(i), dt(i), spol(i), misfit(i)
   enddo


contains
   !============================================================================
   subroutine get_options()
   !============================================================================
      implicit none
      
      narg = command_argument_count()
      ! Incorrect invocation
      if (narg < 1) call usage
      ! Check whether we're reading from stdin
      call get_command_argument(narg,arg)
      if (arg == '-') then
         read_stdin = .true.
         write(0,'(a)') &
            'split_misfit_ecs: Reading density-normalised Cij from stdin...'
      else
         call CIJ_load(arg, C, rho)
         C = C/rho
      endif
      ! Get other options
      iarg = 1
      arg_loop: do while (iarg <= narg - 1)
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
            case('-l')
               call check_nargs(1)
               call get_command_argument(iarg+1,arg)
               read(arg,*,iostat=iostat) t
               call check_read('layer thickness')
               thickness_set = .true.
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
               ! Unknown option
               write(0,'(a)') 'split_misfit: Unknown option "'//trim(arg)//'"'
               call usage
         end select
      enddo arg_loop
      
   end subroutine get_options
   !----------------------------------------------------------------------------

   !============================================================================
   subroutine read_ecs_stdin()
   !============================================================================
      implicit none
      read(*,'(a)',iostat=iostat) arg
      if (iostat /= 0) then
         write(0,'(a)') 'split_misfit_ecs: Can''t read 36 ecs from stdin'
         stop
      endif
      read(arg,*,iostat=iostat) ((C(i,j),j=1,6),i=1,6)
      if (iostat /= 0) then
         write(0,'(a)') 'split_misfit_ecs: Problem reading 36 ecs from stdin ' &
            // 'line: "'//trim(arg)//'"'
         stop
      endif
   end subroutine read_ecs_stdin
   !----------------------------------------------------------------------------

   !============================================================================
   subroutine read_splits_stdin()
   !============================================================================
      implicit none
      i = 0
      do
         read(*,'(a)',iostat=iostat) arg
         if (iostat < 0) exit
         if (iostat > 0) then
            write(0,'(a)') 'split_misfit_ecs: Error reading splits from stdin'
            stop
         endif
         i = i + 1
         if (i > nmax) then
            write(0,'(a)') 'split_misfit_ecs: Number of splits larger than ' &
               // 'maximum compiled limits.  Use fewer values or recompile ' &
               // 'with larger nmax.'
               stop
         endif
         read(arg,*,iostat=iostat) az(i), inc(i), phi(i), dt(i), spol(i)
         if (iostat /= 0) then
            write(0,'(a)') 'split_misfit_ecs: Problem getting az, inc, phi, dt, spol' &
               // ' from line "'//trim(arg)//'"'
            stop
         endif
      enddo
      ! Total number of splits read on stdin
      n = i
   end subroutine read_splits_stdin
   !----------------------------------------------------------------------------

   !============================================================================
   subroutine check_read(str)
   !============================================================================
      implicit none
      character(len=*), intent(in) :: str
      if (iostat > 0) then
         write(0,'(a)') &
            'split_misfit_ecs: Cannot get value of '//trim(str)//' from argument "' &
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
      integer, intent(in) :: n
      if (iarg + n > narg - 1) then
         write(0,'(a)') 'split_misfit_ecs: Not enough arguments supplied for option ' &
         //arg(1:2)
         call usage
      endif
   end subroutine check_nargs
   !----------------------------------------------------------------------------

   !============================================================================
   subroutine usage()
   !============================================================================
      integer :: iout
      if (arg(1:2) == "-h") then
         iout = 5
      else
         iout = 0
      endif
      write(iout,'(a)') &
         "split_misfit_ecs calculates the splitting misfit for a set of", &
         "observations compared with those predicted by an elasticity", &
         "tensor.  Provide density-normalised tensor.  By default, the", &
         "Cij splits are normalised to have the same maximum as the", &
         "observations: specify a layer thickness with -l [thickness].", &
         "Azimuth and inclination are in CIJ_phasevels convention.", &
         "", &
         "Usage: split_misfit_ecs (options) [.ecs file] < (az,inc,phi,dt,spol on stdin)", &
         "   or", &
         "       split_misfit_ecs (options) - < (36 ecs \n az,inc,phi,dt,spol on stdin)", &
         "Options:", &
         "   -d [delta]     : Set sampling interval of wave / s [0.1]", &
         "   -f [freq]      : Set dominant frequency of wave / Hz [0.1]", &
         "   -l [thickness] : Set layer thickness over which splits accrued / km [normalised]", &
         "   -n [noise]     : Set fraction of random noise [0.05]", &
         "   -s [spol]      : Set spol for all pairs of splits [read from stdin or command line]", &
         "   -t ['g'|'r']   : Set wave type: (g)aussian first derivative or (r)icker [g]", &
         "   -h(elp)        : Prints this message"
      stop
   end subroutine usage
   !----------------------------------------------------------------------------

end program split_misfit_ecs
!-------------------------------------------------------------------------------