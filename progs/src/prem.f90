!===============================================================================
! One of Andy Nowacki's Fortran utility modules for dealing with seismic
! anisotropy and other problems.
!
! Andy Nowacki <andy.nowacki@bristol.ac.uk>
!
! See the file LICENCE for licence details.
!===============================================================================
program PREM_vals

use global_1d_models

implicit none

character(len=50) :: arg
character(len=:), parameter :: fmt_iso = '(3f7.3, f6.1, f6.2)', &
                               fmt_iso_depth = '(3f7.3, f6.1, f6.2, f9.3)', &
                               fmt_ani = '(4f7.3, f7.4, f7.3, f6.1, f6.2)', &
                               fmt_ani_depth = '(4f7.3, f7.4, f7.3, f6.1, f6.2, f9.3)'
character(len=:), parameter :: hdr_iso = '#    vp     vs    rho     P     g'
character(len=:), parameter :: hdr_ani = '#   vph    vpv    vsh    vsv    eta    rho     P     g'
character(len=6) :: label = ' depth'
real(8) :: depth, vp, vs, rho, P, g, vph, vpv, vsh, vsv, eta
integer :: iostatus
logical :: aniso = .false., header = .true., radius = .false., read_stdin = .true.

call get_args

!  Use the argument on the command line
if (.not. read_stdin) then
   if (radius) depth = PREM_radius_km - depth
   if (aniso) then
      if (header) write(*,'(a)') hdr_ani
      call prem(depth, vph=vph, vpv=vpv, vsh=vsh, vsv=vsv, eta=eta, rho=rho, g=g)
   else
      if (header) write(*,'(a)') hdr_iso
      call prem(depth, vp=vp, vs=vs, rho=rho, g=g)
   endif
   P = pressure(depth, model='PREM')
   if (aniso) then
      write(*,fmt_ani) vph, vpv, vsh, vsv, eta, rho, P, g
   else
      write(*,fmt_iso) vp, vs, rho, P, g
   endif

else

   if (aniso) then
      if (header) write(*,'(a," ",a)') hdr_ani, trim(label)
   else
      if (header) write(*,'(a," ",a)') hdr_iso, trim(label)
   endif
   iostatus = 0
   do while (iostatus == 0)
      read(*,*,iostat=iostatus) depth
      if (iostatus < 0) stop
      if (iostatus > 0) then
         write(0,'(a)') 'prem: Error: Problem reading value of depth on stdin'
         error stop
      endif
      if (radius) depth = PREM_radius_km - depth
      if (aniso) then
         call prem(depth, vph=vph, vpv=vpv, vsh=vsh, vsv=vsv, eta=eta, rho=rho, g=g)
      else
         call prem(depth, vp=vp, vs=vs, rho=rho, g=g)
      endif
      P = pressure(depth,model='PREM')
      if (aniso) then
         write(*,fmt_ani_depth) vph, vpv, vsh, vsv, eta, rho, P, g, depth
      else
         write(*,fmt_iso_depth) vp, vs, rho, P, g, depth
      endif
   enddo

endif


contains
   subroutine usage
      write(0,'(a)') &
         'Usage: prem (options) [depth / km]', &
         '   or: prem (options) < [depth(s) on stdin]', &
         'Returns the values of PREM (Dziewonski & Anderson, 1981) at a given depth,', &
         'in km/s and g/cm^3, and the pressure and gravity predicted by them.', &
         'If reading values from stdin, the depths are also written out as the final column.', &
         'Options:', &
         '   -a:  Evaluate anisotropic PREM, writing out Vph, Vpv, Vsh, Vsv, eta, rho, P and g', &
         '   -n:  Do not print the header', &
         '   -r:  Read radi(us/i) instead of depth(s), in km'
      error stop
   end subroutine usage

   subroutine get_args
      integer :: iarg, narg
      character(len=250) :: arg
      narg = command_argument_count()
      if (narg == 0) then
         read_stdin = .true.
         return
      endif

      iarg = 1
      do while (iarg <= narg)
         call get_command_argument(iarg, arg)
         select case(arg)
            case ('-a')
               aniso = .true.
               iarg = iarg + 1
            case ('-n')
               header = .false.
               iarg = iarg + 1
            case ('-r')
               radius = .true.
               label = 'radius'
               iarg = iarg + 1
            case default
               if (iarg /= narg) call usage
               read_stdin = .false.
               read(arg,*,iostat=iostatus) depth
               if (iostatus /= 0) then
                  write(0,'(a)') 'prem: Error: Cannot get value of depth from argument "' &
                     // trim(arg) // '"'
                  call usage
               endif
               return
         end select
      enddo
   end subroutine get_args

end program PREM_vals