!===============================================================================
! One of Andy Nowacki's Fortran utility modules for dealing with seismic
! anisotropy and other problems.
!
! Andy Nowacki <andy.nowacki@bristol.ac.uk>
!
! See the file LICENCE for licence details.
!===============================================================================
program lonlat2xyz
!===============================================================================
!  Reads sets of lon,lat(,radius) from stdin and writes out x,y,z

   use spherical_geometry, only: geog2cart

   implicit none

   real(8) :: lon, lat, r, x, y, z
   character(len=250) :: arg
   integer :: iostat
   logical :: fixed_rad = .false.

!  Check input arguments
   call get_args
   
   ! Get coordinates from stdin
   iostat = 0
   do while (iostat == 0)
      if (fixed_rad) then
         read(*,*,iostat=iostat) lon, lat
      else
         read(*,*,iostat=iostat) lon, lat, r
      endif
      if (iostat < 0) exit
      if (iostat > 0) then
         write(0,'(a)') 'lonlat2xyz: problem reading coordinates from stdin.'
         stop 1
      endif
      call geog2cart(lon, lat, r, x, y, z, degrees=.true.)
      write(*,*) x, y, z
   enddo

contains
   !===============================================================================
   subroutine get_args()
   !===============================================================================
      integer :: iarg, narg
      narg = command_argument_count()
      iarg = 1
      do while (iarg <= narg)
         call get_command_argument(iarg, arg)
         select case(arg)
            case('-r')
               if (narg < iarg + 1) call usage
               call get_command_argument(iarg+1, arg)
               read(arg,*,iostat=iostat) r
               if (iostat /= 0) then
                  write(0,'(a)') 'lonlat2xyz: Error: Can''t read radius from ' &
                     // 'argument "' // trim(arg) // '"'
                  stop 1
               endif
               fixed_rad = .true.
               iarg = iarg + 2
            case('-h')
               call usage
            case default
               call usage
         end select
      enddo
   end subroutine get_args
   !-------------------------------------------------------------------------------

   !===============================================================================
   subroutine usage()
   !===============================================================================
      write(0,'(a)') &
         'Usage: lonlat2xyz (options) < [lon] [lat] [radius]', &
         'Reads geographic coordinates from stdin and writes cartesian to stdout.', &
         'Options:', &
         '   -r [radius] : Fix radius.  If not present, r must be on stdin.'
      stop 1
   end subroutine usage
   !-------------------------------------------------------------------------------
   
end program lonlat2xyz
!-------------------------------------------------------------------------------
