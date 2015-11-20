!===============================================================================
! One of Andy Nowacki's Fortran utility modules for dealing with seismic
! anisotropy and other problems.
!
! Andy Nowacki <andy.nowacki@bristol.ac.uk>
!
! See the file LICENCE for licence details.
!===============================================================================
program xyz2lonlat
!===============================================================================
!  Reads sets of x,y,z from stdin and writes out lon,lat,radius

   use spherical_geometry, only: cart2geog

   implicit none

   real(8) :: lon, lat, r, x, y, z, val, norm_rad
   character(len=250) :: arg
   integer :: iostat
   logical :: norm = .false., value = .false.

!  Check input arguments
   call get_args
   
   ! Get coordinates from stdin
   iostat = 0
   do while (iostat == 0)
      if (value) then
         read(*,*,iostat=iostat) x, y, z, val
      else
         read(*,*,iostat=iostat) x, y, z
      endif
      if (iostat < 0) exit
      if (iostat > 0) then
         write(0,'(a)') 'lonlat2xyz: problem reading coordinates from stdin.'
         stop 1
      endif
      call cart2geog(x, y, z, lat, lon, r, degrees=.true.)
      if (norm) r = r*norm_rad
      if (value) then
         write(*,*) lon, lat, r, val
      else
         write(*,*) lon, lat, r
      endif
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
            case('-n')
               norm = .true.
               call get_command_argument(iarg+1, arg)
               read(arg,*) norm_rad
               iarg = iarg + 2
            case('-v')
               value = .true.
               iarg = iarg + 1
            case('-h')
               call usage(.true.)
            case default
               call usage
         end select
      enddo
   end subroutine get_args
   !-------------------------------------------------------------------------------

   !===============================================================================
   subroutine usage(help)
   !===============================================================================
      logical, intent(in), optional :: help
      integer :: unit
      unit = 0
      if (present(help)) then
         if (help) unit = 6
      endif
      write(0,'(a)') &
         'Usage: xyz2lonlat (options) < [x] [y] [z] (value)', &
         'Reads cartesian coordinates from stdin and writes geographic to stdout.', &
         'Options:', &
         '   -n [norm] : Scale the points radius by <norm>', &
         '   -v        : Output the value in the third/fourth column'
      if (present(help)) then
         if (help) stop
      endif
      error stop
   end subroutine usage
   !-------------------------------------------------------------------------------
   
end program xyz2lonlat
!-------------------------------------------------------------------------------
