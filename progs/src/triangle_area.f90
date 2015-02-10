!===============================================================================
! One of Andy Nowacki's Fortran utility modules for dealing with seismic
! anisotropy and other problems.
!
! Andy Nowacki <andy.nowacki@bristol.ac.uk>
!
! See the file LICENCE for licence details.
!===============================================================================
program triangle_area
!===============================================================================
! Calculate the area of a set of triangles on a sphere.

   use spherical_geometry, only: sg_triangle_area

   implicit none

   real(8) :: lon1, lat1, lon2, lat2, lon3, lat3, A
   real(8) :: r = 1._8  ! Default to unit sphere
   character(len=250) :: arg
   integer :: iostat

   call get_args

   iostat = 0
   do while (iostat == 0)
      read(*,*,iostat=iostat) lon1, lat1, lon2, lat2, lon3, lat3
      if (iostat < 0) exit
      if (iostat > 0) then
         write(0,'(a)') 'triangle_area: Error: Cannot read coordinates from stdin'
         stop
      endif
      A = sg_triangle_area(lon1, lat1, lon2, lat2, lon3, lat3, r=r, &
         degrees=.true., quiet=.true.)
      write(*,*) A
   enddo

contains
   subroutine usage
      write(0,'(a)') &
         'Usage: triangle_area (options) < (lon1 lat1 lon2 lat2 lon3 lat3 on stdin)', &
         'Coordinates are in degrees; output is in units of r.', &
         'By default, r is 1 (case for unit sphere)', &
         'Options:', &
         '   -r [radius] : Specify radius of sphere [1]'
      stop
   end subroutine usage

   subroutine get_args()
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
                  write(0,'(a)') 'triangle_area: Error: Can''t read radius from ' &
                     // 'argument "' // trim(arg) // '"'
                  stop 1
               endif
               iarg = iarg + 2
            case('-h')
               call usage
            case default
               call usage
         end select
      enddo
   end subroutine get_args

end program triangle_area
!-------------------------------------------------------------------------------