!===============================================================================
! One of Andy Nowacki's Fortran utility modules for dealing with seismic
! anisotropy and other problems.
!
! Andy Nowacki <andy.nowacki@bristol.ac.uk>
!
! See the file LICENCE for licence details.
!===============================================================================
program gcp_from_points_prog
!===============================================================================
! Wrapper program for sg_gcp_from_s.  Given a point on a sphere and an
! azimuth through that point, return the pole to the great circle along that
! azimuth.  Reads lon, lat and azi from stdin, all in degrees

   use spherical_geometry, only: sg_gcp_from_points

   implicit none

   real(8) :: lon1,lat1,lon2,lat2,long,latg
   integer :: iostat

   if (command_argument_count() /= 0) call usage

   iostat = 0
   do while (iostat == 0)
      read(*,*,iostat=iostat) lon1,lat1,lon2,lat2
      if (iostat > 0) then
         write(0,'(a)') 'gcp_from_points: Error: problem reading lon1,lat1,lon2,lat2 (degrees) from stdin'
         stop
      endif
      if (iostat < 0) exit
      call sg_gcp_from_points(lon1,lat1,lon2,lat2,long,latg,degrees=.true.)
      write(*,'(f11.6,1x,f10.6)') long,latg
   enddo

contains
   !============================================================================
   subroutine usage()
   !============================================================================
      implicit none
      write(*,'(a)') 'Usage: gcp_from_points < [lon1] [lat1] [lon2] [lat2]',&
                     '   Reads two point on sphere from stdin (degrees)',&
                     '   Writes lon and lat of pole to great circle defined by the input to stdout.'
      stop
   end subroutine usage
   !----------------------------------------------------------------------------

end program gcp_from_points_prog
!-------------------------------------------------------------------------------
