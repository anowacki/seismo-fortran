!===============================================================================
! One of Andy Nowacki's Fortran utility modules for dealing with seismic
! anisotropy and other problems.
!
! Andy Nowacki <andy.nowacki@bristol.ac.uk>
!
! See the file LICENCE for licence details.
!===============================================================================
program gcp_from_point_azi_prog
!===============================================================================
! Wrapper program for sg_gcp_from_point_azi.  Given a point on a sphere and an
! azimuth through that point, return the pole to the great circle along that
! azimuth.  Reads lon, lat and azi from stdin, all in degrees

   use spherical_geometry, only: sg_gcp_from_point_azi

   implicit none
   
   real(8) :: lon,lat,azi,long,latg
   integer :: iostat

   if (command_argument_count() /= 0) call usage
   
   iostat = 0
   do while (iostat == 0)
      read(*,*,iostat=iostat) lon,lat,azi
      if (iostat > 0) then
         write(0,'(a)') 'gcp_from_point_azi: Error: problem reading lon,lat,azi (degrees) from stdin'
         stop
      endif
      if (iostat < 0) exit
      call sg_gcp_from_point_azi(lon,lat,azi,long,latg,degrees=.true.)
      write(*,'(f11.6,1x,f10.6)') long,latg
   enddo

contains
   !============================================================================
   subroutine usage()
   !============================================================================
      implicit none
      write(0,'(a)') 'Usage: gcp_from_point_azi < [lon] [lat] [azi]',&
                     '   Reads point on sphere and azimuth from stdin (degrees)',&
                     '   Writes lon and lat of pole to great circle defined by the input to stdout.'
      stop
   end subroutine usage
   !----------------------------------------------------------------------------

end program gcp_from_point_azi_prog
!-------------------------------------------------------------------------------