!===============================================================================
program project_to_gcp_prog
!===============================================================================
!  Wrapper program for sg_project_to_gcp.  This reads the coordinates of the pole
!  to a great circle on the unit sphere and another point and projects that point
!  onto the great circle.
!  Input is in geographic coordinates in degrees.

   use spherical_geometry, only: sg_project_to_gcp

   implicit none
   real(8) :: long, latg, lonp, latp, lon, lat
   integer :: iostat
      
   if (command_argument_count() /= 0) call usage
   
   ! Loop over lines of stdin
   iostat = 0
   do while (iostat == 0)
      read(*,*,iostat=iostat) long,latg,lonp,latp
      if (iostat > 0) then
         write(0,'(a)') 'project_to_gcp: Error: problem reading lon,lat,lon,lat from stdin'
         stop
      endif
      if (iostat < 0) exit
      call sg_project_to_gcp(long,latg,lonp,latp,lon,lat,degrees=.true.)
      write(*,'(f11.6,1x,f10.6)') lon,lat
   enddo

contains
   !============================================================================
   subroutine usage()
   !============================================================================
      implicit none
      write(0,'(a)') 'Usage: project_to_gcp < [lon pole] [lat pole] [lon point] [lat point]',&
                     '   Project a point onto a great circle',&
                     '   First two coordinates are the pole to the great circle, second two are',&
                     '   coordinates of the point to be projected onto the great circle'
      stop
   end subroutine usage
   !----------------------------------------------------------------------------

end program project_to_gcp_prog
!-------------------------------------------------------------------------------