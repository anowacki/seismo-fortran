!===============================================================================
program test_sg_gcp_to_azimuth
!===============================================================================
!  Test the sg_gcp_to_azimuth subroutine.
!  Usage: test_sg_gcp_to_azimuth [long] [latg] [lonp] [latp]

   use spherical_geometry, only: sg_gcp_to_azimuth

   implicit none
   integer, parameter :: rs = 8
   real(rs) :: long,latg,lonp,latp,azi
   character(len=250) :: arg
   
   ! Pole to Greenwich meridian
   long = -90.
   latg = 0.
   
   ! Point on the equator
   lonp = 0.
   latp = 0.
   
   ! If we've supplied some arguments, use them
   if (command_argument_count() == 4) then
      call get_command_argument(1,arg)
      read(arg,*) long
      call get_command_argument(2,arg)
      read(arg,*) latg
      call get_command_argument(3,arg)
      read(arg,*) lonp
      call get_command_argument(4,arg)
      read(arg,*) latp
   endif
   
   write(*,'(f0.6)') sg_gcp_to_azimuth(long,latg,lonp,latp,degrees=.true.)

end program test_sg_gcp_to_azimuth
!-------------------------------------------------------------------------------