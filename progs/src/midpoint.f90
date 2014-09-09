!===============================================================================
! One of Andy Nowacki's Fortran utility modules for dealing with seismic
! anisotropy and other problems.
!
! Andy Nowacki <andy.nowacki@bristol.ac.uk>
!
! See the file LICENCE for licence details.
!===============================================================================
program midpoint
!===============================================================================
!  Compute the midpoint between two points on a sphere.  This use case is common
!  enough to warrant its own program.

   use spherical_geometry, only: azimuth, delta, step, sg_angle_180

   implicit none
   
   integer, parameter :: r8 = selected_real_kind(15,307)
   integer, parameter :: rs = r8
   real(rs) :: lon1,lat1,lon2,lat2,lonm,latm,a,d
   real(rs), parameter :: d_tol = 1e-8
   integer :: iostat = 0
   
   if (command_argument_count() > 0) call usage()
   
   do while (iostat == 0)
      call read_vals()
      if (iostat < 0) exit ! EOF
      d = delta(lon1,lat1,lon2,lat2,degrees=.true.)
      if (abs(180._rs - d) < d_tol) then
         write(0,'(4(a,f0.8),a)') 'midpoint: Warning: points (',lon1,',',lat1, &
            ') and (',lon2,',',lat2,') are antipodal so midpoint isarbitrary'
      endif
      a = azimuth(lon1,lat1,lon2,lat2,degrees=.true.)
      call step(lon1,lat1,a,d/2._rs,lonm,latm,degrees=.true.)
      write(*,'(f13.8,1x,f12.8)') lonm,latm
   enddo

contains
   !===============================================================================
   subroutine read_vals()
   !===============================================================================
      implicit none
      read(*,*,iostat=iostat) lon1,lat1,lon2,lat2
      lon1 = sg_angle_180(lon1)
      lon2 = sg_angle_180(lon2)
      if (iostat > 0) then
         write(0,'(a)') 'midpoint: Error: problem reading lon1,lat1,lon2,lat2 from stdin'
         stop
      endif
      if (abs(lat1) > 90._rs .or. abs(lat2) > 90._rs) then
         write(0,'(a)') 'midpoint: Error: latitude must be in range -90 to 90 degrees.'
         stop
      endif
   end subroutine read_vals
   !-------------------------------------------------------------------------------

   !===============================================================================
   subroutine usage()
   !===============================================================================
      implicit none
      write(0,'(a)') 'Usage: midpoint < [lon1] [lat1] [lon2] [lat2]',&
                     'Compute the midpoint between two points on a sphere'
      stop
   end subroutine usage
   !-------------------------------------------------------------------------------

end program midpoint
!-------------------------------------------------------------------------------