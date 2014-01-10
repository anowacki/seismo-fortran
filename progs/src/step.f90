program step_prog
! Computes the endpoint given a starting point lon,lat, azimuth and angular distance

use spherical_geometry, only: step

implicit none

integer, parameter :: rs = selected_real_kind(15,307) ! DP
real(rs) :: lon1,lat1,lon2,lat2,az,delta
real(rs), parameter :: pi=3.14159265358979323846264338327950288D0
integer :: iostatus

if (command_argument_count() /= 0) then
	write(0,'(a)') 'Usage: step < [lon] [lat] [azi] [delta]',&
                  'Reads a list of values from stdin and writes the endpoint lon,lat to stdout.'
	stop
endif

iostatus = 0
do while (iostatus == 0)
	read(*,*,iostat=iostatus) lon1,lat1,az,delta
	if (iostatus < 0) exit
	if (iostatus > 0) stop 'step: Error: Problem reading lon,lat,azimuth,distance from stdin.'

	if (delta > 180.) then
		write(0,*) 'step: Error: distance must be less than 180 degrees.'
      stop
	else if (lat1 <  -90._rs .or. lat1 >  90._rs) then
		write(0,*) 'step: Error: latitude must be in range -90 - 90.'
      stop
   endif
   call step(lon1,lat1,az,delta,lon2,lat2,degrees=.true.)
	write(*,'(f11.4,f10.4)') lon2,lat2
enddo

end program step_prog

