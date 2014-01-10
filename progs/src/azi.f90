program azi_prog
! Program to take two sets of spherical co-ordinates
! and write out an azimuth from the first to the second.

use spherical_geometry, only: azimuth

implicit none

integer, parameter :: rs = selected_real_kind(15,307)
real(rs) :: lon1,lat1,lon2,lat2
integer :: iostatus

if (command_argument_count() /= 0) then
	write(0,'(a)') 'Usage: azi < [stlon stlat evlon evlat]'
	stop
endif

iostatus = 0
do while (iostatus == 0)
	read(*,*,iostat=iostatus) lon1,lat1,lon2,lat2
	if (iostatus < 0) exit
	if (iostatus > 0) then
		write(0,'(a)') 'azi: Error: problem reading coordinates from stdin.'
		stop
	endif
	write(*,'(f10.6)') azimuth(lon1,lat1,lon2,lat2,degrees=.true.)
enddo

end program azi_prog
