!===============================================================================
! One of Andy Nowacki's Fortran utility modules for dealing with seismic
! anisotropy and other problems.
!
! Andy Nowacki <andy.nowacki@bristol.ac.uk>
!
! See the file LICENCE for licence details.
!===============================================================================
program baz_prog
! Program to take two sets of spherical co-ordinates
! and write out a backazimuth from the sexond to the first

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
		write(0,'(a)') 'baz: Error: problem reading coordinates from stdin.'
		stop
	endif
	write(*,'(f10.6)') azimuth(lon2,lat2,lon1,lat1,degrees=.true.)
enddo

end program baz_prog
