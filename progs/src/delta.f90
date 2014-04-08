program delta
! Program to compute the great circle angle between two points
! on a sphere (e.g. Earth) given latitude and longitude
! Reads from standard input

implicit none

! variables
double precision :: lat1,lon1,lat2,lon2
double precision :: angle,pi
integer          :: iostatus

! write (*,*)"longitude1, latitude1, longitude2, latitude2"
iostatus = 0
do while (iostatus == 0)
	read(*,*,iostat=iostatus) lon1,lat1,lon2,lat2
	if (iostatus < 0) exit
	if (iostatus > 0) stop 'Some problem reading coordinates'
	write(*,'(f10.6," degrees")')angle(lat1,lon1,lat2,lon2)
enddo
end program





function angle(lat1,lon1,lat2,lon2)
implicit none
double precision :: angle,lat1,lon1,lat2,lon2,pi

pi=3.14159265358979323846264338327950288D0

lat1=lat1*pi/1.8D2 ; lon1=lon1*pi/1.8D2
lat2=lat2*pi/1.8D2 ; lon2=lon2*pi/1.8D2

angle=atan2( sqrt( (cos(lat2)*sin(lon2-lon1))**2 + (cos(lat1)*sin(lat2) - &
		sin(lat1)*cos(lat2)*cos(lon2-lon1))**2) , &
		sin(lat1)*sin(lat2) + cos(lat1)*cos(lat2)*cos(lon2-lon1))

angle=angle*1.8D2/pi

return
end function
!=ACOS(SIN(lat1)*SIN(lat2)+COS(lat1)*COS(lat2)*COS(lon2-lon1))*6371