program test

use spherical_geometry

implicit none

real(8) :: x1,y1,d,az,x2,y2

x1 = 0.
y1 = 0.
d = 45.
az = 0.

write(*,*) 'start lon,lat = ',x1,y1
write(*,*) 'azimuth = ',az
write(*,*) 'distance = ',d

call step(x1,y1,az,d,x2,y2,degrees=.true.)

write(*,*) 'expected lon,lat = ',0.,45.
write(*,*) 'calculated lon,lat = ',x2,y2


end program test