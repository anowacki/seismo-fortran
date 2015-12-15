!===============================================================================
! One of Andy Nowacki's Fortran utility modules for dealing with seismic
! anisotropy and other problems.
!
! Andy Nowacki <andy.nowacki@bristol.ac.uk>
!
! See the file LICENCE for licence details.
!===============================================================================
program delta
! Program to compute the great circle angle between two points
! on a sphere (e.g. Earth) given latitude and longitude
! Reads from standard input

implicit none

! variables
double precision :: lat1,lon1,lat2,lon2
double precision :: angle
integer          :: iostatus = 0

call get_args

do while (iostatus == 0)
   read(*,*,iostat=iostatus) lon1,lat1,lon2,lat2
   if (iostatus < 0) exit
   if (iostatus > 0) stop 'Some problem reading coordinates'
   write(*,'(f10.6)') angle(lat1,lon1,lat2,lon2)
enddo

contains

subroutine get_args()
   integer :: narg, iarg
   narg = command_argument_count()
   if (narg > 0) call usage
end subroutine get_args

subroutine usage()
   write(0,'(a)') 'Usage: delta < [(lon1, lat1, lon2, lat2) on stdin]', &
                  'Reads pairs of geographic coordinates in degrees and writes', &
                  'great circle distance on a sphere (in degrees) to stdout.'
   stop
end subroutine usage

end program delta


function angle(lat1i,lon1i,lat2i,lon2i)
implicit none
double precision :: angle
double precision, intent(in) :: lat1i,lon1i,lat2i,lon2i
double precision :: lat1,lon1,lat2,lon2
double precision, parameter :: pi=3.14159265358979323846264338327950288D0

lat1=lat1i*pi/1.8D2 ; lon1=lon1i*pi/1.8D2
lat2=lat2i*pi/1.8D2 ; lon2=lon2i*pi/1.8D2

angle=atan2( sqrt( (cos(lat2)*sin(lon2-lon1))**2 + (cos(lat1)*sin(lat2) - &
                   sin(lat1)*cos(lat2)*cos(lon2-lon1))**2) , &
            sin(lat1)*sin(lat2) + cos(lat1)*cos(lat2)*cos(lon2-lon1))

angle=angle*1.8D2/pi

end function angle

