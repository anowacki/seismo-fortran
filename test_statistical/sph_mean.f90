!===============================================================================
program spherical_mean
!===============================================================================
! Reads in (lon, colat) angles from stdin and computes the mean direction.
! Input and output in degrees

use statistical

implicit none

integer, parameter :: rs = 8
integer, parameter :: nmax = 100000
real(rs), parameter :: pi = 4._rs*atan2(1._rs, 1._rs)
real(rs), parameter :: rad = pi/180._rs, deg = 180._rs/pi
real(rs), dimension(nmax) :: lon, colat
real(rs), allocatable, dimension(:) :: x, y, z
real(rs) :: mean(3), meanlon, meancolat
integer :: n, i, iostat = 0

i = 1
do while (iostat == 0)
   read(*,*,iostat=iostat) lon(i), colat(i)
   if (iostat < 0) then
      n = i - 1
      exit
   else if (iostat > 0) then
      write(0,'(a)') 'sph_mean: Problem reading lon,colat from stdin'
      stop
   endif
   if (colat(i) < 0._rs .or. colat(i) > 180._rs) then
      write(0,'(a)') 'sph_mean: Error: colat is outside range (0, 180) degrees'
      stop
   endif
   i = i + 1
   if (i > nmax) then
      write(0,'(a,i0)') 'sph_mean: Reached compiled nmax limit of ',nmax
      write(0,'(a)') 'Skipping remaining points...'
   endif
enddo

write(*,'(a,i0,a)') 'Got ',n,' pairs from stdin'

allocate(x(n), y(n), z(n))

! Convert to x,y,z
lon(1:n) = rad*lon(1:n)
colat(1:n) = rad*colat(1:n)
x(1:n) = sin(colat(1:n))*cos(lon(1:n))
y(1:n) = sin(colat(1:n))*sin(lon(1:n))
z(1:n) = cos(colat(1:n))

mean = sph_mean(x, y, z)

! Convert back to lon, colat
meancolat = deg*acos(mean(3))
meanlon = deg*atan2(mean(2),mean(1))

write(*,*) meanlon, meancolat

end program spherical_mean
!-------------------------------------------------------------------------------