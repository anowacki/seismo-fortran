!===============================================================================
program test_SHTOOLS
!===============================================================================
!  Test the use of SHTOOLS

   use shtools
   use spherical_geometry, only: sphere_sample
   
   implicit none
   
   integer :: lmax = 360/5
   real(8),parameter :: pi = atan2(1._8,1._8)*4._8
   real(8) :: dd=5. ! approx distance between knot points
   real(8) :: sigma=0.2  ! Gives half-width of 10 degrees
   real(8),allocatable,dimension(:) :: lon,lat,c
   integer :: n  ! number of points
   real(8) :: val,vlon,vlat
   integer :: i,j
   character(len=20) :: arg
   
   real(8),allocatable :: cilm(:,:,:)

!  Get lmax from the command line, if supplied
   if (command_argument_count() > 0) then
      call get_command_argument(1,arg)
      read(arg,*) lmax
   endif
   
!  Create the knot points
   call sphere_sample(dd,lon,lat,n)

!  Make room for the values at each knot and the spherical harmonics
   allocate(c(n))
   allocate(cilm(2,lmax+1,lmax+1))
   
!  Fill the points with a bulge at the equator, width controlled by sigma
   c = exp((-(lon*pi/180._8)**2 - (lat*pi/180._8)**2)/sigma**2)

!  Alternatively, create a checkerboard
!   s%c = sign(1._8,cos(4._8*s%phi*pi/180._8)*cos(4._8*s%theta*pi/180._8))

!  Make the shape sharp sided
   where (c > 0.5)
      c = 1.
   elsewhere
      c = 0.
   endwhere

!  Fit the spherical harmonics
   write(0,'(a,i0.0,a)') 'test_SHTOOLS: Fitting spherical harmonics with lmax = ',lmax,'...'
   call SHExpandLSQ(cilm,c,lat,lon,n,lmax)
   write(0,'(a)') 'test_SHTOOLS: Finished fitting harmonics'
   
   write(0,'(a)') 'test_SHTOOLS: Evaluating interpolated values...'
!  Write out points on a regular lon,lat grid
outer: do j=-45,45,1
          do i=-45,45,1
             vlon = real(i)
             vlat = real(j)
             val = MakeGridPoint(cilm,lmax,vlat,vlon)
             write(*,*) i,j,val
          enddo
       enddo outer


end program