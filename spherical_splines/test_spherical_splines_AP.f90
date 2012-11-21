!===============================================================================
program test_spherical_splines
!===============================================================================
!  Create an array of knots, supply some function, then re-express that using 
!  the spherical_splines module

   use spherical_splines
   use spherical_geometry, only: sphere_sample
   
   real(8),parameter :: pi = atan2(1._8,1._8)*4._8
   
   real(8) :: dd=4. ! approx distance between knot points
   real(8) :: sigma=0.2
   real(8),allocatable,dimension(:) :: lon,lat
   integer :: n  ! number of points
   real(8),dimension(1) :: val,vlon,vcolat
   type(sph_splines_set) :: s
   
!  Create the knot points
   call sphere_sample(dd,lon,lat,n)
   
!  Create the knot points on a regular lon,lat grid
!   n = 360*(180-int(dd))/int(dd**2) + 1
!   allocate(lon(n), lat(n))
!   lon(1) = 0.  ;  lat(1) = 90.
!   vlon(1) = -180.
!   do while(vlon(1) < 180.)
!      vcolat(1) = 90. - dd
!      do while (vcolat(1) > -90.)
!         i = i + 1
!         lon(i) = vlon(1)
!         lat(i) = vcolat(1)
!         vcolat(1) = vcolat(1) - dd
!      enddo
!      vlon(1) = vlon(1) + dd
!   enddo
!   i = i + 1
!   lon(i) = 0.  ;  lat(i) = -90.
!   
!   if (i /= n) then
!      write(0,'(2(a,i0))') 'i /= n !!    i = ',i,' /= n = ',n
!      stop
!   endif
   
!  Make room for the values at each knot
   allocate(s%phi(n),s%theta(n),s%a(n),s%c(n))
   s%phi = lon
   s%theta = 90._8 - lat
   s%n = n
   
!  Fill the points with a bulge at the equator, width controlled by sigma
   s%c = exp((-(s%phi*pi/180._8)**2 - (pi/2._8 - (s%theta*pi/180._8))**2)/sigma**2)
   s%degrees = .true.

!  Alternatively, create a checkerboard
!   s%c = sign(1._8,cos(2._8*s%phi*pi/180._8)*cos(2._8*s%theta*pi/180._8))

!  Make the shape sharp sided
   where (s%c > 0.5)
      s%c = 1.
   elsewhere
      s%c = 0.
   endwhere

!  Choose a desired kernel and scaling for the Abel-Poisson 
!   s%kernel = 'AP'
!   dd = 10.
!   s%h = sph_splines_AP_width2h(dd,degrees=.true.)
   
!  Fit the splines to the data
   call sph_splines_fit_AP(s,half_width=dd)
!   write(*,*) s%a
   
!  Write out points on a regular lon,lat grid
outer:   do j=45,135,3
      do i=-50,50,3
         vlon(1) = real(i)
         vcolat(1) = real(j)
         call sph_splines_eval(vlon,vcolat,s,val,degrees=.true.)
         write(*,*) i,90-j,val(1)
      enddo
   enddo outer
   
!  Write out knot locations
   
end program
