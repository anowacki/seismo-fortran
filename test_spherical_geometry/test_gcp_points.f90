program test
! Check that gcp_points works as expected

use spherical_geometry

integer,parameter :: rs = 8, &
                     npts = 1000 ! Number of points along path
real(rs) :: x1 =  0.,  y1 =  0.,  &
            x2 = 45.,  y2 = 45.  ! Endpoints
real(rs) :: ds, lon(npts), lat(npts) ! Points on path

write(0,'(a,f0.1,1x,f0.1)') 'Start point : ',x1,y1,&
                            'End point :   ',x2,y2

call gcp_points(x1,y1,x2,y2,lon,lat,n=npts,degrees=.true.)

do i=1,npts
   write(*,*) lon(i),lat(i)
enddo

end program test
