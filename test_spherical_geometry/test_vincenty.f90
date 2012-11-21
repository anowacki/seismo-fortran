program test
   use spherical_geometry
   implicit none
   
   real(8) :: x1,x2,y1,y2,d
   
   ! Points
   x1 = 0.d0;  y1 = 0.d0
   x2 = 0.d0;  y2 = 10.d0
   
!  Calcuate distance with Haversie
   d = delta (x1,y1,x2,y2,degrees=.true.)
   write(*,*) 'Haversine distance: ',d
   
!  Calculate distance with vincenty
   call dist_vincenty(x1,y2,x2,y2,d,degrees=.true.)
   write(*,*) 'Vincenty distance:  ',d

end program
