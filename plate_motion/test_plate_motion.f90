program test_plate_motion
!  Find Euler pole for Africa (AF) for NUVEL1A and rotation rate at 0,60
   
   use plate_motion
   
   implicit none
   
   type(plate) :: p1,p2
   character(len=MODEL_NAME_LENGTH) :: model_name
   real(8) :: lon,lat,v(2)
   
   lon = 0._8
   lat = 60._8
   
   p1%name = "EU"
   p2%name = "AF"
   model_name = 'HS3-NUVEL1A'
   
!  Test apm calculation
   call absolute_plate_motion(lon, lat, p1, model_name, v=v)
   
   write(*,'(a,2(i0,a))') 'Absolute plate motion of Africa at ',int(lon),',',int(lat),':'
   write(*,'("   ",a,f0.2,x,f0.2)') 'E, N:    ',v  
   write(*,'("   ",a,f5.1,x,f0.2)') 'az,rate: ', &
      mod(atan2(v(1), v(2))*180._8/3.14159_8+3600._8, 360._8), &
      sqrt(v(1)**2 + v(2)**2)
      
!  Test relative plate motion calculation
   write(*,'(a,2(i0,a))') 'Relative motion of EU to AF at ',int(lon),',',int(lat),':'
   call relative_plate_motion(lon,lat,P1,P2,model_name,v=v)
   write(*,'("   ",f5.1,x,f0.2)')  mod(atan2(v(1), v(2))*180._8/3.14159_8+3600._8, 360._8), &
      sqrt(v(1)**2 + v(2)**2)

   
end program