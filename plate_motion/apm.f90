!===============================================================================
program apm
!===============================================================================
!  Calculates the aboslute plate motion at a given spot for a given plate.
!  Resposibility to pick a plate at that point is up to the user.

   use plate_motion
   
   implicit none
   
   type(plate) :: p
   character(len=MODEL_NAME_LENGTH) :: model_name
   character(len=80) :: arg
   real(8) :: lon,lat,rate,az
   integer :: iostatus
   
!  Check arguments
   if (iargc() /= 0 .and. iargc() /= 3 .and. iargc() /= 4) then
      write(0,'(a)') 'Usage: apm [lon] [lat] [plate] (model)', &
                     '       For full list of plates and models, use command "plate_name_list"'
      stop
   endif
   
!  If requested, get model; otherwise assume HS3-NUVEL1A
   model_name = 'HS3-NUVEL1A'
   if (iargc() == 4) call getarg(4,model_name)

!  Get arguments if we're doing a one-shot calculation
   if (iargc() /= 0) then
      call getarg(1,arg) ;  read(arg,*) lon
      call getarg(2,arg) ;  read(arg,*) lat
      call getarg(3,P%name)
   
!  Calculate relative motion of plate 1 wrt plate 2 and write results
      call absolute_plate_motion(lon,lat,P,model_name=model_name,az=az,rate=rate)
      write(*,'(a)') '# Az(N->E)  rate(mm/a):'
      write(*,'(f0.2,x,f0.3)') az,rate
      stop
      
   else
!  Otherwise we're reading from stdin, so loop over the input.  This version only 
!  accepts HS3-NUVEL1A calculations at the moment.
      write(*,'(3a)') '#   lon    lat plate  azi   rate(mm/a) [',trim(model_name),']'
      iostatus = 0
      do while(iostatus == 0)
         read(*,*,iostat=iostatus) lon,lat,p%name
         if (iostatus < 0) exit  ! EOF
         if (iostatus > 0) then  ! Read error
            write(0,'(a)') 'spread: Problem reading from stdin'
            stop
         endif
         call absolute_plate_motion(lon,lat,P,model_name=model_name,az=az,rate=rate)
         write(*,'(f7.2,x,f6.2,x,a2,x,f0.2,x,f0.3)') lon,lat,p%name,az,rate
      enddo
      
   endif

end program apm