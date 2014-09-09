!===============================================================================
! One of Andy Nowacki's Fortran utility modules for dealing with seismic
! anisotropy and other problems.
!
! Andy Nowacki <andy.nowacki@bristol.ac.uk>
!
! See the file LICENCE for licence details.
!===============================================================================
program spread
!===============================================================================
!  ...sounds rude but calculates the relative plate motion between two plates
!  and is envisaged for calculating spreading rates and directions at divergent
!  plate boundaries.
!
!  Can be used for one-shot calculations like:
!     spread [lon] [lat] [plate1] [plate2] (model)
!
!  or can accept a list of points and plates from stdin like this:
!     echo -e "0 0 AF SA \n 10 10 EU AF" | spread

   use plate_motion
   
   implicit none
   
   type(plate) :: P1,P2
   character(len=MODEL_NAME_LENGTH) :: model_name
   character(len=80) :: arg
   real(8) :: lon,lat,az,rate
   integer :: iostatus
   
!  Check arguments
   if (command_argument_count() /= 0 .and. command_argument_count() /= 4 .and. command_argument_count() /= 5) then
      write(0,'(a)') 'Usage: spread [lon] [lat] [plate1] [plate2] (model)', &
                     '       For full list of plates and models, use command "plate_name_list"'
      stop
   endif
   
!  If requested, get model; otherwise assume HS3-NUVEL1A
   model_name = 'HS3-NUVEL1A'
   if (command_argument_count() == 5) call get_command_argument(5,model_name)

!  Get arguments if we're doing a one-shot calculation
   if (command_argument_count() /= 0) then
      call get_command_argument(1,arg) ;  read(arg,*) lon
      call get_command_argument(2,arg) ;  read(arg,*) lat
      call get_command_argument(3,P1%name)
      call get_command_argument(4,P2%name)
   
!  Calculate relative motion of plate 1 wrt plate 2 and write results
      call relative_plate_motion(lon,lat,P1,P2,model_name=model_name,az=az,rate=rate)
      write(*,'(a)') '# Az(N->E)  rate(mm/a):'
      write(*,'(f0.2,x,f0.3)') az,rate
      stop
      
   else
!  Otherwise we're reading from stdin, so loop over the input.  This version only 
!  accepts HS3-NUVEL1A calculations at the moment.
      write(*,'(3a)') '#   lon    lat P1 P2 azi   rate(mm/a) [',trim(model_name),']'
      iostatus = 0
      do while(iostatus == 0)
         read(*,*,iostat=iostatus) lon,lat,p1%name,p2%name
         if (iostatus < 0) exit  ! EOF
         if (iostatus > 0) then  ! Read error
            write(0,'(a)') 'spread: Problem reading from stdin'
            stop
         endif
         call relative_plate_motion(lon,lat,P1,P2,model_name=model_name,az=az,rate=rate)
         write(*,'(f7.2,x,f6.2,x,a2,x,a2,x,f0.2,x,f0.3)') lon,lat,p1%name,p2%name,az,rate
      enddo
      
   endif

end program spread