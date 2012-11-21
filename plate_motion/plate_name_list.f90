!===============================================================================
program plate_name_list
!===============================================================================
!  Lists all the plates available to the plate_motion module, as well as
!  all the models.

   use plate_motion
   
   implicit none
      
!  Just print everything out, regardless of how the program is run
   call plate_list_all
   
end program