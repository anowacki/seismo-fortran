!===============================================================================
program test_flip
!===============================================================================
! Make sure that the CIJ_flip? subroutines give expected results.  This is done
! by visual comparison using the test_flip.sh which uses CIJ_plot.

   use anisotropy_ajn

   implicit none
   
   real(8), dimension(6,6) :: C,Cf
   real(8) :: rh
   character(len=250) :: file_orig,file_flipped,arg
   
   if (command_argument_count() /= 3) then
      write(*,'(a)') 'Usage: test_flip [x|y|z] [original ecs .ecs] [flipped ecs .ecs]'
      stop
   endif
   
   call get_command_argument(1,arg)
   if (arg /= 'x' .and. arg /= 'y' .and. arg /= 'z') then
      write(*,'(a)') 'Usage: test_flip [x|y|z] [original ecs .ecs] [flipped ecs .ecs]'
      stop
   endif

   call get_command_argument(2,file_orig)
   call get_command_argument(3,file_flipped)
   
   ! An average kimberlite from Mainprice and Silver, PEPI, 1993
   C(1,1) = 230.88e9
   C(1,2) = 85.14e9;  C(2,1) = 85.14e9
   C(1,3) = 88.53e9;  C(3,1) = 88.53e9
   C(1,4) = 0.00e9;   C(4,1) = 0.00e9
   C(1,5) = 0.88e9;   C(5,1) = 0.88e9
   C(1,6) = -1.29e9;  C(6,1) = -1.29e9
   C(2,2) = 252.33e9
   C(2,3) = 85.93e9;  C(3,2) = 85.93e9
   C(2,4) = 0.41e9;   C(4,2) = 0.41e9
   C(2,5) = 0.56e9;   C(5,2) = 0.56e9
   C(2,6) = -1.42e9;  C(6,2) = -1.42e9
   C(3,3) = 237.65e9
   C(3,4) = -1.14e9;  C(4,3) = -1.14e9
   C(3,5) = 1.11e9;   C(5,3) = 1.11e9
   C(3,6) = -0.01e9;  C(6,3) = -0.01e9
   C(4,4) = 76.97e9
   C(4,5) = -0.68e9;  C(5,4) = -0.68e9
   C(4,6) = 0.68e9;   C(6,4) = 0.68e9
   C(5,5) = 72.49e9
   C(5,6) = -0.34e9;  C(6,5) = -0.34e9
   C(6,6) = 74.75e9
   rh = 3333.0
   
   ! Save original
   call CIJ_save(file_orig,C,rh)
   
   ! Flip
   if (arg == 'x') then
      Cf = CIJ_flipx(C)
   else if (arg == 'y') then
      Cf = CIJ_flipy(C)
   else
      Cf = CIJ_flipz(C)
   endif
   
   ! Save flipped version
   call CIJ_save(file_flipped,Cf,rh)

end program test_flip
!-------------------------------------------------------------------------------