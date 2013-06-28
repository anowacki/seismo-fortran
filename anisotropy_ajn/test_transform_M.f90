!===============================================================================
program test_transform_M
!===============================================================================
! Make sure that the CIJ_transform_M routines give the same results.
! Also compare them in terms of speed

   use anisotropy_ajn

   implicit none
   
   integer,parameter :: rs = 8
   real(rs), dimension(6,6) :: C,CT,Cflip,Crot
   real(rs) :: M(3,3),a
   real(rs) :: rh
   real(rs), parameter :: tol = 1.d-5
   character(len=250) :: file_orig,file_flipped,file_rot

   if (command_argument_count() /= 3) then
      write(*,'(a)') 'Usage: test_flip [original ecs .ecs] [flipped ecs .ecs] [rotated ecs .ecs]'
      stop
   endif

   call get_command_argument(1,file_orig)
   call get_command_argument(2,file_flipped)
   call get_command_argument(3,file_rot)

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

   ! Flip across the x=0 plane
   Cflip = CIJ_flipx(C)
   call CIJ_save(trim(file_flipped)//'.ref',Cflip,rh)

   ! Do it the Bower way
   M = 0._rs; M(1,1) = -1._rs; M(2,2) = 1._rs; M(3,3) = 1._rs
   CT = CIJ_transform_M(C,M)
   write(*,'(a)',advance='no') 'Explicit v Bower flip: '
   call CIJ_save(trim(file_flipped)//'.2',CT,rh)
   call check(Cflip,CT,tol)

   ! Rotate about the x axis 30 degrees
   a = 30._rs
   call CIJ_rot3_old(C,a,0._rs,0._rs,Crot)
   call CIJ_save(file_rot,Crot,rh)

   ! Do it the Bower way
   call CIJ_rot3(C,a,0._rs,0._rs,CT)
   call CIJ_save(trim(file_rot)//'.2',CT,rh)
   write(*,'(a)',advance='no') 'Old v Bower rotation: '
   call check(Crot,CT,tol)

end program test_transform_M
!-------------------------------------------------------------------------------

!===============================================================================
subroutine check(A,B,tol)
!===============================================================================
   implicit none
   real(8), intent(in) :: A(6,6),B(6,6),tol
   if (any(2._8*abs(A-B)/(A+B) > tol)) then
      write(*,'(a,e7.1)') 'Two matrices do not agree to +/- ',tol
   else
      write(*,'(a)') 'Two matrices agree'
   endif
end subroutine check
!-------------------------------------------------------------------------------

