!===============================================================================
program test_rot_euler
!===============================================================================
! Make sure that the CIJ_rot_euler routines work as expected.

   use anisotropy_ajn

   implicit none
   
   integer,parameter :: rs = 8
   real(rs), dimension(6,6) :: C, Crot
   real(rs) :: psi1, phi, psi2, rh
   character(len=250) :: file_orig = '/tmp/test_rot_euler.orig', &
                         file_rot  = '/tmp/test_rot_euler.rot', arg

   if (command_argument_count() /= 3) then
      write(*,'(a)') 'Usage: test_rot_euler [psi1] [phi] [psi2]'
      stop
   endif

   call get_command_argument(1,arg)
   read(arg,*) psi1
   call get_command_argument(2,arg)
   read(arg,*) phi
   call get_command_argument(3,arg)
   read(arg,*) psi2
   
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
   call CIJ_save(file_orig, C, rh)

   ! Rotate
   Crot = CIJ_rot_euler(c, psi1, phi, psi2, passive=.false., type='z1x2z3')
   
   ! Write out new one
   call CIJ_save(file_rot, Crot, rh)

end program test_rot_euler
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

