!===============================================================================
! One of Andy Nowacki's Fortran utility modules for dealing with seismic
! anisotropy and other problems.
!
! Andy Nowacki <andy.nowacki@bristol.ac.uk>
!
! See the file LICENCE for licence details.
!===============================================================================
program CIJ_rot_azinc
!===============================================================================
!  Program calling CIJ_rot_azinc

   use anisotropy_ajn
   
   implicit none
   
   real(8) :: C(6,6), CR(6,6), az, inc, phi, rho
   character(len=250) :: file,arg
   integer :: iostatus,i,j
   
!  Check input arguments
   if (command_argument_count() < 3 .or. command_argument_count() > 4) then
      write(0,'(a)') &
         'Usage: CIJ_rot_azinc [az] [inc] [phi] (ecfile)', &
         '  Rotations are applied about vector along (az, inc) by phi, anticlock. looking down axis.', &
         '  Azimuth is from +x1 towards -x2; inc is from x1-x2 towards +x3.', &
         '  If no .ecs file specified, then ECs are read from stdin, 36 constants per line.'
      stop
   endif
   
!  Get input arguments
   call get_command_argument(1,arg) ;  read(arg,*) az
   call get_command_argument(2,arg) ;  read(arg,*) inc
   call get_command_argument(3,arg) ;  read(arg,*) phi
   
!  Get elastic constants
!  if reading from .ecs file, just output in same units: assumes density normalised
   if (command_argument_count() == 4) then
      call get_command_argument(4,file)
      call CIJ_load(file,C,rho)
      CR = CIJ_rotate_az_inc(C,az,inc,phi)
      write(*,*) CR
!  Otherwise, we're reading several constants from stdin
   else
      iostatus = 0
      do while (iostatus == 0)
         read(*,*,iostat=iostatus) ((C(i,j),j=1,6),i=1,6)
         if (iostatus < 0) exit
         if (iostatus > 0) then
            write(0,'(a)') 'CIJ_rot_azinc: problem reading 36 input elastic constants from stdin.'
            stop
         endif
         CR = CIJ_rotate_az_inc(C,az,inc,phi)
         write(*,*) CR
      enddo
   endif

end program CIJ_rot_azinc
!-------------------------------------------------------------------------------
