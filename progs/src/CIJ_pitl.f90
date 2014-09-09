!===============================================================================
! One of Andy Nowacki's Fortran utility modules for dealing with seismic
! anisotropy and other problems.
!
! Andy Nowacki <andy.nowacki@bristol.ac.uk>
!
! See the file LICENCE for licence details.
!===============================================================================
program pitl
!===============================================================================
!  Wrapper to the CIJ_pitl routine.
!  The program computes the effective long-wavelength anisotropic equivalent
!  constants for a periodic, isotropic thin-layered medium, where two alternating
!  layers of constant relative thickness d1/(d1+d2) and d2/(d1+d2) are present.

   use anisotropy_ajn, only: CIJ_pitl

   implicit none

   integer, parameter :: rs = 8
   real(rs) :: C(6,6), vp1, vp2, vs1, vs2, d1, d2, rho1, rho2
   character(len=250) :: arg
   integer :: iostat, narg

   ! Get command line arguments if any
   narg = command_argument_count()
   if (narg /= 0 .and. narg /= 8) call usage

   if (narg == 0) then
      iostat = 0
      do while (iostat == 0)
         read(*,*,iostat=iostat) d1, vp1, vs1, rho1, d2, vp2, vs2, rho2
         if (iostat < 0) exit
         if (iostat > 0) then
            write(0,'(a)') 'CIJ_pitl: Error: Can''t get d1,Vp1,Vs1,rho1,d2,Vp2,Vs2,rho2 from stdin'
            stop 1
         endif
         call get_ptl_vals
         call output
      enddo

   else
      call get_arg(1, 'd1',   d1)
      call get_arg(2, 'vp1',  vp1)
      call get_arg(3, 'vs1',  vs1)
      call get_arg(4, 'rho1', rho1)
      call get_arg(5, 'd2',   d2)
      call get_arg(6, 'vp2',  vp2)
      call get_arg(7, 'vs2',  vs2)
      call get_arg(8, 'rho2', rho2)
      call get_ptl_vals
      call output
   endif

contains
   subroutine get_ptl_vals
      C = CIJ_pitl(d1, vp1, vs1, rho1, d2, vp2, vs2, rho2)
   end subroutine get_ptl_vals

   subroutine usage
      write(0,'(a)') &
         'Usage CIJ_pitl [d1] [Vp1] [Vs1] [rho1] [d2] [Vp2] [Vs2] [rho2]', &
         '      CIJ_pitl < ([d1] [Vp1] [Vs1] [rho1] [d2] [Vp2] [Vs2] [rho2] on stdin)', &
         'Calculate the long-wavelength equivalent elastic constants for a medium', &
         'made up of periodic layering of two isotropic layers, 1 and 2, where the', &
         'relative widths of the layers are d1 and d2.', &
         'Vp, Vs and rho are the P- and S-wave velocities and densities, respectively.', &
         'Uses formulae as given by Postma (1955), Geophysics, 20, 780-806.'
      stop 1
   end subroutine usage

   subroutine get_arg(i, name, var)
      integer, intent(in) :: i
      character(len=*), intent(in) :: name
      real(rs), intent(out) :: var
      call get_command_argument(i, arg)
      read(arg,*,iostat=iostat) var
      if (iostat /= 0) then
         write(0,'(a)') 'CIJ_pitl: Error: Can''t read value of '//trim(name) &
            //'from argument "'//trim(arg)//'"'
         stop 1
      endif
   end subroutine get_arg

   subroutine output
      write(*,*) C
   end subroutine output
end program pitl
!-------------------------------------------------------------------------------