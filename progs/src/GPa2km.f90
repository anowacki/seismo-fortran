!===============================================================================
! One of Andy Nowacki's Fortran utility modules for dealing with seismic
! anisotropy and other problems.
!
! Andy Nowacki <andy.nowacki@bristol.ac.uk>
!
! See the file LICENCE for licence details.
!===============================================================================
program GPa2km
!===============================================================================
!  Find the depth for a given Earth model which is at the requested pressure

use global_1d_models

implicit none

real(8) :: P,P0,P1,P2,depth,depth0,depth1,depth2,dPdz
character(10) :: model
real(8),parameter :: tol=5.  !  Accuracy of the solution / GPa
real(8),parameter :: dz=100.

if (command_argument_count() /= 1 .and. command_argument_count() /= 2) then
   write(0,'(a)') '  Usage: GPa2km [P / GPa] (PREM|AK135)'
   stop
endif

call get_command_argument(1,model); read(model,*) P
model = 'AK135'
if (command_argument_count() == 2) call get_command_argument(2,model)

!  Search for the depth which gives the correct pressure using Newton-Raphson
!  f(x)=0: x_n+1 = x_n - f(x_n) / f'(x_n)
!  f'(x) approx. by (f(x_n)-f(x_n+dx))/dx

!depth1 = 0.
!depth0 = 2890.  ! starting guess of CMB
!do while (abs(depth1 - depth0) > tol)
!   if (depth1 > 0.) depth0 = depth1
!   if (depth1 < 0.) depth0 = 0.
!   dPdz = (pressure(depth0+dz,model=trim(model)) - &
!           pressure(depth0,model=trim(model)))/dz
!   depth1 = depth0 - (pressure(depth,model=trim(model)) - P) / dPdz
!   write(*,*) depth1
!enddo

!  Actually, use the method of halves to get the correct depth:
write(0,'(a)',advance='no') 'GPa2km is thinking'
!  starting guesses
depth0 = 0.      ;  P0 = 0.    ! Top
depth2 = 6371.   ;  P2 = 365.  ! Bottom
depth1 = 6371./2.;  P1 = pressure(depth1,model=trim(model))  ! Middle
do while (abs(P2 - P0) > tol)
   if (P > P0 .and. P < P1) then
      depth0 = depth0
      depth2 = depth1                  ; P2 = pressure(depth2,model=trim(model))
      depth1 = (depth0 + depth2) / 2.  ; P1 = pressure(depth1,model=trim(model))
   else if (P > P1 .and. P < P2) then
      depth0 = depth1                  ; P0 = pressure(depth0,model=trim(model))
      depth2 = depth2
      depth1 = (depth0 + depth2) / 2.  ; P1 = pressure(depth1,model=trim(model))
   else if (P == P1) then
      write(*,*) depth1; stop
   endif
   P0 = pressure(depth0,model=trim(model))
   P1 = pressure(depth1,model=trim(model))
   write(0,'(a)',advance='no') '.'
enddo
write(0,'(a)') ''
write(*,'(i4)') int(depth1)

end program GPa2km
!-------------------------------------------------------------------------------