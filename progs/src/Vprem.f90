program PREM_vals

use global_1d_models

implicit none

character(len=50) :: arg
character(len=50),parameter :: fmt='(3f7.3)'
character(len=50),parameter :: fmt2='(3f7.3,f7.1)'
character(len=50),parameter :: hdr='#    vp     vs    rho'
real(8) :: depth,vp,vs,rho,P,g
integer :: iostatus

if (iargc() > 1) then
   write(0,'(a)') '  Usage: Vprem [depth / km]'
   write(0,'(a)') '     or: reads series of depths from stdin'
   stop
endif

!  Use the argument on the command line
if (iargc() == 1) then
   call getarg(1,arg)
   read(arg,*) depth
   call prem(depth,vp=vp,vs=vs,rho=rho)
   write(*,'(a)') hdr
   write(*,fmt) vp,vs,rho
   stop
endif

iostatus = 0
write(*,'(2a)') trim(hdr),'  depth'
do while (iostatus == 0)
   read(*,*,iostat=iostatus) depth
   if (iostatus /= 0) stop
   call prem(depth,vp=vp,vs=vs,rho=rho)
   write(*,fmt2) vp,vs,rho,depth
enddo

end program PREM_vals