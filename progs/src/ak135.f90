program AK135_vals

use global_1d_models

implicit none

character(len=50) :: arg
character(len=50),parameter :: fmt='(3f7.3,f6.1,f6.1)'
character(len=50),parameter :: fmt2='(3f7.3,f6.1,f6.2,f7.1)'
character(len=50),parameter :: hdr='#    vp     vs    rho     P     g'
real(8) :: depth,vp,vs,rho,P,g
integer :: iostatus

if (iargc() > 1) then
   write(0,'(a)') '  Usage: ak135 [depth / km]'
   write(0,'(a)') '     or: reads series of depths from stdin'
   stop
endif

!  Use the argument on the command line
if (iargc() == 1) then
   call getarg(1,arg)
   read(arg,*) depth
   call ak135(depth,vp=vp,vs=vs,rho=rho)
   P = pressure(depth,model='AK135')
   g = gravity(depth,model='AK135')
   write(*,'(a)') hdr
   write(*,fmt) vp,vs,rho,P,g
   stop
endif

iostatus = 0
write(*,'(2a)') trim(hdr),'  depth'
do while (iostatus == 0)
   read(*,*,iostat=iostatus) depth
   if (iostatus /= 0) stop
   call ak135(depth,vp=vp,vs=vs,rho=rho)
   P = pressure(depth,model='AK135')
   g = gravity(depth,model='AK135')
   write(*,fmt2) vp,vs,rho,P,g,depth
enddo

end program AK135_vals