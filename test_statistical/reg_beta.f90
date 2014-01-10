program reg_beta
use statistical
implicit none
integer,parameter :: rs=8
real(rs) :: z,a,b
character(len=80) :: arg

if (command_argument_count() /= 3) then
   write(*,'(a)') 'Usage: reg_beta [z] [a] [b]'
   stop
endif

call get_command_argument(1,arg) ; read(arg,*) z
call get_command_argument(2,arg) ; read(arg,*) a
call get_command_argument(3,arg) ; read(arg,*) b

write(*,*) reg_beta_func(z,a,b)

end program
