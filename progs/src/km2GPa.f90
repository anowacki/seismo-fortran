program km2GPa

use global_1d_models

implicit none

real(8) :: z
character(len=10) :: model

if (command_argument_count() /= 1 .and. command_argument_count() /= 2) then
   write(0,'(a)') '  Usage: km2GPa [depth /km] (PREM|AK135)'
   stop
endif

call get_command_argument(1,model); read(model,*) z
model = 'AK135'
if (command_argument_count() == 2) call get_command_argument(2,model)

write(*,'(f7.2)') pressure(z,model=trim(model))


end program km2GPa