program plate_name
!  Tell us what the full name of the plate is

use plate_motion

implicit none

character(len=PLATE_NAME_LENGTH) :: ab
character(len=LONG_PLATE_NAME_LENGTH) :: name
character(len=80) :: arg

if (command_argument_count() /= 1) then
   write(0,'(a)') 'Usage: plate_name [plate two-letter abbreviation]'
   stop
endif

call get_command_argument(1,arg)

if (len_trim(arg) /= PLATE_NAME_LENGTH) then
   write(0,'(a)') 'plate_name: Plate abbreviation must be two letters long.'
   stop
endif

ab = arg(1:PLATE_NAME_LENGTH)

name = get_plate_name(ab)


write(*,'(a)') trim(name)

end program plate_name