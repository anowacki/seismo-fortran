!===============================================================================
module get_args
!===============================================================================
!  Andy Nowacki, University of Bristol
!  andy.nowacki@bristol.ac.uk
!
!  get_args is a module for reading the command line options used when running
!  a Fortran program.  It uses the Fortran 2003 standard call get_command_arguments()
!  to read the command line.  If you don't have a modern Fortran compiler, then
!  you can use the call get_command_argument() and iargc() commands.
!  
!  get_args takes a list of argument switch names (e.g., -r, -f) and the expected
!  type  and number of associated arguments and passes back either an unchanged
!  value or that given on the command line.  A module procedure interface allows
!  real (double or single precision), integer or character variables to be called
!  with the same call to the get_arg() subroutine.
!
!  You can specify that the argument must be present to continue running your program
!  with the optional argument required=.true. in the call to get_arg().
!  Supplying a logical variable with 'supplied=value_given' will allow you to tell
!  if a parameter was given on the command line.
!
!  Examples:
!     You want to optionally set the minimum and maximum search range for a variable
!     with the command line flag: '-range [min] [max]'.
!        min = 0.0;  max = 1.0
!        call get_arg('-range',min,max)
!     
!     The user must specify the integer number of something with '-n [number]'.
!        call get_arg('-n',number,required=.true.)
!     
!     An input file name character may be specified to overwrite the default, but
!     we want to know if the default has been changed.
!        input_file = 'input.text'
!        call get_arg('-input',input_file,supplied=input_file_changed)
!
!-------------------------------------------------------------------------------
!  Known limitations and bugs:
!   - There is no way to have a switch turn something on (or off) without giving
!     at least one argument (e.g., you can't "call get_arg('-nobugs')").  You
!     have to supply a dummy argument on the command line.  However, for these 
!     situations you can use "./myprog -bugs 0" or "./myprog -bugs 1" and test
!     for the integer value.
!     (This limitation arises from the module procedure feature of Fortran.  
!     However, if you're using this module at all and would rather things worked 
!     differently, you probably shouldn't be using Fortran to process command 
!     line arguments!)
!   - get_args has no way to check that you aren't supplying a typo of a switch,
!     and will silently ignore the presence of a non-defined switch name.  A
!     possible workaround is a check_args routine supplied with a list of 
!     acceptable switches and their number of arguments.
!   - Arguments to switches have to all be of the same type.
!   - Alternative (e.g., short and long) switches for the same variables have
!     to be set separately in calls.
!   - At the moment, the maxmimum number of arguments to each switch is 9.  This
!     can be easily expanded with some tedious coding.
!
!===============================================================================
   
implicit none
   
   integer,parameter,private :: max_len = 250
   private
   public :: get_arg
   
   interface get_arg
      module procedure get_arg_real8
      module procedure get_arg_real4
      module procedure get_arg_int
      module procedure get_arg_char
   end interface get_arg
   
CONTAINS

!===============================================================================
subroutine get_arg_real8(name,v1,v2,v3,v4,v5,v6,v7,v8,v9,required,supplied)
!===============================================================================
!  Arguments:
!     name : name of switch
!     v1   : value to assign
!     v2,v3,v4,v5,v6,v7,v8,v9 : other optional values
!     required: logical.  If true, program exits if argument isn't supplied
!     supplied: logical.  True if we get a value for this option.

   implicit none
   character(len=*) :: name
   real(8),intent(inout) :: v1
   real(8),intent(inout),optional :: v2,v3,v4,v5,v6,v7,v8,v9
   logical,intent(in),optional :: required
   logical,intent(out),optional :: supplied
   real(8),dimension(9) :: varray
   logical :: strict,p2,p3,p4,p5,p6,p7,p8,p9
   integer :: nargs,iarg,jarg,n,length
   character(max_len) :: arg
   
!  Decide if this option must be present
   strict = .false.
   if (present(required)) strict = required

!  Check how many options required
   p2 = present(v2); p3 = present(v3); p4 = present(v4); p5 = present(v5)
   p6 = present(v6); p7 = present(v7); p8 = present(v8); p9 = present(v9)
   if (.not.p2.and..not.p3.and..not.p4.and..not.p5.and..not.p6.and..not.p7.and..not.p8.and..not.p9) then
      nargs = 1
   else if (p2.and..not.p3.and..not.p4.and..not.p5.and..not.p6.and..not.p7.and..not.p8.and..not.p9) then
      nargs = 2
   else if(p2.and.p3.and..not.p4.and..not.p5.and..not.p6.and..not.p7.and..not.p8.and..not.p9) then
      nargs = 3
   else if(p2.and.p3.and.p4.and..not.p5.and..not.p6.and..not.p7.and..not.p8.and..not.p9) then
      nargs = 4
   else if(p2.and.p3.and.p4.and.p5.and..not.p6.and..not.p7.and..not.p8.and..not.p9) then
      nargs = 5
   else if(p2.and.p3.and.p4.and.p5.and.p6.and..not.p7.and..not.p8.and..not.p9) then
      nargs = 6
   else if(p2.and.p3.and.p4.and.p5.and.p6.and.p7.and..not.p8.and..not.p9) then
      nargs = 7
   else if(p2.and.p3.and.p4.and.p5.and.p6.and.p7.and.p8.and..not.p9) then
      nargs = 8
   else if(p2.and.p3.and.p4.and.p5.and.p6.and.p7.and.p8.and.p9) then
      nargs = 9
   else
      write(0,'(a)') 'get_arg_real: supply values to get sequentially'
      stop 1
   endif
   
!  If no arguments supplied, return
   n = command_argument_count()
   if (n == 0) then
      if (strict) then
         write(0,'(a)') 'Option "'//trim(name)//'" is required.'
         stop 1
      endif
      if (present(supplied)) supplied = .false.
      return
   endif
   
!  Check there are enough arguments
   if (n < nargs + 1) then
      write(0,'(a)') 'Not enough arguments supplied for option "'//trim(name)//'".'
      stop 1
   endif
   
!  Loop over the arguments until we find the one we're after
   iarg = 1
   do while (iarg <= n)
      call get_command_argument(iarg,arg,length)
      if (length > max_len) then
         write(0,'(a)') 'Argument to option "'//trim(name)//'" is too long for get_args module.'
         stop 1
      endif
      if (name == arg) then
         !  Check for enough arguments
         if (n >= iarg + nargs) then
            do jarg=1,nargs
               call get_command_argument(iarg+jarg,arg)
               read(arg,*) varray(jarg)
            enddo
         else
            write(0,'(a)') 'Not enough arguments supplied for option "'//trim(name)//'".'
            stop 1
         endif
         v1 = varray(1)
         if (p2) v2 = varray(2)
         if (p3) v3 = varray(3)
         if (p4) v4 = varray(4)
         if (p5) v5 = varray(5)
         if (p6) v6 = varray(6)
         if (p7) v7 = varray(7)
         if (p8) v8 = varray(8)
         if (p9) v9 = varray(9)
         if (present(supplied)) supplied = .true.
         return
      endif
      iarg = iarg + 1
   enddo
   
!  Getting this far means we haven't found the requested argument
   if (present(supplied)) supplied = .false.
   
!  Exit program is we haven't found the argument we're after
   if (strict) then
      write(0,'(a)') 'Option "'//trim(name)//'" is required.'
      stop 1
   endif

end subroutine get_arg_real8
!-------------------------------------------------------------------------------

!===============================================================================
subroutine get_arg_real4(name,v1,v2,v3,v4,v5,v6,v7,v8,v9,required,supplied)
!===============================================================================
!  Arguments:
!     name : name of switch
!     v1   : value to assign
!     v2,v3,v4,v5,v6,v7,v8,v9 : other optional values
!     required: logical.  If true, program exits if argument isn't supplied
!     supplied: logical.  True if we get a value for this option.

   implicit none
   character(len=*) :: name
   real(4),intent(inout) :: v1
   real(4),intent(inout),optional :: v2,v3,v4,v5,v6,v7,v8,v9
   logical,intent(in),optional :: required
   logical,intent(out),optional :: supplied
   real(4),dimension(9) :: varray
   logical :: strict,p2,p3,p4,p5,p6,p7,p8,p9
   integer :: nargs,iarg,jarg,n,length
   character(max_len) :: arg
   
!  Decide if this option must be present
   strict = .false.
   if (present(required)) strict = required

!  Check how many options required
   p2 = present(v2); p3 = present(v3); p4 = present(v4); p5 = present(v5)
   p6 = present(v6); p7 = present(v7); p8 = present(v8); p9 = present(v9)
   if (.not.p2.and..not.p3.and..not.p4.and..not.p5.and..not.p6.and..not.p7.and..not.p8.and..not.p9) then
      nargs = 1
   else if (p2.and..not.p3.and..not.p4.and..not.p5.and..not.p6.and..not.p7.and..not.p8.and..not.p9) then
      nargs = 2
   else if(p2.and.p3.and..not.p4.and..not.p5.and..not.p6.and..not.p7.and..not.p8.and..not.p9) then
      nargs = 3
   else if(p2.and.p3.and.p4.and..not.p5.and..not.p6.and..not.p7.and..not.p8.and..not.p9) then
      nargs = 4
   else if(p2.and.p3.and.p4.and.p5.and..not.p6.and..not.p7.and..not.p8.and..not.p9) then
      nargs = 5
   else if(p2.and.p3.and.p4.and.p5.and.p6.and..not.p7.and..not.p8.and..not.p9) then
      nargs = 6
   else if(p2.and.p3.and.p4.and.p5.and.p6.and.p7.and..not.p8.and..not.p9) then
      nargs = 7
   else if(p2.and.p3.and.p4.and.p5.and.p6.and.p7.and.p8.and..not.p9) then
      nargs = 8
   else if(p2.and.p3.and.p4.and.p5.and.p6.and.p7.and.p8.and.p9) then
      nargs = 9
   else
      write(0,'(a)') 'get_arg_real: supply values to get sequentially'
      stop 1
   endif
   
!  If no arguments supplied, return
   n = command_argument_count()
   if (n == 0) then
      if (strict) then
         write(0,'(a)') 'Option "'//trim(name)//'" is required.'
         stop 1
      endif
      if (present(supplied)) supplied = .false.
      return
   endif
   
!  Check there are enough arguments
   if (n < nargs + 1) then
      write(0,'(a)') 'Not enough arguments supplied for option "'//trim(name)//'".'
      stop 1
   endif
   
!  Loop over the arguments until we find the one we're after
   iarg = 1
   do while (iarg <= n)
      call get_command_argument(iarg,arg,length)
      if (length > max_len) then
         write(0,'(a)') 'Argument to option "'//trim(name)//'" is too long for get_args module.'
         stop 1
      endif
      if (name == arg) then
         !  Check for enough arguments
         if (n >= iarg + nargs) then
            do jarg=1,nargs
               call get_command_argument(iarg+jarg,arg)
               read(arg,*) varray(jarg)
            enddo
         else
            write(0,'(a)') 'Not enough arguments supplied for option "'//trim(name)//'".'
            stop 1
         endif
         v1 = varray(1)
         if (p2) v2 = varray(2)
         if (p3) v3 = varray(3)
         if (p4) v4 = varray(4)
         if (p5) v5 = varray(5)
         if (p6) v6 = varray(6)
         if (p7) v7 = varray(7)
         if (p8) v8 = varray(8)
         if (p9) v9 = varray(9)
         if (present(supplied)) supplied = .true.
         return
      endif
      iarg = iarg + 1
   enddo
   
!  Getting this far means we haven't found the requested argument
   if (present(supplied)) supplied = .false.
   
!  Exit program is we haven't found the argument we're after
   if (strict) then
      write(0,'(a)') 'Option "'//trim(name)//'" is required.'
      stop 1
   endif

end subroutine get_arg_real4
!-------------------------------------------------------------------------------

!===============================================================================
subroutine get_arg_int(name,v1,v2,v3,v4,v5,v6,v7,v8,v9,required,supplied)
!===============================================================================
!  Arguments:
!     name : name of switch
!     v1   : value to assign
!     v2,v3,v4,v5,v6,v7,v8,v9 : other optional values
!     required: logical.  If true, program exits if argument isn't supplied
!     supplied: logical.  True if we get a value for this option.

   implicit none

   character(len=*) :: name
   integer,intent(inout) :: v1
   integer,intent(inout),optional :: v2,v3,v4,v5,v6,v7,v8,v9
   logical,intent(in),optional :: required
   logical,intent(out),optional :: supplied
   integer,dimension(9) :: varray
   logical :: strict,p2,p3,p4,p5,p6,p7,p8,p9
   integer :: nargs,iarg,jarg,n,length
   character(max_len) :: arg
   
!  Decide if this option must be present
   strict = .false.
   if (present(required)) strict = required

!  Check how many options required
   p2 = present(v2); p3 = present(v3); p4 = present(v4); p5 = present(v5)
   p6 = present(v6); p7 = present(v7); p8 = present(v8); p9 = present(v9)
   if (.not.p2.and..not.p3.and..not.p4.and..not.p5.and..not.p6.and..not.p7.and..not.p8.and..not.p9) then
      nargs = 1
   else if (p2.and..not.p3.and..not.p4.and..not.p5.and..not.p6.and..not.p7.and..not.p8.and..not.p9) then
      nargs = 2
   else if(p2.and.p3.and..not.p4.and..not.p5.and..not.p6.and..not.p7.and..not.p8.and..not.p9) then
      nargs = 3
   else if(p2.and.p3.and.p4.and..not.p5.and..not.p6.and..not.p7.and..not.p8.and..not.p9) then
      nargs = 4
   else if(p2.and.p3.and.p4.and.p5.and..not.p6.and..not.p7.and..not.p8.and..not.p9) then
      nargs = 5
   else if(p2.and.p3.and.p4.and.p5.and.p6.and..not.p7.and..not.p8.and..not.p9) then
      nargs = 6
   else if(p2.and.p3.and.p4.and.p5.and.p6.and.p7.and..not.p8.and..not.p9) then
      nargs = 7
   else if(p2.and.p3.and.p4.and.p5.and.p6.and.p7.and.p8.and..not.p9) then
      nargs = 8
   else if(p2.and.p3.and.p4.and.p5.and.p6.and.p7.and.p8.and.p9) then
      nargs = 9
   else
      write(0,'(a)') 'get_arg_real: supply values to get sequentially'
      stop 1
   endif
   
!  If no arguments supplied, return
   n = command_argument_count()
   if (n == 0) then
      if (strict) then
         write(0,'(a)') 'Option "'//trim(name)//'" is required.'
         stop 1
      endif
      if (present(supplied)) supplied = .false.
      return
   endif
   
!  Check there are enough arguments
   if (n < nargs + 1) then
      write(0,'(a)') 'Not enough arguments supplied for option "'//trim(name)//'".'
      stop 1
   endif
   
!  Loop over the arguments until we find the one we're after
   iarg = 1
   do while (iarg <= n)
      call get_command_argument(iarg,arg,length)
      if (length > max_len) then
         write(0,'(a)') 'Argument to option "'//trim(name)//'" is too long for get_args module.'
         stop 1
      endif
      if (name == arg) then
         !  Check for enough arguments
         if (n >= iarg + nargs) then
            do jarg=1,nargs
               call get_command_argument(iarg+jarg,arg)
               read(arg,*) varray(jarg)
            enddo
         else
            write(0,'(a)') 'Not enough arguments supplied for option "'//trim(name)//'".'
            stop 1
         endif
         v1 = varray(1)
         if (p2) v2 = varray(2)
         if (p3) v3 = varray(3)
         if (p4) v4 = varray(4)
         if (p5) v5 = varray(5)
         if (p6) v6 = varray(6)
         if (p7) v7 = varray(7)
         if (p8) v8 = varray(8)
         if (p9) v9 = varray(9)
         if (present(supplied)) supplied = .true.
         return
      endif
      iarg = iarg + 1
   enddo
   
!  Getting this far means we haven't found the requested argument
   if (present(supplied)) supplied = .false.
   
!  Exit program is we haven't found the argument we're after
   if (strict) then
      write(0,'(a)') 'Option "'//trim(name)//'" is required.'
      stop 1
   endif

end subroutine get_arg_int
!-------------------------------------------------------------------------------

!===============================================================================
subroutine get_arg_char(name,v1,v2,v3,v4,v5,v6,v7,v8,v9,required,supplied)
!===============================================================================
!  Arguments:
!     name : name of switch
!     nargs: how many arguments this switch accepts
!     v1    : value to assign
!     v2,v3,v4,v5,v6,v7,v8,v9 : other optional values
!     required: logical.  If true, program exits if argument isn't supplied
!     supplied: logical.  True if we get a value for this option.

   implicit none
   character(len=*) :: name
   character(len=*),intent(inout),optional :: v1,v2,v3,v4,v5,v6,v7,v8,v9
   logical,intent(in),optional :: required
   logical,intent(out),optional :: supplied
   character(max_len),dimension(9) :: varray
   logical :: strict,p1,p2,p3,p4,p5,p6,p7,p8,p9
   integer :: nargs,iarg,jarg,n,length
   character(max_len) :: arg
   
!  Decide if this option must be present
   strict = .false.
   if (present(required)) strict = required

!  Check how many options required
   p1 = present(v1)
   p2 = present(v2); p3 = present(v3); p4 = present(v4); p5 = present(v5)
   p6 = present(v6); p7 = present(v7); p8 = present(v8); p9 = present(v9)
   if (.not.p1.and..not.p2.and..not.p3.and..not.p4.and..not.p5.and..not.p6.and..not.p7.and..not.p8.and..not.p9) then
      nargs = 0
   else if (p1.and..not.p2.and..not.p3.and..not.p4.and..not.p5.and..not.p6.and..not.p7.and..not.p8.and..not.p9) then
      nargs = 1
   else if (p1.and.p2.and..not.p3.and..not.p4.and..not.p5.and..not.p6.and..not.p7.and..not.p8.and..not.p9) then
      nargs = 2
   else if(p1.and.p2.and.p3.and..not.p4.and..not.p5.and..not.p6.and..not.p7.and..not.p8.and..not.p9) then
      nargs = 3
   else if(p1.and.p2.and.p3.and.p4.and..not.p5.and..not.p6.and..not.p7.and..not.p8.and..not.p9) then
      nargs = 4
   else if(p1.and.p2.and.p3.and.p4.and.p5.and..not.p6.and..not.p7.and..not.p8.and..not.p9) then
      nargs = 5
   else if(p1.and.p2.and.p3.and.p4.and.p5.and.p6.and..not.p7.and..not.p8.and..not.p9) then
      nargs = 6
   else if(p1.and.p2.and.p3.and.p4.and.p5.and.p6.and.p7.and..not.p8.and..not.p9) then
      nargs = 7
   else if(p1.and.p2.and.p3.and.p4.and.p5.and.p6.and.p7.and.p8.and..not.p9) then
      nargs = 8
   else if(p1.and.p2.and.p3.and.p4.and.p5.and.p6.and.p7.and.p8.and.p9) then
      nargs = 9
   else
      write(0,'(a)') 'get_arg_real: supply values to get sequentially'
      stop 1
   endif
   
!  Check we're not supplying character strings longer than our maximum
   if (nargs > 0) then
      if (len(v1) > max_len) &
         write(0,'(a)') 'get_arg_char: Warning: Requesting string longer than max in get_args.'
   endif

!  If no arguments supplied, return
   n = command_argument_count()
   if (n == 0) then
      if (strict) then
         write(0,'(a)') 'Option "'//trim(name)//'" is required.'
         stop 1
      endif
      if (present(supplied)) supplied = .false.
      return
   endif
   
!  Check there are enough arguments
   if (n < nargs + 1) then
      write(0,'(a)') 'Not enough arguments supplied for option "'//trim(name)//'".'
      stop 1
   endif
   
!  Loop over the arguments until we find the one we're after
   iarg = 1
   do while (iarg <= n)
      call get_command_argument(iarg,arg,length)
      if (length > max_len) then
         write(0,'(a)') 'Argument to option "'//trim(name)//'" is too long for get_args module.'
         stop 1
      endif
      if (name == arg) then
         !  Check for enough arguments
         if (n >= iarg + nargs) then
            do jarg=1,nargs
               call get_command_argument(iarg+jarg,varray(jarg))
            enddo
         else
            write(0,'(a)') 'Not enough arguments supplied for option "'//trim(name)//'".'
            stop 1
         endif
!         write(*,*) 'v1 length = ',len(v1)
!         v1 = ''
         if (p1) v1 = varray(1)
         if (p2) v2 = varray(2)
         if (p3) v3 = varray(3)
         if (p4) v4 = varray(4)
         if (p5) v5 = varray(5)
         if (p6) v6 = varray(6)
         if (p7) v7 = varray(7)
         if (p8) v8 = varray(8)
         if (p9) v9 = varray(9)
         if (present(supplied)) supplied = .true.
         return
      endif
      iarg = iarg + 1
   enddo
   
!  Getting this far means we haven't found the requested argument
   if (present(supplied)) supplied = .false.
   
!  Exit program is we haven't found the argument we're after
   if (strict) then
      write(0,'(a)') 'Option "'//trim(name)//'" is required.'
      stop 1
   endif

end subroutine get_arg_char
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
end module get_args