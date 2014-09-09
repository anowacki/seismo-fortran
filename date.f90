!===============================================================================
! One of Andy Nowacki's Fortran utility modules for dealing with seismic
! anisotropy and other problems.
!
! Andy Nowacki <andy.nowacki@bristol.ac.uk>
!
! See the file LICENCE for licence details.
!===============================================================================
module date
!===============================================================================
!  Fortran module containing utility functions for dealing with dates and time
!  
!  Andy Nowacki, University of Bristol
!  2011 - 
!  
!-------------------------------------------------------------------------------
!  Revisions:
!      2011-09-30:  Incept date.  Contains calday and julday only
!-------------------------------------------------------------------------------

!  Declare parameters
   integer,parameter,private :: days_in_month_nonleap(12) = &
      (/31,28,31,30,31,30,31,31,30,31,30,31/)
   integer,parameter,private :: days_in_month_leap(12) = 
      (/31,29,31,30,31,30,31,31,30,31,30,31/)
   character(len=3),parameter,private :: month_name_short(12) = &
      (/'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'/)
   character(len=9),parameter,private :: month_name_long(12) = &
      (/'January','February','March','April','May','June', &
        'July','August','September','October','November','December'/)
   character(len=3),parameter,private :: weekday_name_short(7) = &
      (/'Mon','Tue','Wed','Thu','Fri','Sat','Sun'/)
   character(len=9),parameter,private :: weekday_name_long(7) = &
      (/'Monday','Tuesday','Wednesday','Thursday','Friday','Saturday','Sunday'/)

!  IO units
   integer, parameter, private :: lu_stderr = 0, &
                                  lu_stdin  = 5, &
                                  lu_stdout = 6
!  Declare derived types
   type date
      integer :: yr,mo,d,jd,hr,m,is
      real(4) :: s
      logical :: isleap
   end type date
   
CONTAINS

!===============================================================================
   function date_isleap(year)
!===============================================================================
!  Returns .true. if year is a leap year, false if not
      implicit none
      integer, intent(in) :: year
      logical :: date_isleap
!  Leap years are divisible by 4 or 400, but not 100 (so 2000 was leap, 2100 isn't)
      if ((mod(year,4) == 0 .and. mod(year,100) /= 0) .or. mod(year,400) == 0) then
         date_isleap = .true.
      else
         date_isleap = .false.
      endif
      return
   end function date_isleap
!-------------------------------------------------------------------------------

!===============================================================================
   function date_exists(year,month_or_julday,day)
!===============================================================================
!  Returns true if date exists.
      implicit none
      integer, intent(in) :: year,month
      integer, intent(in), optional :: day
      logical :: date_exists
      
!  If day supplied, then assume this is in format year-month-day
      if (present(day)) then
         if (month < 1 .or. month > 12) then
            write(lu_stderr,'(a,i0,a)') &
               'date: date_exists: ',month,' is not a valid month.'
            stop
         elseif ((date_isleap(year) .and. day /= days_in_month_leap(month)) .or. &
                 (.not.date_isleap(year) .and. day /= days_in_month_nonleap(month))) then
            date_exists = .false.
            return
         else
            date_exists = .true.
            return
!  If no day supplied, assume this is in format year-julday.  This then only makes
!  sure that we're not claiming a false leap year
      else
         if (month < 1 .or. (month > 366 .and. date_isleap(year)) .or. &
                            (month > 365 .and. .not.date_isleap(year))) then
            date_exists = .false.

!===============================================================================
   function date_julday(year,month,day)
!===============================================================================
!  Returns the number of the day in the year
      implicit none
      integer,intent(in) :: year,month,day
      integer :: date_julday
      integer :: i
!  Check this date exists
      if (.not.date_exists(year,month,day=day) then
         write(lu_stderr,'(i0,x,i2.2,x,i2.2,a)') year,month,day,' is not a valid date.'
         stop
      endif
      
      date_julday = 0
      if (.not.isleap(year)) then
         do i=1,month-1
            date_julday = date_julday + days_in_month_nonleap(i)
         enddo
         date_julday = date_julday + day
         return
      else
         do i=1,month-1
            date_julday = date_julday + days_in_month_leap(i)
         enddo
         date_julday = date_julday + day
      endif
      return
   end function date_julday
!-------------------------------------------------------------------------------

!===============================================================================
   function date_day_of_week(year,month,day,long)
!===============================================================================
!  Returns the day of the week.
!  Uses the method of Tomohiko Sakamoto as tken from (3.) at
!  http://c-faq.com/misc/zeller.html
!  This version 
      implicit none
      integer,intent(in) :: year,month,day
      character(len=*) :: date_day_of_week
      logical,optional,intent(in) :: long
      integer :: t(12) = (/0, 3, 2, 5, 0, 3, 5, 1, 4, 6, 2, 4/)
      integer :: iday_of_week
      
      if (month < 3) year = year - 1
      iday_of_week = mod(year + year/4 - year/100 + t(month) + day, 7)
      if (iday_of_week == 0) iday_of_week = 7  ! Sunday is day 7, not day 0
      
      if (present(long)) then
         if (long) then
            date_day_of_week = weekday_name_long(iday_of_week)
            return
         endif
      endif
      
      date_day_of_week = weekday_name_short(iday_of_week)
      return
   end function date_day_of_week
!-------------------------------------------------------------------------------
