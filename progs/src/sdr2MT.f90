!===============================================================================
! One of Andy Nowacki's Fortran utility modules for dealing with seismic
! anisotropy and other problems.
!
! Andy Nowacki <andy.nowacki@bristol.ac.uk>
!
! See the file LICENCE for licence details.
!===============================================================================
program sdr2MT_prog
!===============================================================================
! sdr2MT creates a moment tensor for a double couple earthquake, given a strike,
! dip and rake.
!
! Will output either a list of six moments, or GMT- or SPECFEM3D_*-compatible
! strings.

   use moment_tensor

   implicit none

   integer, parameter :: rs = 8
   real(rs) :: M(6)
   real(rs) :: strike, dip, rake, M0, Mw
   integer :: i, iostat
   logical :: Mw_set = .false., M0_set = .false., sdr_on_stdin = .true.
   character(len=1) :: out_format = 'n'

   call get_args

   if (.not.sdr_on_stdin) then
      if (Mw_set) then
         M = mt_sdr2mt(strike, dip, rake, Mw=Mw)
      else if (M0_set) then
         M = mt_sdr2mt(strike, dip, rake, M0=M0)
      else
         M = mt_sdr2mt(strike, dip, rake)
      endif
      call write_output
   else
      iostat = 0
      do while (iostat == 0)
         read(*,*,iostat=iostat) strike, dip, rake
         if (iostat < 0) exit
         if (iostat > 0) then
            write(0,'(a)') 'sdr2MT: Error: Problem read strike, dip, rake from stdin'
            stop
         endif
         if (Mw_set) then
            M = mt_sdr2mt(strike, dip, rake, Mw=Mw)
         else if (M0_set) then
            M = mt_sdr2mt(strike, dip, rake, M0=M0)
         else
            M = mt_sdr2mt(strike, dip, rake)
         endif
         call write_output
      enddo
   endif

contains
   !============================================================================
   subroutine usage()
   !============================================================================
      write(0,'(a)') &
         'Usage: sdr2MT (options) < [strike] [dip] [rake]', &
         'Output:', &
         '   azi, inc, P, SV, SH, j', &
         'Directions (degrees):', &
         '   strike : Azimuth from north towards east of strike of fault plane', &
         '   dip    : Dip of fault plane from horizontal', &
         '   rake   : Fault motion measured anticlockwise from strike on footwall', &
         'Options:', &
         '   -o [G|S]: Output tensor for (G) plotting with psmeca -Sm option, or', &
         '             (S) use with SPECFEM3D_GLOBE CMTSOLUTION [list of 6 numbers]', &
         '   -sdr [strike dip rake]: Supply orientation on command line [read from stdin]', &
         '   -Mw|-M0 : Specify moment magnitude or moment [Mw = 5]', &
         '             (last -Mw or -M0 switch has effect and previous ones are ignored)', &
         '   -h      : Display this message'
      stop
   end subroutine usage
   !----------------------------------------------------------------------------

   !============================================================================
   subroutine write_output()
   !============================================================================
      select case(out_format)
         case('n')
            write(*,'(6e14.6)') (M(i), i=1,6)
         case('g', 'G')
            write(*,'(a,6e14.6,a)') '0 0 0 ',(M(i), i=1,6),' 23 0 0'
         case('s', 'S')
            write(*,'(a,e13.6)') &
              'Mrr:      ',M(1), &
              'Mtt:      ',M(2), &
              'Mpp:      ',M(3), &
              'Mrt:      ',M(4), &
              'Mrp:      ',M(5), &
              'Mtp:      ',M(6)
        case default
           write(0,'(a)') 'sdr2MT: Error: Output format "'//out_format//'" not recognised'
           stop
      end select
   end subroutine write_output
   !----------------------------------------------------------------------------

   !============================================================================
   subroutine get_args()
   !============================================================================
      implicit none
      integer :: iarg, narg
      character(len=250) :: arg

      narg = command_argument_count()
      iarg = 1
      do while (iarg <= narg)
         call get_command_argument(iarg, arg)
         select case(arg)
            case('-h')
               call usage
            case('-o')
               call get_command_argument(iarg+1, out_format)
               iarg = iarg + 2
            case('-sdr')
               call get_command_argument(iarg+1, arg)
               read(arg,*) strike
               call get_command_argument(iarg+2, arg)
               read(arg,*) dip
               call get_command_argument(iarg+3, arg)
               read(arg,*) rake
               sdr_on_stdin = .false.
               iarg = iarg + 4
            case('-Mw', '-MW', '-mw')
               call get_command_argument(iarg+1, arg)
               read(arg,*) Mw
               Mw_set = .true.
               M0_set = .false.
               iarg = iarg + 2
            case('-M0', '-m0')
               call get_command_argument(iarg+1, arg)
               read(arg,*) M0
               M0_set = .true.
               Mw_set = .false.
               iarg = iarg + 2
            case default
               write(0,'(a)') 'MT2rp: Error: Unrecognised option "'//trim(arg)//'"'
               stop
         end select
      enddo

   end subroutine get_args
   !----------------------------------------------------------------------------

end program sdr2MT_prog
!-------------------------------------------------------------------------------