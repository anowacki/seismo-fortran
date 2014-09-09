!===============================================================================
! One of Andy Nowacki's Fortran utility modules for dealing with seismic
! anisotropy and other problems.
!
! Andy Nowacki <andy.nowacki@bristol.ac.uk>
!
! See the file LICENCE for licence details.
!===============================================================================
program MT2rp
!===============================================================================
!  MT2rp outputs the radiation pattern for a moment tensor, including the source
!  polarisation.

   use moment_tensor

   implicit none
   
   integer, parameter :: rs = 8
   real(rs) :: M(6), azi, inc, P, SV, SH, j
   logical :: dirs_on_stdin = .true., MTs_on_stdin = .true.
   integer :: i, iostat
   
   ! Get command line parameters and options
   call get_args
   
   ! Loop over input on stdin
   if (dirs_on_stdin .or. MTs_on_stdin) then
      iostat = 0
      do while (iostat == 0)
         if (dirs_on_stdin .and. .not.MTs_on_stdin) then
            read(*,*,iostat=iostat) azi, inc
         else if (.not.dirs_on_stdin .and. MTs_on_stdin) then
            read(*,*,iostat=iostat) (M(i), i=1,6)
         else
            read(*,*,iostat=iostat) azi, inc, (M(i), i=1,6)
         endif
         
         ! Check EOF or bad input
         if (iostat < 0) exit
         if (iostat > 0) then
            write(0,'(a)') 'MT2rp: Error: Problem reading values from stdin'
            stop
         endif
         
         call mt_radiation_pattern(M, azi, inc, P, SV, SH, j)
         call write_output
      enddo

   ! If no lines to process on stdin, do that one
   else
      call mt_radiation_pattern(M, azi, inc, P, SV, SH, j)
      call write_output
   endif
   
contains
   !============================================================================
   subroutine usage()
   !============================================================================
      write(0,'(a)') &
         'Usage: MT2rp (options) < [Mrr Mtt Mpp Mrt Mrp Mtp] [azi] [inc]', &
         'Output:', &
         '   azi, inc, P, SV, SH, j', &
         'Directions (degrees):', &
         '   inc : Angle from down towards upwards', &
         '   azi : Azimuth measured from north towards east', &
         '   j   : S wave polarisation looking along ray, measured clokwise from upwards', &
         'Options:', &
         '   -d [azi inc] : Read directions from command line; MTs are on stdin', &
         '                  (can be used with -m)', &
         '   -m [Mrr Mtt Mpp Mrt Mrp Mtp] : Read MT from command line; read dirs', &
         '                  from stdin (can be used with -d)'
      stop      
   end subroutine usage
   !----------------------------------------------------------------------------

   !============================================================================
   subroutine write_output()
   !============================================================================
      write(*,'(2(f7.2,1x),3(e12.6,1x),f7.2)') azi, inc, P, SV, SH, j
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
            case('-d')
               call get_command_argument(iarg+1, arg)
               read(arg,*) azi
               call get_command_argument(iarg+2, arg)
               read(arg,*) inc
               dirs_on_stdin = .false.
               iarg = iarg + 3
            case('-h')
               call usage
            case('-m')
               call get_mt_from_cmd_line(iarg+1)
               MTs_on_stdin = .false.
               iarg = iarg + 7
            case default
               write(0,'(a)') 'MT2rp: Error: Unrecognised option "'//trim(arg)//'"'
               stop
         end select
      enddo
      
   end subroutine get_args
   !----------------------------------------------------------------------------
   
   !===============================================================================
   subroutine get_mt_from_cmd_line(iarg)
   !===============================================================================
      implicit none
      integer, intent(in) :: iarg
      character(len=250) :: arg
      
      do i=iarg, iarg+5
         call get_command_argument(i, arg)
         read(arg,*,iostat=iostat) M(i-iarg+1)
         if (iostat /= 0) then
            write(0,'(a,i0.0,a)') 'MT2rp: Error: Cannot get value ', i-iarg+1, &
               ' of MT from command line argument "'//trim(arg)//'"'
            stop
         endif
      enddo
      
   end subroutine get_mt_from_cmd_line
   !-------------------------------------------------------------------------------
   

end program MT2rp
!-------------------------------------------------------------------------------