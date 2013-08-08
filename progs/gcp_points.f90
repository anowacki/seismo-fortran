!===============================================================================
program gcp_points_prog
!===============================================================================
!  Wrapper program to gcp_points subroutine.  See 
   use spherical_geometry

   implicit none
   
   integer, parameter :: rs = 8
   real(rs) :: lon1,lat1,lon2,lat2,d
   real(rs),allocatable,dimension(:) :: lon,lat
   real(rs) :: ds
   integer,parameter :: np_default = 50  ! Default to 50 points along path
   integer :: i,iostat,npout,np
   logical :: d_supplied = .false., n_supplied = .false., p_supplied = .false.
   ! Path separator, hard-wired for now
   character(len=1), parameter :: sep = ">"
   
   call get_arguments
   
   ! Have asked for constant number of points, or are using the default np
   if (n_supplied .or. (.not.n_supplied .and. .not.d_supplied)) then
      if (.not.n_supplied) np = np_default
      allocate(lon(np),lat(np))
   endif
   
   if (.not.p_supplied) then
      iostat = 0
      do while (iostat == 0)
         call read_points
         if (iostat < 0) exit ! End of file
         ! If we're specifying constant spacing
         if (d_supplied) then
            np = ceiling(delta(lon1,lat1,lon2,lat2,degrees=.true.)/ds) + 1
            allocate(lon(np),lat(np))
            call gcp_points(lon1,lat1,lon2,lat2,lon,lat,ds=ds,npts=npout,degrees=.true.)
            np = npout
            call write_points
            deallocate(lon,lat)
         ! Otherwise we're asking for a constant number of points, which are already
         ! allocated
         else
            call gcp_points(lon1,lat1,lon2,lat2,lon,lat,n=np,degrees=.true.)
            call write_points
         endif
      enddo
   
   ! Using one-shot command line values
   else
      ! If we haven't allocated, then we need to do so for a constant spacing
      if (d_supplied) then
         np = ceiling(delta(lon1,lat1,lon2,lat2,degrees=.true.)/ds) + 1
         allocate(lon(np),lat(np))
         call gcp_points(lon1,lat1,lon2,lat2,lon,lat,ds=ds,npts=npout,degrees=.true.)
         np = npout
         call write_points
         deallocate(lon,lat)
      else
         call gcp_points(lon1,lat1,lon2,lat2,lon,lat,n=np,degrees=.true.)
         call write_points
      endif
   endif
   
   if (allocated(lon)) deallocate(lon,lat)
   
contains
   !============================================================================
   subroutine read_points
   !============================================================================
      implicit none
      read(*,*,iostat=iostat) lon1,lat1,lon2,lat2
      if (iostat > 0) then
         write(0,'(a)') 'gcp_points: Error: problem reading lon1,lat1,lon2,lat2 from stdin'
         stop
      endif
   end subroutine read_points
   !----------------------------------------------------------------------------
   
   !============================================================================
   subroutine write_points
   !============================================================================
      implicit none
      do i=1,np
         write(*,*) lon(i),lat(i)
      enddo
   end subroutine write_points
   !----------------------------------------------------------------------------
   
   !============================================================================
   subroutine get_arguments()
   !============================================================================
      use get_args
      implicit none
      character(len=250) :: arg
      if (command_argument_count() == 0) then
         return
      else
         call get_arg('-d',ds,supplied=d_supplied)
         call get_arg('-n',np,supplied=n_supplied)
         if (d_supplied .and. n_supplied) then
            write(0,'(a)') 'gcp_points: Error: Must supply only one of -n or -d'
            call usage
         endif
         ! Manual search for -p
         i = 1
         arg_loop: do while (i <= command_argument_count())
            call get_command_argument(i,arg)
            if (arg == "-p") then
               if (command_argument_count() < i + 4) then
                  write(0,'(a)') 'gcp_points: Error: not enough parameters supplied for option "-p"'
                  call usage
               endif
               call get_command_argument(i+1,arg); read(arg,*) lon1
               call get_command_argument(i+2,arg); read(arg,*) lat1
               call get_command_argument(i+3,arg); read(arg,*) lon2
               call get_command_argument(i+4,arg); read(arg,*) lat2
               p_supplied = .true.
               exit arg_loop
            endif
            i = i + 1
         enddo arg_loop
      endif
   end subroutine get_arguments
   !----------------------------------------------------------------------------

   !============================================================================
   subroutine usage()
   !============================================================================
   !  Usage subroutine to tell user how to invoke
      write(0,'(a)') 'Usage: gcp_points (options)', &
                     '   -d [ds]                  : Set spacing of points',&
                     '   -n [npts]                : Set number of points',&
                     '   -p [lon1 lat1 lon2 lat2] : Set endpoints on command line',&
                     '',&
                     'By default, gcp_points reads sets of endpoints from stdin as lon1,lat1,lon2,lat2',&
                     'and writes intermediate points to along the great circle path to stdout,',&
                     ' separated by ">"s for use with GMT.'
      stop
   end subroutine usage
   !----------------------------------------------------------------------------
end program gcp_points_prog
!-------------------------------------------------------------------------------
