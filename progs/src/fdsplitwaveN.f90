!===============================================================================
! One of Andy Nowacki's Fortran utility modules for dealing with seismic
! anisotropy and other problems.
!
! Andy Nowacki <andy.nowacki@bristol.ac.uk>
!
! See the file LICENCE for licence details.
!===============================================================================
program fdsplitwaveN
!===============================================================================
! Applies N splitting operators to a set of SAC files in the frequency domain.
! Reads the operators from stdin, or if only one set is desired, from the command
! line.
!
! Usage:
!    echo "phi1 dt1 \n phi2 dt2 \n ... phiN dtN" |
!       fdsplitwaveN [filename]
!    OR:
!       fdsplitwaveN [filename] [phi] [dt]

   use f90sac
   use splitwave

   implicit none
   
   type(SACtrace) :: t1, t2, t3
   integer, parameter :: nsplits_max = 10000
   real, dimension(nsplits_max) :: theta, dt
   real :: theta_temp, dt_temp
   character(len=250) :: file, file1, file2, file3
   integer :: i, iostatus, nsplits

   ! Check arguments
   if (command_argument_count() /= 0 .and. command_argument_count() /= 1) then
      write(0,'(a)') 'Usage:  fdsplitwaveN [filename] < [splitting parameters file]'
      write(0,'(a)') 'Acts on files ''filename.BH[E,N,Z]''.'
      write(0,'(a)') 'Splitting parameters file contains [phi] [dt] on each line.'
      stop
   endif

 ! Get filename base
   file = 'wave'
   if (command_argument_count() == 1) then
      call get_command_argument(1,file)
   endif
  
   file1 = trim(file) // '.BHN'
   file2 = trim(file) // '.BHE'
   file3 = trim(file) // '.BHZ'
  
   ! Load traces
   call f90sac_readtrace(file1,t1)
   call f90sac_readtrace(file2,t2)
   call f90sac_readtrace(file3,t3)
  
   ! Check traces are same length
   if (t1%npts /= t2%npts) then
      write(0,'(a)') 'fdsplitwaveN: Traces must be same length!'
      stop
   endif
   
   ! Read operators from stdin
   iostatus = 0
   i = 0
   do while(iostatus == 0)
      read(*,*,iostat=iostatus) theta_temp, dt_temp
      if (iostatus < 0) exit
      if (iostatus > 0) then
         write(0,'(a)') 'fsplitwaveN: Some problem reading splitting parameters'
         stop
      endif
      i = i + 1
      theta(i) = -theta_temp
      dt(i) = dt_temp
   enddo
   
   nsplits = i
   
   ! Split traces
   call sw_fdsplitN(t1,t2,t3, nsplits, theta(1:nsplits), dt(1:nsplits), quiet=.true.)
   
   ! Write traces out
   call f90sac_writetrace(file1,t1)
   call f90sac_writetrace(file2,t2)
   call f90sac_writetrace(file3,t3)

end program fdsplitwaveN
!-------------------------------------------------------------------------------
