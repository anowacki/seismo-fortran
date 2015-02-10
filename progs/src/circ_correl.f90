!===============================================================================
! One of Andy Nowacki's Fortran utility modules for dealing with seismic
! anisotropy and other problems.
!
! Andy Nowacki <andy.nowacki@bristol.ac.uk>
!
! See the file LICENCE for licence details.
!===============================================================================
program circ_correl_prog
!===============================================================================
!  Return the correlation of two sets of angles read from stdin

   use statistical, only: circ_correl
   
   implicit none
   
   integer,parameter :: nmax = 10000
   real(8) :: atemp,btemp, a(nmax),b(nmax)
   integer :: n, iostatus
   
   if (command_argument_count() /= 0) then
      write(0,'(a)') 'Usage: circ_stats < [list of angles (degrees)]'
      stop
   endif
   
!  Read values from stdin
   iostatus = 0
   n = 0
   do while (iostatus == 0)      
      read(*,*,iostat=iostatus) atemp,btemp
      if (iostatus > 0) then    ! Bad input
         write(0,'(a)') 'circ_stats: Error: problem reading from stdin.'
         stop
      endif
      if (iostatus < 0) exit    ! EOF

      if (n+1 > nmax) then      ! Read past compiled limits of a(nmax)
         write(0,'(a,i0.1,a)') 'circ_stats: Error: Compiled max no. data is ', &
           nmax, '.  Current datum is beyond that.  Recompile with larger nmax.'
         stop
      endif

      a(n+1) = atemp
      b(n+1) = btemp
      n = n + 1
   enddo

!  Compute and output correlation
   write(*,'(f0.6)') circ_correl(a(1:n), b(1:n), degrees=.true.)

end program circ_correl_prog
