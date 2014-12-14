!===============================================================================
program circ_stats
!===============================================================================
!  Returns some circular statistics

   use statistical, only: circ_mean, circ_res_length, circ_variance
   
   implicit none
   
   integer, parameter :: nmax = 100000
   real(8) :: atemp, a(nmax)
   integer :: n, iostatus
   
   if (command_argument_count() /= 0) then
      write(0,'(a)') 'Usage: circ_stats < [list of angles (degrees)]'
      stop
   endif
   
!  Read values from stdin
   iostatus = 0
   n = 0
   do while (iostatus == 0)      
      read(*,*,iostat=iostatus) atemp
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
      n = n + 1
   enddo
   
!  Calculate statistics and write them out
   write(*,'("n = ",i0,"  Mean = ",f0.5,"  R = ",f0.5,"  Var = ",f0.6)') &
         n, circ_mean(a(1:n)), circ_res_length(a(1:n)), circ_variance(a(1:n))
   
end program