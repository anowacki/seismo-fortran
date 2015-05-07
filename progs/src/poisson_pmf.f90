!===============================================================================
program possion_pmf
!===============================================================================
! Return the value of the probability mass function for the Poisson distribution
! with parameter lambda at integer value k.

   use statistical, only: stat_poisson_pmf

   implicit none

   logical :: read_stdin = .false.
   real(8) :: lambda
   integer :: i, k, k2, iline = 1
   integer :: iostat = 0

   call get_args

   if (read_stdin) then
      do while (iostat == 0)
         read(*,*,iostat=iostat) lambda, k
         if (iostat < 0) exit
         if (iostat > 0) then
            write(0,'(a,i0.1,a)') 'poisson_pmf: Error: Problem reading lambda, k ' &
               // 'from line ', iline, ' of stdin'
            error stop
         endif
         write(*,*) k, stat_poisson_pmf(lambda, k)
      enddo
   else
      do i = k, k2
         write(*,*) i, stat_poisson_pmf(lambda, i)
      enddo
   endif

contains
   subroutine usage
      write(0,'(a)') &
         'Usage: possion_pmf [lambda] [k] (k2)', &
         '       OR', &
         '       poisson_pmf < (lambda,k on stdin)', &
         'Return the value of the pmf for the Poisson distribution with parameter', &
         'lambda at integer value k, and optionally between k and k2 inclusive.'
     error stop
   end subroutine usage   

   subroutine get_args
     integer :: iarg, narg
     character(len=250) :: arg
     narg = command_argument_count()
     if (narg == 1 .or. narg > 3) call usage
     if (narg == 0) then
        read_stdin = .true.
     else
        call get_command_argument(1, arg)
        read(arg,*) lambda
        call get_command_argument(2, arg)
        read(arg,*) k
        if (narg == 3) then
           call get_command_argument(3, arg)
           read(arg,*) k2
        else
           k2 = k
        endif
     endif
   end subroutine get_args

end program possion_pmf
!-------------------------------------------------------------------------------