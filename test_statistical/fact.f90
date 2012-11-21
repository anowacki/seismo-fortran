program factorial
!  Test the floating-point factorial approximation
   use statistical
   implicit none
   integer :: n,nfact12
   real(8) :: nfact
   character(len=80) :: arg
   
   if (iargc() /= 1) then
      write(*,'(a)') 'Usage: fact [n]'
      stop
   endif
   call getarg(1,arg) ; read(arg,*) n
   write(*,*) n,'factorial:'
   nfact = fact(n)
   write(*,*) nfact
   if (n <= 12) then
      write(*,*) 'True integer value:'
      nfact12 = fact12(n)
      write(*,*) nfact12
      write(*,*) 'Relative error of approximation =', &
                           abs(real(nfact12)-nfact)/real(nfact12)
   endif
   
end program factorial
