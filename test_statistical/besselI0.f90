program besselI0

   implicit none
   real(8) :: x,dx
   
   x = 0.
   dx = 0.1
   
   do while (x <= 4.)
      write(*,*) x,BessI0(x)
      x = x + dx
   enddo

   contains
      function BessI0(x)
      !  Evaluates the modified Bessel function of the first kind at x.
      !  Polynomial approximation taken from Abramowitz & Stegun, 1964, Handbook
      !   of mathematical functions
         implicit none
         integer,parameter :: r8=8
         integer,parameter :: rs=r8
         real(rs),intent(in) :: x
         real(rs)            :: BessI0
         real(r8)            :: t
         
         if (abs(x) < 3.75_r8) then
            t = abs(x/3.75_r8)
            BessI0 = 1._r8 + &
                     3.5156229_r8*t**2 + 3.0899424_r8*t**4 +  1.2067492_r8*t**6 + &
                     0.2659732_r8*t**8 + 0.0360768_r8*t**10 + 0.0045813_r8*t**12
         else
            t = abs(3.75_r8/x)
            BessI0 = (exp(abs(x))/sqrt(abs(x)))*(0.39894228_rs + &
                     0.01328592_r8*t    + 0.00225319_r8*t**2 - 0.00157565_r8*t**3 + &
                     0.00916281_r8*t**4 - 0.02057706_r8*t**5 + 0.02635537_r8*t**6 - &
                     0.01647633_r8*t**7 + 0.00392377_r8*t**8)
         endif
         return
      end function BessI0

end program