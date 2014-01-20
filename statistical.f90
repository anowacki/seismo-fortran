!===============================================================================
module statistical
!===============================================================================
!  Contains statistical functions
!  Andy Nowacki, University of Bristol
!
!  2011-06-13   + First iteration.  Includes F-distribution, Beta- (B), regularized
!                 Beta- (I) and incomplete Beta-functions
!  2011-06-17   + Added functions from Alan Miller's website to compute percentage
!                 values for the chi2 distribution.
!  2011-08-09   + Added circular statistical functions.

!  Declare public functions
!   public f_dist
   
!  Declare some internal constants
!  ** size constants
      integer, parameter, private :: i4 = selected_int_kind(9) ; ! long int
      integer, parameter, private :: r4 = selected_real_kind(6,37) ; ! SP
      integer, parameter, private :: r8 = selected_real_kind(15,307) ; ! DP
      integer, parameter, private :: c4 = r4
      integer, parameter, private :: c8 = r8
      
!  ** precision selector
      integer, parameter, private :: rs = r8
      integer, parameter, private :: cs = c8
      
!  ** maths constants and other useful things
      real(rs), parameter, private :: pi = 3.141592653589793238462643_rs
      
!  ** IO units
      integer, parameter, private :: lu_stdout = 5
      integer, parameter, private :: lu_stdin  = 6
      integer, parameter, private :: lu_stderr = 0

contains

!===============================================================================
function fact12(n)
!===============================================================================
!  Gives true integer value for n!, n <= 12
   implicit none
   integer,intent(in) :: n
   integer            :: fact12,i
   
   if (n < 0) then
      write(lu_stderr,'(a)') 'statistical: fact: error: n must be >= 0.'
      stop
   else if (n == 0) then
      fact12 = 1
      return
   else if (n > 12) then
      write(lu_stderr,'(a)') 'statistical: fact: error: n cannot be larger than 12 for 4-byte integer calculations.'
      stop
   endif
   fact12 = 1
   do i=1,n
      fact12 = i*fact12
   enddo
   return
end function fact12
!-------------------------------------------------------------------------------

!===============================================================================
function fact(n)
!===============================================================================
!  Gives approximate value of n! up to n<170
!  fact(n) ≈ exp{n.ln(n) - n + ln(n(1 + 4n(1 + 2n)))/6 + ln(pi)/2}
!  Approximation given by Ramanujan (1988)
!  Good to within ~1e-5 of integer value up to n=12.
   implicit none
   integer,intent(in) :: n
   real(rs)           :: fact,rn
   
   if (n == 0) then
      fact = 1._rs
      return
   else if (n < 0) then
      write(lu_stderr,'(a)') 'statistical: fact: error: n must be >= 0.'
      stop
   else if (n > 0 .and. n <= 12) then
      fact = real(fact12(n))
      return
   else if (n > 170) then
      write(lu_stderr,'(a)') 'statistical: fact: error: cannot represent n! in double precision for n > 170.'
      stop
   endif
   
   rn = real(n)
   fact = exp(n*log(rn) - rn + log(rn*(1._rs + 4._rs*rn*(1._rs + 2._rs*rn)))/6._rs &
                                                         + log(pi)/2._rs)
   return
end function fact
!-------------------------------------------------------------------------------

!===============================================================================
function beta_func(p,q)
!===============================================================================
!  Returns the Beta function B(p,q) = Gamma(p)*Gamma(q) / Gamma(p+q)
   implicit none
   real(rs),intent(in) :: p,q
   real(rs)            :: beta_func,rp,rq
   
   rp = real(p)  ;  rq = real(q)
   beta_func = gamma(rp) * gamma(rq) / gamma(rp + rq)
   return
end function beta_func
!-------------------------------------------------------------------------------

!===============================================================================
function incomp_beta_func(z,a,b)
!===============================================================================
!  Returns the incomplete Beta function B(z;a,b) = ∫_0^z(u^(a-1)*(1-u)^(b-1))du
   implicit none
   real(rs),intent(in) :: z,a,b
   real(rs)            :: incomp_beta_func,ra,rb,rn
   real(rs)            :: u,du   ! Dummy variable of integration
   integer             :: n
   
   ra = real(a)
   rb = real(b)
   if (z <= 0. .or. z > 1.) then
      write(lu_stderr,'(a)') 'statistical: incomp_beta_func: error: z must be in range [0,1].'
      stop
   endif
   du = 0.000001_rs
   u = 0._rs
   incomp_beta_func = 0._rs
   do while (u <= z)
      incomp_beta_func = incomp_beta_func + (u**(a-1._rs))*((1._rs-u)**(b-1._rs))*du
      u = u + du
   enddo
   incomp_beta_func = incomp_beta_func / beta_func(a,b)
   return
end function incomp_beta_func
!-------------------------------------------------------------------------------

!===============================================================================
function reg_beta_func(z,a,b)
!===============================================================================
!  Returns the regularised Beta function B(a,b) = B(z;a,b)/B(a,b)
   implicit none
   real(rs),intent(in) :: z,a,b
   real(rs)            :: reg_beta_func
   
   reg_beta_func = incomp_beta_func(z,a,b) / beta_func(a,b)
   return
end function reg_beta_func
!-------------------------------------------------------------------------------

!===============================================================================
function pochhammer(x,n)
!===============================================================================
!  Give the Pochhhammer symbol p(x)_n == gamma(x+n)/gamma(x)
   implicit none
   real(rs),intent(in) :: x
   integer,intent(in)  :: n
   real(rs)            :: pochhammer,rn
   
   rn = real(n)
   pochhammer = gamma(x+rn)/gamma(x)
   return
end function pochhammer
!-------------------------------------------------------------------------------

!===============================================================================
function f_dist(n,m,x)
!===============================================================================
!  Returns the value of the f distribution f(n,m,x)
   implicit none
   integer             :: m,n
   real(rs)            :: x
   real(rs)            :: f_dist
   real(rs)            :: rn,rm,half,one
   
   if (x <= 0.) then
      write(lu_stderr,'(a)') 'statistical: f_dist: error: x must be > 0.'
      stop
   endif
   
!  Convert input integers into real values and set value of half
   rn = real(n)
   rm = real(m)
   half = 0.5_rs
   one = 1.0_rs
   
   f_dist = gamma(half*(rn+rm)) * rn**(half*rn) * rm**(half*rm) * x**(half*rn - one) / &
            (gamma(half*rn) * gamma(half*rm) * (rm+rn*x)**(half*(rn + rm)))
   return
end function f_dist
!-------------------------------------------------------------------------------

!===============================================================================
function f_dist_cum(n,m,x)
!===============================================================================
!  Returns the value of the cumulative F distribution F(n,m,x)
   implicit none
   integer,intent(in)  :: m,n
   real(rs),intent(in) :: x
   real(rs)            :: f_dist_cum,u,du
   real(rs)            :: rm,rn
!  Convert input integers into real values
   rm = real(m)
   rn = real(n)
   
!   f_dist_cum = reg_beta_func(rn*x/(rm+rn*x), rn*0.5_rs, rm*0.5_rs)
   f_dist_cum = 0._rs
   du = 0.0001_rs
   u = 0._rs + du
   do while (u <= x)
      f_dist_cum = f_dist_cum + f_dist(n,m,u)*du
      u = u + du
   enddo
   return
end function f_dist_cum
!-------------------------------------------------------------------------------

!===============================================================================
!===============================================================================
! The following functions are taken from Alan Miller's page at:
! http://jblevins.org/mirror/amiller/
! All are transcriptions into Fortran of algorithms plublished in the Royal Statistical
! Society's Applied Statitstics journal.

!===============================================================================
function ppchi2(p, v, g) RESULT(fn_val)
!===============================================================================
! N.B. Argument IFAULT has been removed.
! This version by Alan Miller
! amiller @ bigpond.net.au
! Latest revision - 27 October 2000
!  Algorithm AS 91   Appl. Statist. (1975) Vol.24, P.35
!
!  To evaluate the percentage points of the chi-squared
!  probability distribution function.
!
!  p must lie in the range 0.000002 to 0.999998,
!  v must be positive,
!  g must be supplied and should be equal to ln(gamma(v/2.0))
!
!  Incorporates the suggested changes in AS R85 (vol.40(1), pp.233-5, 1991)
!  which should eliminate the need for the limited range for p above,
!  though these limits have not been removed from the routine.
!
!  If IFAULT = 4 is returned, the result is probably as accurate as
!  the machine will allow.
!
!  Auxiliary routines required: PPND = AS 111 (or AS 241) and GAMMAD = AS 239.

   IMPLICIT NONE
   INTEGER, PARAMETER    :: dp = rs
   
   REAL (dp), INTENT(IN)  :: p
   REAL (dp), INTENT(IN)  :: v
   REAL (dp), INTENT(IN)  :: g
   REAL (dp)              :: fn_val

! AJN: don't need interface as we're in a module.   
!   INTERFACE
!     FUNCTION gammad(x, p) RESULT(fn_val)
!       IMPLICIT NONE
!       INTEGER, PARAMETER    :: dp = rs
!       REAL (dp), INTENT(IN) :: x, p
!       REAL (dp)             :: fn_val
!     END FUNCTION gammad
!   
!     SUBROUTINE ppnd16 (p, normal_dev, ifault)
!       IMPLICIT NONE
!       INTEGER, PARAMETER      :: dp = rs
!       REAL (dp), INTENT(IN)   :: p
!       INTEGER, INTENT(OUT)    :: ifault
!       REAL (dp), INTENT(OUT)  :: normal_dev
!     END SUBROUTINE ppnd16
!   END INTERFACE
   
! Local variables
   
   REAL (dp)  :: a, b, c, p1, p2, q, s1, s2, s3, s4, s5, s6, t, x, xx
   INTEGER    :: i, if1
   
   INTEGER, PARAMETER    :: maxit = 20
   REAL (dp), PARAMETER  :: aa = 0.6931471806_dp, e = 0.5e-06_dp,         &
                            pmin = 0.000002_dp, pmax = 0.999998_dp,       &
                            zero = 0.0_dp, half = 0.5_dp, one = 1.0_dp,   &
                            two = 2.0_dp, three = 3.0_dp, six = 6.0_dp,   &
                            c1 = 0.01_dp, c2 = 0.222222_dp, c3 = 0.32_dp, &
                            c4 = 0.4_dp, c5 = 1.24_dp, c6 = 2.2_dp,       &
                            c7 = 4.67_dp, c8 = 6.66_dp, c9 = 6.73_dp,     &
                            c10 = 13.32_dp, c11 = 60.0_dp, c12 = 70.0_dp, &
                            c13 = 84.0_dp, c14 = 105.0_dp, c15 = 120.0_dp, &
                            c16 = 127.0_dp, c17 = 140.0_dp, c18 = 175.0_dp, &
                            c19 = 210.0_dp, c20 = 252.0_dp, c21 = 264.0_dp, &
                            c22 = 294.0_dp, c23 = 346.0_dp, c24 = 420.0_dp, &
                            c25 = 462.0_dp, c26 = 606.0_dp, c27 = 672.0_dp, &
                            c28 = 707.0_dp, c29 = 735.0_dp, c30 = 889.0_dp, &
                            c31 = 932.0_dp, c32 = 966.0_dp, c33 = 1141.0_dp, &
                            c34 = 1182.0_dp, c35 = 1278.0_dp, c36 = 1740.0_dp, &
                            c37 = 2520.0_dp, c38 = 5040.0_dp
   
!  Test arguments and initialise
   fn_val = -one
   IF (p < pmin .OR. p > pmax) THEN
      WRITE(lu_stderr,'(a)') 'statistical: PPCHI2: error: p must be between 0.000002 & 0.999998'
      stop !RETURN
   END IF
   IF (v <= zero) THEN
      WRITE(lu_stderr,'(a)') 'statistical: PPCHI2: error: Number of deg. of freedom <= 0'
      stop !RETURN
   END IF
   
   xx = half * v
   c = xx - one
   
!  Starting approximation for small chi-squared
   IF (v < -c5 * LOG(p)) THEN
      fn_val = (p * xx * EXP(g + xx * aa)) ** (one/xx)
      IF (fn_val < e) GO TO 6
      GO TO 4
   END IF
   
!  Starting approximation for v less than or equal to 0.32
   IF (v > c3) GO TO 3
   fn_val = c4
   a = LOG(one-p)
   
   2 q = fn_val
   p1 = one + fn_val * (c7+fn_val)
   p2 = fn_val * (c9 + fn_val * (c8 + fn_val))
   t = -half + (c7 + two * fn_val) / p1 - (c9 + fn_val * (c10 + three * fn_val)) / p2
   fn_val = fn_val - (one - EXP(a + g + half * fn_val + c * aa) * p2 / p1) / t
   IF (ABS(q / fn_val - one) > c1) GO TO 2
   GO TO 4
   
!  Call to algorithm AS 241 - note that p has been tested above.
   3 CALL ppnd16(p, x, if1)
   
!  Starting approximation using Wilson and Hilferty estimate
   p1 = c2 / v
   fn_val = v * (x * SQRT(p1) + one - p1) ** 3
   
!  Starting approximation for p tending to 1
   IF (fn_val > c6 * v + six) fn_val = -two * (LOG(one-p) - c * LOG(half * fn_val) + g)
   
!  Call to algorithm AS 239 and calculation of seven term Taylor series
   4 DO i = 1, maxit
      q = fn_val
      p1 = half * fn_val
      p2 = p - gammad(p1, xx)
    
      t = p2 * EXP(xx * aa + g + p1 - c * LOG(fn_val))
      b = t / fn_val
      a = half * t - b * c
      s1 = (c19 + a * (c17 + a * (c14 + a * (c13 + a * (c12 + c11 * a))))) / c24
      s2 = (c24 + a * (c29 + a * (c32 + a * (c33 + c35 * a)))) / c37
      s3 = (c19 + a * (c25 + a * (c28 + c31 * a))) / c37
      s4 = (c20 + a * (c27 + c34 * a) + c * (c22 + a * (c30 + c36 * a))) / c38
      s5 = (c13 + c21 * a + c * (c18 + c26 * a)) / c37
      s6 = (c15 + c * (c23 + c16 * c)) / c38
      fn_val = fn_val + t * (one + half * t * s1 - b * c * (s1 - b *   &
               (s2 - b * (s3 - b * (s4 - b * (s5 - b * s6))))))
      IF (ABS(q / fn_val - one) > e) RETURN
   END DO
   
   WRITE(lu_stderr,'(a)') 'statistical: PPCHI2: error: Max. number of iterations exceeded'
   stop
   
   6 RETURN
END FUNCTION ppchi2
!-------------------------------------------------------------------------------

!===============================================================================
function gammad(x, p) RESULT(fn_val)
!===============================================================================
!  ALGORITHM AS239  APPL. STATIST. (1988) VOL. 37, NO. 3

!  Computation of the Incomplete Gamma Integral

!  Auxiliary functions required: ALOGAM = logarithm of the gamma
!  function, and ALNORM = algorithm AS66

! ELF90-compatible version by Alan Miller
! Latest revision - 27 October 2000

! N.B. Argument IFAULT has been removed
! AJN: Remove need for alogam by using log_gamma, intrinsic in Fortran 2008 and
!      implemented in most recent compilers.

   IMPLICIT NONE
   INTEGER, PARAMETER    :: dp = rs
   REAL (dp), INTENT(IN) :: x, p
   REAL (dp)             :: fn_val
   
   ! Local variables
   REAL (dp)             :: pn1, pn2, pn3, pn4, pn5, pn6, arg, c, rn, a, b, an
   REAL (dp), PARAMETER  :: zero = 0.d0, one = 1.d0, two = 2.d0, &
                            oflo = 1.d+37, three = 3.d0, nine = 9.d0, &
                            tol = 1.d-14, xbig = 1.d+8, plimit = 1000.d0, &
                            elimit = -88.d0
! EXTERNAL alogam, alnorm
   
   fn_val = zero
   
!  Check that we have valid values for X and P
   IF (p <= zero .OR. x < zero) THEN
     WRITE(lu_stderr,'(a)') 'statistical: gammad(AS239): error: Either p <= 0 or x < 0'
     stop !RETURN
   END IF
   IF (x == zero) RETURN
   
!  Use a normal approximation if P > PLIMIT
   IF (p > plimit) THEN
     pn1 = three * SQRT(p) * ((x / p) ** (one / three) + one /(nine * p) - one)
     fn_val = alnorm(pn1, .false.)
     RETURN
   END IF
   
!  If X is extremely large compared to P then set fn_val = 1
   IF (x > xbig) THEN
     fn_val = one
     RETURN
   END IF
   
   IF (x <= one .OR. x < p) THEN
     
!  Use Pearson's series expansion.
!  (Note that P is not large enough to force overflow in ALOGAM).
!  No need to test IFAULT on exit since P > 0.
!     arg = p * LOG(x) - x - alogam(p + one, ifault)
     arg = p * LOG(x) - x - log_gamma(p + one)
     c = one
     fn_val = one
     a = p
     40   a = a + one
     c = c * x / a
     fn_val = fn_val + c
     IF (c > tol) GO TO 40
     arg = arg + LOG(fn_val)
     fn_val = zero
     IF (arg >= elimit) fn_val = EXP(arg)
     
   ELSE
     
!  Use a continued fraction expansion
!     arg = p * LOG(x) - x - alogam(p, ifault)
     arg = p * LOG(x) - x - log_gamma(p)
     a = one - p
     b = a + x + one
     c = zero
     pn1 = one
     pn2 = x
     pn3 = x + one
     pn4 = x * b
     fn_val = pn3 / pn4
     60   a = a + one
     b = b + two
     c = c + one
     an = a * c
     pn5 = b * pn3 - an * pn1
     pn6 = b * pn4 - an * pn2
     IF (ABS(pn6) > zero) THEN
       rn = pn5 / pn6
       IF (ABS(fn_val - rn) <= MIN(tol, tol * rn)) GO TO 80
       fn_val = rn
     END IF
     
     pn1 = pn3
     pn2 = pn4
     pn3 = pn5
     pn4 = pn6
     IF (ABS(pn5) >= oflo) THEN
       
!  Re-scale terms in continued fraction if terms are large
       pn1 = pn1 / oflo
       pn2 = pn2 / oflo
       pn3 = pn3 / oflo
       pn4 = pn4 / oflo
     END IF
     GO TO 60
     80   arg = arg + LOG(fn_val)
     fn_val = one
     IF (arg >= elimit) fn_val = one - EXP(arg)
   END IF
   
   RETURN
END FUNCTION gammad
!-------------------------------------------------------------------------------

!===============================================================================
subroutine ppnd16 (p, normal_dev, ifault)
!===============================================================================
! ALGORITHM AS241  APPL. STATIST. (1988) VOL. 37, NO. 3

! Produces the normal deviate Z corresponding to a given lower
! tail area of P; Z is accurate to about 1 part in 10**16.

! The hash sums below are the sums of the mantissas of the
! coefficients.   They are included for use in checking
! transcription.

! This ELF90-compatible version by Alan Miller - 20 August 1996
! N.B. The original algorithm is as a function; this is a subroutine

   IMPLICIT NONE
   
   INTEGER, PARAMETER      :: dp = rs
   REAL (dp), INTENT(IN)   :: p
   INTEGER, INTENT(OUT)    :: ifault
   REAL (dp), INTENT(OUT)  :: normal_dev
   
! Local variables
   REAL (dp) :: zero = 0.d0, one = 1.d0, half = 0.5d0,  &
                split1 = 0.425d0, split2 = 5.d0, const1 = 0.180625d0, &
                const2 = 1.6d0, q, r
   
! Coefficients for P close to 0.5
   REAL (dp) :: a0 = 3.3871328727963666080D0, &
                a1 = 1.3314166789178437745D+2, &
                a2 = 1.9715909503065514427D+3, &
                a3 = 1.3731693765509461125D+4, &
                a4 = 4.5921953931549871457D+4, &
                a5 = 6.7265770927008700853D+4, &
                a6 = 3.3430575583588128105D+4, &
                a7 = 2.5090809287301226727D+3, &
                b1 = 4.2313330701600911252D+1, &
                b2 = 6.8718700749205790830D+2, &
                b3 = 5.3941960214247511077D+3, &
                b4 = 2.1213794301586595867D+4, &
                b5 = 3.9307895800092710610D+4, &
                b6 = 2.8729085735721942674D+4, &
                b7 = 5.2264952788528545610D+3
! HASH SUM AB    55.8831928806149014439
   
! Coefficients for P not close to 0, 0.5 or 1.
   REAL (dp) :: c0 = 1.42343711074968357734D0, &
                c1 = 4.63033784615654529590D0, &
                c2 = 5.76949722146069140550D0, &
                c3 = 3.64784832476320460504D0, &
                c4 = 1.27045825245236838258D0, &
                c5 = 2.41780725177450611770D-1, &
                c6 = 2.27238449892691845833D-2, &
                c7 = 7.74545014278341407640D-4, &
                d1 = 2.05319162663775882187D0, &
                d2 = 1.67638483018380384940D0, &
                d3 = 6.89767334985100004550D-1, &
                d4 = 1.48103976427480074590D-1, &
                d5 = 1.51986665636164571966D-2, &
                d6 = 5.47593808499534494600D-4, &
                d7 = 1.05075007164441684324D-9
! HASH SUM CD    49.33206503301610289036
   
! Coefficients for P near 0 or 1.
   REAL (dp) :: e0 = 6.65790464350110377720D0, &
                e1 = 5.46378491116411436990D0, &
                e2 = 1.78482653991729133580D0, &
                e3 = 2.96560571828504891230D-1, &
                e4 = 2.65321895265761230930D-2, &
                e5 = 1.24266094738807843860D-3, &
                e6 = 2.71155556874348757815D-5, &
                e7 = 2.01033439929228813265D-7, &
                f1 = 5.99832206555887937690D-1, &
                f2 = 1.36929880922735805310D-1, &
                f3 = 1.48753612908506148525D-2, &
                f4 = 7.86869131145613259100D-4, &
                f5 = 1.84631831751005468180D-5, &
                f6 = 1.42151175831644588870D-7, &
                f7 = 2.04426310338993978564D-15
! HASH SUM EF    47.52583317549289671629
   
   ifault = 0
   q = p - half
   IF (ABS(q) <= split1) THEN
      r = const1 - q * q
      normal_dev = q * (((((((a7*r + a6)*r + a5)*r + a4)*r + a3)*r + a2)*r + a1)*r + a0) / &
               (((((((b7*r + b6)*r + b5)*r + b4)*r + b3)*r + b2)*r + b1)*r + one)
      RETURN
   ELSE
      IF (q < zero) THEN
         r = p
      ELSE
         r = one - p
      END IF
      IF (r <= zero) THEN
         ifault = 1
         normal_dev = zero
         RETURN
      END IF
      r = SQRT(-LOG(r))
      IF (r <= split2) THEN
         r = r - const2
         normal_dev = (((((((c7*r + c6)*r + c5)*r + c4)*r + c3)*r + c2)*r + c1)*r + c0) / &
                  (((((((d7*r + d6)*r + d5)*r + d4)*r + d3)*r + d2)*r + d1)*r + one)
      ELSE
         r = r - split2
         normal_dev = (((((((e7*r + e6)*r + e5)*r + e4)*r + e3)*r + e2)*r + e1)*r + e0) / &
                  (((((((f7*r + f6)*r + f5)*r + f4)*r + f3)*r + f2)*r + f1)*r + one)
      END IF
      IF (q < zero) normal_dev = - normal_dev
      RETURN
   END IF
END SUBROUTINE ppnd16
!-------------------------------------------------------------------------------

!===============================================================================
function alnorm( x, upper ) RESULT( fn_val )
!===============================================================================
!  Algorithm AS66 Applied Statistics (1973) vol.22, no.3

!  Evaluates the tail area of the standardised normal curve
!  from x to infinity if upper is .true. or
!  from minus infinity to x if upper is .false.

! ELF90-compatible version by Alan Miller
! Latest revision - 29 November 2001

   IMPLICIT NONE
   INTEGER, PARAMETER     ::  dp = rs
   REAL(DP), INTENT(IN)   ::  x
   LOGICAL,   INTENT(IN)  ::  upper
   REAL(DP)               ::  fn_val

!  Local variables
   REAL(DP), PARAMETER   ::  zero=0.0_DP, one=1.0_DP, half=0.5_DP, con=1.28_DP
   REAL(DP)              ::  z, y
   LOGICAL               ::  up

!  Machine dependent constants
   REAL(DP), PARAMETER  ::  ltone = 7.0_DP, utzero = 18.66_DP
   REAL(DP), PARAMETER  ::  p = 0.398942280444_DP, q = 0.39990348504_DP,   &
                            r = 0.398942280385_DP, a1 = 5.75885480458_DP,  &
                            a2 = 2.62433121679_DP, a3 = 5.92885724438_DP,  &
                            b1 = -29.8213557807_DP, b2 = 48.6959930692_DP, &
                            c1 = -3.8052E-8_DP, c2 = 3.98064794E-4_DP,     &
                            c3 = -0.151679116635_DP, c4 = 4.8385912808_DP, &
                            c5 = 0.742380924027_DP, c6 = 3.99019417011_DP, &
                            d1 = 1.00000615302_DP, d2 = 1.98615381364_DP,  &
                            d3 = 5.29330324926_DP, d4 = -15.1508972451_DP, &
                            d5 = 30.789933034_DP

   up = upper
   z = x
   IF( z < zero ) THEN
      up = .NOT. up
      z = -z
   END IF
   IF( z <= ltone  .OR.  (up  .AND.  z <= utzero) ) THEN
      y = half*z*z
      IF( z > con ) THEN
         fn_val = r*EXP( -y )/(z+c1+d1/(z+c2+d2/(z+c3+d3/(z+c4+d4/(z+c5+d5/(z+c6))))))
      ELSE
         fn_val = half - z*(p-q*y/(y+a1+b1/(y+a2+b2/(y+a3))))
      END IF
   ELSE
      fn_val = zero
   END IF

   IF( .NOT. up ) fn_val = one - fn_val
   RETURN
END FUNCTION alnorm
!-------------------------------------------------------------------------------


!===============================================================================
!-------------------------------------------------------------------------------
!  Circular statistical functions
!-------------------------------------------------------------------------------
!===============================================================================

!===============================================================================
function circ_mean(angle,degrees)
!===============================================================================
!  Returns the circular mean of a set of angles.  Input is a column vector
!  of arbitrary length.  Default is for input in degrees.
   implicit none
   real(rs),intent(in) :: angle(:)
   real(rs)            :: circ_mean
   logical,intent(in),optional :: degrees
   logical                     :: radians
   
   radians = .false.
   if (present(degrees)) radians = .not.degrees
   
   if (radians) then
      circ_mean = atan2(sum( sin( angle(1:size(angle)) ) ), &
                        sum( cos( angle(1:size(angle)) ) ) )
   else
      circ_mean = atan2(sum( sin( angle(1:size(angle))*pi/180._rs ) ),  &
                        sum( cos( angle(1:size(angle))*pi/180._rs ) ) )
      circ_mean = circ_mean*180._rs/pi
   endif
   
   return
end function circ_mean
!-------------------------------------------------------------------------------

!===============================================================================
subroutine circ_mean_bootstrap(angle,mu,mu1,mu2,B,degrees)
!===============================================================================
!  Calculate a sample mean using a bootstrap technique.  Assumes a symmetric
!  distribution.  See Fisher, Statistical analysis of circular data, §4.4
!  Optionally specify number of bootstrap samples to take.
   implicit none
   real(rs),intent(in)  :: angle(:)
   real(rs),intent(out) :: mu
   real(rs),intent(out),optional :: mu1,mu2
   integer,intent(in),optional   :: B
   logical,intent(in),optional   :: degrees
   logical :: degrees_in
   
   degrees_in = .true.
   if (present(degrees)) degrees_in = degrees
   
   if (size(angle) <= 9) then  ! For small N, calculate all n**n samples
      call circ_mean_bootstrap_smallN(angle,mu,mu1,mu2,degrees=degrees_in)
   else
      write(lu_stderr,'(a)') 'statistical: circ_mean_bootstrap: only implemented for n<=9 at the moment.'
      stop
   endif
   
   return
end subroutine circ_mean_bootstrap
!-------------------------------------------------------------------------------

!===============================================================================
subroutine circ_mean_bootstrap_smallN(angle,mu,mu1,mu2,degrees,force)
!===============================================================================
!  For a small sample, calculate the bootstrap mean and confidence interval
!  by taking all possible subsamples of the data (= n**n).
!  If n==9, then n**n real*4s takes up ~1.5 GB.  Use this routine with care on
!  lesser or shared machines!  Hence a warning is displayed, unless the force=.true.
!  option is employed on the command line.
   implicit none
   real(rs),intent(in)  :: angle(:)
   real(rs),intent(out) :: mu,mu1,mu2
   logical,intent(in),optional :: degrees,force
   logical :: force_in
   real(rs) :: conversion
   real(r4) :: Csum,Ssum
   real(r4),allocatable :: mean(:)  ! Use single preision for storage of bootstrap samples
   real(rs),allocatable :: rangle(:)
   integer :: n,i,i1,i2,i3,i4,i5,i6,i7,i8,i9
   
!  We're potentially allocating <=1.5 GB of memory, so make sure we know what
!  we're doing.
   force_in = .false.
   if (present(force)) force_in = force
   if (.not.force_in) then
      if (size(angle) >= 8) then
         write(*,'(a,f0.0,a)') 'circ_mean_bootstrap_smallN is about to allocate ',&
                  8.*size(angle)**size(angle)/(1024.**2),' MB of memory.'
         write(*,'(a)') 'Hit return to continue, or ^C to abort.'
         read(*,*)
      endif
   endif
   
!  Check we're only testing with small sample size
   n = size(angle)
   if (n > 9 .or. n == 1) then
      write(lu_stderr,'(a)') &
      'statistical: circ_mean_bootstrap_smallN: number of samples too large for small N subroutine, or n is 1.'
      stop
   endif

!  Allocate memory for all n**n sample means: here we go!
   allocate(mean(n**n))
   
!  Convert angles into radians if necesary
   allocate(rangle(n))
   conversion = pi/180._rs
   if (present(degrees)) then
      if (.not.degrees) conversion = 1._rs
   endif
   rangle = conversion*angle
   
!  Calculate sample mean for all possible bootstrap samples
   i = 1
   if (n == 2) then
      do i1=1,n
       do i2=1,n
               Csum = cos(rangle(i1)) + cos(rangle(i2))
               Ssum = sin(rangle(i1)) + sin(rangle(i2))
               mean(i) = atan2(Ssum,Csum)
               i = i + 1
       enddo
      enddo
      
   else if (n == 3) then
      do i1=1,n
       do i2=1,n
        do i3=1,n
               Csum = cos(rangle(i1)) + cos(rangle(i2)) + cos(rangle(i3))
               Ssum = sin(rangle(i1)) + sin(rangle(i2)) + sin(rangle(i3))
               mean(i) = atan2(Ssum,Csum)
               i = i + 1
        enddo
       enddo
      enddo
      
   else if (n == 4) then
      do i1=1,n
       write(*,'(a,i0,a)') 'Executing ',i1,' of 8 passes.'
       do i2=1,n
        do i3=1,n
         do i4=1,n
                 Csum = cos(rangle(i1)) + cos(rangle(i2)) + cos(rangle(i3)) + &
                        cos(rangle(i4))
                 Ssum = sin(rangle(i1)) + sin(rangle(i2)) + sin(rangle(i3)) + &
                        sin(rangle(i4))
                 mean(i) = atan2(Ssum,Csum)
                 i = i + 1
         enddo
        enddo
       enddo
      enddo

   
   else if (n == 5) then
      do i1=1,n
       write(*,'(a,i0,a)') 'Executing ',i1,' of 8 passes.'
       do i2=1,n
        do i3=1,n
         do i4=1,n
          do i5=1,n
                 Csum = cos(rangle(i1)) + cos(rangle(i2)) + cos(rangle(i3)) + &
                        cos(rangle(i4)) + cos(rangle(i5))
                 Ssum = sin(rangle(i1)) + sin(rangle(i2)) + sin(rangle(i3)) + &
                        sin(rangle(i4)) + sin(rangle(i5))
                 mean(i) = atan2(Ssum,Csum)
                 i = i + 1
          enddo
         enddo
        enddo
       enddo
      enddo

   else if (n == 6) then
      do i1=1,n
       write(*,'(a,i0,a)') 'Executing ',i1,' of 8 passes.'
       do i2=1,n
        do i3=1,n
         do i4=1,n
          do i5=1,n
           do i6=1,n
                 Csum = cos(rangle(i1)) + cos(rangle(i2)) + cos(rangle(i3)) + &
                        cos(rangle(i4)) + cos(rangle(i5)) + cos(rangle(i6))
                 Ssum = sin(rangle(i1)) + sin(rangle(i2)) + sin(rangle(i3)) + &
                        sin(rangle(i4)) + sin(rangle(i5)) + sin(rangle(i6))
                 mean(i) = atan2(Ssum,Csum)
                 i = i + 1
           enddo
          enddo
         enddo
        enddo
       enddo
      enddo

   else if (n == 7) then
      do i1=1,n
       write(*,'(a,i0,a)') 'Executing ',i1,' of 8 passes.'
       do i2=1,n
        do i3=1,n
         do i4=1,n
          do i5=1,n
           do i6=1,n
            do i7=1,n
                 Csum = cos(rangle(i1)) + cos(rangle(i2)) + cos(rangle(i3)) + &
                        cos(rangle(i4)) + cos(rangle(i5)) + cos(rangle(i6)) + &
                        cos(rangle(i7))
                 Ssum = sin(rangle(i1)) + sin(rangle(i2)) + sin(rangle(i3)) + &
                        sin(rangle(i4)) + sin(rangle(i5)) + sin(rangle(i6)) + &
                        sin(rangle(i7))
                 mean(i) = atan2(Ssum,Csum)
                 i = i + 1
            enddo
           enddo
          enddo
         enddo
        enddo
       enddo
      enddo

   else if (n == 8) then
      do i1=1,n
       write(*,'(a,i0,a)') 'Executing ',i1,' of 8 passes.'
       do i2=1,n
        do i3=1,n
         do i4=1,n
          do i5=1,n
           do i6=1,n
            do i7=1,n
             do i8=1,n
                 Csum = cos(rangle(i1)) + cos(rangle(i2)) + cos(rangle(i3)) + &
                        cos(rangle(i4)) + cos(rangle(i5)) + cos(rangle(i6)) + &
                        cos(rangle(i7)) + cos(rangle(i8)) 
                 Ssum = sin(rangle(i1)) + sin(rangle(i2)) + sin(rangle(i3)) + &
                        sin(rangle(i4)) + sin(rangle(i5)) + sin(rangle(i6)) + &
                        sin(rangle(i7)) + sin(rangle(i8)) 
                 mean(i) = atan2(Ssum,Csum)
                 i = i + 1
             enddo
            enddo
           enddo
          enddo
         enddo
        enddo
       enddo
      enddo
         
   else if (n == 9) then
      do i1=1,n
       write(*,'(a,i0,a)') 'Executing ',i1,' of 9 passes.'
       do i2=1,n
        do i3=1,n
         do i4=1,n
          do i5=1,n
           do i6=1,n
            do i7=1,n
             do i8=1,n
              do i9=1,n
                 Csum = cos(rangle(i1)) + cos(rangle(i2)) + cos(rangle(i3)) + &
                        cos(rangle(i4)) + cos(rangle(i5)) + cos(rangle(i6)) + &
                        cos(rangle(i7)) + cos(rangle(i8)) + cos(rangle(i9))
                 Ssum = sin(rangle(i1)) + sin(rangle(i2)) + sin(rangle(i3)) + &
                        sin(rangle(i4)) + sin(rangle(i5)) + sin(rangle(i6)) + &
                        sin(rangle(i7)) + sin(rangle(i8)) + sin(rangle(i9))
                 mean(i) = atan2(Ssum,Csum)
                 i = i + 1
              enddo
             enddo
            enddo
           enddo
          enddo
         enddo
        enddo
       enddo
      enddo
   endif
   
   mu = sum(mean)/(conversion*n**n)
   
   return
end subroutine circ_mean_bootstrap_smallN
!-------------------------------------------------------------------------------

!===============================================================================
function circ_res_length(angle,degrees)
!===============================================================================
!  Returns the resultant length of the angles.  Input is a column vector of
!  arbitrary length.  Default is for input in degrees.
   implicit none
   real(rs),intent(in) :: angle(:)
   real(rs)            :: circ_res_length
   logical,intent(in),optional :: degrees
   logical                     :: radians
   
   radians = .false.
   if (present(degrees)) radians = .not.degrees
   
   if (radians) then
      circ_res_length = sqrt((sum(sin(angle(1:size(angle)))))**2 + &
                             (sum(cos(angle(1:size(angle)))))**2) / real(size(angle))
   else
      circ_res_length = sqrt((sum( sin(angle(1:size(angle))*pi/180._rs)) )**2 + &
                (sum( cos(angle(1:size(angle))*pi/180._rs)) )**2) / real(size(angle))
   endif

   return
end function circ_res_length
!-------------------------------------------------------------------------------

!===============================================================================
function circ_sd(angle,degrees)
!===============================================================================
!  Returns the circular standard deviation.  Input is a column vector of 
!  arbitrary length.  Default is for input in degrees.
   implicit none
   real(rs),intent(in) :: angle(:)
   real(rs)            :: circ_sd
   logical,intent(in),optional :: degrees
   logical                     :: degrees_in
   
   degrees_in = .true.
   if (present(degrees)) degrees_in = degrees
   
   circ_sd = sqrt(-2._rs*log(circ_res_length(angle,degrees=degrees_in)))
   
   return
end function circ_sd
!===============================================================================

!===============================================================================
function circ_variance(angle,degrees)
!===============================================================================
!  Returns the circular variance of a set of angles.  Input is a column vector
!  of arbitrary length.  Default is for input in degrees.
   implicit none
   real(rs),intent(in) :: angle(:)
   real(rs)            :: circ_variance,mean
   logical,intent(in),optional :: degrees
   logical                     :: radians
   
   radians = .false.
   if (present(degrees)) radians = .not.degrees
   
   if (radians) then
      mean = circ_mean(angle,degrees=.false.)
      circ_variance = 1._rs - sum(cos(angle - mean))/real(size(angle))
   else
      mean = circ_mean(angle,degrees=.true.) * pi/180._rs
      circ_variance = 1._rs - sum(cos(angle*pi/180._rs - mean))/real(size(angle))
   endif
   
   return
end function circ_variance
!-------------------------------------------------------------------------------

!===============================================================================
function circ_uncent_trig_moment(angle,p,degrees)
!===============================================================================
!  Returns the uncentered pth trigonometric moment.  Input is a column vector of
!  arbitrary length.  Default is for input in degrees.
   implicit none
   real(rs),intent(in) :: angle(:)
   integer,intent(in)  :: p
   complex(cs)         :: circ_uncent_trig_moment
   real(rs) :: Cp,Sp
   integer  :: n
   logical,intent(in),optional :: degrees
   real(rs) :: conversion
   
!  Check for valid p
   if (p < 1) then
      write(lu_stderr,'(a)') 'statistical: circ_uncent_trig_moment: p must be > 0.'
      stop
   endif
   
!  Convert from degrees to radians unless otherwise
   conversion = pi/180._rs
   if (present(degrees)) then
      if (.not.degrees) conversion = 1._rs
   endif
   
   n = size(angle)
   Cp = sum(cos(real(p)*angle(1:n)*conversion))/real(n)
   Sp = sum(sin(real(p)*angle(1:n)*conversion))/real(n)
   
   circ_uncent_trig_moment = cmplx(Cp, Sp)
   
   return
end function circ_uncent_trig_moment
!-------------------------------------------------------------------------------

!===============================================================================
function circ_cent_trig_moment(angle,p,degrees)
!===============================================================================
!  Returns the centered pth trigonometric moment.  Input is a column vector of
!  arbitrary length.  Default is for input in degrees.
   implicit none
   real(rs),intent(in) :: angle(:)
   integer, intent(in) :: p
   complex(cs)         :: circ_cent_trig_moment
   real(rs) :: Cp,Sp,mean,conversion
   integer  :: n
   logical,intent(in),optional :: degrees
   
!  Check for valid p
   if (p < 1) then
      write(lu_stderr,'(a)') 'statistical: circ_cent_trig_moment: P must be > 0.'
      stop
   endif
   
!  Convert from degrees to radians unless otherwise desired
   conversion = pi/180._rs
   if (present(degrees)) then;  if (.not.degrees) conversion = 1._rs; endif
   
   n = size(angle)
   mean = circ_mean(angle,degrees=degrees)
   Cp = sum(cos(real(p)*conversion*(angle(1:n) - mean)))/real(n)
   Sp = sum(sin(real(p)*conversion*(angle(1:n) - mean)))/real(n)
   
   circ_cent_trig_moment = cmplx(Cp, Sp)
   
   return
end function circ_cent_trig_moment
!-------------------------------------------------------------------------------

!===============================================================================
function circ_dispers(angle,degrees)
!===============================================================================
!  Returns the circular dispersion of a set of angles.  Input is a column vector
!  of arbitrary length.  Default is for input in degrees.
   implicit none
   real(rs),intent(in) :: angle(:)
   real(rs)            :: circ_dispers
   real(rs) :: mean,rho2,R,conversion
   integer :: n
   logical,intent(in),optional :: degrees
   
!  By default, convert from degrees to radians for calculations
   conversion = pi/180._rs
   if (present(degrees)) then
      if (.not.degrees) conversion = 1._rs
   endif
   
   mean = circ_mean(angle,degrees=degrees)
   n = size(angle)
   R = sum(cos(conversion*(angle(1:n) - mean)))/real(n)
   rho2 = sum(cos(2._rs*conversion*(angle(1:n) - mean)))/real(n)
   circ_dispers = (1._rs - rho2)/(2._rs*R**2)
   
   return
end function circ_dispers
!-------------------------------------------------------------------------------

!===============================================================================
function circ_correl(a,b,degrees)
!===============================================================================
!  Returns the circular correlation between two sets of angles.  Input is two
!  column vectors of arbitrary (but the same) length.  Default is for input
!  in dgerees.
   implicit none
   real(rs),intent(in) :: a(:), b(:)
   real(rs)            :: circ_correl, amean, bmean
   logical,intent(in),optional :: degrees
   logical                     :: radians
   
!  Check arrays are same length
   if (size(a) /= size(b)) then
      write(0,'(a)') 'Error: statistical: circ_correl: Input arrays are not same length.'
      stop
   endif
   
   radians = .false.
   if (present(degrees)) radians = .not.degrees
   
   if (radians) then
      amean = circ_mean(a,degrees=.false.)
      bmean = circ_mean(a,degrees=.false.)
      circ_correl = sum( sin(a - amean) * sin(b - bmean) )/ &
            sqrt( sum( sin(a-amean)*sin(a-amean)*sin(b-bmean)*sin(b-bmean) ) )
   else
      amean = circ_mean(a,degrees=.true.) * pi/180._rs
      bmean = circ_mean(b,degrees=.true.) * pi/180._rs
      circ_correl = sum( sin(a*pi/180._rs - amean) * sin(b*pi/180._rs - bmean) )/ &
            sqrt( sum( sin(a*pi/180._rs-amean)*sin(a*pi/180._rs-amean) * &
                  sin(b*pi/180._rs-bmean)*sin(b*pi/180._rs-bmean) ) )   
   endif
   
   return
end function circ_correl
!-------------------------------------------------------------------------------

!===============================================================================
function circ_von_mieses(theta,mu,kappa,degrees)
!===============================================================================
!  Returns the value of the von Mieses distribution for given mean and 'pointiness'
   implicit none

   real(rs),intent(in) :: theta,mu,kappa
   real(rs)            :: circ_von_mieses,conversion
   logical,intent(in),optional :: degrees
   
   conversion = pi/180._rs
   if (present(degrees)) then; if (.not.degrees) conversion = 1._rs; endif
   
   circ_von_mieses = exp(kappa*cos(conversion*(theta - mu)))/(2._rs*pi*BessI0(kappa))

   return
   
   contains
      function BessI0(x)
      !  Evaluates the modified Bessel function of the first kind at x.
      !  Polynomial approximation taken from Abramowitz & Stegun, 1964, Handbook
      !   of mathematical functions.  See also Numerical Recipes §6.6
         implicit none
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

end function circ_von_mieses
!-------------------------------------------------------------------------------

!===============================================================================
function circ_test_random_orient(theta,theta0,alpha,degrees) result(pass)
!===============================================================================
! Test whether a set of orientations (i.e., 180-degrees ambiguous), theta_i,
! point in a certain direction, theta0, statistically differently from random,
! using the V-test as modified by Rayleigh
! (test 95 in: 100 Statistical Tests. G.K. Kanji, SAGE Publications.).
   implicit none
   real(rs), intent(in) :: theta(:), theta0
   real(rs), intent(in), optional :: alpha
   logical, optional, intent(in) :: degrees
   real(rs) :: conversion, alpha_in
   real(rs) :: r, mean, V
   logical :: radians, pass
   integer :: n
   
   n = size(theta)
   if (n < 5) then
      write(0,'(a)') 'statistical: circ_test_random_orient: Error: Cannot ' // &
         'calculate significance for samples < 5'
      stop
   endif
   
   radians = .false.
   alpha_in = 0.05_rs
   if (present(alpha)) alpha_in = alpha
   
   if (present(degrees)) radians = .not.degrees
   conversion = pi/180._rs
   if (radians) conversion = 1._rs

   mean = circ_mean(2._rs*theta, degrees=.not.radians)
   r = circ_res_length(2._rs*theta, degrees=.not.radians)
   V = 2._rs*sqrt(2._rs*size(theta))*r*cos(conversion*mean - 2._rs*conversion*theta0)
   
   ! Lookup values in table and decide on significance
   pass = .false.
   if (V > v_test_table()) pass = .true.
   
   contains
      function v_test_table() result(value)
         integer :: i
         real(rs) :: value
         integer, parameter, dimension(33) :: table_n = &
            (/ 5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26, &
            27,28,29,30,40,50,60,70,100,500,1000 /)
         real(rs), parameter, dimension(6) :: table_alpha = &
            (/ 0.1, 0.05, 0.01, 0.005, 0.001, 0.0001 /)
         real(rs), parameter :: table(33,6) = transpose(reshape( (/ &
            1.3051_rs, 1.6524_rs, 2.2505_rs, 2.4459_rs, 2.7938_rs, 3.0825_rs, & 
            1.3009_rs, 1.6509_rs, 2.2640_rs, 2.4695_rs, 2.8502_rs, 3.2114_rs, & 
            1.2980_rs, 1.6499_rs, 2.2734_rs, 2.4858_rs, 2.8886_rs, 3.2970_rs, & 
            1.2958_rs, 1.6492_rs, 2.2803_rs, 2.4978_rs, 2.9164_rs, 3.3578_rs, & 
            1.2942_rs, 1.6484_rs, 2.2856_rs, 2.5070_rs, 2.9375_rs, 3.4034_rs, & 
            1.2929_rs, 1.6482_rs, 2.2899_rs, 2.5143_rs, 2.9540_rs, 3.4387_rs, & 
            1.2918_rs, 1.6479_rs, 2.2933_rs, 2.5201_rs, 2.9672_rs, 3.4669_rs, & 
            1.2909_rs, 1.6476_rs, 2.2961_rs, 2.5250_rs, 2.9782_rs, 3.4899_rs, & 
            1.2902_rs, 1.6474_rs, 2.2985_rs, 2.5290_rs, 2.9873_rs, 3.5091_rs, & 
            1.2895_rs, 1.6472_rs, 2.3006_rs, 2.5325_rs, 2.9950_rs, 3.5253_rs, & 
            1.2890_rs, 1.6470_rs, 2.3023_rs, 2.5355_rs, 3.0017_rs, 3.5392_rs, & 
            1.2885_rs, 1.6469_rs, 2.3039_rs, 2.5381_rs, 3.0075_rs, 3.5512_rs, & 
            1.2881_rs, 1.6467_rs, 2.3052_rs, 2.5404_rs, 3.0126_rs, 3.5617_rs, & 
            1.2877_rs, 1.6466_rs, 2.3064_rs, 2.5424_rs, 3.0171_rs, 3.5710_rs, & 
            1.2874_rs, 1.6465_rs, 2.3075_rs, 2.5442_rs, 3.0211_rs, 3.5792_rs, & 
            1.2871_rs, 1.6464_rs, 2.3085_rs, 2.5458_rs, 3.0247_rs, 3.5866_rs, & 
            1.2868_rs, 1.6464_rs, 2.3093_rs, 2.5473_rs, 3.0279_rs, 3.5932_rs, & 
            1.2866_rs, 1.6463_rs, 2.3101_rs, 2.5486_rs, 3.0308_rs, 3.5992_rs, & 
            1.2864_rs, 1.6462_rs, 2.3108_rs, 2.5498_rs, 3.0335_rs, 3.6047_rs, & 
            1.2862_rs, 1.6462_rs, 2.3115_rs, 2.5509_rs, 3.0359_rs, 3.6096_rs, & 
            1.2860_rs, 1.6461_rs, 2.3121_rs, 2.5519_rs, 3.0382_rs, 3.6142_rs, & 
            1.2858_rs, 1.6461_rs, 2.3127_rs, 2.5529_rs, 3.0402_rs, 3.6184_rs, & 
            1.2856_rs, 1.6460_rs, 2.3132_rs, 2.5538_rs, 3.0421_rs, 3.6223_rs, & 
            1.2855_rs, 1.6460_rs, 2.3136_rs, 2.5546_rs, 3.0439_rs, 3.6258_rs, & 
            1.2853_rs, 1.6459_rs, 2.3141_rs, 2.5553_rs, 3.0455_rs, 3.6292_rs, & 
            1.2852_rs, 1.6459_rs, 2.3145_rs, 2.5560_rs, 3.0471_rs, 3.6323_rs, & 
            1.2843_rs, 1.6456_rs, 2.3175_rs, 2.5610_rs, 3.0580_rs, 3.6545_rs, & 
            1.2837_rs, 1.6455_rs, 2.3193_rs, 2.5640_rs, 3.0646_rs, 3.6677_rs, & 
            1.2834_rs, 1.6454_rs, 2.3205_rs, 2.5660_rs, 3.0689_rs, 3.6764_rs, & 
            1.2831_rs, 1.6453_rs, 2.3213_rs, 2.5674_rs, 3.0720_rs, 3.6826_rs, & 
            1.2826_rs, 1.6452_rs, 2.3228_rs, 2.5699_rs, 3.0775_rs, 3.6936_rs, & 
            1.2818_rs, 1.6449_rs, 2.3256_rs, 2.5747_rs, 3.0877_rs, 3.7140_rs, & 
            1.2817_rs, 1.6449_rs, 2.3260_rs, 2.5752_rs, 3.0890_rs, 3.7165_rs /), (/6,33/)))
         ! The function just checks that the value you asked for alpha
         ! is near to one of those we offer--if not, you get an error
         do i = 1, size(table_alpha)
            if (abs(table_alpha(i) - alpha_in)/table_alpha(i) < 0.1_rs) then
               value = table(minloc(abs(n-table_n),1), i)
               return
            endif
         enddo
         write(0,'(a,f0.6,a)') 'statistical: circ_test_random_orient: v_test_table: ' &
            //'Error: Cannot find a value for alpha near that requested (',alpha_in,')'
         stop
      end function v_test_table
end function circ_test_random_orient
!-------------------------------------------------------------------------------

!===============================================================================



!_______________________________________________________________________________
end module statistical
