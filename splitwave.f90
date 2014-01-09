!===============================================================================
!  module splitwave defines a set of routines for creating waveforms and applying
!  splitting operators to them, either in the time or frequency domain.
!
!  The forward splitting operator routines require FFFTW3 and are based on code
!  by James Wookey <j.wookey@bristol.ac.uk>.
!
!===============================================================================
module splitwave
!===============================================================================

   implicit none

   ! Kind parameters
   integer, parameter, private :: r4 = selected_real_kind(6,37)
   integer, parameter, private :: r8 = selected_real_kind(15,307)
   integer, parameter, private :: rs = r8  ! Use double precision for reals

   ! Constants
   real(rs), parameter, private :: pi = 4._rs*atan2(1._rs,1._rs)

   ! Default values for routines
   real, parameter, private :: freq_default = 0.1_rs, &
                               delta_default = 0.05_rs, &
                               noise_default = 0.05_rs
   character(len=1), parameter, private :: wavetype_default = 'g'


   ! Input/output units
   integer, parameter, private :: lu_stdin = 5, &
                                  lu_stdout = 6, &
                                  lu_stderr = 0

   ! Hide internal routines
   private :: gaussian1, &
              gaussian2, &
              ricker, &
              sw_error, &
              sw_warning

contains

!===============================================================================
subroutine sw_misfit_ecs(C,az,inc,phi,dt,spol,misfit,t,phi_ecs,dt_ecs,t_scaled, &
   freq,delta,noise,wavetype)
!===============================================================================
!  Calculate the misfits between a set of splitting observations and those
!  predicted by a set of elastic constants.
!  Unless a layer thickness is specified, by default sw_misfit_ecs calculates
!  the splitting for the Cijs then normalises the dt to have the same maximum
!  as the observed splits.  Specifying the layer thickness turns this off.
!  INPUT:
!     C(6,6)       : Density-normalised (Aij) elastic constants / m^2.s^-2
!     az(:)        : Azimuth of rays defined in the CIJ_phasevels frame / deg
!     inc(:)       : Inclination of rays     "   "        "         "
!     phi(:),dt(:) : Splitting parameters of observations (ray frame) / deg, s
!     spol(:)      : Source polarisation of observations (ray frame) / deg
!  OPTIONAL INPUT:
!     t            : Thickness of an assumed layer across which the splits in the
!                    elastic constants are accrued / km
!     freq         : Dominant frequency of wave / Hz
!     noise        : Fractional amount of random noise
!     wavetype     : 'g'aussian (first-derivative) or 'r'icker
!  OUTPUT:
!     misfit(:)    : Array of misfits corresponding to each observation
!  OPTIONAL OUTPUT:
!     phi_ecs(:),dt_ecs(:) : Splits accrued in the elastic constants / deg, s
!     t_scaled     : Thickness used as scaling when not specifying t above / km
   use EmatrixUtils

   implicit none

   real(rs), intent(in) :: C(6,6), az(:), inc(:), phi(:), dt(:), spol(:)
   real(rs), intent(in), optional :: t, freq, delta, noise
   real(rs), intent(out) :: misfit(:)
   real(rs), intent(out), optional :: phi_ecs(:), dt_ecs(:), t_scaled
   real :: freq_in, delta_in, noise_in
   character(len=*), intent(in), optional :: wavetype
   character(len=1) :: wavetype_in
   real(rs), allocatable :: phi2(:), dt2(:)
   real(rs) :: vs1, vs2
   integer :: i, n

   ! Take wave parameters from arguments or set to standard values
   freq_in = freq_default
   noise_in = noise_default
   wavetype_in = wavetype_default
   if (present(freq)) freq_in = real(freq,r4)
   delta_in = 1./(freq_in*200.)
   if (present(noise)) noise_in = real(noise,r4)
   if (present(delta)) delta_in = real(delta,r4)
   if (present(wavetype)) wavetype_in = wavetype(1:1)

   ! Check array sizes
   n = size(az)
   if (size(inc) /= n .or. size(phi) /= n .or. size(dt) /= n .or. size(spol) /= n) &
      call sw_error('sw_misfit_ecs: Arrays for az, inc, phi, dt and spol not the same length')
   if (size(misfit) < n) &
      call sw_error('sw_misfit_ecs: Output array for misfit is too small')
   if (present(phi_ecs)) then
      if (size(phi_ecs) < n) &
         call sw_error('sw_misfit_ecs: Output array for phi_ecs is too small')
   endif
   if (present(dt_ecs)) then
      if (size(dt_ecs) < n) &
         call sw_error('sw_misfit_ecs: Output array for dt_ecs is too small')
   endif

   allocate(phi2(n), dt2(n))

   ! Calculate splitting and misfit for each observation
   do i=1,n
      call CIJ_phasevels(C, 1._rs, az(i), inc(i), pol=phi2(i), vs1=vs1, vs2=vs2)
      if (present(phi_ecs)) phi_ecs(i) = phi2(i)
      dt2(i) = 1._rs/vs2 - 1._rs/vs1
      if (present(dt_ecs)) dt_ecs(i) = dt2(i)
   enddo

   ! If we've specified a layer thickness set dt2, otherwise normalise to
   ! same maximum as input dt (hopefully similar range)
   if (present(t)) then
      dt2 = t*dt2
   else
      if (present(t_scaled)) t_scaled = maxval(dt)/maxval(dt2)
      dt2 = maxval(dt)*dt2/maxval(dt2)
   endif

   if (present(phi_ecs)) phi_ecs = phi2
   if (present(dt_ecs)) dt_ecs = dt2

   do i=1,n
      misfit(i) = sw_misfit(real(phi(i),r4), real(dt(i),r4), real(phi2(i),r4), &
         real(dt2(i),r4), real(spol(i),r4), freq=real(freq_in,r4), &
         delta=real(delta_in,r4), noise=real(noise_in,r4), wavetype=wavetype_in)
   enddo

   deallocate(phi2, dt2)

end subroutine sw_misfit_ecs
!-------------------------------------------------------------------------------

!===============================================================================
function sw_misfit(phi1,dt1,phi2,dt2,spol,freq,delta,noise,wavetype) result(misfit)
!===============================================================================
!  Calculates a misfit between two splitting operators.  This obviously depends
!  on several factors: the source polarisation, the frequency content of the
!  the split waves, and so on.
!  INPUT:
!     phi1,dt1,phi2,dt2 : Two splitting operators to compare
!     spol     : Source polarisation of wave / degrees
!  OPTIONAL INPUT:
!     freq     : Dominant frequency of wave
!     noise    : Fractional amount of random noise
!     wavetype : 'g'aussian (first-derivative) or 'r'icker
   use f90sac

   implicit none

   real :: misfit
   type(SACtrace) :: E,N,Z
   real, intent(in) :: phi1, dt1, phi2, dt2, spol
   real, intent(in), optional :: freq, delta, noise
   character(len=*), intent(in), optional :: wavetype
   real :: f, delta_in, noise_in, c(2,2), eig(2)
   character(len=1) :: wavetype_in
   logical :: debug = .false.  ! Set to .true. to output waves

   f = freq_default
   noise_in = noise_default
   wavetype_in = wavetype_default
   if (present(freq)) f = freq
   delta_in = 1./(f*200.)
   if (present(noise)) noise_in = noise
   if (present(delta)) delta_in = delta
   if (present(wavetype)) wavetype_in = wavetype(1:1)

   ! Make trace
   call sw_create_wave(E,N,Z, freq=f, delta=delta_in, noise=noise_in, &
      wavetype=wavetype_in, spol=spol)
   call debug_writetrace('wave')
   ! Apply the first splitting operator, then remove the second operator
   call sw_splitN(E,N,Z, 1, (/phi1/), (/ dt1/))
   call debug_writetrace('wave_s1')
   call sw_splitN(E,N,Z, 1, (/phi2/), (/-dt2/))
   call debug_writetrace('wave_s1_-s2')
   ! Find covariance matrix and get misfit as ratio of smaller to larger
   ! eigenvalues of covariance matrix
   call f90sac_covar2(N,E,c)
   eig = sw_eig2x2(c)
   misfit = minval(eig)/maxval(eig)

   ! Now do the same but the other way round
   call sw_create_wave(E,N,Z, freq=f, delta=delta_in, noise=noise_in, &
      wavetype=wavetype_in, spol=spol)
   call sw_splitN(E,N,Z, 1, (/phi2/), (/ dt2/))
   call debug_writetrace('wave_s2')
   call sw_splitN(E,N,Z, 1, (/phi1/), (/-dt1/))
   call debug_writetrace('wave_s2_-s1')
   call f90sac_covar2(N,e,c)
   eig = sw_eig2x2(c)
   misfit = (misfit + minval(eig)/maxval(eig))/2.
   call f90sac_deletetrace(E)
   call f90sac_deletetrace(N)
   call f90sac_deletetrace(Z)

contains

   subroutine debug_writetrace(str)
      implicit none
      character(len=*), intent(in) :: str
      if (debug) then
         call f90sac_writetrace(trim(str)//'.BHE', E)
         call f90sac_writetrace(trim(str)//'.BHN', N)
         call f90sac_writetrace(trim(str)//'.BHZ', Z)
      endif
   end subroutine debug_writetrace
end function sw_misfit
!-------------------------------------------------------------------------------

!==============================================================================
subroutine sw_create_wave(E,N,Z,freq,delta,noise,wavetype,spol)
!==============================================================================
!  Creates a synthetic waveform of a number of different kinds.
!  INPUT:
!     E,N,Z:    SACtrace types.  E and N are the East and North components
!     freq:     Desired dominant frequency of waveform / Hz
!     delta:    Sampling frequency of wave / s
!     noise:    Add white noise with max amplitude noise_in of the wave.
!     wavetype: String describing type of waveform desired
!     spol:     Source polarisation measured from N.
!  E,N,Z if not already allocated

   use f90sac

   implicit none

   type(SACtrace), intent(out)   :: E, N, Z
   real, intent(in), optional :: freq, delta, noise, spol
   character(len=*), intent(in), optional :: wavetype
   real(r4), allocatable :: trace(:), random_noise(:)
   real :: delta_in, spol_in, freq_in, noise_in
   integer :: npts
   character(len=1) :: wavetype_in

   ! Set defaults
   freq_in = 0.1
   delta_in = 0.05
   noise_in = 0.0
   wavetype_in = 'g'
   spol_in = 0.0

   ! Check for presence of arguments
   if (present(freq)) freq_in = real(freq)
   if (present(delta)) delta_in = real(delta)
   if (present(noise)) then
      if (noise < 0.) then
         call sw_warning('sw_create_wave: Value for noise must be greater than 0.  Applying no noise.')
         noise_in = 0._r4
      else
         noise_in = real(noise)
      endif
   endif
   if (present(spol)) spol_in = real(spol)
   if (present(delta)) delta_in = real(delta)
   if (present(wavetype)) then
      if (wavetype(1:1) /= 'g' .and. wavetype(1:1) /= 'G' .and. &
          wavetype(1:1) /= 'r' .and. wavetype(1:1) /= 'R') then
         call sw_error('sw_create_wave: Error: Unrecognised wave type: "'//trim(wavetype)//'"' &
            //'  Specify g(aussian1) or r(icker).')
      else
         wavetype_in = wavetype(1:1)
      endif
   endif

   ! Issue sanity check if 1/delta is below the Nyqvist frequency:
   if (delta_in > 1.0/(2.0*freq_in)) then
      call sw_warning('sw_create_wave: 1/delta is below the Nyqvist frequency.  Expect aliasing.')
   endif

   ! The waveform is 5/f s long: npts = 5/(delta * freq)
   npts = int(5.0 / (delta_in * freq_in)) + 1
   allocate(trace(npts))

   ! Create the waveform
   if (wavetype_in == 'g') call gaussian1(freq_in,delta_in,trace,npts)
   if (wavetype_in == 'r') call ricker(freq_in,delta_in,trace,npts)

   ! Now create the trace: starting with the N-component
   call f90sac_newtrace(2*(npts-1)+1,delta_in,N)

   N%kstnm = 'SWAV'
   N%o = 0.0
   N%b = 0.0
   N%e = real(2*(npts-1)+1) * delta_in

   N%cmpaz =  0.0 ; N%cmpinc = 90.0 ; N%kcmpnm ='BHN'

   call f90sac_setdate(N,3000, 1, 0, 0, 0, 0)
   call f90sac_setevent(N, -90.0, 0.0, 100.0)
   call f90sac_setstation(N, 0.0, 0.0, 0.0)
   N%az = 0.0  ;  N%baz = 180.0
   N%gcarc = 90.0

   ! Clone the trace and create the other components
   call f90sac_clonetrace(N,E)
   E%cmpaz = 90.0 ; E%cmpinc = 90.0 ; E%kcmpnm = 'BHE'

   call f90sac_clonetrace(N,Z)
   Z%cmpaz =  0.0 ; Z%cmpinc =  0.0 ; Z%kcmpnm = 'BHZ'

   N%trace = 0.0 ; E%trace = 0.0 ; Z%trace = 0.0
   N%trace(1:npts) = trace(1:npts) * cos(spol_in*real(pi)/180.0)
   E%trace(1:npts) = trace(1:npts) * sin(spol_in*real(pi)/180.0)

   ! Add noise
   if (noise_in > 1.0e-5) then
      allocate(random_noise(2*(npts-1)+1))
      call random_seed()

      call random_number(random_noise)
      random_noise = (random_noise*2.0 - 1.0) * noise_in
      N%trace = N%trace + random_noise

      call random_number(random_noise)
      random_noise = (random_noise*2.0 - 1.0) * noise_in
      E%trace = E%trace + random_noise

      call random_number(random_noise)
      random_noise = (random_noise*2.0 - 1.0) * noise_in
      Z%trace = Z%trace + random_noise

      deallocate(random_noise)
   endif

   ! Add analysis window picks
   N%a = 1.0/freq_in ; N%user0 = N%a ; N%user1 = N%a + delta_in
   E%a = 1.0/freq_in ; E%user0 = N%a ; E%user1 = E%a + delta_in
   Z%a = 1.0/freq_in ; Z%user0 = N%a ; Z%user1 = Z%a + delta_in

   N%f = 4.0/freq_in ; N%user2 = N%f ; N%user3 = N%f + delta_in
   E%f = 4.0/freq_in ; N%user2 = E%f ; N%user3 = E%f + delta_in
   Z%f = 4.0/freq_in ; Z%user2 = N%f ; Z%user3 = Z%f + delta_in

   deallocate(trace)

end subroutine sw_create_wave
!------------------------------------------------------------------------------

!==============================================================================
subroutine gaussian1(freq,delta,trace,npts)
!==============================================================================
!  Creates a wave which is the first derivative of the Gaussian.
!  y(t) = -(4/pi) f {t-(5/2f)} exp -[(4/pi) f {t-(5/2f)}]**2
!  Dominant frequency f, length of trace is 5/f, maxmimum amplitude of 1.
!  No check on delta is performed--one must be aware of this when calling
!  the routine in case it leads to aliasing.

   implicit none

   integer, intent(in) :: npts
   real(r4), intent(in) :: freq, delta
   real(r4), intent(out) :: trace(npts)
   real(r4) :: f, t, rpi
   integer :: i

   f = freq * 10.0_r4 / 3.0_r4

   rpi = real(pi,r4)
   do i=1,npts
      t = real(i-1) * delta
      trace(i) = -(4.0_r4/rpi) * f * (t - 5.0_r4/(2.0_r4*freq)) * &
                 exp( -((4.0_r4/rpi) * f * (t - 5.0_r4/(2.0_r4*freq)))**2 )
   enddo
   ! Normalise amplitude to 1
   trace(1:npts) = trace(1:npts)/(maxval(abs(trace(1:npts))))

end subroutine gaussian1
!------------------------------------------------------------------------------

!==============================================================================
subroutine gaussian2(freq,delta,trace,npts)
!==============================================================================
!  Creates a wave which is the second derivative of the Gaussian.
!  y(t) = (4/pi) [-2 (4/pi) f**2 {t-(5/2f)}**2 + 1] exp -[(4/pi) f {t-(5/2f)]**2
!  Dominant frequency f, length of trace is 5/f, maximum amplitude of 1.
!  No check on delta is performed--one must be aware of this when calling
!  the routine in case it leads to aliasing.

   implicit none

   integer, intent(in) :: npts
   real(r4), intent(in) :: freq, delta
   real(r4), intent(out) :: trace(npts)
   real(r4) :: t, rpi
   integer :: i

   rpi = real(pi,r4)
   do i=1,npts
      t = real(i-1) * delta
      trace(i) = (1._r4 - 2.0_r4*(4.0_r4/rpi) * freq**2 * (t - 5.0_r4/(2.0_r4*freq))**2) &
                 * exp (-((4.0_r4/rpi) * freq * (t - 5.0_r4/(2.0_r4*freq)))**2)
   enddo
   ! Normalise amplitude to 1
   trace(1:npts) = trace(1:npts)/(maxval(abs(trace(1:npts))))

end subroutine gaussian2
!------------------------------------------------------------------------------

!==============================================================================
subroutine ricker(freq,delta,trace,npts)
!==============================================================================
!  Creates a Ricker wave.
!  y(t) = (1 - 2 * t**2 * f**2 * pi**2) exp(-pi**2 * t**2 * f**2)
!  Dominant frequency f, length of trace is 5/f, maximum amplitude of 1.
!  No check on delta is performed--one must be aware of this when calling
!  the routine in case it leads to aliasing.

   implicit none

   integer, intent(in) :: npts
   real(r4), intent(in) :: freq, delta
   real(r4), intent(out) :: trace(npts)
   real(r4) :: t, rpi
   integer :: i

   rpi = real(pi,r4)
   do i = 1,npts
      t = real(i-1) * delta
      trace(i) = (1._r4 - 2._r4*rpi**2 * freq**2 * (t - 5._r4/(2._r4*freq))**2) * &
                 exp(-rpi**2 * freq**2 * (t - 5._r4/(2._r4*freq))**2)
   enddo
   ! Normalise amplitude to 1
   trace(1:npts) = trace(1:npts)/(maxval(abs(trace(1:npts))))

end subroutine ricker
!------------------------------------------------------------------------------

!===============================================================================
subroutine sw_splitN(t1,t2,t3,N,phi_in,dt_in,quiet)
!===============================================================================
!  Split a set of SACtrace waves in the time domain.
!  INPUT:
!     t1,t2,t3:     SACtrace types for the two horizontals (t1,t2) and the
!                   vertical.  These are changed in place.
!     N:            Number of splitting operators to apply.
!     phi_in,dt_in: Arrays of length N with splitting operators.
!     quiet:        If == .false., write out some information to stderr.
!
!  Splitting operators are applied in order

   use f90sac

   implicit none

   type(SACtrace), intent(inout) :: t1, t2, t3
   integer, intent(in) :: N
   real, intent(in) :: phi_in(N), dt_in(N)
   logical, optional, intent(in) :: quiet
   logical :: silent
   real :: dt, theta, sum_dt
   character(len=8) :: kcmpnm1, kcmpnm2
   integer :: i

   ! See if we're being quiet or not: default to silet
   silent = .true.
   if (present(quiet)) silent = quiet

   ! Check traces are same length
   if (t1%npts /= t2%npts) stop 'fdsplitwaveN: Traces must be same length!'

   !  Get old component names
   kcmpnm1 = t1 % kcmpnm   ;   kcmpnm2 = t2 % kcmpnm

   sum_dt = 0.
   ! Loop over the splitting parameters
   do i=1,N
      theta = real(phi_in(i), r4)
      dt = real(dt_in(i), r4)
      call f90sac_unwind(theta)

      ! Rotate traces, shift by dt, and rotate back in FD
      theta = theta - t1 % cmpaz
      call f90sac_rotate2d(t1,t2,theta,0)
      call f90sac_tshift(t2,dt)
      call f90sac_rotate2d(t1,t2,-theta,0)

      if (.not.silent) then
         write(lu_stderr,'(a,f7.2,a1,1x,f9.6)') &
         'Applied splitting with parameters phi,dt: ',theta + t1%cmpaz,',',dt
      endif
      sum_dt = sum_dt + dt
   enddo

   !  Set window markers for SHEBA by moving end window along by the applied tshift
   t1%f = t1%f + sum_dt - modulo(sum_dt,t1%delta)
   t1%user2 = t1%f   ;   t1%user3 = t1%f + t1%delta

   t2%f = t1%f  ;  t2%user2 = t1%user2  ;  t2%user3 = t1%user3
   t3%f = t1%f  ;  t3%user2 = t1%user2  ;  t3%user3 = t1%user3

   ! Replace component names
   t1%kcmpnm = kcmpnm1   ;  t2%kcmpnm = kcmpnm2

end subroutine sw_splitN
!-------------------------------------------------------------------------------

!==============================================================================
subroutine sw_fdsplitN(t1,t2,t3,N,phi_in,dt_in,quiet)
!==============================================================================
!  Split a set of SACtrace waves in the frequency domain.
!  INPUT:
!     t1,t2,t3:     SACtrace types for the two horizontals (t1,t2) and the
!                   vertical.  These are changed in place.
!     N:            Number of splitting operators to apply.
!     phi_in,dt_in: Arrays of length N with splitting operators.
!     quiet:        If == .false., write out some information to stderr.
!
!  Splitting operators are applied in order

   use FFFTW
   use f90sac

   implicit none

   type(SACtrace),intent(inout):: t1,t2,t3     ! The input traces
   integer,intent(in)  :: N                    ! Number of splits
   real,intent(in) :: phi_in(N),dt_in(N)       ! Splitting operators
   logical,optional,intent(in) :: quiet
   complex(r4),allocatable :: f1(:),f2(:)      ! The freq-domain traces
   logical             :: silent
   real                :: dt,theta
   real                :: sum_dt
   character(len=8)    :: kcmpnm1,kcmpnm2
   integer             :: np,npf,i

   ! See if we're being quiet or not: default to silet
   silent = .true.
   if (present(quiet)) silent = quiet

   ! Check traces are same length
   if (t1%npts /= t2%npts) call sw_error('sw_fdsplitN: Input traces must be same length!')

   ! Create fd traces
   np = t1%npts

   ! Calculate FFT: f1 and f2 are allocated
   call FFFTW_fwd(t1%trace,f1)
   call FFFTW_fwd(t2%trace,f2)
   npf = size(f1)
   if (size(f1) /= size(f2)) then
      call sw_error('sw_fdsplitN: Some problem with FFFTW: FFTs are not the same length')
   endif

   !  Get old component names
   kcmpnm1 = t1 % kcmpnm   ;   kcmpnm2 = t2 % kcmpnm

   sum_dt = 0.
   ! Loop over the splitting parameters
   do i=1,N
      theta = phi_in(i)
      dt = dt_in(i)
      call f90sac_unwind(theta)

      ! Rotate traces, shift by dt, and rotate back in FD
      theta = theta - t1 % cmpaz
      call sw_fdrotate2d(f1,f2,npf,theta)
      call sw_fdtshift(f2,npf,t1%delta,dt)
      call sw_fdrotate2d(f1,f2,npf,-theta)

      if (.not.silent) then
         write(lu_stderr,'(a,f7.2,a1,1x,f9.6)') &
         'Applied splitting with parameters phi,dt: ',theta + t1%cmpaz,',',dt
      endif
      sum_dt = sum_dt + dt
   enddo

   !  Set window markers for SHEBA by moving end window along by the applied tshift
   t1%f = t1%f + sum_dt - modulo(sum_dt,t1%delta)
   t1%user2 = t1%f   ;   t1%user3 = t1%f + t1%delta

   t2%f = t1%f  ;  t2%user2 = t1%user2  ;  t2%user3 = t1%user3
   t3%f = t1%f  ;  t3%user2 = t1%user2  ;  t3%user3 = t1%user3

   ! Calculate inverse FFT
   call FFFTW_rev(f1,t1%trace)
   call FFFTW_rev(f2,t2%trace)

   ! Replace component names
   t1%kcmpnm = kcmpnm1   ;  t2%kcmpnm = kcmpnm2

end subroutine sw_fdsplitN
!------------------------------------------------------------------------------

!===============================================================================
subroutine sw_fdrotate2d(f1,f2,n,theta,degrees)
!===============================================================================
!  Rotate a pair of frequency domain traces by theta degrees.  It is assumed that
!  f2 has an azimuth 90 degrees clockwise from f1, e.g., f1 // N, f2 // E.
!  The reference frame is rotated clockwise by theta, hence particle motion is
!  rotated anticlockwise.
!  INPUT:
!     f1,f2     : Complex, frequency-domain traces
!     n         : Length of traces
!     theta     : Rotation angle
!     degrees   : Input angle in degrees if .true., radians if .false.
!                 (Default is .true.)

   implicit none
   complex, intent(inout) :: f1(n),f2(n)
   integer, intent(in) :: n
   real, intent(in) :: theta
   logical, optional, intent(in) :: degrees
   real :: conv, R(2,2), real_rot(2,1), imag_rot(2,1)
   integer :: i

   conv = pi/180._rs ! Default to input in degrees
   if (present(degrees)) then
      if (.not.degrees) conv = 1._rs
   endif

   ! Make rotation matrix
   R(1,:) = (/  cos(theta*conv), -sin(theta*conv) /)
   R(2,:) = (/  sin(theta*conv),  cos(theta*conv) /)

   ! Rotate real and imaginary parts separately
   do i=1,n
      real_rot = matmul(R, reshape((/  real(f1(i)),  real(f2(i)) /), (/2,1/)))
      imag_rot = matmul(R, reshape((/ aimag(f1(i)), aimag(f2(i)) /), (/2,1/)))
      f1(i) = cmplx(real_rot(1,1), imag_rot(1,1))
      f2(i) = cmplx(real_rot(2,1), imag_rot(2,1))
   enddo

end subroutine sw_fdrotate2d
!-------------------------------------------------------------------------------

!===============================================================================
subroutine sw_fdtshift(f,n,delta,tshift)
!===============================================================================
!  Shift a frequency domain trace in time (phase in frequency domain).
!  Positive time shifts move the trace later in time (in the time domain).
!  INPUT:
!     f(n)   : complex, frequency domain trace
!     n      : length of trace
!     delta  : sampling interval / s
!     tshift : amount to shift f2 back by / s

   implicit none
   integer, parameter :: ri = r4
   complex(ri), intent(inout) :: f(n)
   integer, intent(in) :: n
   real(ri), intent(in) :: delta, tshift
   integer :: i
   complex(ri), parameter :: j = (0., -1.)
   complex(ri) :: shift, cdelta, ctshift, cpi

   cpi = cmplx(pi,0., kind=ri)
   cdelta = cmplx(delta,0., kind=ri)
   ctshift = cmplx(tshift,0., kind=ri)
   do i=1,n
      shift = (cmplx(2._ri)*cpi*ctshift/cdelta)*cmplx(i-1,kind=ri)/cmplx(2*(n-1),kind=ri)
      f(i) = f(i) * exp(shift*j)
   enddo

end subroutine sw_fdtshift
!-------------------------------------------------------------------------------

!===============================================================================
function sw_eig2x2(a) result(e)
!===============================================================================
!  Compute the eigenvalues for a 2-by-2 matrix.
   implicit none
   real, intent(in) :: a(2,2)
   real :: e(2)
   real :: trace, det
   real, parameter :: tol = 1.e-7
   ! Check for singularity
   trace = a(1,1) + a(2,2)
   det = a(1,1)*a(2,2) - a(1,2)*a(2,1)
   if (abs(det) < tol) call sw_warning('sw_eig2x2: Matrix is singular')
   e(1) = trace/2. + sqrt(trace**2/4. - det)
   e(2) = trace/2. - sqrt(trace**2/4. - det)
end function sw_eig2x2
!-------------------------------------------------------------------------------

!===============================================================================
subroutine sw_error(str)
!===============================================================================
   implicit none
   character(len=*), intent(in) :: str
   write(lu_stderr,'(a)') 'splitwave: Error: '//trim(str)
   stop
end subroutine sw_error
!-------------------------------------------------------------------------------

!===============================================================================
subroutine sw_warning(str)
!===============================================================================
   implicit none
   character(len=*), intent(in) :: str
   write(lu_stderr,'(a)') 'splitwave: Warning: '//trim(str)
end subroutine sw_warning
!-------------------------------------------------------------------------------

end module splitwave
!-------------------------------------------------------------------------------
