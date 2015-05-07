!===============================================================================
! One of Andy Nowacki's Fortran utility modules for dealing with seismic
! anisotropy and other problems.
!
! Andy Nowacki <andy.nowacki@bristol.ac.uk>
!
! See the file LICENCE for licence details.
!===============================================================================
module global_1d_models
!-------------------------------------------------------------------------------
!  Provides some global 1d isotropic velocity models and a few functions relating
!  to 1d models

      implicit none

!  ** size constants
      integer, parameter, private :: i4 = selected_int_kind(9) ; ! long int
      integer, parameter, private :: r4 = selected_real_kind(6,37) ; ! SP
      integer, parameter, private :: r8 = selected_real_kind(15,307) ; ! DP

!  ** precision selector
      integer, parameter, private :: rs = r8

!  ** maths constants and other useful things
      real(rs), parameter, private :: pi = 3.141592653589793238462643 ;
      real(rs), parameter, private :: to_rad = 1.74532925199433e-002 ;  
      real(rs), parameter, private :: to_deg = 57.2957795130823e0 ;  
      real(rs), parameter, private :: to_km = 111.194926644559 ;      

      real(rs), parameter, private :: big_number = 10.e36 ;     

      real(rs), parameter, private :: bigG = 6.67428e-11  ! Newton's gravitational constant

!  ** Available models:
      character(len=32),parameter,private :: available_models(2) = (/ 'AK135', &
                                                                      'PREM ' /)
   contains

!===============================================================================
   function gravity(depth,model)
!===============================================================================
!  Returns the acceleration due to gravity at a given depth in the specified model.
!  Assumes AK135 unless stated.
!  Input depth in km, converted to m for calculation within the function.
!  g = GM/r**2
!  M = 4π∫ rho(r) r**2 dr

     real(rs),intent(in) :: depth
     character(len=*),optional :: model
     character(len=5) :: model_in
     real(rs) :: gravity,r,r0,r1,M,rho,z
     integer :: ir
     real(rs),parameter :: dr=1.

!  Check if we've specified PREM instead, and stop for bad input
     model_in = 'AK135'
     if (present(model)) model_in = trim(model)
     if (trim(model_in) /= 'PREM' .and. trim(model_in) /= 'AK135') then
       write(0,'(a,a,a)') "gravity: model '",trim(model_in),"' not recognised."
       write(0,'(a)') '   Available models are: ', available_models
       stop
     endif
!  Check our depth isn't silly
     if (depth < 0. .or. depth > 6371.) then
       write(0,'(a)') 'gravity: Depth must be in range 0--6371 km'
       stop
     endif

!  Get radius required
     r = 6371. - depth

     M = 0.
     r0 = 0.
!  Calculate mass up to radius required
     do ir=1,int(r),int(dr)
       r1 = real(ir)                ! Radius to top of shell
       z = 6371. - (r0 + r1)/2.     ! Depth to middle of shell
       if (trim(model_in) == 'AK135') call ak135(z,rho=rho)
       if (trim(model_in) == 'PREM' ) then
          call prem(depth, g=gravity)
          return
       endif
       rho = rho * 1.e3             ! Convert to kg/m^3
       M = M + rho*(4./3.)*pi*((r1 * 1.e3)**3 - (r0 * 1.e3)**3) ! Add mass within shell
       r0 = r1                      ! Next layer upwards
     enddo

!  Calculate g
     gravity = bigG *  M / (r * 1.e3)**2

   end function gravity
!-------------------------------------------------------------------------------

!===============================================================================
   function pressure(depth,model)
!===============================================================================
!  Calculate pressure as a function of depth for a specified model.
!  Assumes AK135 unless otherwise stated.
!  Input depth in km.

      real(rs),intent(in) :: depth
      real(rs)            :: pressure
      character(len=*),optional :: model
      character(len=5)          :: model_in
      real(rs)                  :: rho, r, P, r0, r1, dep
      real(rs),parameter        :: dr=10.   ! integration width

      model_in='AK135'
      if (present(model)) model_in = trim(model)
     if (trim(model_in) /= 'PREM' .and. trim(model_in) /= 'AK135') then
       write(0,'(a,a,a)') "pressure: model '",trim(model_in),"' not recognised."
       write(0,'(a)') '   Available models are: ', available_models
       stop
     endif
!  Check our depth isn't silly
     if (depth < 0. .or. depth > 6371.) then
       write(0,'(a)') 'pressure: Depth must be in range 0--6371 km'
       stop
     endif

!  Get radius required
     r = 6371. - depth

!  Calculate pressure by integrating downwards
      P = 0.
      r1 = 6371.
      do while (r1 >= r)
         dep = 6371. - r1
         if (trim(model_in) == 'AK135') call ak135(dep,rho=rho)
         if (trim(model_in) == 'PREM' ) call prem(dep,rho=rho)
         P = P + gravity(dep,model=model)*(rho*1.e3)*(dr*1.e3)
         r1 = r1 - dr
      enddo

      pressure = P / 1.e9   ! GPa

   end function pressure
!-------------------------------------------------------------------------------

!===============================================================================
subroutine ak135(depth,vp,vs,rho)
!===============================================================================
!  AK135 model

   real(rs),intent(in) :: depth
   real(rs),intent(out),optional :: vp,vs,rho
   integer, parameter :: nlayers = 136
   real(rs) :: den1,den2,vp1,vp2,vs1,vs2,dep1,dep2
   integer :: i
   logical :: twodepths

   real(rs), parameter, dimension(nlayers) ::  dep = (/ &
       0.0,   20.0,   20.0,   35.0,   35.0, &
      77.5,  120.0,  165.0,  210.0,  210.0, &
      260.0,  310.0,  360.0,  410.0,  410.0, &
      460.0,  510.0,  560.0,  610.0,  660.0, &
      660.0,  710.0,  760.0,  809.5,  859.0, &
      908.5,  958.0, 1007.5, 1057.0, 1106.5, &
     1156.0, 1205.5, 1255.0, 1304.5, 1354.0, &
     1403.5, 1453.0, 1502.5, 1552.0, 1601.5, &
     1651.0, 1700.5, 1750.0, 1799.5, 1849.0, &
     1898.5, 1948.0, 1997.5, 2047.0, 2096.5, &
     2146.0, 2195.5, 2245.0, 2294.5, 2344.0, &
     2393.5, 2443.0, 2492.5, 2542.0, 2591.5, &
     2640.0, 2690.0, 2740.0, 2740.0, 2789.7, &
     2839.3, 2891.5, 2891.5, 2939.3, 2989.7, &
     3040.0, 3090.3, 3140.7, 3191.0, 3241.3, &
     3291.7, 3342.0, 3392.3, 3442.6, 3493.0, &
     3543.3, 3593.6, 3644.0, 3694.3, 3744.6, &
     3795.0, 3845.3, 3895.6, 3945.9, 3996.3, &
     4046.6, 4096.9, 4147.3, 4197.6, 4247.9, &
     4298.3, 4348.6, 4398.9, 4449.3, 4499.6, &
     4549.9, 4600.3, 4650.6, 4700.9, 4801.6, &
     4851.9, 4902.2, 4952.6, 5002.9, 5053.2, &
     5103.6, 5153.5, 5153.5, 5204.6, 5255.3, &
     5306.0, 5356.8, 5407.5, 5458.2, 5508.9, &
     5559.6, 5610.3, 5661.0, 5711.7, 5813.2, &
     5863.9, 5914.6, 5965.3, 6016.0, 6066.7, &
     6117.4, 6168.1, 6218.9, 6269.6, 6320.3, &
     6371.0 /)

   real(rs), parameter, dimension(nlayers) :: den = (/ &
      2.44900,  2.44900,  2.71420,  2.71420,  3.29760, &
      3.29940,  3.30130,  3.34870,  3.39600,  3.39600, &
      3.46520,  3.53430,  3.60340,  3.67260,  3.79760, &
      3.86120,  3.92480,  3.98850,  4.05210,  4.11580, &
      4.33930,  4.38960,  4.43990,  4.47010,  4.50290, &
      4.53500,  4.56640,  4.59700,  4.62700,  4.65630, &
      4.68490,  4.71300,  4.74040,  4.76730,  4.79360, &
      4.81950,  4.84480,  4.87000,  4.89420,  4.91820, &
      4.94230,  4.96530,  4.98820,  5.01090,  5.03330, &
      5.05530,  5.07720,  5.09900,  5.12060,  5.14220, &
      5.16380,  5.18480,  5.20610,  5.22700,  5.24810, &
      5.26980,  5.29070,  5.31220,  5.33380,  5.35600, &
      5.37760,  5.39990,  5.42240,  5.42240,  5.42380, &
      5.42510,  5.42650,  9.91450,  9.99420, 10.07220, &
     10.14850, 10.22330, 10.29640, 10.36790, 10.43780, &
     10.50620, 10.57310, 10.63850, 10.70230, 10.76470, &
     10.82570, 10.88520, 10.94340, 11.00010, 11.05550, &
     11.10950, 11.16230, 11.21370, 11.26390, 11.31270, &
     11.36040, 11.40690, 11.45210, 11.49620, 11.53910, &
     11.58090, 11.62160, 11.66120, 11.69980, 11.73730, &
     11.77370, 11.80920, 11.84370, 11.87720, 11.94140, &
     11.97220, 12.00010, 12.03110, 12.05930, 12.08670, &
     12.11330, 12.13910, 12.70370, 12.72890, 12.75300, &
     12.77600, 12.79800, 12.81880, 12.83870, 12.85740, &
     12.87510, 12.89170, 12.90720, 12.92170, 12.94740, &
     12.95860, 12.96880, 12.97790, 12.98590, 12.99290, &
     12.99880, 13.00360, 13.00740, 13.01000, 13.01170, &
     13.01220 /)

   real(rs), parameter, dimension(nlayers) ::  vpm = (/ &
     5.8000,  5.8000,  6.5000,  6.5000,  8.0400,  8.0450,  8.0500, &
     8.1750,  8.3000,  8.3000,  8.4825,  8.6650,  8.8475,  9.0300, &
     9.3600,  9.5280,  9.6960,  9.8640, 10.0320, 10.2000, 10.7900, &
    10.9229, 11.0558, 11.1353, 11.2221, 11.3068, 11.3896, 11.4705, &
    11.5495, 11.6269, 11.7026, 11.7766, 11.8491, 11.9200, 11.9895, &
    12.0577, 12.1245, 12.1912, 12.2550, 12.3185, 12.3819, 12.4426, &
    12.5031, 12.5631, 12.6221, 12.6804, 12.7382, 12.7956, 12.8526, &
    12.9096, 12.9668, 13.0222, 13.0783, 13.1336, 13.1894, 13.2465, &
    13.3018, 13.3585, 13.4156, 13.4741, 13.5312, 13.5900, 13.6494, &
    13.6494, 13.6530, 13.6566, 13.6602,  8.0000,  8.0382,  8.1283, &
     8.2213,  8.3122,  8.4001,  8.4861,  8.5692,  8.6496,  8.7283, &
     8.8036,  8.8761,  8.9461,  9.0138,  9.0792,  9.1426,  9.2042, &
     9.2634,  9.3205,  9.3760,  9.4297,  9.4814,  9.5306,  9.5777, &
     9.6232,  9.6673,  9.7100,  9.7513,  9.7914,  9.8304,  9.8682, &
     9.9051,  9.9410,  9.9761, 10.0103, 10.0439, 10.0768, 10.1415, &
    10.1739, 10.2049, 10.2329, 10.2565, 10.2745, 10.2854, 10.2890, &
    11.0427, 11.0585, 11.0718, 11.0850, 11.0983, 11.1166, 11.1316, &
    11.1457, 11.1590, 11.1715, 11.1832, 11.1941, 11.2134, 11.2219, &
    11.2295, 11.2364, 11.2424, 11.2477, 11.2521, 11.2557, 11.2586, &
    11.2606, 11.2618, 11.2622 /)

   real(rs), parameter, dimension(nlayers) :: vsm = (/ &
     3.4600,  3.4600,  3.8500,  3.8500,  4.4800,  4.4900,  4.5000, &
     4.5090,  4.5180,  4.5230,  4.6090,  4.6960,  4.7830,  4.8700, &
     5.0800,  5.1860,  5.2920,  5.3980,  5.5040,  5.6100,  5.9600, &
     6.0897,  6.2095,  6.2426,  6.2798,  6.3160,  6.3512,  6.3854, &
     6.4187,  6.4510,  6.4828,  6.5138,  6.5439,  6.5727,  6.6008, &
     6.6285,  6.6555,  6.6815,  6.7073,  6.7326,  6.7573,  6.7815, &
     6.8052,  6.8286,  6.8515,  6.8742,  6.8972,  6.9194,  6.9418, &
     6.9627,  6.9855,  7.0063,  7.0281,  7.0500,  7.0720,  7.0931, &
     7.1144,  7.1369,  7.1586,  7.1807,  7.2031,  7.2258,  7.2490, &
     7.2490,  7.2597,  7.2704,  7.2811,  0.0000,  0.0000,  0.0000, &
     0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000, &
     0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000, &
     0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000, &
     0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000, &
     0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000, &
     0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000, &
     3.5043,  3.5187,  3.5314,  3.5435,  3.5551,  3.5661,  3.5765, &
     3.5864,  3.5957,  3.6044,  3.6126,  3.6202,  3.6337,  3.6396, &
     3.6450,  3.6498,  3.6540,  3.6577,  3.6608,  3.6633,  3.6653, &
     3.6667,  3.6675,  3.6678 /)

!  Ensure depth is in km and not too large
   if (depth > 6371.0) then
      write(0,'(a,f9.1,a)') '  ak135: depth ',depth,' is too large.'
      stop
   else if (depth < 0.) then
      write(0,'(a)') '  ak135: Depth cannot be negative.'
      stop
   endif

!  Look for the specific point first
   do i=1,136
      if (depth == dep(i)) then 
         if (present(rho)) rho = den(i)
         if (present(vp)) vp = vpm(i)
         if (present(vs)) vs = vsm(i)
         twodepths = .false.
         exit
      else
!  Otherwise, find the points above and below this depth
         if (depth > dep(i) .and. depth < dep(i+1)) then
            dep1 = dep(i);    dep2 = dep(i+1)
            if (present(rho)) then
               den1 = den(i); den2 = den(i+1)
            endif
            if (present(vp)) then
               vp1 = vpm(i);  vp2 = vpm(i+1)
            endif
            if (present(vs)) then
               vs1 = vsm(i);  vs2 = vsm(i+1)
            endif
            twodepths = .true.
            exit
         endif
      endif
   enddo

   if (twodepths) then
!  Extrapolate linearly to find density
      if (present(rho)) rho = den1 + (depth-dep1)*(den2-den1)/(dep2-dep1)
      if (present(vp))   vp = vp1  + (depth-dep1)*(vp2-vp1) / (dep2-dep1)
      if (present(vs))   vs = vs1  + (depth-dep1)*(vs2-vs1) / (dep2-dep1)
   endif

   return
   end subroutine ak135
!-------------------------------------------------------------------------------

!===============================================================================
subroutine prem_old(depth,vp,vs,rho)
!===============================================================================
!  AK135 model

   real(rs),intent(in) :: depth
   real(rs),intent(out),optional :: vp,vs,rho
   integer, parameter :: nlayers = 88
   real(rs) :: den1,den2,vp1,vp2,vs1,vs2,dep1,dep2
   integer :: i
   logical :: twodepths

   real(rs), parameter, dimension(nlayers) :: dep = (/ &
         0.0,    15.0,    15.0,    24.4,    24.4,    40.0,    60.0, &
        80.0,   115.0,   150.0,   185.0,   220.0,   220.0,   265.0, &
       310.0,   355.0,   400.0,   400.0,   450.0,   500.0,   550.0, &
       600.0,   635.0,   670.0,   670.0,   721.0,   771.0,   871.0, &
       971.0,  1071.0,  1171.0,  1271.0,  1371.0,  1471.0,  1571.0, &
      1671.0,  1771.0,  1871.0,  1971.0,  2071.0,  2171.0,  2271.0, &
      2371.0,  2471.0,  2571.0,  2671.0,  2741.0,  2771.0,  2871.0, &
      2891.0,  2891.0,  2971.0,  3071.0,  3171.0,  3271.0,  3371.0, &
      3471.0,  3571.0,  3671.0,  3771.0,  3871.0,  3971.0,  4017.0, &
      4171.0,  4271.0,  4371.0,  4471.0,  4571.0,  4671.0,  4771.0, &
      4871.0,  4971.0,  5071.0,  5149.5,  5149.5,  5171.0,  5271.0, &
      5371.0,  5471.0,  5571.0,  5671.0,  5771.0,  5871.0,  5971.0, &
      6071.0,  6171.0,  6271.0,  6371.0 /)

   real(rs), parameter, dimension(nlayers) :: den = (/ &
      2.6000,  2.6000,  2.9000,  2.9000,  3.3808,  3.3791,  3.3769, &
      3.3747,  3.3709,  3.3671,  3.3633,  3.3595,  3.4358,  3.4626, &
      3.4895,  3.5164,  3.5433,  3.7238,  3.7868,  3.8498,  3.9128, &
      3.9758,  3.9840,  3.9921,  4.3807,  4.4124,  4.4432,  4.5037, &
      4.5631,  4.6213,  4.6784,  4.7346,  4.7898,  4.8442,  4.8978, &
      4.9507,  5.0030,  5.0547,  5.1059,  5.1567,  5.2071,  5.2573, &
      5.3072,  5.3571,  5.4068,  5.4566,  5.4915,  5.5064,  5.5564, &
      5.5664,  9.9035, 10.0294, 10.1813, 10.3273, 10.4673, 10.6015, &
      10.7301, 10.8532, 10.9709, 11.0833, 11.1907, 11.2930, 11.3904, &
      11.4831, 11.5712, 11.6548, 11.7340, 11.8090, 11.8799, 11.9468, &
      12.0099, 12.0692, 12.1250, 12.1663, 12.7636, 12.7749, 12.8250, &
      12.8707, 12.9121, 12.9491, 12.9818, 13.0101, 13.0340, 13.0536, &
      13.0689, 13.0798, 13.0863, 13.0885 /)

   real(rs), parameter, dimension(nlayers) :: vpm = (/ &
      5.8000,  5.8000,  6.8000,  6.8000,  8.1106,  8.1012,  8.0891, &
      8.0769,  8.0554,  8.0337,  8.0118,  7.9897,  8.5590,  8.6455, &
      8.7321,  8.8187,  8.9052,  9.1340,  9.3899,  9.6459,  9.9018, &
      10.1578, 10.2120, 10.2662, 10.7513, 10.9101, 11.0656, 11.2449, &
      11.4156, 11.5783, 11.7336, 11.8821, 12.0244, 12.1613, 12.2932, &
      12.4207, 12.5447, 12.6655, 12.7839, 12.9004, 13.0158, 13.1305, &
      13.2453, 13.3607, 13.4774, 13.5960, 13.6804, 13.6875, 13.7117, &
      13.7166,  8.0648,  8.1994,  8.3602,  8.5130,  8.6580,  8.7957, &
      8.9263,  9.0502,  9.1675,  9.2787,  9.3842,  9.4841,  9.5788, &
      9.6686,  9.7539,  9.8350,  9.9121,  9.9855, 10.0557, 10.1229, &
      10.1874, 10.2496, 10.3097, 10.3557, 11.0283, 11.0364, 11.0725, &
      11.1054, 11.1352, 11.1619, 11.1854, 11.2058, 11.2230, 11.2371, &
      11.2481, 11.2559, 11.2606, 11.2622 /)

   real(rs), parameter, dimension(nlayers) :: vsm = (/ &
      3.2000,  3.2000,  3.9000,  3.9000,  4.4909,  4.4849,  4.4771, &
      4.4695,  4.4564,  4.4436,  4.4311,  4.4188,  4.6439,  4.6754, &
      4.7069,  4.7384,  4.7699,  4.9326,  5.0784,  5.2243,  5.3701, &
      5.5160,  5.5431,  5.5702,  5.9451,  6.0942,  6.2405,  6.3109, &
      6.3781,  6.4423,  6.5037,  6.5625,  6.6189,  6.6732,  6.7255, &
      6.7761,  6.8251,  6.8729,  6.9196,  6.9654,  7.0105,  7.0553, &
      7.0997,  7.1442,  7.1889,  7.2340,  7.2660,  7.2657,  7.2649, &
      7.2647,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000, &
      0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000, &
      0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000, &
      0.0000,  0.0000,  0.0000,  0.0000,  3.5043,  3.5100,  3.5352, &
      3.5582,  3.5791,  3.5977,  3.6141,  3.6284,  3.6404,  3.6503, &
      3.6579,  3.6634,  3.6667,  3.6678 /)

!  Ensure depth is in km and not too large
   if (depth > 6371.0) then
      write(0,'(a,f9.1,a)') '  PREM: depth ',depth,' is too large.'
      stop
   else if (depth < 0.) then
      write(0,'(a)') '  PREM: Depth cannot be negative.'
      stop
   endif

!  Look for the specific point first
   do i=1,88
      if (depth == dep(i)) then 
         if (present(rho)) rho = den(i)
         if (present(vp)) vp = vpm(i)
         if (present(vs)) vs = vsm(i)
         twodepths = .false.
         exit
      else
!  Otherwise, find the points above and below this depth
         if (depth > dep(i) .and. depth < dep(i+1)) then
            dep1 = dep(i);    dep2 = dep(i+1)
            if (present(rho)) then
               den1 = den(i); den2 = den(i+1)
            endif
            if (present(vp)) then
               vp1 = vpm(i);  vp2 = vpm(i+1)
            endif
            if (present(vs)) then
               vs1 = vsm(i);  vs2 = vsm(i+1)
            endif
            twodepths = .true.
            exit
         endif
      endif
   enddo

   if (twodepths) then
!  Extrapolate linearly to find density
      if (present(rho)) rho = den1 + (depth-dep1)*(den2-den1)/(dep2-dep1)
      if (present(vp))   vp = vp1  + (depth-dep1)*(vp2-vp1) / (dep2-dep1)
      if (present(vs))   vs = vs1  + (depth-dep1)*(vs2-vs1) / (dep2-dep1)
   endif

   return
   end subroutine prem_old
!-------------------------------------------------------------------------------

!===============================================================================
subroutine prem(depth,vp,vs,rho,Qmu,Qkappa,vpv,vph,vsv,vsh,eta,g)
!===============================================================================
! PREM model
! Dziewonski & Anderson, 1981, Preliminary reference Earth model, PEPI, 25,
! 297--356.
! Parameterised in terms of normalised radius x = r/a, a = 6371 km
! Output vp and vs are isotropic equivalent velocities; use v{p,s}{v,h} to get
! the true model anisotropic output.
   integer, parameter :: nlayers = 13
   real(rs), intent(in) :: depth
   real(rs), intent(out), optional :: vp, vs, rho, Qmu, Qkappa, vpv, vph, vsv, &
                                      vsh, eta, g
   real(rs), dimension(nlayers), parameter :: &
      r    = (/1221.5, 3480.0, 3630.0, 5600.0, 5701.0, 5771.0, 5971.0, 6151.0, &
               6291.0, 6346.0, 6356.0, 6368.0, 6371.0/), &
      rho0 = (/13.0885, 12.5815,  7.9565,  7.9565,  7.9565,  5.3197, 11.2494, &
                7.1089,  2.6910,  2.6910,  2.9000,  2.6000,  1.0200/), &
      rho1 = (/ 0.0000, -1.2638, -6.4761, -6.4761, -6.4761, -1.4836, -8.0298, &
               -3.8045,  0.6924,  0.6924,  0.0000,  0.0000,  0.0000/), &
      rho2 = (/-8.8381, -3.6426,  5.5283,  5.5283,  5.5283,  0.0000,  0.0000, &
                0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000/), &
      rho3 = (/ 0.0000, -5.5281, -3.0807, -3.0807, -3.0807,  0.0000,  0.0000, &
                0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000/), &
      vp0  = (/11.2622, 11.0487, 15.3891, 24.9520, 29.2766, 19.0957, 39.7027, &
               20.3926,  4.1875,  4.1875,  6.8000,  5.8000,  1.4500/), &
      vp1 =  (/ 0.0000, -4.0362, -5.3181,-40.4673,-23.6027, -9.8672,-32.6166, &
              -12.2569,  3.9382,  3.9382,  0.0000,  0.0000,  0.0000/), &
      vp2 =  (/-6.3640,  4.8023,  5.5242, 51.4832,  5.5242,  0.0000,  0.0000, &
                0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000/), &
      vp3 =  (/ 0.0000,-13.5732, -2.5514,-26.6419, -2.5514,  0.0000,  0.0000, &
                0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000/), &
      vs0 =  (/ 3.6678,  0.0000,  6.9254, 11.1671, 22.3459,  9.9839, 22.3512, &
                8.9496,  2.1518,  2.1519,  3.9000,  3.2000,  0.0000/), &
      vs1 =  (/ 0.0000,  0.0000,  1.4672,-13.7818,-17.2473, -4.9324,-18.5956, &
               -4.4597,  2.3481,  2.3481,  0.0000,  0.0000,  0.0000/), &
      vs2 =  (/-4.4475,  0.0000, -2.0834, 17.4575, -2.0834,  0.0000,  0.0000, &
                0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000/), &
      vs3 =  (/ 0.0000,  0.0000,  0.9783, -9.2777,  0.9783,  0.0000,  0.0000, &
                0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000/), &
      vpv0 = (/11.2622, 11.0487, 15.3891, 24.9520, 29.2766, 19.0957, 39.7027, &
               20.3926,  0.8317,  0.8317,  6.8000,  5.8000,  1.4500/), &
      vpv1 = (/ 0.0000, -4.0362, -5.3181,-40.4673,-23.6027, -9.8672,-32.6166, &
              -12.2569,  7.2180,  7.2180,  0.0000,  0.0000,  0.0000/), &
      vph0 = (/11.2622, 11.0487, 15.3891, 24.9520, 29.2766, 19.0957, 39.7027, &
               20.3926,  3.5908,  3.5908,  6.8000,  5.8000,  1.4500/), &
      vph1 = (/ 0.0000, -4.0362, -5.3181,-40.4673,-23.6027, -9.8672,-32.6166, &
              -12.2569,  4.6172,  4.6172,  0.0000,  0.0000,  0.0000/), &
      vsv0 = (/ 3.6678,  0.0000,  6.9254, 11.1671, 22.3459,  9.9839, 22.3512, &
                8.9496,  5.8582,  5.8582,  3.9000,  3.2000,  0.0000/), &
      vsv1 = (/ 0.0000,  0.0000,  1.4672,-13.7818,-17.2473, -4.9324,-18.5956, &
               -4.4597, -1.4678, -1.4678,  0.0000,  0.0000,  0.0000/), &
      vsh0 = (/ 3.6678,  0.0000,  6.9254, 11.1671, 22.3459,  9.9839, 22.3512, &
                8.9496, -1.0839, -1.0839,  3.9000,  3.2000,  0.0000/), &
      vsh1 = (/ 0.0000,  0.0000,  1.4672,-13.7818,-17.2473, -4.9324,-18.5956, &
               -4.4597,  5.7176,  5.7176,  0.0000,  0.0000,  0.0000/), &
      Qmu_model =  (/84.6, -1., 312., 312., 312., 143., 143., 143., 80., 600., &
                     600., 600., -1./), &
      Qkappa_model = (/1327.7, 57823., 57823., 57823., 57823., 57823., 57823., &
                       57823., 57823., 57823., 57823., 57823., 57823./), &
      eta0 = (/1., 1., 1., 1., 1., 1., 1., 1., 3.3687, 3.3687, 1., 1., 1./), &
      eta1 = (/1., 1., 1., 1., 1., 1., 1., 1.,-2.4778,-2.4778, 1., 1., 1./)
   real(rs), parameter :: a = r(nlayers)
   real(rs), parameter, dimension(nlayers) :: rho0SI = rho0*1.e3, rho1SI = rho1/a, &
                          rho2SI = rho2*1.e-3/a**2, rho3SI = rho3*1.e-6/a**3
   real(rs) :: x, x0, M, r0, r1
   integer :: i, j

   if (depth > a) then
      write(0,'(a,f9.1,a)') '  PREM: depth ',depth,' is too large.'
      stop
   else if (depth < 0.) then
      write(0,'(a)') '  PREM: Depth cannot be negative.'
      stop
   endif

   x = (a - depth)/a
   i = 1
   do while(i < nlayers)
      if (x < r(i)/a) exit
      i = i + 1
   enddo

   if (present(rho)) rho = rho0(i) + rho1(i)*x + rho2(i)*x**2 + rho3(i)*x**3
   if (present(vp)) vp = vp0(i) + vp1(i)*x + vp2(i)*x**2 + vp3(i)*x**3
   if (present(vs)) vs = vs0(i) + vs1(i)*x + vs2(i)*x**2 + vs3(i)*x**3
   if (present(Qmu)) Qmu = Qmu_model(i)
   if (present(Qkappa)) Qkappa = Qkappa_model(i)
   ! Note for v{p,s}{v,h}, differences are only in x^0 and x^1, so we intentionally
   ! use v{p,s}{2,3}(i)
   if (present(vpv)) vpv = vpv0(i) + vpv1(i)*x + vp2(i)*x**2 + vp3(i)*x**3
   if (present(vph)) vph = vph0(i) + vph1(i)*x + vp2(i)*x**2 + vp3(i)*x**3
   if (present(vsv)) vsv = vsv0(i) + vsv1(i)*x + vs2(i)*x**2 + vs3(i)*x**3
   if (present(vsh)) vsh = vsh0(i) + vsh1(i)*x + vs2(i)*x**2 + vs3(i)*x**3
   if (present(eta)) eta = eta0(i) + eta1(i)*x

   ! Calculate gravity by integrating the density upwards
   if (present(g)) then
      M = 0._rs
      do j = 1, i
         if (j == 1) then
            r0 = 0._rs
         else
            r0 = r(j-1)*1.e3
         endif
         if (j < i) then
            r1 = r(j)*1.e3
         else
            r1 = (a - depth)*1.e3
         endif
         M = M + rho0SI(j)*(r1**3 - r0**3)/3._rs + rho1SI(j)*(r1**4 - r0**4)/4._rs + &
                 rho2SI(j)*(r1**5 - r0**5)/5._rs + rho3SI(j)*(r1**6 - r0**6)/6._rs
      enddo
      g = bigG*(4._rs*pi*M)/r1**2
   endif

end subroutine prem
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
end module global_1d_models