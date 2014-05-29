!===============================================================================
!  Module containing spherical geometry helper functions and subroutines.
!  Andy Nowacki, University of Bristol
!  andy.nowacki@bristol.ac.uk
!
!  History:
!      2011-04-12:  Added sphere_sample subroutine to return an array of points
!                   which evenly sample a sphere.
!      2011-07-18:  Added routines to find Earth radial direction
!      2011-11-08:  Added sph_poly_inout: determines if point is inside or outside
!                   a set of points (ordered) on a sphere.
!
!===============================================================================
module spherical_geometry

!  ** size constants
   integer, parameter, private :: i4 = selected_int_kind(9)       ! long int
   integer, parameter, private :: r4 = selected_real_kind(6,37)   ! SP
   integer, parameter, private :: r8 = selected_real_kind(15,307) ! DP

!  ** precision selector
   integer, parameter, private :: rs = r8

!  ** maths constants and other useful things
   real(rs), parameter, private :: pi = 3.141592653589793238462643_rs
   real(rs), parameter, private :: pi2 = pi/2._rs
   real(rs), parameter, private :: twopi = 2._rs*pi
   real(rs), parameter, private :: to_rad = 1.74532925199433e-002
   real(rs), parameter, private :: to_deg = 57.2957795130823e0
   real(rs), parameter, private :: to_km = 111.194926644559

   real(rs), parameter, private :: big_number = 10.e36
   real(rs), parameter, private :: angle_tol = 1.e-5

   ! Random number sampler state
   logical, save, private :: rng_initialised = .false.

   contains

!------------------------------------------------------------------------------
   function delta(lon1_in,lat1_in,lon2_in,lat2_in,degrees)
!------------------------------------------------------------------------------
!  delta returns the angular distance between two points on a sphere given
!  the lat and lon of each using the Haversine formula
!
     implicit none
     real(rs) :: delta,lat1,lon1,lat2,lon2
     real(rs),intent(in) :: lat1_in,lon1_in,lat2_in,lon2_in
     logical,optional :: degrees

     lat1 = lat1_in;  lon1 = lon1_in
     lat2 = lat2_in;  lon2 = lon2_in
     if (present(degrees)) then
        if (degrees) then
         lat1=lat1*pi/1.8D2 ; lon1=lon1*pi/1.8D2
         lat2=lat2*pi/1.8D2 ; lon2=lon2*pi/1.8D2
        endif
     endif

     delta=atan2( sqrt( (cos(lat2)*sin(lon2-lon1))**2 + (cos(lat1)*sin(lat2) - &
           sin(lat1)*cos(lat2)*cos(lon2-lon1))**2) , &
           sin(lat1)*sin(lat2) + cos(lat1)*cos(lat2)*cos(lon2-lon1))

     if (present(degrees)) then
         if (degrees) delta = delta * 1.8D2/pi
      endif

     return

   end function delta
!==============================================================================

!------------------------------------------------------------------------------
!   function dist_vincenty(lon1_in,lat1_in,lon2_in,lat2_in,R,&
!                           a_in,b_in,azi,baz,degrees)
   subroutine dist_vincenty(lon1_in,lat1_in,lon2_in,lat2_in,dist,R,a_in,b_in,degrees)
!------------------------------------------------------------------------------
!  delta_vincenty uses the Vincenty algorithm to find accurate great circle
!  distances on a flattened sphere.  Input in radians unless specified.
!  Unless R is specified, then distances are on the surface of the Earth.
!  R is fractional radius of points, given the ellipsoid of the Earth.
!  Also computes the azimuth and backazimuth, if asked for
!  Fom: http://www.movable-type.co.uk/scripts/latlong-vincenty.html

      implicit none

      real(rs),intent(in)  :: lon1_in,lat1_in,lon2_in,lat2_in
      real(rs),intent(in),optional :: a_in,b_in,R
      real(rs)              :: conversion
      logical,intent(in),optional :: degrees
      real(rs),intent(out)  :: dist
      real(rs)             :: a,b,f,lon1,lat1,lon2,lat2,L,u1,u2,lambda1,lambda2,&
                              sin_s,cos_s,s,sin_a,cos_a_sq,cos_2sm,C,u_2,&
                              big_A,big_B,ds
      logical              :: isnan
!      real,intent(out),optional :: azi,baz

      write(*,'(a)') 'subroutine dist_vincenty is not working yet.'
      stop

!  Convert from degrees if necessary
      conversion = 1.d0
      if (present(degrees)) then
         if (degrees) conversion = pi/1.8d2
      endif

!  Get the inputs
      lon1 = conversion * (lon1_in)  ;  lat1 = conversion * (lat1_in)
      lon2 = conversion * (lon2_in)  ;  lat2 = conversion * (lat2_in)

!  Default ellipsoidal shape is WGS-84: override these if both a and b are supplied
      a = 6378.137d0  ;  b = 6356.752314245d0
      if (present(a_in).and.present(b_in)) then
         a = a_in  ;  b = b_in
      endif
!  Scale by the fractional radius of the ellipsoid if present
      if (present(R)) then
         a = a*R  ; b = b*R
      endif
      f = (a-b)/a

!  Constants
      L = lon2 - lon1
      u1 = atan((1-f)*tan(lat1))
      u2 = atan((1-f)*tan(lat2))
!  Starting guess of L for lambda
      lambda1 = 0.d0
      lambda2 = L
!  Iterate until convergence
      do while (abs(lambda2-lambda1) < 1.d-12)
         sin_s = sqrt( (cos(u2)*sin(lambda2))**2 + (cos(u1)*sin(u2) - &
                      sin(u1)*cos(u2)*cos(lambda2))**2 )
         if (sin_s==0.d0) then
            dist = 0.d0
            return
         endif
         cos_s = sin(u1)*sin(u2) + cos(u1)*cos(u2)*cos(lambda2)
         s = atan2(sin_s,cos_s)
         sin_a = cos(u1)*cos(u2)*sin(lambda2)/sin_s
         cos_a_sq = 1.d0 - sin_a**2
         cos_2sm = cos_s - 2.d0*sin(u1)*sin(u2)/(cos_a_sq)
         if (isnan(cos_2sm)) cos_2sm = 0.d0
         C = (f/16.d0)*(cos_a_sq)*(4.d0+f*(4.d0-3.d0*cos_a_sq))
         lambda1 = lambda2
         lambda2 = L+(1.d0-C)*f*sin_a*(s+C*sin_s*(cos_2sm+C*cos_s*(-1.d0+2.d0*cos_2sm**2)))
      enddo

      u_2 = cos_a_sq*(a**2-b**2)/(b**2)
      big_A = 1.d0 + (u_2/16384.d0)*(4096.d0+u_2*(-768.d0+u_2*(320.d0-175.d0*u_2)))
      big_B = (u_2/1024.d0)*(256.d0+u_2*(-128.d0+u_2*(74.d0-47.d0*u_2)))
      ds = big_B*sin_s*(cos_2sm+(big_B/4.d0)*(cos_s*(-1.d0+2.d0*cos_2sm**2)-&
           (big_B/6.d0)*cos_2sm*(-3.d0+4.d0*sin_s**2)*(-3.d0+4.d0*cos_2sm**2)))
      dist = b*big_A*(s-ds)
!      if (present(azi)) &
!         azi = atan2(cos(u2)*sin(lambda2), &
!                     cos(u1)*sin(u2)-sin(u1)*cos(u2)*cos(lambda2))
!      if (present(baz)) &
!         baz = atan2(cos(u1)*sin(lambda2), &
!                     -sin(u1)*cos(u2)+cos(u1)*sin(u2)*cos(lambda2))

      write(*,*) 'lambda =',lambda2

   end subroutine dist_vincenty
!==============================================================================

!------------------------------------------------------------------------------
!   function test_dist_vincenty(lon1_in,lat1_in,lon2_in,lat2_in,a_in,b_in,R,degrees)
!------------------------------------------------------------------------------
!   implicit none
!
!   real(rs),intent(in)    :: lon1_in,lon2_in,lat1_in,lat2_in
!   real(rs),intent(in),optional :: a_in,b_in,R
!   real(rs)               :: lon1,lon2,lat1,lat2
!   real(rs)               :: test_dist_vincenty
!   logical,intent(in),optional :: degrees
!   real(rs)               :: a,b,f,L,u1,u2,lambda1,lambda2,sin_s,cos_s,s,sin_a,&
!                             cos_a_sq,cos_2sm,C,u_2,big_A,big_B,Delta_s,dist,&
!                             conversion
!   real(rs),parameter     :: convergence_limit = 1.d-12  !  To within ~6mm
!
!
!   conversion = 1.d0
!   if (present(degrees)) then
!      if (degrees) conversion = pi/1.8d2
!   endif
!
!!  Get the inputs
!   lon1 = conversion * (lon1_in)  ;  lat1 = conversion * (lat1_in)
!   lon2 = conversion * (lon2_in)  ;  lat2 = conversion * (lat2_in)
!
!!  Default ellipsoidal shape is WGS-84: override these if both a and b are supplied
!   a = 6378.137d0 ; b = 6356.752314245d0
!   if (present(a_in).and.present(b_in)) then
!      a = a_in  ;  b = b_in
!   endif
!
!!  Scale by the fractional radius of the ellipsoid if present
!   if (present(R)) then
!      a = a*R  ; b = b*R
!   endif
!
!!  Constants
!   f = (a-b)/a ;
!   L = lon2 - lon1 ;
!   u1 = atan2((1-f)*tan(lat1),1.d0)
!   u2 = atan2((1-f)*tan(lat2),1.d0)
!
!!  Starting guess for the iteration variables
!   lambda1 = big_number
!   lambda2 = L
!
!   do while ( abs(lambda2-lambda1) > convergence_limit )
!     sin_s = sqrt((cos(u2)*sin(lambda2))**2 + (cos(u1)*sin(u2)-sin(u1)*cos(u2)*cos(lambda2))**2)
!     cos_s = sin(u1)*sin(u2) + cos(u1)*cos(u2)*cos(lambda2)
!     s = atan2(sin_s,cos_s)
!     sin_a = cos(u1)*cos(u2)*sin(lambda2)/sin_s
!     cos_a_sq = 1.d0 - sin_a**2
!     cos_2sm = cos_s - s*sin(u1)*sin(u2)/cos_a_sq
!      if (isnan(cos_2sm)) cos_2sm = 0.d0
!     C = (f/16.d0)*cos_a_sq*(4.d0+f*(4.d0-3.d0*cos_a_sq))
!     lambda1 = lambda2
!     lambda2 = L + (1.d0-C)*f*sin_a*(s+C*sin_a*(cos_2sm+C*cos_s*(-1.d0+2.d0*cos_2sm**2)))
!   enddo
!
!   u_2 = cos_a_sq*(a**2-b**2)/b**2
!   big_A = 1.d0+(u_2/16384.d0)*(4096.d0+u_2*(-768.d0+u_2*(320.d0-175.d0*u_2)))
!   big_B = (u_2/1024.d0)*(256.d0+u_2*(-128.d0+u_2*(74.d0-47.d0*u_2)))
!   Delta_s = big_B*sin_s*(cos_2sm+(big_B/4.d0)*(cos_s*(-1.d0+2.d0*cos_2sm**2)- &
!          (big_B/6.d0)*cos_2sm*(-3.d0+4.d0*sin_s**2)*(-3.d0+4.d0*cos_2sm**2)))
!   test_dist_vincenty = b*big_A*(s-Delta_s) ;
!
!
!   return
!
!   end function test_dist_vincenty
!==============================================================================

!------------------------------------------------------------------------------
   subroutine step(lon1_in,lat1_in,az_in,delta_in,lon2,lat2,degrees)
!------------------------------------------------------------------------------
! Computes the endpoint given a starting point lon,lat, azimuth and angular distance

     implicit none

     real(rs),intent(in)  :: lon1_in,lat1_in,az_in,delta_in
     real(rs),intent(out) :: lon2,lat2
     real(rs)             :: lon1,lat1,az,delta
     logical,optional,intent(in) :: degrees
     logical          :: using_deg,deg

     lon1=lon1_in ; lat1=lat1_in ; az=az_in ; delta=delta_in

     using_deg = present(degrees)
     if (using_deg) then
        if (degrees) deg = degrees
     else
        deg = .false.
     endif

     if (using_deg .and. deg) then
       if (delta > 180.) then
          write(*,*)'spherical_geometry: step: Error: distance must be less than 180 degrees.'
          stop
       else if (lon1 <-180 .or. lon1 > 180) then
          write(*,*)'spherical_geometry: step: Error: longitude must be in range -180 - 180.'
          stop
       else if (lat1 <-90 .or. lat1 > 90) then
          write(*,*)'spherical_geometry: step: Error: latitude must be in range -90 - 90.'
          stop
       endif
     else
       if (delta > pi) then
            write(*,*)'spherical_geometry: step: Error: distance must be less than 2pi radians.'
            stop
        else if (lon1 < -pi .or. lon1 > pi) then
            write(*,*)'spherical_geometry: step: Error: longitude must be in range -2pi - 2pi.'
            stop
        else if (lat1 < -pi/2.d0 .or. lat1 > pi/2.d0) then
            write(*,*)'spherical_geometry: step: Error: latitude must be in range -pi - pi.'
          stop
        endif
     endif

     if (using_deg .and. deg) then
!  Convert to radians
        lon1=lon1*pi/1.8D2 ; lat1=lat1*pi/1.8D2
        az=az*pi/1.8D2     ; delta=delta*pi/1.8D2
     endif

!  Calculate point which is delta degrees/radians from lon1,lat1 along az
     lat2 = asin(sin(lat1)*cos(delta) + cos(lat1)*sin(delta)*cos(az))
     lon2 = lon1 + atan2(sin(az)*sin(delta)*cos(lat1),  &
                         cos(delta)-sin(lat1)*sin(lat2) )

     if (using_deg .and. deg) then
!  Convert to degrees
        lat2=1.8D2*lat2/pi  ; lon2=1.8D2*lon2/pi
        if(lon2>1.8D2) lon2=lon2-3.6D2 ; if(lon2<-1.8D2) lon2=lon2+3.6D2
     end if

     return

   end subroutine step
!==============================================================================

!-------------------------------------------------------------------------------
   subroutine gcp_points(lon1,lat1,lon2,lat2,ptslon,ptslat,npts,ds,n,degrees)
!-------------------------------------------------------------------------------
!  Returns arrays of lon and lat points (including the end points)
!  along a great circle path between two endpoints.  The user must specify
!  one of the separation distance, ds (degrees or radians), or number of points
!  (including the end points), n.  geographic coordinates can be in degrees or
!  radians, and the points are returned in the same format.  The array must be
!  at least npts long (optionally returned by the subroutine).

      implicit none
      real(rs),intent(in) :: lon1,lat1,lon2,lat2
      real(rs) :: x1,x2,y1,y2,d,ddelta,azi
      real(rs),dimension(:),intent(out) :: ptslon,ptslat
      integer,intent(out),optional :: npts
      real(rs),intent(in),optional :: ds
      integer,intent(in),optional :: n
      logical,intent(in),optional :: degrees
      real(rs) :: conversion
      integer :: i,npoints

      ! Default to radians; convert from degrees if necessary
      conversion = 1._rs
      if (present(degrees)) then
         if (degrees) conversion = pi/180._rs
      endif
      x1 = lon1*conversion;  x2 = lon2*conversion
      y1 = lat1*conversion;  y2 = lat2*conversion

      ! Check one of ds or n is present
      if (.not.present(ds) .and. .not.present(n)) then
         write(*,'(a)') 'spherical_geometry: gcp_points: Error: one of ds or n must be specified'
         stop
      endif

      ! Get distance and azimuth between two points
      d = delta(x1,y1,x2,y2,degrees=.false.)
      azi = azimuth(x1,y1,x2,y2,degrees=.false.)

      ! Using a fixed distance ds (same units as coordinates)
      if (present(ds)) then
         if (.not.present(npts)) then
            write(*,'(a)') 'spherical_geometry: gcp_points: Error: must supply npts as well as ds if using constant spacing'
            stop
         endif
         ddelta = ds*conversion
         ! Check we haven't asked for too much
         if (ddelta > d) then
            write(*,'(a)') 'spherical_geometry: gcp_points: Error: requested point spacing is larger than distance between points'
            stop
         endif
         ! Calculate number of points
         npoints = ceiling(d/ddelta) + 1
         npts = npoints

      ! Using a fixed number of points
      else
         npoints = n
         ! Have supplied too few points
         if (npoints < 2) then
            write(*,'(a)') 'spherical_geometry: gcp_points: Error: n must be at least 2'
            stop
         ! Have asked for two points--which are the end points we've supplied!
         else if (npoints == 2) then
            write(*,'(a)') &
            'spherical_geometry: gcp_points: Warning: have asked for only two points on path, so returning end points'
            ptslon(2) = x2
            ptslat(2) = y2
            ptslon(1:2) = ptslon(1:2) / conversion
            ptslat(1:2) = ptslat(1:2) / conversion
            return
         endif
         ! Calculate distance between points
         ddelta = d/npoints
      endif

      ! Check we have room for all the points
      if (size(ptslon) < npoints .or. size(ptslat) < npoints) then
         write(*,'(a)') 'spherical_geometry: gcp_points: Error: arrays to hold lon and lat points are not long enough'
         stop
      endif

      ! Fill in points
      ptslon(1) = x1
      ptslat(1) = y1
      do i=2,npoints-1
         call step(x1,y1,azi,real(i-1)*ddelta,ptslon(i),ptslat(i),degrees=.false.)
      enddo
      ptslon(npoints) = x2
      ptslat(npoints) = y2

      ! Convert back
      ptslon(1:npoints) = ptslon(1:npoints) / conversion
      ptslat(1:npoints) = ptslat(1:npoints) / conversion


   end subroutine gcp_points
!===============================================================================

!------------------------------------------------------------------------------
   function azimuth(lon1,lat1,lon2,lat2,degrees)
!  Returns azimuth from point 1 to point 2.
!  From: http://www.movable-type.co.uk/scripts/latlong.html

     implicit none

     real(rs) :: azimuth,lon1,lat1,lon2,lat2
     real(rs) :: rlon1,rlat1,rlon2,rlat2,conversion
     logical,optional :: degrees

     conversion = 1._rs
     if (present(degrees)) then
        if (degrees) conversion = pi/180._rs
     endif

     rlon1 = conversion*lon1  ;  rlon2 = conversion*lon2
     rlat1 = conversion*lat1  ;  rlat2 = conversion*lat2

     azimuth = atan2(sin(rlon2-rlon1)*cos(rlat2) , &
               cos(rlat1)*sin(rlat2) - sin(rlat1)*cos(rlat2)*cos(rlon2-rlon1) )

     if (azimuth < 0) then
        azimuth = azimuth+2._rs*pi
     endif

     azimuth = azimuth / conversion

!     write(*,*)'Azimuth',azimuth

     return

   end function azimuth
!==============================================================================

!------------------------------------------------------------------------------
   subroutine geog2cart(phi_in,theta_in,r,x,y,z,degrees)
!  Returns the cartesian coordinates from geographical ones
!  Theta is latitude, phi is longitude and r is radius

      implicit none

      real(rs),intent(in)  :: theta_in,phi_in,r
      real(rs),intent(out) :: x,y,z
      real(rs)             :: theta,phi,conversion
      logical,optional,intent(in)  :: degrees

      conversion = 1._rs
      if (present(degrees)) then
         if (degrees) conversion = pi/180._rs
      endif
      theta = theta_in * conversion
      phi = phi_in * conversion

      if (theta < -pi/2._rs .or. theta > pi/2._rs) then
         write(*,'(a)') 'Latitude must be in range -pi/2 - pi/2 (-90 - 90 deg).'
         stop
      endif

      x = r * sin(pi/2._rs - theta) * cos(phi)
      y = r * sin(pi/2._rs - theta) * sin(phi)
      z = r * cos(pi/2._rs - theta)

      return

   end subroutine geog2cart
!==============================================================================

!------------------------------------------------------------------------------
   subroutine sph2cart(phi_in,theta_in,r,x,y,z,degrees)
!  Returns the cartesian coordinates from spherical ones
!  Theta is colatitude, phi is longitude and r is radius

      implicit none

      real(rs),intent(in)  :: theta_in,phi_in,r
      real(rs),intent(out) :: x,y,z
      real(rs)             :: theta,phi
      logical,optional,intent(in)  :: degrees

      if (present(degrees)) then
         if (degrees) then
            theta = theta_in * pi / 1.8d2
            phi = phi_in * pi / 1.8d2
         else
            theta = theta_in
            phi = phi_in
         endif
      else
         theta = theta_in
         phi = phi_in
      endif

       if (theta < 0.d0 .or. theta > pi ) then
          write(*,'(a)') 'Colatitude must be in range 0--pi (0--180deg).'
          stop
       endif

       x = r * sin(theta) * cos(phi)
       y = r * sin(theta) * sin(phi)
       z = r * cos(theta)

       return

   end subroutine sph2cart
!==============================================================================

!------------------------------------------------------------------------------
   subroutine cart2geog(x,y,z,theta,phi,r,degrees)
!  Returns the geographic coordinates from cartesian ones.

   implicit none

   real(rs),intent(in)  :: x,y,z
   real(rs),intent(out) :: theta,phi,r
   real(rs)             :: t,p,r_temp
   logical,optional     :: degrees

   r_temp = sqrt(x**2 + y**2 + z**2)

   t = acos(z/r_temp)
   p = atan2(y,x)

   r = r_temp

   if (present(degrees)) then
      if (degrees) then
         theta = 90.d0 - t * 1.8d2/pi
         phi   = p   * 1.8d2/pi
      else
         theta = pi/2.d0 - t
         phi = p
      endif
   else
      theta = pi/2.d0 - t
      phi = p
   endif

   return

   end subroutine cart2geog
!==============================================================================

!------------------------------------------------------------------------------
   subroutine cart2sph(x,y,z,theta,phi,r,degrees)
!  Returns the spherical coordinates from cartesian ones.

   implicit none

   real(rs),intent(in)  :: x,y,z
   real(rs),intent(out) :: theta,phi,r
   real(rs)             :: t,p,r_temp
   logical,optional     :: degrees

   r_temp = sqrt(x**2 + y**2 + z**2)

   t = acos(z/r_temp)
   p = atan2(y,x)

   r = r_temp

   if (present(degrees)) then
      if (degrees) then
         theta = t * 1.8d2/pi
         phi   = p   * 1.8d2/pi
      else
         theta =  t
         phi = p
      endif
   else
      theta = t
      phi = p
   endif

   return

   end subroutine cart2sph
!==============================================================================

!------------------------------------------------------------------------------
   function inclination(r_in,lon_in,lat_in,degrees)
!  Give the inclination of a vector in cartesian coordinates, given
!  the latitude and longitude.
!  Inclination is measured away from the Earth radial direction, hence
!  0 for an upward, vertical ray, 90° for a horizontal ray, 180° for a downward,
!  vertical ray

     implicit none

     real(rs),intent(in)   :: r_in(3),lon_in,lat_in
     real(rs)              :: inclination
     real(rs)              :: r(3),lon,lat,radial(3),conversion,dot
     logical,intent(in),optional :: degrees

!  Convert to radians if necessary
     if (present(degrees)) then
       if (degrees) conversion = pi/1.8d2
       if (.not.degrees) conversion = 1.d0
     else
       conversion = 1.d0
     endif
     lon = conversion * lon_in ; lat = conversion * lat_in

!  Create the (unit) cartesian vector along the Earth radial direction
     radial(1) = cos(lat)*cos(lon)
     radial(2) = cos(lat)*sin(lon)
     radial(3) = sin(lat)

!  Make r into unit vector
     r = r_in / sqrt(r_in(1)**2 + r_in(2)**2 + r_in(3)**2)

!  Compute the dot product and the inclination
     dot = r(1)*radial(1) + r(2)*radial(2) * r(3)*radial(3)
     inclination = abs(acos(dot))

     if (inclination > pi/2.d0) inclination = pi - inclination

     inclination = inclination / conversion
!      write(*,*)'Inclination',inclination

     return

   end function inclination
!==============================================================================

!-------------------------------------------------------------------------------
   function xyz2radial(x,y,z)
!  Given Cartesian coordinates of convention
!     1 goes through (0,0)
!     2 goes through (90E,0)
!     3 goes through N pole,
!  produce the Earth radial direction in Cartesian coordinates

      implicit none
      real(rs),intent(in) :: x,y,z
      real(rs)            :: xyz2radial(3), r

      r = sqrt(x**2 + y**2 + z**2)
      xyz2radial(1) = x / r
      xyz2radial(2) = y / r
      xyz2radial(3) = z / r

      return
   end function xyz2radial
!===============================================================================

!-------------------------------------------------------------------------------
   function lonlat2radial(lon,lat,degrees)
!  Given a longitude and latitude, give the Earth radial direction in the standard
!  Cartesian reference system (see e.g. xyz2radial)
!  Default is input in radians: override with degrees=.true.
      implicit none
      real(rs),intent(in) :: lon,lat
      real(rs)            :: lonlat2radial(3),x,y,z,r
      logical,optional,intent(in) :: degrees
      logical                     :: degrees_in

!  Check for input in degrees and pass on as appropriate
      degrees_in = .false.
      if (present(degrees)) degrees_in = degrees

      r = 1000._rs ! Dummy radius
      call geog2cart(lon, lat, r, x, y, z, degrees=degrees_in)
      lonlat2radial = xyz2radial(x,y,z)

      return
   end function lonlat2radial
!===============================================================================

!-------------------------------------------------------------------------------
   subroutine sg_sphere_sample(d,lon_out,lat_out,n_out)
!  Evenly sample a sphere given an input distance d between adjacent points.
!  Points are in longitude range -180 to 180.
!  lon and lat are column vectors which are assigned within the subroutine.

      implicit none

      real(rs),intent(in)  :: d
      real(rs),allocatable,intent(out) :: lon_out(:), lat_out(:)
      integer,intent(out)  :: n_out
      integer,parameter    :: nmax=50000
      real(rs)             :: lon(nmax),lat(nmax)
      real(rs)             :: dlon,dlat,dlon_i,lon_i,lat_i
      integer              :: i,n,n_i

      n = 1

      dlat = d
      dlon = dlat  ! At the equator

      lat(n)=90.; lon(n)=0.

      lat_i = lat(1) - dlat
      do while (lat_i > -90.)
         dlon_i = dlon/sin((90.-lat_i)*pi/180.)
         n_i = nint(360./dlon_i)
         do i=1,n_i
            n = n + 1
            if (n > nmax) then
               write(0,'(a)') 'sphere_sample: number of points greater than nmax',&
                              '  Change compiled limits or increase point spacing d.'
               stop
            endif
            lat(n) = lat_i
            lon(n) = lon_i
            lon_i = modulo(lon_i + dlon_i, 360.)
         enddo
         lon_i = modulo(lon_i + dlon_i, 360.)
         lat_i = lat_i - dlat
      enddo

      n = n + 1
      lat(n)=-90. ; lon(n) = 0.

      if (allocated(lon_out)) then
         if (size(lon_out) /= n) then
            deallocate(lon_out)
            allocate(lon_out(n))
         endif
      else
         allocate(lon_out(n))
      endif

      if (allocated(lat_out)) then
         if (size(lat_out) /= n) then
            deallocate(lat_out)
            allocate(lat_out(n))
         endif
      else
         allocate(lat_out(n))
      endif

      lon_out(1:n) = mod(lon(1:n) + 180., 360.) - 180.
      lat_out(1:n) = lat(1:n)
      n_out = n

      return

    end subroutine sg_sphere_sample
!===============================================================================

!-------------------------------------------------------------------------------
   subroutine sphere_sample(d,lon_out,lat_out,n_out)
! Deprecated synonym for sg_sphere_sample: we provide this wrapper subroutine to
! for backwards compatibility.
      implicit none
      real(rs), intent(in) :: d
      real(rs), intent(out), allocatable, dimension(:) :: lon_out, lat_out
      integer, intent(out) :: n_out
      call sg_sphere_sample(d,lon_out,lat_out,n_out)
   end subroutine sphere_sample
!===============================================================================

!-------------------------------------------------------------------------------
   function sph_poly_inout(x,y,px,py,degrees)
!  Takes in assumed-shape arrays (vectors) for points on a sphere, which must be
!  ordered either clockwise or anticlockwise.  The function assumes that the first
!  and last points are not the same, but this doesn't matter anyway.
!  x,y:         trial point in lon,lat
!  px(:),py(:): polygon vertices in lon,lat
!
!  NOTE: This algorithm won't work for sample points on the north or south poles,
!        because the azimuths will always be 0 or 180, and hence the total will
!        always be zero.  This can be alleviated by implementing the algorithm
!        described in:
!                  Schettino (1999). Polygon intersections in spherical
!        topology: applications to plate tectonics.  Computers & Geosciences, 25
!        (1) 61-69. doi:10.1016/S0098-3004(98)00081-8

      implicit none

      real(rs),intent(in) :: x,y
      real(rs),intent(in),dimension(:) :: px,py
      logical,intent(in),optional :: degrees
      logical :: sph_poly_inout
      logical :: deg
      real(rs) :: conversion, s, a0, a1, da, tx, ty, tpx0, tpy0, tpx1, tpy1
      integer :: i,n
      real(rs),parameter :: tol = 1._rs  ! Tolerance in *degrees*

!  Check for same size arrays
      if (size(px) /= size(py)) then
         write(0,'(a)') &
'spherical_geometry: sph_poly_inout: Error: polygon coordinate vectors must be same length.'
         stop
      endif

!  Check for degrees/radians
      deg = .false.
      conversion = 1._rs
      if (present(degrees)) then
         deg = degrees
         if (degrees) conversion = pi/180._rs
      endif

      tx = conversion*x
      ty = conversion*y

!  Check for point on vertex
      if (any(x == px .and. y == py)) then
         write(0,'(a)') 'spherical_geometry: sph_poly_inout: point is on vertex.'
         stop
      endif

!  Check for point on poles
      if (y == 90. .or. y == -90.) then
         write(0,'(a)') 'spherical_geometry: sph_poly_inout: point is on one of the poles.'
         stop
      endif

!  Loop over sides and calculate sum of angles.  If ~360, inside.  If ~0, outside
      n = size(px)
      s = 0.
      do i = 1,n-1
         tpx0 = conversion*px(i)
         tpy0 = conversion*py(i)
         tpx1 = conversion*px(i+1)
         tpy1 = conversion*py(i+1)
         a0 = azimuth(tx,ty,tpx0,tpy0,degrees=.true.)
         a1 = azimuth(tx,ty,tpx1,tpy1,degrees=.true.)
         da = a1 - a0
         do while (da > 180._rs)
            da = da - 360._rs
         enddo
         do while (da < -180._rs)
            da = da + 360._rs
         enddo
         s = s + da
      enddo
!  Calculate difference between last and first.  da == 0 if given a closed set of points.
      tpx0 = conversion*px(n)
      tpy0 = conversion*py(n)
      tpx1 = conversion*px(1)
      tpy1 = conversion*py(1)
      a0 = azimuth(tx,ty,tpx0,tpy0,degrees=.true.)
      a1 = azimuth(tx,ty,tpx1,tpy1,degrees=.true.)
      da = a1 - a0
      do while (da > 180._rs)
         da = da - 360._rs
      enddo
      do while (da < -180._rs)
         da = da + 360._rs
      enddo
      s = s + da
!      write(*,*) s

!  Test for in or out
      if (360._rs - abs(s) <= tol) then
         sph_poly_inout = .true.
      else
         sph_poly_inout = .false.
      endif

      return
   end function sph_poly_inout
!===============================================================================

!===============================================================================
function sg_torad(a)
!===============================================================================
!  Convert from degrees to radians
   implicit none
   real(rs), intent(in) :: a
   real(rs) :: sg_torad
   sg_torad = a*to_rad
end function sg_torad
!-------------------------------------------------------------------------------

!===============================================================================
function sg_todeg(a)
!===============================================================================
!  Convert from radians to degrees
   implicit none
   real(rs), intent(in) :: a
   real(rs) :: sg_todeg
   sg_todeg = a*to_deg
end function sg_todeg
!-------------------------------------------------------------------------------

!===============================================================================
subroutine sg_project_to_gcp(long,latg,lonp,latp,lon,lat,degrees)
!===============================================================================
!  Using the pole to a great circle path about a sphere with coordinates long,
!  latg, project a point onto that path (point at lonp,latp)
!  Input is assumed to be in radians.
   implicit none
   real(rs), intent(in) :: long,latg,lonp,latp
   real(rs), intent(out) :: lon,lat
   logical, optional, intent(in) :: degrees
   real(rs), dimension(3) :: g,p,gp,pp
   real(rs) :: r
   logical :: deg

   deg = .false.
   if (present(degrees)) deg = degrees
   ! Convert to vectors
   g = sg_lonlat2vec(long,latg,degrees=deg)
   p = sg_lonlat2vec(lonp,latp,degrees=deg)

   ! Check whether the point already lies on the gcp
   if (abs(dot_product(g,p)) < angle_tol) then
      call cart2geog(p(1),p(2),p(3),lat,lon,r,degrees=deg)
      return
   endif

   ! Compute the pole to the gcp containing g and the point
   gp = sg_cross_prod(g,p)

   ! Check for the point and the pole being the same--we can't handle this
   if (sqrt(gp(1)**2 + gp(2)**2 + gp(3)**2) < angle_tol) then
      write(0,'(a)') 'spherical_geometry: sg_project_to_gcp: Error: point and pole to plane are the same'
      stop
   endif

   gp = gp/sqrt(sum(gp**2))  ! Normalise to unit vector
   pp = sg_cross_prod(gp,g)
   ! Have to swap the sign of the point if the angle between the starting and
   ! projected point is more than 90 deg.  Everything is a unit vector.
   if (acos(dot_product(p,pp)) > pi2) pp = -pp

   call cart2geog(pp(1),pp(2),pp(3),lat,lon,r,degrees=deg)

end subroutine sg_project_to_gcp
!-------------------------------------------------------------------------------

!===============================================================================
function sg_cross_prod(a,b) result(c)
!===============================================================================
!  Compute the cross product for two three-vectors.
   implicit none
   real(rs), dimension(3), intent(in) :: a,b
   real(rs) :: c(3)
   c = (/a(2)*b(3) - a(3)*b(2),  a(3)*b(1) - a(1)*b(3),  a(1)*b(2) - a(2)*b(1)/)
end function sg_cross_prod
!-------------------------------------------------------------------------------

!===============================================================================
function sg_lonlat2vec(lon,lat,degrees) result(v)
!===============================================================================
!  Convert a longitude and latitude on a unit sphere to a vector.
   implicit none
   real(rs), intent(in) :: lon,lat
   logical, optional, intent(in) :: degrees
   real(rs) :: r,x,y,z,v(3)
   logical :: deg
   deg = .false.
   if (present(degrees)) deg = degrees
   r = 1._rs
   call geog2cart(lon,lat,r,x,y,z,degrees=deg)
   v = (/x, y, z/)
end function sg_lonlat2vec
!-------------------------------------------------------------------------------

!===============================================================================
subroutine sg_gcp_from_point_azi(lon,lat,azi,lonp,latp,degrees)
!===============================================================================
!  Given a point on a unit sphere and and azimuth from that point, return the
!  pole to the great cricle path created by that point and azimuth
   implicit none
   real(rs), intent(in) :: lon,lat,azi
   real(rs), intent(out) :: lonp,latp
   logical, optional, intent(in) :: degrees
   logical :: deg
   real(rs), dimension(3) :: pstart, pstep, gp
   real(rs) :: r,conversion
   real(rs), parameter :: step_dist = 45._rs

   deg = .false.
   if (present(degrees)) deg = degrees
   conversion = 1._rs
   if (deg) conversion = to_rad
   ! Convert to a vector
   pstart = sg_lonlat2vec(lon, lat, degrees=deg)
   ! Find a point along the gcp
   call step(lon, lat, azi, step_dist*conversion, lonp, latp, degrees=deg)
   pstep = sg_lonlat2vec(lonp,latp,degrees=deg)
   ! Find the pole to the two points
   gp = sg_cross_prod(pstart,pstep)
   call cart2geog(gp(1),gp(2),gp(3),latp,lonp,r,degrees=deg)

end subroutine sg_gcp_from_point_azi
!-------------------------------------------------------------------------------

!===============================================================================
subroutine sg_gcp_from_points(lon1,lat1,lon2,lat2,lonp,latp,degrees)
!===============================================================================
!  Given two points on a sphere (geographic coordinates), return the coordinates
!  of the pole to the great circle containing the two points.
   implicit none
   real(rs), intent(in) :: lon1,lat1,lon2,lat2
   real(rs), intent(out) :: lonp,latp
   logical, optional, intent(in) :: degrees
   logical :: deg
   real(rs), dimension(3) :: p1,p2,g
   real(rs) :: r,conversion

   deg = .false.
   if (present(degrees)) deg = degrees
   conversion = 1._rs
   if (deg) conversion = to_rad
   ! Check the points don't overlap
   if (conversion*abs(lon1-lon2) < angle_tol .and. conversion*abs(lat1-lat2) < angle_tol) then
      write(0,'(a)') 'spherical_geometry: sg_gcp_from_points: Error: two points overlap'
      stop
   endif
   ! Check the points aren't antipodal
   if (abs(delta(lon1,lat1,lon2,lat2,degrees=deg) - pi/conversion) < angle_tol) then
      write(0,'(a)') 'spherical_geometry: sg_gcp_from_points: Error: two points are antipodal'
      stop
   endif
   ! Convert the coordinates to vectors and compute the cross product, which
   ! gives the pole to the great circle
   p1 = sg_lonlat2vec(lon1,lat1,degrees=deg)
   p2 = sg_lonlat2vec(lon2,lat2,degrees=deg)
   g = sg_cross_prod(p1,p2)
   g = g/sqrt(sum(g**2))
   call cart2geog(g(1),g(2),g(3),latp,lonp,r,degrees=deg)

end subroutine sg_gcp_from_points
!-------------------------------------------------------------------------------

!===============================================================================
function sg_gcp_to_azimuth(long,latg,lonp,latp,degrees) result(azi)
!===============================================================================
!  Calculate the local azimuth of the great circle defined by pole (long,latg)
!  at a point (lonp,latp).  Default is for input/output in radians; use
!  degrees=.true. for azimuth in degrees.
   implicit none
   real(rs), intent(in) :: long,latg,lonp,latp
   real(rs) :: azi
   logical, optional, intent(in) :: degrees
   real(rs) :: lon,lat,conversion

   conversion = 1._rs
   if (present(degrees)) then
      if (degrees) conversion = to_rad
   endif
   ! Check values
   if (latg < -pi2/conversion .or. latg > pi2/conversion .or. &
       latp < -pi2/conversion .or. latp > pi2/conversion) then
      write(0,'(a)') 'spherical_geometry: sg_gcp_to_azimuth: Error: '// &
         'latitude must be in range -pi/2 to pi/2 (-90 to 90 deg)'
      stop
   endif
   ! Project the point onto the great circle
   call sg_project_to_gcp(long,latg,lonp,latp,lon,lat,degrees=degrees)
   ! The local azimuth is pi/2 (90 deg) from the azimuth to the pole, given
   ! in the range 0-pi or 0-360
   azi = modulo(azimuth(lon,lat,long,latg,degrees=degrees) + pi2/conversion + &
         twopi/conversion, twopi/conversion)
end function sg_gcp_to_azimuth
!-------------------------------------------------------------------------------

!===============================================================================
function sg_angle_180(a) result(angle)
!===============================================================================
!  Make an angle be in the range -180 to 180 deg
   implicit none
   real(rs), intent(in) :: a
   real(rs) :: angle
   angle = modulo(a+180._rs,360._rs) - 180._rs
end function sg_angle_180
!-------------------------------------------------------------------------------

!===============================================================================
function sg_angle_360(a) result(angle)
!===============================================================================
!  Make an angle be in the range 0 to 360 deg
   implicit none
   real(rs), intent(in) :: a
   real(rs) :: angle
   angle = modulo(a,360._rs)
end function sg_angle_360
!-------------------------------------------------------------------------------

!===============================================================================
function sg_angle_diff(a,b) result(diff)
!===============================================================================
!  Provide the difference between two angles, handling cases where they are
!  either side of 0/360 or -180/180
   implicit none
   real(rs), intent(in) :: a,b
   real(rs) :: diff
   diff = abs(sg_angle_360(b) - sg_angle_360(a))
   if (diff > 180._rs) diff = 360._rs - diff
end function sg_angle_diff
!-------------------------------------------------------------------------------

!===============================================================================
subroutine sg_initialise_rng(s)
!===============================================================================
! Initialise the random number generator.  This involves seeding the RNG with
! an array of integers of a minimum length.  This is 12 with gfortran 4.6.
! OPTIONAL: Supply an integer to mix in to the seed.  This may be useful when
! calling in parallel to make sure the PRNG state is initialised differently.
   implicit none
   integer, optional, intent(in) :: s
   integer :: values(8), n, i
   integer, dimension(:), allocatable :: seed
   if (.not. rng_initialised) then
      call random_seed(size=n)  ! Find minimum size of array used to seed
      allocate(seed(n))
      call date_and_time(values=values)
      if (present(s)) values(7) = s*values(7)
      seed = values(8) + values(7)*(/ (i-1, i=1,n) /)
      call random_seed(put=seed)
      deallocate(seed)
      rng_initialised = .true.
   endif
end subroutine sg_initialise_rng
!-------------------------------------------------------------------------------

!===============================================================================
subroutine sg_random_points_sph(lon, colat, degrees)
!===============================================================================
! Return a set of points on a sphere randomly distributed, in spherical
! coordinates.
! From: http://mathworld.wolfram.com/SpherePointPicking.html
! (NB: They use (to me) non-standard theta=lon and phi=colat!)
   implicit none
   real(rs), intent(out), dimension(:) :: lon, colat
   logical, optional, intent(in) :: degrees
   logical :: radians
   ! Check array lengths
   if (size(lon) /= size(colat)) then
      write(0,'(a)') 'spherical_geometry: sg_random_points_sph: Error: arrays ' &
         // 'for lon and colat are not the same length'
      stop
   endif
   radians = .true.
   if (present(degrees)) radians = .not.degrees
   call sg_initialise_rng
   call random_number(colat)
   colat = acos(2._rs*colat - 1._rs)
   call random_number(lon)
   lon = twopi*lon

   if (.not.radians) then
      lon = to_deg*lon
      colat = to_deg*colat
   endif
end subroutine sg_random_points_sph
!-------------------------------------------------------------------------------

!===============================================================================
subroutine sg_random_points_geog(lon, lat, degrees)
!===============================================================================
! Return a set of points on a sphere randomly distributed, in geographical
! coordinates.
   implicit none
   real(rs), intent(out), dimension(:) :: lon, lat
   logical, optional :: degrees
   logical :: radians
   radians = .true.
   if (present(degrees)) radians = .not.degrees
   call sg_random_points_sph(lon, lat)
   lat = pi2 - lat
   if (.not.radians) then
      lon = to_deg*lon
      lat = to_deg*lat
   endif
end subroutine sg_random_points_geog
!-------------------------------------------------------------------------------

!===============================================================================
subroutine sg_random_point_sph(lon, colat, degrees)
!===============================================================================
   implicit none
   real(rs), intent(out) :: lon, colat
   logical, optional, intent(in) :: degrees
   logical :: radians
   real(rs), dimension(1) :: alon, acolat
   radians = .true.
   if (present(degrees)) radians = .not.degrees
   call sg_random_points_sph(alon, acolat, degrees=.not.radians)
   lon = alon(1)
   colat = acolat(1)
end subroutine sg_random_point_sph
!-------------------------------------------------------------------------------

!===============================================================================
subroutine sg_random_point_geog(lon, lat, degrees)
!===============================================================================
   implicit none
   real(rs), intent(out) :: lon, lat
   logical, optional, intent(in) :: degrees
   logical :: radians
   real(rs), dimension(1) :: alon, alat
   radians = .true.
   if (present(degrees)) radians = .not.degrees
   call sg_random_points_geog(alon, alat, degrees=.not.radians)
   lon = alon(1)
   lat = alat(1)
end subroutine sg_random_point_geog
!-------------------------------------------------------------------------------

!______________________________________________________________________________
end module spherical_geometry
