!===============================================================================
module moment_tensor
!===============================================================================
! The moment_tensor module provides routines for dealing with moment tensors,
! such as creating them from double-couple parameters (strike, dip, rake, moment)
! or evaluating source radiation patterns from them.
!
! This module uses the Harvard/Global CMT convention for all input/output:
!    x // r     (local radial = upwards)
!    y // theta (local south = colatitude)
!    z // phi   (local east = (co)longitude)
!
! (Individual routines may use different conventions internally.)
!
! MTs should be passed around as 6-element arrays with the following construction:
!    MT = [Mrr, Mtt, Mpp, Mrt, Mrp, Mtp]
! They can then be accessed using the vector index parameters like so:
!    MT(tt) = 1.5.

   implicit none
   
   ! All variables and routines private by default
   private
   
   ! Precision selectors
   integer, parameter :: r4 = selected_real_kind(6, 37), &
                         r8 = selected_real_kind(15, 307)
   integer, parameter :: rs = r8

   ! Constants
   real(rs), parameter :: pi = 4._rs*atan2(1._rs, 1._rs)
   real(rs), parameter :: deg = 180._rs/pi, rad = pi/180._rs
   real(rs), parameter :: ZERO = 0._rs, &
                          ONE = 1._rs, &
                          HALF = 1._rs/2._rs, &
                          TWOPI = 2._rs*pi, &
                          PI2 = pi/2._rs

   ! IO
   integer, parameter :: lu_stderr = 0, &
                         lu_stdin  = 5, &
                         lu_stdout = 6

   ! Conversion from 6-vector to array elements
   integer, parameter :: rr = 1, tt = 2, pp = 3, rt = 4, rp = 5, tp = 6

   public :: mt_radiation_pattern

contains

!===============================================================================
subroutine mt_radiation_pattern(M, azi, inc, P, SV, SH, j)
!===============================================================================
! Calculate the radiation pattern from a general moment tensor (or list of).
! Uses formula given in pp. 70 ff. of
!   The Seismic Wavefield, Volume 1.  Kennett, B.L.N., Cambridge University Press.
! In Kennett's notation:
!   x // N
!   y // E
!   z // down
!   i is angle from z towards x-y plane (incidence angle away from down)
!   phi is angle from x towards y (i.e., azimuth clocwise from N)
!   V is direction upwards (radial) when looking along ray at source
!   H is direction to the right when looking along ray at source:
!
!	^ V (SV)
!   |
!   | j   /       View looking along ray
!   |-   /
!   | \ /
!   |  /
!   | /
!   |/
!   x----------> H (SH)
!
! INPUT:
!   M(6) : Moment tensor
!   azi  : Azimuth from north towards east of ray (deg)
!   inc  : Inclination from downwards towards the radial of ray (deg)
! OUTPUT:
!   P, SV, SH : Amplitudes of ray in P, SV and SH direction (defined above)
!   j    : Source polarisation measured from SV towards SH (i.e., like phi' in splitting)
   implicit none
   real(rs), intent(in) :: M(6), azi, inc
   real(rs), intent(out) :: P, SV, SH, j
   real(rs) :: Mxx, Myy, Mzz, Mxy, Mxz, Myz, a, i, sini, cosi, sinphi, cosphi, &
               sintwophi, costwophi

   ! Convert to Kennett's convention
   Mxx =  M(tt)
   Myy =  M(pp)
   Mzz =  M(rr)
   Mxy = -M(tp)
   Mxz =  M(rt)
   Myz = -M(rp)

   ! Convert to radians
   a = rad*azi
   i = rad*inc

   ! Some shortcuts
	sini = sin(i)
	cosi = cos(i)
	sinphi = sin(a)
	cosphi = cos(a)
	sintwophi = sin(2._rs*a)
	costwophi = cos(2._rs*a)
	! Calculate radiation pattern for P, SV and SH
	P = (sini**2)*(Mxx*cosphi**2 + Mxy*sintwophi + Myy*sinphi**2 - Mzz) &
		+ 2._rs*sini*cosi*(Mxz*cosphi + Myz*sinphi) + Mzz
	SV = sini*cosi*(Mxx*cosphi**2 + Mxy*sintwophi + Myy*sinphi**2 - Mzz) &
		+ cos(2._rs*i)*(Mxz*cosphi + Myz*sinphi)
	SH = sini*((Myy - Mxx)*sinphi*cosphi + Mxy*costwophi) &
		+ cosi*(Myz*cosphi - Mxz*sinphi)
	! Source polarisation in ray frame, measured from upwards towards the right
	j = deg*atan2(SV,SH)

end subroutine mt_radiation_pattern
!-------------------------------------------------------------------------------

end module moment_tensor
!-------------------------------------------------------------------------------
