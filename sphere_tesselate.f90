!===============================================================================
! One of Andy Nowacki's Fortran utility modules for dealing with seismic
! anisotropy and other problems.
!
! Andy Nowacki <andy.nowacki@bristol.ac.uk>
!
! See the file LICENCE for licence details.
!===============================================================================
module sphere_tesselate
!===============================================================================
! Create a tesselation of the unit sphere

! use spherical_geometry

implicit none

private

! Size constants
integer, parameter, private :: i4 = selected_int_kind(9)
integer, parameter, private :: r4 = selected_real_kind(6,37)
integer, parameter, private :: r8 = selected_real_kind(15,307)
! Precision selector
integer, parameter, private :: rs = r8

! Useful constants
real(rs), parameter, private :: pi = 4._rs*atan2(1._rs,1._rs)

! Debugging variables
logical, parameter, private :: debug = .false.

! I/O unit for files
integer, parameter, private :: lu = 10
! File format version descriptor
integer, parameter, private :: st_version = 1

! Tolerance for duplicate point searching; this must be less than the closest
! expected point spacing.  It is the 3D distance between the points
real(rs), parameter, private :: point_tol = 1.e-4_rs

! Location for cache files.  The directory is pre-set to a default location,
! and routines look in the environment variable SEISMO_FORTRAN_DATA for the
! data directory.  Routines can set st_cache_dir to fix where the files go,
! but can't override the environment variable.
character(len=250), parameter, private :: st_cache_env_var = 'SEISMO_FORTRAN_DATA'
character(len=250), public :: st_cache_dir = '$HOME/Applications/modules/data'

! Public types exposed to calling routines
! A point in 3d space
type, public :: point
   real(rs) :: x, y, z
end type point

! A triangle on the unit sphere.  a, b and c are the integer indices of another
! array which contains the coordinates of the vertices of the triangle
type, public :: triangle
   integer :: a, b, c
end type triangle

! A tesselation of a sphere, with np vertices and nt triangles
! p contains all unique points in the tesselation
! t contains the indices of p which correspond to the vertices of each triangle
! b contains the barycentric coordinates of each triangle
! Hence to get the corners of triangle i, you can say:
!     pt = tess%p(tess%t(i))
!     x = pt%x;  y = pt%y;  z = pt%z
type, public :: tesselation
   integer :: level
   type(point), allocatable, dimension(:) :: p
   type(triangle), allocatable, dimension(:) :: t
   integer :: np, nt
   type(point), allocatable, dimension(:) :: b
end type tesselation


! Declare the public-facing routines
interface st_norm_p
   module procedure :: st_norm_p_array, st_norm_p_single
end interface st_norm_p

public :: &
   st_dump_points, &
   st_dump_triangles, &
   st_icosahedron, &
   st_iterate_level, &
   st_level_from_points, &
   st_load_cache, &
   st_load_tesselation, &
   st_new, &
   st_num_faces, &
   st_num_points, &
   st_norm_p, &
   st_save_cache, &
   st_save_tesselation, &
   st_rotate


contains

!===============================================================================
subroutine st_save_cache(t)
!===============================================================================
!  Save a tesselation to the default cache location.  This is read from the
!  environment variable $SEISMO_FORTRAN_DATA, or defaults to st_cache_dir
!  if that is not available.
   type(tesselation), intent(inout) :: t
   character(len=250) :: dir, file
   logical :: dir_exists

   call get_environment_variable(st_cache_env_var, dir)
   if (dir == '') dir = st_cache_dir
   ! Check that the directory exists
   inquire(file=dir, exist=dir_exists)
   if (.not.dir_exists) then
      write(0,'(a)') 'st_save_cache: Error: Directory "' // trim(dir) // &
         '" does not exist'
      error stop
   endif
   write(file,'("ico_pole_",i0.1,".tess")') t%level
   file = trim(dir) // '/' // trim(file)
   call st_save_tesselation(t, file)
end subroutine st_save_cache
!-------------------------------------------------------------------------------

!===============================================================================
subroutine st_load_cache(t, exists)
!===============================================================================
!  Load a tesselation from the default cache location.  The t%level value must
!  be filled with the desired tesselation level.  At the moment, only polar,
!  icosahedral tesselations are supported.
!  Optionally, the exists argument specified whether the cache file exists.
!  If this is supplied, the routine will not exit with an error when the file
!  does not exist.
   type(tesselation), intent(inout) :: t
   logical, optional, intent(out) :: exists
   logical :: file_exists
   character(len=250) :: dir, file

   call get_environment_variable(st_cache_env_var, dir)
   if (dir == '') dir = st_cache_dir
   if (t%level < 0 .or. t%level > 15) then
      write(0,'(a,i0.1)') 'st_load_cache: Error: tesselation level is ', t%level
      error stop
   endif
   write(file,'("ico_pole_",i0.1,".tess")') t%level
   file = trim(dir) // '/' // trim(file)
   inquire(file=file, exist=file_exists)
   if (present(exists)) then
      exists = file_exists
      if (.not.file_exists) return
   else
      if (.not.file_exists) then
         write(0,'(a)') 'st_load_cache: Error: Cache file "' // trim(file) &
            // '" does not exist'
         error stop
      endif
   endif
   call st_load_tesselation(file, t)
end subroutine st_load_cache
!-------------------------------------------------------------------------------

!===============================================================================
subroutine st_save_tesselation(t, file, ascii)
!===============================================================================
!  Save a tesselation to an output file.  This can be ASCII or Fortran binary;
!  use ascii=.true. for the former.  These can then be re-read in later.
!  Binary format is the default
   type(tesselation), intent(inout) :: t
   character(len=*), intent(in) :: file
   logical, optional, intent(in) :: ascii
   logical :: ascii_in
   integer :: i, iostat
   character(len=30) :: string

   ascii_in = .false.
   if (present(ascii)) ascii_in = ascii
   if (.not.allocated(t%b)) call st_compute_barycentres(t)
   ! ASCII file
   if (ascii_in) then
      open(lu, file=file, iostat=iostat)
      call check_open
      write(lu,'(a,i0.1)') '# sphere_tesselate version ', st_version
      write(lu,'(i0.1," ",i0.1," ",i0.1)') t%level, t%np, t%nt ! Check on level and num points
      do i = 1, t%np
         write(lu,*) t%p(i)
      enddo
      do i = 1, t%nt
         write(lu,*) t%t(i)
      enddo
      do i = 1, t%nt
         write(lu,*) t%b(i)
      enddo
      close(lu)

   ! Fortran binary
   else
      open(lu, file=file, form='unformatted', iostat=iostat)
      call check_open
      write(string,'(a,i0.1)') '# sphere_tesselate version ', st_version
      write(lu) string
      write(lu) st_version
      write(lu) t%level
      write(lu) t%np
      write(lu) t%nt
      write(lu) t%p
      write(lu) t%t
      write(lu) t%b
      close(lu)
   endif

   contains
      subroutine check_open
         if (iostat /= 0) then
            write(0,'(a)') 'st_save_tesselation: Error: Cannot open file "', &
               trim(file), '" for reading'
            error stop
         endif
      end subroutine check_open
end subroutine st_save_tesselation
!-------------------------------------------------------------------------------

!===============================================================================
subroutine st_load_tesselation(file, t, ascii)
!===============================================================================
!  Read a previously-saved tesselation.  This can be ASCII or Fortran, the latter
!  being the default.
   character(len=*), intent(in) :: file
   type(tesselation), intent(out) :: t
   logical, intent(in), optional :: ascii
   logical :: ascii_in
   character(len=250) :: dummy1, dummy2, dummy3
   character(len=30) :: version_string
   integer :: i, iostat, version

   ascii_in = .false.
   if (present(ascii)) ascii_in = ascii
   ! ASCII file
   if (ascii_in) then
      open(lu, file=file, iostat=iostat)
      call check_open
      read(lu,*) dummy1, dummy2, dummy3, version
      call check_version
      read(lu,*) t%level, t%np, t%nt
      call check_num_points
      call reallocate_arrays
      do i = 1, t%np
         read(lu,*) t%p(i)
      enddo
      do i = 1, t%nt
         read(lu,*) t%t(i)
      enddo
      do i = 1, t%nt
         read(lu,*) t%b(i)
      enddo
      close(lu)

   ! Fortran binary
   else
      open(lu, file=file, form='unformatted', iostat=iostat)
      call check_open
      read(lu) version_string
      read(lu) version
      call check_version
      read(lu) t%level
      read(lu) t%np
      read(lu) t%nt
      call check_num_points
      call reallocate_arrays
      read(lu) t%p
      read(lu) t%t
      read(lu) t%b
      close(lu)
   endif

   contains
      subroutine reallocate_arrays
         if (allocated(t%p)) then
            if (size(t%p) /= t%np) then
               deallocate(t%p)
               allocate(t%p(t%np))
            endif
         else
            allocate(t%p(t%np))
         endif
         if (allocated(t%t)) then
            if (size(t%t) /= t%nt) then
               deallocate(t%t)
               allocate(t%t(t%nt))
            endif
         else
            allocate(t%t(t%nt))
         endif
         if (allocated(t%b)) then
            if (size(t%b) /= t%nt) then
               deallocate(t%b)
               allocate(t%b(t%nt))
            endif
         else
            allocate(t%b(t%nt))
         endif
      end subroutine reallocate_arrays

      subroutine check_open
         if (iostat /= 0) then
            write(0,'(a)') 'st_load_tesselation: Error: Cannot open file "', &
               trim(file), '" for reading'
            error stop
         endif
      end subroutine check_open

      subroutine check_version
         if (version /= st_version) then
            write(0,'(2(a,i0.1),a)') 'st_load_tesselation: Error: Input file "' // &
               trim(file) // '" has a different version (', version, &
               ') to that expected (', st_version ,')'
            error stop
         endif
      end subroutine check_version

      subroutine check_num_points
         if (t%np /= st_num_points(t%level) .or. t%nt /= st_num_faces(t%level)) then
            write(0,'(5(a,i0.1),a)') 'st_load_tesselation: Error: Input file "' // &
               trim(file) // '" has (np, nt) = (', t%np, ', ', t%nt, '), which is not ' &
               // 'as expected for level ', t%level, '(', t%np, ', ', t%nt, ')'
            error stop
         endif
      end subroutine check_num_points
end subroutine st_load_tesselation
!-------------------------------------------------------------------------------

!===============================================================================
subroutine st_dump_points(t, geog)
!===============================================================================
!  Write the points of a tesselation to stdout.
!  If geog == .true., write out lon, lat of points instead.
   use spherical_geometry, only: cart2geog
   type(tesselation), intent(in) :: t
   logical, intent(in), optional :: geog
   integer :: i
   real(rs) :: lon, lat, r

   do i = 1, t%np
      if (present(geog)) then
         if (geog) then
            call cart2geog(t%p(i)%x, t%p(i)%y, t%p(i)%z, lat, lon, r, &
               degrees=.true.)
            write(*,*) lon, lat
         else
            write(*,*) t%p(i)
         endif
      else
         write(*,*) t%p(i)
      endif
   enddo
end subroutine st_dump_points
!-------------------------------------------------------------------------------

!===============================================================================
subroutine st_dump_triangles(t, geog)
!===============================================================================
!  Write the triangles to stdout, separated by a line beginning '>' (which
!  is suitable for GMT multisegment output).
!  If geog == .true., write out lon, lat instead.
   use spherical_geometry, only: cart2geog
   type(tesselation), intent(in) :: t
   logical, intent(in), optional :: geog
   integer :: i
   real(rs) :: lon, lat, r, lona, lata

   do i = 1, t%nt
      write(*,'(a)') '>'
      if (present(geog)) then
         if (geog) then
            call cart2geog(t%p(t%t(i)%a)%x, t%p(t%t(i)%a)%y, t%p(t%t(i)%a)%z, &
               lata, lona, r, degrees=.true.)
            write(*,*) lona, lata
            call cart2geog(t%p(t%t(i)%b)%x, t%p(t%t(i)%b)%y, t%p(t%t(i)%b)%z, &
               lat, lon, r, degrees=.true.)
            write(*,*) lon, lat
            call cart2geog(t%p(t%t(i)%c)%x, t%p(t%t(i)%c)%y, t%p(t%t(i)%c)%z, &
               lat, lon, r, degrees=.true.)
            write(*,*) lon, lat
            write(*,*) lona, lata
         else
            write(*,*) t%p(t%t(i)%a)
            write(*,*) t%p(t%t(i)%b)
            write(*,*) t%p(t%t(i)%c)
            write(*,*) t%p(t%t(i)%a)
         endif
      else
         write(*,*) t%p(t%t(i)%a)
         write(*,*) t%p(t%t(i)%b)
         write(*,*) t%p(t%t(i)%c)
         write(*,*) t%p(t%t(i)%a)
      endif
   enddo
end subroutine st_dump_triangles
!-------------------------------------------------------------------------------

!===============================================================================
function st_new(level, shape, pole) result(t)
!===============================================================================
!  Create a new tesselation at one's desired level, ranging from 0 upwards.
!  Optionally, supply the name of a starting shape.  Current options are:
!     'icosahedron' ('ico' or 'i') [default]
!  Optionally, ask for a tesselation with points aligned on the poles with the
!  pole=.true. argument.
   implicit none
   integer, intent(in) :: level
   logical, intent(in), optional :: pole
   character(len=*), intent(in), optional :: shape
   character(len=50) :: shape_in
   type(tesselation) :: t
   logical :: pole_in
   integer :: i

   if (level < 0) then
      write(0,'(a)') 'st_new: Error: level must be >= 0'
      error stop
   endif
   pole_in = .false.
   if (present(pole)) pole_in = pole
   shape_in = 'icosahedron'
   if (present(shape)) shape_in = tolower(shape)
   ! Choose starting shape
   select case(shape_in)
      case('i', 'ico', 'icosahedron')
         t = st_icosahedron(pole=pole_in)
      case default
         write(0,'(a)') 'st_new: Error: shape "'//trim(shape_in)//'" not recognised'
         error stop
   end select
   ! Iterate to the desired level
   do i = 1, level
      call st_iterate_level(t)
   enddo
end function st_new
!-------------------------------------------------------------------------------

!===============================================================================
function st_icosahedron(pole) result(t)
!===============================================================================
!  Return the coordinates and triangles of an icosahedron (12-sided).
!  By default, the original icosahedron is arranged so that no point lies on the
!  pole.  This is optimal when looking at data at the poles (see Teanby, C&G, 2006).
!  Pass pole=.true. to instead use an icosahedron where the poles have points;
!  this is the arrangement more common in global simulations (see Sadourny,
!  Arakawa and Mintz, Monthly Weather Review, 1968).
   logical, intent(in), optional :: pole
   logical :: pole_in
   type(tesselation) :: t
   real(rs), parameter :: phi = 2._rs*cos(pi/5._rs), &
      a = 1._rs/sqrt(5._rs), & ! Height of points for polar orientation
      b = 2._rs/sqrt(5._rs), & ! Radius of points for polar orientation
      pi5 = pi/5._rs

   ! Number of points
   t%level = 0
   t%np = 12
   t%nt = 20
   if (allocated(t%p)) deallocate(t%p)
   allocate(t%p(t%np))
   if (allocated(t%t)) deallocate(t%t)
   allocate(t%t(t%nt))

   pole_in = .false.
   if (present(pole)) pole_in = pole

   ! First orientation is the division of the octahedron into segments with the golden
   ! ratio, phi
   if (.not.pole_in) then
      ! Vertices of triangles, unnormalised
      t%p(1)  = point( 0._rs,    phi,  1._rs)
      t%p(2)  = point( 0._rs,   -phi,  1._rs)
      t%p(3)  = point( 0._rs,    phi, -1._rs)
      t%p(4)  = point( 0._rs,   -phi, -1._rs)
      t%p(5)  = point( 1._rs,  0._rs,    phi)
      t%p(6)  = point(-1._rs,  0._rs,    phi)
      t%p(7)  = point( 1._rs,  0._rs,   -phi)
      t%p(8)  = point(-1._rs,  0._rs,   -phi)
      t%p(9)  = point(   phi,  1._rs,  0._rs)
      t%p(10) = point(  -phi,  1._rs,  0._rs)
      t%p(11) = point(   phi, -1._rs,  0._rs)
      t%p(12) = point(  -phi, -1._rs,  0._rs)

      ! Indices describing location of triangle vertices
      t%t( 1)%a =  2 ;  t%t( 1)%b =  4 ;  t%t( 1)%c = 11
      t%t( 2)%a =  5 ;  t%t( 2)%b =  2 ;  t%t( 2)%c = 11
      t%t( 3)%a =  9 ;  t%t( 3)%b =  5 ;  t%t( 3)%c = 11
      t%t( 4)%a =  7 ;  t%t( 4)%b =  9 ;  t%t( 4)%c = 11
      t%t( 5)%a = 11 ;  t%t( 5)%b =  7 ;  t%t( 5)%c =  4
      t%t( 6)%a =  4 ;  t%t( 6)%b =  2 ;  t%t( 6)%c = 12
      t%t( 7)%a =  6 ;  t%t( 7)%b = 12 ;  t%t( 7)%c =  2
      t%t( 8)%a =  2 ;  t%t( 8)%b =  5 ;  t%t( 8)%c =  6
      t%t( 9)%a =  1 ;  t%t( 9)%b =  6 ;  t%t( 9)%c =  5
      t%t(10)%a =  5 ;  t%t(10)%b =  9 ;  t%t(10)%c =  1
      t%t(11)%a =  3 ;  t%t(11)%b =  1 ;  t%t(11)%c =  9
      t%t(12)%a =  9 ;  t%t(12)%b =  7 ;  t%t(12)%c =  3
      t%t(13)%a =  8 ;  t%t(13)%b =  3 ;  t%t(13)%c =  7
      t%t(14)%a =  7 ;  t%t(14)%b =  4 ;  t%t(14)%c =  8
      t%t(15)%a = 12 ;  t%t(15)%b =  8 ;  t%t(15)%c =  4
      t%t(16)%a = 12 ;  t%t(16)%b =  6 ;  t%t(16)%c = 10
      t%t(17)%a =  6 ;  t%t(17)%b =  1 ;  t%t(17)%c = 10
      t%t(18)%a =  1 ;  t%t(18)%b =  3 ;  t%t(18)%c = 10
      t%t(19)%a =  3 ;  t%t(19)%b =  8 ;  t%t(19)%c = 10
      t%t(20)%a = 10 ;  t%t(20)%b = 12 ;  t%t(20)%c =  8

      ! Normalise onto unit sphere
      call st_norm_p(t%p)

   ! Second orientation places two points at +/- z, then places points around small
   ! circles at z = +/-a = 1/sqrt(5).
   ! Coordinates are from http://mathworld.wolfram.com/Icosahedron.html
   else
      t%p(1)  = point(           0._rs,            0._rs,  1._rs)
      t%p(2)  = point(      b*cos(pi5),       b*sin(pi5),      a)
      t%p(3)  = point(b*cos(3._rs*pi5), b*sin(3._rs*pi5),      a)
      t%p(4)  = point(b*cos(5._rs*pi5), b*sin(5._rs*pi5),      a)
      t%p(5)  = point(b*cos(7._rs*pi5), b*sin(7._rs*pi5),      a)
      t%p(6)  = point(b*cos(9._rs*pi5), b*sin(9._rs*pi5),      a)
      t%p(7)  = point(               b,            0._rs,     -a)
      t%p(8)  = point(b*cos(2._rs*pi5), b*sin(2._rs*pi5),     -a)
      t%p(9)  = point(b*cos(4._rs*pi5), b*sin(4._rs*pi5),     -a)
      t%p(10) = point(b*cos(6._rs*pi5), b*sin(6._rs*pi5),     -a)
      t%p(11) = point(b*cos(8._rs*pi5), b*sin(8._rs*pi5),     -a)
      t%p(12) = point(           0._rs,            0._rs, -1._rs)

      ! Indices of triangle vertices
      t%t(1)  = triangle( 6,  2,  1)
      t%t(2)  = triangle( 2,  3,  1)
      t%t(3)  = triangle( 3,  4,  1)
      t%t(4)  = triangle( 4,  5,  1)
      t%t(5)  = triangle( 5,  6,  1)
      t%t(6)  = triangle(11,  7,  6)
      t%t(7)  = triangle( 7,  2,  6)
      t%t(8)  = triangle( 7,  8,  2)
      t%t(9)  = triangle( 8,  3,  2)
      t%t(10) = triangle( 8,  9,  3)
      t%t(11) = triangle( 9,  4,  3)
      t%t(12) = triangle( 9, 10,  4)
      t%t(13) = triangle(10,  5,  4)
      t%t(14) = triangle(10, 11,  5)
      t%t(15) = triangle(11,  6,  5)
      t%t(16) = triangle(12,  7, 11)
      t%t(17) = triangle(12,  8,  7)
      t%t(18) = triangle(12,  9,  8)
      t%t(19) = triangle(12, 10,  9)
      t%t(20) = triangle(12, 11, 10)
   endif

end function st_icosahedron
!-------------------------------------------------------------------------------

!===============================================================================
subroutine st_iterate_level(u)
!===============================================================================
! Increase the tesselation level by one.
! In the scheme I use, we take the old triangle and divide into four new ones.
! Point a takes new T1, b T2, c T3 and T4 is in the middle.
! The old triangle is replaced by T1.
! New points d, e and f are added to the end of the existing point list.
!             c
!             /\
!           /    \
!         /   3    \
!     f /------------\ e
!     /  \    4    /   \
!   /   1  \     /   2   \
! a _________\_/__________ b
!             d

   ! The current tesselation; this will be reallocated and refilled
   type(tesselation), intent(inout) :: u
   ! Arrays for the next level up
   type(triangle), allocatable, dimension(:) :: t
   type(point), allocatable, dimension(:) :: p
   type(point) :: a, b, c, d, e, f
   integer :: it, nt, np, level, ipa, ipb, ipc, ipd, ipe, ipf, it_last, ip_last

   ! Make room for the next level in temporary tesselation held in p and t
   level = u%level + 1
   if (debug) write(0,'(a,i0.1,a,i0.1,a)') 'Moving from level ', u%level, ' to level ', level
   nt = st_num_faces(level)
   np = st_num_points(level)
   if (debug) write(0,'(2(a,i0.1))') 'nt = ', nt, '   np = ', np
   allocate(t(nt))
   allocate(p(np))
   ! Fill up the new tesselation with the existing np points and nt faces
   t(1:u%nt) = u%t
   p(1:u%np) = u%p
   ! Go through each triangle; in each case, we create a new point in the middle
   ! of each of the edges and add these points to the end of the list.
   ! Then replace the original triangle with the one closest to the original point
   ! a, and add the others onto the end of the list.
   ip_last = u%np ! Current last point counter
   it_last = u%nt ! Current last triangle counter
   do it = 1, u%nt
      if (debug) write(0,'(a,i0.1)') 'Triangle ', it
      ! The existing points' indices
      ipa = t(it)%a
      ipb = t(it)%b
      ipc = t(it)%c
      ! The existing points
      a = p(ipa)
      b = p(ipb)
      c = p(ipc)
      if (debug) then
         write(0,'(a)') '   Corners:'
         write(0,'(a,3f9.5,1x,i0.1)') '      a: ', a, ipa
         write(0,'(a,3f9.5,1x,i0.1)') '      b: ', b, ipb
         write(0,'(a,3f9.5,1x,i0.1)') '      c: ', c, ipc
         write(0,'(a,i0.1,1x,i0.1)') '   Current # triangles, points: ', it_last, ip_last
      endif
      ! The halfway points along each side
      d = st_halfway_pt(a, b)
      e = st_halfway_pt(b, c)
      f = st_halfway_pt(c, a)
      ! If this is the first go round for this increment, then add to the point list
      if (it == 1) then
         ipd = ip_last + 1
         ipe = ip_last + 2
         ipf = ip_last + 3
         ip_last = ip_last + 3
      ! Otherwise, search for these points already existing and add if necessary
      else
         if (.not.st_point_exists(d, p(u%np+1:ip_last), ipd)) ip_last = ip_last + 1
         ipd = u%np + ipd
         if (.not.st_point_exists(e, p(u%np+1:ip_last), ipe)) ip_last = ip_last + 1
         ipe = u%np + ipe
         if (.not.st_point_exists(f, p(u%np+1:ip_last), ipf)) ip_last = ip_last + 1
         ipf = u%np + ipf
      endif
      if (debug) then
         write(0,'(a)') '   New corners:'
         write(0,'(a,3f9.5,1x,i0.1)') '      d: ', d, ipd
         write(0,'(a,3f9.5,1x,i0.1)') '      e: ', e, ipe
         write(0,'(a,3f9.5,1x,i0.1)') '      f: ', f, ipf
         write(0,'(a,i0.1,1x,i0.1)') '   Current # triangles, points: ', it_last, ip_last
      endif
      p(ipd) = d
      p(ipe) = e
      p(ipf) = f
      ! Update the old triangle (T1)
      t(it)%b = ipd
      t(it)%c = ipf
      ! Add the new triangles
      t(it_last + 1) = triangle(ipb, ipe, ipd) ! T2
      t(it_last + 2) = triangle(ipc, ipf, ipe) ! T3
      t(it_last + 3) = triangle(ipd, ipe, ipf) ! T4
      ! Update counter
      it_last = it_last + 3
   enddo
   if (ip_last /= np) then
      write(0,'(2(a,i0.1),a)') 'st_iterate_level: Error: Finished iteration but ' // &
         'point counter is not as expected (ip_last = ',ip_last,' not ',np,')'
      stop
   endif
   if (it_last /= nt) then
      write(0,'(2(a,i0.1),a)') 'st_iterate_level: Error: Finished iteration but ' // &
         'point counter is not as expected (it_last = ',it_last,' not ',nt,')'
      stop
   endif

   ! Reallocate old tesselation and fill in with new one
   deallocate(u%t)
   deallocate(u%p)
   allocate(u%t(nt))
   allocate(u%p(np))
   u%level = level
   u%nt = nt
   u%np = np
   u%t = t
   u%p = p
   ! Clear out temporary tesselation
   deallocate(t, p)

end subroutine st_iterate_level
!-------------------------------------------------------------------------------

!===============================================================================
subroutine st_compute_barycentres(t)
!===============================================================================
!  Assuming the barycentring calculation has not been done, compute the barycentres
!  of each triangle, filling in the t%b set of points
   type(tesselation), intent(inout) :: t
   integer :: i
   if (allocated(t%b)) then
      if (size(t%b) /= t%nt) then
         deallocate(t%b)
         allocate(t%b(t%nt))
      endif
   else
      allocate(t%b(t%nt))
   endif
   do i = 1, t%nt
      t%b(i) = st_barycentre([t%p(t%t(i)%a), t%p(t%t(i)%b), t%p(t%t(i)%c)])
   enddo
end subroutine st_compute_barycentres
!-------------------------------------------------------------------------------

!===============================================================================
function st_point_exists(p, a, ip) result(exists)
!===============================================================================
!  Return .true. if point p exists in the array of points a.
!  Optionally, return the index of the point in a in ip, if supplied.
!  If the point does not exist and ip is supplied, ip is set to one more than
!  the size of a.
   implicit none
   type(point), intent(in) :: p, a(:)
   integer, intent(out), optional :: ip
   logical :: exists
   integer :: i
   exists = .false.
   do i = 1, size(a)
      if (sqrt((p%x - a(i)%x)**2 + (p%y - a(i)%y)**2 + (p%z - a(i)%z)**2) &
          <= point_tol) then
         ! Set ip to the index of the matching point if it exists
         if (present(ip)) ip = i
         exists = .true.
         return
      endif
   enddo
   ! Set the index to the the next one up if it doesn't
   if (present(ip)) ip = size(a) + 1
end function st_point_exists
!-------------------------------------------------------------------------------

!===============================================================================
function st_halfway_pt(a, b) result(c)
!===============================================================================
! Return the halfway point between two other points
   type(point), intent(in) :: a, b
   type(point) :: c
   c = point((a%x+b%x)/2._rs, (a%y+b%y)/2._rs, (a%z+b%z)/2._rs)
   call st_norm_p(c)
end function st_halfway_pt
!-------------------------------------------------------------------------------

!===============================================================================
function st_num_faces(level) result(ntriangles)
!===============================================================================
!  Return the number of faces/triangles for a given tesselation level
   integer, intent(in) :: level
   integer :: ntriangles
   ntriangles = 20*2**(2*level)
end function st_num_faces
!-------------------------------------------------------------------------------

!===============================================================================
function st_num_points(level) result(npoints)
!===============================================================================
!  Return the number of points for a given tesselation level
   integer, intent(in) :: level
   integer :: npoints
   npoints = 10*2**(2*level) + 2
end function st_num_points
!-------------------------------------------------------------------------------

!===============================================================================
function st_level_from_points(n, warn) result(level)
!===============================================================================
!  Return the level of the tesselation from the number of points, and return
!  -1 if the number of points does not match an expected tesselation level.
   integer, intent(in) :: n
   integer :: level
   logical, intent(in), optional :: warn
   logical :: warn_in
   warn_in = .false.
   if (present(warn)) warn_in = warn
   level = 0
   do while (level <= 10)
      if (st_num_points(level) == n) return
      level = level + 1
   enddo
   level = -1
   if (warn_in) write(0,'(a,i0.1)') &
      'st_level_from_points: Warning: Not an expected number of points: ', n
end function st_level_from_points
!-------------------------------------------------------------------------------

!===============================================================================
subroutine st_norm_p_array(p)
!===============================================================================
!  Normalise an array of points onto the unit sphere
!  This is overloaded with the st_norm_p interface
   type(point), intent(inout) :: p(:)
   integer :: i
   real(rs) :: r
   do i = 1, size(p)
      r = sqrt(p(i)%x**2 + p(i)%y**2 + p(i)%z**2)
      p(i)%x = p(i)%x/r
      p(i)%y = p(i)%y/r
      p(i)%z = p(i)%z/r
   enddo
end subroutine st_norm_p_array
!-------------------------------------------------------------------------------

!===============================================================================
subroutine st_norm_p_single(p)
!===============================================================================
!  Normalise a point onto the unit sphere
!  This is overloaded with the st_norm_p interface
   type(point), intent(inout) :: p
   real(rs) :: r
   r = sqrt(p%x**2 + p%y**2 + p%z**2)
   p%x = p%x/r
   p%y = p%y/r
   p%z = p%z/r
end subroutine st_norm_p_single
!-------------------------------------------------------------------------------

!===============================================================================
subroutine st_rotate(t, alpha, beta, gamma, degrees)
!===============================================================================
!  Rotate the points in a tesselation about the x, y and z axes by a, b and c.
!  Angles are in radians, unless degrees=.true. is passed in.
!  The rotations are applied, in order, about the x, y and z axes, anticlockwise
!  when looking down these axes (right-hand rule).

   use spherical_geometry, only: sg_torad
   type(tesselation), intent(inout) :: t
   real(rs), intent(in) :: alpha, beta, gamma
   logical, intent(in), optional :: degrees
   real(rs) :: a, b, c, R(3,3), v(3)
   integer :: i

   a = alpha
   b = beta
   c = gamma
   if (present(degrees)) then
      if (degrees) then
         a = sg_torad(a)
         b = sg_torad(b)
         c = sg_torad(c)
      endif
   endif

   ! Make rotation matrix and apply to each point
   R = st_rotmat(a, b, c)
   do i = 1, t%np
      v = matmul(R, [t%p(i)%x, t%p(i)%y, t%p(i)%z])
      t%p(i) = point(v(1), v(2), v(3))
   enddo
end subroutine st_rotate
!-------------------------------------------------------------------------------

!===============================================================================
function st_rotmat(a, b, c) result(R)
!===============================================================================
!  Return a 3x3 rotation matrix given three angles a, b and c.
!  They represent, in turn, an anticlockwise rotation about
!  the x, y and z axes when looking down towards the origin.
!  Angles are in radians.
   real(rs), intent(in) :: a, b, c
   real(rs), dimension(3,3) :: Rx, Ry, Rz, R
   real(rs) :: sina, sinb, sinc, cosa, cosb, cosc

   sina = sin(a); sinb = sin(b); sinc = sin(c)
   cosa = cos(a); cosb = cos(b); cosc = cos(c)
   Rx = transpose(reshape([1._rs, 0._rs, 0._rs, &
                           0._rs,  cosa, -sina, &
                           0._rs,  sina,  cosa], [3,3]))
   Ry = transpose(reshape([ cosb, 0._rs,  sinb, &
                           0._rs, 1._rs, 0._rs, &
                           -sinb, 0._rs,  cosb], [3,3]))
   Rz = transpose(reshape([ cosc, -sinc, 0._rs, &
                            sinc,  cosc, 0._rs, &
                           0._rs, 0._rs, 1._rs], [3,3]))
   R = matmul(Rz, matmul(Ry, Rx))
end function st_rotmat
!-------------------------------------------------------------------------------

!===============================================================================
function st_barycentre(p)
!===============================================================================
!  Return the barycentre of an array of points, normalised to the radius of 1.
   type(point), intent(in), dimension(:) :: p
   type(point) :: st_barycentre
   integer :: n
   n = size(p)
   st_barycentre = point(sum(p%x), sum(p%y), sum(p%z))
   call st_norm_p(st_barycentre)
end function st_barycentre
!-------------------------------------------------------------------------------

!===============================================================================
function tolower(a) result(b)
!===============================================================================
!  Convert the uppercase characters in a string to lowercase
   implicit none
   character(len=*), intent(in) :: a
   character(len=len(a)) :: b
   integer, parameter :: capA = iachar('A'), capZ = iachar('Z')
   integer :: i, ia
   b = a
   do i = 1, len_trim(a)
      ia = iachar(a(i:i))
      if (ia >= capA .and. ia <= capZ) &
         b(i:i) = achar(iachar(a(i:i)) + 32)
   enddo
end function tolower
!-------------------------------------------------------------------------------
end module sphere_tesselate
!-------------------------------------------------------------------------------
