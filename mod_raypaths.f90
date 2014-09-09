!===============================================================================
! One of Andy Nowacki's Fortran utility modules for dealing with seismic
! anisotropy and other problems.
!
! Andy Nowacki <andy.nowacki@bristol.ac.uk>
!
! See the file LICENCE for licence details.
!===============================================================================
module raypaths
!===============================================================================
!  Module which contains definitions of raypath structure and subroutines to
!  manipulate them

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
      
!===============================================================================
!  Define the raypath type
   type raypath
      real(rs),allocatable,dimension(:) :: x,y,z
      integer                           :: N
   end type raypath
!===============================================================================

   contains
   
!===============================================================================
   subroutine raypaths_read(file,ray,taup)
!===============================================================================
!  Read in the coordinates of a ray from file.
!  These ray files have x,y,z coordinates of ray on each line, with comment 
!  lines containing #, % or > as the first character on the line.
!
!  An alternative input format is for output from taup_path.  These have format:
!     > name info
!     distance/deg radius/km lat/deg lon/deg
!  Specifying taup=.true. as an argument enables this conversion
   
   implicit none
   
   character(len=*),intent(in)  :: file
   type(raypath),intent(inout)  :: ray
   character(len=1)             :: comment
   integer                      :: npts,iostatus
   logical, intent(in), optional :: taup
   logical :: taup_input
   real(rs) :: delta,r,lon,lat,colon,colat
   
!  Check if we have taup_input; default to false
   taup_input = .false.
   if (present(taup)) taup_input = taup
   
!  Check that this is a 'new' ray
   if (allocated(ray%x) .or. allocated(ray%y) .or. allocated(ray%z)) then
      write(0,'(2a)') 'raypaths_read: require ''new'' raypath structure, but ',&
                     '               memory already allocated.'
      stop
   endif
   
!  Open file
   open(20,file=trim(file),status='old',iostat=iostatus)
   if (iostatus /=0) then
      write(0,'(3a)') 'raypaths_read: file "',trim(file),'" could not be opened for reading.'
      stop
   endif
   
!  Work out number of points along ray
   npts = 0
   do while (iostatus == 0)
      read(20,fmt='(1a)',iostat=iostatus) comment
      if (iostatus < 0) exit
      if (iostatus > 0) then
         write(0,'(a,a)') 'raypaths_read: Problem reading from file ',trim(file)
         stop
      endif
      if (comment == '#' .or. comment == '%' .or. comment == '>') cycle
      npts = npts + 1
   enddo
   
!  Allocate space for memory
   allocate(ray%x(npts), ray%y(npts), ray%z(npts))
   ray%N = npts
   
!  Read in data
   npts = 1
   rewind(20)
   do while (npts <= ray%N)
      read(20,'(1a)') comment
      if (comment == '#' .or. comment == '%' .or. comment == '>') cycle
      backspace(20)
      if (taup_input) then
         read(20,*) delta, r, lat, lon
         colat = to_rad*(90._rs - lat)
         colon = to_rad*lon
         ray%x(npts) = r * sin(colat) * cos(colon)
         ray%y(npts) = r * sin(colat) * sin(colon)
         ray%z(npts) = r * cos(colat)
      else
         read(20,*) ray%x(npts), ray%y(npts), ray%z(npts)
      endif
      npts = npts + 1
   enddo
   
   close(20)

   return
   
   end subroutine raypaths_read
!-------------------------------------------------------------------------------

!===============================================================================
   subroutine raypaths_write(file,ray)
!===============================================================================
!  Write out a ray path to a file
      implicit none
      type(raypath),intent(in) :: ray
      character(len=*),intent(in) :: file
      integer :: i,iostatus
      
      open(20,file=trim(file),iostat=iostatus)
      if (iostatus /= 0) then
         write(0,'(a)') 'raypaths_write: Problem opening file "'//trim(file)//'"'
         stop 1
      endif
      do i=1,ray%N
         write(20,*) ray%x(i), ray%y(i), ray%z(i)
      enddo
      close(20)
   end subroutine raypaths_write
!-------------------------------------------------------------------------------

!===============================================================================
   subroutine raypaths_clone(in,out)
!===============================================================================
!  Create a new raypath from an old one
     implicit none
     type(raypath),intent(in)  :: in
     type(raypath),intent(out) :: out
     
     if (allocated(out%x) .or. allocated(out%y) .or. allocated(out%z)) then
       write(0,'(2a)') 'raypaths_clone: expecting ''new'' raypath structure, ',&
                   'but memory already allocated for it.'
       stop
     endif
     
     allocate(out%x(in%N), out%y(in%N), out%z(in%N))
     out % N = in % N
     out % x = in % x
     out % y = in % y
     out % z = in % z
     
     return
   end subroutine raypaths_clone
!-------------------------------------------------------------------------------

!===============================================================================
   subroutine raypaths_delete(ray)
!===============================================================================
!  Delete a raypath and deallocate memory
     implicit none
     type(raypath),intent(inout) :: ray
     
     deallocate(ray%x,ray%y,ray%z)
     ray%N=0
     
     return
   end subroutine raypaths_delete
!-------------------------------------------------------------------------------

!===============================================================================
   subroutine raypaths_translate(ray,x,y,z)
!===============================================================================
!  Translates a raypath structure to a different location in Cartesian space.
     implicit none
     type(raypath),intent(inout) :: ray
     real(rs),intent(in)         :: x,y,z
     
     ray%x = ray%x + x
     ray%y = ray%y + y
     ray%z = ray%z + z
     
     return
   end subroutine raypaths_translate
!-------------------------------------------------------------------------------

!===============================================================================
   subroutine raypaths_rotate(ray,a_in,b_in,c_in,x,y,z,degrees)
!===============================================================================
!  Rotate a raypath about a certain point in the grid (optional; default is origin)
!  Rotation a is about x1 axis, +ve from 3 -> 2
!  Rotation b is about x2 axis, +ve from 1 -> 3
!  Rotation c is about x3 axis, +ve from 1 -> 2 (anticlock)
!  Rotations are applied in this order.
!  Angles are in degrees by default, but can be supplied by stating degrees=.false.

   implicit none
   
   type(raypath),intent(inout) :: ray
   real(rs),intent(in)         :: a_in,b_in,c_in
   real(rs),optional           :: x,y,z
   real(rs)                    :: a,b,c
   real(rs)                    :: v(3)
   real(rs),dimension(3,3)     :: R,Ra,Rb,Rc
   logical,optional            :: degrees
   integer                     :: i
   
!  Get inputs and convert to radians unless otherwise
   a = a_in*pi/180.
   b = b_in*pi/180.
   c = c_in*pi/180.
   if (present(degrees)) then
      if (.not.degrees) then
         a = a_in
         b = b_in
         c = c_in
      endif
   endif
   
!  Build the rotation matrices
   Ra(1,1) = 1.       ; Ra(1,2) = 0.       ; Ra(1,3) = 0.
   Ra(2,1) = 0.       ; Ra(2,2) = cos(a)   ; Ra(2,3) = -sin(a)
   Ra(3,1) = 0.       ; Ra(3,2) = sin(a)   ; Ra(3,3) = cos(a)
   
   Rb(1,1) = cos(b)   ; Rb(1,2) = 0.       ; Rb(1,3) = sin(b)
   Rb(2,1) = 0.       ; Rb(2,2) = 1.       ; Rb(2,3) = 0.
   Rb(3,1) = -sin(b)  ; Rb(3,2) = 0.       ; Rb(3,3) = cos(b)
   
   Rc(1,1) = cos(c)   ; Rc(1,2) = -sin(c)  ; Rc(1,3) = 0.
   Rc(2,1) = sin(c)   ; Rc(2,2) = cos(c)   ; Rc(2,3) = 0.
   Rc(3,1) = 0.       ; Rc(3,2) = 0.       ; Rc(3,3) = 1.
   
   R = matmul(Ra,Rb)
   R = matmul(Rc,R)
   
!  To simplify the rotation, shift so the rotation centre is at the origin
   if (present(x) .and. present(y) .and. present(z)) &
      call raypaths_translate(ray,-x,-y,-z)
   
!  Apply the rotations
   do i=1,ray%N
      v(1) = ray%x(i)
      v(2) = ray%y(i)
      v(3) = ray%z(i)
      v = matmul(R,v)
      ray%x(i) = v(1)
      ray%y(i) = v(2)
      ray%z(i) = v(3)
   enddo
   
!  Shift back so the centre of rotation is where it started
   if (present(x) .and. present(y) .and. present(z)) &
      call raypaths_translate(ray,x,y,z)
   
   return
   
   end subroutine raypaths_rotate
!-------------------------------------------------------------------------------
  
   

!_______________________________________________________________________________
end module raypaths