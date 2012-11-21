!==============================================================================
module EC_grid
!==============================================================================
!  Module containing routines for manipulating regularly-spaced grids of elastic
!  constants.  Defines a type for them and can do various things such as load
!  these structures into memory, write out to and read from files.
!
!  Andy Nowacki, University of Bristol
!  andy.nowacki@bristol.ac.uk
!  2011/02
!
!  UPDATES:
!     * 2011-03-16   Adding a type for grids with real coordinates and spacing.
!                    This will require a stricter file format, as the limits can't
!                    be determined just by reading the file and comparing
!                    coordinates (STILL TO DO).
!     * 2011-07-26   Implemented and checked the binary EC_grid file functions
!                    EC_grid_load_bin and EC_grid_write_bin.  These aren't very
!                    quick, so should consider C routines instead...
!     * 2012-03-01   Begun implementation of new EC_grid types which include
!                    density.  This breaks backward compatibility with old binary
!                    and ASCII grid files, but these can be converted fairly
!                    easily.
!     * 2012-11-20   Begin testing of EC_gridr routines, with only binary files 
!                    allowed at this stage.  Introduce a 'nver' variable at start
!                    of binary files to distinguish from assumed-integer.
!-------------------------------------------------------------------------------
!  FORTRAN BINARY FORMAT FILES
!  
!  These are written in one of two ways.  The correct way is as follows:
!   [1] nver,npts,nx,ny,nz,x1,x2,y1,y2,z1,z2,dx,dy,dz
!   [2+] c11,c12...c16,c22...c66 (loops over z slowest, x fastest)
!  
!  nver is 1 for explicit integer coordinates, and 2 for real coordinates.
!  (Real coordinates are the default.)
!  
!  If nver is not given (i.e., the first field of the first line is npts), then
!  it is assumed that the binary file describes an assumed-integer coordinate frame.
!  The format is:
!   [1] npts,nx,ny,nz,ix1,ix2,iy1,iy2,iz1,iz2,idx,idy,idz
!   [2+] (as above)
!
!  This will be deprecated in the future.
!-------------------------------------------------------------------------------  

   implicit none
   
!  ** size constants
   integer, parameter, private :: i4 = selected_int_kind(9)       ! long int
   integer, parameter, private :: r4 = selected_real_kind(6,37)   ! SP
   integer, parameter, private :: r8 = selected_real_kind(15,307) ! DP

!  ** precision selector
   integer, parameter, private :: rs = r8
   
!  ** maths constants and other useful things
   real(rs), parameter, private :: pi = 3.141592653589793238462643_rs
   integer, parameter, private :: big_p_integer = 2147483647
   integer, parameter, private :: big_n_integer = -2147483647
   
!  ** IO constants
   integer, parameter, private :: lu_stdin  = 5
   integer, parameter, private :: lu_stdout = 6
   integer, parameter, private :: lu_stderr = 0
   
!  ** Version constants for explicit real and binary grid types (used internally)
   integer, parameter, private :: nver_assumed_i = 0
   integer, parameter, private :: nver_i = 1
   integer, parameter, private :: nver_r = 2
!  ** Version constants for external use
   integer, parameter, public :: EC_grid_nver_assumed_i = nver_assumed_i
   integer, parameter, public :: EC_grid_nver_i = nver_i
   integer, parameter, public :: EC_grid_nver_r = nver_r
   
!  ** Limit on maximum size of arrays: currently 2 GB
   real, parameter, private :: MAX_GRID_SIZE_IN_BYTES = 2.*1024.**3
   
!------------------------------------------------------------------------------
!  Define the ECgrid type: assumed integer coordinates
   type :: ECgrid
      integer(i4) :: npts,nx,ny,nz,ix1,ix2,iy1,iy2,iz1,iz2,idx,idy,idz
      real(r8),allocatable ::  ecs(:,:,:,:,:)
      real(r8),allocatable ::  x(:),y(:),z(:)
   end type ECgrid
   
!  ECgridi: explicitly integer coordinates
   type :: ECgridi
      integer(i4) :: npts,nx,ny,nz,ix1,ix2,iy1,iy2,iz1,iz2,idx,idy,idz
      real(r8),allocatable ::  ecs(:,:,:,:,:)
      real(r8),allocatable ::  x(:),y(:),z(:)
   end type ECgridi
   
!  ECgridr type: real coordinates
   type :: ECgridr
      integer(i4) :: npts,nx,ny,nz
      real(r8)    :: x1,x2,y1,y2,z1,z2,dx,dy,dz
      real(r8),allocatable :: ecs(:,:,:,:,:)
      real(r8),allocatable :: x(:),y(:),z(:)
   end type ECgridr
!------------------------------------------------------------------------------
!  Define the ECRgrid types: These also include density
!  ECRgrid: assumed real coordinates
   type :: ECRgrid
      integer(i4) :: npts,nx,ny,nz
      real(r8)    :: x1,x2,y1,y2,z1,z2,dx,dy,dz
      real(r8),allocatable :: ecs(:,:,:,:,:)
      real(r8),allocatable :: rho(:,:,:)
      real(r8),allocatable :: x(:),y(:),z(:)
   end type ECRgrid

!  ECRgridi: explicit integer coordinates
   type :: ECRgridi
      integer(i4) :: npts,nx,ny,nz,ix1,ix2,iy1,iy2,iz1,iz2,idx,idy,idz
      real(r8),allocatable ::  ecs(:,:,:,:,:)
      real(r8),allocatable ::  rho(:,:,:)
      real(r8),allocatable ::  x(:),y(:),z(:)
   end type ECRgridi
!-------------------------------------------------------------------------------   
   
   CONTAINS
   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!                         ASSMUED INTEGER COORDINATES
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!===============================================================================
   subroutine EC_grid_new(ix1,ix2,idx,iy1,iy2,idy,iz1,iz2,idz,grid)
!===============================================================================
!  Create a new grid

   implicit none
   
   type(ECgrid),intent(out)  :: grid
   integer,intent(in)        :: ix1,ix2,idx,iy1,iy2,idy,iz1,iz2,idz
   integer                   :: k
   
   if ((ix1>=ix2).or.(iy1>=iy2).or.(iz1>=iz2).or.(idx<=0).or.(idy<=0).or.(idz<=0)) then
      write(lu_stderr,'(a)') 'EC_grid_new: problem with array limits or spacing.'
      write(lu_stderr,'(a)') '    ix1<ix2, iy1<iy2, iz1<iz1, (idx,idy,idz)>0'
      stop
   endif
   
   grid % ix1 = ix1  ;  grid % ix2 = ix2
   grid % iy1 = iy1  ;  grid % iy2 = iy2
   grid % iz1 = iz1  ;  grid % iz2 = iz2
   grid % idx = idx
   grid % idy = idy
   grid % idz = idz
   
!  Calcuate array dimensions
   grid % nx = (ix2 - ix1) / idx + 1
   grid % ny = (iy2 - iy1) / idy + 1
   grid % nz = (iz2 - iz1) / idz + 1
   grid % npts = grid%nx * grid%ny * grid%nz
   
   if (allocated(grid%ecs)) deallocate(grid%ecs)
   if (allocated(grid%x)) deallocate(grid%x)
   if (allocated(grid%y)) deallocate(grid%y)
   if (allocated(grid%z)) deallocate(grid%z)

!  Allocate memory   
   allocate(grid%ecs(grid%nx,grid%ny,grid%nz,6,6))
   allocate(grid%x(grid%nx), grid%y(grid%ny), grid%z(grid%nz))
   
!  Fill in dimensions
   grid%ecs = 0.
   do k=0,grid%nx-1;  grid % x(k+1) = real(ix1) + real(idx) * real(k);   enddo
   do k=0,grid%ny-1;  grid % y(k+1) = real(iy1) + real(idy) * real(k);   enddo
   do k=0,grid%nz-1;  grid % z(k+1) = real(iz1) + real(idz) * real(k);   enddo
   
   return
   end subroutine EC_grid_new
!-------------------------------------------------------------------------------

!==============================================================================
   subroutine EC_grid_load(ecs_file,grid,necs,quiet)
!==============================================================================
!  Load a 3D grid of elastic constants into memory.
!  Elastic constants should be oriented according to the grid, where
!     x // 1   ,   y // 2   ,   z // 3, 
!  assuming a right-handed coordinate system.  Later tranformations can be performed.

   implicit none

   type(ECgrid),intent(out)    :: grid
   character(len=*),intent(in) :: ecs_file
   integer   :: ix,iy,iz,ix1,ix2,iy1,iy2,iz1,iz2,idx,idy,idz,&
                itempx,itempy,itempz,nx,ny,nz,npts,i,j,k,nec
   integer,optional,intent(in) :: necs
   logical,optional,intent(in) :: quiet
   logical                     :: silent
   real(rs)   :: ecs_in(6,6),multiplier
   integer    :: iostatus
   
!  Check optional arguments
   if (present(necs)) then
      if (necs /=21 .and. necs /= 36) stop 'EC_grid_load: nec must be 21 or 36'
      nec = necs
   else
      nec = 36
   endif
   
   silent = .false.
   if (present(quiet)) then
      silent = quiet
   endif
   
!  Open file
   open(20,file=ecs_file,status='old',iostat=iostatus)
   if (iostatus /= 0) then
      write(0,'(a,a)') 'Problem opening ECs file ',ecs_file
      stop
   endif  
   
   iostatus = 0
   ix1 = big_p_integer ; ix2 = big_n_integer
   iy1 = big_p_integer ; iy2 = big_n_integer
   iz1 = big_p_integer ; iz2 = big_n_integer
   npts = 0

!  Read limits of data and find dimensions
   do while (iostatus == 0)
      read(20,fmt=*,iostat=iostatus) ix,iy,iz
      if (iostatus /= 0) exit
      if (ix < ix1) ix1 = ix
      if (ix > ix2) ix2 = ix
      if (iy < iy1) iy1 = iy
      if (iy > iy2) iy2 = iy
      if (iz < iz1) iz1 = iz
      if (iz > iz2) iz2 = iz
      npts = npts + 1
   enddo
   
   if (.not.silent) then
     write(0,'(a)') '=================================='
     write(0,'(a)') 'Dimensions of box are:'
     write(0,'(a,2i5)') '  x:  ',ix1,ix2
     write(0,'(a,2i5)') '  y:  ',iy1,iy2
     write(0,'(a,2i5)') '  z:  ',iz1,iz2
     write(0,'(a)') '=================================='
   endif
   
!  Start again and work out the dimensions
   rewind(20)
!  Count the number of points with the same x and z as the first line.
   read(20,*) itempx,itempy,itempz
   ny = 1
   do i=2,npts
      read(20,fmt=*,iostat=iostatus) ix,iy,iz
      if (ix == itempx .and. iz == itempz) ny = ny + 1
   enddo
   
!  Count the number of x at this y and z
   rewind(20)
   read(20,*) itempx,itempy,itempz
   nx = 1
   do i=2,npts
      read(20,fmt=*,iostat=iostatus) ix,iy,iz
      if (iy == itempy .and. iz == itempz) nx = nx + 1
   enddo
   
!  Count the number of z at this x and y: just for now, to check!
   rewind(20)
   read(20,*) itempx,itempy,itempz
   nz = 1
   do i=2,npts
      read(20,fmt=*,iostat=iostatus) ix,iy,iz
      if (ix == itempx .and. iy == itempy) nz = nz + 1
   enddo
   
!  Sanity check
   if (nx*ny*nz /= npts) then
      write(0,'(a)') 'load_ecs_cart: Number of points in EC file does not match dimensions!'
      stop
   endif
   
   idx = (ix2 - ix1) / (nx-1)
   idy = (iy2 - iy1) / (ny-1)
   idz = (iz2 - iz1) / (nz-1)
   if (.not.silent) then
     write(0,'(3(a,i4))') 'Grid size is    ',nx,' x ',ny,' x ',nz
     write(0,'(3(a,i3))') 'Grid step size is ',idx,' x ',idy,' x ',idz
     write(0,'(a,i6)')    'Total number of ECs is ',npts
     write(0,'(a/)')      '=================================='
   endif
   
!  Allocate space for the ECs
   allocate(grid%ecs(nx,ny,nz,6,6),stat=iostatus)
   allocate(grid%x(nx), grid%y(ny), grid%z(nz))
   if (iostatus /= 0) stop 'EC_grid_load: error allocating space for ECs'
   
!  Read each line and fill the array
   rewind(20)
   do k=1,npts
      if (nec == 36) then
         read(20,*) itempx,itempy,itempz,ecs_in
      else if (nec == 21) then
         read(20,*) itempx,itempy,itempz,((ecs_in(i,j),j=i,6),i=1,6)
         do i=1,6
            do j=1,6
               ecs_in(j,i) = ecs_in(i,j)
            enddo
         enddo
      endif
!     For the first layer, determine whether we're in Pa or GPa; issue a warning regardless.
      if (k == 1) then
         if (ecs_in(1,1) > 1.d5) then
            multiplier = 1.d0
            if (.not.silent) &
            write(0,'(a)') '  EC_grid_load: Input Cij file appears to be in Pa.'
         else
            multiplier = 1.d9
            if (.not.silent) &
            write(0,'(a)') '  EC_grid_load: Input Cij file appears to be in GPa.'
         endif
      endif

!  Fill the array with the constants in Pa
      grid%ecs((itempx-ix1+idx)/idx,&
               (itempy-iy1+idy)/idy,&
               (itempz-iz1+idz)/idz,:,:) = ecs_in(:,:) * multiplier
      
   enddo
   
   close(20)

   grid%npts = nx*ny*nz
   grid % nx = nx  ;   grid % ny = ny  ;   grid % nz = nz
   grid % idx = idx  ;  grid % idy = idy  ;   grid % idz = idz
   grid % ix1 = ix1  ;  grid % ix2 = ix2
   grid % iy1 = iy1  ;  grid % iy2 = iy2
   grid % iz1 = iz1  ;  grid % iz2 = iz2
   
!  Fill in dimensions
   do k=0,nx-1;  grid % x(k+1) = real(ix1) + real(idx) * real(k);   enddo
   do k=0,ny-1;  grid % y(k+1) = real(iy1) + real(idy) * real(k);   enddo
   do k=0,nz-1;  grid % z(k+1) = real(iz1) + real(idz) * real(k);   enddo
   
   end subroutine EC_grid_load
!-------------------------------------------------------------------------------

!===============================================================================
   subroutine EC_grid_inquire(ecs_file,grid,quiet)
!===============================================================================
!  Find out about a grid and return its dimensions.  Does not allocate memory or
!  read in the array or coordinates, however
   
   implicit none
   type(ECgrid),intent(out)    :: grid
   character(len=*),intent(in) :: ecs_file
   integer   :: ix,iy,iz,ix1,ix2,iy1,iy2,iz1,iz2,idx,idy,idz,&
                itempx,itempy,itempz,nx,ny,nz,npts,i
   logical,optional,intent(in) :: quiet
   logical                     :: silent
   integer    :: iostatus
   
!  Check optional arguments
   silent = .false.
   if (present(quiet)) then
      silent = quiet
   endif
   
!  Open file
   open(20,file=ecs_file,status='old',iostat=iostatus)
   if (iostatus /= 0) then
      write(0,'(a,a)') 'EC_grid_inquire: Problem opening ECs file ',ecs_file
      stop
   endif  
   
   iostatus = 0
   ix1 = big_p_integer ; ix2 = big_n_integer
   iy1 = big_p_integer ; iy2 = big_n_integer
   iz1 = big_p_integer ; iz2 = big_n_integer
   npts = 0

!  Read limits of data and find dimensions
   do while (iostatus == 0)
      read(20,fmt=*,iostat=iostatus) ix,iy,iz
      if (iostatus /= 0) exit
      if (ix < ix1) ix1 = ix
      if (ix > ix2) ix2 = ix
      if (iy < iy1) iy1 = iy
      if (iy > iy2) iy2 = iy
      if (iz < iz1) iz1 = iz
      if (iz > iz2) iz2 = iz
      npts = npts + 1
   enddo
   
   if (.not.silent) then
     write(lu_stdout,'(a)')         '=================================='
     write(lu_stdout,'(a)')         'Dimensions of box are:'
     write(lu_stdout,'(a,i0,x,i0)') '  x:  ',ix1,ix2
     write(lu_stdout,'(a,i0,x,i0)') '  y:  ',iy1,iy2
     write(lu_stdout,'(a,i0,x,i0)') '  z:  ',iz1,iz2
     write(lu_stdout,'(a)')         '=================================='
   endif
   
!  Start again and work out the dimensions
   rewind(20)
!  Count the number of points with the same x and z as the first line.
   read(20,*) itempx,itempy,itempz
   ny = 1
   do i=2,npts
      read(20,fmt=*,iostat=iostatus) ix,iy,iz
      if (ix == itempx .and. iz == itempz) ny = ny + 1
   enddo
   
!  Count the number of x at this y and z
   rewind(20)
   read(20,*) itempx,itempy,itempz
   nx = 1
   do i=2,npts
      read(20,fmt=*,iostat=iostatus) ix,iy,iz
      if (iy == itempy .and. iz == itempz) nx = nx + 1
   enddo
   
!  Count the number of z at this x and y: just for now, to check!
   rewind(20)
   read(20,*) itempx,itempy,itempz
   nz = 1
   do i=2,npts
      read(20,fmt=*,iostat=iostatus) ix,iy,iz
      if (ix == itempx .and. iy == itempy) nz = nz + 1
   enddo
   
!  Sanity check
   if (nx*ny*nz /= npts) then
      write(0,'(a)') 'EC_grid_inquire: Number of points in EC file does not match dimensions!'
      stop
   endif
   
   idx = (ix2 - ix1) / (nx-1)
   idy = (iy2 - iy1) / (ny-1)
   idz = (iz2 - iz1) / (nz-1)
   if (.not.silent) then
     write(lu_stdout,'(3(a,i0))') 'Grid size is    ',nx, ' x ',ny, ' x ',nz
     write(lu_stdout,'(3(a,i0))') 'Grid spacing is ',idx,' x ',idy,' x ',idz
     write(lu_stdout,'(a,i0)')    'Total number of ECs is ',npts
     write(lu_stdout,'(a/)')      '=================================='
   endif
   
!  Fill in return values for grid, apart from coordinates and ecs
   grid % nx = nx  ;   grid % ny = ny  ;   grid % nz = nz
   grid % idx = idx  ;  grid % idy = idy  ;   grid % idz = idz
   grid % ix1 = ix1  ;  grid % ix2 = ix2
   grid % iy1 = iy1  ;  grid % iy2 = iy2
   grid % iz1 = iz1  ;  grid % iz2 = iz2
   

   end subroutine EC_grid_inquire
!-------------------------------------------------------------------------------

!===============================================================================
   subroutine EC_grid_inquire_bin(ecs_file,grid,quiet)
!===============================================================================
!  Get information about a binary grid file
   
   implicit none
   
   type(ECgrid),intent(inout)  :: grid
   character(len=*),intent(in) :: ecs_file
   logical,intent(in),optional :: quiet
   logical                     :: silent
   integer                     :: iostatus
   
   silent = .true.
   if (present(quiet)) silent = quiet
   
   iostatus = 0
   open(20,file=trim(ecs_file),iostat=iostatus,form='unformatted',status='old')
   if (iostatus /= 0) then
      write(lu_stderr,'(a,a)') &
         'EC_grid_inquire_bin: cannot open file for reading: ', trim(ecs_file)
      stop
   endif
   
!  Read header line
   read(20,iostat=iostatus) grid%npts, grid%nx, grid%ny, grid%nz, grid%ix1, grid%ix2, &
             grid%iy1, grid%iy2, grid%iz1, grid%iz2, grid%idx, grid%idy, grid%idz
   close(20)
   if (iostatus /= 0) then
      write(lu_stderr,'(a,a)') &
         'EC_grid_inquire_bin: cannot read header line for binary file ',trim(ecs_file)
      stop
   endif
   
! If desired, print out information
   if (.not.silent) then
     write(lu_stdout,'(a)')         '=================================='
     write(lu_stdout,'(a)')         'Dimensions of box are:'
     write(lu_stdout,'(a,i0,x,i0)') '  x:  ',grid%ix1,grid%ix2
     write(lu_stdout,'(a,i0,x,i0)') '  y:  ',grid%iy1,grid%iy2
     write(lu_stdout,'(a,i0,x,i0)') '  z:  ',grid%iz1,grid%iz2
     write(lu_stdout,'(a)')         '=================================='
     write(lu_stdout,'(3(a,i0))') 'Grid size is    ',grid%nx, ' x ',grid%ny, ' x ',grid%nz
     write(lu_stdout,'(3(a,i0))') 'Grid spacing is ',grid%idx,' x ',grid%idy,' x ',grid%idz
     write(lu_stdout,'(a,i0)')    'Total number of ECs is ',grid%npts
     write(lu_stdout,'(a/)')      '=================================='
   endif
   
   end subroutine EC_grid_inquire_bin
!-------------------------------------------------------------------------------

!===============================================================================
   subroutine EC_grid_write(fname,grid,necs)
!===============================================================================
!  Write the EC_grid type to an ASCII text file.
!  Assumes output is all 36 elastic constants, otherwise can specify 21.

   implicit none
   
   type(ECgrid),intent(in)       :: grid
   character(len=*),intent(in)   :: fname
   character(len=80)             :: fmt
   integer,optional,intent(in)   :: necs
   integer                       :: iostatus,ix,iy,iz,kx,ky,kz,i,j,nec
   integer                       :: lenx,leny,lenz
   
!  Check for requested output format; default to all 36 constants to be written out.
   if (present(necs)) then
      if (necs /= 21 .and. necs /= 36) stop 'EC_grid_write: nec must be 21 or 36'
      nec = necs
   else
      nec = 36
   endif
   
   open(20,file=trim(fname),iostat=iostatus)
   if (iostatus /= 0) stop 'EC_grid_write: error opening file for writing'
   
!  Determine output format; minimum length is 5
   lenx = 5   ;   leny = lenx   ;   lenz = lenx
   if (int(log10(real(grid%ix1))+1.) > lenx) lenx = int(log10(real(grid%ix1))+1.)
   if (int(log10(real(grid%ix2))+1.) > lenx) lenx = int(log10(real(grid%ix2))+1.)
   if (int(log10(real(grid%iy1))+1.) > leny) leny = int(log10(real(grid%iy1))+1.)
   if (int(log10(real(grid%iy2))+1.) > leny) leny = int(log10(real(grid%iy2))+1.)
   if (int(log10(real(grid%iz1))+1.) > lenz) lenz = int(log10(real(grid%iz1))+1.)
   if (int(log10(real(grid%iz2))+1.) > lenz) lenz = int(log10(real(grid%iz2))+1.)
   
   write(fmt,'(3(a4,i3),a1,i2,a)') '  (i',lenx,',x,i',leny,',x,i',lenz,',',nec,'e24.16)'  

   do kx=1,grid%nx
      ix = grid%ix1 + grid%idx*(kx-1)
      do ky=1,grid%ny
         iy = grid%iy1 + grid%idy*(ky-1)
         do kz=1,grid%nz
            iz = grid%iz1 + grid%idz*(kz-1)
            
            if (nec == 21) then
               write(20,fmt) ix,iy,iz,((grid%ecs(kx,ky,kz,i,j),j=i,6),i=1,6)
            else
               write(20,fmt) ix,iy,iz,grid%ecs(kx,ky,kz,:,:)
            endif
         enddo
      enddo
   enddo
   
   close(20)
   
   end subroutine EC_grid_write
!-------------------------------------------------------------------------------

!===============================================================================
   subroutine EC_grid_load_bin(fname,grid,quiet)
!===============================================================================
!  Load the ECgrid type from a binary file

   implicit none
   
   character(len=*),intent(in) :: fname
   type(ECgrid),intent(inout)  :: grid
   logical,intent(in),optional :: quiet
   logical                     :: silent
   integer :: iostatus,ix,iy,iz,i,j,k
   
   silent=.true.
   if (present(quiet)) silent=quiet
   
!  open file   
   open(20,file=fname,iostat=iostatus,form='unformatted',status='old')
   if (iostatus /= 0) then
      write(lu_stderr,'(a,a)') &
                     'EC_grid_load_bin: cannot open file for reading: ',fname
      stop
   endif
   
!  Read header line
   read(20) grid%npts, grid%nx, grid%ny, grid%nz, grid%ix1, grid%ix2, &
             grid%iy1, grid%iy2, grid%iz1, grid%iz2, grid%idx, grid%idy, grid%idz
   
!  Deallocate memory if reusing existing grid object and dimensions are different
   if (allocated(grid%ecs)) then
      if (size(grid%ecs,1) /= grid%nx .or. &
          size(grid%ecs,2) /= grid%ny .or. &
          size(grid%ecs,3) /= grid%nz)  deallocate(grid%ecs)
   endif
   if (allocated(grid%x)) then
      if (size(grid%x) /= grid%nx) deallocate(grid%x)
   endif
   if (allocated(grid%y)) then
      if (size(grid%y) /= grid%ny) deallocate(grid%y)
   endif
   if (allocated(grid%z)) then
      if (size(grid%z) /= grid%nz) deallocate(grid%z)
   endif
   
!  Allocate memory and first check this grid isn't too big
   if (real(grid%nx)*real(grid%ny)*real(grid%nz)*real(r8) > MAX_GRID_SIZE_IN_BYTES) then
      write(lu_stderr,'(a,i0,a)') 'EC_grid_load_bin: Allocated size of grid is ',&
         real(grid%nx)*real(grid%ny)*real(grid%nz)*real(r8)/real(1024)**2,' MB: stopping.'
      stop
   endif
   if (.not.allocated(grid%ecs)) allocate(grid%ecs(grid%nx,grid%ny,grid%nz,6,6))
   if (.not.allocated(grid%x)) allocate(grid%x(grid%nx))
   if (.not.allocated(grid%y)) allocate(grid%y(grid%ny))
   if (.not.allocated(grid%z)) allocate(grid%z(grid%nz))
   
!  Load constants
   do iz=1,grid%nz
      do iy=1,grid%ny
         do ix=1,grid%nx
            read(20) ((grid%ecs(ix,iy,iz,i,j),j=i,6),i=1,6)
            do i=1,6; do j=i,6
               grid%ecs(ix,iy,iz,j,i) = grid%ecs(ix,iy,iz,i,j)
            enddo; enddo
         enddo
      enddo
   enddo
   
   close(20)
   
!  Fill in dimensions
   do k=0,grid%nx-1; grid % x(k+1) = real(grid%ix1) + real(grid%idx) * real(k); enddo
   do k=0,grid%ny-1; grid % y(k+1) = real(grid%iy1) + real(grid%idy) * real(k); enddo
   do k=0,grid%nz-1; grid % z(k+1) = real(grid%iz1) + real(grid%idz) * real(k); enddo

   if (.not.silent) then
     write(lu_stdout,'(a)')         '=================================='
     write(lu_stdout,'(a)')         'Dimensions of box are:'
     write(lu_stdout,'(a,i0,x,i0)') '  x:  ',grid%ix1,grid%ix2
     write(lu_stdout,'(a,i0,x,i0)') '  y:  ',grid%iy1,grid%iy2
     write(lu_stdout,'(a,i0,x,i0)') '  z:  ',grid%iz1,grid%iz2
     write(lu_stdout,'(a)')         '=================================='
     write(lu_stdout,'(3(a,i0))')   'Grid size is    ',grid%nx, ' x ',grid%ny, ' x ',grid%nz
     write(lu_stdout,'(3(a,i0))')   'Grid spacing is ',grid%idx,' x ',grid%idy,' x ',grid%idz
     write(lu_stdout,'(a,i0)')      'Total number of ECs is ',grid%npts
     write(lu_stdout,'(a/)')        '=================================='
   endif

   end subroutine EC_grid_load_bin
!------------------------------------------------------------------------------   

!===============================================================================
   subroutine EC_grid_write_bin(fname,grid)
!===============================================================================
!  Write the ECgrid type to a binary, direct-access file.
!  Format is:
!   [1] npts,nx,ny,nz,ix1,ix2,iy1,iy2,iz1,iz2,idx,idy,idz
!   [2+] c11,c12...c16,c22...c66 (loops over z slowest, x fastest)

   implicit none
   
   type(ECgrid),intent(in)  :: grid
   character(len=*),intent(in) :: fname
   integer  :: iostatus,ix,iy,iz,i,j
    
   open(20,file=fname,iostat=iostatus,form='unformatted')
   if (iostatus /= 0) then
      write(lu_stderr,'(a,a)') &
                     'EC_grid_write_bin: Cannot open file for writing: ',trim(fname)
      stop
   endif
   
!  Write the header line
   write(20) grid%npts, grid%nx, grid%ny, grid%nz, grid%ix1, grid%ix2, &
             grid%iy1, grid%iy2, grid%iz1, grid%iz2, grid%idx, grid%idy, grid%idz
             
!  Write the elastic constants
   do iz=1,grid%nz
      do iy=1,grid%ny
         do ix=1,grid%nx
            write(20) ((grid%ecs(ix,iy,iz,i,j),j=i,6),i=1,6)
         enddo
      enddo
   enddo
   
   close(20)
   
   end subroutine EC_grid_write_bin
!------------------------------------------------------------------------------
   
!==============================================================================
   subroutine EC_grid_dump(grid,necs)
!===============================================================================
!  Write the EC_grid type to stdout.
!  Assumes output is all 36 elastic constants, otherwise can specify 21.

   implicit none
   
   type(ECgrid),intent(in)       :: grid
   character(len=80)             :: fmt
   integer,optional,intent(in)   :: necs
   integer                       :: ix,iy,iz,kx,ky,kz,i,j,nec
   integer                       :: lenx,leny,lenz
   
!  Check for requested output format; default to all 36 constants to be written out.
   if (present(necs)) then
      if (necs /= 21 .and. necs /= 36) stop 'EC_grid_dump: necs must be 21 or 36'
      nec = necs
   else
      nec = 36
   endif
   
!  Determine output format; minimum length is 5
   lenx = 5   ;   leny = lenx   ;   lenz = lenx
   if (int(log10(real(grid%ix1))+2.) > lenx) lenx = int(log10(real(grid%ix1))+2.)
   if (int(log10(real(grid%ix2))+2.) > lenx) lenx = int(log10(real(grid%ix2))+2.)
   if (int(log10(real(grid%iy1))+2.) > leny) leny = int(log10(real(grid%iy1))+2.)
   if (int(log10(real(grid%iy2))+2.) > leny) leny = int(log10(real(grid%iy2))+2.)
   if (int(log10(real(grid%iz1))+2.) > lenz) lenz = int(log10(real(grid%iz1))+2.)
   if (int(log10(real(grid%iz2))+2.) > lenz) lenz = int(log10(real(grid%iz2))+2.)
   
   write(fmt,'(3(a4,i3),a1,i2,a)') '  (i',lenx,',x,i',leny,',x,i',lenz,',',nec,'e24.16)'  
   write(fmt,'(a,i2,a)') '(3(i0,x),',nec,'e24.16)'

   do kx=1,grid%nx
      ix = grid%ix1 + grid%idx*(kx-1)
      do ky=1,grid%ny
         iy = grid%iy1 + grid%idy*(ky-1)
         do kz=1,grid%nz
            iz = grid%iz1 + grid%idz*(kz-1)
            
            if (nec == 21) then
               write(*,fmt) ix,iy,iz,((grid%ecs(kx,ky,kz,i,j),j=i,6),i=1,6)
            else
               write(*,fmt) ix,iy,iz,grid%ecs(kx,ky,kz,:,:)
            endif
         enddo
      enddo
   enddo
   
   end subroutine EC_grid_dump
!------------------------------------------------------------------------------

!==============================================================================
   subroutine EC_grid_clone(in,out)
!  Create a new grid of ECs from another

   implicit none
   
     type(ECgrid),intent(in)   :: in
     type(ECgrid),intent(out)  :: out
    
     out % npts = in % npts
     out % nx = in % nx
     out % ny = in % ny
     out % nz = in % nz
     out % ix1 = in % ix1
     out % ix2 = in % ix2
     out % iy1 = in % iy1
     out % iy2 = in % iy2
     out % iz1 = in % iz1
     out % iz2 = in % iz2
     out % idx = in % idx
     out % idy = in % idy
     out % idz = in % idz
     
     allocate(out%ecs(out%nx, out%ny, out%nz, 6, 6))
     allocate(out%x(out%nx), out%y(out%ny), out%z(out%nz))
     
     out % ecs = in % ecs
     out % x = in % x
     out % y = in % y
     out % z = in % z
      
   end subroutine EC_grid_clone
!------------------------------------------------------------------------------

!==============================================================================
   subroutine EC_grid_delete(grid)
!==============================================================================
!  Zeros out all values and deallocates the space for the arrays
   
   implicit none
   
   type(ECgrid),intent(inout) :: grid
   
     grid % npts = 0
     grid % nx = 0
     grid % ny = 0
     grid % nz = 0
     grid % ix1 = 0
     grid % ix2 = 0
     grid % iy1 = 0
     grid % iy2 = 0
     grid % iz1 = 0
     grid % iz2 = 0
     grid % idx = 0
     grid % idy = 0
     grid % idz = 0  
     deallocate(grid%ecs)
     deallocate(grid%x)
     deallocate(grid%y)
     deallocate(grid%z)
   
   end subroutine EC_grid_delete
!------------------------------------------------------------------------------

!==============================================================================
   function EC_grid_check_necs(fname)
!==============================================================================
!  Determines the number of ECs per line in an EC_grid-type file

   implicit none
   
   character(len=*),intent(in) :: fname
   integer           :: EC_grid_check_necs
   real(rs)          :: x,y,z, test(21)
   real(rs),parameter:: test_comparison=1.d0
   integer           :: iostatus
   
!  Test for which way round the input file is
   open(20,file=fname,status='old',iostat=iostatus)
   if (iostatus /= 0) then
      write(0,'(a,a)') 'EC_grid_check_necs: can''t open file:',trim(fname)
      stop
   endif
  
!  Compare the values in the 4th column with that in the 10th (i.e. the second and
!  seventh columns of the ECs):
!  if n=36, they should be the same (have read c12 and c21).  Otherwise, n=21.  
!  If this is true by accident, the constants are in completely the wrong fromat anyway.
   test(1) = 0.
   do while (abs(test(1)) < test_comparison)
      read(20,*) x,y,z,test
   enddo
   if (abs(test(2) - test(7)) < test_comparison) then
      EC_grid_check_necs = 36
   else
      EC_grid_check_necs = 21
   endif
  
   close(20)
      
   end function EC_grid_check_necs
!-------------------------------------------------------------------------------

!===============================================================================
   subroutine EC_grid_check_type_bin(file,nver,quiet)
!===============================================================================
!  Check whether a file represents a binary (1) implicit integer, (2) explicit
!  integer or (3) explicit-real ECgrid type file.

   implicit none
   
   character(len=*),intent(in) :: file
   integer,intent(out)         :: nver
   logical,optional,intent(in) :: quiet
   logical                     :: silent
   integer                     :: iostatus
   
!  Check for verbosity level
   if (present(quiet)) then
      silent = quiet
   else
      silent = .true.
   endif
   
   open(20,file=file,iostat=iostatus,form='unformatted',status='old')
   if (iostatus /= 0) then
      write(0,'(a,a)') 'EC_grid_check_type_bin: can''t open file: ',trim(file)
      stop
   endif
   
!  Read first item in file header
   read(20,iostat=iostatus) nver
   if (iostatus /= 0) then
      write(0,'(a,a,a)') 'EC_grid_check_type_bin: can''t read first item in header of ',&
         'binary file ',trim(file)
      stop
   endif
   close(20)
   
!  Check version
   if (nver /= nver_i .and. nver /= nver_r) nver = nver_assumed_i
   
!  Write this out if we want
   if (.not.silent) then
      select case (nver)
         case(nver_assumed_i)
            write(lu_stdout,'(a,a,a)') 'EC_grid_check_type_bin: file "',trim(file),&
               ' is type "assumed-integer-coordinates"'
         case(nver_i)
            write(lu_stdout,'(a,a,a)') 'EC_grid_check_type_bin: file "',trim(file),&
               ' is type "explicit-integer-coordinates"'
         case(nver_r)
            write(lu_stdout,'(a,a,a)') 'EC_grid_check_type_bin: file "',trim(file),&
               ' is type "explicit-real-coordinates"'
      end select
   endif
   
   end subroutine EC_grid_check_type_bin
!-------------------------------------------------------------------------------

!===============================================================================
   subroutine EC_grid_to_gridr(gridr,grid,gridi)
!===============================================================================
!  Converts assumed- or explicit-integer dimensioned grids into ones with real
!  dimensions.  Supply either grid or gridi via optional arguments.

   implicit none
   
   type(ECgridr),intent(out) :: gridr
   type(ECgrid),intent(in),optional :: grid
   type(ECgridi),intent(in),optional :: gridi
   real(r8) :: x1,x2,y1,y2,z1,z2,dx,dy,dz
   
!  Check for correct options
   if (.not.present(grid).and..not.present(gridi)) then
      write(lu_stderr,'(a)') 'EC_grid_to_gridr: Supply either assumed-int or explicit-int input grids.'
      stop
   else if (present(grid).and.present(gridi)) then
      write(lu_stderr,'(a)') 'EC_grid_to_gridr: Only one of assumed-ont or explicit-int grids can be supplied.'
      stop
   endif
   
!  Get limits for new grid, create new one, and copy constants across
   if (present(grid)) then
      x1 = real(grid%ix1);  x2 = real(grid%ix2)
      y1 = real(grid%iy1);  y2 = real(grid%iy2)
      z1 = real(grid%iz1);  z2 = real(grid%iz2)
      dx = real(grid%idx);  dy = real(grid%idy);  dz = real(grid%idz)
      call EC_gridr_new(x1,x2,dx,y1,y2,dy,z1,z2,dz,gridr)
      gridr%ecs = grid%ecs
   else
      x1 = real(gridi%ix1);  x2 = real(gridi%ix2)
      y1 = real(gridi%iy1);  y2 = real(gridi%iy2)
      z1 = real(gridi%iz1);  z2 = real(gridi%iz2)
      dx = real(gridi%idx);  dy = real(gridi%idy);  dz = real(gridi%idz)
      call EC_gridr_new(x1,x2,dx,y1,y2,dy,z1,z2,dz,gridr)
      gridr%ecs = gridi%ecs
   endif
   
   end subroutine EC_grid_to_gridr
!-------------------------------------------------------------------------------


!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!                         EXPLICIT INTEGER COORDINATES
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!===============================================================================
   subroutine EC_gridi_new(ix1,ix2,idx,iy1,iy2,idy,iz1,iz2,idz,grid)
!===============================================================================
!  Create a new grid

   implicit none
   
   type(ECgrid),intent(out)  :: grid
   integer,intent(in)        :: ix1,ix2,idx,iy1,iy2,idy,iz1,iz2,idz
   integer                   :: k
   
   if ((ix1>=ix2).or.(iy1>=iy2).or.(iz1>=iz2).or.(idx<=0).or.(idy<=0).or.(idz<=0)) then
      write(lu_stderr,'(a)') 'EC_grid_new: problem with array limits or spacing.'
      write(lu_stderr,'(a)') '    ix1<ix2, iy1<iy2, iz1<iz1, (idx,idy,idz)>0'
      stop
   endif
   
   grid % ix1 = ix1  ;  grid % ix2 = ix2
   grid % iy1 = iy1  ;  grid % iy2 = iy2
   grid % iz1 = iz1  ;  grid % iz2 = iz2
   grid % idx = idx
   grid % idy = idy
   grid % idz = idz
   
!  Calcuate array dimensions
   grid % nx = (ix2 - ix1) / idx + 1
   grid % ny = (iy2 - iy1) / idy + 1
   grid % nz = (iz2 - iz1) / idz + 1
   grid % npts = grid%nx * grid%ny * grid%nz
   
   if (allocated(grid%ecs)) deallocate(grid%ecs)
   if (allocated(grid%x)) deallocate(grid%x)
   if (allocated(grid%y)) deallocate(grid%y)
   if (allocated(grid%z)) deallocate(grid%z)

!  Allocate memory   
   allocate(grid%ecs(grid%nx,grid%ny,grid%nz,6,6))
   allocate(grid%x(grid%nx), grid%y(grid%ny), grid%z(grid%nz))
   grid%ecs = 0.
!  Fill in dimensions
   do k=0,grid%nx-1;  grid % x(k+1) = real(ix1) + real(idx) * real(k);   enddo
   do k=0,grid%ny-1;  grid % y(k+1) = real(iy1) + real(idy) * real(k);   enddo
   do k=0,grid%nz-1;  grid % z(k+1) = real(iz1) + real(idz) * real(k);   enddo
   
   end subroutine EC_gridi_new
!-------------------------------------------------------------------------------

!==============================================================================
   subroutine EC_gridi_load(ecs_file,grid,necs,quiet)
!==============================================================================
!  Load a 3D grid of elastic constants into memory.
!  Elastic constants should be oriented according to the grid, where
!     x // 1   ,   y // 2   ,   z // 3, 
!  assuming a right-handed coordinate system.  Later tranformations can be performed.

   implicit none

   type(ECgridi),intent(out)    :: grid
   character(len=*),intent(in) :: ecs_file
   integer   :: ix,iy,iz,ix1,ix2,iy1,iy2,iz1,iz2,idx,idy,idz,&
                itempx,itempy,itempz,nx,ny,nz,npts,i,j,k,nec
   integer,optional,intent(in) :: necs
   logical,optional,intent(in) :: quiet
   logical                     :: silent
   real(rs)   :: ecs_in(6,6),multiplier
   integer    :: iostatus
   
!  Check optional arguments
   if (present(necs)) then
      if (necs /=21 .and. necs /= 36) stop 'EC_grid_load: nec must be 21 or 36'
      nec = necs
   else
      nec = 36
   endif
   
   silent = .false.
   if (present(quiet)) then
      silent = quiet
   endif
   
!  Open file
   open(20,file=ecs_file,status='old',iostat=iostatus)
   if (iostatus /= 0) then
      write(*,'(a,a)') 'Problem opening ECs file ',ecs_file
      stop
   endif  
   
   iostatus = 0
   ix1 = big_p_integer ; ix2 = big_n_integer
   iy1 = big_p_integer ; iy2 = big_n_integer
   iz1 = big_p_integer ; iz2 = big_n_integer
   npts = 0

!  Read limits of data and find dimensions
   do while (iostatus == 0)
      read(20,fmt=*,iostat=iostatus) ix,iy,iz
      if (iostatus /= 0) exit
      if (ix < ix1) ix1 = ix
      if (ix > ix2) ix2 = ix
      if (iy < iy1) iy1 = iy
      if (iy > iy2) iy2 = iy
      if (iz < iz1) iz1 = iz
      if (iz > iz2) iz2 = iz
      npts = npts + 1
   enddo
   
   if (.not.silent) then
     write(0,'(a)') '=================================='
     write(0,'(a)') 'Dimensions of box are:'
     write(0,'(a,2i5)') '  x:  ',ix1,ix2
     write(0,'(a,2i5)') '  y:  ',iy1,iy2
     write(0,'(a,2i5)') '  z:  ',iz1,iz2
     write(0,'(a)') '=================================='
   endif
   
!  Start again and work out the dimensions
   rewind(20)
!  Count the number of points with the same x and z as the first line.
   read(20,*) itempx,itempy,itempz
   ny = 1
   do i=2,npts
      read(20,fmt=*,iostat=iostatus) ix,iy,iz
      if (ix == itempx .and. iz == itempz) ny = ny + 1
   enddo
   
!  Count the number of x at this y and z
   rewind(20)
   read(20,*) itempx,itempy,itempz
   nx = 1
   do i=2,npts
      read(20,fmt=*,iostat=iostatus) ix,iy,iz
      if (iy == itempy .and. iz == itempz) nx = nx + 1
   enddo
   
!  Count the number of z at this x and y: just for now, to check!
   rewind(20)
   read(20,*) itempx,itempy,itempz
   nz = 1
   do i=2,npts
      read(20,fmt=*,iostat=iostatus) ix,iy,iz
      if (ix == itempx .and. iy == itempy) nz = nz + 1
   enddo
   
!  Sanity check
   if (nx*ny*nz /= npts) then
      write(*,'(a)') 'load_ecs_cart: Number of points in EC file does not match dimensions!'
      stop
   endif
   
   idx = (ix2 - ix1) / (nx-1)
   idy = (iy2 - iy1) / (ny-1)
   idz = (iz2 - iz1) / (nz-1)
   if (.not.silent) then
     write(0,'(3(a,i4))') 'Grid size is    ',nx,' x ',ny,' x ',nz
     write(0,'(3(a,i3))') 'Grid step size is ',idx,' x ',idy,' x ',idz
     write(0,'(a,i6)')    'Total number of ECs is ',npts
     write(0,'(a/)')      '=================================='
   endif
   
!  Allocate space for the ECs
   allocate(grid%ecs(nx,ny,nz,6,6),stat=iostatus)
   allocate(grid%x(nx), grid%y(ny), grid%z(nz))
   if (iostatus /= 0) stop 'EC_grid_load: error allocating space for ECs'
   
!  Read each line and fill the array
   rewind(20)
   do k=1,npts
      if (nec == 36) then
         read(20,*) itempx,itempy,itempz,ecs_in
      else if (nec == 21) then
         read(20,*) itempx,itempy,itempz,((ecs_in(i,j),j=i,6),i=1,6)
         do i=1,6
            do j=1,6
               ecs_in(j,i) = ecs_in(i,j)
            enddo
         enddo
      endif
!     For the first layer, determine whether we're in Pa or GPa; issue a warning regardless.
      if (k == 1) then
         if (ecs_in(1,1) > 1.d5) then
            multiplier = 1.d0
            if (.not.silent) &
            write(0,'(a)') '  read_ecs_cart: Input Cij file appears to be in Pa.'
         else
            multiplier = 1.d9
            if (.not.silent) &
            write(0,'(a)') '  read_ecs_cart: Input Cij file appears to be in GPa.'
         endif
      endif

!  Fill the array with the constants in Pa
      grid%ecs((itempx-ix1+idx)/idx,&
               (itempy-iy1+idy)/idy,&
               (itempz-iz1+idz)/idz,:,:) = ecs_in(:,:) * multiplier
      
   enddo
   
   close(20)

   grid % nx = nx  ;   grid % ny = ny  ;   grid % nz = nz
   grid % idx = idx  ;  grid % idy = idy  ;   grid % idz = idz
   grid % ix1 = ix1  ;  grid % ix2 = ix2
   grid % iy1 = iy1  ;  grid % iy2 = iy2
   grid % iz1 = iz1  ;  grid % iz2 = iz2
   
!  Fill in dimensions
   do k=0,nx-1;  grid % x(k+1) = real(ix1) + real(idx) * real(k);   enddo
   do k=0,ny-1;  grid % y(k+1) = real(iy1) + real(idy) * real(k);   enddo
   do k=0,nz-1;  grid % z(k+1) = real(iz1) + real(idz) * real(k);   enddo
   
   end subroutine EC_gridi_load
!------------------------------------------------------------------------------

!===============================================================================
   subroutine EC_gridi_inquire_bin(ecs_file,grid,quiet)
!===============================================================================
!  Get information about a binary grid file with real coordinates
   
   implicit none
   
   type(ECgridi),intent(inout)  :: grid
   character(len=*),intent(in) :: ecs_file
   logical,intent(in),optional :: quiet
   logical                     :: silent
   integer                     :: iostatus,nver
   
   silent = .true.
   if (present(quiet)) silent = quiet
   
   iostatus = 0
   open(20,file=trim(ecs_file),iostat=iostatus,form='unformatted',status='old')
   if (iostatus /= 0) then
      write(lu_stderr,'(a,a)') &
         'EC_grid_inquire_bin: cannot open file for reading: ', trim(ecs_file)
      stop
   endif
   
!  Read header line
   read(20,iostat=iostatus) nver, grid%npts, grid%nx, grid%ny, grid%nz, grid%ix1, grid%ix2, &
             grid%iy1, grid%iy2, grid%iz1, grid%iz2, grid%idx, grid%idy, grid%idz
   close(20)
   if (iostatus /= 0) then
      write(lu_stderr,'(a,a)') &
         'EC_gridr_inquire_bin: cannot read header line for binary file ',trim(ecs_file)
      stop
   endif
!  Check correct file version
   if (nver /= 1) then
      write(lu_stderr,'(a)') &
         'EC_gridi_inquire_bin: file appears to be incorrect format for explicit integer coordinates.'
      stop
   endif
   
!  If desired, print out information
   if (.not.silent) then
     write(lu_stdout,'(a)')         '=================================='
     write(lu_stdout,'(a)')         'Dimensions of box are:'
     write(lu_stdout,'(a,i0,x,i0)') '  x:  ',grid%ix1,grid%ix2
     write(lu_stdout,'(a,i0,x,i0)') '  y:  ',grid%iy1,grid%iy2
     write(lu_stdout,'(a,i0,x,i0)') '  z:  ',grid%iz1,grid%iz2
     write(lu_stdout,'(a)')         '=================================='
     write(lu_stdout,'(3(a,i0))') 'Grid size is    ',grid%nx, ' x ',grid%ny, ' x ',grid%nz
     write(lu_stdout,'(3(a,i0))') 'Grid spacing is ',grid%idx,' x ',grid%idy,' x ',grid%idz
     write(lu_stdout,'(a,i0)')    'Total number of ECs is ',grid%npts
     write(lu_stdout,'(a/)')      '=================================='
   endif
   
   end subroutine EC_gridi_inquire_bin
!-------------------------------------------------------------------------------

!==============================================================================
   subroutine EC_gridi_write(fname,grid,necs)
!===============================================================================
!  Write the EC_grid type to an ASCII text file.
!  Assumes output is all 36 elastic constants, otherwise can specify 21.

   implicit none
   
   type(ECgridi),intent(in)       :: grid
   character(len=*),intent(in)   :: fname
   character(len=80)             :: fmt
   integer,optional,intent(in)   :: necs
   integer                       :: iostatus,ix,iy,iz,kx,ky,kz,i,j,nec
   integer                       :: lenx,leny,lenz
   
!  Check for requested output format; default to all 36 constants to be written out.
   if (present(necs)) then
      if (necs /= 21 .and. necs /= 36) stop 'EC_grid_write: nec must be 21 or 36'
      nec = necs
   else
      nec = 36
   endif
   
   open(20,file=trim(fname),iostat=iostatus)
   if (iostatus /= 0) stop 'EC_grid_write: error opening file for writing'
   
!  Determine output format; minimum length is 5
   lenx = 5   ;   leny = lenx   ;   lenz = lenx
   if (int(log10(real(grid%ix1))+1.) > lenx) lenx = int(log10(real(grid%ix1))+1.)
   if (int(log10(real(grid%ix2))+1.) > lenx) lenx = int(log10(real(grid%ix2))+1.)
   if (int(log10(real(grid%iy1))+1.) > leny) leny = int(log10(real(grid%iy1))+1.)
   if (int(log10(real(grid%iy2))+1.) > leny) leny = int(log10(real(grid%iy2))+1.)
   if (int(log10(real(grid%iz1))+1.) > lenz) lenz = int(log10(real(grid%iz1))+1.)
   if (int(log10(real(grid%iz2))+1.) > lenz) lenz = int(log10(real(grid%iz2))+1.)
   
   write(fmt,'(3(a2,i3),a1,i2,a)') '(i',lenx,',i',leny,',i',lenz,',',nec,'e24.16)'  

   do kx=1,grid%nx
      ix = grid%ix1 + grid%idx*(kx-1)
      do ky=1,grid%ny
         iy = grid%iy1 + grid%idy*(ky-1)
         do kz=1,grid%nz
            iz = grid%iz1 + grid%idz*(kz-1)
            
            if (nec == 21) then
               write(20,fmt) ix,iy,iz,((grid%ecs(kx,ky,kz,i,j),j=i,6),i=1,6)
            else
               write(20,fmt) ix,iy,iz,grid%ecs(kx,ky,kz,:,:)
            endif
         enddo
      enddo
   enddo
   
   close(20)
   
   end subroutine EC_gridi_write
!------------------------------------------------------------------------------

!==============================================================================
   subroutine EC_gridi_load_bin(fname,grid,quiet)
!  Load the ECgrid type from a binary file

   implicit none
   
   character(len=*),intent(in) :: fname
   type(ECgridi),intent(out) :: grid
   logical,intent(in),optional :: quiet
   integer :: iostatus,ix,iy,iz,i,j,k,nver
   
   write(lu_stderr,'(a)') 'EC_grid_load_bin hasn''t been tested yet!'
   stop
   
   
   open(20,file=fname,iostat=iostatus,form='unformatted',status='old')
   if (iostatus /= 0) then
      write(lu_stderr,'(a,a)') &
                     'EC_grid_load_bin: cannot open file for reading: ',fname
      stop
   endif
   
   read(20) nver, grid%npts, grid%nx, grid%ny, grid%nz, grid%ix1, grid%ix2, &
             grid%iy1, grid%iy2, grid%iz1, grid%iz2, grid%idx, grid%idy, grid%idz
   
   ! Check for version of grid
   if (nver /= nver_i) then
      write(0,'(a)') &
         'EC_gridi_load_bin: Bad version number for integer-dimensioned ECgrid binary file.',&
         '  Check grid is of correct format.'
      stop
   endif

!  Deallocate memory if reusing existing grid object and dimensions are different
   if (allocated(grid%ecs)) then
      if (size(grid%ecs,1) /= grid%nx .or. &
          size(grid%ecs,2) /= grid%ny .or. &
          size(grid%ecs,3) /= grid%nz)  deallocate(grid%ecs)
   endif
   if (allocated(grid%x)) then
      if (size(grid%x) /= grid%nx) deallocate(grid%x)
   endif
   if (allocated(grid%y)) then
      if (size(grid%y) /= grid%ny) deallocate(grid%y)
   endif
   if (allocated(grid%z)) then
      if (size(grid%z) /= grid%nz) deallocate(grid%z)
   endif
   
!  Allocate memory
   if (.not.allocated(grid%ecs)) allocate(grid%ecs(grid%nx,grid%ny,grid%nz,6,6))
   if (.not.allocated(grid%x)) allocate(grid%x(grid%nx))
   if (.not.allocated(grid%y)) allocate(grid%y(grid%ny))
   if (.not.allocated(grid%z)) allocate(grid%z(grid%nz))

!  Read in constants
   do iz=1,grid%nz
      do iy=1,grid%ny
         do ix=1,grid%nx
            read(20) ((grid%ecs(ix,iy,iz,i,j),j=i,6),i=1,6)
            do i=1,6; do j=i,6
               grid%ecs(ix,iy,iz,j,i) = grid%ecs(ix,iy,iz,i,j)
            enddo; enddo
         enddo
      enddo
   enddo
   
   close(20)
   
!  Fill in dimensions
   do k=0,grid%nx-1; grid % x(k+1) = real(grid%ix1) + real(grid%idx) * real(k); enddo
   do k=0,grid%ny-1; grid % y(k+1) = real(grid%iy1) + real(grid%idy) * real(k); enddo
   do k=0,grid%nz-1; grid % z(k+1) = real(grid%iz1) + real(grid%idz) * real(k); enddo

   if (present(quiet)) then
   if (.not.quiet) then
     write(lu_stdout,'(a)')         '=================================='
     write(lu_stdout,'(a)')         'Dimensions of box are:'
     write(lu_stdout,'(a,i0,x,i0)') '  x:  ',grid%ix1,grid%ix2
     write(lu_stdout,'(a,i0,x,i0)') '  y:  ',grid%iy1,grid%iy2
     write(lu_stdout,'(a,i0,x,i0)') '  z:  ',grid%iz1,grid%iz2
     write(lu_stdout,'(a)')         '=================================='
     write(lu_stdout,'(3(a,i0))') 'Grid size is    ',grid%nx, ' x ',grid%ny, ' x ',grid%nz
     write(lu_stdout,'(3(a,i0))') 'Grid spacing is ',grid%idx,' x ',grid%idy,' x ',grid%idz
     write(lu_stdout,'(a,i0)')    'Total number of ECs is ',grid%npts
     write(lu_stdout,'(a/)')      '=================================='
   endif
   endif
   
   end subroutine EC_gridi_load_bin
!------------------------------------------------------------------------------   

!==============================================================================
   subroutine EC_gridi_write_bin(fname,grid)
!  Write the ECgrid type to a binary, direct-access file.
!  Format is:
!   [1] npts,nx,ny,nz,ix1,ix2,iy1,iy2,iz1,iz2,idx,idy,idz
!   [2+] c11,c12...c16,c21...c66 (loops over x slowest, z fastest)

   implicit none
   
   type(ECgridi),intent(in)  :: grid
   character(len=*),intent(in) :: fname
   integer  :: iostatus,ix,iy,iz
   
   write(lu_stderr,'(a)') 'EC_grid_write_bin hasn''t been tested yet!'
   stop
   
   open(20,file=fname,iostat=iostatus,form='unformatted')
   if (iostatus /= 0) then
      write(lu_stderr,'(a,a)') &
                     'EC_grid_write_bin: Cannot open file for reading: ',fname
      stop
   endif
   
!  Write the header line
   write(20) grid%npts, grid%nx, grid%ny, grid%nz, grid%ix1, grid%ix2, &
             grid%iy1, grid%iy2, grid%iz1, grid%iz2, grid%idx, grid%idy, grid%idz
             
!  Write the elastic constants
   do ix=1,grid%nx
      do iy=1,grid%ny
         do iz=1,grid%nz
            write(20) grid%ecs(ix,iy,iz,:,:)
         enddo
      enddo
   enddo
   
   close(20)
   
   end subroutine EC_gridi_write_bin
!------------------------------------------------------------------------------
   
!==============================================================================
   subroutine EC_gridi_dump(grid,necs)
!===============================================================================
!  Write the EC_grid type to stdout.
!  Assumes output is all 36 elastic constants, otherwise can specify 21.

   implicit none
   
   type(ECgridi),intent(in)      :: grid
   character(len=80)             :: fmt
   integer,optional,intent(in)   :: necs
   integer                       :: ix,iy,iz,kx,ky,kz,i,j,nec
   integer                       :: lenx,leny,lenz
   
!  Check for requested output format; default to all 36 constants to be written out.
   if (present(necs)) then
      if (necs /= 21 .and. necs /= 36) stop 'EC_grid_write: necs must be 21 or 36'
      nec = necs
   else
      nec = 36
   endif
   
!  Determine output format; minimum length is 5
   lenx = 5   ;   leny = lenx   ;   lenz = lenx
   if (int(log10(real(grid%ix1))+1.) > lenx) lenx = int(log10(real(grid%ix1))+1.)
   if (int(log10(real(grid%ix2))+1.) > lenx) lenx = int(log10(real(grid%ix2))+1.)
   if (int(log10(real(grid%iy1))+1.) > leny) leny = int(log10(real(grid%iy1))+1.)
   if (int(log10(real(grid%iy2))+1.) > leny) leny = int(log10(real(grid%iy2))+1.)
   if (int(log10(real(grid%iz1))+1.) > lenz) lenz = int(log10(real(grid%iz1))+1.)
   if (int(log10(real(grid%iz2))+1.) > lenz) lenz = int(log10(real(grid%iz2))+1.)
   
   write(fmt,'(3(a2,i3),a1,i2,a)') '(i',lenx,',i',leny,',i',lenz,',',nec,'e24.16)'  

   do kx=1,grid%nx
      ix = grid%ix1 + grid%idx*(kx-1)
      do ky=1,grid%ny
         iy = grid%iy1 + grid%idy*(ky-1)
         do kz=1,grid%nz
            iz = grid%iz1 + grid%idz*(kz-1)
            
            if (nec == 21) then
               write(*,fmt) ix,iy,iz,((grid%ecs(kx,ky,kz,i,j),j=i,6),i=1,6)
            else
               write(*,fmt) ix,iy,iz,grid%ecs(kx,ky,kz,:,:)
            endif
         enddo
      enddo
   enddo
   
   end subroutine EC_gridi_dump
!------------------------------------------------------------------------------

!==============================================================================
   subroutine EC_gridi_clone(in,out)
!  Create a new grid of ECs from another

   implicit none
   
     type(ECgridi),intent(in)   :: in
     type(ECgridi),intent(out)  :: out
    
     out % npts = in % npts
     out % nx = in % nx
     out % ny = in % ny
     out % nz = in % nz
     out % ix1 = in % ix1
     out % ix2 = in % ix2
     out % iy1 = in % iy1
     out % iy2 = in % iy2
     out % iz1 = in % iz1
     out % iz2 = in % iz2
     out % idx = in % idx
     out % idy = in % idy
     out % idz = in % idz
     
     allocate(out%ecs(out%nx, out%ny, out%nz, 6, 6))
     allocate(out%x(out%nx), out%y(out%ny), out%z(out%nz))
     
     out % ecs = in % ecs
     out % x = in % x
     out % y = in % y
     out % z = in % z
      
   end subroutine EC_gridi_clone
!------------------------------------------------------------------------------

!==============================================================================
   subroutine EC_gridi_delete(grid)
!==============================================================================
!  Zeros out all values and deallocates the space for the arrays
   
   implicit none
   
   type(ECgridi),intent(inout) :: grid
   
     grid % npts = 0
     grid % nx = 0
     grid % ny = 0
     grid % nz = 0
     grid % ix1 = 0
     grid % ix2 = 0
     grid % iy1 = 0
     grid % iy2 = 0
     grid % iz1 = 0
     grid % iz2 = 0
     grid % idx = 0
     grid % idy = 0
     grid % idz = 0  
     deallocate(grid%ecs)
     deallocate(grid%x)
     deallocate(grid%y)
     deallocate(grid%z)
   
   end subroutine EC_gridi_delete
!------------------------------------------------------------------------------
   

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!                         EXPLICIT REAL COORDINATES
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!===============================================================================
   subroutine EC_gridr_new(x1,x2,dx,y1,y2,dy,z1,z2,dz,grid)
!===============================================================================
!  Create a new grid with real coordinates.
!TODO: Make arguments optional, so one can specify either coordinate bounds
!      and number of points in each dimnension, or (as now), the grid spacing.

   implicit none
   
   type(ECgridr),intent(out)  :: grid
   real(r8),intent(in)        :: x1,x2,dx,y1,y2,dy,z1,z2,dz
   integer                    :: k
   
   if ((x1>=x2).or.(y1>=y2).or.(z1>=z2).or.(dx<=0.).or.(dy<=0.).or.(dz<=0.)) then
      write(lu_stderr,'(a)') 'EC_gridr_new: problem with array limits or spacing.'
      write(lu_stderr,'(a)') '    x1<x2, y1<y2, z1<z1, (dx,dy,dz)>0'
      stop
   endif
   
   grid % x1 = x1  ;  grid % x2 = x2
   grid % y1 = y1  ;  grid % y2 = y2
   grid % z1 = z1  ;  grid % z2 = z2
   grid % dx = dx
   grid % dy = dy
   grid % dz = dz
   
!  Calcuate array dimensions
   grid % nx = nint((x2 - x1) / dx) + 1
   grid % ny = nint((y2 - y1) / dy) + 1
   grid % nz = nint((z2 - z1) / dz) + 1
   grid % npts = grid%nx * grid%ny * grid%nz
   
   if (allocated(grid%ecs)) deallocate(grid%ecs)
   if (allocated(grid%x)) deallocate(grid%x)
   if (allocated(grid%y)) deallocate(grid%y)
   if (allocated(grid%z)) deallocate(grid%z)

!  Allocate memory
   allocate(grid%ecs(grid%nx,grid%ny,grid%nz,6,6))
   allocate(grid%x(grid%nx), grid%y(grid%ny), grid%z(grid%nz))
   
!  Fill in dimensions
   grid%ecs = 0.
   do k=0,grid%nx-1;  grid % x(k+1) = x1 + dx * real(k);   enddo
   do k=0,grid%ny-1;  grid % y(k+1) = y1 + dy * real(k);   enddo
   do k=0,grid%nz-1;  grid % z(k+1) = z1 + dz * real(k);   enddo
   
   end subroutine EC_gridr_new
!-------------------------------------------------------------------------------

!!==============================================================================
!   subroutine EC_gridr_load(ecs_file,grid,necs,quiet)
!!==============================================================================
!!  Load a 3D grid of elastic constants into memory.
!!  Elastic constants should be oriented according to the grid, where
!!     x // 1   ,   y // 2   ,   z // 3, 
!!  assuming a right-handed coordinate system.  Later tranformations can be performed.
!
!   implicit none
!
!   type(ECgridr),intent(out)    :: grid
!   character(len=*),intent(in) :: ecs_file
!   character(len=80)           :: fmt
!   integer   :: ix,iy,iz,ix1,ix2,iy1,iy2,iz1,iz2,idx,idy,idz,&
!                itempx,itempy,itempz,nx,ny,nz,npts,i,j,k,nec
!   integer,optional,intent(in) :: necs
!   logical,optional,intent(in) :: quiet
!   logical                     :: silent
!   real(rs)   :: ecs_in(6,6),multiplier
!   real(rs)   :: x,y,z,  dx,dy,dz
!   integer    :: iostatus
!   
!   write(*,*) 'EC_gridr_load has not been tested.'; stop
!   
!!  Check optional arguments
!   if (present(necs)) then
!      if (necs /=21 .and. necs /= 36) stop 'EC_grid_load: nec must be 21 or 36'
!      nec = necs
!   else
!      nec = 36
!   endif
!   
!   silent = .false.
!   if (present(quiet)) then
!      silent = quiet
!   endif
!   
!!  Open file
!   open(20,file=ecs_file,status='old',iostat=iostatus)
!   if (iostatus /= 0) then
!      write(*,'(a,a)') 'Problem opening ECs file ',ecs_file
!      stop
!   endif  
!   
!   iostatus = 0
!   ix1 = big_p_integer ; ix2 = big_n_integer
!   iy1 = big_p_integer ; iy2 = big_n_integer
!   iz1 = big_p_integer ; iz2 = big_n_integer
!   npts = 0
!
!!  Read limits of data and find dimensions
!   do while (iostatus == 0)
!      read(20,fmt=*,iostat=iostatus) x,y,z
!      if (iostatus /= 0) exit
!      if (ix < ix1) ix1 = ix
!      if (ix > ix2) ix2 = ix
!      if (iy < iy1) iy1 = iy
!      if (iy > iy2) iy2 = iy
!      if (iz < iz1) iz1 = iz
!      if (iz > iz2) iz2 = iz
!      npts = npts + 1
!   enddo
!   
!   if (.not.silent) then
!     write(0,'(a)') '=================================='
!     write(0,'(a)') 'Dimensions of box are:'
!     write(0,'(a,2i5)') '  x:  ',ix1,ix2
!     write(0,'(a,2i5)') '  y:  ',iy1,iy2
!     write(0,'(a,2i5)') '  z:  ',iz1,iz2
!     write(0,'(a)') '=================================='
!   endif
!   
!!  Start again and work out the dimensions
!   rewind(20)
!!  Count the number of points with the same x and z as the first line.
!   read(20,*) itempx,itempy,itempz
!   ny = 1
!   do i=2,npts
!      read(20,fmt=*,iostat=iostatus) ix,iy,iz
!      if (ix == itempx .and. iz == itempz) ny = ny + 1
!   enddo
!   
!!  Count the number of x at this y and z
!   rewind(20)
!   read(20,*) itempx,itempy,itempz
!   nx = 1
!   do i=2,npts
!      read(20,fmt=*,iostat=iostatus) ix,iy,iz
!      if (iy == itempy .and. iz == itempz) nx = nx + 1
!   enddo
!   
!!  Count the number of z at this x and y: just for now, to check!
!   rewind(20)
!   read(20,*) itempx,itempy,itempz
!   nz = 1
!   do i=2,npts
!      read(20,fmt=*,iostat=iostatus) ix,iy,iz
!      if (ix == itempx .and. iy == itempy) nz = nz + 1
!   enddo
!   
!!  Sanity check
!   if (nx*ny*nz /= npts) then
!      write(*,'(a)') 'load_ecs_cart: Number of points in EC file does not match dimensions!'
!      stop
!   endif
!   
!   idx = (ix2 - ix1) / (nx-1)
!   idy = (iy2 - iy1) / (ny-1)
!   idz = (iz2 - iz1) / (nz-1)
!   if (.not.silent) then
!     write(0,'(3(a,i4))') 'Grid size is    ',nx,' x ',ny,' x ',nz
!     write(0,'(3(a,i3))') 'Grid step size is ',idx,' x ',idy,' x ',idz
!     write(0,'(a,i6)')    'Total number of ECs is ',npts
!     write(0,'(a/)')      '=================================='
!   endif
!   
!!  Allocate space for the ECs
!   allocate(grid%ecs(nx,ny,nz,6,6),stat=iostatus)
!   allocate(grid%x(nx), grid%y(ny), grid%z(nz))
!   if (iostatus /= 0) stop 'EC_grid_load: error allocating space for ECs'
!   
!!  Read each line and fill the array
!   rewind(20)
!   do k=1,npts
!      if (nec == 36) then
!         read(20,*) itempx,itempy,itempz,ecs_in
!      else if (nec == 21) then
!         read(20,*) itempx,itempy,itempz,((ecs_in(i,j),j=i,6),i=1,6)
!         do i=1,6
!            do j=1,6
!               ecs_in(j,i) = ecs_in(i,j)
!            enddo
!         enddo
!      endif
!!     For the first layer, determine whether we're in Pa or GPa; issue a warning regardless.
!      if (k == 1) then
!         if (ecs_in(1,1) > 1.d5) then
!            multiplier = 1.d0
!            if (.not.silent) &
!            write(0,'(a)') '  read_ecs_cart: Input Cij file appears to be in Pa.'
!         else
!            multiplier = 1.d9
!            if (.not.silent) &
!            write(0,'(a)') '  read_ecs_cart: Input Cij file appears to be in GPa.'
!         endif
!      endif
!
!!  Fill the array with the constants in Pa
!      grid%ecs((itempx-ix1+idx)/idx,&
!               (itempy-iy1+idy)/idy,&
!               (itempz-iz1+idz)/idz,:,:) = ecs_in(:,:) * multiplier
!      
!   enddo
!   
!   close(20)
!
!   grid % nx = nx  ;   grid % ny = ny  ;   grid % nz = nz
!   grid % idx = idx  ;  grid % idy = idy  ;   grid % idz = idz
!   grid % ix1 = ix1  ;  grid % ix2 = ix2
!   grid % iy1 = iy1  ;  grid % iy2 = iy2
!   grid % iz1 = iz1  ;  grid % iz2 = iz2
!   
!!  Fill in dimensions
!   do k=0,nx-1;  grid % x(k+1) = real(ix1) + real(idx) * real(k);   enddo
!   do k=0,ny-1;  grid % y(k+1) = real(iy1) + real(idy) * real(k);   enddo
!   do k=0,nz-1;  grid % z(k+1) = real(iz1) + real(idz) * real(k);   enddo
!   
!   return
!   
!   end subroutine EC_gridr_load
!!------------------------------------------------------------------------------
!
!!==============================================================================
!   subroutine EC_gridr_write(fname,grid,necs)
!!===============================================================================
!!  Write the EC_grid type to an ASCII text file.
!!  Assumes output is all 36 elastic constants, otherwise can specify 21.
!
!   implicit none
!   
!   type(ECgridr),intent(in)       :: grid
!   character(len=*),intent(in)   :: fname
!   character(len=80)             :: fmt
!   integer,optional,intent(in)   :: necs
!   integer                       :: iostatus,ix,iy,iz,kx,ky,kz,i,j,nec
!   integer                       :: lenx,leny,lenz
!
!   write(*,*) 'EC_gridr_write has not been tested.'; stop
!   
!!  Check for requested output format; default to all 36 constants to be written out.
!   if (present(necs)) then
!      if (necs /= 21 .and. necs /= 36) stop 'EC_grid_write: nec must be 21 or 36'
!      nec = necs
!   else
!      nec = 36
!   endif
!   
!   open(20,file=trim(fname),iostat=iostatus)
!   if (iostatus /= 0) stop 'EC_grid_write: error opening file for writing'
!   
!!  Determine output format; minimum length is 5
!   lenx = 5   ;   leny = lenx   ;   lenz = lenx
!   if (int(log10(real(grid%ix1))+1.) > lenx) lenx = int(log10(real(grid%ix1))+1.)
!   if (int(log10(real(grid%ix2))+1.) > lenx) lenx = int(log10(real(grid%ix2))+1.)
!   if (int(log10(real(grid%iy1))+1.) > leny) leny = int(log10(real(grid%iy1))+1.)
!   if (int(log10(real(grid%iy2))+1.) > leny) leny = int(log10(real(grid%iy2))+1.)
!   if (int(log10(real(grid%iz1))+1.) > lenz) lenz = int(log10(real(grid%iz1))+1.)
!   if (int(log10(real(grid%iz2))+1.) > lenz) lenz = int(log10(real(grid%iz2))+1.)
!   
!   write(fmt,'(3(a2,i3),a1,i2,a)') '(i',lenx,',i',leny,',i',lenz,',',nec,'e24.16)'  
!
!   do kx=1,grid%nx
!      ix = grid%ix1 + grid%idx*(kx-1)
!      do ky=1,grid%ny
!         iy = grid%iy1 + grid%idy*(ky-1)
!         do kz=1,grid%nz
!            iz = grid%iz1 + grid%idz*(kz-1)
!            
!            if (nec == 21) then
!               write(20,fmt) ix,iy,iz,((grid%ecs(kx,ky,kz,i,j),j=i,6),i=1,6)
!            else
!               write(20,fmt) ix,iy,iz,grid%ecs(kx,ky,kz,:,:)
!            endif
!         enddo
!      enddo
!   enddo
!   
!   close(20)
!   
!   return
!   
!   end subroutine EC_gridr_write
!!------------------------------------------------------------------------------

!===============================================================================
   subroutine EC_gridr_load_bin(fname,grid,quiet)
!===============================================================================
!  Load the ECgridr type from a binary file
!  Format is:
!   [1] nver,npts,nx,ny,nz,x1,x2,y1,y2,z1,z2,dx,dy,dz
!   [2+] c11,c12,...,c16,c22,c23,...,c66 (loops over x fastest, z slowest)

   implicit none
   
   character(len=*),intent(in) :: fname
   type(ECgridr),intent(out) :: grid
   logical,intent(in),optional :: quiet
   integer :: iostatus,ix,iy,iz,i,j,k,nver
   
   open(20,file=fname,iostat=iostatus,form='unformatted',status='old')
   if (iostatus /= 0) then
      write(lu_stderr,'(a,a)') &
                     'EC_gridr_load_bin: cannot open file for reading: ',fname
      stop
   endif
   
   read(20) nver,grid%npts, grid%nx, grid%ny, grid%nz, grid%x1, grid%x2, &
             grid%y1, grid%y2, grid%z1, grid%z2, grid%dx, grid%dy, grid%dz
   
   ! Check for version of grid
   if (nver /= nver_r) then
      write(0,'(a)') &
         'EC_gridr_load_bin: Bad version number for real-dimensioned ECgrid binary file.',&
         '  Check grid is of correct format.'
      stop
   endif
   
!  Deallocate memory if reusing existing grid object and dimensions are different
   if (allocated(grid%ecs)) then
      if (size(grid%ecs,1) /= grid%nx .or. &
          size(grid%ecs,2) /= grid%ny .or. &
          size(grid%ecs,3) /= grid%nz)  deallocate(grid%ecs)
   endif
   if (allocated(grid%x)) then
      if (size(grid%x) /= grid%nx) deallocate(grid%x)
   endif
   if (allocated(grid%y)) then
      if (size(grid%y) /= grid%ny) deallocate(grid%y)
   endif
   if (allocated(grid%z)) then
      if (size(grid%z) /= grid%nz) deallocate(grid%z)
   endif
   
!  Allocate memory
   if (.not.allocated(grid%ecs)) allocate(grid%ecs(grid%nx,grid%ny,grid%nz,6,6))
   if (.not.allocated(grid%x)) allocate(grid%x(grid%nx))
   if (.not.allocated(grid%y)) allocate(grid%y(grid%ny))
   if (.not.allocated(grid%z)) allocate(grid%z(grid%nz))

   ! Read in constants, looping over z slowest
   k=0
   do iz=1,grid%nz
      do iy=1,grid%ny
         do ix=1,grid%nx
            read(20) ((grid%ecs(ix,iy,iz,i,j),j=i,6),i=1,6)
            do i=1,6; do j=i,6
               grid%ecs(ix,iy,iz,j,i) = grid%ecs(ix,iy,iz,i,j)
            enddo; enddo
         enddo
      enddo
   enddo
   
   close(20)
   
!  Fill in dimensions
   do k=0,grid%nx-1; grid % x(k+1) = grid%x1 + grid%dx * real(k); enddo
   do k=0,grid%ny-1; grid % y(k+1) = grid%y1 + grid%dy * real(k); enddo
   do k=0,grid%nz-1; grid % z(k+1) = grid%z1 + grid%dz * real(k); enddo

   if (present(quiet)) then
   if (.not.quiet) then
     write(lu_stdout,'(a)')         '=================================='
     write(lu_stdout,'(a)')         'Dimensions of box are:'
     write(lu_stdout,'(a,es13.6,x,es13.6)') '  x:  ',grid%x1,grid%x2
     write(lu_stdout,'(a,es13.6,x,es13.6)') '  y:  ',grid%y1,grid%y2
     write(lu_stdout,'(a,es13.6,x,es13.6)') '  z:  ',grid%z1,grid%z2
     write(lu_stdout,'(a)')         '=================================='
     write(lu_stdout,'(3(a,i0))') 'Grid size is    ',grid%nx, ' x ',grid%ny, ' x ',grid%nz
     write(lu_stdout,'(3(a,es13.6))') 'Grid spacing is ',grid%dx,' x ',grid%dy,' x ',grid%dz
     write(lu_stdout,'(a,i0)')    'Total number of ECs is ',grid%npts
     write(lu_stdout,'(a/)')      '=================================='
   endif
   endif
   
   end subroutine EC_gridr_load_bin
!------------------------------------------------------------------------------   

!===============================================================================
   subroutine EC_gridr_inquire_bin(ecs_file,grid,quiet)
!===============================================================================
!  Get information about a binary grid file with real coordinates
   
   implicit none
   
   type(ECgridr),intent(inout)  :: grid
   character(len=*),intent(in) :: ecs_file
   logical,intent(in),optional :: quiet
   logical                     :: silent
   integer                     :: iostatus,nver
   
   silent = .true.
   if (present(quiet)) silent = quiet
   
   iostatus = 0
   open(20,file=trim(ecs_file),iostat=iostatus,form='unformatted',status='old')
   if (iostatus /= 0) then
      write(lu_stderr,'(a,a)') &
         'EC_grid_inquire_bin: cannot open file for reading: ', trim(ecs_file)
      stop
   endif
   
!  Read header line
   read(20,iostat=iostatus) nver, grid%npts, grid%nx, grid%ny, grid%nz, grid%x1, grid%x2, &
             grid%y1, grid%y2, grid%z1, grid%z2, grid%dx, grid%dy, grid%dz
   close(20)
   if (iostatus /= 0) then
      write(lu_stderr,'(a,a)') &
         'EC_gridr_inquire_bin: cannot read header line for binary file ',trim(ecs_file)
      stop
   endif
!  Check correct file version
   if (nver /= 2) then
      write(lu_stderr,'(a)') &
         'EC_gridr_inquire_bin: file appears to be incorrect format for explicit real coordinates.'
      stop
   endif
   
! If desired, print out information
   if (.not.silent) then
     write(lu_stdout,'(a)')         '=================================='
     write(lu_stdout,'(a)')         'Dimensions of box are:'
     write(lu_stdout,'(a,es13.6,x,es13.6)') '  x:  ',grid%x1,grid%x2
     write(lu_stdout,'(a,es13.6,x,es13.6)') '  y:  ',grid%y1,grid%y2
     write(lu_stdout,'(a,es13.6,x,es13.6)') '  z:  ',grid%z1,grid%z2
     write(lu_stdout,'(a)')         '=================================='
     write(lu_stdout,'(3(a,i0))') 'Grid size is    ',grid%nx, ' x ',grid%ny, ' x ',grid%nz
     write(lu_stdout,'(3(a,es13.6))') 'Grid spacing is ',grid%dx,' x ',grid%dy,' x ',grid%dz
     write(lu_stdout,'(a,i0)')    'Total number of ECs is ',grid%npts
     write(lu_stdout,'(a/)')      '=================================='
   endif
   
   end subroutine EC_gridr_inquire_bin
!-------------------------------------------------------------------------------

!==============================================================================
   subroutine EC_gridr_write_bin(fname,grid)
!  Write the ECgridr type to a binary, direct-access file.
!  Format is:
!   [1] npts,nx,ny,nz,x1,x2,y1,y2,z1,z2,dx,dy,dz
!   [2+] c11,c12...c16,c21...c66 (loops over x slowest, z fastest)

   implicit none
   
   type(ECgridr),intent(in)  :: grid
   character(len=*),intent(in) :: fname
   integer  :: iostatus,ix,iy,iz,i,j
   
!   write(lu_stderr,'(a)') 'EC_gridr_write_bin hasn''t been tested yet!'
!   stop
   
   open(20,file=fname,iostat=iostatus,form='unformatted')
   if (iostatus /= 0) then
      write(lu_stderr,'(a,a)') &
                     'EC_grid_write_bin: Cannot open file for reading: ',fname
      stop
   endif
   
!  Write the header line
   write(20) nver_r, grid%npts, grid%nx, grid%ny, grid%nz, grid%x1, grid%x2, &
             grid%y1, grid%y2, grid%z1, grid%z2, grid%dx, grid%dy, grid%dz
             
!  Write the elastic constants
   do iz=1,grid%nz
      do iy=1,grid%ny
         do ix=1,grid%nx
            write(20) ((grid%ecs(ix,iy,iz,i,j),j=i,6),i=1,6)
         enddo
      enddo
   enddo
   
   close(20)
   
   end subroutine EC_gridr_write_bin
!------------------------------------------------------------------------------
   
!==============================================================================
   subroutine EC_gridr_dump(grid,necs)
!===============================================================================
!  Write the EC_grid type to stdout.
!  Assumes output is all 36 elastic constants, otherwise can specify 21.

   implicit none
   
   type(ECgridr),intent(in)      :: grid
   character(len=80)             :: fmt
   integer,optional,intent(in)   :: necs
   real(r8)                      :: x,y,z
   integer                       :: kx,ky,kz,i,j,nec
   
   write(*,*) 'EC_gridr_load has not been tested.'; stop

!  Check for requested output format; default to all 36 constants to be written out.
   if (present(necs)) then
      if (necs /= 21 .and. necs /= 36) stop 'EC_grid_write: necs must be 21 or 36'
      nec = necs
   else
      nec = 36
   endif
   
   write(fmt,'(a,i2,a)') '(3e18.10,',nec,'e24.16)'  

   do kx=1,grid%nx
      x = grid%x(kx)
      do ky=1,grid%ny
         y = grid%y(ky)
         do kz=1,grid%nz
            z = grid%z(kz)
            
            if (nec == 21) then
               write(*,fmt) x,y,z,((grid%ecs(kx,ky,kz,i,j),j=i,6),i=1,6)
            else
               write(*,fmt) x,y,z,grid%ecs(kx,ky,kz,:,:)
            endif
         enddo
      enddo
   enddo
   
   end subroutine EC_gridr_dump
!------------------------------------------------------------------------------

!==============================================================================
   subroutine EC_gridr_clone(in,out)
!===============================================================================
!  Create a new grid of ECs from another

   implicit none
   
     type(ECgridr),intent(in)   :: in
     type(ECgridr),intent(out)  :: out
    
     out % npts = in % npts
     out % nx = in % nx
     out % ny = in % ny
     out % nz = in % nz
     out % x1 = in % x1
     out % x2 = in % x2
     out % y1 = in % y1
     out % y2 = in % y2
     out % z1 = in % z1
     out % z2 = in % z2
     out % dx = in % dx
     out % dy = in % dy
     out % dz = in % dz
     
     allocate(out%ecs(out%nx, out%ny, out%nz, 6, 6))
     allocate(out%x(out%nx), out%y(out%ny), out%z(out%nz))
     
     out % ecs = in % ecs
     out % x = in % x
     out % y = in % y
     out % z = in % z
      
   end subroutine EC_gridr_clone
!------------------------------------------------------------------------------

!==============================================================================
   subroutine EC_gridr_delete(grid)
!==============================================================================
!  Zeros out all values and deallocates the space for the arrays
   
   implicit none
   
   type(ECgridr),intent(inout) :: grid
   
     grid % npts = 0
     grid % nx = 0
     grid % ny = 0
     grid % nz = 0
     grid % x1 = 0.
     grid % x2 = 0.
     grid % y1 = 0.
     grid % y2 = 0.
     grid % z1 = 0.
     grid % z2 = 0.
     grid % dx = 0.
     grid % dy = 0.
     grid % dz = 0.
     deallocate(grid%ecs)
     deallocate(grid%x)
     deallocate(grid%y)
     deallocate(grid%z)
   
   end subroutine EC_gridr_delete
!------------------------------------------------------------------------------
  

!______________________________________________________________________________
end module EC_grid
!==============================================================================