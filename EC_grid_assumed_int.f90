!==============================================================================
module EC_grid
!==============================================================================
!  Module containing routines for manipulating regularly-spaced grids of elastic
!  constants.  Defines a type for them and can do various things such as load
!  these structures into memory.
!
!  AJN 2011/02
!
!  This compatibility version always assumes that the grid coordinates are
!  in integers, which simplifies reading and writing the grid.

   implicit none
   
!  Declare functions
   public :: EC_grid_load
   public :: EC_grid_write
!   public :: EC_grid_load_bin
!   public :: EC_grid_write_bin
   public :: EC_grid_clone
   public :: EC_grid_dump
   public :: EC_grid_delete
   
   
!  ** size constants
   integer, parameter, private :: i4 = selected_int_kind(9) ; ! long int
   integer, parameter, private :: r4 = selected_real_kind(6,37) ; ! SP
   integer, parameter, private :: r8 = selected_real_kind(15,307) ; ! DP

!  ** precision selector
   integer, parameter, private :: rs = r8
   
!  ** maths constants and other useful things
   real(r8), parameter, private :: pi = 3.141592653589793238462643d0
   integer, parameter, private :: big_p_integer = 2147483647
   integer, parameter, private :: big_n_integer = -2147483647
   
!  ** IO constants
   integer, parameter, private :: lu_stdin  = 5
   integer, parameter, private :: lu_stdout = 6
   integer, parameter, private :: lu_stderr = 0
   
!------------------------------------------------------------------------------
!  Define the ECgrid type: assumed integer coordinates
   type :: ECgrid
      integer(i4) :: npts,nx,ny,nz,ix1,ix2,iy1,iy2,iz1,iz2,idx,idy,idz
      real(r8),allocatable ::  ecs(:,:,:,:,:)
      real(r8),allocatable ::  x(:),y(:),z(:)
   end type ECgrid
!------------------------------------------------------------------------------

   CONTAINS
   
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
   character(len=80)           :: fmt
   integer   :: ix,iy,iz,ix1,ix2,iy1,iy2,iz1,iz2,idx,idy,idz,&
                itempx,itempy,itempz,nx,ny,nz,npts,i,j,k,nec
   integer,optional,intent(in) :: necs
   logical,optional,intent(in) :: quiet
   logical                     :: silent
   real(rs)   :: ecs_in(6,6),multiplier
   real(rs)   :: x,y,z,  dx,dy,dz
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
   
   return
   
   end subroutine EC_grid_load
!------------------------------------------------------------------------------

!==============================================================================
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
   
   return
   
   end subroutine EC_grid_write
!------------------------------------------------------------------------------

!==============================================================================
   subroutine EC_grid_load_bin(fname,grid)
!  Load the ECgrid type from a binary file

   implicit none
   
   character(len=*),intent(in) :: fname
   type(ECgrid),intent(out) :: grid
   integer :: iostatus,ix,iy,iz,k
   
   write(lu_stderr,'(a)') 'EC_grid_load_bin hasn''t been tested yet!'
   stop
   
   open(20,file=fname,iostat=iostatus,form='unformatted',status='old')
   if (iostatus /= 0) then
      write(lu_stderr,'(a,a)') &
                     'EC_grid_load_bin: cannot open file for reading: ',fname
      stop
   endif
   
   read(20) grid%npts, grid%nx, grid%ny, grid%nz, grid%ix1, grid%ix2, &
             grid%iy1, grid%iy2, grid%iz1, grid%iz2, grid%idx, grid%idy, grid%idz
   
   do ix=1,grid%nx
      do iy=1,grid%ny
         do iz=1,grid%nz
            read(20) grid%ecs(ix,iy,iz,:,:)
         enddo
      enddo
   enddo
   
   close(20)
   
!  Fill in dimensions
   do k=0,grid%nx-1; grid % x(k+1) = real(grid%ix1) + real(grid%idx) * real(k); enddo
   do k=0,grid%ny-1; grid % y(k+1) = real(grid%iy1) + real(grid%idy) * real(k); enddo
   do k=0,grid%nz-1; grid % z(k+1) = real(grid%iz1) + real(grid%idz) * real(k); enddo

   return
   
   end subroutine EC_grid_load_bin
!------------------------------------------------------------------------------   

!==============================================================================
   subroutine EC_grid_write_bin(fname,grid)
!  Write the ECgrid type to a binary, direct-access file.
!  Format is:
!	[1] npts,nx,ny,nz,ix1,ix2,iy1,iy2,iz1,iz2,idx,idy,idz
!   [2+] c11,c12...c16,c21...c66 (loops over x slowest, z fastest)

   implicit none
   
   type(ECgrid),intent(in)  :: grid
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
   
   return
   
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
   integer                       :: iostatus,ix,iy,iz,kx,ky,kz,i,j,nec
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
   
   return
   
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
      
      return
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
   
   return
   
   end subroutine EC_grid_delete
!------------------------------------------------------------------------------

!==============================================================================
   function EC_grid_check_necs(fname)
!==============================================================================
!  Determines the number of ECs per line in an EC_grid-type file

   implicit none
   
   character(len=*),intent(in) :: fname
   integer           :: EC_grid_check_necs, necs
   type (ECgrid)     :: grid
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
   read(20,*) x,y,z,test
   if (abs(test(2) - test(7)) < test_comparison) then
      EC_grid_check_necs = 36
   else
      EC_grid_check_necs = 21
   endif
  
   close(20)
      
   return

   end function EC_grid_check_necs
   
   

!______________________________________________________________________________
end module EC_grid
!==============================================================================