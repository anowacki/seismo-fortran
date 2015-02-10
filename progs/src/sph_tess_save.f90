!===============================================================================
! One of Andy Nowacki's Fortran utility modules for dealing with seismic
! anisotropy and other problems.
!
! Andy Nowacki <andy.nowacki@bristol.ac.uk>
!
! See the file LICENCE for licence details.
!===============================================================================
program sph_tess_save
!===============================================================================
! Create and save a sph_tess file for later use by a program which in turn
! uses the sphere_tesselate module to do stuff with in.

   use sphere_tesselate, only: &
      tesselation, point, triangle, st_icosahedron, st_new, st_save_tesselation, &
      st_load_tesselation

   implicit none
   type(tesselation) :: t, u
   character(len=250) :: file
   logical :: ascii = .false.
   integer :: level
   real(8), parameter :: point_tol = 1e-5
   
   call get_args
   
   write(*,'(a)') 'Creating tesselation...'
   t = st_new(level, 'icosahedron', pole=.true.)
   write(*,'(a)') 'Saving tesselation to disk...'
   call st_save_tesselation(t, file, ascii=ascii)
   write(*,'(a)') 'Reading tesselation from disk'
   call st_load_tesselation(file, u, ascii=ascii)
   write(*,'(a)') 'Checking input and output match...'
   call assert_points_equal(t%p, u%p)
   call assert_triangles_equal(t%t, u%t)
contains
   
   subroutine usage
      write(0,'(a)') &
         'Usage: sph_tess_save (options) [level] [file]', &
         'Create and save a file containing a tesselation of the unit sphere', &
         'for later use by other sphere_tesselate programs.', &
         'Specify the tesselation level >= 0.', &
         'Options:', &
         '   -a : Write an ASCII file, not a Fortran binary file'
      error stop
   end subroutine usage
   
   subroutine get_args
      character(len=250) :: arg
      integer :: iarg, narg
      narg = command_argument_count()
      if (narg < 2) call usage
      iarg = 1
      do while (iarg <= narg - 2)
         call get_command_argument(iarg, arg)
         select case(arg)
            case('-a');    ascii = .true.; iarg = iarg + 1
            case default;  call usage
         end select
      enddo
      call get_command_argument(narg - 1, arg)
      read(arg,*) level
      call get_command_argument(narg, file)
   end subroutine get_args
   
   subroutine assert_points_equal(a, b)
      type(point), intent(in), dimension(:) :: a, b
      if (any(abs(a%x - b%x) > point_tol) .or. any(abs(a%y - b%y) > point_tol) &
          .or. any(abs(a%z - b%z) > point_tol)) &
         write(0,'(a,3f11.8,a,f411.8)') 'Two points are not equal! ', a, ' != ', b
   end subroutine assert_points_equal
   
   subroutine assert_triangles_equal(t, u)
      type(triangle), intent(in), dimension(:) :: t, u
      if (any(t%a /= u%a) .or. any(t%b /= u%b) .or. any(t%c /= u%c)) &
         write(0,'(a)') 'Two triangles are not equal!'
   end subroutine assert_triangles_equal
end program sph_tess_save
!-------------------------------------------------------------------------------