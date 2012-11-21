!===============================================================================
program test_sphere_sample
!===============================================================================
!  Ensure the sphere_sample subroutine works as expected

   use spherical_geometry
   
   implicit none
   
   real(8),allocatable, dimension(:) :: lat,lon
   integer :: n  ! Number of points
   integer :: i
   real(8) :: d  ! Distance between points
   
   d = 15._8
   
   call sphere_sample(d,lon,lat,n)
   
   do i=1,n
      write(*,'(i4.4,x,2f6.1)') i,lon(i),lat(i)
   enddo
   
end program test_sphere_sample
!-------------------------------------------------------------------------------