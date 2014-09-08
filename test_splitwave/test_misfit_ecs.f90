program test
   ! Test that sw_misfit_ecs gives the correct answers
   use spherical_geometry, only: sphere_sample
   use splitwave
   use anisotropy_ajn, only: CIJ_phase_vels, CIJ_rot3
   implicit none

   ! Constants for olivine (Abramson)
   integer, parameter :: rs = 8
   real(rs) :: Col(6,6) = 1.e7_8 * reshape( (/ &
        9.5529_8,  2.0298_8,  2.1341_8,      0._8,      0._8,      0._8, &
        2.0298_8,  5.8569_8,  2.2891_8,      0._8,      0._8,      0._8, &
        2.1341_8,  2.2891_8,  6.9598_8,      0._8,      0._8,      0._8, &
            0._8,      0._8,      0._8,  1.9076_8,      0._8,      0._8, &
            0._8,      0._8,      0._8,      0._8,  2.2951_8,      0._8, &
            0._8,      0._8,      0._8,      0._8,      0._8,  2.3457_8/), (/6,6/))
   real(rs) :: rho_ol = 3355._rs
   real(rs) :: vs1, vs2, a, b, c
   real(rs), allocatable :: lon(:), lat(:), phi(:), dt(:), spol(:), &
                            phi_ecs(:), dt_ecs(:), misfit(:)
   integer :: i, n, time(8), size_rand
   integer, allocatable :: seed(:)
   real(rs) :: d = 10._rs, &  ! Spacing of points on sphere
               t = 500._rs,&  ! Layer thickness
               t_scaled       ! Scaled layer thickness from sw_misfit_ecs

   ! Sample the sphere evenly
   call sphere_sample(d, lon, lat, n)
   write(*,'(a,i0.0,a)') 'Testing with ', n, ' orientations'
   allocate(phi(n), dt(n), spol(n), phi_ecs(n), dt_ecs(n), misfit(n))

   ! Compute splits for these orientations
   do i=1,n
      call CIJ_phase_vels(Col, lon(i), lat(i), pol=phi(i), vs1=vs1, vs2=vs2)
      dt(i) = t*(1._rs/vs2 - 1._rs/vs1)
   enddo

   ! Assign random spol
   call random_seed(size=size_rand)
   allocate(seed(size_rand))
   call date_and_time(values=time)  ! Values are y,m,d,zone,h,m,s,msec
   seed = time(8)
   call random_seed(put=seed)      ! Seed with total of values
   call random_number(spol)
   spol = 180._rs * spol

   ! Compare these splits with those predicted in sw_misfit_ecs
   call sw_misfit_ecs(Col, lon, lat, phi, dt, spol, misfit, phi_ecs=phi_ecs, &
      dt_ecs=dt_ecs, t_scaled=t_scaled)
   write(*,'(a,f10.8)') '   Mean misfit v. self: ', sum(misfit)/size(misfit)
   write(*,'(a,f10.2)') '   Scaled thickness: ',t_scaled

   ! Now compare with olivine tensor rotated randomly
   call random_seed()
   call random_number(a)
   a = 360._rs*a
   call random_number(b)
   b = 360._rs*b
   call random_number(c)
   c = 360._rs*c
   call CIJ_rot3(Col, a, b, c, Col)
   call sw_misfit_ecs(Col, lon, lat, phi, dt, spol, misfit, phi_ecs=phi_ecs, &
      dt_ecs=dt_ecs, t_scaled=t_scaled)
   write(*,'(a,f10.8)') '   Mean misfit v. rand: ',sum(misfit)/size(misfit)
   write(*,'(a,f10.2)') '   Scaled thickness: ',t_scaled

   deallocate(lon, lat, phi, dt, spol, phi_ecs, dt_ecs, misfit)

end program test
