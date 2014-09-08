!===============================================================================
program phasevels
!===============================================================================
!  Read in ecs from stdin or from a file, and output information for a given
!  inclination and azimuth
!
!  Relies on CIJ_phase_vels in module anisotropy_ajn.  See the source
!  for information on the coordinate system and angle conventions.

   use anisotropy_ajn, only: CIJ_phase_vels, CIJ_load

   implicit none
   
   real(8) :: ecs(6,6),inc,azi,rho,pol,avs,vp,vs1,vs2,vsmean
   character(len=250) :: file,arg
   character(len=1000) :: line
   integer :: iostatus = 0, narg
   
   narg = command_argument_count()
   if (narg /= 1 .and. narg /= 2 .and. narg /= 3) then
      write(0,'(a)') 'Usage: CIJ_phase_vels [inc] [azi] (ecfile)',&
                     'or:    CIJ_phase_vels [ecfile] < (inc,azi on stdin)',&
                     '   Inc is angle from 1-2 plane towards 3',&
                     '   Azi is angle from 1 towards 2 in 1-2 plane',&
                     '   If no input file provided, ecs (density-normalised) are read from stdin, c11, c12, etc. (36)'
      stop
   endif
   
   ! Looping over sets of orientations on stdin
   if (narg == 1) then
      call get_command_argument(1,file)
      call CIJ_load(file,ecs,rho)
      call write_header_line
      iostatus = 0
      do while (iostatus == 0)
         read(*,'(a)',iostat=iostatus) arg
         if (iostatus < 0) exit
         if (iostatus > 0) then
            write(0,'(a)') 'CIJ_phase_vels: Error: Problem reading line from stdin'
            stop
         endif
         read(arg,*,iostat=iostatus) inc,azi
         azi = modulo(azi,360._8)
         if (iostatus /= 0) then
            write(0,'(a)') &
               'CIJ_phase_vels: Error: Problem reading inc,azi from line "' &
                  //trim(arg)//'"'
            stop
         endif
         call CIJ_phase_vels(ecs/rho,azi,inc,pol=pol,avs=avs,vp=vp,vs1=vs1,vs2=vs2)
         call write_output
      enddo

   else
      call get_command_argument(1,arg) ;  read(arg,*) inc
      call get_command_argument(2,arg) ;  read(arg,*) azi

!  Get elastic constants
!  If reading from an .ecs file, MUST NOT BE DENSITY-NORMALISED!!!
      if (narg == 3) then  ! One set from input file
         call get_command_argument(3,file)
         call CIJ_load(file,ecs,rho)
!  Check whether we're in GPa, not Pa
         if (ecs(1,1) < 5000.) ecs = ecs * 1.e9
         call CIJ_phase_vels(ecs/rho,azi,inc,pol=pol,avs=avs,vp=vp,vs1=vs1,vs2=vs2)
         call write_header_line
         call write_output

!  If reading ecs from stdin, MUST BE DENSITY-NORMALISED!!
      else if (narg == 2) then  ! Many sets from stdin
         call write_header_line
         do while (iostatus == 0)
            call read_ecs_stdin
            if (iostatus < 0) exit
!  Check whether we're in GPa, not Pa
            if (ecs(1,1) < 5000.) ecs = ecs * 1.e9
            call CIJ_phase_vels(ecs,azi,inc,pol=pol,avs=avs,vp=vp,vs1=vs1,vs2=vs2)
            call write_output
         enddo
      endif
   endif

contains
   subroutine write_header_line
      if (narg == 1) then
         write(*,'(a)') '   pol      avs        vp       vs1       vs2    inc    azi'
      else
         write(*,'(a)') '   pol      avs        vp       vs1       vs2'
      endif
   end subroutine write_header_line

   subroutine read_ecs_stdin
      read(*,'(a)',iostat=iostatus) line
      if (iostatus < 0) return
      if (iostatus > 0) then
         write(0,'(a)') 'CIJ_phase_vels: Error: Problem reading line from stdin'
         stop
      endif
      read(line,*,iostat=iostatus) ecs
      if (iostatus /= 0) then
         write(0,'(a)') 'CIJ_phase_vels: Error: Problem reading 36 ecs from line "' &
            //trim(line)//'"'
         stop
      endif
   end subroutine read_ecs_stdin

   subroutine write_output
      call ms2kms
      if (narg == 1) then
         write(*,'(f6.1,f9.4,3f10.4,2(f7.1))') pol,avs,vp,vs1,vs2,inc,azi
      else
         write(*,'(f6.1,f9.4,3f10.4)') pol,avs,vp,vs1,vs2
      endif
   end subroutine write_output

   subroutine ms2kms
      ! Convert velocities to km/s from m/s
      vp = vp/1000._8
      vs1 = vs1/1000._8
      vs2 = vs2/1000._8
   end subroutine ms2kms
end program phasevels
