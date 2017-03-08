program load_neighbours_from_cache

! use the sphere_tesselate module 
! (this is what we want to test after all)
use sphere_tesselate

! implicit none as always
implicit none

! declare a tesselation which defines our icosahedron
! (kind of like a mesh)
! type(tesselation) :: tess
type(point) :: pt
type(triangle) :: tr
type(list_neighbours) :: lnbrs
integer :: pt_id, level


! build the icosadron with the following arguments:
level = 2
! -- shape = icosahedron
! -- point at poles = true
! tess = st_new(level,'ico',.true.)


!!!! Read from cache
write(*,*) 'Testing Fortran binary formatted neighbours files'
! save list to ascii file for human inspection
call st_load_neighbours_cache(level,lnbrs)

! Check basic info
write(*,*) 'lnbrs level, nt, np:',lnbrs%level,lnbrs%nt,lnbrs%np

! check it's read in properly
pt_id = 98
write(*,*) 'pt_id:',pt_id
write(*,*) 'neighbours:',lnbrs%list(pt_id)%idxs


end program
