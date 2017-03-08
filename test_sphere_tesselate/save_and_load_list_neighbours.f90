program save_and_load_list_neighbours
! given a point find all neighbouring points

! use the sphere_tesselate module 
! (this is what we want to test after all)
use sphere_tesselate

! implicit none as always
implicit none

! declare a tesselation which defines our icosahedron
! (kind of like a mesh)
type(tesselation) :: tess
type(point) :: pt
type(triangle) :: tr
type(list_neighbours) :: lnbrs
integer :: pt_id


! build the icosadron with the following arguments:
! -- level = 0
! -- shape = icosahedron
! -- point at poles = true
tess = st_new(2,'ico',.true.)

! generate list of neighbours
call st_generate_list_neighbours(tess,lnbrs)

!!!! Fortran Binary
write(*,*) 'Testing Fortran binary formatted neighbours files'
! save list to ascii file for human inspection
call st_save_neighbours(lnbrs,'test_list_of_neighbours.bin')
! now load it back in
call st_load_neighbours('test_list_of_neighbours.bin',lnbrs)
! check it's read in properly
pt_id = 98
write(*,*) 'pt_id:',pt_id
write(*,*) 'neighbours:',lnbrs%list(pt_id)%idxs


!!!! ASCII
write(*,*) 'Testing ASCII formatted neighbours files'

! save list to ascii file for human inspection

call st_save_neighbours(lnbrs,'test_list_of_neighbours.txt',.true.)

! now load it back in
call st_load_neighbours('test_list_of_neighbours.txt',lnbrs,.true.)
! check it's read in properly
pt_id = 98
write(*,*) 'pt_id:',pt_id
write(*,*) 'neighbours:',lnbrs%list(pt_id)%idxs






end program
