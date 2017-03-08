program generate_list_neighbours
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


write (*,*) 'Calling st_generate_list_neighbours'
call st_generate_list_neighbours(tess,lnbrs)

! show some elements in the list
pt_id = 98
write(*,*) 'pt_id:',pt_id
write(*,*) 'neighbours:',lnbrs%list(pt_id)%idxs

! show whole list
! do pt_id = 1,13
! write(*,*) 'whole list:',lnbrs%list(pt_id)%idxs
! end do

end program
