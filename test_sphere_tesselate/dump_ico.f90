program dump_ico
! print out triangles from icosahedron

! use the sphere_tesselate module 
! (this is what we want to test after all)
use sphere_tesselate

! implicit none as always
implicit none



! declare a tesselation which defines our icosahedron
! (kind of like a mesh)
type(tesselation) :: t

! build the icosadron with the following arguments:
! -- level = 0
! -- shape = icosahedron
! -- point at poles = true
t = st_new(0,'ico',.true.)

! now spit out the triangles
! -- t is our icosahedron (in the tesselation structure)
! -- .true. says we want geographical co-ordinates (lat,lon)
call st_dump_triangles(t,.true.)

! print *,"Test"

end program
