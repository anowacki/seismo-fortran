program test

use get_args

real(8) :: dt
real(4) :: r
character(50) :: f1,f2,f3,f4
integer :: i

f1='xxx'
f2=f1
f3=f1

!call get_arg('-dt',dt)
!write(*,*) 'dt = ',dt
!
!call get_arg('-n',i)
!write(*,*) 'i = ',i
!
call get_arg('-f',f1,f2,f3)
write(*,*) 'files: ',trim(f1)//' '//trim(f2)//' '//trim(f3)

end program test
