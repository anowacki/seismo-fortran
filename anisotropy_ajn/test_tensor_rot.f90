program test_tensor_rot
!  Test that my own implementation of CIJ_rotate3() is correct
!  It is, but it's much slower than Mainprice's
   
   use anisotropy_ajn   
   
   implicit none
   
   integer,parameter :: rs=8
   real,parameter :: pi = 3.141592653589793238462643_rs
   real(rs) :: T(3,3,3,3), Tr(3,3,3,3), CIJ(6,6), CIJr(6,6), CIJ1(6,6)
   real(rs) :: alp,bet,gam
   integer  :: i,j,k,l,a,b,c,d, n
   integer  :: seconds0=0, seconds1=0, rate
   real(rs),dimension(3,3) :: R,Ra,Rb,Rc,g   ! Rotation matrices
   integer :: nrot=5000
   character(len=80) :: arg
   
   if (iargc() == 1) then
      call getarg(1,arg)
      read(arg,*) nrot
   endif   
   
!************************************************
!  Rotate tensor the explicit way

!  Make a symmetrical matrix
   do i=1,6
      do j=1,6
         CIJ(i,j) = real(i*j) !*1.d9
      enddo
   enddo

!  Start clock
   seconds0 = 0
   call system_clock(seconds0,rate)
   
!  Fill in T from C
   T = Cij2cijkl(CIJ)
   
!  Rotation angles
   alp = 25._rs * pi / 180._rs
   bet = 0. !30._rs * pi / 180._rs
   gam = 0. !45._rs * pi / 180._rs
   
!  Create rotation matrices
   Ra(1,1) =  1.     ; Ra(1,2) =  0.     ; Ra(1,3) =  0.
   Ra(2,1) =  0.     ; Ra(2,2) =  cos(alp) ; Ra(2,3) =  sin(alp)
   Ra(3,1) =  0.     ; Ra(3,2) = -sin(alp) ; Ra(3,3) =  cos(alp)
   
   Rb(1,1) =  cos(bet) ; Rb(1,2) =  0.     ; Rb(1,3) = -sin(bet)
   Rb(2,1) =  0.     ; Rb(2,2) =  1.     ; Rb(2,3) =  0.
   Rb(3,1) =  sin(bet) ; Rb(3,2) =  0.     ; Rb(3,3) =  cos(bet)
   
   Rc(1,1) =  cos(gam) ; Rc(1,2) =  sin(gam) ; Rc(1,3) =  0.
   Rc(2,1) = -sin(gam) ; Rc(2,2) =  cos(gam) ; Rc(2,3) =  0.
   Rc(3,1) =  0.     ; Rc(3,2) =  0.     ; Rc(3,3) =  1.
   
   R = matmul(Rb,Ra)
   g = matmul(Rc,R)
   
!  Rotate tensor 500 times:
   do n=1,nrot
   Tr = 0._rs
   do l=1,3
      do k=1,3
         do j=1,3
            do i=1,3
               do d=1,3
                  do c=1,3
                     do b=1,3
                        do a=1,3
                           Tr(i,j,k,l) = &
                              Tr(i,j,k,l) + g(i,a)*g(j,b)*g(k,c)*g(l,d)*T(a,b,c,d)
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
   enddo
   enddo

!  Stop clock
   call system_clock(seconds1,rate)
   
!  Show results
   CIJr = cijkl2Cij(Tr)
   CIJ1 = CIJr
   write(*,'(a,f0.4,a)') 'Output from inbuilt rotation took ', &
                  (seconds1-seconds0)/real(rate),' seconds:'
   write(*,'(6f12.4)') CIJr
   write(*,*)
   
!************************************************
!  Call CIJ_rot3
   
!  Make a symmetrical matrix
   do i=1,6
      do j=1,6
         CIJ(i,j) = real(i*j) !*1.d9
      enddo
   enddo
   
!  Compare to CIJ_rot3
   call system_clock(seconds0,rate)
   do n=1,nrot
      call CIJ_rot3(CIJ,alp*180./pi,bet*180./pi,gam*180./pi,CIJr)
   enddo
   call system_clock(seconds1,rate)
   write(*,'(a,f6.4,a)') 'Output from CIJ_rot3 took: ', &
                  (seconds1-seconds0)/real(rate),' seconds:'
   write(*,'(6f12.4)') CIJr
   write(*,*)
   
   write(*,'(a)') 'Percentage error:'
   write(*,'(6f7.0)') 100._rs*(CIJ1-CIJr)/CIJr
   
end program