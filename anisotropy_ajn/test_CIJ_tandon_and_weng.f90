program test_CIJ_tandon_and_weng

use anisotropy_ajn

implicit none

      integer, parameter :: i4 = selected_int_kind(9) ; ! long int
      integer, parameter :: r4 = selected_real_kind(6,37) ; ! SP
      integer, parameter :: r8 = selected_real_kind(15,307) ; ! DP
      
!  ** precision selector
      integer, parameter :: rs = r8
      
!  ** maths constants and other useful things
      real(rs), parameter :: pi = 3.141592653589793238462643 ;
      real(rs), parameter :: to_rad = 1.74532925199433e-002 ;  
      real(rs), parameter :: to_deg = 57.2957795130823e0 ;  
      real(rs), parameter :: to_km = 111.194926644559 ;      

      real(rs), parameter :: big_number = 10.e36 ;      

   real(rs)   :: CC(6,6),CCexp(6,6),vp,vs,rho,vpi,vsi,rhoi,c,del,rhoout
   
   vp = 5.8d3
   vs = 3.46d3
   rho = 2.92d3
   
   del = 0.01
   c = 0.02
   
   vpi = 4.d3
   vsi = 0.d0
   rhoi = 2.2d3

   call CIJ_tandon_and_weng(vp,vs,rho,del,c,vpi,vsi,rhoi,CC,rhoout)
   
   write(*,'(a)')  'Output should look like: (x10^7)'
   write(*,'(a)') &
'    3.2634    0.9878    0.9878    0.0000    0.0000    0.0000',&
'    0.9878    3.3358    0.9784    0.0000    0.0000    0.0000',&
'    0.9878    0.9784    3.3358    0.0000    0.0000    0.0000',&
'    0.0000    0.0000    0.0000    1.1787    0.0000    0.0000',&
'    0.0000    0.0000    0.0000    0.0000    0.5574    0.0000',&
'    0.0000    0.0000    0.0000    0.0000    0.0000    0.5574'

   CCexp = 1.d7 * reshape( &
   (/  3.2634,  0.9878,  0.9878,  0.0000,  0.0000,  0.0000, &
       0.9878,  3.3358,  0.9784,  0.0000,  0.0000,  0.0000, &
       0.9878,  0.9784,  3.3358,  0.0000,  0.0000,  0.0000, &
       0.0000,  0.0000,  0.0000,  1.1787,  0.0000,  0.0000, &
       0.0000,  0.0000,  0.0000,  0.0000,  0.5574,  0.0000, &
       0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.5574  /), (/6,6/))

   write(*,'(/a)') 'Output: (x10^7)'
   write(*,'(6f10.4)') CC/1.d7
   
   ! Different if any value differs by more than 0.1 %
   if (any(abs(CC - CCexp) > abs(CCexp/1.d4))) then
      write(*,'(/a,e10.4)') 'TEST FAILED: total difference between sum of matrices = ',sum(abs(CC-CCexp))
   else
      write(*,'(/a)') 'Test passed'
   endif

end program test_CIJ_tandon_and_weng