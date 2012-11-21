!
!  Simple test program for f90sac library
!
   program test_f90sac
   use f90sac  ! use the f90sac module
      implicit none
      type (SACTrace) :: t1,t2,t3,tc
      character (len=80) :: fn
       

      fn = 'test.sac' ;
      call f90sac_newtrace(101,0.05,tc)

      tc % trace(51) = 1.0
      tc % evla = 30.0
      call f90sac_writetrace(fn,tc)

      
      if (f90sac_isBigEndian()) then
         print*,'Machine is Big-Endian'
      else
         print*,'Machine is Little-Endian'
      endif
         
      if (f90sac_force_byteswap) then
         print*,'Forced byte-swapping is ON'
      else
         print*,'Forced byte-swapping is OFF'
      endif   
      
      print*,'Now test that you can read TEST.SAC into your version of SAC'
!
!     ** now check that we can edit the header
!      
      tc % evla = 60.
      call f90sac_writeheader(fn,tc)
            
   end program test_f90sac
