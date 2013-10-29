! Test the splitwave module.  Create a wave, split it in the time and frequency
! domain, then write out the split traces.
program test
   
   use splitwave
   use f90sac
   
   type(SACtrace) :: E,N,Z, Eorig,Norig,Zorig, Efd,Nfd,Zfd
   ! Wave parameters
   real :: delta = 0.1, freq = 0.1, spol = 45., noise = 0.01
   ! Splitting parameters
   real :: phi(1) = (/10./), dt(1) = (/2./)
   
   call sw_create_wave(E,N,Z,freq,delta,noise=noise,wavetype='r',spol=spol)
   Eorig = E
   Efd = E
   Norig = N
   Nfd = N
   Zorig = Z
   Zfd = Z
   
   ! Apply splitting of -phi
   call sw_fdsplitN(Nfd, Efd, Zfd, size(phi), -phi, dt, quiet=.false.)
   call sw_splitN(  N,   E,   Z,   size(phi), -phi, dt, quiet=.false.)
   
   call f90sac_writetrace('wave.BHE', Eorig)
   call f90sac_writetrace('wave.BHN', Norig)
   call f90sac_writetrace('wave.BHZ', Zorig)

   call f90sac_writetrace('wave_fdsplit.BHE', Efd)
   call f90sac_writetrace('wave_fdsplit.BHN', Nfd)
   call f90sac_writetrace('wave_fdsplit.BHZ', Zfd)

   call f90sac_writetrace('wave_tdsplit.BHE', E)
   call f90sac_writetrace('wave_tdsplit.BHN', N)
   call f90sac_writetrace('wave_tdsplit.BHZ', Z)
   
   
end program test