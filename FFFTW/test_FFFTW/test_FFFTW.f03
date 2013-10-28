! Test the FFFTW module: create a wave with specified frequencies and
! write the input, output time-domain traces, plus the amplitude and phase.
program test
   use FFFTW
   implicit none
   integer, parameter :: rs = 8
   integer, parameter :: npts = 2**11-1
   real(rs), parameter :: pi = 4._rs*atan2(1._rs, 1._rs)
   real(rs) :: t(npts), tval(npts), t_orig(npts)
   real(rs) :: dt, df, freq
   complex(rs) :: tc(npts), tc_orig(npts)
   complex(rs), allocatable :: f(:), fc(:)
   real(rs), allocatable :: fval(:), fvalc(:)
   integer :: i
   
   ! Sampling interval
   dt = 0.01_rs ! / seconds
   ! Frequency of sine wave
   freq = 2._rs ! / Hz
   
   ! Frequency sampling interval
   df = 1._rs/(dt*real(npts, kind=rs))
   
   ! Create signal for real and complex traces
   do i=1,npts
      t(i) = sin(2._rs*pi*real(i-1, kind=rs)*dt*freq) &
             + sin(2._rs*pi*real(i-1,rs)*dt*freq*2._rs)
   enddo
   ! Normalise amplitdue to ±1
   t = t/maxval(abs(t))
   tc = cmplx(t,-t, kind=rs)
   t_orig = t
   tc_orig = tc
   
   ! Write out original signal: real
   open(10,file="input_real.xy")
   write(10,'(e12.6,1x,e12.6)') (real(i-1,rs)*dt,t(i), i=1,size(t))
   close(10)
   ! Complex
   open(10,file="input_complex.xy")
   write(10,'(e12.6,1x,e12.6,1x,e12.6)') (real(i-1,rs)*dt,real(tc(i)),aimag(tc(i)), i=1,size(t))
   close(10)
   
   ! Create Fourier-transformed trace
   call FFFTW_fwd(t,f)
   call FFFTW_fwd(tc,fc)
   
   ! Allocate memory for frequency values
   allocate(fval(size(f)), fvalc(size(fc)))
   
   ! Fill in time and frequency values
   tval = [(real(i-1,rs)*dt, i=1,size(t))]
   fval = [(real(i-1,rs)*df, i=1,size(fval))]
   fvalc = [(real(i-1,rs)*df, i=1,size(fvalc))]
   
   ! Low-pass inf.-order filter in FD
!    where (fval > 3.0_rs)
!       f = complex(0._rs,0_rs)
!    end where
   
   ! Write out amplitude spectrum of FFT
   ! Real-to-complex
   open(10,file="amplitude_real.xy")
   write(10,'(e12.6,1x,e12.6)') (fval(i),abs(f(i)), i=1,size(f))
   close(10)
   ! Phase
   open(10,file="phase_real.xy")
   write(10,'(e12.6,1x,e12.6)') (fval(i),atan2(real(f(i)),aimag(f(i))), i=1,size(f))
   close(10)
   ! Complex-to-complex
   open(10,file="amplitude_complex.xy")
   write(10,'(e12.6,1x,e12.6)') (fvalc(i),abs(fc(i)), i=1,size(fc))
   close(10)
   ! Phase
   open(10,file="phase_complex.xy")
   write(10,'(e12.6,1x,e12.6)') (fvalc(i),atan2(real(fc(i)),aimag(fc(i))), i=1,size(fc))
   close(10)
   
   ! Create inverse trace
   call FFFTW_rev(f,t)
   call FFFTW_rev(fc,tc)
   
   ! Write out FFT-iFFT trace: real
   open(10,file="output_real.xy")
   write(10,'(e12.6,1x,e12.6)') (tval(i),t(i), i=1,size(t))
   close(10)
   
   ! Write out difference between input and output: real
   open(10,file="output-input_real.xy")
   write(10,'(e12.6,1x,e12.6)') (tval(i),abs(t(i)-t_orig(i)), i=1,size(t))
   close(10)
   
   write(*,'(a,e12.6)') 'Real: Average residual between input and output traces: ', &
      sum(abs(t-t_orig))/size(t)
   write(*,'(a,e12.6)') 'Real: Maximum residual between input and output traces: ', &
      maxval(abs(t-t_orig))
   
   ! Write out FFT-iFFT trace: complex
   open(10,file="output_complex.xy")
   write(10,'(e12.6,1x,e12.6,1x,e12.6)') (tval(i),real(tc(i)),aimag(tc(i)), i=1,size(tc))
   close(10)

   ! Write out difference between input and output: complex
   open(10,file="output-input_complex.xy")
   write(10,'(e12.6,1x,e12.6)') (tval(i),abs(tc(i)-tc_orig(i)), i=1,size(tc))
   close(10)

   write(*,'(a,e12.6)') 'Complex: Average residual between input and output traces: ', &
      sum(abs(tc-tc_orig))/size(tc)
   write(*,'(a,e12.6)') 'Complex: Maximum residual between input and output traces: ', &
      maxval(abs(tc-tc_orig))
   
   
end program
