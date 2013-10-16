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
   complex(rs), allocatable :: f(:)
   real(rs), allocatable :: fval(:)
   integer :: i
   
   ! Sampling interval
   dt = 0.01_rs ! / seconds
   ! Frequency of sine wave
   freq = 2 ! / Hz
   
   ! Frequency sampling interval
   df = 1._rs/(dt*real(npts, kind=rs))
   
   ! Create signal
   do i=1,npts
      t(i) = sin(2._rs*pi*real(i-1, kind=rs)*dt*freq) &
             + sin(2._rs*pi*real(i-1,rs)*dt*freq*2._rs)
   enddo
   t = t/maxval(abs(t))
   t_orig = t
   
   ! Write out original signal
   open(10,file="input.xy")
   write(10,'(e12.6,1x,e12.6)') (real(i-1,rs)*dt,t(i), i=1,size(t))
   close(10)
   
   ! Create Fourier-transformed trace
   call FFFTW_fwd(t,f)
   
   ! Allocate memory for frequency values
   allocate(fval(size(f)))
   
   ! Fill in time and frequency values
   tval = [(real(i-1,rs)*dt, i=1,size(t))]
   fval = [(real(i-1,rs)*df, i=1,size(fval))]
   
   ! Low-pass inf.-order filter in FD
!    where (fval > 3.0_rs)
!       f = complex(0._rs,0_rs)
!    end where
   
   ! Write out amplitude spectrum of FFT
   open(10,file="amplitude.xy")
   write(10,'(e12.6,1x,e12.6)') (fval(i),abs(f(i)), i=1,size(f))
   close(10)
   ! Phase
   open(10,file="phase.xy")
   write(10,'(e12.6,1x,e12.6)') (fval(i),atan2(real(f(i)),aimag(f(i))), i=1,size(f))
   close(10)
   
   ! Create inverse trace
   call FFFTW_rev(f,t)
   
   ! Write out FFT-iFFT trace
   open(10,file="output.xy")
   write(10,'(e12.6,1x,e12.6)') (tval(i),t(i), i=1,size(t))
   close(10)
   
   ! Write out difference between input and output
   open(10,file="output-input.xy")
   write(10,'(e12.6,1x,e12.6)') (tval(i),abs(t(i)-t_orig(i)), i=1,size(t))
   close(10)
   
   write(*,'(a,e12.6)') 'Average residual between input and output traces: ', &
      sum(abs(t-t_orig))/size(t)
   write(*,'(a,e12.6)') 'Maximum residual between input and output traces: ', &
      maxval(abs(t-t_orig))
   
   
end program
