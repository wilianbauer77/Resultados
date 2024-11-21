c
      subroutine rfftcv(x,nx,icont,ix,work,iwk)
c
c     Real FFT using CONVEX VECLIB
c
c      programed by Y. Omura, 12/14/93
c
c     nx = number of elements of x,  must be a power of 2
c     ix = nx+2
c     iwk = 3*nx/2
c     icont = 1    FFT
c     icont = -1   Inverse FFT
c
c     Functions are basically same with the "realft" of LIBKEMPO1 
c     except for the sign of the sin( ) coefficients after FFT.
c
      real*4 x(ix),work(iwk)
      if(icont==1) then
        call src1ft(x,nx,work,-3,ier)
        call src1ft(x,nx,work,1,ier)
        x(2) = x(nx+1)
      else if(icont==-1) then
        x(nx+1) = x(2)
        call src1ft(x,nx,work,-2,ier)
        call sscal(nx,0.5,x,1)
      endif
      return
      end
