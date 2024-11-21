c**********************************************************************
      subroutine ecrrct
      include "paramt.h"
      common /constc/ tcs, bx0, rho0, slx, nx, nxp1, nxp2, npt, ns
      common /fieldc/ ex(ix), ey(ix), ez(ix), by(ix), bz(ix),
     &                ajx(ix), ajy(ix), ajz(ix), rho(ix)
      common /work1c/ work1(ix),work2(ix)
      common /ecrctc/ rkfact(ix)
c
      nxh=nx/2
      do 100 i=2,nxp1
        work1(i-1) = rho(i) - ex(i) + ex(i-1)
  100 continue
      call realft(work1,nx,1)
      do 200 i=1,nx
        work1(i) = work1(i)*rkfact(i)
  200 continue
      call realft(work1,nx,-1)
      work1(nxp1) = work1(1)
      do 300 i=2,nxp1
        ex(i) = ex(i) + work1(i-1) - work1(i)
  300 continue
        ex(1) = ex(nxp1)
      return
      end
