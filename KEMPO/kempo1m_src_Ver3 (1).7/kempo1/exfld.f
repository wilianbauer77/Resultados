c*******************************************************************
      subroutine exfld
      include "paramt.h"
      common /constc/ tcs, bx0, rho0, slx, nx, nxp1, nxp2, npt, ns
      common /fieldc/ ex(ix), ey(ix), ez(ix), by(ix), bz(ix),
     &                ajx(ix), ajy(ix), ajz(ix), rho(ix)
      common /exfldc/ w(im),cf(im),cb(im),pw(im)
      common /excnst/ w0,wp0,dep,det,dex,sgn
c
      do 100 i=2,nxp1
        ex(i) = ex(i) - 2.*ajx(i)
        ey(i) = ey(i) - tcs*( bz(i) - bz(i-1) ) - 2.*ajy(i)
        ez(i) = ez(i) + tcs*( by(i+1) - by(i) ) - 2.*ajz(i)
  100 continue
      ex(1)    = ex(nxp1)
      ey(nxp2) = ey(2)
      ez(1)    = ez(nxp1)
      return
      end
