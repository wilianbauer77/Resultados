c*******************************************************************
      subroutine bfield
      include "paramt.h"
      common /constc/ tcs, bx0, rho0, slx, nx, nxp1, nxp2, npt, ns
      common /fieldc/ ex(ix), ey(ix), ez(ix), by(ix), bz(ix),
     &                ajx(ix), ajy(ix), ajz(ix), rho(ix)
c
      do 200 i=2,nxp1
        by(i) = by(i) + ez(i) - ez(i-1)
        bz(i) = bz(i) - ey(i+1) + ey(i)
  200 continue
      by(nxp2) = by(2)
      bz(nxp2) = bz(2)
      bz(1)    = bz(nxp1)
      return
      end
