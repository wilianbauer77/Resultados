c**********************************************************************
      subroutine positn
      include "paramt.h"
      common /constc/ tcs, bx0, rho0, slx, nx, nxp1, nxp2, npt, ns
      common /prtclc/ x(in), vx(in), vy(in),vz(in)
c
      do 100 i = 1, npt
        x(i) = x(i)+vx(i)
        if(x(i)<0.0) x(i) = x(i)+slx
        if(x(i).ge.slx) x(i) = x(i)-slx
  100 continue
      return
      end
