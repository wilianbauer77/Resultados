c**********************************************************************
      subroutine charge
      include "paramt.h"
      common /constc/ tcs, bx0, rho0, slx, nx, nxp1, nxp2, npt, ns
      common /fieldc/ ex(ix), ey(ix), ez(ix), by(ix), bz(ix),
     &                ajx(ix), ajy(ix), ajz(ix), rho(ix)
      common /prtclc/ x(in), vx(in), vy(in),vz(in)
      common /ptprmc/ wp(is), qm(is), q(is), vpe(is), vpa(is),
     &                vd(is), pch(is), np(is)
c
      do 100 i=1,nxp2
        rho(i) = rho0
  100 continue
c
          n2 = 0
      do 210 k=1,ns 
        n1 = n2
        n2 = n1 + np(k)
        do 200 m = n1+1, n2
           i  =  x(m)+ 2.0
           s2 = (x(m)+ 2.0 - i)*q(k)
           s1 = q(k) - s2
           rho(i)   = rho(i)   + s1
           rho(i+1) = rho(i+1) + s2
  200   continue
  210 continue
      rho(2)    = rho(2) + rho(nxp2) - rho0
      rho(1)    = rho(nxp1)
      rho(nxp2) = rho(2)
      return
      end
