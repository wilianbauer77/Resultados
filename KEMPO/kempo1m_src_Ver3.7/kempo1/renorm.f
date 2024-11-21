c****************************************************************
      subroutine renorm
      include "paramt.h"
      common /constc/ tcs, bx0, rho0, slx, nx, nxp1, nxp2, npt, ns
      common /ptprmc/ wp(is), qm(is), q(is), vpe(is), vpa(is),
     &                vd(is), pch(is), np(is)
      common /otherc/ vmin,vmax
      common /resclc/ rex, ret, rev, ree, reb, rej, rer, res
      common /inputc/ dx, dt, cv, wc, angle
c
c-- distance
      rex = 1.0/dx
c-- time
      ret = 2.0/dt
c-- velocity
      rev = rex/ret
c-- electric field, charge, and mass
      ree = rex/(ret*ret)
c-- magnetic field
      reb = 1.0/ret
c-- current density
      rej = rex/(ret*ret*ret)
c-- charge density
      rer = 1.0/(ret*ret)
c-- energy density
      res = (rex*rex)/(ret*ret*ret*ret)
c
      vmin = vmin*rev
      vmax = vmax*rev
      cv   = cv  *rev
      wc   = wc  /ret
c
      do 10 k=1,ns
        wp(k)  = wp(k) /ret
        vpe(k) = vpe(k)*rev
        vpa(k) = vpa(k)*rev
        vd(k)  = vd(k) *rev
   10 continue
      return
      end
