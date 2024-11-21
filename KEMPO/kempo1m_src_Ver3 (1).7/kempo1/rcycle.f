******************************************************************
      subroutine rcycle
      include "paramt.h"
      common /inputc/ dx, dt, cv, wc, angle
      common /ptprmc/ wp(is), qm(is), q(is), vpe(is), vpa(is),
     &                vd(is), pch(is), np(is)
      common /constc/ tcs, bx0, rho0, slx, nx, nxp1, nxp2, npt, ns
      common /rotatc/ sinth, costh
      common /prtclc/ x(in), vx(in), vy(in),vz(in)
      common /ngyroc/ ngs, jperp, dpsi, ncycle
c
      data l,m /0,0/
c
      if(ngs<0.or.ngs>ns) return
      if(ncycle<=0) return
      twopi = 6.283185308
c--- angle of the static magnetic field
      theta = twopi/360.0*angle
c--- range of angle for nongyrotropy
      dpsir = twopi/360.0*dpsi
c--- pitch angle of the particle drift
      phi   = twopi/360.0*pch(ngs)
      vdpa  = vd(ngs)*cos(phi)
      vdpe  = vd(ngs)*sin(phi)
c
c ------------- Particle Reinitialization for Vperp -------------
c
      rnps = np(ngs)  - 0.001
      do 100 k=1,ncycle
        i = 1 + rnps*unrndm(m)
c       vxi = costh*vx(i) + sinth*vy(i)
        vxi   = vpa(ngs)*strndm(l) + vdpa
        phase = twopi*unrndm(m)
        vyi   = vpe(ngs)*strndm(l) + vdpe*cos(phase)
        vzi   = vpe(ngs)*strndm(l) + vdpe*sin(phase)
        vri   = sqrt(vzi*vzi+vyi*vyi)
        psi   = twopi*0.25+dpsir*(unrndm(m) -0.5)
        if(jperp==0.and.mod(k,2)==0) vri = -vri
        vyi   = vri*cos(psi)
        vz(i) = vri*sin(psi)
        vx(i) = costh*vxi - sinth*vyi
        vy(i) = sinth*vxi + costh*vyi
  100 continue
      return
      end
