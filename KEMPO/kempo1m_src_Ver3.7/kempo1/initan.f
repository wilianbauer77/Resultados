c*******************************************************************
      subroutine initan
      include "paramt.h"
      common /inputc/ dx, dt, cv, wc, angle
      common /ptprmc/ wp(is), qm(is), q(is), vpe(is), vpa(is),
     &                vd(is), pch(is), np(is)
      common /constc/ tcs, bx0, rho0, slx, nx, nxp1, nxp2, npt, ns
      common /rotatc/ sinth, costh
      common /ecrctc/ rkfact(ix)
      common /prtclc/ x(in), vx(in), vy(in),vz(in)
      common /fieldc/ ex(ix), ey(ix), ez(ix), by(ix), bz(ix),
     &                ajx(ix), ajy(ix), ajz(ix), rho(ix)
      common /resclc/ rex, ret, rev, ree, reb, rej, rer, res
c
      dimension xs(is),xl(is)
      dimension unr1(1024), str1(1024), str2(1024), str3(1024)
c
      twopi = 6.283185308
      theta = twopi/360.0*angle
      sinth = sin(theta)
      costh = cos(theta)
      bx0  = wc/qm(1)*costh
      by0  = wc/qm(1)*sinth
      tcs = 2.0*cv*cv
      slx = nx
      nxp1 = nx + 1
      nxp2 = nx + 2
c
      npt=0
      rho0 = 0.0
      do 10 k = 1,ns
        xs(k) = 0.0
        xl(k) = slx
        npt   = npt + np(k)
        q(k)  = (slx/float(np(k))) * (wp(k)**2) / qm(k)
        rho0  = rho0 + q(k)*np(k)
   10 continue
      rho0 = -rho0/slx
c
      rkmin = twopi/slx
      nxh = nx/2
      fft = 1.0/float(nxh)
      do 300 i=1,nxh-1
        rk  = sin(rkmin*i*0.5)*2.0
        rkfact(2*i+1) = (1.0/rk**2) * fft
        rkfact(2*i+2) = rkfact(2*i+1)
  300 continue
      rkfact(1) = 0.0
      rk  = sin(rkmin*nxh*0.5)*2.0
      rkfact(2) = (1.0/rk**2) * fft
c ------------- particle initialization -----------------
      l = 0
      m = 0
      n2=0
      do 200 k=1,ns
        n1=n2
        n2=n1+np(k)
        phi = twopi/360.0*pch(k)
        vdpa = vd(k)*cos(phi)
        vdpe = vd(k)*sin(phi)
        rkk = rkmin*2
c       vmod = 2.0*rev
c       xmod = vmod/rkk
        mm = (n2 - n1)/1024
        left = n2-n1 - 1024*mm
        if(left>0) mm = mm + 1
        i1 = n1
        do ii = 1, mm
        i2 = i1 + 1024
c       call strnda(str1, 1024)
c       call strnda(str2, 1024)
c       call strnda(str3, 1024)
c       call unrnda(unr1, 1024)
        do i = i1 + 1, min(i2,n2)  
c         l = i - i1
          x(i)  = xs(k) + xl(k)*(i - n1 -1)/float(np(k))
          vxi   = vpa(k)*strndm(l) + vdpa
c         vxi  =  vxi +  vmod*cos(rkk*x(i))
          phase = twopi*unrndm(m)
          vyi   = vpe(k)*strndm(l) + vdpe*cos(phase)
          vz(i) = vpe(k)*strndm(l) + vdpe*sin(phase)
c         if(k==ns) then
c---- modification for non-gyrotropy
c         vri = sqrt(vyi*vyi + vz(i)*vz(i))
c         vyi = 0.
c--- non-gyrotropy without current
c         vz(i) = vri*(mod(i,2)*2-1)
c--- non-gyrotropy with current
c           vz(i) = vri
c-----------------------------------------------------
c         endif
c   rotation to the direction of the magnetic field
          vx(i) = costh*vxi - sinth*vyi
          vy(i) = sinth*vxi + costh*vyi
c         x(i) =    x(i) +  xmod*sin(rkk*x(i))
        end do
        i1 = i1 + 1024
        end do 
  200 continue
c ------------- field initialization -----------------
      do 20 i = 1,nxp2
        ex(i) = 0.0
        ey(i) = 0.0
        ez(i) = 0.0
        by(i) = by0
        bz(i) = 0.0
   20 continue
      return
      end
