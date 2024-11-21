c*******************************************************************
      subroutine bvlcty
      include "paramt.h"
      common /constc/ tcs, bx0, rho0, slx, nx, nxp1, nxp2, npt, ns
      common /inputc/ dx, dt, cv, wc, angle
      common /fieldc/ ex(ix), ey(ix), ez(ix), by(ix), bz(ix),
     #                ajx(ix), ajy(ix), ajz(ix), rho(ix)
      common /prtclc/ x(in), vx(in), vy(in),vz(in)
      common /ptprmc/ wp(is), qm(is), q(is), vpe(is), vpa(is),
     #                vd(is), pch(is), np(is)
      common /work1c/ work1(ix),work2(ix)
      real*8 ux,uy,uz,uxt,uyt,uzt,cvm,cs,ex1,ey1,ez1,bx1,by1,bz1
      real*8 g,boris
c
      cvm = cv
      cs = cvm*cvm
c
      do 100 i = 2, nxp1
        work1(i) = 0.5 * ( ex(i-1) + ex(i) )
  100 continue      
c
      work1(nxp2) = work1(2)
c
      do 110 i = 2, nxp1
        work2(i) = 0.5 * ( by(i+1) + by(i) )
  110 continue      
c
      work2(1) = work2(nxp1)
c
      n2=0
      do 210 k=1,ns 
         n1  = n2
         n2  = n1 + np(k)
c 
      do 200 m = n1+1, n2
c
           i   =  x(m) + 2.0
           sf2 = (x(m) + 2.0 - i)*qm(k)
           sf1 = qm(k) - sf2
           ih  =  x(m) + 1.5
           sh2 = (x(m) + 1.5 - ih)*qm(k)
           sh1 = qm(k) - sh2
           i1  = i + 1
           ih1 = ih+ 1
        ex1 = sf1*work1(i)  + sf2*work1(i1)
        ey1 = sf1*ey(i)     + sf2*ey(i1)
        ez1 = sh1*ez(ih)    + sh2*ez(ih1)
        bx1 = bx0*qm(k)
        by1 = sh1*work2(ih) + sh2*work2(ih1)
        bz1 = sh1*bz(ih)    + sh2*bz(ih1)
c
c --- Relativistic Modification  by Buneman's Method ----
c       g= 1.
        g = cvm/sqrt(cs-vx(m)*vx(m)-vy(m)*vy(m)-vz(m)*vz(m))
c       g1 = g
        ux = vx(m)*g + ex1
        uy = vy(m)*g + ey1
        uz = vz(m)*g + ez1
        g = cvm/sqrt(cs + ux*ux + uy*uy + uz*uz)
c       g2 = g
        bx1 = bx1*g
        by1 = by1*g
        bz1 = bz1*g
        boris = 2./(1. + bx1*bx1 + by1*by1 + bz1*bz1)
        uxt = ux + uy*bz1 - uz*by1
        uyt = uy + uz*bx1 - ux*bz1
        uzt = uz + ux*by1 - uy*bx1
        ux = ux + boris*(uyt*bz1 - uzt*by1) + ex1
        uy = uy + boris*(uzt*bx1 - uxt*bz1) + ey1
        uz = uz + boris*(uxt*by1 - uyt*bx1) + ez1
        g = cvm/sqrt(cs + ux*ux + uy*uy + uz*uz)
c       print*, g1, g2, g
        vx(m) = ux*g  
        vy(m) = uy*g  
        vz(m) = uz*g  
c
  200 continue
  210 continue
      return
      end
