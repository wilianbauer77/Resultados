c*******************************************************************
      subroutine rvlcty
      include "paramt.h"
      common /constc/ tcs, bx0, rho0, slx, nx, nxp1, nxp2, npt, ns
      common /fieldc/ ex(ix), ey(ix), ez(ix), by(ix), bz(ix),
     #                ajx(ix), ajy(ix), ajz(ix), rho(ix)
      common /prtclc/ x(in), vx(in), vy(in),vz(in)
      common /ptprmc/ wp(is), qm(is), q(is), vpe(is), vpa(is),
     #                vd(is), pch(is), np(is)
      common /work1c/ work1(ix),work2(ix)
c
      cs = tcs*0.5
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
         bx1 = bx0*qm(k)
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
        by1 = sh1*work2(ih) + sh2*work2(ih1)
        bz1 = sh1*bz(ih)    + sh2*bz(ih1)
c --- Relativistic Modification  by Romashin and Omura's Method ---
      vs = vx(m)*vx(m) + vy(m)*vy(m) + vz(m)*vz(m)
      gamma1 = 1.0 / sqrt( 1.0 - vs/cs )
      gamma2 = gamma1 
     &       + (ex1*vx(m) + ey1*vy(m) + ez1*vz(m))/ cs
      vxt = (vx(m) * gamma1 + ex1 + vy(m)*bz1 - vz(m)*by1)
     &       / gamma2
      vyt = (vy(m) * gamma1 + ey1 + vz(m)*bx1 - vx(m)*bz1)
     &       / gamma2
      vzt = (vz(m) * gamma1 + ez1 + vx(m)*by1 - vy(m)*bx1)
     &       / gamma2
c--- modification ---
c     vs = vxt*vxt + vyt*vyt + vzt*vzt
c     gamma2 = 1.0 / sqrt( 1.0 - vs/cs )
      gamma3 = gamma2**3
c-------------------
      bx1 = (bx1 - vyt*ez1 + vzt*ey1) / gamma2
      by1 = (by1 - vzt*ex1 + vxt*ez1) / gamma2
      bz1 = (bz1 - vxt*ey1 + vyt*ex1) / gamma2
      ex1 = ex1 / gamma3
      ey1 = ey1 / gamma3
      ez1 = ez1 / gamma3
c----------------------------------------------------------------- 
c
               boris = 2./(1. + bx1*bx1 + by1*by1 + bz1*bz1)
               vx(m) = vx(m) + ex1
               vy(m) = vy(m) + ey1
               vz(m) = vz(m) + ez1
               vxt    = vx(m) + vy(m)*bz1 - vz(m)*by1
               vyt    = vy(m) + vz(m)*bx1 - vx(m)*bz1
               vzt    = vz(m) + vx(m)*by1 - vy(m)*bx1
               vx(m) = vx(m) + boris*(vyt*bz1 - vzt*by1)
               vy(m) = vy(m) + boris*(vzt*bx1 - vxt*bz1)
               vz(m) = vz(m) + boris*(vxt*by1 - vyt*bx1)
               vx(m) = vx(m) + ex1
               vy(m) = vy(m) + ey1
               vz(m) = vz(m) + ez1
c
  200 continue
  210 continue
      return
      end
