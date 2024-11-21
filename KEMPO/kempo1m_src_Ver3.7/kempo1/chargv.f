c**********************************************************************
      subroutine chargv
      include "paramt.h"
      parameter(lvec=128)
      dimension wrk1(lvec,ix)
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
      do 150 i=1,nxp2
      do 150 k = 1,lvec
          wrk1(k,i) = 0.0
  150 continue
c     call szero(lvec*nxp2,wrk1,1) 
c
      n2 = 0
      do 210 k=1,ns 
        n1 = n2
        n2 = n1 + np(k)
        do 200 ik = n1+1,n2,lvec
c$dir no_recurrence
        do 200 m  = ik,min(ik+lvec-1,n2)
           l  = m - ik + 1
           i  =  x(m)+ 2.0
           s2 = (x(m)+ 2.0 - i)*q(k)
           s1 = q(k) - s2
           wrk1(l,i )    = wrk1(l,i )   + s1
           wrk1(l,i+1 )  = wrk1(l,i+1 ) + s2
  200   continue
  210 continue
      do 300 l=1,lvec
      do 300 i=1,nxp2
        rho(i) = rho(i) + wrk1(l,i)
  300 continue
      rho(2)    = rho(2) + rho(nxp2) - rho0
      rho(1)    = rho(nxp1)
      rho(nxp2) = rho(2)
      return
      end
