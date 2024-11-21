c**********************************************************************
      subroutine curntv
      include "paramt.h"
      parameter(lvec=128)
      common /inputc/ dx, dt, cv, wc, angle
      common /constc/ tcs, bx0, rho0, slx, nx, nxp1, nxp2, npt, ns
      common /fieldc/ ex(ix), ey(ix), ez(ix), by(ix), bz(ix),
     &                ajx(ix), ajy(ix), ajz(ix), rho(ix)
      common /prtclc/ x(in), vx(in), vy(in),vz(in)
      common /ptprmc/ wp(is), qm(is), q(is), vpe(is), vpa(is),
     &                vd(is), pch(is), np(is)
      common /timecm/ itime,ntime,iecrct,iwrite,jobnum
      common /digcom/ omega
      common /resclc/ rex, ret, rev, ree, reb, rej, rer, res
      dimension wrk1(lvec,ix),wrk2(lvec,ix),wrk3(lvec,ix)
c
      do 100 i=1,nxp2
        ajx(i) = 0.0
        ajy(i) = 0.0
        ajz(i) = 0.0
  100 continue
      do 150 i=1,nxp2
      do 150 l=1,lvec
         wrk1(l,i) = 0.0
         wrk2(l,i) = 0.0
         wrk3(l,i) = 0.0
  150 continue
c
      n2 = 0
      do 210 k=1,ns 
        n1 = n2
        n2 = n1 + np(k)
        qh = q(k)*0.5
        do 200 ik = n1+1,n2,lvec
c$dir no_recurrence
        do 200 m  = ik,min(ik+lvec-1,n2)
          l = m - ik + 1
          ih = x(m) + 1.5
          s2 = (x(m) + 1.5 - ih)*q(k)
          s1 = q(k) - s2
          ih1 = ih + 1
          wrk2(l,ih )  = wrk2(l,ih ) + vy(m)*s1
          wrk2(l,ih1)  = wrk2(l,ih1) + vy(m)*s2
          wrk3(l,ih )  = wrk3(l,ih ) + vz(m)*s1
          wrk3(l,ih1)  = wrk3(l,ih1) + vz(m)*s2
c--------- charge conservation method ---------
          qhs = qh * sign(1.0, vx(m))
          avx=abs(vx(m))
          x1 = x(m) + 2.0 - avx
          x2 = x(m) + 2.0 + avx
          i1 = x1
          i2 = x2
          wrk1(l,i1) = wrk1(l,i1) + (i2 - x1)*qhs
          wrk1(l,i2) = wrk1(l,i2) + (x2 - i2)*qhs
c----------------------------------------------
  200   continue
  210 continue
c
      do 300 i=1,nxp2
      do 300 l=1,lvec
        ajx(i) = ajx(i) + wrk1(l,i)
        ajy(i) = ajy(i) + wrk2(l,i)
        ajz(i) = ajz(i) + wrk3(l,i)
  300 continue
c
      ajx(nxp1) = ajx(1) + ajx(nxp1)
      ajx(2)    = ajx(2) + ajx(nxp2)
      ajy(nxp1) = ajy(1) + ajy(nxp1)
      ajy(2)    = ajy(2) + ajy(nxp2)
      ajy(1)    = ajy(nxp1)
      ajz(nxp1) = ajz(1) + ajz(nxp1)
      ajz(2)    = ajz(2) + ajz(nxp2)
c
      do 350 i = nxp1, 2,-1
        ajy(i) = (ajy(i) + ajy(i-1))*0.5
  350 continue
c-------- cancel the uniform component ---------------
      juncan = 1
      if(juncan==1) then
      ajxu = 0.0
      ajyu = 0.0
      ajzu = 0.0
      do 400 i = 2,nxp1
        ajxu = ajxu + ajx(i)
        ajyu = ajyu + ajy(i)
        ajzu = ajzu + ajz(i)
  400 continue
      ajxu = ajxu/float(nx)
      ajyu = ajyu/float(nx)
      ajzu = ajzu/float(nx)
      do 500 i = 2,nxp1
        ajx(i) = ajx(i) - ajxu
        ajy(i) = ajy(i) - ajyu
        ajz(i) = ajz(i) - ajzu
  500 continue
      endif
c------ external current source ------
c     if(itime<1024) then
c       omega = 0.5
c     else 
c       omega = 0.5+0.001*(itime - 1024)
c     endif
c     i = nx/2+1
c     ajy(i) = ajy(i) + 80.0*rej*sin( omega*dt*itime )
c---addition  by libo liu
c     wf = 2.0
c     cus = 20.*rej
c     ajz(512) = ajz(512) + cus*sin(wf*itime*dt)
c
      return
      end
