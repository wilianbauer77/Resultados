c**********************************************************************
      subroutine energy
      include "paramt.h"
      parameter(iw=8193)
      common /constc/ tcs, bx0, rho0, slx, nx, nxp1, nxp2, npt, ns
      common /fieldc/ ex(ix), ey(ix), ez(ix), by(ix), bz(ix),
     &                ajx(ix), ajy(ix), ajz(ix), rho(ix)
      common /prtclc/ x(in), vx(in), vy(in),vz(in)
      common /ptprmc/ wp(is), qm(is), q(is), vpe(is), vpa(is),
     &                vd(is), pch(is), np(is)
      common /timecm/ itime,ntime,iecrct,iwrite,jobnum
      common /diagcm/ iediag, ifdiag, ikdiag, ipdiag, isdiag, ivdiag,
     &                ieplot, ifplot, ikplot, ipplot, isplot, ivplot
      common /resclc/ rex, ret, rev, ree, reb, rej, rer, res
      common /inputc/ dx, dt, cv, wc, angle
      common /work3c/ wkx(iw,is), wky(iw,is), wkz(iw,is),
     &                wdx(iw,is), wdy(iw,is), wdz(iw,is),
     &                wk1(iw), wk2(iw), wk3(iw), wk4(iw)
      common /rotatc/ sinth, costh
      save ic
      data ic/0/
c
      if(ic==0) t1=dt*itime
      ic=ic+1
      te=0.0
      tb=0.0 
      n2=0
      do 60 k=1,ns
      n1=n2
      n2=n1+np(k)
      rm=q(k)/qm(k)
      tkx=0.0
      tky=0.0
      tkz=0.0
      tdx=0.0
      tdy=0.0
      tdz=0.0
      do 10 i=n1+1,n2
        tkx = tkx + vx(i)*vx(i)
        tky = tky + vy(i)*vy(i)
        tkz = tkz + vz(i)*vz(i)
        tdx = tdx + vx(i)
        tdy = tdy + vy(i)
        tdz = tdz + vz(i)
   10 continue
      wkx(ic,k) = 0.5*rm*tkx/slx/res
      wky(ic,k) = 0.5*rm*tky/slx/res
      wkz(ic,k) = 0.5*rm*tkz/slx/res
      wdx(ic,k) = 0.5*rm*tdx*tdx/float(np(k))/slx/res
      wdy(ic,k) = 0.5*rm*tdy*tdy/float(np(k))/slx/res
      wdz(ic,k) = 0.5*rm*tdz*tdz/float(np(k))/slx/res
   60 continue
      do 20 i=2,nxp1
        te = te + ex(i)*ex(i) + ey(i)*ey(i) + ez(i)*ez(i)
   20 continue
      by0  = wc/qm(1)*sinth
      do 30 i=2,nxp1
        tb = tb + (by(i) - by0)**2 + bz(i)**2
   30 continue
      wk1(ic) = 0.5*te/float(nx) /res
      wk2(ic) = 0.25*tcs*tb/float(nx) /res
      wk3(ic) = 0.
      do 40 k=1,ns
        wk3(ic) = wk3(ic) + wkx(ic,k) + wky(ic,k) + wkz(ic,k)
   40 continue
      wk4(ic) = wk1(ic) + wk2(ic) + wk3(ic)
c
      if(ic==ieplot.or.ic==iw) then
      t2=dt*itime
      nt = 20004
      call newpen(5)
      call symbol(0.5,25.8,0.7,'energy',0.,6)
      call qlook(wk1,ic, 7., 2.,10.,10.,t1,t2,'time',nt,'electric',8)
      call qlook(wk2,ic,22., 2.,10.,10.,t1,t2,'time',nt,'magnetic',8)
      call qlook(wk3,ic, 7.,15.,10.,10.,t1,t2,'time',nt,'kinetic',7)
      call qlook(wk4,ic,22.,15.,10.,10.,t1,t2,'time',nt,'total',5)
      call chart
      wk11=wk1(1)
      wk21=wk2(1)
      wk31=wk3(1)
      wk41=wk4(1)
c
      open(50,file='efeng.dat',status='unknown')
      do 91 i = 1, ic
      write(50,*) wk1(i)
   91 continue
      close(50)
      open(51,file='bfeng.dat',status='unknown')
      do 92 i = 1, ic
      write(51,*) wk2(i)
   92 continue
      close(51)
      open(52,file='kieng.dat',status='unknown')
      do 93 i = 1, ic
      write(52,*) wk3(i)
   93 continue
      close(52)
c

      do 50 i=1,ic
        wk1(i) = wk1(i) - wk11
        wk2(i) = wk2(i) - wk21
        wk3(i) = wk3(i) - wk31
        wk4(i) = wk4(i) - wk41
   50 continue
      nt = 4
      call newpen(5)
      call symbol(0.5,25.8,0.7,'energy',0.,6)
      call symbol(1.,12.,0.7,'variation',90.,9)
      call qlook(wk1,ic, 7., 2.,10.,10.,t1,t2,'time',nt,'electric',8)
      call prmplt(12.,12.3,0.45,0.,'E0',2,wk11,3)
      call qlook(wk2,ic,22., 2.,10.,10.,t1,t2,'time',nt,'magnetic',8)
      call prmplt(27.,12.3,0.45,0.,'M0',2,wk21,3)
      call qlook(wk3,ic, 7.,15.,10.,10.,t1,t2,'time',nt,'kinetic',7)
      call prmplt(12.,25.3,0.45,0.,'K0',2,wk31,3)
      call qlook(wk4,ic,22.,15.,10.,10.,t1,t2,'time',nt,'total',5)
      call prmplt(27.,25.3,0.45,0.,'T0',2,wk41,3)
      call chart
      do 70 k=1,ns
      do 80 i=1,ic
        wk1(i) = wkx(i,k) + wky(i,k) + wkz(i,k)
        wk2(i) = wdx(i,k) + wdy(i,k) + wdz(i,k)
        wk3(i) = wk1(i) - wk2(i)
        tpara  = wkx(i,k) - wdx(i,k)
        if(tpara<1.e-7) tpara = 1.e-7
        wk4(i) = 0.5*(wky(i,k)+wkz(i,k)-wdy(i,k)-wdz(i,k))/tpara
   80 continue
      call newpen(5)
      call symbol(0.5,25.8,0.7,'energy',0.,6)
      call newpen(3)
      call symbol(24.,25.7,0.8,'species',0.0,7)
      call number(31.,25.7,0.8,float(k),0.0,-1)
      call qlook(wk3,ic, 7., 2.,10.,10.,t1,t2,'time',nt,'thermal',7)
      call qlook(wk4,ic,22., 2.,10.,10.,t1,t2,
     &          'time',nt,'anisotropy',10)
      call qlook(wk2,ic, 7.,15.,10.,10.,t1,t2,'time',nt,'drift',5)
      call qlook(wk1,ic,22.,15.,10.,10.,t1,t2,'time',nt,'total',5)
      call chart
   70 continue
      ic=0
      end if
c
      return
      end
