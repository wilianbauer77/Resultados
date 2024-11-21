***********************************************************************
      subroutine renrgy
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
      cs=tcs*0.5
      te=0.0
      tb=0.0 
      n2=0
      do 60 k=1,ns
      n1=n2
      n2=n1+np(k)
      rm=q(k)/qm(k)
      gm   = 0.0
      gmvx =0.0
      gmvy =0.0
      gmvz =0.0
      rn = np(k)
      do 10 i=n1+1,n2
        vs = vx(i)*vx(i) + vy(i)*vy(i) + vz(i)*vz(i)
        gamma = 1.0 / sqrt(1.0 - vs /cs)
        gm = gm + gamma
        gmvx = gmvx + gamma*vx(i)
        gmvy = gmvy + gamma*vy(i)
        gmvz = gmvz + gamma*vz(i)
   10 continue
      wkx(ic,k) = rm*cs*(gm - rn)/slx/res
      wdx(ic,k) = gmvx/gm
      wdy(ic,k) = gmvy/gm
      wdz(ic,k) = gmvz/gm
      vds = (gmvx*gmvx + gmvy*gmvy + gmvz*gmvz)/(gm*gm)
      gammad = 1.0 / sqrt(1.0 - vds /cs)
      wky(ic,k) = rm*cs*(gammad - 1.0)*rn/slx/res
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
        wk3(ic) = wk3(ic) + wkx(ic,k)
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
        wk1(i) = wkx(i,k)
        wk2(i) = wky(i,k)
        wk3(i) = wk1(i) - wk2(i)
        tpara  = wdx(i,k)**2
        if(tpara<1.e-7) tpara = 1.e-7
        wk4(i) = 0.5*(wdy(i,k)**2 + wdz(i,k)**2)/tpara
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
