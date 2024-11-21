c**********************************************************************
      subroutine vdsplt
      include "paramt.h"
      common /constc/ tcs, bx0, rho0, slx, nx, nxp1, nxp2, npt, ns
      common /prtclc/ x(in), vx(in), vy(in),vz(in)
      common /ptprmc/ wp(is), qm(is), q(is), vpe(is), vpa(is),
     &                vd(is), pch(is), np(is)
      common /diagcm/ iediag, ifdiag, ikdiag, ipdiag, isdiag, ivdiag,
     &                ieplot, ifplot, ikplot, ipplot, isplot, ivplot
      common /timecm/ itime,ntime,iecrct,iwrite,jobnum
      common /inputc/ dx, dt, cv, wc, angle
      common /otherc/ vmin,vmax
      common /resclc/ rex, ret, rev, ree, reb, rej, rer, res
      common /work1c/ work1(ix),work2(ix)
c
      nv = 101
      if(nv>ix) nv = ix
      n2 = 0
      do 10 k=1,ns 
        n1 = n2
        n2 = n1 + np(k)
        v1 = vmin
        v2 = vmax
        if(vmin==vmax) call maxmin(vx(n1+1),np(k),v1,v2)
        dvi = float(nv-1)/(v2-v1)
      do 20 i = 1,nv
         work1(i) = 0.
   20    continue
      do 30 m = n1+1, n2
        if(vx(m)<v1.or.vx(m).ge.v2) go to 30
          rv = (vx(m)-v1)*dvi + 1.0
          i  = rv
          s2 = rv - i
          s1 = 1.0 - s2
          work1(i)   = work1(i)   + s1
          work1(i+1) = work1(i+1) + s2
   30 continue
      do 40 i = 1,nv
          work1(i) = work1(i)*dvi*rev/float(np(k))
   40 continue
      call qlook(work1,nv,7.,15.5,10.,10.,v1/rev,v2/rev,
     &                  'vx',20002,'f(vx)',5)
      if(vmin==vmax) call maxmin(vy(n1+1),np(k),v1,v2)
      dvi = float(nv-1)/(v2-v1)
      do 22 i = 1,nv
         work1(i) = 0.
   22    continue
      do 32 m = n1+1, n2
        if(vy(m)<v1.or.vy(m).ge.v2) go to 32
          rv = (vy(m)-v1)*dvi + 1.0
          i  = rv
          s2 = rv - i
          s1 = 1.0 - s2
          work1(i)   = work1(i)   + s1
          work1(i+1) = work1(i+1) + s2
   32 continue
      do 42 i = 1,nv
          work1(i) = work1(i)*dvi*rev/float(np(k))
   42 continue
      call qlook(work1,nv,7.,2.,10.,10.,v1/rev,v2/rev,
     &                  'vy',20002,'f(vy)',5)
      if(vmin==vmax) call maxmin(vz(n1+1),np(k),v1,v2)
      dvi = float(nv-1)/(v2-v1)
      do 24 i = 1,nv
         work1(i) = 0.
   24    continue
      do 34 m = n1+1, n2
        if(vz(m)<v1.or.vz(m).ge.v2) go to 34
          rv = (vz(m)-v1)*dvi + 1.0
          i  = rv
          s2 = rv - i
          s1 = 1.0 - s2
          work1(i)   = work1(i)   + s1
          work1(i+1) = work1(i+1) + s2
   34 continue
      do 44 i = 1,nv
          work1(i) = work1(i)*dvi*rev/float(np(k))
   44 continue
      call qlook(work1,nv,22.,2.,10.,10.,v1/rev,v2/rev,
     &                  'vz',20002,'f(vz)',5)
      call newpen(5)
      call symbol(24.,24.,0.8,'species',0.,7)
      call number(30.4,24.,0.8,float(k),0.,-1)
      call prmplt(24.,22.,0.8,0.,'time',4,itime*dt,2)
      call chart
   10 continue
        return
        end
