c**********************************************************************
      subroutine fldplt
      include "paramt.h"
      common /constc/ tcs, bx0, rho0, slx, nx, nxp1, nxp2, npt, ns
      common /fieldc/ ex(ix), ey(ix), ez(ix), by(ix), bz(ix),
     &                ajx(ix), ajy(ix), ajz(ix), rho(ix)
      common /timecm/ itime,ntime,iecrct,iwrite,jobnum
      common /resclc/ rex, ret, rev, ree, reb, rej, rer, res
      common /inputc/ dx, dt, cv, wc, angle
      common /diagcm/ iediag, ifdiag, ikdiag, ipdiag, isdiag, ivdiag,
     &                ieplot, ifplot, ikplot, ipplot, isplot, ivplot
      common /work1c/ work1(ix),work2(ix)
c
      if(mod(ifplot,2)==1)  then
      do 10 i=2,nxp1
        work1(i)=ex(i)/ree
   10 continue
      call qlook(work1(2),nx,6.,19.,10.,5.,0.,slx/rex,'x',1,'ex',2)
c     call qlook2(work2(2),nx,0)
      do 20 i=2,nxp1
        work1(i)=ey(i)/ree
   20 continue
      call qlook(work1(2),nx,6.,11.,10.,5.,0.,slx/rex,'x',1,'ey',2)
      do 30 i=2,nxp1
        work1(i)=ez(i)/ree
   30 continue
      call qlook(work1(2),nx,6.,3.,10.,5.,0.,slx/rex,'x',1,'ez',2)
      do 40 i=2,nxp1
        work1(i)=by(i)/reb
   40 continue
      call qlook(work1(2),nx,23.,11.,10.,5.,0.,slx/rex,'x',1,'by',2)
      do 50 i=2,nxp1
        work1(i)=bz(i)/reb
   50 continue
      call qlook(work1(2),nx,23.,3.,10.,5.,0.,slx/rex,'x',1,'bz',2)
      do 60 i=2,nxp1
        work1(i)=rho(i)/rer
   60 continue
      call qlook(work1(2),nx,23.,19.,10.,5.,0.,slx/rex,'x',1,'rho',3)
      work2(1)=0.0
      work2(2)=0.0
      call newpen(1)
      call qlook2(work2,2,0)
      call newpen(3)
      call prmplt(25.,25.,0.7,0.,'time',4,itime*dt,2)
      call chart
      endif
c
      if(mod(ifplot,4).ge.2)  then
      call qlkmd2(0.2,0.5)
      do 65 i=2,nxp1
        work1(i)=rho(i)/rer
   65 continue
c     work2(1)= 1.0
c     work2(2)= -1.0
c     call qlook(work2(1),2,8.,3.,20.,5.,0.,slx/rex,'x',1+11000,'rho',3)
c     call qlook2(work1(2),nx,1)
      call qlook(work1(2),nx,8.,3.,20.,5.,0.,slx/rex,'x',1,'rho',3)
      work2(1)=0.0
      work2(2)=0.0
      call newpen(1)
      call qlook2(work2,2,0)
      do 15 i=2,nxp1
        work1(i)=ex(i)/ree
   15 continue
c     work2(1)= 2.0
c     work2(2)= -2.0
c     call qlook(work2(1),2,8.,11.,20.,5.,0.,slx/rex,'x',1+11000,'ex',2)
c     call qlook2(work1(2),nx,1)
      call qlook(work1(2),nx,8.,11.,20.,5.,0.,slx/rex,'x',1,'ex',2)
      work2(2) = 0.0
      do 70 i=2,nx
        work2(i+1)=work2(i) - ex(i)
   70 continue
      phi0 = work2(2)
      do 80 i=3,nxp1
        phi0 = phi0 + work2(i)
   80 continue
      phi0 = phi0/float(nx)
      do 90 i=2,nxp1
        work1(i) = (work2(i) - phi0)/(ree*rex)
   90 continue
c     work2(1)= 8.0
c     work2(2)= -8.0
c     call qlook(work2(1),2,8.,19.,20.,5.,0.,slx/rex,'x',1+11000,
c    &                    'potential',9)
c     call qlook2(work1(2),nx,1)
      call qlook(work1(2),nx,8.,19.,20.,5.,0.,slx/rex,'x',1,
     &                    'potential',9)
      work2(1)=0.0
      work2(2)=0.0
      call newpen(1)
      call qlook2(work2,2,0)
      call newpen(3)
      call prmplt(25.,25.,0.7,0.,'time',4,itime*dt,2)
      call chart
      call qlkmd2(0.0,0.0)
      endif
      return
      end
