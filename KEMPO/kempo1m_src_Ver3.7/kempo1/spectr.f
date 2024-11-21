c**********************************************************************
      subroutine spectr
      include "paramt.h"
      parameter(imax=64,jmax=2048,icomp=5)
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
      common /rotatc/ sinth, costh
      common /work1c/ work1(ix),work2(ix)
      common /work2c/ wk2d(imax,jmax,icomp)
      dimension workt(jmax)
      character*2 comp(5)
      dimension xp(4),yp(4)
      save comp,xp,yp,ic 
      data comp/'ex','ey','ez','by','bz'/
      data xp / 7.,22., 7.,22./
      data yp /16.,16., 2., 2./
      data ic/0/
c
      ncomp = 5
      minmod = 1 
      maxmod = 16
      ikp = 32
      iwp = 128
      ifb = 0
c
      if(maxmod.ge.imax/2) maxmod = imax/2-1
      if(ic==0) t1=dt*itime
      nxrd=imax
      if(nx<nxrd) nxrd=nx
      ic=ic+1
      ncp=1
      do 10 i=1,nx
        work1(i)=ex(i+1) /ree
  10  continue
      call skfft(work1,nx,wk2d(1,ic,ncp),nxrd)
      if(ncp==ncomp) goto 60
      ncp=ncp+1
      do 20 i=1,nx
        work1(i)=ey(i+1) /ree
  20  continue
      call skfft(work1,nx,wk2d(1,ic,ncp),nxrd)
      if(ncp==ncomp) goto 60
      ncp=ncp+1
      do 30 i=1,nx
        work1(i)=ez(i+1) /ree
  30  continue
      call skfft(work1,nx,wk2d(1,ic,ncp),nxrd)
      if(ncp==ncomp) goto 60
      ncp=ncp+1
      by0  = wc/qm(1)*sinth
      do 40 i=1,nx
        work1(i)=(by(i+1)-by0) /reb
  40  continue
      call skfft(work1,nx,wk2d(1,ic,ncp),nxrd)
      if(ncp==ncomp) goto 60
      ncp=ncp+1
      do 50 i=1,nx
        work1(i)=bz(i+1) /reb
  50  continue
      call skfft(work1,nx,wk2d(1,ic,ncp),nxrd)
  60  if(ic/=isplot.and.ic/=jmax) return
      dtime=dt*isdiag*ic
      t2=dt*itime
c
      do 80 ncp = 1, ncomp
        i = (minmod-1)*2 + 1
        do 75 k = minmod, maxmod
          i = i + 2
          m = mod(k-minmod,4) + 1
          do 70 j=1,ic
            workt(j) = sqrt(wk2d(i,j,ncp)**2 + wk2d(i+1,j,ncp)**2)
   70     continue
          call qlook(workt,ic,xp(m),yp(m),10.,10.,
     &                t1,t2,'time',4,comp(ncp),2)
          call newpen(5)
          call prmplt(xp(m)+6.,yp(m)+10.5,0.5,0.,
     &                'mode',4,float(k),0)
          if(m==4) call chart
   75   continue
   80 continue
      if(m/=4) call chart
c
      if(ikp>nxrd/2) ikp = nxrd/2
      if(iwp>ic/2) iwp = ic/2
      do 90 ncp=1,ncomp
        call newpen(5)
        call symbol(1.,25.5,1.0,comp(ncp),0.,2)
        call newpen(3)
        call wkfft(wk2d(1,1,ncp),imax,jmax,nxrd,ic,work1,workt,1)
        call wkplot(wk2d(1,1,ncp),imax,jmax,nxrd,ic,work1,workt,
     &               8.,5.,20.,20.,slx/rex,dtime,ikp,iwp,1)
        call symbol(29.,24.,0.8,'forward',0.,7)
        call prmplt(29.,22.,0.8,0.,'t1',2,t1,2)
        call prmplt(29.,20.,0.8,0.,'t2',2,t2,2)
        call chart
        call wkplot(wk2d(1,1,ncp),imax,jmax,nxrd,ic,work1,workt,
     &                 8.,5.,20.,20.,slx/rex,dtime,ikp,iwp,-1)
        call prmplt(29.,22.,0.8,0.,'t1',2,t1,2)
        call prmplt(29.,20.,0.8,0.,'t2',2,t2,2)
        call symbol(29.,24.,0.8,'backward',0.,8)
        call chart
        call wkplot(wk2d(1,1,ncp),imax,jmax,nxrd,ic,work1,workt,
     &                 8.,5.,20.,20.,slx/rex,dtime,ikp,iwp,0)
        call prmplt(29.,22.,0.8,0.,'t1',2,t1,2)
        call prmplt(29.,20.,0.8,0.,'t2',2,t2,2)
        call symbol(29.,24.,0.8,'total',0.,5)
        call chart
  90  continue
      ic = 0
      return
      end
