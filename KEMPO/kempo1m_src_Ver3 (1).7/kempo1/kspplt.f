c**********************************************************************
      subroutine kspplt
      include "paramt.h"
      common /constc/ tcs, bx0, rho0, slx, nx, nxp1, nxp2, npt, ns
      common /fieldc/ ex(ix), ey(ix), ez(ix), by(ix), bz(ix),
     &                ajx(ix), ajy(ix), ajz(ix), rho(ix)
      common /timecm/ itime,ntime,iecrct,iwrite,jobnum
      common /diagcm/ iediag, ifdiag, ikdiag, ipdiag, isdiag, ivdiag,
     &                ieplot, ifplot, ikplot, ipplot, isplot, ivplot
      common /resclc/ rex, ret, rev, ree, reb, rej, rer, res
      common /inputc/ dx, dt, cv, wc, angle
      common /work1c/ work1(ix),work2(ix)
c
      rk=6.283185/slx*ikplot *rex
      do 10 i = 1, nx
          work1(i) = ex(i+1) /ree
   10 continue
      call realft(work1,nx,1)
      fact = 2.0/float(nx) 
      j=2
      do 70 i = 3,nx-1,2
        work2(j) = sqrt( work1(i)**2 + work1(i+1)**2 )*fact
        j = j + 1
   70 continue
      work2(1) = abs(work1(1))*fact*0.5
      work2(j) = abs(work1(2))*fact*0.5
      call qlook(work2,ikplot+1,7.,15.,10.,10.,0.,rk,'k',1,'ex',2)
      do 12 i = 1, nx
          work1(i) = ey(i+1) /ree
   12 continue
      call realft(work1,nx,1)
      fact = 2.0/float(nx) 
      j=2
      do 72 i = 3,nx-1,2
        work2(j) = sqrt( work1(i)**2 + work1(i+1)**2 )*fact
        j = j + 1
   72 continue
      work2(1) = abs(work1(1))*fact*0.5
      work2(j) = abs(work1(2))*fact*0.5
      call qlook(work2,ikplot+1,7.,2.,10.,10.,0.,rk,'k',1,'ey',2)
      do 14 i = 1, nx
          work1(i) = ez(i+1) /ree
   14 continue
      call realft(work1,nx,1)
      fact = 2.0/float(nx) 
      j=2
      do 74 i = 3,nx-1,2
        work2(j) = sqrt( work1(i)**2 + work1(i+1)**2 )*fact
        j = j + 1
   74 continue
      work2(1) = abs(work1(1))*fact*0.5
      work2(j) = abs(work1(2))*fact*0.5
      call qlook(work2,ikplot+1,22.,15.,10.,10.,0.,rk,'k',1,'ez',2)
c
      do 20 i = 1, nx
          work1(i) = bz(i+1) /reb
   20 continue
      call realft(work1,nx,1)
      fact = 2.0/float(nx) 
      j=2
      do 80 i = 3,nx-1,2
        work2(j) = sqrt( work1(i)**2 + work1(i+1)**2 )*fact
        j = j + 1
   80 continue
      work2(1) = abs(work1(1))*fact*0.5
      work2(j) = abs(work1(2))*fact*0.5
      rk=6.283185/slx*ikplot *rex
      call qlook(work2,ikplot+1,22.,2.,10.,10.,0.,rk,'k',1,'bz',2)
      call prmplt(25.,25.8,0.6,0.,'time',4,dt*itime,2)
      call chart
      return
      end
