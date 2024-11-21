c**********************************************************************
      subroutine ksppl2
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
      dimension work3(ix) 
c
      if(itime==0) then
        open(80,file='kspec32.dat',status='unknown')
      endif
      rk=6.283185/slx*ikplot *rex
      do 10 i = 1, nx
          work1(i) = by(i+1) /reb
          work1(i) = bz(i+1) /reb
   10 continue
      call realft(work1,nx,1)
      call realft(work2,nx,1)
      fact = 2.0/float(nx) 
      j=2
      do 70 i = 3,nx-1,2
        work3(j) = sqrt( work1(i)**2 + work1(i+1)**2 +
     &                   work2(i)**2 + work2(i+1)**2 )*fact
        j = j + 1
   70 continue
      work3(1) = sqrt(work1(1)**2 + work2(1)**2)*fact*0.5
      work3(j) = sqrt(work1(2)**2 + work2(2)**2)*fact*0.5

      do 12 i = 2, 33
          write(80, *) work3(i)
   12 continue

      return
      end
