c*******************************************************************
      subroutine writer
      include "paramt.h"
      common /constc/ tcs, bx0, rho0, slx, nx, nxp1, nxp2, npt, ns
      common /fieldc/ ex(ix), ey(ix), ez(ix), by(ix), bz(ix),
     &                ajx(ix), ajy(ix), ajz(ix), rho(ix)
      common /prtclc/ x(in), vx(in), vy(in),vz(in)
      common /ptprmc/ wp(is), qm(is), q(is), vpe(is), vpa(is),
     &                vd(is), pch(is), np(is)
      common /ecrctc/ rkfact(ix)
      common /resclc/ rex, ret, rev, ree, reb, rej, rer, res
      common /inputc/ dx, dt, cv, wc, angle
      common /rotatc/ sinth, costh
      common /timecm/ itime,ntime,iecrct,iwrite,jobnum
      common /diagcm/ iediag, ifdiag, ikdiag, ipdiag, isdiag, ivdiag,
     &                ieplot, ifplot, ikplot, ipplot, isplot, ivplot
      common /otherc/ vmin,vmax
c
      if(jobnum==0) return
      ind=90
      jx = ix
      js = is
      jn = in
      open(ind,file='kempo1.cont',status='unknown',
     &       form='unformatted',access='sequential')
      write(ind) jx,js,jn,itime,ntime,iecrct,iwrite,jobnum
      write(ind) iediag,ifdiag,ikdiag,ipdiag,isdiag,ivdiag
      write(ind) ieplot,ifplot,ikplot,ipplot,isplot,ivplot
      write(ind) tcs,bx0,rho0,slx,nx,nxp1,nxp2,npt,ns
      write(ind) dx,dt,cv,wc,angle,sinth,costh,vmin,vmax
      write(ind) rex,ret,rev,ree,reb,rej,rer,res
      write(ind) wp,qm,q,vpe,vpa,vd,pch,np
      write(ind) (ex(i),     i = 1, nxp2)
      write(ind) (ey(i),     i = 1, nxp2)
      write(ind) (ez(i),     i = 1, nxp2)
      write(ind) (by(i),     i = 1, nxp2)
      write(ind) (bz(i),     i = 1, nxp2)
      write(ind) (ajx(i),    i = 1, nxp2)
      write(ind) (ajy(i),    i = 1, nxp2)
      write(ind) (ajz(i),    i = 1, nxp2)
      write(ind) (rho(i),    i = 1, nxp2)
      write(ind) (rkfact(i), i = 1, nxp2)
      i2 = 0
      do 10 j = 1, npt/nx
        i1 = i2 + 1
        i2 = i2 + nx
        write(ind) (x(i),  i = i1, i2)
        write(ind) (vx(i), i = i1, i2)
        write(ind) (vy(i), i = i1, i2)
        write(ind) (vz(i), i = i1, i2)
   10 continue
      if(i2<npt) then
        i1 = i2 + 1
        write(ind) ( x(i), i = i1, npt)
        write(ind) (vx(i), i = i1, npt)
        write(ind) (vy(i), i = i1, npt)
        write(ind) (vz(i), i = i1, npt)
      endif
      close(ind)
      return
      end
