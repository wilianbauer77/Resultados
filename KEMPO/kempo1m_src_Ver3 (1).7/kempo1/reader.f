c*******************************************************************
      subroutine reader
      Include "paramt.h"
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
      ind=90
      open(ind,file='kempo1.cont',status='old',
     &       form='unformatted',access='sequential')
      read(ind) jx,js,jn,itime,jtime,jecrct,jwrite,jobnum
      read(ind) jediag,jfdiag,jkdiag,jpdiag,jsdiag,jvdiag
      read(ind) jeplot,jfplot,jkplot,jpplot,jsplot,jvplot
      read(ind) tcs,bx0,rho0,slx,nx,nxp1,nxp2,npt,ns
      read(ind) dx,dt,cv,wc,angle,sinth,costh,vmin,vmax
      read(ind) rex,ret,rev,ree,reb,rej,rer,res
      read(ind) wp,qm,q,vpe,vpa,vd,pch,np
      if(nxp2>ix)  go to 20
      if(ns  >is)  go to 20
      if(npt >in)  go to 20
      read(ind) (ex(i),     i = 1, nxp2)
      read(ind) (ey(i),     i = 1, nxp2)
      read(ind) (ez(i),     i = 1, nxp2)
      read(ind) (by(i),     i = 1, nxp2)
      read(ind) (bz(i),     i = 1, nxp2)
      read(ind) (ajx(i),    i = 1, nxp2)
      read(ind) (ajy(i),    i = 1, nxp2)
      read(ind) (ajz(i),    i = 1, nxp2)
      read(ind) (rho(i),    i = 1, nxp2)
      read(ind) (rkfact(i), i = 1, nxp2)
      i2 = 0
      do 10 j = 1, npt/nx
        i1 = i2 + 1
        i2 = i2 + nx
        read(ind) (x(i),  i = i1, i2)
        read(ind) (vx(i), i = i1, i2)
        read(ind) (vy(i), i = i1, i2)
        read(ind) (vz(i), i = i1, i2)
   10 continue
      if(i2<npt) then
        i1 = i2 + 1
        read(ind) ( x(i), i = i1, npt)
        read(ind) (vx(i), i = i1, npt)
        read(ind) (vy(i), i = i1, npt)
        read(ind) (vz(i), i = i1, npt)
      endif
      jobnum=jobnum+1
      close(ind)
      return
  20  continue
      close(ind)
      write(0,*) 'inconsistent parameters : ix,is,in'
      write(0,*) 'ix=',nxp2,'is=',ns,'in=',npt
      stop
      end
