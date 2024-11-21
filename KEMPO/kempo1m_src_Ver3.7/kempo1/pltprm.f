c****************************************************************
      subroutine pltprm
      include "paramt.h"
      common /constc/ tcs, bx0, rho0, slx, nx, nxp1, nxp2, npt, ns
      common /ptprmc/ wp(is), qm(is), q(is), vpe(is), vpa(is),
     &                vd(is), pch(is), np(is)
      common /timecm/ itime,ntime,iecrct,iwrite,jobnum
      common /diagcm/ iediag, ifdiag, ikdiag, ipdiag, isdiag, ivdiag,
     &                ieplot, ifplot, ikplot, ipplot, isplot, ivplot
      common /otherc/ vmin,vmax
      common /inputc/ dx, dt, cv, wc, angle
      common /ngyroc/ ngs, jperp, dpsi, ncycle
c----------------------------------------------------------
      h=0.6
      call newpen(3)
      call symbol(1.,25.,1.,'kempo1 parameters',0.,17)
      call prmplt(2.,   23.,h,0.,'nx    ',6,float(nx    ),0)
      call prmplt(-999.,23.,h,0.,'ntime ',6,float(ntime ),0)
      call prmplt(-999.,23.,h,0.,'iecrct',6,float(iecrct),0)
      call prmplt(-999.,23.,h,0.,'cv    ',6,cv           ,2)
      call prmplt(-999.,23.,h,0.,'dx    ',6,dx           ,2)
      call prmplt(-999.,23.,h,0.,'dt    ',6,dt           ,2)
      call prmplt(-999.,23.,h,0.,'wc    ',6,wc           ,2)
      call prmplt(-999.,23.,h,0.,'angle ',6,angle        ,2)
      call prmplt(  13.,23.,h,0.,'iediag',6,float(iediag),0)
      call prmplt(-999.,23.,h,0.,'isdiag',6,float(isdiag),0)
      call prmplt(-999.,23.,h,0.,'ifdiag',6,float(ifdiag),0)
      call prmplt(-999.,23.,h,0.,'ikdiag',6,float(ikdiag),0)
      call prmplt(-999.,23.,h,0.,'ivdiag',6,float(ivdiag),0)
      call prmplt(-999.,23.,h,0.,'ipdiag',6,float(ipdiag),0)
      call prmplt(-999.,23.,h,0.,'iwrite',6,float(iwrite),0)
      call prmplt(-999.,23.,h,0.,'jobnum',6,float(jobnum),0)
      call prmplt(  24.,23.,h,0.,'ieplot',6,float(ieplot),0)
      call prmplt(-999.,23.,h,0.,'isplot',6,float(isplot),0)
      call prmplt(-999.,23.,h,0.,'ifplot',6,float(ifplot),0)
      call prmplt(-999.,23.,h,0.,'ikplot',6,float(ikplot),0)
      call prmplt(-999.,23.,h,0.,'ivplot',6,float(ivplot),0)
      call prmplt(-999.,23.,h,0.,'ipplot',6,float(ipplot),0)
      call prmplt(-999.,23.,h,0.,'ngs   ',6,float(ngs   ),0)
      call prmplt(-999.,23.,h,0.,'jperp ',6,float(jperp ),0)
      call prmplt(-999.,23.,h,0.,'dpsi  ',6,dpsi         ,2)
      call prmplt(-999.,23.,h,0.,'ncycle',6,float(ncycle),0)
      do 10 i=1,ns
        x=2.+10.*(i-1)
        y=12.
        call symbol(x,y,h,'species',0.,7)
        call number(x+0.8*8,y,h,float(i),0.,-1)
        y=10.
        call prmplt(x    ,y,h,0.,'np ',3,float(np(i)),0)
        call prmplt(-999.,y,h,0.,'qm ',3,qm(i)       ,2)
        call prmplt(-999.,y,h,0.,'wp ',3,wp(i)       ,2)
        call prmplt(-999.,y,h,0.,'vd ',3,vd(i)       ,2)
        call prmplt(-999.,y,h,0.,'pch',3,pch(i)      ,2)
        call prmplt(-999.,y,h,0.,'vpa',3,vpa(i)      ,2)
        call prmplt(-999.,y,h,0.,'vpe',3,vpe(i)      ,2)
   10 continue
      call chart
      return
      end
