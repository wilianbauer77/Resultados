c**********************************************************************
      subroutine phsplt
      include "paramt.h"
c
      common /constc/ tcs, bx0, rho0, slx, nx, nxp1, nxp2, npt, ns
      common /prtclc/ x(in), vx(in), vy(in),vz(in)
      common /ptprmc/ wp(is), qm(is), q(is), vpe(is), vpa(is),
     &                vd(is), pch(is), np(is)
      common /timecm/ itime,ntime,iecrct,iwrite,jobnum
      common /diagcm/ iediag, ifdiag, ikdiag, ipdiag, isdiag, ivdiag,
     &                ieplot, ifplot, ikplot, ipplot, isplot, ivplot
      common /otherc/ vmin,vmax
      common /resclc/ rex, ret, rev, ree, reb, rej, rer, res
      common /inputc/ dx, dt, cv, wc, angle
      common /rotatc/ sinth, costh
      dimension ipen(7)
      data ipen/5,3,6,2,4,7,1/
c
      ipmax = 20000
c
      v1=vmin
      v2=vmax
      if(vmin==vmax)  call maxmin(vx,npt,v1,v2)
      if(mod(ipplot,2)==1) then
        xfact=20.0/slx
        yfact=20.0/(v2-v1)
        call newpen(3)
        call prmplt(16.,25.,0.8,0.,'time',4,itime*dt,2)
        call newpen(4)
        call xaxis1(8.,4.,20.,10,2,0.5,0.8,0.,slx/rex,3,'x',1)
        call xaxis1(8.,24.,20.,10,2,-0.5,0.0,0.,slx/rex,2,'x',1)
        call yaxis1(8.,4.,20.,10,2,0.5,0.8,v1/rev,v2/rev,3,'vx',2)
        call yaxis1(28.,4.,20.,10,2,-0.5,0.0,v1/rev,v2/rev,3,'vx',2)
        n2=0
        do 10 k=1,ns
          n1 = n2
          n2 = n1 + np(k)
          n3 = np(k)/ipmax + 1
          call newpen(ipen(k))
          do 20 i=n1+1,n2,n3
            xx = 8.0 + x(i)*xfact
            yy = 4.0 + (vx(i)-v1)*yfact
            call plot(xx,yy,3)
            call plot(xx+0.05,yy,2)
   20     continue
   10   continue
      call chart
      end if
      if(mod(ipplot,4).ge.2) then
        xfact=20.0/(v2-v1)
        yfact=20.0/(v2-v1)
        call newpen(5)
        call prmplt(26.,25.,0.8,0.,'time',4,itime*dt,2)
        call newpen(4)
        call xaxis1(10.,4.,20.,10,2,0.5,0.7,v1/rev,v2/rev,3,
     &              'v-perp/xy',9)
        call xaxis1(10.,24.,20.,10,2,-0.5,0.0,v1/rev,v2/rev,3,'p',1)
        call yaxis1(10.,4.,20.,10,2,0.5,0.7,v1/rev,v2/rev,3,'vz',2)
        call yaxis1(30.,4.,20.,10,2,-0.5,0.0,v1/rev,v2/rev,3,'vz',2)
        n2=0
        do 30 k=1,ns
          n1 = n2
          n2 = n1 + np(k)
          n3 = np(k)/ipmax + 1
          call newpen(ipen(k))
          do 40 i=n1+1,n2,n3
            vpx = costh*vx(i) + sinth*vy(i)
            vpy =-sinth*vx(i) + costh*vy(i) 
            xx = 10.0 + (vpy-v1)*xfact
            yy = 4.0 + (vz(i)-v1)*yfact
            call plot(xx,yy,3)
            call plot(xx+0.05,yy,2)
   40     continue
   30   continue
      call chart
      end if
      if(mod(ipplot,8).ge.4) then
        xfact=28.0/(v2-v1)
        yfact=14.0/v2
        call newpen(5)
        call prmplt(24.,23.,0.8,0.,'time',4,itime*dt,2)
        call newpen(4)
        call xaxis1(5.,8.,28.,20,4,0.5,0.7,v1/rev,v2/rev,3,'v-para',6)
        call xaxis1(5.,22.,28.,20,4,-0.5,0.0,v1/rev,v2/rev,3,'v-para',6)
        call yaxis1(5.,8.,14.,10,2,0.5,0.7,0.,v2/rev,3,'v-perp',6)
        call yaxis1(33.,8.,14.,10,2,-0.5,0.0,0.,v2/rev,3,'v-perp',6)
        n2=0
        do 50 k=1,ns
          n1 = n2
          n2 = n1 + np(k)
          n3 = np(k)/ipmax + 1
          call newpen(ipen(k))
          do 60 i=n1+1,n2,n3
            vpx = costh*vx(i) + sinth*vy(i)
            vpy =-sinth*vx(i) + costh*vy(i) 
            vperp = sqrt(vpy*vpy+vz(i)*vz(i))
            xx = 5.0 + (vpx-v1)*xfact
            yy = 8.0 + vperp*yfact
            call plot(xx,yy,3)
            call plot(xx+0.05,yy,2)
   60     continue
   50   continue
      call chart
      end if
      return
      end
