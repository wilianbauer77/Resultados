C****************************************************************
      subroutine input
      include "paramt.h"
      common /constc/ tcs, bx0, rho0, slx, nx, nxp1, nxp2, npt, ns
      common /ptprmc/ wp(is), qm(is), q(is), vpe(is), vpa(is),
     #                vd(is), pch(is), np(is)
      common /timecm/ itime,ntime,iecrct,iwrite,jobnum
      common /diagcm/ iediag, ifdiag, ikdiag, ipdiag, isdiag, ivdiag,
     &                ieplot, ifplot, ikplot, ipplot, isplot, ivplot      
      common /otherc/ vmin,vmax
      common /inputc/ dx, dt, cv, wc, angle
c----------------------------------------------------------
      dx = 1.0
      dt = 0.025
      nx = 1024
      ntime = 1024
      iediag = 32
      isdiag = 32
      ifdiag = 99999
      ipdiag = ntime/2048
      ivdiag = 99999
      ikdiag = 99999
      ieplot = ntime/iediag
      ifplot = 2
      ikplot = nx/8
      isplot = 512
      ipplot = 1
      ivplot = 1
      vmin   = -40.0
      vmax   =  40.0
      cv     = 40
      wc     = -1.0
      angle  = 0.
      iecrct = 9999999
      iwrite = 8192*32
      jobnum = 2
c
      ns=3
      wp(1)  = sqrt(0.5)
      qm(1)  = -1.0
      vpe(1) = 1.0
      vpa(1) = 1.0
      vd(1)  = 0.0
      pch(1) = 0.0
      np(1)  = 4096
c
      wp(2)  = sqrt(0.5)
      qm(2)  = -1.0
      vpe(2) = 1.0
      vpa(2) = 1.0
      vd(2)  = 20.0
      pch(2) = 0.0
      np(2)  = 4096
c
      wp(3)  = 0.1
      qm(3)  = 0.01
      vpe(3) = 2.0
      vpa(3) = 2.0
      vd(3)  = 0.0
      pch(3) = 0.0
      np(3)  = 4096 
c---parameters for nongyrotropy
      ngs    = 0
      jperp  = 0
      dpsi   = 0.
      ncycle = 0
c
c----------------------------------------------------------
      return
      end
