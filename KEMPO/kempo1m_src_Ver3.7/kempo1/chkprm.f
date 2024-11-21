c****************************************************************
      subroutine chkprm
      include "paramt.h"
      common /constc/ tcs, bx0, rho0, slx, nx, nxp1, nxp2, npt, ns
      common /ptprmc/ wp(is), qm(is), q(is), vpe(is), vpa(is),
     &                vd(is), pch(is), np(is)
      common /timecm/ itime,ntime,iecrct,iwrite,jobnum
      common /diagcm/ iediag, ifdiag, ikdiag, ipdiag, isdiag, ivdiag,
     &                ieplot, ifplot, ikplot, ipplot, isplot, ivplot
      common /otherc/ vmin,vmax
      common /resclc/ rex, ret, rev, ree, reb, rej, rer, res
      common /inputc/ dx, dt, cv, wc, angle
c
c-- size of arrays
      if(nx+2>ix) stop 'number of grids (nx) is too large'
      if(ns>is) stop 'number of species (ns) is too large'
      npa = 0
      do 10 i=1,ns
      npa = npa + np(i)
   10 continue
      if(npa>in) stop 'number of particles is too large'
c-- courant condition
      if(dx/dt<=cv) then
        print*,'courant condition is not satisfied'
        print*,'  make dt less than ', dx/cv
        stop
      end if
c-- paramters for diagnostics, etc.
      if(iediag==0) iediag = 99999999
      if(isdiag==0) isdiag = 99999999
      if(ifdiag==0) ifdiag = 99999999
      if(ipdiag==0) ipdiag = 99999999
      if(ikdiag==0) ikdiag = 99999999
      if(iecrct==0) iecrct = 99999999
      if(iwrite==0) iwrite = 99999999
      return
      end
