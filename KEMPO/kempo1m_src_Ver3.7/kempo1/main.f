c****************************************************************
c
c      1D Electromagnetic Full Particle Code  : KEMPO1
c
c           by Yoshiharu Omura
c           Radio Atmospheric Science Center, Kyoto University
c           Uji, Kyoto, 611, Japan   
c               E-mail: omura@kurasc.kyoto-u.ac.jp
c               FAX:  +81-774-31-8463
c
c                   Version 4.4    November 27, 1992
c
c****************************************************************
      program main
      common /timecm/ itime,ntime,iecrct,iwrite,jobnum
      common /diagcm/ iediag, ifdiag, ikdiag, ipdiag, isdiag, ivdiag,
     &                ieplot, ifplot, ikplot, ipplot, isplot, ivplot
c
      call plots
      call factor(0.9)
      call input
c     call chkprm
      call pltprm
      call renorm
      if(jobnum<=1) then
        itime = 0
        call inital
        call positn
        call chargv
        call ecrrct
      else
        call reader
      endif
      call fldplt
      call phsplt
      call vdsplt
      ist = itime
      print*, 'main loop'
c
      do 100 j = ist+1, ist+ntime
        itime = j
c       call bfield
        call velcty
c       call bvlcty
        call positn
c       call currnt
c       call curntv
c       call bfield
c       call efield
        call positn
c       call rcycle
c       if( mod(j,iecrct)==0) then
            call charge
            call ecrrct
c       endif
        if( mod(j,ifdiag)==0 ) then
c           if( mod(ifdiag,iecrct)/=0 ) call charge
            call fldplt
        endif
        if( mod(j,ikdiag)==0 ) call kspplt
        if( mod(j,ipdiag)==0 ) call phsplt
        if( mod(j,ivdiag)==0 ) call vdsplt
        if( mod(j,isdiag)==0 ) call spectr
        if( mod(j,iediag)==0 ) call energy
        if( mod(j,iwrite)==0 ) call writer
 100  continue
      call writer
      call plot(0.,0.,999)
      stop
      end
