      subroutine vinit(ar,na,v0,vt,xlos)
      implicit real*8 (a-h,o-z)
c    ***************************
c    *                         *
c    *    subroutine vinit     *
c    *                         *
c    ***************************
c
c   ..................................................
c   .                                                .
c   .  this subroutine gives a maxwell distribution  .
c   .  (xlos==0), or a loss cone distriburion      .
c   .  (xlos/=0). xlos is the loss cone factor.    .
c   .                                                .
c   ..................................................
c
      real*4 ar(na),v0,vt,xlos

      ppii = 3.14159265358979
      sq2  = 1.414213562

      if ( xlos/=0. ) goto 20

      vd  = sq2*vt
      nph = na/2
      ds  = sqrt(ppii)/float(na)

      n = 0
      x = 0.0
      do 10 i=1,nph
         dx = ds/exp(-x**2)
         x = x+dx
         dvp = vd*(x-0.5*dx)
         n = n+1
         ar(n) = v0+dvp
         n = n+1
         ar(n) = v0-dvp
   10 continue
      call rexch(na,ar)
      return

   20 continue
      axlos = abs(xlos)

      if ( xlos>0. ) then
         vt = v0*sq2/sqrt(axlos)
      else
         v0 = vt*sqrt(axlos)/sq2
      end if

      dv   =  vt
      a    =  1.0/(dv*dv)
      xx   =  0.5*(axlos+1.0)
      xxm  = -xx
      ts   =  0.5*(a**xxm)*dgamma(xx)
      ds   =  ts/float(na)
      sqvc =  axlos*0.5/a
      v1 = sqrt(sqvc)
      v2 = v1
      j  = 1

   60 continue
         avm = -a*v1*v1
         fv1 = v1**axlos*exp(avm)
         dvd = ds/fv1
          v1 = v1-dvd
         if ( v1<=0.0 ) goto 70
         ar(j) = v1+0.5*dvd
         j = j+1
      goto 60

   70 continue
          avm = -a*v2*v2
          fv2 = v2**axlos*exp(avm)
          dvd = ds/fv2
           v2 = v2+dvd
        ar(j) = v2-0.5*dvd
            j = j+1
      if ( j<=na ) goto 70

      call rexch(na,ar)
      return
      end

      real function dgamma(xx)
      implicit real*8(a-h,o-z)

      dimension cof(6)

      data cof,stp /76.18009173,-86.50532033,24.01409822,
     &     -1.231739516,0.120858003e-2,-0.536382e-5,2.50662827465/
      data half,one,fpf/0.5,1.0,5.5/

        x = xx-one
      tmp = x+fpf
      tmp = (x+half)*log(tmp)-tmp
      ser = one
      do 10 j=1,6
           x = x+one
         ser = ser+cof(j)/x
   10 continue
      dgamma = exp(tmp)*stp*ser

      return
      end
