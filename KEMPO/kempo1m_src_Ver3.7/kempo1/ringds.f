c    *********************************
      subroutine ringds(ar,na,vt,bbb)
c    *********************************
c
      dimension ar(na)
      ppii=3.14159265358979d0
      sq2=1.41421356
      dv=vt
      ds=vt**2/2./float(na)
      a=1./vt**2
      cc=1./(1.-bbb)
      v1=vt
      v2=v1
      j=1
   60 continue
      avm=-a*v1*v1
      fv1=cc*v1*(exp(avm)-exp(avm/bbb))
      dvd=ds/fv1
      v1=v1-dvd
      if(v1<=0.0) go to 70
      ar(j)=v1+0.5*dvd
      j=j+1
      go to 60
   70 continue
      avm=-a*v2*v2
      fv2=cc*v2*(exp(avm)-exp(avm/bbb))
      dvd=ds/fv2
      v2=v2+dvd
      ar(j)=v2-0.5*dvd
      j=j+1
      if(j<=na) go to 70
      call rexch(na,ar)
      return
      end
