C     *************************
      REAL FUNCTION RNDOFF(X,N)
C     *************************
C     PROGRAMED BY Y. OMURA
C
      IS=1
      IF(X<0.) IS=-1
      XA=ABS(X)
      NA=N
      IF(NA<0) NA=N+1
      FACT=10.**NA
      XA=XA*FACT
      IXA=INT(XA)
      DXA=XA-FLOAT(IXA)
      IF(DXA.GE.0.5) IXA=IXA+1
      ADD=0.001
      IF(IXA==0) ADD=0.
      RNDOFF=(FLOAT(IXA)+ADD)/FACT*FLOAT(IS)
      RETURN
      END
