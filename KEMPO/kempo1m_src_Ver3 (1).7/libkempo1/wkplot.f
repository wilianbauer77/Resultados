C     ************************************************************
      SUBROUTINE WKPLOT(AR,N1,M1,N,M,WK1,WK2,X0,Y0,XL,YL,SL,ST,NQ,
     #                  MQ,ID)
C     ************************************************************
C
C         PROGRAMED BY Y. OMURA 
C
      CHARACTER XTITLE*43,YTITLE*43
      DIMENSION AR(N1,M1),WK1(N1),WK2(M1)
C
        KAKUDO=0
      PI=3.1415926
      IC=ID
      IF(WK1(1)==999.) THEN
        IPMAX=1
        EMAX=WK1(2)
      ELSE
        IPMAX=0
      ENDIF
C
      NP=NQ
      MP=MQ
      N2=N/2
      M2=M/2
      IF((NP<1).OR.(NP>N2)) NP=N2
      IF((MP<1).OR.(MP>M2)) MP=M2
      MP1=MP+1
      WMAX=2.*PI*FLOAT(NP)/SL
      YMAX=2.*PI*FLOAT(MP)/ST
      M21=M2+1
      DX=XL/FLOAT(NP)
      DY=YL/FLOAT(MP)
C
C--- PLOT THE AXIS ----
      WMAX=RNDOFF(WMAX,4)
      YMAX=RNDOFF(YMAX,4)
      SIGN=1.
         LOGPLT=0
C
      NWMAX=3
      NYMAX=3
      IF(ABS(IC).GE.10000) THEN
        IC=IC/ABS(IC)*(ABS(IC)-10000)
        DANG=WK1(3)
        KAKUDO=1
      END IF
      IF(ABS(IC).GE.1000) THEN
        IC=IC/ABS(IC)*(ABS(IC)-1000)
        LOGPLT=1
      ENDIF
      IF(IC<0) SIGN=-1.
      DYD=YL/YMAX*ASCALE(YMAX,-0.4)/8.
      IF((YMAX<10.).AND.(YMAX.GE.2.)) THEN
        DYD=YL/YMAX*ASCALE(YMAX,-2.)/8.
      ENDIF
      IF((YMAX<8.).AND.(YMAX.GE.6.)) THEN
        DYD=YL/YMAX*ASCALE(YMAX,-2.)/6.
      ENDIF
      IF((YMAX<2.).AND.(YMAX.GE.1.0)) THEN
        DYD=YL/YMAX*ASCALE(YMAX,-1.)/4.
      ENDIF
      IF((YMAX<1.).AND.(YMAX.GE.0.5)) THEN
        DYD=YL/YMAX*ASCALE(YMAX,-5.)/2.
      ENDIF
C
C      IF(WMAX<0.1.OR.YMAX>10.0 ) NWMAX=1001
C      IF(YMAX<1.0.OR.YMAX>10.0 ) NYMAX=1001
      IF(DYD==0.) DYD=YL/(YMAX*10.)*INT(YMAX*10.)/8.
        XTITLE='K'
        NKT=1
      YTITLE='OMEGA'
      NWT=5
      CALL XAXIS1(X0,Y0,XL,4,2,XL*0.02,XL*0.05,
     &            0.,WMAX*SIGN,NWMAX,XTITLE,NKT)
      CALL XAXIS1(X0,Y0+YL,XL,4,2,-XL*0.02,0.,
     &            0.,WMAX*SIGN,NWMAX,XTITLE,NKT)
      CALL YAXIS2(X0,Y0,YL,DYD,2,XL*0.02,XL*0.05,
     &            0.,YMAX,NYMAX,YTITLE,NWT)
      CALL YAXIS2(X0+XL,Y0,YL,DYD,2,-XL*0.02,0.,
     &            0.,YMAX,NYMAX,YTITLE,NWT)
C
      IF(KAKUDO==1) THEN
        CALL PRMPLT(X0+XL*0.5,Y0+YL*1.15,XL*0.04,0.,'ANGLE',
     $               5,DANG,-2)
      END IF
C--- GET THE MAXIMUM COMPONENT OF AR ----------
      FMAX=1.E-7
      FMIN=1.E+9
      DO 15 L=1,N
   15 WK1(L)=0.
      I0=3
      IF(IC==-1) I0=4
      J=1
      DO 14 K=1,MP
      WK1(1)=AR(1,J)
      I=I0
      DO 12 L=2,NP
      WK1(L)=AR(I,J)
      IF(IC==0) WK1(L)=WK1(L)+AR(I+1,J)
      I=I+2
   12 CONTINUE
      IF(NP==N2) WK1(N2+1)=AR(2,J)
      WK1(N2+1)=AR(2,J)
      CALL MAXMIN(WK1,N,VMIN,VMAX)
      IF(VMIN<FMIN) FMIN=VMIN
      IF(VMAX>FMAX) FMAX=VMAX
      J=J+2
   14 CONTINUE
      I=I0
      DO 99 L=2,NP
      WK1(L)=AR(I,2)
      IF(IC==0) WK1(L)=WK1(L)+AR(I+1,2)
      I=I+2
   99 CONTINUE
      IF(NP==N2) WK1(N2+1)=AR(2,2)
      WK1(N2+1)=AR(2,2)
      CALL MAXMIN(WK1,N,VMIN,VMAX)
      IF(VMIN<FMIN) FMIN=VMIN
      IF(VMAX>FMAX) FMAX=VMAX
      IF(IPMAX==1) FMAX=EMAX
      CALL PRMPLT(X0+XL*0.5,Y0+YL*1.05,XL*0.04,0.,'MAX',3,FMAX,3)
      FACT=DX/FMAX
      IF(LOGPLT/=0) THEN
      IF(FMIN<=0) FMIN=FMAX*1.E-2
      AFMAX=ALOG10(FMAX)
      AFMIN=ALOG10(FMIN)
      FACT=DX/(AFMAX-AFMIN)
      ENDIF
      J=1
      DO 22 K=1,M2
        WK2(K)=AR(1,J)
        J=J+2
   22 CONTINUE
      WK2(M21)=AR(1,2)
      PX=X0
      PY=Y0
      DD=WK2(1)*FACT
      IF(LOGPLT/=0) DD=(ALOG10(WK2(1))-AFMIN)*FACT
      IF(DD>DX) DD=DX
       IF(DD<0.) DD=0.
      CALL PLOT(PX,PY,3)
      CALL PLOT(PX,PY+YL,2)
      CALL PLOT(PX+DD,PY,3)
      DO 32 K=2,MP1
        PY=PY+DY
        DD=WK2(K)*FACT
      IF(LOGPLT/=0) DD=(ALOG10(WK2(K))-AFMIN)*FACT
      IF(DD>DX) DD=DX
       IF(DD<0.) DD=0.
        CALL PLOT(PX+DD,PY,2)
   32 CONTINUE
      PY=Y0
      DO 42 K=1,MP1
        DD=WK2(K)*FACT
      IF(LOGPLT/=0) DD=(ALOG10(WK2(K))-AFMIN)*FACT
      IF(DD>DX) DD=DX
       IF(DD<0.) DD=0.
        CALL PLOT(PX,PY,3)
        CALL PLOT(PX+DD,PY,2)
        PY=PY+DY
   42 CONTINUE
      I=I0
      X=X0+DX
      DO 10 L=2,NP
      J=1
      DO 20 K=1,M2
        WK2(K)=AR(I,J)
      IF(IC==0) WK2(K)=WK2(K)+AR(I+1,J)
        J=J+2
   20 CONTINUE
      WK2(M21)=AR(I,2)
      PX=X
      PY=Y0
      DD=WK2(1)*FACT
      IF(LOGPLT/=0) DD=(ALOG10(WK2(1))-AFMIN)*FACT
      IF(DD>DX) DD=DX
       IF(DD<0.) DD=0.
      CALL PLOT(PX,PY,3)
      CALL PLOT(PX,PY+YL,2)
      CALL PLOT(PX+DD,PY,3)
      DO 30 K=2,MP1
        PY=PY+DY
        DD=WK2(K)*FACT
      IF(LOGPLT/=0) DD=(ALOG10(WK2(K))-AFMIN)*FACT
      IF(DD>DX) DD=DX
       IF(DD<0.) DD=0.
        CALL PLOT(PX+DD,PY,2)
   30 CONTINUE
      PY=Y0
      DO 40 K=1,MP1
        DD=WK2(K)*FACT
      IF(LOGPLT/=0) DD=(ALOG10(WK2(K))-AFMIN)*FACT
      IF(DD>DX) DD=DX
       IF(DD<0.) DD=0.
        CALL PLOT(PX,PY,3)
        CALL PLOT(PX+DD,PY,2)
        PY=PY+DY
   40 CONTINUE
      X=X+DX
      I=I+2
   10 CONTINUE
      IF(NP==N2) I=2
      J=1
      DO 21 K=1,M2
        WK2(K)=AR(I,J)
        J=J+2
   21 CONTINUE
      WK2(M21)=AR(I,2)
      PX=X
      PY=Y0
      DD=WK2(1)*FACT
      IF(LOGPLT/=0) DD=(ALOG10(WK2(1))-AFMIN)*FACT
      IF(DD>DX) DD=DX
       IF(DD<0.) DD=0.
      CALL PLOT(PX,PY,3)
      CALL PLOT(PX,PY+YL,2)
      CALL PLOT(PX+DD,PY,3)
      DO 31 K=2,MP1
        PY=PY+DY
        DD=WK2(K)*FACT
      IF(LOGPLT/=0) DD=(ALOG10(WK2(K))-AFMIN)*FACT
      IF(DD>DX) DD=DX
       IF(DD<0.) DD=0.
        CALL PLOT(PX+DD,PY,2)
   31 CONTINUE
      PY=Y0
      DO 41 K=1,MP1
        DD=WK2(K)*FACT
      IF(LOGPLT/=0) DD=(ALOG10(WK2(K))-AFMIN)*FACT
      IF(DD>DX) DD=DX
       IF(DD<0.) DD=0.
        CALL PLOT(PX,PY,3)
        CALL PLOT(PX+DD,PY,2)
        PY=PY+DY
   41 CONTINUE
      WK1(2)=FMAX
      IF(IPMAX==1) WK1(1)=999.
      WK1(3)=FMIN
      WK1(4)=YMAX
      RETURN
      END
