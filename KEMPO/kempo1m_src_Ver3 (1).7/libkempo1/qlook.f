C     ********************************************
      SUBROUTINE QLOOK(AR,N,X,Y,XAL,YAL,XMIN,XMAX,
     &                   NXTEXT,NXI,NYTEXT,NY)
C     ********************************************
C     PROGRAMED BY Y. OMURA, RASC/KYOTO UNIV.
C
      DIMENSION AR(N)
      CHARACTER NXTEXT*(*), NYTEXT*(*)
      COMMON /QLKCM1/IC,VMIN,VMAX,CRMIN,CRMAX,CX,CY,CXAL,CYAL
      DATA NFX,NFY,IPPEN,ITPEN /3,3,3,4/
      DATA DTMOD, HMOD /0.,0./
C
      NX=NXI
      IPLOT=0
      JPLOT=0
      KPLOT=0
      CALL NEWPEN(ITPEN)
      IF(IABS(NX).GE.20000) THEN
         KPLOT=1
         NX=NX/IABS(NX)*(IABS(NX)-20000)
      ENDIF
      IF(IABS(NX).GE.10000) THEN
        JPLOT=1
        NX=NX/IABS(NX)*(IABS(NX)-10000)
      ENDIF
      IF(IABS(NX).GE.1000.AND.IABS(NX)<10000) THEN
        IPLOT=1
        NX=NX/IABS(NX)*(IABS(NX)-1000)
      ENDIF
      CALL MAXMIN(AR,N,RMIN,RMAX)
      IF(RMIN/=0.) THEN
        A=(RMAX-RMIN)/RMIN
        IF(ABS(A)<1.E-5) THEN
      H=MAX(XAL,YAL)/22.
      CALL SYMBOL(X+XAL*0.2,Y+YAL*0.5,H,'NO VARIATION',0.,12)
      CALL PRMPLT(X+XAL*0.2,Y+YAL*0.4,H*0.8,0.,'VALUE',5,RMIN,3)
C       PRINT *,'+QLOOK+  NO VARIATION IN DATA : Y=',RMIN
      IPLOT=1
      ENDIF
      ENDIF
      CRMIN=RMIN
      VMIN=RMIN
      VMAX=RMAX
      CRMAX=RMAX
      CX=X
      CY=Y
      CXAL=XAL
      CYAL=YAL
      IC=0
      IF(JPLOT==1) GOTO 30
      IF(RMIN>0..AND.NX.GE.0) THEN
        A=(RMAX-RMIN)/RMIN
        IF(A>10..AND.NY.GE.0) THEN
          A=0.
          DO 10 I=1,N
   10     A=A+AR(I)
          A=A/FLOAT(N)
          DL=ABS(A-(RMIN+RMAX)*0.5)
          A=0.
          DO 20 I=1,N
   20     A=A+ALOG10(AR(I))
          A=A/FLOAT(N)
          DG=ABS(A-(ALOG10(RMIN)+ALOG10(RMAX))*0.5)
          IF(DG<DL) IC=1
        ELSE
          IF(NY<0) IC=1
        ENDIF
      END IF
      NXA=IABS(NX)
      NYA=IABS(NY)
      IF(IC==1) GO TO 30
      CALL ETRANS(RMIN,A,NP1)
      CALL ETRANS(RMAX,A,NP2)
      IF((RMIN<0.).AND.(RMAX>0.)) THEN
        AMIN=ABS(RMIN)
        AMAX=ABS(RMAX)
        IF(AMAX.GE.AMIN) THEN
        VMAX=ASCALE(RMAX,2.)
        A=2.*10.**(NP2-NP1)
        VMIN=ASCALE(RMIN,A)
        ELSE
        VMIN=ASCALE(RMIN,2.)
        A=2.*10.**(NP1-NP2)
        VMAX=ASCALE(RMAX,A)
        END IF
      ELSE IF(RMIN==0.) THEN
      VMIN=RMIN
      VMAX=ASCALE(RMAX,2.)
      ELSE IF(RMAX==0.) THEN
      VMIN=RMAX
      VMAX=ASCALE(RMIN,2.)
      ELSE IF(RMIN>0.) THEN
        A=(RMAX-RMIN)/RMAX
        A=ASCALE(A,2.)*0.5
        VMAX=ASCALE(RMAX,A)
        A=A*10.**(NP2-NP1)
        VMIN=ASCALE(RMIN,-A)
      ELSE
        A=(RMIN-RMAX)/RMIN
      A=ASCALE(A,2.)*0.5
      VMAX=ASCALE(RMIN,A)
        A=A*10.**(NP1-NP2)
        VMIN=ASCALE(RMAX,-A)
      END IF
      IF(KPLOT==1) THEN
        VMIN=0.0
        IF(ABS(RMAX).GE.ABS(RMIN)) THEN
          VMAX=ASCALE(RMAX,2.)
         ELSE
          VMAX=ASCALE(RMIN,2.)
         ENDIF
        ENDIF
   30 CONTINUE
      NXA=IABS(NX)
      NYA=IABS(NY)
      INX=INT(XAL)
      INX=INT(FLOAT(INX)/2.)*2
      INY=INT(YAL)
      INY=INT(FLOAT(INY)/2.)*2
      IF(DTMOD==0.) THEN
        DT=MAX(XAL,YAL)/30.
      ELSE
        DT = DTMOD
      ENDIF
      IF(HMOD==0.) THEN
        H=MAX(XAL,YAL)/22.
      ELSE
        H = HMOD
      ENDIF
      IF(XMIN==0..AND.XMAX>0.) THEN
      IAXIS=1
      CALL ETRANS(XMAX,A,NP)
      A=RNDOFF(A,3)
      IF(A.GE.8.) THEN
        DAX=1.
        LDX=4
      ELSE IF(A.GE.6.) THEN
        DAX=0.5
        LDX=6
      ELSE IF(A.GE.5.) THEN
        DAX=0.5
        LDX=5
      ELSE IF(A.GE.4.) THEN
        DAX=0.5
        LDX=4
      ELSE IF(A.GE.3.) THEN
        DAX=0.5
        LDX=3
      ELSE IF(A.GE.2.) THEN
        DAX=0.2
        LDX=5
      ELSE
        DAX=0.1
        LDX=5
      ENDIF
        DAX=DAX*10.**NP/XMAX*XAL
        CALL XAXIS2(X,Y,XAL,DAX,LDX,DT,H,XMIN,XMAX,NFX,NXTEXT,NXA)
      ELSE
      CALL XAXIS1(X,Y,XAL,INX,2,DT,H,XMIN,XMAX,NFX,NXTEXT,NXA)
      IAXIS=0
      ENDIF
      IF(IC==0) THEN
      IF(IAXIS==1) THEN
        CALL XAXIS2(X,Y+YAL,XAL,DAX,LDX,-DT,0.,0.,0.,0,0,0)
      ELSE
      CALL XAXIS1(X,Y+YAL,XAL,INX,2,-DT,0.,0.,0.,0,0,0)
      ENDIF
      IF(VMIN==0..AND.VMAX>0.) THEN
      CALL ETRANS(VMAX,B,NQ)
      B=RNDOFF(B,3)
      IF(B.GE.8.) THEN
        DBY=1.
        LDY=4
      ELSE IF(B.GE.6.) THEN
        DBY=0.5
        LDY=6
      ELSE IF(B.GE.5.) THEN
        DBY=0.5
        LDY=5
      ELSE IF(B.GE.4.) THEN
        DBY=0.5
        LDY=4
      ELSE IF(B.GE.3.) THEN
        DBY=0.5
        LDY=3
      ELSE IF(B.GE.2.) THEN
        DBY=0.2
        LDY=5
      ELSE
        DBY=0.1
        LDY=5
      ENDIF
        DBY=DBY*10.**NQ/(VMAX-VMIN)*YAL
      CALL YAXIS2(X,Y,YAL,DBY,LDY,DT,H,VMIN,VMAX,NFY,NYTEXT,NYA)
      CALL YAXIS2(X+XAL,Y,YAL,DBY,LDY,-DT,0.,0.,0.,0,0,0)
      ELSE IF(ABS(VMIN)==ABS(VMAX)) THEN
      CALL ETRANS(VMAX,B,NQ)
      B=RNDOFF(B,3)
      IF(B.GE.8.) THEN
        DBY=2.
        LDY=4
      ELSE IF(B.GE.6.) THEN
        DBY=1.
        LDY=6
      ELSE IF(B.GE.5.) THEN
        DBY=1.
        LDY=5
      ELSE IF(B.GE.4.) THEN
        DBY=1.0
        LDY=4
      ELSE IF(B.GE.3.) THEN
        DBY=1.0
        LDY=3
      ELSE IF(B.GE.2.) THEN
        DBY=0.5
        LDY=4
      ELSE
        DBY=0.2
        LDY=5
      ENDIF
        DBY=DBY*10.**NQ/(VMAX-VMIN)*YAL
      CALL YAXIS2(X,Y,YAL,DBY,LDY,DT,H,VMIN,VMAX,NFY,NYTEXT,NYA)
      CALL YAXIS2(X+XAL,Y,YAL,DBY,LDY,-DT,0.,0.,0.,0,0,0)
      ELSE
      CALL YAXIS1(X,Y,YAL,INY,2,DT,H,VMIN,VMAX,NFY,NYTEXT,NYA)
      CALL YAXIS1(X+XAL,Y,YAL,INY,2,-DT,0.,0.,0.,0,0,0)
      ENDIF
      ELSE
      CALL PLOT(X+XAL,Y,3)
      CALL PLOT(X+XAL,Y+YAL,2)
      CALL PLOT(X,Y+YAL,2)
      RMIN0=RMAX*1.E-10
      IF(RMIN<RMIN0) THEN
        RMIN=RMIN0
        CRMIN=RMIN
      ELSE
        CALL ETRANS(RMAX,RR1,NN1)
        CALL ETRANS(RMIN,RR2,NN2)
        IF(NN1==NN2) THEN
          RMIN=9.5*(10.**(NN2-1))
          CRMIN=RMIN
          IF(RR1<2.1) THEN
            RMAX=2.1*(10.**NN1)
            CRMAX=RMAX
          END IF
        END IF
      END IF
      CALL LYAXIS(X,Y,YAL,-DT,H,RMIN,RMAX,NYTEXT,NYA)
      END IF
      CALL NEWPEN(IPPEN)
      IF(IPLOT/=1) CALL QLOOK2(AR,N,1)
      RETURN
C     ********************
      ENTRY QLKMOD(NEX,NEY,JPPEN,JTPEN)
C     *********************
C
      NFX=NEX
      NFY=NEY
      IPPEN=JPPEN
      ITPEN=JTPEN
C
      RETURN
C
C     ********************
      ENTRY QLKMD2(DTIN,HIN)
C     *********************
C
      DTMOD = DTIN
      HMOD = HIN
C
      RETURN
      END
