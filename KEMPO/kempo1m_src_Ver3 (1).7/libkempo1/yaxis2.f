C     ************************************************
      SUBROUTINE YAXIS2(X,Y,DY,DY1,MD1,DT,H,VMIN,VMAX,
     &           N,NTEXT,NT)
C     ************************************************
C     PROGRAMED BY Y. OMURA
C
      DIMENSION VM(50)
      CHARACTER NTEXT*(*)
      MD0=MD1
      IF(MD1.GE.1000) MD0=MD1-1000
      DY2=DY1*FLOAT(MD0)
      CALL PLOT(X,Y,3)
      CALL PLOT(X,Y+DY,2)
      DDY=DY1
      YY=Y
      SDT=0.6*DT
      X1=X
      IF(MD1.GE.1000) X1=X-SDT
      X2=X+SDT
      ND1M=RNDOFF(DY/DY1,3)+1.
      DO 10 I=1,ND1M
        CALL PLOT(X1,YY,3)
        CALL PLOT(X2,YY,2)
        YY=YY+DDY
   10 CONTINUE
      X1=X
      IF(MD1.GE.1000) X1=X-DT
      X2=X+DT
      DDY=DY2
      YY=Y
      ND3=RNDOFF(DY/DY2,3)+1.
      ND2=ND3-1
      DO 30 I=1,ND3
        CALL PLOT(X1,YY,3)
        CALL PLOT(X2,YY,2)
        YY=YY+DDY
   30 CONTINUE
      IF(H==0.0) RETURN
      ZMIN=ABS(VMIN)
      ZMAX=ABS(VMAX)
      IF(ZMIN>ZMAX) ZMAX=ZMIN
      CALL ETRANS(ZMAX,ZMAX,NPW)
      GFACT=10.**NPW
      NPN=N-1
      IF((NPW<1).OR.(NPW.GE.N)) GO TO 80
      GFACT=1.
      NPN=N-NPW-1
   80 CONTINUE
      IF(N<=0) GFACT=1.
      IF(N==1) NPN=-1
      IF(N<=-1) NPN=-N
      GMIN=VMIN/GFACT
      GMAX=VMAX/GFACT
      NCMAX=1
      CD=0.0
      HA=H
      IF(H.GE.0.0) GO TO 40
      HA=-H
   40 DSY=0.5*HA
      IF(N==0) GO TO 110
      IF(N==999) GO TO 555
  120 DG=(GMAX-GMIN)/10.**(-NPN)
      DG=ABS(DG)
      IF(DG>2.) GO TO 110
      NPN=NPN+1
      GO TO 120
  110 CONTINUE
      IF(NPN==0) NPN=-1
      DVM=(GMAX-GMIN)/(DY/DY2)
      VMM=GMIN
      DO 70 I=1,ND3
        VM(I)=VMM
        VMM=VMM+DVM
   70 CONTINUE
      DTH=DT*H
      IF(MD1.GE.1000) CD=ABS(DT)
      IF(DTH<0.0)  CD=ABS(DT)
      NCMAX=0
      DO 50 J=1,ND3
        Z=VM(J)
        Z=RNDOFF(Z,NPN)
        ZABS=ABS(Z)
        I=0
   15   I=I+1
        ZMAX=10.**I
        IF(ZABS.GE.ZMAX) GO TO 15
        IF(Z<0) I=I+1
        NC=I+NPN+1
      IF(NC>NCMAX) NCMAX=NC
        DSX=H*(FLOAT(NC)+0.7)+CD
        IF(H<0.0) DSX=-HA*0.7-CD
        XX=X-DSX
        YY=Y-DSY+DDY*(FLOAT(J)-1.)
        CALL NUMBER(XX,YY,HA,Z,0.,NPN)
   50 CONTINUE
C ===================================
  555 CONTINUE
C     CALL MSYMCT(NTEXT,NT,ANT,LINES)
      ANT=FLOAT(NT)
      LINES=1
C ===================================
      IF(H.GE.0.) DSX=H*(FLOAT(NCMAX)+1.5)+CD+(LINES-1)*1.2*HA*12./7.
      IF(H<0.0) DSX=-HA*(FLOAT(NCMAX)+2.7)-CD
      XX=X-DSX
      YY=Y+0.5*DY-0.6*HA*ANT
      IF(NT==0) GO TO 77
      CALL SYMBOL(XX,YY,1.2*HA,NTEXT,90.,NT)
   77 XX=X-3.*HA
      IF(N==999) RETURN
      IF(H<0.0) XX=X
      YY=Y+DY+HA*1.2
      CALL ENUMBR(XX,YY,1.2*HA,GFACT,0.,-2)
      RETURN
      END
