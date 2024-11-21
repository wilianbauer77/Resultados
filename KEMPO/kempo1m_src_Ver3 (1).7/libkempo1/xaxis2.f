C     ************************************************
      SUBROUTINE XAXIS2(X,Y,DX,DX1,MD1,DT,H,VMIN,VMAX,
     &                    N,NTEXT,NT)
C     ************************************************
C     PROGRAMED BY Y. OMURA
C
      DIMENSION VM(50)
      CHARACTER*1 NTEXT(NT)
      MD0=MD1
      IF(MD1.GE.1000) MD0=MD1-1000
      DX2=DX1*FLOAT(MD0)
      CALL PLOT(X,Y,3)
      CALL PLOT(X+DX,Y,2)
      DDX=DX1
      XX=X
      SDT=0.6*DT
      Y1=Y
      IF(MD1.GE.1000) Y1=Y-SDT
      Y2=Y+SDT
      ND1M=RNDOFF(DX/DX1,3)+1.
      DO 10 I=1,ND1M
        CALL PLOT(XX,Y1,3)
        CALL PLOT(XX,Y2,2)
        XX=XX+DDX
   10 CONTINUE
      Y1=Y
      IF(MD1.GE.1000) Y1=Y-DT
      Y2=Y+DT
      DDX=DX2
      XX=X
      ND3=RNDOFF(DX/DX2,3) + 1.
      ND2=ND3-1
      DO 30 I=1,ND3
        CALL PLOT(XX,Y1,3)
        CALL PLOT(XX,Y2,2)
        XX=XX+DDX
   30 CONTINUE
      IF(H==0.0) RETURN
      NPN=N-1
      IF(N<=1) NPN=-1
      ZMAX=ABS(VMAX)
      ZMIN=ABS(VMIN)
      IF(ZMIN>ZMAX) ZMAX=ZMIN
      CALL ETRANS(ZMAX,ZMAX,NPW)
      GFACT=10.**NPW
      IF((NPW<1).OR.(NPW.GE.N)) GO TO 80
      GFACT=1.
      NPN=N-NPW-1
   80 CONTINUE
      IF(N<=-1) NPN=-N
      IF(N<=0) GFACT=1.
      GMAX=VMAX/GFACT
      GMIN=VMIN/GFACT
      CD=0.0
      DTH=DT*H
      IF(MD1.GE.1000) CD=ABS(DT)
      IF(DTH<0.0)  CD=ABS(DT)
      HA=H
      DSY=1.9*H+CD
      IF(H.GE.0.0) GO TO 40
      HA=-H
      DSY=-0.9*HA-CD
   40 CONTINUE
      IF(N==999) GO TO 555
      IF(N==0) GO TO 110
  120 DG=(GMAX-GMIN)/10.**(-NPN)
      DG=ABS(DG)
      IF(DG>2.) GO TO 110
      NPN=NPN+1
      GO TO 120
  110 CONTINUE
      IF(NPN==0) NPN=-1
      DVM=(GMAX-GMIN)/(DX/DX2)
      VMM=GMIN
      DO 70 I=1,ND3
        VM(I)=VMM
        VMM=VMM+DVM
  70  CONTINUE
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
        DSX=HA*FLOAT(NC)*0.5
        XX=X-DSX+DDX*(FLOAT(J)-1.)
        YY=Y-DSY
        CALL NUMBER(XX,YY,HA,Z,0.,NPN)
   50 CONTINUE
C ===================================
C 555 CALL MSYMCT(NTEXT,NT,ANT,LINES)
  555 ANT=FLOAT(NT)
      LINES=1
C ===================================
      DSY=3.9*H+CD
      IF(H<0.0) DSY=-2.8*HA-CD-(LINES-1)*1.2*HA*12./7.
      XX=X+0.5*DX-0.6*HA*ANT
      YY=Y-DSY
      IF(NT==0) GO TO 77
      CALL SYMBOL(XX,YY,1.2*HA,NTEXT,0.,NT)
      IF(N==999) RETURN
   77 XX=X+DX-3.6*HA
      XT=X+0.5*DX+0.6*HA*FLOAT(NT)
      IF(XT.GE.XX) XX=XT+HA
      CALL ENUMBR(XX,YY,1.2*HA,GFACT,0.,-2)
      RETURN
      END
