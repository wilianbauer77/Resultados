C     **************************
      SUBROUTINE QLOOK2(BR,N,IP)
C     **************************
C       PROGRAMED BY Y. OMURA
C
      DIMENSION BR(N)
      COMMON /QLKCM1/IC,VMIN,VMAX,RMIN,RMAX,X,Y,XAL,YAL
C
      DD=XAL/40.
      IF(IP<=-1) DD=XAL/10.
      BD=XAL/80.
      PD=XAL/100.
      IF(IC==0) THEN
      YFACT=YAL/(VMAX-VMIN)
      DDX=XAL/FLOAT(N-1)
      YMAX=Y+YAL
      XX=X
      YY=Y+(BR(1)-VMIN)*YFACT
      IF(YY>YMAX) YY=YMAX
      IF(YY<Y)  YY=Y
      CALL DPLOT6(XX,YY,IP,DD,BD,PD)
      IPEN=3
      DO 10 I=2,N
      XX=XX+DDX
      IF(BR(I)==999.) GO TO 10
      YY=Y+(BR(I)-VMIN)*YFACT
      IF(BR(I-1)==999.) THEN
        IP2=3
      ELSE
        IP2=2
      END IF
      IF(YY>YMAX) THEN
        YY=YMAX
        IF(IPEN==2) THEN
          CALL DPLOT(XX,YY,IP2)
          IPEN=3
        ELSE
          CALL DPLOT(XX,YY,3)
        ENDIF
      ELSE IF(YY<Y) THEN
        YY=Y
        IF(IPEN==2) THEN
          CALL DPLOT(XX,YY,IP2)
          IPEN=3
        ELSE
          CALL DPLOT(XX,YY,3)
        ENDIF
      ELSE
        CALL DPLOT(XX,YY,IP2)
        IPEN=2
      ENDIF
  10  CONTINUE
      ELSE
      ALR=ALOG10(RMIN)
      YFACT=YAL/(ALOG10(RMAX)-ALR)
      DDX=XAL/FLOAT(N-1)
      YMAX=Y+YAL
      XX=X
      YY=Y+(ALOG10(BR(1))-ALR)*YFACT
      IF(YY>YMAX) YY=YMAX
      IF(YY<Y)  YY=Y
      CALL DPLOT6(XX,YY,IP,DD,BD,PD)
      IPEN=3
      DO 20 I=2,N
      IF(BR(I)==999.) GO TO 20
      XX=XX+DDX
      YY=Y+(ALOG10(BR(I))-ALR)*YFACT
      IF(BR(I-1)==999.) THEN
        IP2=3
      ELSE
        IP2=2
      END IF
      IF(YY>YMAX) THEN
        YY=YMAX
        IF(IPEN==2) THEN
          CALL DPLOT(XX,YY,IP2)
          IPEN=3
        ELSE
          CALL DPLOT(XX,YY,3)
        ENDIF
      ELSE IF(YY<Y) THEN
        YY=Y
        IF(IPEN==2) THEN
          CALL DPLOT(XX,YY,IP2)
          IPEN=3
        ELSE
          CALL DPLOT(XX,YY,3)
        ENDIF
      ELSE
        CALL DPLOT(XX,YY,IP2)
        IPEN=2
      ENDIF
  20  CONTINUE
      ENDIF
      RETURN
      END
