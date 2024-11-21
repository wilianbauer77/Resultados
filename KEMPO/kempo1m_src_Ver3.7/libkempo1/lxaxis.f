C     ****************************************************
      SUBROUTINE LXAXIS(XI,YI,XLI,DT,H,VMIN,VMAX,NTEXT,NT)
C     ****************************************************
C     PROGRAMED BY Y. OMURA
C
      CHARACTER NTEXT*(*)
C
      AMAX=ABS(VMAX)
      AMIN=ABS(VMIN)
      IF(AMIN>AMAX) THEN
        Y=YI
        X=XI+XLI
        XL=-XLI
        DUMMY=AMAX
        AMAX=AMIN
        AMIN=DUMMY
      ELSE
        X=XI
        Y=YI
        XL=XLI
      END IF
      CALL PLOT(X,Y,3)
      CALL PLOT(X+XL,Y,2)
      X1=ALOG10(AMIN)
      X2=ALOG10(AMAX)
      D=XL/(X2-X1)
      DTY=0.0
      DTH=DT*H
      HA=ABS(H)
      DF=ABS(D*0.3)
      IF(HA>DF) HA=DF
      IF(DTH<0.) DTY=ABS(DT)
      DSY=1.8*HA+DTY
      IF(H<0.) DSY=-0.5*HA-DTY
      VS=ASCALE(AMIN,1.0)
      CCC=(VS-AMIN)/AMIN
      IF((CCC>0.99).AND.(CCC<1.01)) VS=AMIN
      CALL ETRANS(VS,XT,NP)
      I=INT(RNDOFF(XT,1))
      YY=Y-DSY
      DDX=X2-X1
      MAXNP=1
      IF(DDX<1.5) GO TO 50
   10 IF(I.GE.10) THEN
        I=1
       NP=NP+1
      END IF
      V=FLOAT(I)*(10.**NP)
      IF(V>AMAX) GO TO 20
      DD=DT*0.5
      XX=X+(ALOG10(V)-X1)*D
      IF(I==1) THEN
      ISFT=1
      CALL ETRANS(V,XT,MP)
      IF(MP<0) ISFT=ISFT+1
      IAMP=IABS(MP)
      IF(IAMP.GE.10) ISFT=ISFT+1
      IF(ISFT>MAXNP) MAXNP=ISFT
      XXP=XX-0.5*HA*(FLOAT(ISFT)*0.8+2.)
         DD=DT
        CALL ENUMBR(XXP,YY,HA,V,0.,-3)
      END IF
      CALL PLOT(XX,Y,3)
      CALL PLOT(XX,Y+DD,2)
      I=I+1
      GO TO 10
   50 CONTINUE
      DDS=0.05*D*0.66666
      HAH=HA*1.2
      IMAX=9
      IF(DDS<HAH) THEN
        IMAX=8
      HAS=DDS
      ELSE
        HAS=HA*0.7
      END IF
      YYS=Y-0.5*HA-DTY-HAS
      IF(H<0.) YYS=Y-DSY
   30 IF(I.GE.10) THEN
        I=1
       NP=NP+1
      END IF
      V=FLOAT(I)*(10.**NP)
      IF(V>AMAX) GO TO 20
      DD=DT*0.5
      XX=X+(ALOG10(V)-X1)*D
      IF(I==1) THEN
         DD=DT
      ISFT=1
      CALL ETRANS(V,XT,MP)
      IF(MP<0) ISFT=ISFT+1
      IAMP=IABS(MP)
      IF(IAMP.GE.10) ISFT=ISFT+1
      XXP=XX-0.5*HA*(FLOAT(ISFT)*0.8+2.)
        CALL ENUMBR(XXP,YY,HA,V,0.,-3)
      ELSE IF(I<=IMAX) THEN
      CALL NUMBER(XX-0.5*HAS,YYS,HAS,FLOAT(I),0.,-1)
      END IF
      CALL PLOT(XX,Y,3)
      CALL PLOT(XX,Y+DD,2)
      I=I+1
      GO TO 30
   20 CONTINUE
C      CALL MSYMCT(NTEXT,NT,ANT,LINES)
      ANT=FLOAT(NT)
      LINES=1
C
      YY=Y-3.7*HA-DTY
      IF(H<0.) YY=Y+2.5*HA+DTY+(LINES-1)*1.2*HA*12./7.
      XX=X+XL*0.5-0.6*HA*ANT
      IF(NT==0) RETURN
      CALL SYMBOL(XX,YY,HA*1.2,NTEXT,0.,NT)
      RETURN
      END
