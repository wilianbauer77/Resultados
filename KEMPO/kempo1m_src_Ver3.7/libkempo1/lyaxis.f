C     ****************************************************
      SUBROUTINE LYAXIS(XI,YI,YLI,DT,H,VMIN,VMAX,NTEXT,NT)
C     ****************************************************
C
      CHARACTER NTEXT*(*)
C
      AMAX=ABS(VMAX)
      AMIN=ABS(VMIN)
      IF(AMIN>AMAX) THEN
        X=XI
        Y=YI+YLI
        YL=-YLI
        DUMMY=AMAX
        AMAX=AMIN
        AMIN=DUMMY
      ELSE
        X=XI
        Y=YI
        YL=YLI
      END IF
      CALL PLOT(X,Y,3)
      CALL PLOT(X,Y+YL,2)
      Y1=ALOG10(AMIN)
      Y2=ALOG10(AMAX)
      D=YL/(Y2-Y1)
      DTX=0.0
      DTH=DT*H
      HA=ABS(H)
      IF(HA>D) HA=ABS(DT)
      IF(DTH<0.) DTX=ABS(DT)
      DSX=2.5*HA+DTX
      IF(H<0.) DSX=-0.5*HA-DTX
      VS=ASCALE(AMIN,1.0)
      CCC=(VS-AMIN)/AMIN
      IF((CCC>0.99).AND.(CCC<1.01)) VS=AMIN
      CALL ETRANS(VS,XT,NP)
      I=INT(RNDOFF(XT,1))
      XX=X-DSX
      DDY=Y2-Y1
      MAXNP=1
      SGN=1.0
      IF(H<0.) SGN=0.
      IF(DDY<1.5) GO TO 50
   10 IF(I.GE.10) THEN
        I=1
       NP=NP+1
      END IF
      V=FLOAT(I)*(10.**NP)
      IF(V>AMAX) GO TO 20
      DD=DT*0.5
      YY=Y+(ALOG10(V)-Y1)*D
      IF(I==1) THEN
      ISFT=1
      CALL ETRANS(V,XT,MP)
      IF(MP<0) ISFT=ISFT+1
      IAMP=IABS(MP)
      IF(IAMP.GE.10) ISFT=ISFT+1
      IF(ISFT>MAXNP) MAXNP=ISFT
      XXP=XX-0.8*HA*FLOAT(ISFT)*SGN
         DD=DT
        CALL ENUMBR(XXP,YY-0.5*HA,HA,V,0.,-3)
      END IF
      CALL PLOT(X,YY,3)
      CALL PLOT(X+DD,YY,2)
      I=I+1
      GO TO 10
   50 CONTINUE
      DDS=0.05*D*0.66666
      HAH=HA*0.9
      IMAX=9
      IF(DDS<HAH) THEN
        IMAX=8
      HAS=DDS
      ELSE
        HAS=HA*0.7
      END IF
      XXS=X-0.5*HA-DTX-HAS
      IF(H<0.) XXS=X-DSX
   30 IF(I.GE.10) THEN
        I=1
       NP=NP+1
      END IF
      V=FLOAT(I)*(10.**NP)
      IF(V>AMAX) GO TO 20
      DD=DT*0.5
      YY=Y+(ALOG10(V)-Y1)*D
      IF(I==1) THEN
         DD=DT
      ISFT=1
      CALL ETRANS(V,XT,MP)
      IF(MP<0) ISFT=ISFT+1
      IAMP=IABS(MP)
      IF(IAMP.GE.10) ISFT=ISFT+1
      IF(ISFT>MAXNP) MAXNP=ISFT
      XXP=XX-0.8*HA*FLOAT(ISFT)*SGN
        CALL ENUMBR(XXP,YY-0.5*HA,HA,V,0.,-3)
      ELSE IF(I<=IMAX) THEN
      CALL NUMBER(XXS,YY-0.5*HAS,HAS,FLOAT(I),0.,-1)
      END IF
      CALL PLOT(X,YY,3)
      CALL PLOT(X+DD,YY,2)
      I=I+1
      GO TO 30
   20 CONTINUE
C      CALL MSYMCT(NTEXT,NT,ANT,LINES)
      ANT=FLOAT(NT)
      LINES=1
C
      XX=X-(FLOAT(MAXNP)*0.8+3.0)*HA-DTX-(LINES-1)*1.2*HA*12./7.
      IF(H<0.) XX=X+(FLOAT(MAXNP)*0.8+4.4)*HA+DTX
      YY=Y+YL*0.5-0.6*HA*ANT
      IF(NT==0) RETURN
      CALL SYMBOL(XX,YY,HA*1.2,NTEXT,90.,NT)
      RETURN
      END
