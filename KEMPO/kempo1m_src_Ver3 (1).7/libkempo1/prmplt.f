C     *************************************************
      SUBROUTINE PRMPLT(XI,YI,HI,ANGLI,NTEXT,NT,RN,NEN)
C     *************************************************
C     PROGRAMED BY Y. OMURA
C
      CHARACTER NTEXT*(*)
      DATA X,Y,H,ANGL/0.,0.,1.,0./
      IF((XI==999).OR.(YI==999.)) THEN
        X=X-H*SIN(ANGL/180.*3.14159)*2.0
        Y=Y+H*COS(ANGL/180.*3.14159)*2.0
      ELSE IF(XI==-999.) THEN
        X=X+H*SIN(ANGL/180.*3.14159)*2.0
        Y=Y-H*COS(ANGL/180.*3.14159)*2.0
      ELSE
        X=XI
        Y=YI
        H=HI
        ANGL=ANGLI
      ENDIF
      IAB=IABS(NEN)
      IF(IAB.GE.1000) THEN
      ID=INT(FLOAT(IAB)/100.)*100
      FACT=FLOAT(ID)/1000.
      NE=NEN/IAB*(IAB-ID)
      ELSE
      FACT=1.
      NE=NEN
      END IF
      CALL SYMBOL(X,Y,H*FACT,NTEXT,ANGL,NT)
C      CALL MSYMCT(NTEXT,NT,AT,LINES)
      AT=FLOAT(NT)
      LINES=1
      XX=X+H*FACT*AT+H
      YY=Y
      CALL XYROT(XX,YY,X,Y,ANGL)
      CALL SYMBOL(XX,YY,H,'=',ANGL,1)
      XX=X+H*FACT*AT+H*3.
      YY=Y
      CALL XYROT(XX,YY,X,Y,ANGL)
      IF(NE>0) THEN
      CALL GNUMBR(XX,YY,H,RN,ANGL,NE)
      ELSE IF( NE==0) THEN
      CALL NUMBER(XX,YY,H,RNDOFF(RN,-1),ANGL,-1)
      ELSE
      CALL NUMBER(XX,YY,H,RNDOFF(RN,-NE),ANGL,-NE)
      ENDIF
      RETURN
      END
