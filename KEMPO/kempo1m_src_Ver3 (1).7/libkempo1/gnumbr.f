C     ***********************************
      SUBROUTINE GNUMBR(X,Y,H,RNB,ANGL,N)
C     ***********************************
C     PROGRAMED BY Y. OMURA
C
      RN=RNB
      AR=ABS(RN)
      IF(AR/=0.) GO TO 50
      CALL NUMBER(X,Y,H,RN,ANGL,1)
      RETURN
   50 CONTINUE
      ND=N-1
      IF(ND<=0) ND=-1
      NE=N+5
      IF(N<=1) NE=NE-1
      CALL ETRANS(AR,A,NP)
      NP1=NP+1
      IF(N<0) GO TO 70
      IF(NP<=-1) THEN
        NC=N+1-NP
        IF(NC<NE) THEN
      RN=RNDOFF(RN,NC-2)
        CALL NUMBER(X,Y,H,RN,ANGL,NC-2)
        ELSE
        CALL ENUMBR(X,Y,H,RN,ANGL,ND)
        END IF
      ELSE IF(NP1.GE.N) THEN
        IF(NP1<NE) THEN
        CALL NUMBER(X,Y,H,RN,ANGL,-1)
        ELSE
        CALL ENUMBR(X,Y,H,RN,ANGL,ND)
        END IF
      ELSE IF(NP1<N) THEN
      RN=RNDOFF(RN,N-NP1)
        CALL NUMBER(X,Y,H,RN,ANGL,N-NP1)
      END IF
      RETURN
   70 ND=-N-1
      IF(ND==0) ND=-1
      IF(NP1>6) THEN
        RN=RNDOFF(RN,ND)
        CALL ENUMBR(X,Y,H,RN,ANGL,ND)
      ELSE
        CALL NUMBER(X,Y,H,RN,ANGL,-1)
      END IF
      RETURN
      END
