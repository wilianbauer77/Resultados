C     *********************************
      SUBROUTINE ENUMBR(X,Y,H,R,ANGL,N)
C     *********************************
C     PROGRAMED BY Y. OMURA
C
      R0=ABS(R)
      IS=1
      IF(R<0.) IS=-1
      I=0
      IF(R0==0.) GO TO 40
   10 IF(R0.GE.1.) GO TO 20
      I=I-1
      R0=R0*10.
      GO TO 10
   20 IF(R0<10.) GO TO 30
      I=I+1
      R0=R0/10.
      GO TO 20
   30 NP=2+N
      NRD=N
      IF(NRD<1) NRD=1
      R0=RNDOFF(R0,NRD)
      IF(R0<10.) GO TO 40
      I=I+1
      R0=R0/10.
   40 CONTINUE
      R0=R0*FLOAT(IS)
      EI=FLOAT(I)
      HF=H*0.8
      IF(N==-3) GO TO 70
      IF(N<=-2) GO TO 50
      CALL NUMBER(X,Y,H,R0,ANGL,N)
   50 IF(I==0) RETURN
      NP=2+N
      IF(IS==-1) NP=NP+1
      XP=X+H*FLOAT(NP)
      YP=Y
      CALL XYROT(XP,YP,X,Y,ANGL)
      CALL SYMBOL(XP,YP,HF,'X',ANGL,1)
      XP=X+H*FLOAT(NP)+HF
      YP=Y
      CALL XYROT(XP,YP,X,Y,ANGL)
      CALL SYMBOL(XP,YP,H,'10',ANGL,2)
      XP=X+H*FLOAT(NP+2)+HF
      YP=Y+0.5*H
      CALL XYROT(XP,YP,X,Y,ANGL)
      CALL NUMBER(XP,YP,HF,EI,ANGL,-1)
      RETURN
   70 CONTINUE
      XP=X
      YP=Y
      CALL SYMBOL(XP,YP,H,'10',ANGL,2)
      XP=X+H*2.
      YP=Y+0.5*H
      CALL XYROT(XP,YP,X,Y,ANGL)
      CALL NUMBER(XP,YP,HF,EI,ANGL,-1)
      RETURN
      END
