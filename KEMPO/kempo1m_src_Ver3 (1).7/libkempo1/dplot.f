C     ****************************
      SUBROUTINE DPLOT(XP,YP,IPEN)
	save
      ENTRY DPLOT6(XP,YP,IPEN,DA,BA,PA)
C     ****************************
C     COPIED FROM XYGRAPH BY T. SATO, DEPT. OF E.E., KYOTO UNIV.
C     MODIFIED BY Y. OMURA
C
      DIMENSION A(102),M(102)
      DATA MOVE /-1/ 
      DATA MODE/-1/
      DATA XPM,YPM,D,B,P /5*0./
C
      IF(IPEN<=1) GO TO 100
      IF(MODE.GE.0) GO TO 200
      CALL PLOT(XP,YP,IPEN)
      RETURN
  100 N=6
      IF(IPEN==-999) THEN
      XP=XPM
      YP=YPM
      IPEN=-MODE
      DA=D
      BA=B
      PA=P
      RETURN
      END IF
      MODE=-IPEN
      CALL PLOT(XP,YP,3)
      IF(MODE==-1) RETURN
      IF(MODE>50) MODE=50
      D=0.5
      IF(N.GE.4) D=DA
      B=D*0.5
      IF(MODE.GE.1) B=D*0.2
      IF(N.GE.5) B=BA
      P=B
      IF(N.GE.6) P=PA
      XPM=XP
      YPM=YP
      T=0.
      NA=(MODE+1)*2
      A(1)=D
      A(2)=D+B
      C=A(2)
      M(1)=2
      M(2)=3
      IF(MODE==0) RETURN
      DO 110 I=3,NA,2
      A(I)=A(I-1)+P
      A(I+1)=A(I)+B
      M(I)=2
      M(I+1)=3
  110 CONTINUE
      C=A(NA)
      RETURN
  200 IF(XP==XPM.AND.YP==YPM) RETURN
      S=SQRT((XP-XPM)**2+(YP-YPM)**2)
      IF(IPEN==2) GO TO 300
      IF(MOVE==1) T=AMOD(S+T,C)
      IF(MOVE==-1) T=0.
      XPM=XP
      YPM=YP
      CALL PLOT(XP,YP,3)
      RETURN
  300 DX=(XP-XPM)/S
      DY=(YP-YPM)/S
      DO 10 I=1,NA
      IF(A(I).GE.T) GO TO 20
   10 CONTINUE
      I=NA
   20 IF(T+S>C) GO TO 400
      DO 210 J=I,NA
      IF(T+S<=A(J)) GO TO 220
      CALL PLOT(XPM+DX*(A(J)-T),YPM+DY*(A(J)-T),M(J))
  210 CONTINUE
      J=NA
  220 CALL PLOT(XP,YP,M(J))
      T=T+S
      XPM=XP
      YPM=YP
      RETURN
  400 DO 30 J=I,NA
   30 CALL PLOT(XPM+DX*(A(J)-T),YPM+DY*(A(J)-T),M(J))
      XPM=XPM+DX*(C-T)
      YPM=YPM+DY*(C-T)
      S=S+T-C
      T=0.
      L=S/C
      IF(L==0) GO TO 500
      S=S-L*C
      DO 40 I=1,L
      DO 50 J=1,NA
   50 CALL PLOT(XPM+DX*A(J),YPM+DY*A(J),M(J))
      XPM=XPM+DX*C
      YPM=YPM+DY*C
   40 CONTINUE
  500 DO 60 I=1,NA
      IF(A(I).GE.S) GO TO 70
   60 CONTINUE
      I=NA
   70 IF(I==1) GO TO 80
      DO 90 J=1,I-1
   90 CALL PLOT(XPM+DX*A(J),YPM+DY*A(J),M(J))
   80 CALL PLOT(XP,YP,M(I))
      T=S
      XPM=XP
      YPM=YP
      RETURN
      END
