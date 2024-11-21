C     *************************
      SUBROUTINE ETRANS(X,A,NP)
C     *************************
C     PROGRAMED BY Y. OMURA
C
      X0=ABS(X)
      IS=1
      IF(X<0.) IS=-1
      A=0.
      NP=0
      I=0
      IF(X0==0.) RETURN
   10 IF(X0.GE.1.) GO TO 20
      I=I-1
      X0=X0*10.
      GO TO 10
   20 IF(X0<9.999) GO TO 30
      I=I+1
      X0=X0/10.
      GO TO 20
   30 A=X0*FLOAT(IS)
      NP=I
      RETURN
      END
