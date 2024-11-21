C     ***********************************
      SUBROUTINE MAXMIN (XXX,N,AMIN,AMAX)
C     ***********************************
      DIMENSION XXX(N)
C
      AMAX=XXX(1)
      AMIN=XXX(1)
      DO 10 I=1,N
      IF(XXX(I).GE.AMAX) AMAX=XXX(I)
      IF(XXX(I)<=AMIN) AMIN=XXX(I)
   10 CONTINUE
C
      RETURN
      END
