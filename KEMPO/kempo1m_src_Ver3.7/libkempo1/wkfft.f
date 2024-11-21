C     *******************************************
      SUBROUTINE WKFFT(AR,N1,M1,N,M,WK1,WK2,ICNT)
C     *******************************************
C
C-------------BY Y.OMURA   RASC, KYOTO UNIV.  -----
C--- SUBROUTINE TO FOURIER TRANSFORM IN SPACE AND TIME ----
C---  ICNT=0  FFT IN BOTH X(Z) AND Y(T) COMPORNENTS ----
C---  ICNT=1  FFT IN Y(T) COMPORNENT         ----
C---  ICNT=2   FFT IN X(Z) COMPORNENT    ---
C
      DIMENSION AR(N1,M1)
      DIMENSION WK1(N1),WK2(M1)
      DATA PI/3.1415926/
C
      RNI=2./FLOAT(N)
      RMI=2./FLOAT(M)
      IF(ICNT==1) GO TO 35
      DO 30 J=1,M
        DO 10 I=1,N
   10   WK1(I)=AR(I,J)
        CALL REALFT(WK1,N,1)
        DO 20 I=1,N
   20   AR(I,J)=WK1(I)*RNI
        AR(1,J)=0.5*AR(1,J)
        AR(2,J)=0.5*AR(2,J)
   30 CONTINUE
   35 CONTINUE
      IF(ICNT==2) GO TO 45
      DO 40 I=1,N
        DO 50 J=1,M
   50   WK2(J)=AR(I,J)
        CALL REALFT(WK2,M,1)
        DO 60 J=1,M
   60   AR(I,J)=WK2(J)*RMI
   40 CONTINUE
   45 CONTINUE
      DO 91 I=1,N
        AR(I,1)=ABS(AR(I,1))*0.5
        AR(I,2)=ABS(AR(I,2))*0.5
   91 CONTINUE
      N2=N/2
      M2=M/2
      DO 64 I=1,2
      J=3
      DO 65 L=2,M2
      AR1=AR(I,J)
      AR2=AR(I,J+1)
      SQ=AR1*AR1+AR2*AR2
      ARA=SQRT(SQ)
      IF(ARA==0.) ARA=0.0001
      T1=ACOS(AR2/ARA)
      IF(AR1<0.) T1=T1+PI
      AR(I,J)=ARA
      AR(I,J+1)=T1
      J=J+2
   65 CONTINUE
   64 CONTINUE
      DO 66 J=1,2
      I=3
      DO 67 L=2,N2
      AR1=AR(I,J)
      AR2=AR(I+1,J)
      SQ=AR1*AR1+AR2*AR2
      ARA=SQRT(SQ)
      AR(I,J)=ARA
      AR(I+1,J)=ARA
      I=I+2
   67 CONTINUE
   66 CONTINUE
      J=3
      DO 70 L=2,M2
      I=3
      DO 80 K=2,N2
        CC=AR(I,J)
        CS=AR(I,J+1)
        SC=AR(I+1,J)
        SS=AR(I+1,J+1)
        SQ=(CS-SC)**2+(CC+SS)**2
        AR(I,J)=0.5*SQRT(SQ)
        SQ=(CS+SC)**2+(CC-SS)**2
        AR(I+1,J)=0.5*SQRT(SQ)
        AR1=AR(I,J)
        IF(AR1==0.) AR1=0.0001
        AR2=AR(I+1,J)
        IF(AR2==0.) AR2=0.0001
        TC1=0.5*(CS-SC)/AR1
        TC2=0.5*(CS+SC)/AR2
        T1=ACOS(TC1)
        TSIGN=CC+SS
        IF(TSIGN<0.) T1=T1+PI
        T2=ACOS(TC2)
        TSIGN=CC-SS
        IF(TSIGN<0.) T2=T2+PI
        AR(I,J+1)=T1
        AR(I+1,J+1)=T2
        I=I+2
   80 CONTINUE
      J=J+2
   70 CONTINUE
      RETURN
      END
