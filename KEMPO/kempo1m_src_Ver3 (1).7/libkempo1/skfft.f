C     ***************************
      SUBROUTINE SKFFT(AR,N,BR,M)
C     ***************************
C     PROGRAMED BY Y. OMURA
C
      DIMENSION AR(N),BR(M)
      CALL REALFT(AR,N,1)
      RNI=2./FLOAT(N)
      DO 10 I=1,N
   10 AR(I)=AR(I)*RNI
      DO 20 I=1,M
        BR(I)=AR(I)
   20 CONTINUE
      IF(N>M) BR(2)=AR(M+1)
      BR(1) = 0.5*BR(1)
      RETURN
      END
