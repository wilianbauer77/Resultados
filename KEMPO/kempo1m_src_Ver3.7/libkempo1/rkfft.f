C     *******************************
      SUBROUTINE RKFFT(AR,N,BR,M,INC)
C     *******************************
C     PROGRAMED BY Y. OMURA
C
      DIMENSION AR(N),BR(M)
      NM=N/M
      M2=M/2
      CALL REALFT(AR,N,1)
      RNI=2./FLOAT(N)
      DO 10 I=1,N
   10 AR(I)=AR(I)*RNI
      BR(1)=AR(1)
      BR(2)=AR(2)
      NM2=NM*2
      IF((INC<1).OR.(INC>NM)) INC=NM
      I=1+INC*2
      J=3
      DO 20 L=2,M2
        BR(J)=AR(I)
        BR(J+1)=AR(I+1)
        J=J+2
        I=I+NM2
   20 CONTINUE
      RETURN
      END
