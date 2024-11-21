C***************************************************************
C FFT OF SINGLE REAL FUNCTION
C   INPUT DATA HAVE 2N ELEMENTS
C     ISIGN = 1 FOR FOURIER TRANSFORM
C          TRANSFORMD DATA SHOULD BE MULTIPLIED BY 1/N
C          WITH SEQUANCE OF
C          C(0),C(N),C(1),S(1),C(2),S(2),......C(N-1),S(N-1)
C       WHILE X(J) : J=1,2,....,2N IS EXPRESSED AS
C          X(J)= 0.5*C(0) + SUM(C(K)COS(..)+S(K)SIN(..)) + 0.5*C(N)
C
C     ISIGN = -1 FOR INVERSE FOURIER TRANSFORM
C     
C     REFERENCE: NUMERICAL RECIPES BY W.H. PRESS ET AL., CAMBRIDGE 1986
C        MODIFIED BY Y. OMURA, SEPTEMBER, 1989
C
      SUBROUTINE REALFT(DATA,N2,ISIGN)
      REAL*8 WR,WI,WPR,WPI,WTEMP,THETA
      DIMENSION DATA(N2)
      N=N2/2
      THETA=3.141592653589793D0/DBLE(N)
      C1=0.5
      IF(ISIGN==1) THEN
        C2=-0.5
        CALL FOUR1(DATA,N2,1)
      ELSE
        C2=0.5
        THETA=-THETA
      ENDIF
      WPR=-2.D0*DSIN(0.5D0*THETA)**2
      WPI=DSIN(THETA)
      WR=1.D0+WPR
      WI=WPI
      N2P3=2*N+3
      DO 11 I1=3,N-1,2
        WRS=SNGL(WR)
        WIS=SNGL(WI)
        H1R=C1*(DATA(I1)+DATA(N2P3-I1-1))
        H1I=C1*(DATA(I1+1)-DATA(N2P3-I1))
        H2R=-C2*(DATA(I1+1)+DATA(N2P3-I1))
        H2I=C2*(DATA(I1)-DATA(N2P3-I1-1))
        DATA(I1)=H1R+WRS*H2R-WIS*H2I
        DATA(I1+1)=H1I+WRS*H2I+WIS*H2R
        DATA(N2P3-I1-1)=H1R-WRS*H2R+WIS*H2I
        DATA(N2P3-I1)=-H1I+WRS*H2I+WIS*H2R
        WTEMP=WR
        WR=WR*WPR-WI*WPI+WR
        WI=WI*WPR+WTEMP*WPI+WI
   11 CONTINUE
      IF(ISIGN==1) THEN
        H1R=DATA(1)
        DATA(1)=H1R+DATA(2)
        DATA(2)=H1R-DATA(2)
      ELSE
        H1R=DATA(1)
        DATA(1)=C1*(H1R+DATA(2))
        DATA(2)=C1*(H1R-DATA(2))
        CALL FOUR1(DATA,N2,-1)
      ENDIF
      RETURN
      END

