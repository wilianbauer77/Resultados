C     ***************************
      REAL FUNCTION ASCALE(X,RNI)
C     ***************************
C     PROGRAMED BY Y. OMURA
C
      IF(ABS(RNI).GE.1000) THEN
        RN=INT(RNI/1000.)*(ABS(RNI)-1000.)
        ARN=ABS(RN)
        ICONT=1
      ELSE
        RN=RNI
        ARN=ABS(RN)
        ICONT=0
      ENDIF
      IF(RN/=0.) GO TO 10
      ASCALE=X
      RETURN
   10 CALL ETRANS(X,A,NP)
      AA=ABS(A)
      IF(ICONT==1.AND.A<0.) THEN
      IS=1
      IF(X<0.) IS=-1
      OPT=INT(AA/ARN)*ARN
      IF(RN<0.) OPT=OPT+ARN
      ASCALE=(OPT*10.**NP)*FLOAT(IS)
      ELSE
      IS=1
      IF(X<0.) IS=-1
      OPT=INT(AA/ARN)*ARN+ARN
      IF(RN<0.) OPT=OPT-ARN
      ASCALE=(OPT*10.**NP)*FLOAT(IS)
      ENDIF
      RETURN
      END
