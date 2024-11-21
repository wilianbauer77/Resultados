      real FUNCTION BITRVS(IC,N)
c     implicit real*4 (a-h,o-z)
C     ********************
C--- BIT-REVERS-NUMBER IS GENERATED --
C---  SET IC=1 TO INITILIZE THE FUNCTION --
C---  N=2**I  : NUMBER OF THE NUMBERS TO BE GENERATED --
C
      DATA I/1/
      IF(IC==1) I=N
      XS=0.
      F=0.5
      IX=I
   10 XS=XS+(IX-INT(IX*0.5)*2)*F
      IX=IX/2
      IF(IX<=0) GO TO 15
      F=F*0.5
      GO TO 10
   15 IC=IC+1
      I=I+1
      IF(IC>N) IC=1
      BITRVS=XS
      RETURN
      END
