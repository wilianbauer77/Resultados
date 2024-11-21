C**********************************************************************
      REAL FUNCTION STRNDM(IY)
      XX=0.
      DO 10 K=1,12
       XX = XX + ran(IY)
   10 CONTINUE
      STRNDM = XX - 6.0
      RETURN
      END
