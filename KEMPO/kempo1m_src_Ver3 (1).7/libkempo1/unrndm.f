C**********************************************************************
      REAL FUNCTION UNRNDM(IY)
c      JK=1073741824
c      IY=IY*153+7391
c      IY=IY-JK*(IY/JK)
c      Y=IY
c      X=JK
c      UNRNDM=Y/X
c-- calling CONVEX VECLIB: RAN
      UNRNDM=RAN(IY)
      RETURN
      END
