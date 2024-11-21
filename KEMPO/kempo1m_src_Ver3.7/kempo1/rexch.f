c     ***********************
      subroutine rexch(n,xxx)
c     ***********************
c --- exchange at random ---
      dimension xxx(n)
      data ii/0/
c
      nexch=2*n
      do 10 i=1,nexch
        n1=int(n*unrndm(ii))
        n2=int(n*unrndm(ii))
        if(n1<=0) go to 10
        if(n2<=0) go to 10
        if(n1>n) go to 10
        if(n2>n) go to 10
c
        a1=xxx(n1)
        xxx(n1)=xxx(n2)
        xxx(n2)=a1
   10 continue
      return
      end
