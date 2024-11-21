      dimension wk1(10), wk2(10)
c
      do 10  i = 1, 10
       wk1(i) = i
       wk2(i) = -i
   10 continue      
c
      open(50,file='eng.dat',status='unknown')
      write(50,*) wk1, wk2
      close(50)
      stop
      end

