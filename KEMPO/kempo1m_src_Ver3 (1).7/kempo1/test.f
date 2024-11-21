      dimension str(1024)
      do i = 1, 40
         call strnda(str, 1024)
         print*,str(1000:1024)
      end do
      stop
      end
