c     real*8 a(100)
c     do 10 i=100,1,-1
c      a(i) = i
c 10  continue
c     call qsort(a,100)
c     print*,a
c     stop
c     end
c     
c     Copyright (C) AKIRA SAWADA, 1989
c     programed by sawada@kimura4.kuee.kyoto-u.ac.jp
c     Everyone is permitted to do anything on this program
c     including copying, transplanting, debugging and modifying.
c     
      
************************************************************************
      subroutine qsort(x,n)
************************************************************************
      real*4  x(n)
      parameter (ipmax=20)
      dimension il(ipmax),ir(ipmax)
      
      ip     = 1
      il(ip) = 1
      ir(ip) = n

   10 continue
      if(ip>0) then
        
        if(il(ip).ge.ir(ip)) then
          ip = ip-1
        else
          i  = il(ip)
          j  = ir(ip)
          k  = (i+j)/2
          p = max(min(x(i),x(j)),min(x(j),x(k)),min(x(k),x(i)))
          
   20     continue
          if(i<=j) then
            
   30       continue
            if(x(i)<p) then
              i    = i+1
              goto 30
            end if

   40       continue
            if(x(j)>p) then
              j    = j-1
              goto 40
            end if
            
            if(i<=j) then
              t    = x(i)
              x(i) = x(j)
              x(j) = t
              i    = i+1
              j    = j-1
            end if
            
            goto 20
          end if
          
          if((j-il(ip))<(ir(ip)-i)) then
            il(ip+1) = il(ip)
            ir(ip+1) = j
            il(ip  ) = i
          else
            il(ip+1) = i
            ir(ip+1) = ir(ip)
            ir(ip  ) = j
          end if
          ip = ip+1
        end if
        
        goto 10
      end if
      
      return
      end
