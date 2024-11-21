c     ii = 1
c     open(6,file='temp.txt')
cdo i = 1, 102400
c      t = strnd1(ii)
c   write(6,'(f12.6)') t
c     end do
cclose(6)
c     stop
c     end
c**********************************************************************
      REAL FUNCTION STRND1(iy)
	save is
      parameter(in = 1024)
      dimension ra(in)
      data is/0/
	if(iy==0) then
	   rewind 998
	   iy = 1
	   is = 0
      endif
      if(is==0) then
        call strnda(ra,in)
      endif
	is = is + 1
	strnd1=ra(is)
      if(is.ge.in ) is = 0
      return
      END
