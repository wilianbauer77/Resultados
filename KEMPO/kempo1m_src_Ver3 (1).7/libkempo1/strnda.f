c      dimension ra(128)
c      do i = 1, 20
c       call strnda(ra,128)
c	  print*,ra(1), ra(2)
c      end do 
c      stop
c      end
c 
c     ========================
      subroutine strnda(ra, n)
      save is, irec
c     ========================
c
      parameter (in = 1024,  nrec = in*in + 1 )
      dimension ra(n)
      data is, irec /0,0/
      ind = 998
c
      if(n>in) then
         print*, 'too large record length is specified ', n
      endif
      if(is==0) then
	open(ind, file='/dk1/dat/strndm.dat',
     &           action='read',status='old',form='unformatted')
	is = 1
      end if
      irec = irec + 1
      read(ind) ra
      if(irec==nrec) then
        rewind ind
        irec = 0
      endif
      return
      end
