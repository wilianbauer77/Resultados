      program main
      parameter (in=1024, nrec = in*in+1)
      dimension ra(in)
      ind = 98
      open(ind, file='/dk1/dat/unrndm.dat',
     &            status='unknown',form='unformatted')
      is = 123456789
      do m = 1, nrec
      do k = 1, in
	ra(k) = ran(is)
      end do
      write(ind) ra
      end do 
      close(ind)
      stop
      end
