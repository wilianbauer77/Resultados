      program main
      parameter (in=1024, nrec = in*in + 1)
      dimension ra(in)
      ind = 98
      open(ind, file='/dk1/dat/strndm.dat',
     &       status='unknown',form='unformatted')
      is = 123456789
      do m = 1, nrec
      do k = 1, in
	ra(k) = strndm(is)
      end do
      write(ind) ra
      end do 
      close(ind)
      stop
      end

C**********************************************************************
      REAL FUNCTION STRNDM(IY)
      X=0.
      DO  K=1,12
        X = X + ran(IY)
      end do
      STRNDM = X - 6.0
      RETURN
      END
