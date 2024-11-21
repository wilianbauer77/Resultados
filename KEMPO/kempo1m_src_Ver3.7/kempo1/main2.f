        program main2
        include "paramt.h"
        common /prtclc/ x(in), vx(in), vy(in),vz(in)
        common /constc/ tcs, bx0, rho0, slx, nx, nxp1, nxp2, npt, ns
        common /ptprmc/ wp(is), qm(is), q(is), vpe(is), vpa(is),
     #                  vd(is), pch(is), np(is)
        common /fieldc/ ex(ix), ey(ix), ez(ix), by(ix), bz(ix),
     #                  ajx(ix), ajy(ix), ajz(ix), rho(ix)
        common /work1c/ work1(ix),work2(ix)

        call input
        call renorm
        call inital
        call positn
        call charge
        call ecrrct

        do i=1,10
          call bfield
          call velcty
          call positn
          call currnt
          call positn
          call bfield
          call efield
	  call charge
	  call ecrrct
        enddo
        
        do i=1,nx
          write(*,10)ex(i),ey(i),ez(i),rho(i)
10        format(' ',4f14.10)
        enddo

        stop
        end
