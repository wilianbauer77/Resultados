!> Scans through the requested wavenumber interval, computes corresponding frequencies and writes them to output file
!program main
program main
  use param_mod
  implicit none
  character(len= *), parameter :: dsha= 'arq/dshark'
  character(len= 3) :: file_num, set_num
  character(len= 18) :: input_file, log_file, dat_file

  !complex :: omega_start, increment
  real :: start, finish, start2, finish2
  integer :: nk, ik, iw, iarb, narb, comml, datu, nroots, logu

  real :: kstart, kend, dk
  real, allocatable, dimension (:) :: krange, splcoeff1, splcoeff2
  complex, allocatable, dimension (:,:) :: solution, omega_start
  complex, allocatable, dimension(:) :: omega_init, increment, solsk

  ! Get data file number and set up the file names.
   if(command_argument_count() == 0)then
      print*, '   ERROR: You have to provide the number of the imput data file.'
      print*, '   Command line: $>./deshark *** (where * = 0 - 9)'
      print*, '   Format of input file: dhark***.set (where * = 0 - 9)'
      stop
   end if
   call get_command_argument(number= 1, value= file_num, length= comml)
   if(comml /= 3)then
      print*, 'ERROR: The data file number must have 3 digits: *** (* = 0 - 9)'
      stop
   end if
   input_file= dsha//file_num//'.set'
   log_file=   dsha//file_num//'.log'
   dat_file=   dsha//file_num//'.dat'

   open(newunit= datu, file= dat_file, status= 'unknown')
   open(newunit= logu, file= log_file)

   call cpu_time(start)

   write(logu, *) 'Read input data'
   call read_data(input_file, set_num, nroots, omega_start, increment, kstart, kend, nk)
   write(logu, *) '...done'
   allocate(omega_start(nroots)) ; omega_start= omega_init
   allocate(solsk(nroots))
   allocate(solution(nk, nroots))
   allocate(krange(nk))

   dk=(kend-kstart)/(1.0*nk)
 
   do ik=1,nk
      krange(ik)=kstart+(ik-1)*dk
   enddo

  !open(unit=7,status='unknown',file='omega.dat')

  !call cpu_time(start)

  !read parameters from input.dat
   !allocate(splcoeff1(npara_max-1,nperp_max-1,4,3,narb))
   !allocate(splcoeff2(npara_max-1,nperp_max-1,4,3,narb))

   !do iarb=1, narb
   !   call get_splinecoeff(iarb,splcoeff1(:,:,:,:,iarb),splcoeff2(:,:,:,:,iarb))
   !enddo


   do ik=1, nk
      write(logu, *) ' '
      write(logu, fmt= '(A7,I6,A10,F12.8)') '--------',ik,'--------k=', krange(ik)

      call cpu_time(start2)

      call muller(logu, nroots, omega_start, krange(ik), solsk)
      solution(ik, :)= solsk
      call cpu_time(finish2)

      write(logu,'(A9,E20.10,A9,E20.10)')  '   omega:', real(solution(ik)), '   gamma:',aimag(solution(ik))
      write(logu, *) 'time elapsed:', finish2-start2


      if ((ik.ge.3).and.(ik<nk))  then

        !if three subsequent solutions omega(k) are found, use quadratic polynomial fit 
        !to guess next starting frequency for Muller iteration
        
         do iw= 1, nroots
            call polyfit(krange(ik-2:ik+1), solution(ik-2:ik, iw), omega_start(iw))
         end do
         !call polyfit(krange(ik-2:ik+1),solution(ik-2:ik),omega_start)

     else

        !for the first two solution omega(k) guess next starting frequency for Muller iteration
        !by raising the computed omega by an increment which is provided by the user
        omega_start=solution(ik)+increment

     end if
     write(datu, fmt= '(F12.8,E20.10,E20.10)') krange(ik), real(solution(ik)), aimag(solution(ik))
   enddo


  call cpu_time(finish)
 
  write(logu,*) 'Total time elapsed: ', finish-start
 
  close(logu) ; close(datu)


  deallocate(krange,solution, omega_start, solsk)
  deallocate(mu,q,kappa,dens)
  deallocate(beta_para,beta_perp,beta_ratio)
  !deallocate(splcoeff1,splcoeff2)


end program main
