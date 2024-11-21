program main
    use param_mod
    use leopardv3_mod
    implicit none
    character(len= *), parameter :: dsha= 'arq/dshark'
    character(len= 3) :: file_num, set_num
    character(len= 18) :: input_file, log_file, dat_file

    real :: start, finish, start2, finish2
    integer :: nk, ik, iw, iarb, narb, comml, datu, nroots, logu
    real :: kstart, kend, dk
    complex :: omega_init, increment
    real, allocatable, dimension(:) :: krange
    complex, allocatable, dimension(:,:) :: solution
    complex, allocatable, dimension(:) :: omega_start, solsk

    ! Get data file number and set up the file names
    if (command_argument_count() == 0) then
        print*, '   ERROR: You have to provide the number of the input data file.'
        print*, '   Command line: $> ./deshark *** (where * = 0 - 9)'
        print*, '   Format of input file: dshark***.set (where * = 0 - 9)'
        stop
    end if

    call get_command_argument(number= 1, value= file_num, length= comml)
    if (comml /= 3) then
        print*, 'ERROR: The data file number must have 3 digits: *** (* = 0 - 9)'
        stop
    end if

    input_file = dsha // file_num // '.set'
    log_file = dsha // file_num // '.log'
    dat_file = dsha // file_num // '.dat'

    open(newunit= datu, file= dat_file, status= 'unknown')
    open(newunit= logu, file= log_file)

    call cpu_time(start)

    write(logu, *) 'Read input data'

    ! Allocate 'omega_init' and 'increment' to avoid runtime error
    ! allocate(omega_init(nroots), increment(nroots))

    call read_data(input_file, set_num, nroots, omega_init, increment, kstart, kend, nk)
    !PRINT*, nk, nroots
    write(logu, *) '...done'

    ! Allocate 'omega_start' after reading 'nroots'
    allocate(omega_start(nroots)); omega_start = omega_init
    allocate(solsk(nroots))
    allocate(solution(nk, nroots))
    allocate(krange(nk))

    dk = (kend - kstart) / (1.0 * nk)

    do ik = 1, nk
        krange(ik) = kstart + (ik - 1) * dk
    end do

    do ik = 1, nk
        write(logu, *) ' '
        write(logu, fmt= '(A7,I6,A10,F12.8)') '--------', ik, '--------k=', krange(ik)

        call cpu_time(start2)

        call muller(logu, nroots, omega_start, krange(ik), solsk)
        !PRINT*, nroots
        solution(ik, :) = solsk
        call cpu_time(finish2)

        ! Write real and imaginary parts of the solution
        ! do iw = 1, nroots
        !     write(logu, '(A9,E20.10,A9,E20.10)') '   omega:', real(solution(ik, iw)), '   gamma:', aimag(solution(ik, iw))
        ! end do
        write(logu, *) 'time elapsed:', finish2 - start2

        if ((ik >= 3) .and. (ik < nk)) then
            ! Ajustar chamada de polyfit para rank correto
            do iw = 1, nroots
                call polyfit(krange(ik-2:ik+1), solution(ik-2:ik, iw), omega_start(iw))
            end do
        else
            ! Para as duas primeiras soluções omega(k), guess próximo valor inicial para a iteração de Muller
            omega_start = solsk + increment
        end if

        ! Write results to output file if within nk limit
        ! if (ik <= nk) then
        !     do iw = 1, nroots
        !         write(datu, fmt= '(F12.8,E20.10,E20.10)') krange(ik), real(solution(ik, iw)), aimag(solution(ik, iw))
        !     end do
            ! end if
        
        ! write(datu, fmt= '(F12.8,E20.10,E20.10)') krange(ik), real(solution(ik, iw)), aimag(solution(ik, iw))
        write(datu, fmt= '(F8.5,*(x, :, es12.5))') krange(ik), solsk

    end do

    call cpu_time(finish)
    write(logu, fmt= '(/, 2(a,g0))') 'Total time elapsed(s): ', finish-start

    close(logu)
    close(datu)

    deallocate(krange, solution, omega_start, solsk)
    ! deallocate(mu,q)
    ! deallocate(beta_para,beta_perp,beta_ratio)
    ! deallocate(dens, drift)

end program main
