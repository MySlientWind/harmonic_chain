include 'atoms.f90'
include 'interaction.f90'
include 'integration.f90'
include 'randomlize.f90'

program main!markov interalated
    use atoms; use integrate; use pair; use energy_calculate; use randomlize; use omp_lib
    implicit none
    integer(kind=4) :: j, i, point_number, ensemble_number = 100, point_number_max = 11!0
    real(kind=8) :: x_init(natom + 3), x(natom + 3), v(natom + 3)
    integer :: n_threads = 10
    real(kind=8) :: e_w_whole = 0d0
    real(kind=8) :: e_w_ave
    real(kind=8) :: e_ini
    real(kind=8) :: deltae
    real(kind=8) :: e_fin
    real(kind=8) :: start, finish
    real(kind=8) :: length_point(21), rho(21), delta_f(21)
    call cpu_time(start)
    call main_init
    
    ! x_init
    call random_number(x_init)
    do i = 1, natom + 3
        x_init(i) = dble(i) - 1.5d0 + x_init(i)
    end do
    x_init(1) = 0d0
    x_init(natom_reservior + 2) = length_reservior
    x_init(natom + 3) = length_system
    x_init(natom_reservior + 3:natom + 2) = x_init(natom_reservior + 3:natom + 2) - length_reservior
    do i = 1, 10000
        x_init = random_x(x_init)
    end do

    ! open file
    open (unit=17, file='./data/time.txt')
    open (unit=15, file='./data/deltaf.txt')

    do point_number = 1, point_number_max
        write (*, *) point_number
        length_point(point_number) = 4d0 - dble(point_number - 1)*0.1d0
        rho(point_number) = dble(natom_system)/length_point(point_number)
        e_w_whole = 0d0
        x = x_init
        !$OMP PARALLEL DO DEFAULT(PRIVATE) SCHEDULE(DYNAMIC) SHARED(point_number)&
        !$OMP PRIVATE(i, j, x, v, x_init, e_ini, e_fin, deltae)&
        !$OMP NUM_THREADS(n_threads) REDUCTION(+:e_w_whole)
        do i = 1, ensemble_number
            x_init = random_x(x_init)
            x = x_init
            v = random_v()
            e_ini = energy(x, v) !before rk, compute each initial energy of every process
            do j = 1, 10000*(point_number - 1)
                call rk4(x, v)
            end do
            e_fin = energy(x, v) !before rk, compute each initial energy of every process
            deltae = e_fin - e_ini
            e_w_whole = e_w_whole + exp(-beta*deltae)
        end do
        !$OMP END PARALLEL DO
        e_w_ave = e_w_whole/dble(ensemble_number)
        delta_f(point_number) = -log(e_w_ave)/beta/natom_system
        write (15, *) rho(point_number), delta_f(point_number)
    end do
    call cpu_time(finish)
    write (17, *) finish - start

contains
subroutine main_init()
    call init_random_seed()
    call atoms_init()
end subroutine main_init

subroutine init_random_seed()
    implicit none
    integer :: i, n, clock
    integer, dimension(:), allocatable :: seed
    call random_seed(size=n)
    allocate (seed(n))
    call system_clock(count=clock)
    seed = clock + 37*(/(i - 1, i=1, n)/)
    call random_seed(PUT=seed)
    deallocate (seed)
end


end program