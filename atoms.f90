module atoms
    !--------------------------------------------------------------------
    !purpose: construct system
    !--------------------------------------------------------------------
    implicit none
    ! atoms
    integer(kind=4), parameter :: natom = 10, natom_reservior = 7, natom_system = 3
    integer(kind=8), parameter :: interact_atom_reservior = 3, interact_atom_system = 2
    real(kind=8), parameter :: length_reservior = dble(natom_reservior + 1), length_system = dble(natom_system + 1)
    real(kind=8), parameter :: k = 1d0, k_int = 0.1d0 !the strength of interaction
    real(kind=8) :: m(natom + 3) ! the mass of particle
    ! the wall
    real(kind=8), parameter :: v_wall = -0.01d0
    ! temperature
    real(kind=8), parameter :: beta = 50d0

contains
    subroutine atoms_init 
        m(1:natom + 3:2) = 1d0
        m(2:natom + 3:2) = sqrt(2d0)
    end subroutine atoms_init

end module atoms