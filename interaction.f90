module pair
    !--------------------------------------------------------------------
    !purpose: interaction
    !--------------------------------------------------------------------
    implicit none
contains
    !============================================================================
    function pair_potential(x, k)
        !--------------------------------------------------------------------
        !purpose:
        !--------------------------------------------------------------------
        implicit none
        real(kind=8) :: x, pair_potential, k
        pair_potential = k*x**2/2
        return
    end function
    !============================================================================
    function pair_force(x, k)
        !--------------------------------------------------------------------
        !purpose:
        !--------------------------------------------------------------------
        implicit none
        real(kind=8) :: x, pair_force, k
        pair_force = -k*x
        return
    end function
end module pair

module energy_calculate
    !--------------------------------------------------------------------
    !purpose:
    !--------------------------------------------------------------------
    implicit none
contains
    function kinetic_energy(v_e)
        !--------------------------------------------------------------------
        !purpose:
        !--------------------------------------------------------------------
        use atoms; 
        implicit none
        real(kind=8) :: v_e(:)
        real(kind=8) :: kinetic_energy
        kinetic_energy = 0.5d0*sum(v_e**2*m)
        return
    end function
    !============================================================================
    !============================================================================
    function potential_energy(x_e)
        !--------------------------------------------------------------------
        !purpose:
        !--------------------------------------------------------------------
        use atoms; use pair

        implicit none
        integer(kind=4) :: iatom
        real(kind=8) :: x_e(:)
        real(kind=8) :: potential_energy

        potential_energy = pair_potential(x_e(2) - x_e(1) - 1d0, k)
        do iatom = 2, natom_reservior + 1
            potential_energy = potential_energy + pair_potential(x_e(iatom + 1) - x_e(iatom) - 1d0, k)
        end do
        potential_energy = potential_energy + pair_potential(x_e(natom_reservior + 3) - 1d0, k)
        do iatom = natom_reservior + 3, natom + 2
            potential_energy = potential_energy + pair_potential(x_e(iatom + 1) - x_e(iatom) - 1d0, k)
        end do
        potential_energy = potential_energy + pair_potential(x_e(interact_atom_system + natom_reservior + 2)&
        & - x_e(interact_atom_reservior + 1) &
        & - dble(interact_atom_system - interact_atom_reservior), k_int)
        return
    end function

    function energy(x_e, v_e)
        !--------------------------------------------------------------------
        !purpose: compute average energy
        !--------------------------------------------------------------------
        use atoms; use pair
        implicit none
        real(kind=8) :: x_e(:), v_e(:)
        real(kind=8) :: energy
        energy = potential_energy(x_e) + kinetic_energy(v_e)
        return
    end function
end module energy_calculate
