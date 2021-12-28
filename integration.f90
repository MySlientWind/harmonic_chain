module integrate
    !--------------------------------------------------------------------
    !purpose:
    !--------------------------------------------------------------------
    use atoms
    implicit none
    real(kind=8) :: dt = 0.001d0
    real(kind=8) :: dh = 0.0005d0
    !============================================================================
contains
    subroutine rk4(x, v)
        !--------------------------------------------------------------------
        !purpose: 4 order Runge Kutta method
        !--------------------------------------------------------------------
        implicit none
        real(kind=8) :: dotx1(natom + 3), dotv1(natom + 3)
        real(kind=8) :: dotx2(natom + 3), dotv2(natom + 3)
        real(kind=8) :: dotx3(natom + 3), dotv3(natom + 3)
        real(kind=8) :: dotx4(natom + 3), dotv4(natom + 3)
        real(kind=8) :: x(natom + 3), v(natom + 3), x0(natom + 3), v0(natom + 3)
        x0 = x; v0 = v
        call ode(x, v, dotx1, dotv1)
        x = x0 + dh*dotx1; v = v0 + dh*dotv1
        call ode(x, v, dotx2, dotv2)
        x = x0 + dh*dotx2; v = v0 + dh*dotv2
        call ode(x, v, dotx3, dotv3)
        x = x0 + dt*dotx3; v = v0 + dt*dotv3
        call ode(x, v, dotx4, dotv4)
        x = x0 + dt*(dotx1 + 2d0*dotx2 + 2d0*dotx3 + dotx4)/6d0
        v = v0 + dt*(dotv1 + 2d0*dotv2 + 2d0*dotv3 + dotv4)/6d0
        return
    end subroutine rk4
    

    subroutine ode(x, v, dotx, dotv)
        !--------------------------------------------------------------------
        !purpose: ordinary differential equations
        !--------------------------------------------------------------------
        use pair; use atoms
        implicit none
        real(kind=8) :: x(natom + 3), v(natom + 3), dotx(natom + 3), dotv(natom + 3)
        integer(kind=4) :: iatom
        
        dotx = v

        dotv(1) = 0d0
        dotv(natom + 3) = 0d0
        do iatom = 2, natom + 2
            dotv(iatom) = pair_force(x(iatom) - x(iatom - 1) - 1d0, k)/m(iatom)&
            &- pair_force(x(iatom + 1) - x(iatom) - 1d0, k)/m(iatom)
        end do
        dotv(natom_reservior + 2) = 0d0
        dotv(natom_reservior + 3) = pair_force(x(natom_reservior + 3) - x(1) - 1d0, k)/m(natom_reservior + 3)&
        &- pair_force(x(natom_reservior + 4) - x(natom_reservior + 3) - 1d0, k)/m(natom_reservior + 3)

        dotv(interact_atom_reservior + 1) = dotv(interact_atom_reservior + 1)&
        &- pair_force(x(interact_atom_system + natom_reservior + 2) - x(interact_atom_reservior + 1)&
        &- dble(interact_atom_system - interact_atom_reservior), k_int)
        dotv(interact_atom_system + natom_reservior + 2) = dotv(interact_atom_system + natom_reservior + 2)&
        &+ pair_force(x(interact_atom_system + natom_reservior + 2) - x(interact_atom_reservior + 1)&
        &- dble(interact_atom_system - interact_atom_reservior), k_int)
        return 
    end subroutine ode
end module integrate