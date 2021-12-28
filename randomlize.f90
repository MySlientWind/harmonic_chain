module randomlize
contains
    function random_v()
        use atoms
        implicit none
        integer(kind=4) :: iatom
        real(kind=8) :: random_v(natom + 3)
        do iatom = 2, natom + 2
            random_v(iatom) = gauss(1d0/sqrt(m(iatom)*beta))
        end do
        random_v(1) = 0d0
        random_v(natom_reservior + 2) = 0d0
        random_v(natom + 3) = v_wall
        return
    end function

    function random_x(x_old)
        !--------------------------------------------------------------------
        !purpose:generate distribution
        !--------------------------------------------------------------------
        use pair; use atoms; use energy_calculate
        implicit none
        real(kind=8) :: delta_x = 0.005d0
        real(kind=8) :: random_x(natom + 3), x_old(natom + 3), x_new(natom + 3)
        real(kind=8) :: r_r, r_x, r_atom, h1, h2
        integer(kind=4) :: atom

        call random_number(r_x)
        call random_number(r_atom)
        atom = floor(r_atom*dble(natom) + 2)
        if (atom > natom_reservior + 1) then
            atom = atom + 1
        end if
        x_new = x_old
        x_new(atom) = 2*delta_x*r_x + x_old(atom) - delta_x

        h1 = potential_energy(x_old) !old
        h2 = potential_energy(x_new) !new

        if (h2 < h1) then
            random_x = x_new
        else
            call random_number(r_r)
            if (r_r < exp(beta*(h1 - h2))) then
                random_x = x_new
            else
                random_x = x_old
            end if
        end if
    end function

    function gauss(sigma)
        !--------------------------------------------------------------------
        !purpose: gaussian distribution
        !--------------------------------------------------------------------
        implicit none
        real(kind=8) :: sigma
        real(kind=8) :: gauss
        real(kind=8), parameter :: pi = 4d0*datan(1d0)

        real(kind=8) :: phi
        real(kind=8) :: gamma

        call random_number(phi)
        call random_number(gamma)
        phi = 2d0*pi*phi
        gamma = -dlog(gamma)
        gauss = sigma*dsqrt(2.0*gamma)*dcos(phi)

        return
    end function

end module randomlize
