module Parameters
    ! Contains constants and parameters for the calculation.
    ! Constants are in atomic units (atu/AU) unless specified.
    integer, parameter :: xr = kind(0.d0)
    real(kind = xr), parameter :: m_ele = 1.0_xr ! mass of a electron
    real(kind = xr), parameter :: h_bar = 1.0_xr ! Planck constant
    real(kind = xr), parameter :: pi = 3.141592653587932385_xr
    real(kind = xr), parameter :: ev_to_hartree = 27.211_xr
    real(kind = xr), parameter :: ev_to_erg = 6.242e+11_xr
    real(kind = xr), parameter :: a0_sqr_to_mb = 100._xr * 0.529177_xr ** 2.0_xr ! 1 Bohr^2 to MegaBarn
    real(kind = xr), parameter :: atu_to_as = 24.188843265857_xr ! Atomic unit time to atto-second
    complex(kind = xr), parameter :: i = (0.0_xr, 1.0_xr)

    abstract interface
        function function_template(r)
            real(kind = 8) :: r
            real(kind = 8) :: function_template
        end function
    end interface

end module Parameters