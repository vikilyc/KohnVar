! Demonstration of the Complex Kohn Variational method for solving a sample one electron s-wave scattering 
! problem from a short-range potential that essentially reproduces some of the resulte from the following 
! reference, only the 1D case is modeled by this script
!
! McCurdy, C. W., Rescigno, T. N., & Schneider, B. I. (1987). 
! Interrelation between variational principles for scattering amplitudes and generalizedr-matrix theory. 
! Physical Review A, 36(5), 2061–2066. https://doi.org/10.1103/physreva.36.2061 
!
! Author: Yuchen Liu
! Date: 13 June 2022

program scatter
    use Parameters
    use Integrators
    implicit none
    integer ::n_grid_pts, n_basis, idx1, idx2, pos, info, counter, iter
    real(kind = xr), allocatable, dimension(:) :: x, weights, j, j_d2dx2, w, rwork
    complex(kind = xr), allocatable, dimension(:) :: h_plus, h_plus_d2dx2, c, s, e_eigvals, work, swork
    complex(kind = xr), allocatable, dimension(:, :) :: psi, psi_d2dx2, M, hamiltonian, vs, M_copy
    real(kind = xr) :: slater_lambda, k, r_max, energy, r_min, slater_density
    complex(kind = xr) :: lambda_result
    integer, allocatable, dimension(:) :: ipiv


    procedure(function_template), pointer :: v_pot => null ()
    n_grid_pts = 103
    ! Number of discrete grid points used for integration by Simpson's method
    ! Must be odd number

    k = 0.15_xr 
    ! Sample electron momentum for elastic scattering

    slater_lambda = 2.5_xr   
    ! Decay rate of Slater type basis functions of our choice in this demonstration
    
    r_max = 30.0_xr 
    r_min = 1e-8_xr
    ! Upper and lower bound for the range of our model, r_min != 0 to avoid zero-division issues

    energy = (k ** 2._xr) / 2._xr
    ! Energy in Hartree (Atomic unit) for the electron with momentum k

    n_basis = 8
    ! Number of Slater basis functions used in our basis set 

    v_pot => v_pot_sr

    allocate (x(n_grid_pts))
    allocate (weights(n_grid_pts))
    call SimpsonWeights(x, weights, r_min, r_max, n_grid_pts)
        ! Initialization of integration points and their associated weight in our 1D model

    allocate (j(n_grid_pts))            ! Incoming boundary condition wave function on the grids
    allocate (j_d2dx2(n_grid_pts))      ! The second derivative of j's
    allocate (h_plus(n_grid_pts))       ! Outgoing boundary condition wave function on the grids
    allocate (h_plus_d2dx2(n_grid_pts)) ! The second derivative h+'s on the grids

    do pos = 1, n_grid_pts
        ! Generating j, h+ and their second derivatives from analytic functions
        ! Possible to use finite differences here for the derivatives
        j(pos) = sin(k * x(pos))
        j_d2dx2(pos) = - sin(k * x(pos)) * (k ** 2.0_xr)
        h_plus(pos) = zexp(i * k * x(pos)) * g(x(pos))
        h_plus_d2dx2(pos) = zexp(i * k * x(pos)) * (g_d2dx2(x(pos)) &
                            & + i * 2._xr * k * g_ddx(x(pos))  &
                            & - (k ** 2._xr) * g(x(pos)))
    end do

    allocate (psi(0 : n_basis, n_grid_pts))
    allocate (psi_d2dx2(0 : n_basis, n_grid_pts))

    psi(0, :) = h_plus(:)
    psi_d2dx2(0, :) = h_plus_d2dx2(:)
    
    ! generating the trial functions (Slater basis functions)
    if (n_basis > 0) then
        do idx1 = 1, n_basis 
            slater_density = 0.0_xr 
            ! The integration values of each slater trial functions used for normalization
            do pos = 1, n_grid_pts
                psi(idx1, pos) = slater(x(pos), idx1, slater_lambda)
                slater_density = slater_density + real(psi(idx1, pos)) ** 2.0_xr * weights(pos)
                psi_d2dx2(idx1, pos) = slater_d2dx2(x(pos), idx1 , slater_lambda)
            end do

            do pos = 1, n_grid_pts
                ! Normalization of trial functions
                psi(idx1, pos) = psi(idx1, pos) / sqrt(slater_density)
                psi_d2dx2(idx1, pos) = psi_d2dx2(idx1, pos) / sqrt(slater_density)
            end do
        end do
    end if

    allocate (hamiltonian(0 : n_basis , 0 : n_basis))
    allocate (vs(0 : n_basis, 0 : n_basis))
    allocate (w(0 : n_basis))
    allocate (e_eigvals(0 : n_basis))
    allocate (rwork(0 : n_basis))

    ! Generating the Hamiltonian matrix with the trial function psi
    do idx1 = 0, n_basis 
        do idx2 = 0, n_basis
            hamiltonian(idx1, idx2) =  complex(0.0_xr, 0.0_xr)
            do pos = 1, n_grid_pts
                ! H_ij = <psi_i^* | H | psi_j> 
                hamiltonian(idx1, idx2) = hamiltonian(idx1, idx2) +  weights(pos) * &
                                        & dconjg(psi(idx1, pos)) * (- psi_d2dx2(idx2, pos) / 2.0_xr  &
                                        & + v_pot(x(pos)) * psi(idx2, pos))
            end do
        end do
    end do

    allocate (s(0 : n_basis))
    ! The S matrix
    do idx1 = 0, n_basis
        s(idx1) = complex(0.0_xr, 0.0_xr)
        do pos = 1, n_grid_pts
            ! s_i = <psi_i^* | L | j_l> L = H - KE
            ! with l = 0 so the centifugal term in H can be ignored in the current implementation.
            s(idx1) = s(idx1) + weights(pos) * (- psi(idx1, pos) * j(pos) * energy &
                        & + psi(idx1, pos) * (- j_d2dx2(pos) / 2.0_xr  + v_pot(x(pos)) * j(pos)))                    
        end do
    end do

    allocate (M(0 : n_basis, 0 : n_basis))
    allocate (M_copy(0 : n_basis, 0 : n_basis))
    ! The M matrix
    counter = 0
    do idx1 = 0, n_basis
        do idx2 = 0, n_basis
            M(idx1, idx2) =  complex(0.0_xr, 0.0_xr)
            do pos = 1, n_grid_pts
                ! M_ij = <psi_i^* | E - H | psi_j> 
                M(idx1, idx2) = M(idx1, idx2) +  weights(pos) * &
                                        & (- psi(idx1, pos) * energy * psi(idx2, pos)&
                                       & + psi(idx1, pos) * (- psi_d2dx2(idx2, pos) / 2.0_xr + v_pot(x(pos)) * psi(idx2, pos)))    
            end do
        end do
    end do
    M_copy(:, :) = M(:, :)

    allocate (ipiv(0 : n_basis))
    allocate (c(0 : n_basis))
    allocate (work(0 : n_basis))
    allocate (swork(0 : (n_basis+1)*(n_basis+1)))
    c(:) = s(:) 

    call zcgesv(n_basis+1, 1, M_copy, n_basis+1, ipiv, s, n_basis+1, c, n_basis+1, work, swork, rwork, iter, info)
    ! LAPACK routine for solving Ax=b with complex symmetric matrices.
    ! modifies the rhs of Ax=b in-place.
 
    lambda_result = - dot_product(conjg(s), c)
    ! Calculate the phase shift as what the Kohn's method trys to predict

    do pos = 1, n_grid_pts
        lambda_result = lambda_result + weights(pos) *  (j(pos) ** 2.0_xr) *  v_pot(x(pos)) 
    end do

    lambda_result =  - lambda_result * 2.0_xr / k
    print *, "lambda = "
    print *, lambda_result
    print *, "tan(lambda) = "
    print *, imag(lambda_result) / real(lambda_result)
    ! Calculate the tangent of the phase shift as what the Kohn's method trys to predict

    contains
        function v_pot_sr(r) result(pot)
            ! the short range potential
            integer, parameter :: xr = kind(0.d0)
            real(kind = xr) :: r, pot
            pot = - dexp(-r)
        end function v_pot_sr

        !function v_pot_lr(r) result(pot)
        !    ! the short range potential
        !    integer, parameter :: xr = kind(0.d0)
        !    real(kind = xr) :: r, pot
        !    pot = ((1 - dexp(-r)) ** 2.0_xr) / (r ** 3.0_xr)
        !end function v_pot_lr

        function slater(r, n, lmda) result(f)
            integer, parameter :: xr = kind(0.d0)
            real(kind = xr) :: r, lmda, f
            integer :: n
            f = (r ** n) * dexp(- lmda * r)
        end function slater

        function slater_d2dx2(r, n, lmda) result(f)
            integer, parameter :: xr = kind(0.d0)
            real(kind = xr) :: r, lmda, f
            integer :: n
            f = (n ** 2  + (lmda ** 2) * (r ** 2) - 2 * lmda * n * r - n) * r ** (n - 2) * dexp(- lmda * r) 
        end function slater_d2dx2
        
        function g(x) result(f)
            integer, parameter :: xr = kind(0.d0)
            real(kind = xr), intent(in) :: x
            real(kind = xr) :: f
            real(kind = xr) :: alpha = 1.0_xr
            f = 1.0_xr - dexp(- alpha * x)
        end function g

        function g_ddx(x) result(f)
            integer, parameter :: xr = kind(0.d0)
            real(kind = xr), intent(in) :: x
            real(kind = xr) :: f
            real(kind = xr) :: alpha = 1.0_xr
            f = alpha * dexp(- alpha * x)
        end function g_ddx

        function g_d2dx2(x) result(f)
            integer, parameter :: xr = kind(0.d0)
            real(kind = xr), intent(in) :: x
            real(kind = xr) :: f
            real(kind = xr) :: alpha = 1.0_xr
            f = - (alpha ** 2.0_xr) * dexp(- alpha * x)
        end function g_d2dx2
end program
