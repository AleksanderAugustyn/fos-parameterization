!> Shared reference machinery for the fos-parameterization test suites.
!!
!! Validates the FoS -> R(theta) conversion against references computed from
!! the cylindrical representation (independent of the conversion under test):
!!   - Volume: V = pi * int rho^2 dz = 4 pi / 3 exactly (a2 volume constraint)
!!   - Surface: S = 2 pi c * int sqrt(f(u)/c + f'(u)^2/(4 c^4)) du
!! Ported from WMMM tests/fos_geometry_validation_test.f08 (module-level parts).
module fos_test_reference_mod

    use precision_utilities_mod, only: ik, rk
    use mathematical_and_physical_constants_mod, only: PI_C
    use mathematical_utilities_mod, only: compute_gauss_legendre_quadrature_s
    use fos_parameterization_mod, only: compute_radius_at_theta_s, compute_rho_at_z_s, &
            compute_fos_f_and_derivatives_s

    implicit none

    private

    public :: init_quadrature_s
    public :: compute_reference_surface_f
    public :: evaluate_shape_quality_s
    public :: spheroid_surface_area_f
    public :: find_neck_in_radii_s
    public :: N_GL_THETA, N_GL_REF, V_SPHERE
    public :: gl_ref_x, gl_ref_w

    ! Quadrature sizes
    integer(kind = ik), parameter :: N_GL_THETA = 4096_ik  !! R(theta)-side V/S integrals
    integer(kind = ik), parameter :: N_GL_REF = 512_ik    !! cylindrical S_ref and u-integrals

    real(kind = rk), parameter :: V_SPHERE = 4.0_rk * PI_C / 3.0_rk

    ! Gauss-Legendre nodes/weights on [-1, 1], filled by init_quadrature_s
    real(kind = rk) :: gl_theta_x(N_GL_THETA), gl_theta_w(N_GL_THETA)
    real(kind = rk) :: gl_ref_x(N_GL_REF), gl_ref_w(N_GL_REF)

contains

    subroutine init_quadrature_s()
        call compute_gauss_legendre_quadrature_s(N_GL_THETA, gl_theta_x, gl_theta_w)
        call compute_gauss_legendre_quadrature_s(N_GL_REF, gl_ref_x, gl_ref_w)
    end subroutine init_quadrature_s

    !> Reference surface area from the cylindrical FoS representation:
    !! S = 2 pi c * int_{-1}^{1} sqrt( f(u)/c + f'(u)^2 / (4 c^4) ) du
    !! Independent of the R(theta) conversion under test. The integrand stays
    !! bounded at the poles because rho * drho/dz = f'/(2 c^2).
    function compute_reference_surface_f(params) result(s_ref)
        real(kind = rk), intent(in) :: params(:)
        real(kind = rk) :: s_ref

        integer(kind = ik) :: i
        real(kind = rk) :: c, f_val, fp_val, arg

        c = params(1)
        s_ref = 0.0_rk

        do i = 1_ik, N_GL_REF
            call compute_fos_f_and_derivatives_s(params, gl_ref_x(i), f_val, fp_val)
            arg = max(f_val, 0.0_rk) / c + fp_val**2 / (4.0_rk * c**4)
            s_ref = s_ref + gl_ref_w(i) * sqrt(arg)
        end do

        s_ref = 2.0_rk * PI_C * c * s_ref
    end function compute_reference_surface_f

    !> Probes the R(theta) conversion at Gauss-Legendre nodes in x = cos(theta).
    !!
    !! Volume:  V = (2 pi / 3) int R^3 dx
    !! Surface: S = 2 pi int R sqrt(R^2 + (dR/dtheta)^2) dx
    !! dR/dtheta comes analytically from the implicit surface relation
    !! F(R, theta) = R sin(theta) - rho(R cos(theta)) = 0:
    !!   dR/dtheta = (R cos(theta) + rho' R sin(theta)) / (rho' cos(theta) - sin(theta))
    !! using the exact rho' from compute_rho_at_z_s - no finite differences.
    !! Round-trip residual |R sin(theta) - rho(R cos(theta))| checks each node directly.
    subroutine evaluate_shape_quality_s(params, z_shift, s_ref, dv_rel, ds_rel, rt_max)
        real(kind = rk), intent(in) :: params(:)
        real(kind = rk), intent(in) :: z_shift
        real(kind = rk), intent(in) :: s_ref
        real(kind = rk), intent(out) :: dv_rel
        real(kind = rk), intent(out) :: ds_rel
        real(kind = rk), intent(out) :: rt_max

        integer(kind = ik) :: i
        real(kind = rk) :: x, sin_theta, theta, r, z, rho, drho_dz
        real(kind = rk) :: denom, dr_dtheta, vol, surf

        vol = 0.0_rk
        surf = 0.0_rk
        rt_max = 0.0_rk

        do i = 1_ik, N_GL_THETA
            x = gl_theta_x(i)
            sin_theta = sqrt(max(1.0_rk - x**2, 0.0_rk))
            theta = acos(x)

            call compute_radius_at_theta_s(params, theta, z_shift, r)
            vol = vol + gl_theta_w(i) * r**3

            z = r * x
            call compute_rho_at_z_s(params, z, z_shift, rho, drho_dz)
            rt_max = max(rt_max, abs(r * sin_theta - rho))

            denom = drho_dz * x - sin_theta
            if (abs(denom) > 1.0e-12_rk) then
                dr_dtheta = (r * x + drho_dz * r * sin_theta) / denom
            else
                ! Ray nearly tangent to the surface (star-convexity margin).
                ! Any resulting inaccuracy shows up in the surface comparison.
                dr_dtheta = 0.0_rk
            end if
            surf = surf + gl_theta_w(i) * r * sqrt(r**2 + dr_dtheta**2)
        end do

        vol = vol * 2.0_rk * PI_C / 3.0_rk
        surf = surf * 2.0_rk * PI_C
        dv_rel = vol / V_SPHERE - 1.0_rk
        ds_rel = surf / s_ref - 1.0_rk
    end subroutine evaluate_shape_quality_s

    !> Closed-form spheroid surface area for FoS with all a_k = 0:
    !! f = 1 - u^2, so semi-axes are polar cp = c and equatorial ae = 1/sqrt(c).
    function spheroid_surface_area_f(c) result(s)
        real(kind = rk), intent(in) :: c
        real(kind = rk) :: s

        real(kind = rk) :: ae, cp, e

        cp = c
        ae = 1.0_rk / sqrt(c)

        if (abs(c - 1.0_rk) < 1.0e-12_rk) then
            s = 4.0_rk * PI_C
        else if (c > 1.0_rk) then
            ! Prolate: e^2 = 1 - (ae/cp)^2 = 1 - 1/c^3
            e = sqrt(1.0_rk - 1.0_rk / c**3)
            s = 2.0_rk * PI_C * ae**2 * (1.0_rk + (cp / (ae * e)) * asin(e))
        else
            ! Oblate: e^2 = 1 - (cp/ae)^2 = 1 - c^3
            e = sqrt(1.0_rk - c**3)
            s = 2.0_rk * PI_C * ae**2 * (1.0_rk + ((1.0_rk - e**2) / e) * atanh(e))
        end if
    end function spheroid_surface_area_f

    !> Neck scan over an R(theta) radii array (uniform theta grid on [0, pi]).
    !! Test-local clone of WMMM radius_grid_mod::find_neck_in_radii_s — used only
    !! to mirror WMMM's Tier-A physics filter (exclude pronounced-neck shapes at
    !! low elongation from conversion-accuracy statistics). Policy thresholds
    !! live at the call site, not here.
    pure subroutine find_neck_in_radii_s(radii, has_neck, neck_radius, neck_depth)
        real(kind = rk), intent(in) :: radii(:)
        logical, intent(out) :: has_neck
        real(kind = rk), intent(out) :: neck_radius
        real(kind = rk), intent(out) :: neck_depth

        integer(kind = ik) :: i, n, max1_idx, max2_idx, left_idx, right_idx
        real(kind = rk) :: theta_spacing, max1_rho, max2_rho
        real(kind = rk) :: rho(size(radii))

        has_neck = .false.
        neck_radius = 0.0_rk
        neck_depth = 0.0_rk
        n = size(radii, kind = ik)
        if (n < 5_ik) return

        theta_spacing = PI_C / real(n - 1_ik, rk)
        do i = 1_ik, n
            rho(i) = radii(i) * sin(real(i - 1_ik, rk) * theta_spacing)
        end do

        ! Two largest interior local maxima of rho(theta)
        max1_idx = 0_ik
        max2_idx = 0_ik
        max1_rho = -1.0_rk
        max2_rho = -1.0_rk
        do i = 2_ik, n - 1_ik
            if (rho(i) > rho(i - 1) .and. rho(i) > rho(i + 1)) then
                if (rho(i) > max1_rho) then
                    max2_rho = max1_rho
                    max2_idx = max1_idx
                    max1_rho = rho(i)
                    max1_idx = i
                else if (rho(i) > max2_rho) then
                    max2_rho = rho(i)
                    max2_idx = i
                end if
            end if
        end do

        if (max1_idx == 0_ik .or. max2_idx == 0_ik) return

        left_idx = min(max1_idx, max2_idx)
        right_idx = max(max1_idx, max2_idx)

        neck_radius = huge(1.0_rk)
        do i = left_idx, right_idx
            if (rho(i) < neck_radius) neck_radius = rho(i)
        end do
        neck_depth = min(max1_rho, max2_rho) - neck_radius
        has_neck = .true.

    end subroutine find_neck_in_radii_s

end module fos_test_reference_mod
