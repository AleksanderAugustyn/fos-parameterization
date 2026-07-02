!> Breakdown-boundary tests (umbrella design §5): the validity boundary the
!! module reports must coincide with the mathematical breakdown (interior
!! f(u) <= 0), and everywhere the module reports valid, volume conservation
!! and the R(theta) round-trip must hold to quadrature precision.
program fos_param_breakdown_test

    use precision_utilities_mod, only: ik, rk
    use mathematical_and_physical_constants_mod, only: PI_C
    use fos_parameterization_mod, only: compute_fos_radius_grid_s, &
            compute_fos_f_and_derivatives_s, compute_rho_z_grid_s, &
            validate_rho_grid_s, rho_z_grid_t, &
            FOS_VALID, FOS_ERROR_RHO_NEGATIVE, FOS_ERROR_BEAK_SINGULARITY
    use fos_test_reference_mod, only: init_quadrature_s, compute_reference_surface_f, &
            evaluate_shape_quality_s, find_neck_in_radii_s, gl_ref_x, gl_ref_w, N_GL_REF
    use test_utils_mod, only: assert_true, assert_int_eq, test_summary

    implicit none

    integer(kind = ik), parameter :: N_GRID_SMALL = 101_ik
    integer(kind = ik), parameter :: N_RHO_INTERNAL = 1001_ik
    real(kind = rk), parameter :: TOL_VOLUME = 1.0e-11_rk
    real(kind = rk), parameter :: TOL_SURFACE_REL = 1.0e-3_rk
    real(kind = rk), parameter :: TOL_ROUND_TRIP = 1.0e-9_rk
    ! Test-local mirrors of WMMM radius_grid_mod physics policy (Tier-A filter)
    real(kind = rk), parameter :: NECK_MIN_DEPTH_MIRROR = 0.25_rk
    real(kind = rk), parameter :: NECK_MIN_ELONGATION_MIRROR = 1.2_rk
    ! Symmetric-family analytic boundaries: f(0) = 1 - 4 a4 / 3
    real(kind = rk), parameter :: A4_BEAK_BOUNDARY = 0.999_rk * 0.75_rk   ! f(0) = F_MIN_THRESHOLD
    real(kind = rk), parameter :: A4_RHO_BOUNDARY = 0.75_rk               ! f(0) = 0

    call init_quadrature_s()
    call test_symmetric_family_code_boundaries()
    call test_beak_detector_two_sided()
    call test_error_fires_before_volume_breaks()
    call test_volume_conserved_everywhere_valid()
    call test_summary()

contains

    !> Exact code boundaries on the symmetric (c=2, a4) family. Check order in
    !! validate_rho_grid_s is rho-check then beak, both before star-convexity,
    !! so the code at each analytic boundary is deterministic.
    subroutine test_symmetric_family_code_boundaries()
        real(kind = rk) :: params(7), radii(N_GRID_SMALL), z_shift
        logical :: is_valid
        character(len = 256) :: message
        integer(kind = ik) :: code

        write(*, '(A)') '=== Symmetric-family code boundaries ==='

        ! Moderate neck, far from breakdown: accepted.
        params = 0.0_rk
        params(1) = 2.0_rk
        params(3) = 0.30_rk
        call compute_fos_radius_grid_s(params, N_GRID_SMALL, radii, z_shift, is_valid, &
                message, error_code = code)
        call assert_int_eq(code, FOS_VALID, 'boundary: c=2 a4=0.30 valid')

        ! Just past the beak boundary (f(0) ~ 1e-6 > 0, f_min < 1e-3): beak code.
        params(3) = A4_BEAK_BOUNDARY * (1.0_rk + 1.0e-3_rk)
        call compute_fos_radius_grid_s(params, N_GRID_SMALL, radii, z_shift, is_valid, &
                message, error_code = code)
        call assert_int_eq(code, FOS_ERROR_BEAK_SINGULARITY, &
                'boundary: just past beak threshold -> FOS_ERROR_BEAK_SINGULARITY')

        ! Past f(0) = 0 (grid point at u=0 has rho = 0): rho-negative code.
        params(3) = A4_RHO_BOUNDARY * (1.0_rk + 1.0e-3_rk)
        call compute_fos_radius_grid_s(params, N_GRID_SMALL, radii, z_shift, is_valid, &
                message, error_code = code)
        call assert_int_eq(code, FOS_ERROR_RHO_NEGATIVE, &
                'boundary: f(0) < 0 -> FOS_ERROR_RHO_NEGATIVE')
    end subroutine test_symmetric_family_code_boundaries

    !> Two-sided beak boundary at the layer that owns the detector:
    !! validate_rho_grid_s has no star-convexity gate, so both sides of the
    !! analytic boundary a4* = 0.999 * 0.75 are exactly testable here —
    !! impossible through compute_fos_radius_grid_s, where the star-convexity
    !! margin rejects these necked shapes regardless.
    subroutine test_beak_detector_two_sided()
        real(kind = rk) :: params(7)
        type(rho_z_grid_t) :: grid
        logical :: is_valid
        character(len = 256) :: message
        integer(kind = ik) :: code

        write(*, '(A)') '=== Beak detector: two-sided boundary (unit level) ==='

        params = 0.0_rk
        params(1) = 2.0_rk

        ! Just below the boundary: f(0) = 1 - (4/3) a4 ~ 2e-3 > F_MIN_THRESHOLD.
        params(3) = A4_BEAK_BOUNDARY * (1.0_rk - 1.0e-3_rk)
        call compute_rho_z_grid_s(params, N_RHO_INTERNAL, grid, code, message)
        call assert_int_eq(code, FOS_VALID, 'unit: grid built below beak boundary')
        call validate_rho_grid_s(grid, params, is_valid, code, message)
        call assert_int_eq(code, FOS_VALID, 'unit: 0.999 a4* passes beak detection')

        ! Just above: f(0) ~ 1e-6, in (0, F_MIN_THRESHOLD) -> beak code.
        params(3) = A4_BEAK_BOUNDARY * (1.0_rk + 1.0e-3_rk)
        call compute_rho_z_grid_s(params, N_RHO_INTERNAL, grid, code, message)
        call assert_int_eq(code, FOS_VALID, 'unit: grid built above beak boundary')
        call validate_rho_grid_s(grid, params, is_valid, code, message)
        call assert_int_eq(code, FOS_ERROR_BEAK_SINGULARITY, &
                'unit: 1.001 a4* -> FOS_ERROR_BEAK_SINGULARITY')

        ! Past f(0) = 0: interior grid point at u = 0 (N_RHO_INTERNAL odd) has
        ! rho = 0 -> the rho check fires before the beak check.
        params(3) = A4_RHO_BOUNDARY * (1.0_rk + 1.0e-3_rk)
        call compute_rho_z_grid_s(params, N_RHO_INTERNAL, grid, code, message)
        call assert_int_eq(code, FOS_VALID, 'unit: grid built past f(0) = 0')
        call validate_rho_grid_s(grid, params, is_valid, code, message)
        call assert_int_eq(code, FOS_ERROR_RHO_NEGATIVE, &
                'unit: f(0) < 0 -> FOS_ERROR_RHO_NEGATIVE')
    end subroutine test_beak_detector_two_sided

    !> Scan a4 across the breakdown at c=2: wherever the geometric cylindrical
    !! volume pi int max(f,0) du deviates from 4 pi / 3, the module must have
    !! rejected the shape (error fires no later than conservation breaks).
    subroutine test_error_fires_before_volume_breaks()
        real(kind = rk) :: params(7), radii(N_GRID_SMALL), z_shift
        real(kind = rk) :: a4, f_val, v_geom, v_err
        logical :: is_valid
        character(len = 256) :: message
        character(len = 96) :: label
        integer(kind = ik) :: code, i, k, n_broken, n_rejected

        write(*, '(A)') '=== Direction test: error no later than volume break ==='

        n_broken = 0_ik
        n_rejected = 0_ik

        do k = 0_ik, 30_ik
            a4 = 0.70_rk + 0.002_rk * real(k, rk)
            params = 0.0_rk
            params(1) = 2.0_rk
            params(3) = a4

            ! Geometric volume of the actual body in reduced units:
            ! V_geom = pi int rho^2 dz with rho^2 = max(f,0)/c and dz = c du,
            ! so the c factors cancel: V_geom = pi int max(f(u), 0) du.
            v_geom = 0.0_rk
            do i = 1_ik, N_GL_REF
                call compute_fos_f_and_derivatives_s(params, gl_ref_x(i), f_val)
                v_geom = v_geom + gl_ref_w(i) * max(f_val, 0.0_rk)
            end do
            v_geom = PI_C * v_geom
            v_err = abs(v_geom - 4.0_rk * PI_C / 3.0_rk)

            call compute_fos_radius_grid_s(params, N_GRID_SMALL, radii, z_shift, is_valid, &
                    message, n_rho_grid = N_RHO_INTERNAL, error_code = code)

            if (v_err > TOL_VOLUME) then
                n_broken = n_broken + 1_ik
                write(label, '(A,F6.4,A,ES10.3)') 'direction: a4=', a4, ' vol broken by ', v_err
                call assert_true(.not. is_valid, trim(label) // ' -> must be rejected')
            end if
            if (.not. is_valid) n_rejected = n_rejected + 1_ik
        end do

        ! The scan must actually cross the boundary for the test to mean anything.
        call assert_true(n_broken > 0_ik, 'direction: scan reaches volume-broken territory')
        call assert_true(n_rejected > 0_ik, 'direction: scan contains rejected points')
    end subroutine test_error_fires_before_volume_breaks

    !> Coarse (c, a3, a4) scan: every shape the module reports valid must
    !! conserve volume, match the cylindrical reference surface, and round-trip
    !! R(theta) <-> rho(z) to quadrature precision. Pronounced-neck shapes at
    !! low elongation are excluded from accuracy stats (mirrors WMMM Tier A's
    !! physics filter; such shapes are rejected downstream by physics policy).
    subroutine test_volume_conserved_everywhere_valid()
        real(kind = rk), parameter :: C_VALUES(6) = &
                [0.8_rk, 1.0_rk, 1.5_rk, 2.0_rk, 2.5_rk, 3.0_rk]
        real(kind = rk), parameter :: A3_VALUES(4) = [0.0_rk, 0.15_rk, 0.30_rk, 0.45_rk]

        real(kind = rk) :: params(7), radii(N_GRID_SMALL), z_shift
        real(kind = rk) :: a4, s_ref, dv_rel, ds_rel, rt_max
        real(kind = rk) :: neck_radius, neck_depth, elongation
        logical :: is_valid, has_neck
        character(len = 256) :: message
        character(len = 96) :: label
        integer(kind = ik) :: code, ic, i3, k, n_checked, n_rejected, n_filtered

        write(*, '(A)') '=== Volume/surface/round-trip conserved everywhere valid ==='

        n_checked = 0_ik
        n_rejected = 0_ik
        n_filtered = 0_ik

        do ic = 1_ik, size(C_VALUES, kind = ik)
            do i3 = 1_ik, size(A3_VALUES, kind = ik)
                do k = 0_ik, 9_ik
                    a4 = -0.2_rk + 0.1_rk * real(k, rk)
                    params = 0.0_rk
                    params(1) = C_VALUES(ic)
                    params(2) = A3_VALUES(i3)
                    params(3) = a4

                    call compute_fos_radius_grid_s(params, N_GRID_SMALL, radii, z_shift, &
                            is_valid, message, n_rho_grid = N_RHO_INTERNAL, error_code = code)
                    if (.not. is_valid) then
                        n_rejected = n_rejected + 1_ik
                        cycle
                    end if

                    ! Physics filter mirroring WMMM Tier A (see module docstring).
                    call find_neck_in_radii_s(radii, has_neck, neck_radius, neck_depth)
                    elongation = 0.5_rk * (radii(1) + radii(N_GRID_SMALL))
                    if (has_neck .and. neck_depth > NECK_MIN_DEPTH_MIRROR &
                            .and. elongation < NECK_MIN_ELONGATION_MIRROR) then
                        n_filtered = n_filtered + 1_ik
                        cycle
                    end if
                    n_checked = n_checked + 1_ik

                    s_ref = compute_reference_surface_f(params)
                    call evaluate_shape_quality_s(params, z_shift, s_ref, dv_rel, ds_rel, rt_max)

                    write(label, '(A,F4.2,A,F4.2,A,F5.2)') 'scan c=', C_VALUES(ic), &
                            ' a3=', A3_VALUES(i3), ' a4=', a4
                    call assert_true(abs(dv_rel) < TOL_VOLUME, trim(label) // ': volume conserved')
                    call assert_true(abs(ds_rel) < TOL_SURFACE_REL, trim(label) // ': surface ok')
                    call assert_true(rt_max < TOL_ROUND_TRIP, trim(label) // ': round trip ok')
                end do
            end do
        end do

        write(*, '(A,I0,A,I0,A,I0)') 'Scan: checked ', n_checked, &
                ', rejected ', n_rejected, ', physics-filtered ', n_filtered
        call assert_true(n_checked >= 50_ik, 'scan: enough valid shapes checked')
        call assert_true(n_rejected >= 1_ik, 'scan: crosses the validity boundary')
    end subroutine test_volume_conserved_everywhere_valid

end program fos_param_breakdown_test
