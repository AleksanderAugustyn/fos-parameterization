!> Status-code tests: FoS codes are >= 100, shared shape_core codes re-export,
!! and every code maps to a non-empty message.
program fos_param_status_test
    use precision_utilities_mod, only: ik
    use fos_parameterization_mod, only: FOS_VALID, FOS_ERROR_RHO_NEGATIVE, &
            FOS_ERROR_NOT_STAR_CONVEX, FOS_ERROR_INVALID_C, &
            FOS_ERROR_BEAK_SINGULARITY, FOS_ERROR_CONVERGENCE, &
            FOS_ERROR_BUFFER_MISMATCH, status_message
    use shape_core_mod, only: SHAPE_VALID, SHAPE_ERROR_TOO_MANY_PARAMS
    use test_utils_mod, only: assert_int_eq, assert_true, test_summary
    implicit none

    call assert_int_eq(int(FOS_VALID), int(SHAPE_VALID), 'FOS_VALID aliases SHAPE_VALID')
    call assert_int_eq(int(FOS_ERROR_RHO_NEGATIVE), 100, 'rho_negative = 100')
    call assert_int_eq(int(FOS_ERROR_NOT_STAR_CONVEX), 101, 'not_star_convex = 101')
    call assert_int_eq(int(FOS_ERROR_INVALID_C), 102, 'invalid_c = 102')
    call assert_int_eq(int(FOS_ERROR_BEAK_SINGULARITY), 103, 'beak = 103')
    call assert_int_eq(int(FOS_ERROR_CONVERGENCE), 104, 'convergence = 104')
    call assert_int_eq(int(FOS_ERROR_BUFFER_MISMATCH), 105, 'buffer = 105')
    call assert_true(len_trim(status_message(102_ik)) > 0, 'message for 102 nonempty')
    call assert_true(len_trim(status_message(SHAPE_ERROR_TOO_MANY_PARAMS)) > 0, 'shared code has message')
    call assert_true(len_trim(status_message(999_ik)) > 0, 'unknown code has fallback message')
    call test_summary()
end program fos_param_status_test
