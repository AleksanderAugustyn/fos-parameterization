!> Minimal assertion helpers for the fos-parameterization test suites.
module test_utils_mod

    use precision_utilities_mod, only: ik, rk

    implicit none

    private

    public :: assert_true, assert_int_eq, assert_close, assert_abs_close, test_summary

    integer(kind = ik) :: n_pass = 0_ik
    integer(kind = ik) :: n_fail = 0_ik

contains

    subroutine assert_true(cond, label)
        logical,            intent(in) :: cond
        character(len = *), intent(in) :: label
        if (cond) then
            n_pass = n_pass + 1_ik
        else
            n_fail = n_fail + 1_ik
            write(*, '(A,A)') 'FAIL: ', label
        end if
    end subroutine assert_true

    subroutine assert_int_eq(got, want, label)
        integer(kind = ik), intent(in) :: got, want
        character(len = *), intent(in) :: label
        if (got == want) then
            n_pass = n_pass + 1_ik
        else
            n_fail = n_fail + 1_ik
            write(*, '(A,A,A,I0,A,I0)') 'FAIL: ', label, ' — got ', got, ', want ', want
        end if
    end subroutine assert_int_eq

    !> Mixed absolute/relative closeness: |got - want| <= tol * max(1, |want|).
    subroutine assert_close(got, want, tol, label)
        real(kind = rk),    intent(in) :: got, want, tol
        character(len = *), intent(in) :: label
        if (abs(got - want) <= tol * max(1.0_rk, abs(want))) then
            n_pass = n_pass + 1_ik
        else
            n_fail = n_fail + 1_ik
            write(*, '(A,A,A,ES23.16,A,ES23.16)') 'FAIL: ', label, ' — got ', got, ', want ', want
        end if
    end subroutine assert_close

    !> Pure absolute closeness — tolerances ported from WMMM's FoS suite are absolute.
    subroutine assert_abs_close(got, want, tol, label)
        real(kind = rk),    intent(in) :: got, want, tol
        character(len = *), intent(in) :: label
        if (abs(got - want) <= tol) then
            n_pass = n_pass + 1_ik
        else
            n_fail = n_fail + 1_ik
            write(*, '(A,A,A,ES23.16,A,ES23.16,A,ES10.3)') 'FAIL: ', label, &
                    ' — got ', got, ', want ', want, ', |diff| = ', abs(got - want)
        end if
    end subroutine assert_abs_close

    subroutine test_summary()
        write(*, '(A,I0,A,I0,A)') 'Tests: ', n_pass, ' passed, ', n_fail, ' failed.'
        if (n_fail > 0_ik) error stop 1
    end subroutine test_summary

end module test_utils_mod
