!> Minimal assertion helpers for the fos-parameterization test suites.
module test_utils_mod

    use precision_utilities_mod, only: ik, ikl, rk

    implicit none

    private

    public :: assert_true, assert_int_eq, assert_close, assert_abs_close, &
            assert_bits_eq, test_summary

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

    !> IEEE754 bit equality — the zero-tolerance assert of the contract's
    !! bitwise family. `==` on reals is banned by -Wcompare-reals, and a value
    !! compare would in any case call +0.0 and -0.0 equal; only the bit patterns
    !! answer the question the contract asks.
    !!
    !! `volatile` on the two locals is load-bearing: under Release's -ffast-math
    !! (-fno-signed-zeros) the compiler may elide a store of -0.0 over a location
    !! it knows holds +0.0, so the pattern under test would never reach memory.
    subroutine assert_bits_eq(a, b, label)
        real(kind = rk),    intent(in) :: a, b
        character(len = *), intent(in) :: label
        real(kind = rk), volatile :: av, bv
        av = a
        bv = b
        call assert_true(transfer(av, 0_ikl) == transfer(bv, 0_ikl), label)
    end subroutine assert_bits_eq

    subroutine test_summary()
        write(*, '(A,I0,A,I0,A)') 'Tests: ', n_pass, ' passed, ', n_fail, ' failed.'
        if (n_fail > 0_ik) error stop 1
    end subroutine test_summary

end module test_utils_mod
