module elgamal
    implicit none
    integer(8) c1, c2, d

    interface operator(.nxor.)
        module procedure nxor
    end interface

    interface operator(.lsh.)
        module procedure l_shift
    end interface

    interface operator(.rsh.)
        module procedure r_shift
    end interface
contains
    integer(8) function XorShift(n)
        implicit none
        integer(8), intent(inout) :: n
        integer(8), save :: x = 123456789, y = 362436069, z = 521288629, w
        integer(8) t, i
        n = n ** (n .rsh. 7_8) ** (n .lsh. 7_8)
        w = n
        do i = 0, 18
            t = x .nxor. (x .lsh. 11_8)
            x = y; y = z; z = w
            w = (w .nxor. (w .rsh. 19_8)) .nxor. (t .nxor. (t .rsh. 8_8))
        end do
        XorShift = w
    end function XorShift

    pure integer(8) function modPow(a, k, n)
        implicit none
        integer(8), intent(in) :: a, k, n
        integer(8) i, va, t
        t = mod(a, n)
        if (a .eq. 0 .or. n .eq. 0) then
            modPow = 0
        else if (k .eq. 0) then
            modPow = mod(1, n)
        end if
        va = 1
        do i = 0, k-1
            va = va * t
            if (va >= n) va = mod(va, n)
        end do
        modPow = va
    end function modPow

    integer(8) function nxor(a, b) result(c)
        implicit none
        integer(8), intent(in) :: a, b
        c = xor(a, b)
    end function nxor

    integer(8) function l_shift(a, b) result(c)
        implicit none
        integer(8), intent(in) :: a, b
        c = lshift(a, b)
    end function l_shift

    integer(8) function r_shift(a, b) result(c)
        implicit none
        integer(8), intent(in) :: a, b
        c = rshift(a, b)
    end function r_shift
end module elgamal

subroutine elgamal_encrypt(m, p, g, y)
    use elgamal
    implicit none
    integer(8), intent(inout) :: m, p, g, y
    integer(8) r, tmp
    if (0 <= m .and. m < p) then
        tmp = time()
        do
            tmp = tmp + 1
            r = mod(XorShift(tmp), 1000000)
            if (r < 2 .and. r > (p-1)) cycle
            exit
        end do
        c1 = modPow(g, r, p)
        c2 = mod((m * modPow(y, r, p)), p)
    else
        stop
    end if
contains
end subroutine elgamal_encrypt

subroutine elgamal_decrypt(p, sk)
    use elgamal
    implicit none
    integer(8), intent(inout) :: p, sk
    d = mod((c2 * modPow(c1, p - 1 - sk, p)), p)
end subroutine elgamal_decrypt

program main
    use elgamal
    implicit none
    integer(8) p, q, g, x, y, tmp, m
    tmp = time()
    do
        tmp = tmp + 1
        q = mod(XorShift(tmp), 1000000)
        do
            if (.not. is_prime(q)) then
                q = q + 1
            else
                exit
            end if
        end do
        p = 2 * q + 1
        if (is_prime(p)) exit
    end do
    g = primitive_root(p, q)
    tmp = time()
    do
        tmp = tmp + 1
        x = mod(XorShift(tmp), 1000000)
        if (x < 2 .and. x > (p-1)) cycle
        exit
    end do
    y = modPow(g, x, p)

    print *, "pk:", p, g, y
    print *, "sk:", x

    m = 1919
    print *, "m:", m
    call elgamal_encrypt(m, p, g, y)
    print *, "c:", c1, c2

    call elgamal_decrypt(p, x)
    print *, "d:", d
contains
    logical function is_prime(n)
        implicit none
        integer(8), intent(inout) :: n
        integer(8) i
        if (mod(n, 2) .eq. 0 .or. n <= 2) then
            is_prime = .false.
            return
        end if
        do i = 3, int(sqrt(real(n))), 2
            if (mod(n, i) .eq. 0) then
                is_prime = .false.
                return
            end if
        end do
        is_prime = .true.
    end function is_prime

    integer(8) function primitive_root(p, q)
        use elgamal
        implicit none
        integer(8), intent(inout) :: p, q
        integer(8) g, tmp
        tmp = time()
        do
            tmp = tmp + 1
            g = mod(XorShift(tmp), 1000000)
            if (g < 3 .and. g > p) cycle
            if (modPow(g, 2_8, p) .eq. 1) cycle
            if (modPow(g, q, p) .eq. 1) cycle
            exit
        end do
        primitive_root = g
    end function primitive_root
end program main
