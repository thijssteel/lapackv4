program test

    integer, parameter :: n = 3

    integer :: i, j, l

    real :: a(n, n), b(n, n), c(n, n)

    do i = 1, n
        do j = 1, n
            a(i, j) = i + j
            b(i, j) = i - j
        end do
    end do

    do i = 1, n
        do j = 1, n
            C(i, j) = 0
            do l = 1, n
                C(i, j) = C(i, j) + a(i, l) * b(l, j)
            end do
        end do
    end do

    call sgemm( 'N', 'N', n, n, n, 1.0, a, n, b, n, -1.0, c, n )

    do i = 1, n
        do j = 1, n
            write(*, '(ES16.6)', advance='no') c(i, j)
        end do
        print *
    end do


end program