program test

    integer, parameter :: n = 3

    integer :: i, j

    real :: a(n, n), b(n, n), c(n, n)

    do i = 1, n
        do j = 1, n
            a(i, j) = i + j
            b(i, j) = i - j
            c(i, j) = 0.0
        end do
    end do

    call sgemm( 'N', 'N', n, n, n, 1.0, a, n, b, n, 0.0, c, n )

    do i = 1, n
        do j = 1, n
            print *, c(i, j)
        end do
    end do


end program