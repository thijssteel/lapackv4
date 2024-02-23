program test

    integer, parameter :: n = 3
    integer, parameter :: m = 2

    integer :: i, j

    real :: a(m, n), x(n), y(m)

    do i = 1, m
        do j = 1, n
            a(i, j) = i + j
        end do
    end do

    do i = 1, n
        x(i) = i
    end do

    ! Calculate y = a * x using loops
    do i = 1, m
        y(i) = 0.0
        do j = 1, n
            y(i) = y(i) + a(i, j) * x(j)
        end do
    end do

    call sgemv( 'N', m, n, 1.0, a, m, x, 1, -1.0, y, 1 )

    print *, y


end program