module kinds
    use, intrinsic :: iso_fortran_env, only: int32
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none

    ! Define the kind of the indices
    ! Make sure this is consistent with the kind of the indices in the
    ! C code (lapack_idx)
    integer, parameter :: idx_kind = c_int

end module kinds