module kinds
    use, intrinsic :: iso_fortran_env, only: int32

    ! Define the kind of the indices
    ! Make sure this is consistent with the kind of the indices in the
    ! C code (lapack_idx)
    integer, parameter :: idx_kind = int32

end module kinds