#include <iostream>

#include "lapack_c/mangling.h"
#include "lapack_c/rot.h"

// Macro to define both the Fortran interface and its C wrapper
#define INSTANTIATE_C_ROT(T, TC, name_lower, name_upper)                  \
    void FORTRAN_TO_C_NAME(name_lower, name_upper)(                       \
        const lapack_idx* n, T* x, const lapack_idx* incx, T* y,          \
        const lapack_idx* incy, const TC* c, const T* s);                 \
    void lapack_c_##name_lower(lapack_idx n, T* x, lapack_idx incx, T* y, \
                               lapack_idx incy, TC c, T s)                \
    {                                                                     \
        FORTRAN_TO_C_NAME(name_lower, name_upper)                         \
        (&n, x, &incx, y, &incy, &c, &s);                                 \
    }

// Actual definitions
extern "C" {

INSTANTIATE_C_ROT(float, float, srot, SROT)
INSTANTIATE_C_ROT(double, double, drot, DROT)
INSTANTIATE_C_ROT(lapack_float_complex, float, crot, CROT)
INSTANTIATE_C_ROT(lapack_double_complex, double, zrot, ZROT)

}  // extern "C"

#undef INSTANTIATE_C_ROT