#ifndef LAPACK_C_MANGLING_H
#define LAPACK_C_MANGLING_H

// NOTE: this file is a modification from CBLAS and BLAS++

#include "lapack_c/defines.h"

// -----------------------------------------------------------------------------
// Fortran name mangling depends on compiler.
// Define FORTRAN_UPPER for uppercase,
// define FORTRAN_LOWER for lowercase (IBM xlf),
// else the default is lowercase with appended underscore
// (GNU gcc, Intel icc, PGI pgfortan, Cray ftn).
#ifndef FORTRAN_TO_C_NAME
    #if defined(FORTRAN_TO_C_UPPER)
        #define FORTRAN_TO_C_NAME( lower, UPPER ) UPPER
    #elif defined(FORTRAN_TO_C_LOWER)
        #define FORTRAN_TO_C_NAME( lower, UPPER ) lower
    #else
        #define FORTRAN_TO_C_NAME( lower, UPPER ) lower##_
    #endif
#endif

#endif        //  #ifndef BLAS_MANGLING_H