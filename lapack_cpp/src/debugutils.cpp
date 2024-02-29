
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

#include "lapack_cpp/base.hpp"

namespace lapack_cpp {

template<typename T, Layout layout, typename idx_t>
std::string visualize_matrix(const Matrix<T, layout, idx_t>& A){

    const idx_t m = std::min<idx_t>(100, A.num_rows());
    const idx_t n = std::min<idx_t>(100, A.num_columns());


    std::stringstream stream;
    stream << "{ \"kind\":{ \"plotly\": true },\"data\":[{";

    // Add the header (matrix index)
    stream << "\"header\":{\"values\":[";
    for (idx_t j = 0; j < n; ++j) {
        stream << j;
        if (j + 1 < n) stream << ", ";
    }
    stream << "]},";

    // Add the matrix values
    stream << "\"cells\":{\"values\":[";
    for (idx_t j = 0; j < n; ++j) {
        stream << "[";
        for (idx_t i = 0; i < m; ++i) {
            stream << "\"" << std::setprecision(3) << A(i, j) << "\"";
            if (i + 1 < m) stream << ", ";
        }
        stream << "]";
        if (j + 1 < n) stream << ", ";
    }

    stream << "]},";
    stream << "\"type\": \"table\"}],\"layout\": {}}";

    return stream.str();
}

// Explicitly instantiate so they are available in gdb
std::string visualize_matrix_r( const Matrix<float, Layout::ColMajor, lapack_idx_t>& A){
    return visualize_matrix(A);
}
std::string visualize_matrix_d( const Matrix<double, Layout::ColMajor, lapack_idx_t>& A){
    return visualize_matrix(A);
}
std::string visualize_matrix_c( const Matrix<std::complex<float>, Layout::ColMajor, lapack_idx_t>& A){
    return visualize_matrix(A);
}
std::string visualize_matrix_z( const Matrix<std::complex<double>, Layout::ColMajor, lapack_idx_t>& A){
    return visualize_matrix(A);
}

std::string visualize_matrix_r_rm( const Matrix<float, Layout::RowMajor, lapack_idx_t>& A){
    return visualize_matrix(A);
}
std::string visualize_matrix_d_rm( const Matrix<double, Layout::RowMajor, lapack_idx_t>& A){
    return visualize_matrix(A);
}
std::string visualize_matrix_c_rm( const Matrix<std::complex<float>, Layout::RowMajor, lapack_idx_t>& A){
    return visualize_matrix(A);
}
std::string visualize_matrix_z_rm( const Matrix<std::complex<double>, Layout::RowMajor, lapack_idx_t>& A){
    return visualize_matrix(A);
}


} // namespace lapack_cpp