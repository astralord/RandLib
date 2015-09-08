#ifndef MATRIX_H
#define MATRIX_H

#include "randlib_global.h"
#include <vector>

/**
 * @brief The Matrix class
 *
 */
class RANDLIBSHARED_EXPORT Matrix
{
    size_t n, m;
    std::vector<double> data;
public:
    Matrix(const size_t height, const size_t width,
           const double initial_value = 0.0);
    Matrix(const Matrix & other);
    Matrix & operator=(const Matrix & other);
    
    inline size_t height() const { return n; }
    inline size_t width() const { return m; }

    double &operator()(const size_t i, const size_t j);
    double operator()(const size_t i, const size_t j) const;

    // TODO:
    // addition and multiplication
};

#endif // MATRIX_H
