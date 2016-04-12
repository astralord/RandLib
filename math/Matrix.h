#ifndef MATRIX_H
#define MATRIX_H

#include "randlib_global.h"
#include <vector>

/**
 * @brief The Matrix class
 *
 * Allows to work with matrices (maximum size 20x20)
 */
template <size_t n, size_t m>
class RANDLIBSHARED_EXPORT Matrix {
protected:
    std::vector<double> data;
public:
    Matrix(const double initialValue = 0.0);
    Matrix(const Matrix<n, m> & other);
    Matrix<n, m> & operator=(const Matrix<n, m> & other);
    ~Matrix() {}

    inline size_t height() const { return n; }
    inline size_t width() const { return m; }

    double &operator()(const size_t i, const size_t j);
    double operator()(const size_t i, const size_t j) const;

    double &operator()(const size_t i);
    double operator()(const size_t i) const;

    Matrix<n, m> &operator+=(const Matrix<n, m> &right);

    /// multiplication
    template <size_t q>
    friend const Matrix<n, q> operator*(const Matrix<n, m>& left, const Matrix<m, q>& right) {
        Matrix<n, q> product;
        /// pretty straightforward
        for (size_t i = 0; i != n; ++i) {
            for (size_t j = 0; j != q; ++j) {
                for (size_t k = 0; k != m; ++k)
                    product(i, j) += left(i, k) * right(k, j);
            }
        }
        return product;
    }

    bool getSquare(Matrix<n, n> & squaredMatrix) const;
    void fill(double value);
    void clear();

    bool getTransposed(Matrix<m, n> & transposedMatrix) const;

    void toIdentity();
    /// Gauss-Jordan method (not fast and non-optimized)
    bool getInverse(Matrix<m, n> &invertedMatrix) const;
};

template <size_t n>
using Vector = Matrix<n, 1>;

template <size_t n>
using VectorT = Matrix<1, n>;

template <size_t n>
using SquareMatrix = Matrix<n, n>;

#endif // MATRIX_H
