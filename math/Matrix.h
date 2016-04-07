#ifndef MATRIX_H
#define MATRIX_H

#include "randlib_global.h"
#include <vector>

template <size_t n>
class SquareMatrix;

/**
 * @brief The Matrix class
 *
 */
template <size_t n, size_t m>
class RANDLIBSHARED_EXPORT Matrix {
protected:
    std::vector<double> data;
public:
    Matrix(double initialValue);
    Matrix(const Matrix<n, m> & other);
    Matrix<n, m> & operator=(const Matrix<n, m> & other);
    ~Matrix() {}

    inline size_t height() const { return n; }
    inline size_t width() const { return m; }

    double &operator()(const size_t i, const size_t j) {
        return data[i * m + j];
    }

    double operator()(const size_t i, const size_t j) const {
        return data[i * m + j];
    }

    /// addition
    friend const Matrix operator+(const Matrix& left, const Matrix& right) {
        if (left.height() != right.height() || left.width() != right.width())
            return left; // we should throw exception here
        Matrix sum(left.height(), left.width());
        for (size_t i = 0; i != left.data.size(); ++i)
            sum.data[i] = left.data[i] + right.data[i];
        return sum;
    }

    Matrix<n, m> &operator+=(const Matrix<n, m> &right);

    /// multiplication
    friend const Matrix operator*(const Matrix& left, const Matrix& right) {
        if (left.width() != right.height())
            return left; // we should throw exception here
        Matrix product(left.height(), right.width());
        /// pretty straightforward
        for (size_t i = 0; i != product.height(); ++i) {
            for (size_t j = 0; j != product.width(); ++j) {
                for (size_t k = 0; k != left.width(); ++k)
                    product(i, j) += left(i, k) * right(k, j);
            }
        }
        return product;
    }

    bool getSquare(SquareMatrix<n> & squaredMatrix);
    void fill(double value);
    void clear();
};


/**
 * @brief The SquareMatrix class
 */
template <size_t n>
class RANDLIBSHARED_EXPORT SquareMatrix : public Matrix<n, n> {

public:
    SquareMatrix(const double initial_value = 0.0);
    ~SquareMatrix() {}

    inline size_t size() { return n; }

    void toIdentity();

    /// Gauss-Jordan method (not fast and non-optimized)
    bool invert(SquareMatrix & invertedMatrix);
};

#endif // MATRIX_H
