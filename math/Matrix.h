#ifndef MATRIX_H
#define MATRIX_H

#include "randlib_global.h"
#include <vector>

class SquareMatrix;

/**
 * @brief The Matrix class
 *
 */
class RANDLIBSHARED_EXPORT Matrix {
protected:
    size_t n, m;
    std::vector<double> data;
public:
    Matrix(const size_t height, const size_t width,
           const double initial_value = 0.0);
    Matrix(const Matrix & other);
    Matrix & operator=(const Matrix & other);
    ~Matrix() {}

    inline size_t height() const { return n; }
    inline size_t width() const { return m; }

    double &operator()(const size_t i, const size_t j);
    double operator()(const size_t i, const size_t j) const;

    /// addition
    friend const Matrix operator+(const Matrix& left, const Matrix& right) {
        if (left.height() != right.height() || left.width() != right.width())
            return left; // we should throw exception here
        Matrix sum(left.height(), left.width());
        for (size_t i = 0; i != left.data.size(); ++i)
            sum.data[i] = left.data[i] + right.data[i];
        return sum;
    }

    Matrix &operator+=(const Matrix& right);

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

    bool getSquare(SquareMatrix & squaredMatrix);
    void clear();
};


/**
 * @brief The SquareMatrix class
 */
class RANDLIBSHARED_EXPORT SquareMatrix : public Matrix{

public:
    SquareMatrix(const size_t height, const double initial_value = 0.0);
    ~SquareMatrix() {}

    inline size_t size() { return n; }

    void toIdentity();

    /// Gauss-Jordan method (not fast and non-optimized)
    bool invert(SquareMatrix & invertedMatrix);
};

#endif // MATRIX_H
