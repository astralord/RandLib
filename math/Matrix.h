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
    bool isTransposed;
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

    bool getSquare(SquareMatrix & squaredMatrix) const;
    void clear();
    void transpose();
};

/**
 * @brief The Vector class
 */
class RANDLIBSHARED_EXPORT Vector : public Matrix {
public:
    Vector(const size_t height, const double initial_value = 0.0);
    ~Vector() {}

    inline size_t size() const { return n; }

    /// addition
    friend const Vector operator+(const Vector& left, const Vector& right) {
        if (left.size() != right.size())
            return left; // we should throw exception here
        Vector sum(left.size());
        for (size_t i = 0; i != left.data.size(); ++i)
            sum.data[i] = left.data[i] + right.data[i];
        return sum;
    }

    Vector &operator+=(const Vector& right);

    /// multiplication
    friend const Vector operator*(const Matrix& left, const Vector& right) {
        if (left.width() != right.height())
            return right; // we should throw exception here
        Vector product(left.height());
        /// pretty straightforward
        for (size_t i = 0; i != product.size(); ++i) {
            for (size_t k = 0; k != left.width(); ++k)
                product(i, 0) += left(i, k) * right(k, 0); // make an operator() with one parameter!
        }
        return product;
    }

    friend const Vector operator*(const Vector& left, const Matrix& right) {
        if (left.width() != right.height())
            return left; // we should throw exception here
        Vector product(right.width());
        product.transpose();
        /// pretty straightforward
        for (size_t i = 0; i != product.size(); ++i) {
            for (size_t k = 0; k != left.width(); ++k)
                product(0, i) += left(0, k) * right(k, i); // make an operator() with one parameter!
        }
        return product;
    }
};


/**
 * @brief The SquareMatrix class
 */
class RANDLIBSHARED_EXPORT SquareMatrix : public Matrix {
public:
    SquareMatrix(const size_t height, const double initial_value = 0.0);
    ~SquareMatrix() {}

    inline size_t size() const { return n; }

    void toIdentity();

    /// Gauss-Jordan method (not fast and non-optimized)
    bool getInverse(SquareMatrix & invertedMatrix);
};

#endif // MATRIX_H
