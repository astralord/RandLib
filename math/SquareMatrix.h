#ifndef SQUAREMATRIX_H
#define SQUAREMATRIX_H

#include "Matrix.h"

template <size_t n>
class RANDLIBSHARED_EXPORT SquareMatrix : public Matrix<n, n>
{
public:
    SquareMatrix(const double initialValue = 0.0) : Matrix<n, n>(initialValue) {}
    SquareMatrix(const SquareMatrix<n> & other) : Matrix<n, n> (other) {}
    ~SquareMatrix() {}

    void toIdentity();
    /// Gauss-Jordan method (not fast and non-optimized)
    bool getInverse(SquareMatrix<n> &invertedMatrix) const;
};

#endif // SQUAREMATRIX_H
