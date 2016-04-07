#include "EconometricsEngineer.h"

EconometricsEngineer::EconometricsEngineer()
{

}

template <size_t n, size_t m>
bool EconometricsEngineer::OLS(const Matrix<n, m> &X, const Vector<n> &Y, Vector<m> &B)
{
    SquareMatrix<n> A;
    X.getSquare(A);
    SquareMatrix<n> C;
    A.getInverse(C);

    // AVOID COPYING JUST FOR TRANSPOSITION
    Matrix<n, m> Z = X;
    Z.transpose();

    Vector<m> W = Z * Y;
    B = C * W;
    return true;
}

