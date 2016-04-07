#include "EconometricsEngineer.h"

EconometricsEngineer::EconometricsEngineer()
{

}

bool EconometricsEngineer::OLS(const Matrix &X, const Vector &Y, Vector &B)
{
    /// Sanity checks
    if (B.width() != X.height() || X.height() != Y.height())
        return false;

    SquareMatrix A(X.height());
    X.getSquare(A);
    SquareMatrix C(A.size());
    A.getInverse(C);

    // AVOID COPYING JUST FOR TRANSPOSITION
    Matrix Z = X;
    Z.transpose();

    Vector W = Z * Y;
    B = C * W;
    return true;
}

