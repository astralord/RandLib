#ifndef ECONOMETRICSENGINEER_H
#define ECONOMETRICSENGINEER_H

#include <cmath>

#include "../math/RandMath.h"
#include "../math/Matrix.h"
#include "randlib_global.h"

/**
 * @brief The EconometricsEngineer class
 */
class RANDLIBSHARED_EXPORT EconometricsEngineer
{
public:
    EconometricsEngineer();

    /**
     * @brief OLS
     * Ordinary least squares
     * @param X
     * @param Y
     * @param B
     * @return (X'X)^{-1} (X'Y)
     */
    template <size_t n, size_t m>
    static bool OLS(const Matrix<n, m> &X, const Vector<n> &Y, Vector<m> &B);
};

#endif // ECONOMETRICSENGINEER_H
