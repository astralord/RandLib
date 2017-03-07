#ifndef BETAMATH
#define BETAMATH

#include "RandMath.h"

/// Beta-related functions

namespace RandMath
{
/**
 * @brief logBeta
 * Calculate logarithm of beta function
 * @param a
 * @param b
 * @return log(B(a, b)) = log(Γ(a)) + log(Γ(b)) - log(Γ(a + b))
 */
double logBeta(double a, double b);

/**
 * @brief beta
 * Calculate beta function
 * @param a
 * @param b
 * @return B(a, b) = Γ(a) * Γ(b) / Γ(a + b)
 */
double beta(double a, double b);

/**
 * @brief ibeta
 * Fast calculation of regularized beta function, using precalculated values
 * @param x
 * @param a
 * @param b
 * @param logBetaFun log(B(a, b))
 * @param logX log(x)
 * @param log1mX log(1-x)
 * @return I(x, a, b) = B(x, a, b) / B(a, b)
 */
double ibeta(double x, double a, double b, double logBetaFun, double logX, double log1mX);

/**
 * @brief ibeta
 * Calculate regularized beta function
 * @param x
 * @param a
 * @param b
 * @return I(x, a, b) = B(x, a, b) / B(a, b)
 */
double ibeta(double x, double a, double b);

}

#endif // BETAMATH

