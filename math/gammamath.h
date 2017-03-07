#ifndef GAMMAMATH
#define GAMMAMATH

#include "RandMath.h"

/// Gamma-related functions

namespace RandMath
{

/**
 * @brief digamma
 * @param x
 * @return digamma(x) = ln(Γ(x))' = d(Γ(x))/Γ(x)
 */
double digamma(double x);

/**
 * @brief logmdigamma
 * @param x
 * @return digamma(x) - log(x)
 */
double digammamLog(double x);

/**
 * @brief trigamma
 * @param x
 * @return trigamma(x) = digamma'(x)
 */
double trigamma(double x);

/**
 * @brief lpgamma
 * Calculate logarithm of lower regularised incomplete gamma function log(P(a, x))
 * @param a
 * @param x
 * @return log(P(a, x))
 */
double lpgamma(double a, double x);

 /**
 * @brief pgamma
 * Calculate lower regularised incomplete gamma function P(a, x)
 * @param a
 * @param x
 * @return P(a, x)
 */
double pgamma(double a, double x);

 /**
 * @brief lqgamma
 * Calculate logarithm of upper incomplete gamma function
 * @param a
 * @param x
 * @return log(Q(a, x))
 */
double lqgamma(double a, double x);

/**
 * @brief qgamma
 * Calculate upper incomplete gamma function
 * @param a
 * @param x
 * @return Q(a, x)
 */
double qgamma(double a, double x);

}

#endif // GAMMAMATH

