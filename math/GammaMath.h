#ifndef GAMMAMATH
#define GAMMAMATH

#include "RandMath.h"

/// Gamma-related functions

namespace RandMath
{

/**
 * @fn factorial
 * Calculate n! using table values for small n <= 255
 * and built-in gamma function for large n > 255
 * @param n non-negative integer number
 * @return n!
 */
long double factorial(double n);

/**
 * @fn doubleFactorial
 * Calculate n!!
 * @param n non-negative integer number
 * @return n!!
 */
long double doubleFactorial(int n);

/**
 * @fn lfact
 * Calculate log(n!) using table values for small n <= 255
 * and built-in lgamma function for large n > 255
 * @param n non-negative integer number
 * @return log(n!)
 */
long double lfact(size_t n);

/**
 * @brief ldoublefact
 * @param n
 * @return log(n!!)
 */
long double ldfact(size_t n);

/**
 * @fn binom
 * Calculate binomial coefficient
 * @param n non-negative integer
 * @param k non-negative integer, k < n
 * @return C(n, k)
 */
long double binom(size_t n, size_t k);

/**
 * @fn digamma
 * @param x input parameter
 * @return digamma(x) = ln(Γ(x))' = d(Γ(x))/Γ(x)
 */
double digamma(double x);

/**
 * @fn digammamLog
 * @param x input parameter
 * @return digamma(x) - log(x)
 */
double digammamLog(double x);

/**
 * @fn trigamma
 * @param x input parameter
 * @return trigamma(x) = digamma'(x)
 */
double trigamma(double x);

/**
 * @fn lpgamma
 * Calculate logarithm of lower incomplete gamma function,
 * accelerated by using precalculated value of log(a) and log(Γ(a))
 * @param a non-negative parameter
 * @param x non-negative parameter
 * @param logA log(a)
 * @param lgammaA log(Γ(a))
 * @return log(P(a, x))
 */
double lpgamma(double a, double x, double logA, double lgammaA);

/**
 * @fn lpgamma
 * Calculate logarithm of lower incomplete gamma function,
 * accelerated by using precalculated value of log(x)
 * @param a non-negative parameter
 * @param x non-negative parameter
 * @param logX log(x)
 * @return log(P(a, x))
 */
double lpgamma(double a, double x, double logX);

/**
 * @fn lpgamma
 * Calculate logarithm of lower incomplete gamma function
 * @param a non-negative parameter
 * @param x non-negative parameter
 * @return log(P(a, x))
 */
double lpgamma(double a, double x);

/**
 * @fn pgamma
 * Calculate lower regularised incomplete gamma function,
 * accelerated by using precalculated value of log(a) and log(Γ(a))
 * @param a non-negative parameter
 * @param x non-negative parameter
 * @param logA log(a)
 * @param lgammaA log(Γ(a))
 * @return P(a, x)
 */
double pgamma(double a, double x, double logA, double lgammaA);

/**
 * @fn pgamma
 * Calculate lower regularised incomplete gamma function,
 * accelerated by using precalculated value of log(x)
 * @param a non-negative parameter
 * @param x non-negative parameter
 * @param logX log(x)
 * @return P(a, x)
 */
double pgamma(double a, double x, double logX);

 /**
 * @fn pgamma
 * Calculate lower regularised incomplete gamma function
 * @param a non-negative parameter
 * @param x non-negative parameter
 * @return P(a, x)
 */
double pgamma(double a, double x);

/**
 * @fn lqgamma
 * Calculate logarithm of upper incomplete gamma function,
 * accelerated by using precalculated values of log(a) and log(Γ(a))
 * @param a non-negative parameter
 * @param x non-negative parameter
 * @param logA log(a)
 * @param lgammaA log(Γ(a))
 * @return log(Q(a, x))
 */
double lqgamma(double a, double x, double logA, double lgammaA);

/**
 * @fn lqgamma
 * Calculate logarithm of upper incomplete gamma function,
 * accelerated by using precalculated value of log(x)
 * @param a non-negative parameter
 * @param x non-negative parameter
 * @param logX log(x)
 * @return log(Q(a, x))
 */
double lqgamma(double a, double x, double logX);

 /**
 * @fn lqgamma
 * Calculate logarithm of upper incomplete gamma function
 * @param a non-negative parameter
 * @param x non-negative parameter
 * @return log(Q(a, x))
 */
double lqgamma(double a, double x);

/**
 * @fn lqgamma
 * Calculate upper incomplete gamma function,
 * accelerated by using precalculated values of log(a) and log(Γ(a))
 * @param a non-negative parameter
 * @param x non-negative parameter
 * @param logA log(a)
 * @param lgammaA log(Γ(a))
 * @return Q(a, x)
 */
double qgamma(double a, double x, double logA, double lgammaA);

/**
 * @fn qgamma
 * Calculate upper incomplete gamma function,
 * accelerated by using precalculated value of log(x)
 * @param a non-negative parameter
 * @param x non-negative parameter
 * @param logX log(x)
 * @return Q(a, x)
 */
double qgamma(double a, double x, double logX);

/**
 * @fn qgamma
 * Calculate upper incomplete gamma function
 * @param a non-negative parameter
 * @param x non-negative parameter
 * @return Q(a, x)
 */
double qgamma(double a, double x);

}

#endif // GAMMAMATH

