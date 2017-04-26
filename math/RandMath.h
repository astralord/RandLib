#ifndef RANDMATH_H
#define RANDMATH_H

#include <cmath>
#include <functional>
#include <complex>
#include <vector>
#include <utility>
#include <algorithm>

#include "Constants.h"
#include "GammaMath.h"
#include "BetaMath.h"
#include "NumericMath.h"

#include "randlib_global.h"

namespace RandMath
{
/**
 * @brief areClose
 * @param a
 * @param b
 * @param eps
 * @return |a - b| < eps * max(a, b)
 */
bool areClose(double a, double b, double eps = 1e-6);

/**
 * @brief sign
 * @param x
 * @return sign of x
 */
int sign(double x);

/**
 * @brief atan
 * @param x
 * @return safe atan(x)
 */
double atan(double x);

/**
 * @brief log1pexp
 * @param x
 * @return log(1 + exp(x))
 */
double log1pexp(double x);

/**
 * @brief log1mexp
 * @param x
 * @return log(1 - exp(x))
 */
double log1mexp(double x);

/**
 * @brief logexpm1
 * @param x
 * @return log(exp(x) - 1)
 */
double logexpm1(double x);

/**
 * @brief log2mexp
 * @param x
 * @return log(2 - exp(x))
 */
double log2mexp(double x);

/**
 * @brief factorial
 * Calculate n! using table values for small n <= 255
 * and built-in gamma function for large n > 255
 * @param n non-negative integer number
 * @return n!
 */
long double factorial(double n);

/**
 * @brief doubleFactorial
 * Calculate n!!
 * @param n
 * @return
 */
long double doubleFactorial(int n);

/**
 * @brief binomialCoef
 * Calculate binomial coefficient C(n,k)
 * @param n
 * @param k
 * @return C(n,k) = n! / (k! * (n - k)!)
 */
long double binomialCoef(int n, int k);

/**
 * @brief erfinv
 * @param p
 * @return inverse error function: such x that erf(x) = p
 */
double erfinv(double p);

/**
 * @brief erfcinv
 * @param p
 * @return inverse complementary error function: such x that erfc(x) = p
 */
double erfcinv(double p);

/**
 * @brief linearInterpolation
 * @param a left boundary
 * @param b right boundary
 * @param fa is equal to f(a)
 * @param fb is qual to f(b)
 * @param x the interpolated function parameter
 * @return f(x)
 */
double linearInterpolation(double a, double b, double fa, double fb, double x);

/**
 * @brief harmonicNumber
 * @param exponent
 * @param number
 * @return sum_{i=1}^{number} i^{-exponent}
 */
double harmonicNumber(double exponent, int number);

/**
 * @brief modifiedBessel
 * @param x
 * @param n
 * @return I_n(x)
 */
double logModifiedBesselFirstKind(double x, double nu);

/**
 * @brief modifiedBesselSecondKind
 * @param x
 * @param n
 * @return K_n(x)
 */
double logModifiedBesselSecondKind(double x, double nu);

/**
 * @brief zetaRiemann
 * @param s
 * @return Riemann zeta function
 */
double zetaRiemann(double s);

/**
 * @brief W0Lambert
 * @param x
 * @param epsilon
 * @return W0 branch of Lambert W function
 */
double W0Lambert(double x, double epsilon = 1e-11);

/**
 * @brief W1Lambert
 * @param x
 * @param epsilon
 * @return W-1 branch of Lambert W function
 */
double Wm1Lambert(double x, double epsilon = 1e-11);

/**
 * @brief MarcumP
 * @param mu
 * @param x
 * @param y
 * @return 1 - Marcum Q-function
 */
double MarcumP(double mu, double x, double y);

/**
 * @brief MarcumQ
 * @param mu
 * @param x
 * @param y
 * @return Marcum Q-function
 */
double MarcumQ(double mu, double x, double y);
}

#endif // RANDMATH_H
