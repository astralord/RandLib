#ifndef RANDMATH_H
#define RANDMATH_H

#include <cmath>
#include <functional>
#include <complex>
#include <vector>
#include <utility>
#include <algorithm>
#include <numeric>

#include "Constants.h"
#include "GammaMath.h"
#include "BetaMath.h"
#include "NumericMath.h"

#include "RandLib_global.h"

namespace RandMath
{
/**
 * @fn areClose
 * @param a
 * @param b
 * @param eps
 * @return |a - b| < eps * max(a, b)
 */
bool areClose(double a, double b, double eps = 1e-6);

/**
 * @fn sign
 * @param x
 * @return sign of x
 */
int sign(double x);

/**
 * @fn atan
 * @param x
 * @return safe atan(x)
 */
double atan(double x);

/**
 * @fn log1pexp
 * @param x
 * @return log(1 + exp(x))
 */
double log1pexp(double x);

/**
 * @fn log1mexp
 * @param x
 * @return log(1 - exp(x))
 */
double log1mexp(double x);

/**
 * @fn logexpm1
 * @param x
 * @return log(exp(x) - 1)
 */
double logexpm1(double x);

/**
 * @fn log2mexp
 * @param x
 * @return log(2 - exp(x))
 */
double log2mexp(double x);

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
 * @param n
 * @return n!!
 */
long double doubleFactorial(int n);

/**
 * @fn erfinv
 * @param p
 * @return inverse error function: such x that erf(x) = p
 */
double erfinv(double p);

/**
 * @fn erfcinv
 * @param p
 * @return inverse complementary error function: such x that erfc(x) = p
 */
double erfcinv(double p);

/**
 * @fn harmonicNumber
 * @param exponent
 * @param number
 * @return sum_{i=1}^{number} i^{-exponent}
 */
double harmonicNumber(double exponent, int number);

/**
 * @fn modifiedBessel
 * @param x
 * @param n
 * @return I_n(x)
 */
double logModifiedBesselFirstKind(double x, double nu);

/**
 * @fn modifiedBesselSecondKind
 * @param x
 * @param n
 * @return K_n(x)
 */
double logModifiedBesselSecondKind(double x, double nu);

/**
 * @fn zetaRiemann
 * @param s
 * @return Riemann zeta function
 */
double zetaRiemann(double s);

/**
 * @fn W0Lambert
 * @param x
 * @param epsilon
 * @return W0 branch of Lambert W function
 */
double W0Lambert(double x, double epsilon = 1e-11);

/**
 * @fn W1Lambert
 * @param x
 * @param epsilon
 * @return W-1 branch of Lambert W function
 */
double Wm1Lambert(double x, double epsilon = 1e-11);

/**
 * @fn MarcumP
 * @param mu
 * @param x
 * @param y
 * @param sqrtX √x
 * @param sqrtY √y
 * @param logX log(x)
 * @param logY log(y)
 * @return 1 - Marcum Q-function
 */
double MarcumP(double mu, double x, double y, double sqrtX, double sqrtY, double logX, double logY);

/**
 * @fn MarcumP
 * @param mu
 * @param x
 * @param y
 * @return 1 - Marcum Q-function
 */
double MarcumP(double mu, double x, double y);

/**
 * @fn MarcumQ
 * @param mu
 * @param x
 * @param y
 * @param sqrtX √x
 * @param sqrtY √y
 * @param logX log(x)
 * @param logY log(y)
 * @return Marcum Q-function
 */
double MarcumQ(double mu, double x, double y, double sqrtX, double sqrtY, double logX, double logY);

/**
 * @fn MarcumQ
 * @param mu
 * @param x
 * @param y
 * @return Marcum Q-function
 */
double MarcumQ(double mu, double x, double y);
}

#endif // RANDMATH_H
