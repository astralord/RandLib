#ifndef RANDMATH_H
#define RANDMATH_H

#include <cmath>
#include <functional>
#include <complex>
#include <vector>
#include <cmath>
#include <utility>

#include "Constants.h"

#include "randlib_global.h"

#include <QDebug>

constexpr long double MIN_POSITIVE = 1e-21l;

typedef std::pair <double, double> DoublePair;
typedef std::pair <int, int> IntPair;

typedef std::tuple <double, double, double> DoubleTriplet;

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
 * @brief digamma
 * @param x
 * @return digamma(x) = d(ln(Gamma(x)))/dx = d(Gamma(x))/Gamma(x)
 */
double digamma(double x);

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

/**
 * @brief betaFun
 * Calculate Beta function
 * @param a
 * @param b
 * @return Gamma(a) * Gamma(b) / Gamma(a + b)
 */
double betaFun(double a, double b);

/**
 * @brief regularizedBetaFun
 * @param x
 * @param a
 * @param b
 * @return
 */
double regularizedBetaFun(double x, double a, double b);

/**
 * @brief incompleteBetaFun
 * @param x
 * @param a
 * @param b
 * @return
 */
double incompleteBetaFun(double x, double a, double b);

/**
 * @brief integral
 * @param funPtr integrand
 * @param a lower boundary
 * @param b upper boundary
 * @param epsilon tolerance
 * @param maxRecursionDepth how deep should the algorithm go
 * @return
 */
long double integral(const std::function<double (double)> &funPtr, double a, double b,
                            double epsilon = 1e-11, int maxRecursionDepth = 11);

/**
 * @brief findRoot
 * Newton's root-finding procedure,
 * using first and second derivatives
 * @param funPtr mapping x |-> (f(x), f'(x), f''(x))
 * @param root starting point and such x that f(x) = 0
 * @param epsilon tolerance
 * @return true if success, false otherwise
 */
bool findRoot(const std::function<DoubleTriplet (double)> &funPtr, double & root, double funTol = 1e-10, double stepTol = 1e-6);

/**
 * @brief findRoot
 * Newton's root-finding procedure,
 * using first derivative
 * @param funPtr mapping x |-> (f(x), f'(x))
 * @param root starting point and such x that f(x) = 0
 * @param epsilon tolerance
 * @return true if success, false otherwise
 */
bool findRoot(const std::function<DoublePair (double)> &funPtr, double & root, double funTol = 1e-10, double stepTol = 1e-6);

/**
 * @brief findRoot
 * Brent's root-finding procedure
 * @param funPtr mapping x |-> f(x)
 * @param a lower boundary
 * @param b upper boundary
 * @param root starting point and such x that f(x) = 0
 * @param epsilon tolerance
 * @return true if success, false otherwise
 */
bool findRoot(const std::function<double (double)> &funPtr, double a, double b, double & root, double epsilon = 1e-8);

/**
 * @brief findMin
 * Combined Brent's method
 * @param funPtr
 * @param a lower boundary
 * @param c upper boundary
 * @param root such x that funPtr(x) is min
 * @param epsilon tolerance
 * @return
 */
bool findMin(const std::function<double (double)> &funPtr, const DoubleTriplet &abc, const DoubleTriplet &fabc, double &root, double epsilon = 1e-8);

/**
 * @brief findMin
 * Combined Brent's method
 * @param funPtr
 * @param closePoint point that is nearby minimum
 * @param root such x that funPtr(x) is min
 * @param epsilon tolerance
 * @return
 */
bool findMin(const std::function<double (double)> &funPtr, double closePoint, double &root, double epsilon = 1e-8);

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
