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

constexpr double MIN_POSITIVE = 1e-21;

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
 * @brief factorial
 * Calculate n! using table values for small n <= 255
 * and Moivre-Stirling formula for large n > 255
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
 * @return psi(x) = d(ln(Gamma(x)))/dx = d(Gamma(x))/Gamma(x)
 */
 double digamma(double x);

/**
 * @brief trigamma
 * @param x
 * @return psi'(x) = d(psi(x))/dx
 */
 double trigamma(double x);

/**
 * @brief lowerIncGamma
 * Calculate lower incomplete gamma function
 * @param a
 * @param x
 * @return
 */
 long double lowerIncGamma(double a, double x);

/**
 * @brief logLowerIncGamma
 * Calculate logarithm of lower incomplete gamma function
 * @param a
 * @param x
 * @return
 */
 long double logLowerIncGamma(double a, double x);

/**
 * @brief upperIncGamma
 * Calculate upper incomplete gamma function
 * @param a
 * @param x
 * @return
 */
 long double upperIncGamma(double a, double x);

/**
 * @brief upperIncGamma
 * Calculate logarithm of upper incomplete gamma function
 * @param a
 * @param x
 * @return
 */
long double logUpperIncGamma(double a, double x);

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
 * @brief gammaHalf
 * Calculate Gamma function of half integer: gamma(k/2)
 * @param k
 * @return
 */
long double gammaHalf(size_t k);

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
 * Newton's root-finding procedure
 * @param funPtr
 * @param root
 * @param epsilon
 * @return true if success, false otherwise
 */
bool findRoot(const std::function<DoubleTriplet (double)> &funPtr, double & root, double epsilon = 1e-10);

/**
 * @brief findRoot
 * Newton's root-finding procedure
 * @param funPtr
 * @param root starting point
 * @param epsilon
 * @return true if success, false otherwise
 */
bool findRoot(const std::function<DoublePair (double)> &funPtr, double & root, double epsilon = 1e-10);

/**
 * @brief findRoot
 * Brent's root-finding procedure
 * @param funPtr
 * @param a lower boundary
 * @param b upper boundary
 * @param root such x that funPtr(x) = 0
 * @param epsilon tolerance
 * @return true if success, false otherwise
 */
bool findRoot(const std::function<double (double)> &funPtr, double a, double b, double & root,
                     double epsilon = 1e-10);

/**
 * @brief findMin
 * Combined Brent's method
 * @param funPtr
 * @param a lower boundary
 * @param b upper boundary
 * @param root such x that funPtr(x) is min
 * @param epsilon tolerance
 * @return
 */
bool findMin(const std::function<double (double)> &funPtr, double a, double b, double & root,
                    double epsilon = 1e-10);

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
double modifiedBesselFirstKind(double x, double n);

/**
 * @brief modifiedBesselSecondKind
 * @param x
 * @param n
 * @return K_n(x)
 */
double modifiedBesselSecondKind(double x, double n);

/**
 * @brief BernoulliNumber
 * @param n
 * @return Bernoulli number, calculated by Akiyama–Tanigawa algorithm
 */
double BernoulliNumber(int n);

/**
 * @brief zetaRiemann
 * @param s
 * @return Riemann zeta function
 */
double zetaRiemann(double s);
}

#endif // RANDMATH_H
