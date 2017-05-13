#ifndef GAMMAMATH
#define GAMMAMATH

#include "RandMath.h"

/// Gamma-related functions

namespace RandMath
{

/**
 * @fn lfact
 * @param n
 * @return log(n!)
 */
double lfact(size_t n);

/**
 * @fn digamma
 * @param x
 * @return digamma(x) = ln(Γ(x))' = d(Γ(x))/Γ(x)
 */
double digamma(double x);

/**
 * @fn logmdigamma
 * @param x
 * @return digamma(x) - log(x)
 */
double digammamLog(double x);

/**
 * @fn trigamma
 * @param x
 * @return trigamma(x) = digamma'(x)
 */
double trigamma(double x);

/**
 * @fn lpgamma
 * Calculate logarithm of lower incomplete gamma function,
 * accelerated by using precalculated value of log(a) and log(Γ(a))
 * @param a
 * @param x
 * @param logA log(a)
 * @param lgammaA log(Γ(a))
 * @return log(P(a, x))
 */
double lpgamma(double a, double x, double logA, double lgammaA);

/**
 * @fn lpgamma
 * Calculate logarithm of lower incomplete gamma function,
 * accelerated by using precalculated value of log(x)
 * @param a
 * @param x
 * @param logX log(x)
 * @return log(P(a, x))
 */
double lpgamma(double a, double x, double logX);

/**
 * @fn lpgamma
 * Calculate logarithm of lower incomplete gamma function
 * @param a
 * @param x
 * @return log(P(a, x))
 */
double lpgamma(double a, double x);

/**
 * @fn pgamma
 * Calculate lower regularised incomplete gamma function,
 * accelerated by using precalculated value of log(a) and log(Γ(a))
 * @param a
 * @param x
 * @param logA log(a)
 * @param lgammaA log(Γ(a))
 * @return P(a, x)
 */
double pgamma(double a, double x, double logA, double lgammaA);

/**
 * @fn pgamma
 * Calculate lower regularised incomplete gamma function,
 * accelerated by using precalculated value of log(x)
 * @param a
 * @param x
 * @param logX log(x)
 * @return P(a, x)
 */
double pgamma(double a, double x, double logX);

 /**
 * @fn pgamma
 * Calculate lower regularised incomplete gamma function
 * @param a
 * @param x
 * @return P(a, x)
 */
double pgamma(double a, double x);

/**
 * @fn lqgamma
 * Calculate logarithm of upper incomplete gamma function,
 * accelerated by using precalculated values of log(a) and log(Γ(a))
 * @param a
 * @param x
 * @param logA log(a)
 * @param lgammaA log(Γ(a))
 * @return log(Q(a, x))
 */
double lqgamma(double a, double x, double logA, double lgammaA);

/**
 * @fn lqgamma
 * Calculate logarithm of upper incomplete gamma function,
 * accelerated by using precalculated value of log(x)
 * @param a
 * @param x
 * @param logX log(x)
 * @return log(Q(a, x))
 */
double lqgamma(double a, double x, double logX);

 /**
 * @fn lqgamma
 * Calculate logarithm of upper incomplete gamma function
 * @param a
 * @param x
 * @return log(Q(a, x))
 */
double lqgamma(double a, double x);

/**
 * @fn lqgamma
 * Calculate upper incomplete gamma function,
 * accelerated by using precalculated values of log(a) and log(Γ(a))
 * @param a
 * @param x
 * @param logA log(a)
 * @param lgammaA log(Γ(a))
 * @return Q(a, x)
 */
double qgamma(double a, double x, double logA, double lgammaA);

/**
 * @fn qgamma
 * Calculate upper incomplete gamma function,
 * accelerated by using precalculated value of log(x)
 * @param a
 * @param x
 * @param logX log(x)
 * @return Q(a, x)
 */
double qgamma(double a, double x, double logX);

/**
 * @fn qgamma
 * Calculate upper incomplete gamma function
 * @param a
 * @param x
 * @return Q(a, x)
 */
double qgamma(double a, double x);

}

#endif // GAMMAMATH

