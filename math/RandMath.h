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

namespace RandMath
{

template <int TABLE_SIZE>
constexpr std::array<LongDoublePair, TABLE_SIZE> createZiggurat(const std::function<LongDoublePair (const LongDoublePair &)> &funPtr, LongDoublePair first, LongDoublePair second);

/**
 * @fn areClose
 * @param a
 * @param b
 * @param eps
 * @return |a - b| < eps * max(a, b)
 */
template<typename RealType>
bool areClose(RealType a, RealType b, RealType eps = 1e-6)
{
    if (a == b)
        return true;
    RealType fa = std::fabs(a);
    RealType fb = std::fabs(b);
    return std::fabs(b - a) < eps * std::max(fa, fb);
}

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
double softplus(double x);

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
double logexpm1l(double x);

/**
 * @fn log2mexp
 * @param x
 * @return log(2 - exp(x))
 */
double log2mexp(double x);

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
 * @brief xexpxsqerfc
 * @param x
 * @return x * exp(x^2) * erfc(x)
 */
double xexpxsqerfc(double x);

/**
 * @fn harmonicNumber
 * @param exponent
 * @param number
 * @return sum_{i=1}^{number} i^{-exponent}
 */
double harmonicNumber(double exponent, int number);

/**
 * @fn logBesselI
 * Calculates logarithm of modified Bessel function of the 1st kind
 * @param nu
 * @param x
 * @return log(I_ν(x))
 */
long double logBesselI(double nu, double x);

/**
 * @fn logBesselK
 * Calculates logarithm of modified Bessel function of the 2nd kind
 * @param nu
 * @param x
 * @return log(K_ν(x))
 */
long double logBesselK(double nu, double x);

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
