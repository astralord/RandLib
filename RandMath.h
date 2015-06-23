#ifndef RANDMATH_H
#define RANDMATH_H

#include <QtMath>
#include <cmath>
#include "randlib_global.h"

static constexpr double M_SQRT3     = 1.73205080756887729353;
static constexpr double M_SQRTPI    = 1.77245385090551602730;
static constexpr double M_SQRT2PI   = 2.50662827463100050242;
static constexpr double M_1_SQRTPI  = 0.56418958354775628695;
static constexpr double M_1_SQRT2PI = 0.39894228040143267794;
static constexpr double M_1_E       = 0.36787944117144232160;

static constexpr double MIN_POSITIVE = 1e-21;

class RANDLIBSHARED_EXPORT RandMath
{
    /**
     * @brief getTenthFactorial
     * @param n = 10 * k <= 250
     * @return exact value of n!
     */
    static long double getTenthFactorial(int n);

    /**
     * @brief stirlingFactorial
     * @param n integer number
     * @return n! according to Moivre-Stirling formula
     */
    static long double stirlingFactorial(int n);

public:
    RandMath();

    /**
     * @brief fastFactorial
     * Calculate n! using table values for small n <= 250
     * and Moivre-Stirling formula for large n > 250
     * @param n integer number
     * @return n!
     */
    static long double fastFactorial(int n);

    /**
     * @brief doubleFactorial
     * Calculate n!!
     * @param n
     * @return
     */
    static long double doubleFactorial(int n);

    /**
     * @brief lowerIncGamma
     * Calculate lower incomplete gamma function
     * @param a
     * @param x
     * @return
     */
    static long double lowerIncGamma(double a, double x);

    /**
     * @brief upperIncGamma
     * Calculate upper incomplete gamma function
     * @param a
     * @param x
     * @return
     */
    static long double upperIncGamma(double a, double x);

    /**
     * @brief betaFun
     * Calculate Beta function
     * @param x
     * @param y
     * @return Gamma(x) * Gamma(y) / Gamma(x + y)
     */
    static long double betaFun(double x, double y);

    /**
     * @brief gammaHalf
     * Calculate Gamma function of half integer: gamma(k/2)
     * @param k
     * @return
     */
    static long double gammaHalf(int k);
};

#endif // RANDMATH_H
