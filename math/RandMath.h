#ifndef RANDMATH_H
#define RANDMATH_H

#include <QtMath>
#include <cmath>
#include "randlib_global.h"
#include <functional>
#include <QDebug>

#ifndef INFINITY
#include <limits>
static double INFINITY = std::numeric_limits<double>::infinity();
#endif

#ifndef NAN
#include <limits>
static double NAN = std::numeric_limits<double>::quiet_NaN();
#endif

constexpr double M_1_E       = 0.36787944117144232160;
constexpr double M_SQRT3     = 1.73205080756887729353;
constexpr double M_SQRTPI    = 1.77245385090551602730;
constexpr double M_SQRT2PI   = 2.50662827463100050242;
constexpr double M_1_SQRTPI  = 0.56418958354775628695;
constexpr double M_1_SQRT2PI = 0.39894228040143267794;

constexpr double MIN_POSITIVE = 1e-21;

#define SWAP(a, b) (((a) += (b)), ((b) -= (a)), ((a) += (b)), ((b) = -(b)))

/**
 * @brief The RandMath class
 */
class RANDLIBSHARED_EXPORT RandMath
{
public:
    RandMath();

private:
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
    /**
     * @brief fastFactorial
     * Calculate n! using table values for small n <= 255
     * and Moivre-Stirling formula for large n > 255
     * @param n integer number
     * @return n!
     */
    static long double fastFactorial(int n);

    /**
     * @brief fastLogFactorial
     * Calculate log(n!) using table values for small n <= 254
     * and approximate formula for large n > 254
     * @param n integer number
     * @return
     */
    static long double fastLogFactorial(int n);

    /**
     * @brief doubleFactorial
     * Calculate n!!
     * @param n
     * @return
     */
    static long double doubleFactorial(int n);

    /**
     * @brief binomialCoef
     * Calculate binomial coefficient C(n,k)
     * @param n
     * @param k
     * @return C(n,k) = n! / (k! * (n - k)!)
     */
    static long double binomialCoef(int n, int k);

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

    /**
     * @brief erfInv
     * @param p
     * @return inversed error function
     */
    static long double erfInv(double p);

    /**
     * @brief erfcinv
     * @param p
     * @return inversed additional error function
     */
    static long double erfcinv(double p);

private:

    /**
     * @brief adaptiveSimpsonsAux
     * auxiliary function for calculation of integral
     * @param funPtr
     * @param a lower boundary
     * @param b upper boundary
     * @param epsilon
     * @param S
     * @param fa
     * @param fb
     * @param fc
     * @param bottom
     * @return
     */
    static long double adaptiveSimpsonsAux(const std::function<double (double)> &funPtr, double a, double b,
                                           double epsilon, double S, double fa, double fb, double fc, int bottom);

public:

    /**
     * @brief integral
     * @param funPtr integrand
     * @param a lower boundary
     * @param b upper boundary
     * @param epsilon tolerance
     * @param maxRecursionDepth how deep should the algorithm go
     * @return
     */
    static long double integral(const std::function<double (double)> funPtr, double a, double b,
                                double epsilon = 1e-10, int maxRecursionDepth = 10);


    /**
     * @brief findRoot
     * Brent's root-finding procedure
     * @param funPtr
     * @param a lower boundary
     * @param b upper boundary
     * @param root such var as funPtr(var) = 0
     * @param epsilon tolerance
     * @return true if success, false otherwise
     */
    static bool findRoot(const std::function<double (double)> &funPtr, double a, double b, double & root,
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
    static double linearInterpolation(double a, double b, double fa, double fb, double x);
};

#endif // RANDMATH_H
