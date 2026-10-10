#ifndef SPECIALFUNCTIONS_H
#define SPECIALFUNCTIONS_H

/// Special functions that the C++17 standard library provides only optionally.
///
/// libstdc++ and MSVC implement the C++17 mathematical special functions, libc++
/// does not. These wrappers forward to std:: when it is available and fall back to
/// RandLib's own implementations otherwise (or always, with
/// RANDLIB_PORTABLE_SPECIAL_FUNCTIONS defined).

namespace RandMath
{

/**
 * @fn besselI
 * @param nu order, ν ≥ 0
 * @param x argument, x ≥ 0
 * @return modified Bessel function of the first kind I_ν(x)
 */
long double besselI(double nu, double x);

/**
 * @fn besselK
 * @param nu order, ν ≥ 0
 * @param x argument, x > 0
 * @return modified Bessel function of the second kind K_ν(x)
 */
long double besselK(double nu, double x);

/**
 * @fn riemannZeta
 * @param s argument, s > 1
 * @return Riemann zeta function ζ(s)
 */
long double riemannZeta(double s);

/**
 * @fn lgammaCorrection
 * @param x argument, x ≥ 10
 * @return log Γ(x) - ((x - 1/2) log x - x + log √(2π)), the remainder of Stirling's formula
 */
double lgammaCorrection(double x);

namespace portable
{
long double besselI(double nu, double x);
long double besselK(double nu, double x);
long double riemannZeta(double s);
}

}

#endif // SPECIALFUNCTIONS_H
