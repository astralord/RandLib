#include "BetaMath.h"

namespace RandMath {

long double logBeta(double a, double b)
{
    if (a <= 0 || b <= 0)
        return NAN;
    double apb = a + b;
    int roundA = std::round(a), roundB = std::round(b);
    int roundApB = std::round(apb); 
    long double lgammaA = areClose(a, roundA) ? lfact(roundA - 1) : std::lgammal(a);
    long double lgammaB = lgammaA;
    if (a != b)
        lgammaB = areClose(b, roundB) ? lfact(roundB - 1) : std::lgammal(b);
    long double lgammaApB = areClose(apb, roundApB) ? lfact(roundApB - 1) : std::lgammal(apb);
    return lgammaA + lgammaB - lgammaApB;
}

long double beta(double a, double b)
{
    return std::exp(logBeta(a, b));
}

double ibetaPowerSeries1(double x, double a, double b, double logBetaFun, double logX, double log1mX)
{
    double addon = 1.0, sum = 0.0;
    /// avoid underflow
    double u = (x > 1e-5) ? -x / (1.0 - x) : -std::exp(logX - log1mX);
    int n = 1;
    do {
        addon *= (n - b) * u / (a + n);
        sum += addon;
        ++n;
    } while (std::fabs(addon) > MIN_POSITIVE * std::fabs(sum));
    double y = a * logX;
    y += (b - 1) * log1mX;
    y -= logBetaFun;
    y += std::log1p(sum);
    return std::exp(y) / a;
}

double ibetacontinuedFraction1(double x, double a, double b, int number, double logBetaFun, double logX, double log1mX)
{
    int n = number;
    /// avoid underflow
    double u = (x > 1e-5) ? x / (1.0 - x) : std::exp(logX - log1mX);
    double frac = 0.0;
    while (n > 1) {
        double cn = (b - n + 1) / (a + n - 1) * u;
        frac = cn / (1 + cn - frac);
        --n;
    }
    double c1 = b / a;
    double F = c1 / (1 + c1 - frac);
    double y = a * logX;
    y += (b - 1) * log1mX;
    y -= logBetaFun;
    y = std::exp(y) / b;
    y *= F / (1.0 - F);
    return y;
}

double ibetaPowerSeries2(double x, double a, double b, double logBetaFun, double logX, double log1mX)
{
    double addon = 1.0, sum = 0.0;
    int n = 0;
    do {
        addon *= (n + a + b) * x / (a + n + 1);
        sum += addon;
        ++n;
    } while (std::fabs(addon) > MIN_POSITIVE * std::fabs(sum));
    double y = a * logX;
    y += b * log1mX;
    y -= logBetaFun;
    y += std::log1p(sum);
    return std::exp(y) / a;
}

double ibetaContinuedFraction2(double x, double a, double b, int number, double logBetaFun, double logX, double log1mX)
{
    int n = number;
    double frac = 0.0;
    while (n > 1) {
        double cn = (a + b + n - 2) / (a + n - 1) * x;
        frac = cn / (1 + cn - frac);
        --n;
    }
    double c1 = (a + b - 1) / a;
    double F = c1 / (1 + c1 - frac);
    double y = a * logX;
    y += b * log1mX;
    y -= logBetaFun;
    y = std::exp(y) / (a + b - 1);
    y *= F / (1.0 - F);
    return y;
}

double ibeta(double x, double a, double b, double logBetaFun, double logX, double log1mX)
{
    /// Check special values
    if (a <= 0 || b <= 0 || x < 0.0 || x > 1.0)
        return NAN;
    if (x == 0.0)
        return 0.0;
    if (x == 1.0)
        return 1.0;
    if (b == 1.0)
        return std::exp(a * logX);
    if (a == 1.0)
        return -std::expm1(b * log1mX);

    /// If x is larger than mean of Beta distribution,
    /// convergence of complementary distribution function is faster
    if (x > a / (a + b))
        return 1.0 - ibeta(1.0 - x, b, a);

    /// We truncate after 0.45 instead of 0.5 here
    /// in order to avoid x being too close to 0.5
    if (x < 0.45) {
        if (b < 1)
            return ibetaPowerSeries1(x, a, b, logBetaFun, logX, log1mX);
        double number = std::floor(b);
        /// equivalent continued fraction #1
        double ecd = ibetacontinuedFraction1(x, a, b, number, logBetaFun, logX, log1mX);
        double apn = a + number, bmn = b - number;
        double residue = (b == number) ? 0.0 : ibetaPowerSeries1(x, apn, bmn, logBeta(apn, bmn), logX, log1mX);
        return ecd + residue;
    }

    if (b < 1)
        return ibetaPowerSeries2(x, a, b, logBetaFun, logX, log1mX);
    double number = std::floor(b);
    /// equivalent continued fraction #2
    double ecd = ibetaContinuedFraction2(x, a, b, number, logBetaFun, logX, log1mX);
    double apn = a + number;
    double residue = ibetaPowerSeries2(x, apn, b, logBeta(apn, b), logX, log1mX);
    return ecd + residue;
}

double ibeta(double x, double a, double b)
{
    /// Check special values
    if (a <= 0 || b <= 0 || x < 0.0 || x > 1.0)
        return NAN;
    if (x == 0.0)
        return 0.0;
    if (x == 1.0)
        return 1.0;
    if (b == 1.0)
        return std::pow(x, a);
    if (a == 1.0) {
        double y = b * std::log1p(-x);
        return -std::expm1(y);
    }
    double logBetaFun = logBeta(a, b);
    double logX = std::log(x), log1mX = std::log1p(-x);
    return ibeta(x, a, b, logBetaFun, logX, log1mX);
}

}
