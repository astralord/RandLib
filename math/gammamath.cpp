#include "GammaMath.h"

namespace RandMath
{

/**
 * @brief digammamLogForLargeX
 * @param x
 * @return digamma(x) - log(x) for x > 7
 */
double digammamLogForLargeX(double x)
{
    /// set up tables
    static constexpr long double taylorCoef[] = {
         0.00833333333333333333l,
        -0.00396825396825396825l,
         0.00416666666666666666l,
        -0.00757575757575757576l,
         0.02115090296908478727l,
        -0.08333333333333333333l
    };
    static constexpr int bounds[] = {
        6000, 320, 75, 30, 20, 10, 7
    };

    /// choose the amount of terms in series
    int degree = 0;
    while (x < bounds[degree])
        ++degree;

    /// apply
    long double firstTerm = taylorCoef[5] / x - 0.5;
    long double y = firstTerm / x;
    for (int i = 0; i < degree; ++i)
        y += taylorCoef[i] / std::pow(x, 2 * i + 4);
    return y;
}

double digamma(double x)
{
    /// Negative argument
    if (x < 0.0) {
        double y = 1.0 - x;
        double z = M_PI / std::tan(M_PI * y);
        return digamma(y) + z;
    }
    /// shift to minimum value,
    /// for which series expansion is applicable
    long double sum = 0;
    while (x < 7.0) {
        sum -= 1.0 / x;
        ++x;
    }
    return sum + digammamLogForLargeX(x) + std::log(x);
}

double digammamLog(double x)
{
    if (x < 0.0)
        return NAN;
    if (x == 0.0)
        return -INFINITY;
    return (x < 7.0) ? digamma(x) - std::log(x) : digammamLogForLargeX(x);
}

double trigamma(double x)
{
    /// Negative argument
    if (x < 0.0)
    {
        double z = M_PI / std::sin(M_PI * x);
        return z * z - trigamma(1.0 - x);
    }
    /// set up tables
    static constexpr long double taylorCoef[] = {
         0.16666666666666666666l,
        -0.03333333333333333333l,
         0.02380952380952380952l,
        -0.03333333333333333333l,
         0.07575757575757575757l,
        -0.25311355311355311355l
    };
    static constexpr int bounds[] = {
        1000, 140, 47, 25, 15, 10
    };

    /// shift to minimum value,
    /// for which series expansion is applicable
    long double y = 0;
    while (x < 10.0) {
        y += 1.0 / (x * x);
        ++x;
    }

    /// choose the amount of terms in series
    int degree = 1;
    while (x < bounds[degree - 1])
        ++degree;

    /// apply
    long double firstTerm = 1.0 / x;
    firstTerm += 0.5 / (x * x);
    firstTerm += taylorCoef[0] / std::pow(x, 3);
    y += firstTerm;
    for (int i = 1; i < degree; ++i)
        y += taylorCoef[i] / std::pow(x, 2 * i + 3);
    return y;
}

enum REGULARISED_GAMMA_METHOD_ID {
    PT,
    QT,
    PUA,
    QUA,
    CF,
    UNDEFINED
};

double pgammaRaw(double a, double x, double logX, REGULARISED_GAMMA_METHOD_ID mId);
double lpgammaRaw(double a, double x, double logX, REGULARISED_GAMMA_METHOD_ID mId);
double qgammaRaw(double a, double x, double logX, REGULARISED_GAMMA_METHOD_ID mId);
double lqgammaRaw(double a, double x, double logX, REGULARISED_GAMMA_METHOD_ID mId);

REGULARISED_GAMMA_METHOD_ID getRegularizedGammaMethodId(double a, double x, double logX)
{
    if (x < 0.0 || a < 0.0)
        return UNDEFINED;
    double alpha = (x < 0.5) ? M_LN2 / (M_LN2 - logX) : x;
    if ((a <= 12 && a >= alpha) || 0.3 * a >= x)
        return PT;
    if (x <= 1.5 && a <= alpha)
        return QT;
    if (a <= 12 || 2.35 * a <= x)
        return CF;
    return (a > alpha) ? PUA : QUA;
}

double incompleteGammaUniformExpansion(double a, double x, double logX, bool isP)
{
    /// Uniform asymptotic expansion for P(a, x) if isP == true, or Q(a, x) otherwise
    static constexpr long double d[] = {-0.33333333333333333333l, 0.08333333333333333333l, -0.01481481481481481481l, 0.00115740740740740741l,
                                         0.00035273368606701940l, -0.000178755144032922l, 0.0000391926317852244l, -0.00000218544851067999l,
                                        -0.00000185406221071516l, 0.829671134095309e-6l, -0.176659527368261e-6l, 0.670785354340150e-8l,
                                         0.102618097842403e-7l, -0.438203601845335e-8l, 0.914769958223678e-9l, -0.255141939949460e-10l,
                                        -0.583077213255043e-10l, 0.243619480206674e-10l, -0.502766928011417e-11l, 0.110043920319559e-12l,
                                         0.337176326240099e-12l, -0.139238872241816e-12l, 0.285348938070474e-13l, -0.513911183424242e-15l,
                                        -0.197522882943494428e-16l, 0.809952115670456133e-17};
    static constexpr int N = 25;
    double lambda = x / a;
    double logA = std::log(a);
    double logLambda = logX - logA;
    double aux = x - a - a * logLambda;
    double etaSq = 0.0, base = 0.5;
    if (aux > 0.0) { /// otherwise, x ~ a and aux ~ 0.0
        etaSq = 2 * (lambda - 1.0 - logLambda);
        base = 0.5 * std::erfc(std::sqrt(aux));
    }
    long double sum = 0.0l;
    double betanp2 = d[N], betanp1 = d[N - 1];
    for (int n = N - 2; n >= 0; --n) {
        double beta = (n + 2) * betanp2 / a;
        beta += (x < a && (n & 1)) ? -d[n] : d[n];
        sum += beta * std::pow(etaSq, 0.5 * n);
        betanp2 = betanp1;
        betanp1 = beta;
    }
    sum *= a / (a + betanp2);
    double z = aux + 0.5 * (M_LN2 + M_LNPI + logA);
    double y = std::exp(-z) * sum;
    return base + (isP ? -y : y);
}

double lpgammaRaw(double a, double x, double logX, REGULARISED_GAMMA_METHOD_ID mId)
{
    if (mId == PT)
    {
        /// Taylor expansion of P(a,x)
        int n0 = 70.0 * x / a + 7; /// ~ from 7 to 28
        long double sum = 0.0;
        double lgammaAp1 = std::lgamma(a + 1);
        for (int n = n0; n > 0; --n) {
            double addon = n * logX - std::lgamma(a + n + 1) + lgammaAp1;
            addon = std::exp(addon);
            sum += addon;
        }
        return a * logX - x + std::log1p(sum) - lgammaAp1;
    }
    return (mId == PUA) ? std::log(pgammaRaw(a, x, logX, mId)) : std::log1p(-qgammaRaw(a, x, logX, mId));
}

double lpgamma(double a, double x, double logX)
{
    if (x < 0.0 || a < 0.0)
        return NAN;
    if (x == 0.0)
        return 0.0;
    if (a == 1.0)
        return std::log1p(-std::exp(-x));
    REGULARISED_GAMMA_METHOD_ID mId = getRegularizedGammaMethodId(a, x, logX);
    return lpgammaRaw(a, x, logX, mId);
}

double lpgamma(double a, double x)
{
    if (x < 0.0 || a < 0.0)
        return NAN;
    if (x == 0.0)
        return 0.0;
    if (a == 1.0)
        return std::log1p(-std::exp(-x));
    return lpgamma(a, x, std::log(x));
}

double pgammaRaw(double a, double x, double logX, REGULARISED_GAMMA_METHOD_ID mId)
{
    if (mId == PT)
        return std::exp(lpgammaRaw(a, x, logX, mId));
    if (mId == PUA)
        return incompleteGammaUniformExpansion(a, x, logX, true);
    return (mId == QUA) ? 1.0 - qgammaRaw(a, x, logX, mId) : -std::expm1(lqgammaRaw(a, x, logX, mId));
}

double pgamma(double a, double x, double logX)
{
    if (x < 0.0 || a < 0.0)
        return NAN;
    if (x == 0.0)
        return 1.0;
    if (a == 1.0)
        return -std::expm1(-x);
    REGULARISED_GAMMA_METHOD_ID mId = getRegularizedGammaMethodId(a, x, logX);
    return pgammaRaw(a, x, logX, mId);
}

double pgamma(double a, double x)
{
    if (x < 0.0 || a < 0.0)
        return NAN;
    if (x == 0.0)
        return 1.0;
    if (a == 1.0)
        return -std::expm1(-x);
    return pgamma(a, x, std::log(x));
}

double qtGammaExpansionAux(double a, double logX)
{
    /// auxilary function for Taylor expansion of Q(a, x)
    long double sum = 0.0;
    for (int n = 1; n != 20; ++n) {
        double addon = n * logX - std::lgamma(n + 1);
        addon = std::exp(addon);
        addon /= (a + n);
        sum += (n & 1) ? -addon : addon;
    }
    sum *= a;
    double y = std::log1p(sum);
    y += a * logX;
    y -= std::lgamma(a + 1);
    return y;
}

double lqgammaRaw(double a, double x, double logX, REGULARISED_GAMMA_METHOD_ID mId)
{
    if (mId == QT)
    {
        double y = qtGammaExpansionAux(a, logX);
        y = -std::exp(y);
        return std::log1p(y);
    }
    if (mId == CF)
    {
        /// Continued fraction
        int k0 = std::min(40.0 / (x - 1) + 5.0, 60.0);
        long double sum = 0.0;
        double rhok = 0.0, tk = 1.0;
        for (int k = 1; k <= k0; ++k) {
            /// Calculate a(k)
            double ak = k * (a - k);
            double temp = x + 2 * k - a;
            ak /= temp * temp - 1;
            /// Calculate rho(k)
            ++rhok;
            rhok *= ak;
            rhok /= -(1 + rhok);
            /// Calculate t(k) and add it to the sum
            tk *= rhok;
            sum += tk;
        }
        double y = std::log1p(sum);
        y += a * logX;
        y -= x + std::lgamma(a);
        y -= std::log1p(x - a);
        return y;
    }
    return (mId == QUA) ? std::log(qgammaRaw(a, x, logX, mId)) : std::log1p(-pgammaRaw(a, x, logX, mId));
}

double lqgamma(double a, double x, double logX)
{
    if (x < 0.0 || a < 0.0)
        return NAN;
    if (x == 0.0)
        return -INFINITY;
    if (a == 1.0)
        return -x;
    REGULARISED_GAMMA_METHOD_ID mId = getRegularizedGammaMethodId(a, x, logX);
    return lqgammaRaw(a, x, logX, mId);
}

double lqgamma(double a, double x)
{
    if (x < 0.0 || a < 0.0)
        return NAN;
    if (x == 0.0)
        return -INFINITY;
    if (a == 1.0)
        return -x;
    return lqgamma(a, x, std::log(x));
}

double qgammaRaw(double a, double x, double logX, REGULARISED_GAMMA_METHOD_ID mId)
{
    if (mId == CF)
        return std::exp(lqgammaRaw(a, x, logX, mId));
    if (mId == QT)
        return -std::expm1(qtGammaExpansionAux(a, logX));
    if (mId == QUA)
        return incompleteGammaUniformExpansion(a, x, logX, false);
    return (mId == PUA) ? 1.0 - pgammaRaw(a, x, logX, mId) : -std::expm1(lpgammaRaw(a, x, logX, mId));
}

double qgamma(double a, double x, double logX)
{
    if (x < 0.0 || a < 0.0)
        return NAN;
    if (x == 0.0)
        return 0.0;
    if (a == 1.0)
        return std::exp(-x);
    REGULARISED_GAMMA_METHOD_ID mId = getRegularizedGammaMethodId(a, x, logX);
    return qgammaRaw(a, x, logX, mId);
}

double qgamma(double a, double x)
{
    if (x < 0.0 || a < 0.0)
        return NAN;
    if (x == 0.0)
        return 0.0;
    if (a == 1.0)
        return std::exp(-x);
    return qgamma(a, x, std::log(x));
}

}
