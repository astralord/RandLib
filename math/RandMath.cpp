#include "RandMath.h"
#include <functional>

namespace RandMath
{

bool areClose(double a, double b, double eps)
{
    if (a == b)
        return true;
    double fa = std::fabs(a);
    double fb = std::fabs(b);
    return std::fabs(b - a) < eps * std::max(fa, fb);
}

int sign(double x)
{
    return (x > 0) ? 1 : ((x < 0) ? -1 : 0);
}

double atan(double x)
{
    /// For small absolute values we use standard technique
    /// Otherwise we use relation
    /// atan(x) = +/-π/2 - atan(1/x)
    /// to avoid numeric problems
    if (x == 0.0)
        return 0.0;
    if (x > 1.0)
        return M_PI_2 - std::atan(1.0 / x);
    return (x < -1.0) ? -M_PI_2 - std::atan(1.0 / x) : std::atan(x);
}

double log1pexp(double x)
{
    if (x < 20.0)
        return std::log1p(std::exp(x));
    return (x < 35.0) ? x + std::exp(-x) : x;
}

double log1mexp(double x)
{
    return (x < -M_LN2) ? std::log1p(-std::exp(x)) : std::log(-std::expm1(x));
}

double logexpm1(double x)
{
    if (x < 20.0)
        return std::log(std::expm1(x));
    return (x < 35.0) ? x - std::exp(-x) : x;
}

double log2mexp(double x)
{
    return std::log1p(-std::expm1(x));
}

double erfinvChebyshevSeries(double x, long double t, const long double *array, int size)
{
    /// We approximate inverse erf via Chebyshev polynomials
    long double Tn = t, Tnm1 = 1, sum = 0.0;
    for (int i = 1; i != size; ++i) {
        sum += array[i] * Tn;
        /// calculate next Chebyshev polynomial
        double temp = Tn;
        Tn *= 2 * t;
        Tn -= Tnm1;
        Tnm1 = temp;
    }
    return x * (array[0] + sum);
}

double erfinvAux1(double beta)
{
    /// |1 - p| < 5e-16
    static constexpr long double D5 = -9.199992358830151031278420l, D6 = 2.794990820124599493768426l;
    static constexpr int muSize = 26;
    static constexpr long double mu[muSize] = {.9885750640661893136460358l, .0108577051845994776160281l, -.0017511651027627952594825l, .0000211969932065633437984l,
                                               .0000156648714042435087911l, -.05190416869103124261e-5l, -.00371357897426717780e-5l, .00012174308662357429e-5l,
                                              -.00001768115526613442e-5l, -.119372182556161e-10l, .003802505358299e-10l, -.000660188322362e-10l, -.000087917055170e-10l,
                                              -.3506869329e-15l, -.0697221497e-15l, -.0109567941e-15l, -.0011536390e-15l, -.0000263938e-15l, .05341e-20l, -.22610e-20l,
                                               .09552e-20l, -.05250e-20l, .02487e-20l, -.01134e-20l, .00420e-20l};
    return erfinvChebyshevSeries(beta, D5 / std::sqrt(beta) + D6, mu, muSize);
}

double erfinvAux2(double beta)
{
    /// 5e-16 < |1 - p| < 0.025
    static constexpr long double D3 = -0.5594576313298323225436913l, D4 = 2.287915716263357638965891l;
    static constexpr int deltaSize = 38;
    static constexpr long double delta[deltaSize] = {.9566797090204925274526373l, -.0231070043090649036999908l, -.0043742360975084077333218l, -.0005765034226511854809364l,
                                                    -.0000109610223070923931242l, .0000251085470246442787982l, .0000105623360679477511955l, .27544123300306391503e-5l,
                                                     .04324844983283380689e-5l, -.00205303366552086916e-5l, -.00438915366654316784e-5l, -.00176840095080881795e-5l,
                                                    -.00039912890280463420e-5l, -.00001869324124559212e-5l, .00002729227396746077e-5l, .00001328172131565497e-5l,
                                                     .318342484482286e-10l, .016700607751926e-10l, -.020364649611537e-10l, -.009648468127965e-10l, -.002195672778128e-10l,
                                                    -.000095689813014e-10l, .000137032572230e-10l, .000062538505417e-10l, .000014584615266e-10l, .1078123993e-15l,
                                                    -.0709229988e-15l, -.0391411775e-15l, -.0111659209e-15l, -.0015770366e-15l, .0002853149e-15l, .0002716662e-15l,
                                                     .0000176835e-15l, .09828e-20l, .20464e-20l, .08020e-20l, .01650e-20l};
    return erfinvChebyshevSeries(beta, D3 * beta + D4, delta, deltaSize);
}

double erfinvAux3(double beta)
{
    /// 0.8 < p < 0.975
    static constexpr long double D1 = -1.548813042373261659512742l, D2 = 2.565490123147816151928163l;
    static constexpr int lambdaSize = 27;
    static constexpr long double lambda[lambdaSize] = {.9121588034175537733059200l, -.0162662818676636958546661l, .0004335564729494453650589l, .0002144385700744592065205l,
                                                       .26257510757648130176e-5l, -.30210910501037969912e-5l, -.00124060618367572157e-5l, .00624066092999917380e-5l,
                                                      -.00005401247900957858e-5l, -.00014232078975315910e-5l, .343840281955305e-10l, .335848703900138e-10l, -.014584288516512e-10l,
                                                      -.008102174258833e-10l, .000525324085874e-10l, .000197115408612e-10l, -.000017494333828e-10l, -.4800596619e-15l, .0557302987e-15l,
                                                       .0116326054e-15l, -.0017262489e-15l, -.0002784973e-15l, .0000524481e-15l, .65270e-20l, -.15707e-20l, -.01475e-20l, .00450e-20l};
    return erfinvChebyshevSeries(beta, D1 * beta + D2, lambda, lambdaSize);
}

double erfinvAux4(double p)
{
    /// 0 < p < 0.8
    static constexpr int xiSize = 39;
    static constexpr long double xi[xiSize] = {.9928853766189408231495800l, .1204675161431044864647846l, .0160781993420999447267039l, .0026867044371623158279591l,
                                               .0004996347302357262947170l, .0000988982185991204409911l, .0000203918127639944337340l, .43272716177354218758e-5l,
                                               .09380814128593406758e-5l, .02067347208683427411e-5l, .00461596991054300078e-5l, .00104166797027146217e-5l,
                                               .00023715009995921222e-5l, .00005439284068471390e-5l, .00001255489864097987e-5l, .291381803663201e-10l,
                                               .067949421808797e-10l, .015912343331569e-10l, .003740250585245e-10l, .000882087762421e-10l, .000208650897725e-10l,
                                               .000049488041039e-10l, .000011766394740e-10l, .2803855725e-15l, .0669506638e-15l, .0160165495e-15l, .0038382583e-15l,
                                               .0009212851e-15l, .0002214615e-15l, .0000533091e-15l, .0000128488e-15l, .31006e-20l, .07491e-20l, .01812e-20l,
                                               .00439e-20l, .00106e-20l, .00026e-20l, .00006e-20l, .00002e-20l};
    return erfinvChebyshevSeries(p, p * p / 0.32 - 1.0, xi, xiSize);
}

double erfinv(double p)
{
    /// Consider special cases
    if (p < 0.0)
        return -erfinv(-p);
    if (p > 1.0)
        return NAN;
    if (p == 1.0)
        return INFINITY;
    if (p == 0.0)
        return 0.0;
    if (p < 0.8)
        return erfinvAux4(p);
    /// Handle tails
    double beta = std::sqrt(-std::log1p(-p * p));
    if (p < 0.9975)
        return erfinvAux3(beta);
    return (1.0 - p < 5e-16) ? erfinvAux1(beta) : erfinvAux2(beta);
}

double erfcinv(double p)
{
    /// Consider special cases
    if (p > 1.0)
        return -erfcinv(2.0 - p);
    if (p == 0.0)
        return INFINITY;
    if (p == 1.0)
        return 0.0;
    if (p > 0.2)
        return erfinvAux4(1.0 - p);
    double pSq = p * p, p2 = 2 * p;
    double beta = std::sqrt(-std::log(p2 - pSq));
    if (p > 0.0025)
        return erfinvAux3(beta);
    return (p > 5e-16) ? erfinvAux2(beta) : erfinvAux1(beta);
}

double xexpxsqerfc(double x)
{
    static constexpr int MAX_X = 10;
    static constexpr int N = 10;
    if (x < MAX_X) {
        double y = x * x;
        y += std::log(std::erfc(x));
        return x * std::exp(y);
    }
    double log2xSq = M_LN2 + 2 * std::log(x);
    double sum = 0.0;
    for (int n = 1; n != N; ++n) {
        double add = RandMath::ldfact(2 * n - 1);
        add -= n * log2xSq;
        add = std::exp(add);
        sum += (n & 1) ? -add : add;
    }
    return (1.0 + sum) / M_SQRTPI;
}

double harmonicNumber(double exponent, int number)
{
    if (number < 1)
        return 0;
    if (exponent == 1)
        return M_EULER + digamma(number + 1);
    if (exponent == 2)
        return M_PI_SQ / 6.0 - trigamma(number + 1);
    double res = 1.0;
    for (int i = 2; i <= number; ++i)
        res += std::pow(i, -exponent);
    return res;
}

long double logBesselI(double nu, double x)
{
    if (x < 0) {
        double roundNu = std::round(nu);
        bool nuIsInt = areClose(nu, roundNu);
        if (nuIsInt) {
            int nuInt = roundNu;
            return (nuInt % 2) ? NAN : logBesselI(nu, -x);
        }
        return -INFINITY;
    }

    if (x == 0) {
        if (nu == 0)
            return 0.0;
        double roundNu = std::round(nu);
        bool nuIsInt = areClose(nu, roundNu);
        return (nu > 0 || nuIsInt) ? -INFINITY : INFINITY;
    }

    if (std::fabs(nu) == 0.5) {
        /// log(sinh(x)) or log(cosh(x))
        long double y = 0.5 * (M_LN2 - M_LNPI - std::log(x));
        y -= M_LN2;
        y += x;
        y += (nu > 0) ? RandMath::log1pexp(-2 * x) : RandMath::log1mexp(-2 * x);
        return y;
    }

    if (nu < 0) {
        /// I(−ν, x) = I(ν, x) + 2 / π sin(πν) K(ν, x)
        long double besseli = std::cyl_bessel_il(-nu, x);
        long double sinPiNu = std::sin(M_PI * nu);
        long double y = (sinPiNu == 0) ? besseli : besseli - M_2_PI * sinPiNu * std::cyl_bessel_kl(-nu, x);
        return std::log(y);
    }

    long double besseli = std::cyl_bessel_il(nu, x);
    return std::isfinite(besseli) ? std::log(besseli) : x - 0.5 * (M_LN2 + M_LNPI + std::log(x));
}

long double logBesselK(double nu, double x)
{
    if (nu < 0.0)
        return NAN; /// K(-ν, x) = -K(ν, x) < 0

    if (x == 0.0)
        return INFINITY;

    long double besselk = 0;
    if (nu == 0.5 || (besselk = std::cyl_bessel_kl(nu, x)) == 0)
        return 0.5 * (M_LNPI - M_LN2 - std::log(x)) - x;

    if (!std::isfinite(besselk))
        return (nu == 0) ? std::log(-std::log(x)) : std::lgamma(nu) - M_LN2 - nu * std::log(0.5 * x);

    return std::log(besselk);
}

/**
 * @fn WLambert
 * @param x
 * @param w0
 * @param epsilon
 * @return
 */
double WLambert(double x, double w0, double epsilon)
{
    double w = w0;
    double step = 0;
    do {
        double ew = std::exp(w);
        double wew = w * ew;
        double numerator1 = wew - x;
        double wp1 = w + 1;
        double denominator1 = ew * wp1;
        double numerator2 = (w + 2) * numerator1;
        double denominator2 = 2 * wp1;
        step = numerator2 / denominator2;
        step = numerator1 / (denominator1 - step);
        w -= step;
    } while (std::fabs(step) > epsilon);
    return w;
}

double W0Lambert(double x, double epsilon)
{
    double w = 0;
    if (x < -M_1_E)
        return NAN;
    if (x > 10) {
        double logX = std::log(x);
        double loglogX = std::log(logX);
        w = logX - loglogX;
    }
    return WLambert(x, w, epsilon);
}

double Wm1Lambert(double x, double epsilon)
{
    double w = -2;
    if (x < -M_1_E || x > 0)
        return NAN;
    if (x > -0.1) {
        double logmX = std::log(-x);
        double logmlogmX = std::log(-logmX);
        w = logmX - logmlogmX;
    }
    return WLambert(x, w, epsilon);
}

/**
 * @fn MarcumPSeries
 * @param mu
 * @param x
 * @param y
 * @param logX log(x)
 * @param logY log(y)
 * @return series expansion for Marcum-P function
 */
double MarcumPSeries(double mu, double x, double y, double logX, double logY)
{
    /// ~log(2πε) for ε = 1e-16
    static constexpr double ln2piEps = -35.0;
    double lgammamu = std::lgamma(mu);
    double C = lgammamu - ln2piEps + mu;

    /// solving equation f(n) = 0
    /// to find first negleted term
    double root = std::max(0.5 * (mu * mu + 4 * x * y - mu), 1.0);
    double logXY = logX + logY;
    if (!RandMath::findRoot([C, mu, logXY] (double n)
    {
        double npmu = n + mu;
        double logn = std::log(n), lognpmu = std::log(npmu);
        double first = logn - 2 - logXY;
        first *= n;
        first += npmu * lognpmu;
        first -= C;
        double second = logn + lognpmu - logXY;
        double third = 1.0 / n + 1.0 / npmu;
        return DoubleTriplet(first, second, third);
    }, root))
        /// unexpected return
        return NAN;

    /// series expansion
    double sum = 0.0;
    /// sanity check
    int n0 = std::max(std::ceil(root), 5.0);
    double mpn0 = mu + n0;
    double P = pgamma(mpn0, y, logY);
    double diffP = (mpn0 - 1) * logY - y - std::lgamma(mpn0);
    diffP = std::exp(diffP);
    for (int n = n0; n > 0; --n) {
        double term = n * logX - x - lfact(n);
        double mupnm1 = mu + n - 1;
        term = std::exp(term) * P;
        sum += term;
        if (n % 5 == 0) {
            /// every 5 iterations we recalculate P and diffP
            /// in order to achieve enough accuracy
            P = pgamma(mupnm1, y, logY);
            diffP = (mupnm1 - 1) * logY - y - std::lgamma(mupnm1);
            diffP = std::exp(diffP);
        }
        else {
            /// otherwise we use recurrent relations
            P += diffP;
            diffP *= mupnm1 / y;
        }
    }
    /// add the last 0-term
    double lastTerm = lpgamma(mu, y, logY) - x;
    lastTerm = std::exp(lastTerm);
    sum += lastTerm;
    return sum;
}

/**
 * @fn MarcumPAsymptoticForLargeXY
 * @param mu
 * @param x
 * @param y
 * @param sqrtX √x
 * @param sqrtY √y
 * @return asymptotic expansion for Marcum-P function for large x*y
 */
double MarcumPAsymptoticForLargeXY(double mu, double x, double y, double sqrtX, double sqrtY)
{
    double xi = 2 * sqrtX * sqrtY;
    double sigma = (y + x) / xi - 1;
    double sigmaXi = sigma * xi;
    double rho = sqrtY / sqrtX;
    double aux = std::erfc(sqrtX - sqrtY);
    double Phi = (sigma == 0) ? 0.0 : sign(x - y) * std::sqrt(M_PI / sigma) * aux;
    double Psi0 = aux / std::sqrt(rho);
    double logXi = M_LN2 + 0.5 * std::log(x * y);
    /// sanity check
    int n0 = std::max(std::ceil(sigmaXi), 7.0);
    double sum = 0.0;
    double A1 = 1.0, A2 = 1.0;
    for (int n = 1; n <= n0; ++n) {
        /// change φ
        double nmHalf = n - 0.5;
        Phi *= -sigma;
        double temp = -sigmaXi - nmHalf * logXi;
        Phi += std::exp(temp);
        Phi /= nmHalf;
        /// calculate A(μ-1) and A(μ)
        double coef1 = (nmHalf - mu) * (nmHalf + mu);
        double coef2 = (nmHalf - mu + 1) * (nmHalf + mu - 1);
        double denom = -2 * n;
        A1 *= coef1 / denom;
        A2 *= coef2 / denom;
        /// compute term ψ and add it to the sum
        double Psi = Phi * (A2 - A1 / rho);
        sum += (n & 1) ? Psi : -Psi;
    }
    sum /= M_SQRT2PI;
    sum += Psi0;
    sum *= 0.5 * std::pow(rho, mu);
    return sum;
}

/**
 * @fn MarcumPForMuLessThanOne
 * @param mu
 * @param x
 * @param y
 * @return
 */
double MarcumPForMuLessThanOne(double mu, double x, double y, double logX, double logY)
{
    // TODO: check Krishnamoorthy paper for alternative representation

    /// in this case we use numerical integration,
    /// however we have singularity point at 0,
    /// so we get rid of it by subtracting the function
    /// which has the same behaviour at this point
    double aux = x + mu * M_LN2 + std::lgamma(mu);
    double log2x = M_LN2 + logX;
    double I = (M_LN2 + logY) * mu - aux;
    I = std::exp(I) / mu;

    double mum1 = mu - 1.0;
    I += RandMath::integral([x, log2x, mum1, aux] (double t)
    {
        if (t <= 0)
            return 0.0;
        /// Calculate log of leveling factor
        double log2T = std::log(2 * t);
        double exponent = mum1 * log2T;
        exponent -= aux;

        /// Calculate log(2f(t))
        double logBessel = RandMath::logBesselI(mum1, 2 * std::sqrt(x * t));
        double z = mum1 * (log2T - log2x);
        double log2F = 0.5 * z - t - x + logBessel;

        /// Return difference f(t) - factor
        return std::exp(log2F) - 2 * std::exp(exponent);
    }, 0, y);
    return I;
}

double MarcumQIntegrand(double theta, double xi, double sqrt1pXiSq, double mu, double y)
{
    double sinTheta = std::sin(theta);
    double theta_sinTheta = theta / sinTheta;
    double rho = std::hypot(theta_sinTheta, xi);
    double theta_sinThetapRho = theta_sinTheta + rho;
    double r = 0.5 * theta_sinThetapRho / y;
    double cosTheta = std::cos(theta);
    double psi = cosTheta * rho - sqrt1pXiSq;
    double frac = theta_sinThetapRho / (1.0 + sqrt1pXiSq);
    psi -= std::log(frac);
    double numerator = (sinTheta - theta * cosTheta) / (sinTheta * rho);
    numerator += cosTheta - r;
    numerator *= r;
    double denominator = r * r - 2.0 * r * cosTheta + 1.0;
    double f = numerator / denominator;
    return std::exp(mu * psi) * f;
}

double MarcumQIntergralRepresentation(double mu, double x, double y, double sqrtX, double sqrtY)
{
    double xi = 2.0 * sqrtX * sqrtY / mu;
    double sqrt1pXiSq = std::sqrt(1.0 + xi * xi);
    double yPrime = y / mu;
    double s0 = 0.5 * (1.0 + sqrt1pXiSq) / yPrime;
    double phi = x / s0 + y * s0 - std::log(s0) * mu;
    std::function<double (double)> integrandPtr = std::bind(&MarcumQIntegrand, std::placeholders::_1, xi, sqrt1pXiSq, mu, yPrime);
    double integral = RandMath::integral(integrandPtr, -M_PI, M_PI);
    return 0.5 * std::exp(-x - y + phi) / M_PI * integral;
}

double MarcumP(double mu, double x, double y, double sqrtX, double sqrtY, double logX, double logY)
{
    /* 1 - ok
     * 2 - ok
     * 3 - no
     * 4 - no
     * 5 - not yet
     * */
    if (x < 0.0 || y <= 0.0)
        return 0.0;

    if (x < 30)
        return MarcumPSeries(mu, x, y, logX, logY);

    double xi = 2 * sqrtX * sqrtY;
    if (xi > 30 && mu * mu < 2 * xi)
        return MarcumPAsymptoticForLargeXY(mu, x, y, sqrtX, sqrtY);

    double temp = std::sqrt(4 * x + 2 * mu);
    double f1 = x + mu - temp, f2 = x + mu + temp;
    if (y > f1 && y < f2)
        // IF mu > 135
        return 1.0 - MarcumQIntergralRepresentation(mu, x, y, sqrtX, sqrtY);

    // TODO: implement the rest techniques

    double mum1 = mu - 1;
    return RandMath::integral([mum1, logX, x](double t) {
        if (t < 0.0)
            return 0.0;
        if (t == 0.0)
            return (mum1 == 0) ? 0.5 * std::exp(-x) : 0.0;
        double logBesseli = RandMath::logBesselI(mum1, 2 * std::sqrt(x * t));
        double z = 0.5 * mum1 * (std::log(t) - logX);
        double h = t + x;
        return std::exp(logBesseli + z - h);
    }, 0, y);
}

double MarcumP(double mu, double x, double y)
{
    double sqrtX = std::sqrt(x), sqrtY = std::sqrt(y);
    double logX = std::log(x), logY = std::log(y);
    return MarcumP(mu, x, y, sqrtX, sqrtY, logX, logY);
}

double MarcumQ(double mu, double x, double y, double sqrtX, double sqrtY, double logX, double logY)
{
    // TODO: implement and use, when mu + x > y
    return 1.0 - MarcumP(mu, x, y, sqrtX, sqrtY, logX, logY);
}

double MarcumQ(double mu, double x, double y)
{
    double sqrtX = std::sqrt(x), sqrtY = std::sqrt(y);
    double logX = std::log(x), logY = std::log(y);
    return MarcumQ(mu, x, y, sqrtX, sqrtY, logX, logY);
}
}

