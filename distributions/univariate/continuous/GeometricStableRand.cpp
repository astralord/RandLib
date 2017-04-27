#include "GeometricStableRand.h"
#include "LaplaceRand.h"
#include "LevyRand.h"
#include "UniformRand.h"
#include "CauchyRand.h"

GeometricStableRand::GeometricStableRand(double exponent, double skewness, double scale, double location)
    : LimitingDistribution(exponent, skewness, scale, location),
      Z(exponent, skewness)
{
    SetParameters(exponent, skewness, scale, location);
}

std::string GeometricStableRand::Name() const
{
    return "Geometric Stable("
            + toStringWithPrecision(GetExponent()) + ", "
            + toStringWithPrecision(GetSkewness()) + ", "
            + toStringWithPrecision(GetScale()) + ", "
            + toStringWithPrecision(GetLocation()) + ")";
}

void GeometricStableRand::SetParameters(double exponent, double skewness, double scale, double location)
{
    LimitingDistribution::SetParameters(exponent, skewness);
    LimitingDistribution::SetScale(scale);
    LimitingDistribution::SetLocation(location);

    Z.SetParameters(alpha, beta);

    if (alpha == 2) {
        double sigma2 = gamma + gamma;
        k = mu * mu + sigma2 * sigma2;
        k = std::sqrt(k) - mu;
        k /= sigma2;
        kInv = 1.0 / k;
        kSq = k * k;
        pdfCoef = std::log(gamma * (k + kInv));
        cdfCoef = -std::log1p(kSq);
    }
}

double GeometricStableRand::pdfLaplace(double x) const
{
    return std::exp(logpdfLaplace(x));
}

double GeometricStableRand::logpdfLaplace(double x) const
{
    double y = x / gamma;
    y *= (x < 0) ? kInv : -k;
    return y - pdfCoef;
}

double GeometricStableRand::cdfLaplace(double x) const
{
    double y = x / gamma;
    if (x < 0) {
        y *= kInv;
        y = std::exp(cdfCoef + y);
        return kSq * y;
    }
    return -std::expm1(cdfCoef - k * y);
}

double GeometricStableRand::cdfLaplaceCompl(double x) const
{
    double y = x / gamma;
    if (x < 0) {
        y *= kInv;
        y = std::exp(cdfCoef + y);
        return 1.0 - kSq * y;
    }
    return std::exp(cdfCoef - k * y);
}

double GeometricStableRand::pdfByLevy(double x) const
{
    if (mu == 0) {
        /// Invert parameter for case of -Levy
        x *= beta;
        double halfSigmaInv = 0.5 / gamma;
        double z = halfSigmaInv * x;
        double sqrtZ = std::sqrt(z);
        double y = -std::erfc(sqrtZ);
        y *= std::exp(z);
        y += M_1_SQRTPI / sqrtZ;
        y *= halfSigmaInv;
        return y;
    }

    /// If mu != 0
    double a = x / gamma, c = mu / gamma;
    if (beta < 0) {
        a = -a;
        c = -c;
    }
    if (a < 0 && c >= 0)
        return 0.0;
    if (c < 0.5) {
        double h = std::sqrt(1 - 2 * c);
        double h0 = 1 + h;
        double h1 = 1 - h;
        double r = x / (mu * c);
        if (a < 0 && c < 0)
            return -2 / (gamma * h * h1) * std::exp(r * (h0 - c));
        double K = h1 - c;
        K = std::exp(r * K);
        K /= 2 * gamma * c * h;
        double p =  std::exp(r * h);
        double gamma = std::sqrt(0.5 * a) / c;
        double g = h0 * gamma;
        double q = h1 * gamma;
        double y = p;
        if (c > 0) /// c is in (0, 0.5)
            y *= std::erfc(g); /// g can be too small while p is too big
        else
            y *= -std::erfc(-g);
        y *= h0 * p;
        y -= h1 * std::erfc(q);
        y *= K;
        return y;
    }

    /// we do numerical integration in this case
    // (as I wasn't able to handle complex integrals)
    double adivc = a / c;
    double y = RandMath::integral([this, a, c, adivc] (double t)
    {
        if (t <= 0 || t >= adivc)
            return 0.0;
        double amct = a - c * t;
        double integrand = 0.5 * t / amct;
        ++integrand;
        integrand *= t;
        integrand = std::exp(-integrand);
        integrand *= t;
        integrand /= amct * std::sqrt(amct);
        return integrand;
    }, 0, adivc);
    return y / (gamma * M_SQRT2PI);
}

double GeometricStableRand::pdfByCauchy(double x) const
{
    // TODO: find analytical solution, maybe using residue technique
    // otherwise, use peak, which is t = exp(x / sqrt(mu^2 + sigma^2))
    // to split the integral on two
    if (x == 0)
        return INFINITY;
    return -gamma * M_1_PI * RandMath::integral([this, x] (double t)
    {
        if (t <= 0.0 || t >= 1.0)
            return 0.0;
        double logt = std::log(t);
        double a = gamma * logt;
        a *= a;
        double b = x + mu * logt;
        b *= b;
        return logt / (a + b);
    }, 0, 1);
}

double GeometricStableRand::f(const double & x) const
{
    /// Laplace case
    if (alpha == 2)
        return pdfLaplace(x);

    /// Cut zeros for half-infinite distribution
    if (alpha < 1 && std::fabs(beta) == 1) {
        if (x < 0 && mu >= 0 && beta > 0)
            return 0.0;
        if (x > 0 && mu <= 0 && beta < 0)
            return 0.0;
    }

    /// Cauchy case
    if (alpha == 1.0 && beta == 0.0)
        return pdfByCauchy(x);

    /// Levy case
    if (alpha == 0.5 && std::fabs(beta) == 1)
        return pdfByLevy(x);

    if (x == 0) {
        if (alpha == 1)
            return INFINITY;
        if (mu == 0) {
            double coef = std::tgamma(1.0 - alphaInv);
            if (!std::isfinite(coef))
                return INFINITY;
            return Z.f(0) * coef / gamma;
        }

        // can be hashed
        double c = std::lgamma(alpha);
        double signMu = RandMath::sign(mu);
        c += (alpha - 1) * std::log(gamma);
        c -= (alpha + 1) * std::log(signMu * mu / gamma);
        c = std::exp(c);
        c *= std::sin(M_PI_2 * alpha);
        c /= M_PI;
        c *= alpha * (1 - signMu * beta);

        double int1 = RandMath::integral([this, c] (double z)
        {
            if (z <= 0 || z >= 1)
                return 0.0;
            double denominator = 1.0 - z;
            double t = z / denominator;
            double temp = gamma * std::pow(t, alphaInv);
            double par = -mu * t / temp;
            double y = Z.f(par) / temp;
            y -= c / std::pow(t, alpha);
            y *= std::exp(-t);
            return y / (denominator * denominator);
        },
        0, 1);

        double int2 = c * std::tgamma(1 - alpha); // hash
        return int1 + int2;
    }

    return RandMath::integral([this, x] (double z)
    {
        if (z <= 0 || z >= 1)
            return 0.0;
        double logz = std::log(z);
        double tau = 1.0 / (gamma * std::pow(-logz, alphaInv));
        double par = tau * (x + mu * logz);
        return tau * Z.f(par);
    },
    0, 1);
}

double GeometricStableRand::logf(const double & x) const
{
    return std::log(f(x));
}

double GeometricStableRand::F(const double & x) const
{
    if (alpha == 2)
        return cdfLaplace(x);
    return RandMath::integral([this, x] (double z)
    {
        if (z <= 0)
            return 1.0;
        if (z >= 1) {
            if (mu == 0)
                return Z.F(0);
            if (alpha > 1)
                return (mu > 0) ? 0.0 : 1.0;
            double par = (alpha < 1) ? 0 : -mu / gamma;
            return Z.F(par);
        }
        double denominator = 1.0 - z;
        double t = z / denominator;
        double par = x - mu * t;
        par /= std::pow(t, alphaInv) * gamma;
        double y = std::exp(-t) * Z.F(par);
        return y / (denominator * denominator);
    },
    0, 1);
}

double GeometricStableRand::variateForUnityExponent() const
{
    double W = ExponentialRand::StandardVariate();
    double X = Z.Variate();
    return W * (mu + gamma * X + 2 * gamma * std::log(gamma * W) / M_PI);
}

double GeometricStableRand::variateForCommonExponent() const
{
    double W = ExponentialRand::StandardVariate();
    double X = Z.Variate();
    return mu * W + std::pow(W, alphaInv) * gamma * X;
}

double GeometricStableRand::variateByLevy(bool positive) const
{
    double W = ExponentialRand::StandardVariate();
    double X = LevyRand::StandardVariate();
    if (!positive)
        X = -X;
    X *= gamma * W;
    X += mu;
    return X * W;
}

double GeometricStableRand::variateByCauchy() const
{
    double W = ExponentialRand::StandardVariate();
    double X = CauchyRand::Variate(mu, gamma);
    return X * W;
}

double GeometricStableRand::Variate() const
{
    if (alpha == 2)
        return (mu == 0) ? LaplaceRand::Variate(0, gamma) : LaplaceRand::Variate(0, gamma, k);
    if (alpha == 0.5) {
        if (beta == 1)
            return variateByLevy(true);
        if (beta == -1)
            return variateByLevy(false);
    }
    if (alpha == 1)
        return (beta == 0) ? variateByCauchy() : variateForUnityExponent();
    return variateForCommonExponent();
}

void GeometricStableRand::Sample(std::vector<double> &outputData) const
{
    if (alpha == 2) {
        if (mu == 0) {
            for (double &var : outputData)
                var = LaplaceRand::Variate(0, gamma);
        }
        else {
            for (double &var : outputData)
                var = LaplaceRand::Variate(0, gamma, k);
        }
    }
    else if (alpha == 0.5 && std::fabs(beta) == 1) {
        if (beta > 0) {
            for (double & var : outputData)
                var = variateByLevy(true);
        }
        else {
            for (double & var : outputData)
                var = variateByLevy(false);
        }
    }
    else if (alpha == 1) {
        if (beta == 0) {
            for (double &var : outputData)
                var = variateByCauchy();
        }
        else {
            for (double &var : outputData)
                var = variateForUnityExponent();
        }
    }
    else {
        for (double &var : outputData)
            var = variateForCommonExponent();
    }
}

double GeometricStableRand::Variance() const
{
    if (alpha == 2) {
        double y = gamma * gamma / kSq;
        return (1.0 + kSq * kSq) * y;
    }
    return INFINITY;
}

std::complex<double> GeometricStableRand::CFImpl(double t) const
{
    return 1.0 / (1.0 + psi(t));
}

double GeometricStableRand::Median() const
{
    if (alpha == 2) {
        if (k > 1) {
            double y = std::log1p(1.0 / kSq) - M_LN2;
            return gamma * k * y;
        }
        else {
            double y = std::log1p(kSq) - M_LN2;
            return -gamma / k * y;
        }
    }
    return ContinuousDistribution::Median();
}

double GeometricStableRand::Mode() const
{
    // TODO: calculate mode for more cases
    if (alpha == 1 || alpha == 2)
        return 0.0;
    return ContinuousDistribution::Mode();
}

double GeometricStableRand::Skewness() const
{
    if (alpha == 2) {
        double kSq = k * k;
        double k4 = kSq * kSq;
        double k6 = k4 * kSq;
        double z = (k4 + 1);
        z *= std::sqrt(z);
        double y = 2 * (1 - k6);
        return y / z;
    }
    return NAN;
}

double GeometricStableRand::ExcessKurtosis() const
{
    if (alpha == 2)
    {
        if (k == 1)
            return 3.0;
        return NAN; // TODO!
    }
    return NAN;
}
