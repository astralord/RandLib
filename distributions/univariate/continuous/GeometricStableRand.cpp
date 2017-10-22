#include "GeometricStableRand.h"
#include "LaplaceRand.h"
#include "LevyRand.h"
#include "UniformRand.h"
#include "CauchyRand.h"

ShiftedGeometricStableDistribution::ShiftedGeometricStableDistribution(double exponent, double skewness, double scale, double location, double shift)
{
    SetParameters(exponent, skewness, scale, location, shift);
}

void ShiftedGeometricStableDistribution::SetParameters(double exponent, double skewness, double scale, double location, double shift)
{
    Z.SetParameters(exponent, skewness);
    alpha = Z.GetExponent();
    beta = Z.GetSkewness();
    gamma = (scale > 0.0) ? scale : M_SQRT2;
    logGamma = std::log(gamma);
    mu = location;
    m = shift;
}

void ShiftedGeometricStableDistribution::SetLocation(double location)
{
    mu = location;
}

void ShiftedGeometricStableDistribution::SetShift(double shift)
{
    m = shift;
}

void ShiftedGeometricStableDistribution::SetScale(double scale)
{
    gamma = (scale > 0.0) ? scale : M_SQRT2;
    logGamma = std::log(gamma);
}

void ShiftedGeometricStableDistribution::SetAsymmetry(double asymmetry)
{
    kappa = asymmetry;
    if (kappa <= 0)
        kappa = 1.0;
    kappaInv = 1.0 / kappa;
    kappaSq = kappa * kappa;
    log1pKappaSq = std::log1p(kappaSq);
    double logK = std::log(kappa);
    pdfCoef = logGamma + log1pKappaSq - logK;
    cdfCoef = 2 * logK - log1pKappaSq;
}

SUPPORT_TYPE ShiftedGeometricStableDistribution::SupportType() const
{
    if (alpha < 1) {
        if (beta == 1 && mu >= 0)
            return RIGHTSEMIFINITE_T;
        if (beta == -1 && mu <= 0)
            return LEFTSEMIFINITE_T;
    }
    return INFINITE_T;
}

double ShiftedGeometricStableDistribution::MinValue() const
{
    if (alpha < 1 && beta == 1 && mu >= 0)
        return m;
    return -INFINITY;
}

double ShiftedGeometricStableDistribution::MaxValue() const
{
    if (alpha < 1 && beta == -1 && mu <= 0)
        return m;
    return INFINITY;
}

double ShiftedGeometricStableDistribution::pdfLaplace(double x) const
{
    return std::exp(logpdfLaplace(x));
}

double ShiftedGeometricStableDistribution::logpdfLaplace(double x) const
{
    double y = x / gamma;
    y *= (x < 0) ? kappaInv : -kappa;
    return y - pdfCoef;
}

double ShiftedGeometricStableDistribution::cdfLaplace(double x) const
{
    double y = x / gamma;
    if (x < 0) {
        y *= kappaInv;
        y += cdfCoef;
        return std::exp(y);
    }
    return -std::expm1(-log1pKappaSq - kappa * y);
}

double ShiftedGeometricStableDistribution::cdfLaplaceCompl(double x) const
{
    double y = x / gamma;
    if (x < 0) {
        y *= kappaInv;
        y += cdfCoef;
        return -std::expm1(y);
    }
    return std::exp(-log1pKappaSq - kappa * y);
}

double ShiftedGeometricStableDistribution::pdfByLevy(double x) const
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
            return -2.0 / (gamma * h * h1) * std::exp(r * (h0 - c));
        double K = h1 - c;
        K = std::exp(r * K);
        K /= 2 * gamma * c * h;
        double p =  std::exp(r * h);
        double temp = std::sqrt(0.5 * a) / c;
        double g = h0 * temp;
        double q = h1 * temp;
        double y = p;
        /// if c is in (0, 0.5)
        if (c > 0) {
            /// g can be too small while p is too big
            y *= std::erfc(g);
        }
        else {
            y *= -std::erfc(-g);
        }
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

double ShiftedGeometricStableDistribution::pdfByCauchy(double x) const
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

double ShiftedGeometricStableDistribution::f(const double & x) const
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

double ShiftedGeometricStableDistribution::logf(const double & x) const
{
    return (alpha == 2) ? logpdfLaplace(x) : std::log(f(x));
}

double ShiftedGeometricStableDistribution::F(const double & x) const
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

double ShiftedGeometricStableDistribution::variateForUnityExponent() const
{
    double W = ExponentialRand::StandardVariate();
    double X = Z.Variate();
    return W * (mu + gamma * X + 2 * gamma * std::log(gamma * W) / M_PI);
}

double ShiftedGeometricStableDistribution::variateForCommonExponent() const
{
    double W = ExponentialRand::StandardVariate();
    double X = Z.Variate();
    return mu * W + std::pow(W, alphaInv) * gamma * X;
}

double ShiftedGeometricStableDistribution::variateByLevy(bool positive) const
{
    double W = ExponentialRand::StandardVariate();
    double X = LevyRand::StandardVariate();
    if (!positive)
        X = -X;
    X *= gamma * W;
    X += mu;
    return X * W;
}

double ShiftedGeometricStableDistribution::variateByCauchy() const
{
    double W = ExponentialRand::StandardVariate();
    double X = mu + gamma * CauchyRand::StandardVariate();
    return X * W;
}

double ShiftedGeometricStableDistribution::Variate() const
{
    if (alpha == 2) {
        double X = (kappa == 1.0) ? LaplaceRand::StandardVariate() : AsymmetricLaplaceRand::StandardVariate(kappa);
        return gamma * X;
    }
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

void ShiftedGeometricStableDistribution::Sample(std::vector<double> &outputData) const
{
    if (alpha == 2) {
        if (kappa == 1.0) {
            for (double &var : outputData)
                var = gamma * LaplaceRand::StandardVariate();
        }
        else {
            for (double &var : outputData)
                var = gamma * AsymmetricLaplaceRand::StandardVariate(kappa);
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

double ShiftedGeometricStableDistribution::Mean() const
{
    if (alpha > 1)
        return m + mu;
    if (beta == 1)
        return INFINITY;
    return (beta == -1) ? -INFINITY : NAN;
}

double ShiftedGeometricStableDistribution::Variance() const
{
    if (alpha == 2) {
        return mu * mu + 2 * gamma * gamma;
    }
    return INFINITY;
}

std::complex<double> ShiftedGeometricStableDistribution::CFImpl(double t) const
{
    double x = 0;
    if (alpha != 2 && beta != 0) {
        x = (alpha == 1) ? M_2_PI * std::log(t) : std::tan(M_PI_2 * alpha);
        x *= beta;
    }
    double re = std::pow(gamma * t, alpha);
    std::complex<double> psi = std::complex<double>(re, re * x - mu * t);
    return 1.0 / (1.0 + psi);
}

double ShiftedGeometricStableDistribution::Median() const
{
    if (alpha == 2) {
        double y = std::log1p(1.0 / kappaSq) - M_LN2;
        return m + gamma * kappa * y;
    }
    return ContinuousDistribution::Median();
}

double ShiftedGeometricStableDistribution::Mode() const
{
    // TODO: calculate mode for more cases
    if (alpha == 1 || alpha == 2)
        return m;
    return ContinuousDistribution::Mode();
}

double ShiftedGeometricStableDistribution::Skewness() const
{
    if (alpha == 2) {
        double k4 = kappaSq * kappaSq;
        double k6 = k4 * kappaSq;
        double z = (k4 + 1);
        z *= std::sqrt(z);
        double y = 2 * (1 - k6);
        return y / z;
    }
    return NAN;
}

double ShiftedGeometricStableDistribution::ExcessKurtosis() const
{
    if (alpha == 2)
    {
        double denominator = kappaSq + kappaInv * kappaInv;
        denominator *= denominator;
        return 6.0 - 12.0 / denominator;
    }
    return NAN;
}

std::string GeometricStableRand::Name() const
{
    return "Geometric Stable("
            + toStringWithPrecision(GetExponent()) + ", "
            + toStringWithPrecision(GetSkewness()) + ", "
            + toStringWithPrecision(GetScale()) + ", "
            + toStringWithPrecision(GetLocation()) + ")";
}

void GeometricStableRand::ChangeAsymmetry()
{
    if (alpha != 2)
        return;
    double gamma2 = 2 * gamma;
    double asymmetry = mu * mu + gamma2 * gamma2;
    asymmetry = std::sqrt(asymmetry) - mu;
    asymmetry /= gamma2;
    SetAsymmetry(asymmetry);
}

void GeometricStableRand::SetParameters(double exponent, double skewness)
{
    ShiftedGeometricStableDistribution::SetParameters(exponent, skewness);
    ChangeAsymmetry();
}

void GeometricStableRand::SetLocation(double location)
{
    ShiftedGeometricStableDistribution::SetLocation(location);
    ChangeAsymmetry();
}

void GeometricStableRand::SetScale(double scale)
{
    ShiftedGeometricStableDistribution::SetScale(scale);
    ChangeAsymmetry();
}
