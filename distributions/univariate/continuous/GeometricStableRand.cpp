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
    if (exponent < 0.1 || exponent > 2.0)
        throw std::invalid_argument("Geometric-Stable distribution: exponent should be in the interval [0.1, 2]");
    if (std::fabs(skewness) > 1.0)
        throw std::invalid_argument("Geometric-Stable distribution: skewness should be in the interval [-1, 1]");

    alpha = exponent;
    alphaInv = 1.0 / alpha;
    beta = skewness;
    Z.SetParameters(alpha, beta);

    SetScale(scale);
    SetLocation(location);
    SetShift(shift);

    if (alpha == 2.0)
        distributionType = (mu == 0.0) ? LAPLACE : ASYMMETRIC_LAPLACE;
    else if (alpha == 0.5)
        distributionType = (std::fabs(beta) == 1.0) ? LEVY : ONEHALF_EXPONENT;
    else if (alpha == 1.0)
        distributionType = (beta == 0.0) ? CAUCHY : UNITY_EXPONENT;
    else
        distributionType = GENERAL;
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
    if (scale <= 0.0)
        throw std::invalid_argument("Geometric-Stable distribution: scale should be positive");
    gamma = scale;
    logGamma = std::log(gamma);
}

void ShiftedGeometricStableDistribution::SetAsymmetry(double asymmetry)
{
    if (asymmetry <= 0.0)
        throw std::invalid_argument("Asymmetric-Laplace distribution: asymmetry parameter should be positive");
    kappa = asymmetry;
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
        if (x < 0)
            return 0.0;
        if (x == 0)
            return INFINITY;
        double halfSigmaInv = 0.5 / gamma;
        double z = std::sqrt(halfSigmaInv * x);
        double y = M_1_SQRTPI - RandMath::xexpxsqerfc(z);
        return 0.5 * y / (gamma * z);
    }

    if (x == 0) {
        return std::fabs(1.0 / mu);
    }

    double a = x / gamma, c = mu / gamma;
    if (beta < 0) {
        a = -a;
        c = -c;
    }

    if (a < 0) {
        if (c >= 0.0)
            return 0.0;
        double h = std::sqrt(1 - 2 * c); // can be hashed
        double hm1 = h - 1;
        double y = 2 * a / (hm1 * hm1);
        y -= std::log(gamma * hm1 * h);
        y += M_LN2;
        return std::exp(y);
    }

    double sqrt2a = std::sqrt(2 * a);

    if (c < 0.5) {
        double h = std::sqrt(1 - 2 * c); // can be hashed
        double term1 = RandMath::sign(c) * sqrt2a / (1 - h);
        double term2 = sqrt2a / (1 + h);
        double y = RandMath::xexpxsqerfc(term1);
        y -= RandMath::xexpxsqerfc(term2);
        y /= h * gamma * sqrt2a;
        return y;
    }

    if (c == 0.5) {
        double y = RandMath::xexpxsqerfc(sqrt2a);
        y *= (8.0 * a + 2.0) / sqrt2a;
        y -= 4 * sqrt2a / M_SQRTPI;
        return y / gamma;
    }

    /// we do numerical integration in the case of c > 0.5
    // (as I wasn't able to handle complex integrals)
    double adivc = x / mu;
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
    double x0 = x - m;
    /// Cut zeros for half-infinite distribution
    if (alpha < 1 && std::fabs(beta) == 1) {
        if (x0 < 0 && mu >= 0 && beta > 0)
            return 0.0;
        if (x0 > 0 && mu <= 0 && beta < 0)
            return 0.0;
    }

    /// Check if analytical expression is known
    if (distributionType == LAPLACE || distributionType == ASYMMETRIC_LAPLACE)
        return pdfLaplace(x0);
    if (distributionType == CAUCHY)
        return pdfByCauchy(x0);
    if (distributionType == LEVY)
        return pdfByLevy(x0);

    if (x0 == 0) {
        if (alpha == 1)
            return INFINITY;
        if (mu == 0) {
            double coef = std::tgamma(1.0 - alphaInv);
            if (!std::isfinite(coef))
                return INFINITY;
            return Z.f(0) * coef / gamma;
        }
    }

    return RandMath::integral([this, x0] (double z)
    {
        if (z <= 0 || z >= 1)
            return 0.0;
        double logz = std::log(z);
        double tau = 1.0 / (gamma * std::pow(-logz, alphaInv));
        double par = tau * (x0 + mu * logz);
        return tau * Z.f(par);
    },
    0, 1);
}

double ShiftedGeometricStableDistribution::logf(const double & x) const
{
    return (distributionType == LAPLACE || distributionType == ASYMMETRIC_LAPLACE) ? logpdfLaplace(x - m) : std::log(f(x));
}

double ShiftedGeometricStableDistribution::F(const double & x) const
{
    double x0 = x - m;
    if (distributionType == LAPLACE || distributionType == ASYMMETRIC_LAPLACE)
        return cdfLaplace(x0);
    return RandMath::integral([this, x0] (double z)
    {
        if (z <= 0) {
            if (x0 != 0.0)
                return (x0 < 0.0) ? 0.0 : 1.0;
            if (mu == 0)
                return Z.F(0);
            if (distributionType == UNITY_EXPONENT)
                return Z.F(-mu / gamma);
            if (alpha < 1.0)
                return Z.F(0);
            return (mu < 0.0) ? 1.0 : 0.0;
        }

        if (z >= 1) {
            if (mu != 0.0) {
                if (distributionType == UNITY_EXPONENT)
                    return Z.F(-mu / gamma);
                if (alpha > 1.0)
                    return Z.F(0);
                return (mu < 0.0) ? 1.0 : 0.0;
            }
            return Z.F(0);
        }

        double denominator = 1.0 - z;
        double t = z / denominator;
        double par = x0 - mu * t;
        par /= std::pow(t, alphaInv) * gamma;
        double y = std::exp(-t) * Z.F(par);
        return y / (denominator * denominator);
    },
    0, 1);
}

double ShiftedGeometricStableDistribution::variateForUnityExponent(double z) const
{
    double W = ExponentialRand::StandardVariate();
    double X = z + 2 * beta * std::log(gamma * W) / M_PI;
    X *= gamma;
    X += mu;
    return X * W;
}

double ShiftedGeometricStableDistribution::variateForGeneralExponent(double z) const
{
    double W = ExponentialRand::StandardVariate();
    double X = std::pow(W, alphaInv) * gamma * z;
    return mu * W + X;
}

double ShiftedGeometricStableDistribution::variateForOneHalfExponent(double z) const
{
    double W = ExponentialRand::StandardVariate();
    double X = mu + gamma * W * z;
    return X * W;
}

double ShiftedGeometricStableDistribution::variateByCauchy(double z) const
{
    double W = ExponentialRand::StandardVariate();
    double X = mu + gamma * z;
    return X * W;
}

double ShiftedGeometricStableDistribution::Variate() const
{
    switch (distributionType) {
    case LAPLACE:
        return gamma * LaplaceRand::StandardVariate();
    case ASYMMETRIC_LAPLACE:
        return gamma * AsymmetricLaplaceRand::StandardVariate(kappa);
    case ONEHALF_EXPONENT:
    case LEVY:
        return variateForOneHalfExponent(Z.Variate());
    case CAUCHY:
        return variateByCauchy(Z.Variate());
    case UNITY_EXPONENT:
        return variateForUnityExponent(Z.Variate());
    case GENERAL:
        return variateForGeneralExponent(Z.Variate());
    }
    return NAN;
}

void ShiftedGeometricStableDistribution::Sample(std::vector<double> &outputData) const
{
    switch (distributionType) {
    case LAPLACE: {
        for (double & var : outputData)
            var = gamma * LaplaceRand::StandardVariate();
    }
        break;
    case ASYMMETRIC_LAPLACE: {
        for (double & var : outputData)
            var = gamma * AsymmetricLaplaceRand::StandardVariate(kappa);
    }
        break;
    case ONEHALF_EXPONENT:
    case LEVY: {
        Z.Sample(outputData);
        for (double & var : outputData)
            var = variateForOneHalfExponent(var);
    }
        break;
    case CAUCHY: {
        Z.Sample(outputData);
        for (double & var : outputData)
            var = variateByCauchy(var);
    }
        break;
    case UNITY_EXPONENT: {
        Z.Sample(outputData);
        for (double & var : outputData)
            var = variateForUnityExponent(var);
    }
        break;
    case GENERAL:
    default: {
        Z.Sample(outputData);
        for (double & var : outputData)
            var = variateForGeneralExponent(var);
    }
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
    if (distributionType == LAPLACE || distributionType == ASYMMETRIC_LAPLACE) {
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
    if (distributionType == LAPLACE || distributionType == ASYMMETRIC_LAPLACE) {
        return quantileLaplace(0.5);
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
    if (distributionType == LAPLACE || distributionType == ASYMMETRIC_LAPLACE) {
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
    if (distributionType == LAPLACE || distributionType == ASYMMETRIC_LAPLACE)
    {
        double denominator = kappaSq + kappaInv * kappaInv;
        denominator *= denominator;
        return 6.0 - 12.0 / denominator;
    }
    return NAN;
}

double ShiftedGeometricStableDistribution::quantileLaplace(double p) const
{
    if (p < kappaSq / (1 + kappaSq)) {
        double q = p * (1.0 / kappaSq + 1.0);
        q = std::log(q);
        q *= kappa * gamma;
        return m + q;
    }
    else {
        double q = (kappaSq + 1) * (1 - p);
        q = std::log(q);
        q *= gamma / kappa;
        return m - q;
    }
}

double ShiftedGeometricStableDistribution::quantileLaplace1m(double p) const
{
    if (p > 1.0 / (1 + kappaSq)) {
        double q = (1.0 - p) * (1.0 / kappaSq + 1.0);
        q = std::log(q);
        q *= kappa * gamma;
        return m + q;
    }
    else {
        double q = (kappaSq + 1) * p;
        q = std::log(q);
        q *= gamma / kappa;
        return m - q;
    }
}

GeometricStableRand::GeometricStableRand(double exponent, double skewness, double scale, double location)
    : ShiftedGeometricStableDistribution(exponent, skewness, scale, location)
{
    ChangeAsymmetry();
}

String GeometricStableRand::Name() const
{
    return "Geometric Stable("
            + toStringWithPrecision(GetExponent()) + ", "
            + toStringWithPrecision(GetSkewness()) + ", "
            + toStringWithPrecision(GetScale()) + ", "
            + toStringWithPrecision(GetLocation()) + ")";
}

void GeometricStableRand::ChangeAsymmetry()
{
    if (distributionType != LAPLACE && distributionType != ASYMMETRIC_LAPLACE)
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
