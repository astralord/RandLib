#include "LaplaceRand.h"
#include "../discrete/BernoulliRand.h"

AsymmetricLaplaceDistribution::AsymmetricLaplaceDistribution(double shift, double scale, double asymmetry)
    : ShiftedGeometricStableDistribution(2.0, 0.0, scale, 0.0, shift)
{
    ShiftedGeometricStableDistribution::SetAsymmetry(asymmetry);
    ChangeLocation();
}

void AsymmetricLaplaceDistribution::ChangeLocation()
{
    SetLocation((1.0 - kappaSq) * gamma * kappaInv);
}

void AsymmetricLaplaceDistribution::SetScale(double scale)
{
    if (scale <= 0.0)
        throw std::invalid_argument("Laplace distribution: scale should be positive");
    ShiftedGeometricStableDistribution::SetScale(scale);
    ChangeLocation();
}

double AsymmetricLaplaceDistribution::f(const double & x) const
{
    return pdfLaplace(x - m);
}

double AsymmetricLaplaceDistribution::logf(const double & x) const
{
    return logpdfLaplace(x - m);
}

double AsymmetricLaplaceDistribution::F(const double & x) const
{
    return cdfLaplace(x - m);
}

double AsymmetricLaplaceDistribution::S(const double & x) const
{
    return cdfLaplaceCompl(x - m);
}

double AsymmetricLaplaceDistribution::Variate() const
{
    double X = (kappa == 1) ? LaplaceRand::StandardVariate() : AsymmetricLaplaceRand::StandardVariate(kappa);
    return m + gamma * X;
}

void AsymmetricLaplaceDistribution::Sample(std::vector<double> &outputData) const
{
    if (kappa == 1) {
        for (double & var : outputData)
            var = m + gamma * LaplaceRand::StandardVariate();
    }
    else {
        for (double & var : outputData)
            var = m + gamma * AsymmetricLaplaceRand::StandardVariate(kappa);
    }
}

std::complex<double> AsymmetricLaplaceDistribution::CFImpl(double t) const
{
    double bt = gamma * t;
    double btSq = bt * bt;
    double denominator = (1 + kappaSq * btSq) * (1 + btSq / kappaSq);
    std::complex<double> y(std::cos(m * t), std::sin(m * t));
    std::complex<double> x(1, -kappa * bt), z(1, bt * kappaInv);
    return x * y * z / denominator;
}

double AsymmetricLaplaceDistribution::quantileImpl(double p) const
{
    return quantileLaplace(p);
}

double AsymmetricLaplaceDistribution::quantileImpl1m(double p) const
{
    return quantileLaplace1m(p);
}

double AsymmetricLaplaceDistribution::Entropy() const
{
    double y = kappaInv + kappa;
    return std::log1p(gamma * y);
}

void AsymmetricLaplaceDistribution::FitShift(const std::vector<double> &sample)
{
    /// Calculate median (considering asymmetry)
    /// we use root-finding algorithm for median search
    double minVar = *std::min_element(sample.begin(), sample.end());
    double maxVar = *std::max_element(sample.begin(), sample.end());
    double median = 0.5 * (minVar + maxVar);

    if (!RandMath::findRoot([this, sample] (double med)
    {
        double y = 0.0;
        for (const double & x : sample) {
            if (x > med)
                y -= kappaSq;
            else if (x < med)
                ++y;
        }
        return y;
    },
    minVar, maxVar, median
    ))
        throw std::runtime_error(fitErrorDescription(UNDEFINED_ERROR, "Error in root-finding procedure"));

    SetShift(median);
}

void AsymmetricLaplaceDistribution::FitScale(const std::vector<double> &sample)
{
    double deviation = 0.0;
    for (const double & x : sample) {
        if (x > m)
            deviation += kappaSq * (x - m);
        else
            deviation -= (x - m);
    }
    deviation /= (kappa * sample.size());

    SetScale(deviation);
}

void AsymmetricLaplaceDistribution::FitShiftAndScale(const std::vector<double> &sample)
{
    FitShift(sample);
    FitScale(sample);
}

String AsymmetricLaplaceRand::Name() const
{
    return "Asymmetric-Laplace(" + toStringWithPrecision(GetShift()) + ", "
                                 + toStringWithPrecision(GetScale()) + ", "
            + toStringWithPrecision(GetAsymmetry()) + ")";
}

void AsymmetricLaplaceRand::SetAsymmetry(double asymmetry)
{
    ShiftedGeometricStableDistribution::SetAsymmetry(asymmetry);
    ChangeLocation();
}

double AsymmetricLaplaceRand::StandardVariate(double asymmetry)
{
    double x = ExponentialRand::StandardVariate() / asymmetry;
    double y = ExponentialRand::StandardVariate() * asymmetry;
    return x - y;
}

void AsymmetricLaplaceRand::FitAsymmetry(const std::vector<double> &sample)
{
    double xPlus = 0.0, xMinus = 0.0;
    for (const double & x : sample) {
        if (x < m)
            xMinus -= (x - m);
        else
            xPlus += (x - m);
    }

    if (xPlus == xMinus) {
        SetAsymmetry(1.0);
        return;
    }

    double gammaN = gamma * sample.size();
    double root = 1.0;
    double minBound, maxBound;
    if (xPlus < -xMinus) {
        minBound = 1.0;
        maxBound = std::sqrt(xMinus / xPlus);
    }
    else {
        minBound = std::sqrt(xMinus / xPlus);
        maxBound = 1.0;
    }

    if (!RandMath::findRoot([sample, xPlus, xMinus, gammaN] (double t)
    {
        double tSq = t * t;
        double y = 1.0 - tSq;
        y /= (t * (tSq + 1.0));
        y *= gammaN;
        y += xMinus / tSq - xPlus;
        return y;
    }, minBound, maxBound, root))
        throw std::runtime_error(fitErrorDescription(UNDEFINED_ERROR, "Error in root-finding procedure"));

    SetAsymmetry(root);
}

void AsymmetricLaplaceRand::FitShiftAndAsymmetry(const std::vector<double> &sample)
{
    FitShift(sample);
    FitAsymmetry(sample);
}

void AsymmetricLaplaceRand::FitScaleAndAsymmetry(const std::vector<double> &sample)
{
    int n = sample.size();
    double xPlus = 0.0, xMinus = 0.0;
    for (const double & x : sample) {
        if (x < m)
            xMinus -= (x - m);
        else
            xPlus += (x - m);
    }
    xPlus /= n;
    xMinus /= n;

    if (xMinus == 0) {
        /// X ~ Exp(1 / xPlus)
        throw std::runtime_error(fitErrorDescription(UNDEFINED_ERROR, "Distribution might be exponentially distributed"));
    }
    if (xPlus == 0) {
        /// -X ~ Exp(1 / xMinus)
        throw std::runtime_error(fitErrorDescription(UNDEFINED_ERROR, "Distribution might be exponentially distributed"));
    }

    double xPlusSqrt = std::sqrt(xPlus), xMinusSqrt = std::sqrt(xMinus);
    double scale = xPlusSqrt + xMinusSqrt;
    scale *= std::sqrt(xPlusSqrt * xMinusSqrt);

    SetScale(scale);
    SetAsymmetry(std::pow(xMinus / xPlus, 0.25));
}

void AsymmetricLaplaceRand::Fit(const std::vector<double> &sample)
{
    FitShift(sample);
    FitScaleAndAsymmetry(sample);
}


String LaplaceRand::Name() const
{
    return "Laplace(" + toStringWithPrecision(GetShift()) + ", "
                      + toStringWithPrecision(GetScale()) + ")";
}

double LaplaceRand::StandardVariate()
{
    double W = ExponentialRand::StandardVariate();
    return BernoulliRand::StandardVariate() ? W : -W;
}
