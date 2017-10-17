#include "LaplaceRand.h"
#include "../discrete/BernoulliRand.h"

LaplaceRand::LaplaceRand(double shift, double scale, double asymmetry)
    : ShiftedGeometricStableDistribution(2.0, 0.0, scale, 0.0, shift)
{
    SetAsymmetry(asymmetry);
}

std::string LaplaceRand::Name() const
{
    return "Laplace(" + toStringWithPrecision(GetShift()) + ", "
                      + toStringWithPrecision(GetScale()) + ", "
                      + toStringWithPrecision(GetAsymmetry()) + ")";
}

void LaplaceRand::ChangeLocation()
{
    SetLocation((1.0 - kappaSq) * gamma * kappaInv);
}

void LaplaceRand::SetScale(double scale)
{
    ShiftedGeometricStableDistribution::SetScale(scale);
    ChangeLocation();
}

void LaplaceRand::SetAsymmetry(double asymmetry)
{
    ShiftedGeometricStableDistribution::SetAsymmetry(asymmetry);
    ChangeLocation();
}

double LaplaceRand::f(const double & x) const
{
    return pdfLaplace(x - m);
}

double LaplaceRand::logf(const double & x) const
{
    return logpdfLaplace(x - m);
}

double LaplaceRand::F(const double & x) const
{
    return cdfLaplace(x - m);
}

double LaplaceRand::S(const double & x) const
{
    return cdfLaplaceCompl(x - m);
}

double LaplaceRand::Variate() const
{
    return (kappa == 1) ? m + gamma * LaplaceRand::StandardVariate() : LaplaceRand::Variate(m, gamma, kappa);
}

double LaplaceRand::StandardVariate()
{
    double W = ExponentialRand::StandardVariate();
    return BernoulliRand::StandardVariate() ? W : -W;
}

double LaplaceRand::Variate(double location, double scale, double asymmetry)
{
    // TODO: set sanity check here or make this function private
    double x = ExponentialRand::StandardVariate() / asymmetry;
    double y = ExponentialRand::StandardVariate() * asymmetry;;
    return location + scale * (x - y);
}

void LaplaceRand::Sample(std::vector<double> &outputData) const
{
    if (kappa == 1) {
        for (double & var : outputData)
            var = m + gamma * LaplaceRand::StandardVariate();
    }
    else {
        for (double & var : outputData)
            var = LaplaceRand::Variate(m, gamma, kappa);
    }
}

std::complex<double> LaplaceRand::CFImpl(double t) const
{
    double bt = gamma * t;
    double btSq = bt * bt;
    double denominator = (1 + kappaSq * btSq) * (1 + btSq / kappaSq);
    std::complex<double> y(std::cos(m * t), std::sin(m * t));
    std::complex<double> x(1, -kappa * bt), z(1, bt * kappaInv);
    return x * y * z / denominator;
}

double LaplaceRand::quantileImpl(double p) const
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

double LaplaceRand::quantileImpl1m(double p) const
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

double LaplaceRand::Entropy() const
{
    double y = kappaInv + kappa;
    return std::log1p(gamma * y);
}

void LaplaceRand::FitLocation(const std::vector<double> &sample)
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
        throw std::runtime_error(fitError(UNDEFINED_ERROR, "Error in root-finding procedure"));

    SetShift(median);
}

void LaplaceRand::FitScale(const std::vector<double> &sample)
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

void LaplaceRand::FitAsymmetry(const std::vector<double> &sample)
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
        throw std::runtime_error(fitError(UNDEFINED_ERROR, "Error in root-finding procedure"));

    SetAsymmetry(root);
}

void LaplaceRand::FitLocationAndScale(const std::vector<double> &sample)
{
    FitLocation(sample);
    FitScale(sample);
}

void LaplaceRand::FitLocationAndAsymmetry(const std::vector<double> &sample)
{
    FitLocation(sample);
    FitAsymmetry(sample);
}

void LaplaceRand::FitScaleAndAsymmetry(const std::vector<double> &sample)
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
        throw std::runtime_error(fitError(UNDEFINED_ERROR, "Distribution might be exponentially distributed"));
    }
    if (xPlus == 0) {
        /// -X ~ Exp(1 / xMinus)
        throw std::runtime_error(fitError(UNDEFINED_ERROR, "Distribution might be exponentially distributed"));
    }

    double xPlusSqrt = std::sqrt(xPlus), xMinusSqrt = std::sqrt(xMinus);
    double scale = xPlusSqrt + xMinusSqrt;
    scale *= std::sqrt(xPlusSqrt * xMinusSqrt);

    SetScale(scale);
    SetAsymmetry(std::pow(xMinus / xPlus, 0.25));
}

void LaplaceRand::Fit(const std::vector<double> &sample)
{
    FitLocation(sample);
    FitScaleAndAsymmetry(sample);
}
