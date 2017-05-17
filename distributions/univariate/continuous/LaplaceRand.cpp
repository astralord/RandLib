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
    return (kappa == 1) ? LaplaceRand::Variate(m, gamma) : LaplaceRand::Variate(m, gamma, kappa);
}

double LaplaceRand::Variate(double location, double scale)
{
    bool sign = BernoulliRand::StandardVariate();
    double W = scale * ExponentialRand::StandardVariate();
    return location + (sign ? W : -W);
}

double LaplaceRand::Variate(double location, double scale, double asymmetry)
{
    double x = ExponentialRand::StandardVariate() / asymmetry;
    double y = ExponentialRand::StandardVariate() * asymmetry;;
    return location + scale * (x - y);
}

void LaplaceRand::Sample(std::vector<double> &outputData) const
{
    if (kappa == 1) {
        for (double & var : outputData)
            var = LaplaceRand::Variate(m, gamma);
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
    double y = 1 + kappaSq;
    y *= kappaInv * gamma;
    return log1p(y);
}

void LaplaceRand::FitLocationMLE(const std::vector<double> &sample)
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

void LaplaceRand::FitScaleMLE(const std::vector<double> &sample)
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

void LaplaceRand::FitAsymmetryMLE(const std::vector<double> &sample)
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

    double sigmaN = gamma * sample.size();
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

    if (!RandMath::findRoot([sample, xPlus, xMinus, sigmaN] (double t)
    {
        double tSq = t * t;
        double y = 1.0 - tSq;
        y /= (t * (tSq + 1.0));
        y *= sigmaN;
        y += xMinus / tSq - xPlus;
        return y;
    }, minBound, maxBound, root))
        throw std::runtime_error(fitError(UNDEFINED_ERROR, "Error in root-finding procedure"));

    SetAsymmetry(root);
}

void LaplaceRand::FitLocationAndScaleMLE(const std::vector<double> &sample)
{
    FitLocationMLE(sample);
    FitScaleMLE(sample);
}

void LaplaceRand::FitLocationAndAsymmetryMLE(const std::vector<double> &sample)
{
    FitLocationMLE(sample);
    FitAsymmetryMLE(sample);
}

void LaplaceRand::FitScaleAndAsymmetryMLE(const std::vector<double> &sample)
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

void LaplaceRand::FitLocationScaleAndAsymmetryMLE(const std::vector<double> &sample)
{
    FitLocationMLE(sample);
    FitScaleAndAsymmetryMLE(sample);
}

double LaplaceRand::GetAsymmetryFromSkewness(double skewness)
{
    double root = 1.0;
    if (!RandMath::findRoot([skewness] (double x)
    {
        /// f(x)
        double x2 = x * x, x3 = x2 * x, x4 = x2 * x2;
        double x4p1 = x4 + 1.0;
        double aux = x4p1 * std::sqrt(x4p1);
        double f1 = 2 * (1 - x4 * x2) / aux - skewness;
        /// f'(x)
        double y = x4p1 * aux;
        double f2 = -12 * x3 * (1 + x2) / y;
        return DoublePair(f1, f2);
    }, root))
        return NAN;
    return root;
}

void LaplaceRand::FitLocationMM(const std::vector<double> &sample)
{
    double mean = sampleMean(sample);
    SetShift(mean - gamma * (1.0 / kappa - kappa));
}

void LaplaceRand::FitScaleMM(const std::vector<double> &sample)
{
    if (kappa == 1) {
        /// can't derive scale from mean
        double var = sampleVariance(sample, m);
        SetScale(std::sqrt(0.5 * var));
    }
    else {
        double y = sampleMean(sample);
        SetScale((y - m) * kappa / (1.0 - kappaSq));
    }
}

void LaplaceRand::FitAsymmetryMM(const std::vector<double> &sample)
{
    double mean = sampleMean(sample);
    double aux = (m - mean) / gamma;
    double asymmetry = aux * aux + 4;
    asymmetry = std::sqrt(asymmetry);
    asymmetry += aux;
    SetAsymmetry(0.5 * asymmetry);
}

void LaplaceRand::FitLocationAndScaleMM(const std::vector<double> &sample)
{
    double mean = sampleMean(sample);
    double var = sampleVariance(sample, mean);
    double scale = (var * kappaSq) / (1 + kappaSq * kappaSq);
    SetScale(std::sqrt(scale));
    SetShift(mean - gamma * (1.0 / kappa - kappa));
}

void LaplaceRand::FitLocationAndAsymmetryMM(const std::vector<double> &sample)
{
    double mean = sampleMean(sample);
    double skewness = sampleSkewness(sample, mean);
    /// Getting asymmetry from skewness (because from variance result is ambiguous)
    double root = GetAsymmetryFromSkewness(skewness);
    if (!std::isfinite(root))
        throw std::runtime_error(fitError(UNDEFINED_ERROR, "Error in root-finding procedure"));
    SetAsymmetry(root);
    SetShift(mean - gamma * (1.0 / kappa - kappa));
}

void LaplaceRand::FitScaleAndAsymmetryMM(const std::vector<double> &sample)
{
    double mean = sampleMean(sample);
    double var = sampleVariance(sample, mean);
    double z = mean - m;
    double a = var / (z * z);
    double am1 = a - 1;
    double asymmetrySq = (a - std::sqrt(a + am1)) / am1;
    SetAsymmetry(std::sqrt(asymmetrySq));
    SetScale(z * kappa / (1.0 - kappaSq));
}

void LaplaceRand::FitLocationScaleAndAsymmetryMM(const std::vector<double> &sample)
{
    double mean = sampleMean(sample);
    double var = sampleVariance(sample, mean);
    double skewness = sampleSkewness(sample, mean, std::sqrt(var));
    double root = GetAsymmetryFromSkewness(skewness);
    if (!std::isfinite(root))
        throw std::runtime_error(fitError(UNDEFINED_ERROR, "Error in root-finding procedure"));
    SetAsymmetry(root);
    double scale = (var * kappaSq) / (1 + kappaSq * kappaSq);
    SetScale(std::sqrt(scale));
    SetShift(mean - gamma * (1.0 / kappa - kappa));
}
