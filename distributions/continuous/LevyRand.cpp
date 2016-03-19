#include "LevyRand.h"
#include "NormalRand.h"

LevyRand::LevyRand(double location, double scale) : StableRand(0.5, 1)
{
    setParameters(location, scale);
}

std::string LevyRand::name()
{
    return "Levy(" + toStringWithPrecision(getLocation()) + ", " + toStringWithPrecision(getScale()) + ")";
}

void LevyRand::setParameters(double location, double scale)
{
    StableRand::setParameters(0.5, 1.0, scale, location);
    pdfCoef = M_1_SQRT2PI * std::sqrt(sigma);
}

double LevyRand::f(double x) const
{
    if (x <= mu)
        return 0;
    double xInv = 1.0 / (x - mu);
    double y = -0.5 * sigma * xInv;
    y = std::exp(y);
    y *= xInv;
    y *= std::sqrt(xInv);
    return pdfCoef * y;
}

double LevyRand::F(double x) const
{
    if (x <= mu)
        return 0;
    double y = x - mu;
    y += y;
    y = sigma / y;
    y = std::sqrt(y);
    return std::erfc(y);
}

double LevyRand::variate() const
{
    double rv = NormalRand::standardVariate();
    rv *= rv;
    rv = sigma / rv;
    return mu + rv;
}

double LevyRand::variate(double location, double scale)
{
    return location + scale * standardVariate();
}

double LevyRand::standardVariate()
{
    double rv = NormalRand::standardVariate();
    return 1.0 / (rv * rv);
}

std::complex<double> LevyRand::CF(double t) const
{
    std::complex<double> y(0.0, -2 * sigma * t);
    y = -std::sqrt(y);
    y += std::complex<double>(0.0, mu * t);
    return std::exp(y);
}

double LevyRand::Quantile(double p) const
{
    if (p == 0.0)
        return mu;
    double x = pdfCoef * M_SQRT2PI / NormalRand::standardQuantile(0.5 * p);
    return mu + x * x;
}

double LevyRand::Mode() const
{
    return sigma / 3.0 + mu;
}

bool LevyRand::checkValidity(const QVector<double> &sample)
{
    for (double var : sample) {
        if (var <= mu)
            return false;
    }
    return true;
}

bool LevyRand::fitScale_MLE(const QVector<double> &sample)
{
    double n = sample.size();
    if (n <= 0 || !checkValidity(sample))
        return false;
    long double invSum = 0.0;
    for (double var : sample)
        invSum += 1.0 / (var - mu);
    setParameters(mu, n / invSum);
    return true;
}
