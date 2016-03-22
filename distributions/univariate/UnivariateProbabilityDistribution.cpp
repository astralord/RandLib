#include "UnivariateProbabilityDistribution.h"

UnivariateProbabilityDistribution::UnivariateProbabilityDistribution()
{
}

std::complex<double> UnivariateProbabilityDistribution::CF(double t) const
{
    // TODO:
    return std::complex<double>(t);
}

void UnivariateProbabilityDistribution::cf(const QVector<double> &t, QVector<std::complex<double> > &y) const
{
    int size = std::min(t.size(), y.size());
    for (int i = 0; i != size; ++i)
        y[i] = CF(t[i]);
}

double UnivariateProbabilityDistribution::Median() const
{
    return Quantile(0.5);
}

double UnivariateProbabilityDistribution::Skewness() const
{
    double mu = Mean();
    if (std::isnan(mu) || std::isinf(mu))
        return NAN;

    double var = Variance();
    if (std::isnan(var) || std::isinf(var))
        return NAN;

    double sum = ExpectedValue([this, mu] (double x)
    {
        double skew = x - mu;
        return skew * skew * skew;
    }, mu);

    return sum / (var * std::sqrt(var));
}

double UnivariateProbabilityDistribution::ExcessKurtosis() const
{
    double mu = Mean();
    if (std::isnan(mu) || std::isinf(mu))
        return NAN;

    double var = Variance();
    if (std::isnan(var) || std::isinf(var))
        return NAN;

    double sum = ExpectedValue([this, mu] (double x)
    {
        double kurtosis = x - mu;
        kurtosis *= kurtosis;
        return kurtosis * kurtosis;
    }, mu);

    return sum / (var * var) - 3;
}

double UnivariateProbabilityDistribution::Kurtosis() const
{
    return ExcessKurtosis() + 3.0;
}
