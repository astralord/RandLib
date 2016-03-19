#include "ProbabilityDistribution.h"
#include <sstream>      // std::ostringstream
#include <iomanip>      // std::setprecision

ProbabilityDistribution::ProbabilityDistribution()
{
}

std::string ProbabilityDistribution::toStringWithPrecision(const double a_value, const int n)
{
    std::ostringstream out;
    out << std::setprecision(n) << a_value;
    return out.str();
}

void ProbabilityDistribution::sample(QVector<double> &outputData) const
{
    for (double &var : outputData)
        var = variate();
}

void ProbabilityDistribution::cdf(const QVector<double> &x, QVector<double> &y) const
{
    int size = std::min(x.size(), y.size());
    for (int i = 0; i != size; ++i)
        y[i] = F(x[i]);
}

std::complex<double> ProbabilityDistribution::CF(double t) const
{
    // TODO:
    return std::complex<double>(t);
}

void ProbabilityDistribution::cf(const QVector<double> &t, QVector<std::complex<double> > &y) const
{
    int size = std::min(t.size(), y.size());
    for (int i = 0; i != size; ++i)
        y[i] = CF(t[i]);
}

double ProbabilityDistribution::Median() const
{
    return Quantile(0.5);
}

double ProbabilityDistribution::Skewness() const
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

double ProbabilityDistribution::ExcessKurtosis() const
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

double ProbabilityDistribution::Kurtosis() const
{
    return ExcessKurtosis() + 3.0;
}
