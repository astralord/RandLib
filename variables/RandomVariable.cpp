#include "RandomVariable.h"
#include <sstream>      // std::ostringstream
#include <iomanip>      // std::setprecision

RandomVariable::RandomVariable()
{
}

std::string RandomVariable::toStringWithPrecision(const double a_value, const int n)
{
    std::ostringstream out;
    out << std::setprecision(n) << a_value;
    return out.str();
}

void RandomVariable::sample(QVector<double> &outputData)
{
    for (double &var : outputData)
        var = variate();
}

void RandomVariable::cdf(const QVector<double> &x, QVector<double> &y)
{
    int size = std::min(x.size(), y.size());
    for (int i = 0; i != size; ++i)
        y[i] = F(x[i]);
}

void RandomVariable::cf(const QVector<double> &t, QVector<std::complex<double> > &y)
{
    int size = std::min(t.size(), y.size());
    for (int i = 0; i != size; ++i)
        y[i] = CF(t[i]);
}

double RandomVariable::Skewness() const
{
    double mu = E();
    if (std::isnan(mu) || std::isinf(mu))
        return NAN;

    double var = Var();
    if (std::isnan(var) || std::isinf(var))
        return NAN;

    double sum = ExpectedValue([this, mu] (double x)
    {
        double skew = x - mu;
        return skew * skew * skew;
    }, mu);

    double sigma = std::sqrt(var);

    return sum / (var * sigma);
}

double RandomVariable::ExcessKurtosis() const
{
    double mu = E();
    if (std::isnan(mu) || std::isinf(mu))
        return NAN;

    double var = Var();
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
