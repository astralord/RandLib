#include "ContinuousRand.h"

void ContinuousRand::pdf(const QVector<double> &x, QVector<double> &y) const
{
    size_t size = std::min(x.size(), y.size());
    for (size_t i = 0; i != size; ++i)
        y[i] = f(x[i]);
}

double ContinuousRand::Skewness() const
{
    /// WARNING: attempt to calculate skewness by numerical method
    /// (use for those distributions which have no explicit formula)
    double mu = E();
    double var = Var();
    double sigma = std::sqrt(var);
    double sigma3 = 3.0 * sigma;
    static constexpr double epsilon = 1e-10;
    /// get such low boundary 'x' that |integrand(x)| < epsilon
    double lowBoundary = mu;
    double integrand = 0;
    do {
       lowBoundary -= sigma3;
       double aux = (lowBoundary - mu) / sigma;
       integrand = aux * aux * aux;
       integrand *= f(lowBoundary);
    } while (-integrand > epsilon);


    /// get such upper boundary 'x' that |integrand(x)| < epsilon
    double upperBoundary = mu;
    do {
       upperBoundary += sigma3;
       double aux = (upperBoundary - mu) / sigma;
       integrand = aux * aux * aux;
       integrand *= f(upperBoundary);
    } while (integrand > epsilon);

    double integral = RandMath::integral([this, mu] (double x)
    {
        double aux = x - mu;
        return aux * aux * aux * f(x);
    },
    lowBoundary, upperBoundary);

    return integral / (sigma * var);
}

double ContinuousRand::ExcessKurtosis() const
{
    /// WARNING: attempt to calculate kurtosis by numerical method
    /// (use for those distributions which have no explicit formula)
    double mu = E();
    double var = Var();
    double var3 = 3.0 * var;
    static constexpr double epsilon = 1e-10;
    /// get such low boundary 'x' that |integrand(x)| < epsilon
    double lowBoundary = mu;
    double integrand = 0;
    do {
       lowBoundary -= var3;
       integrand = lowBoundary - mu;
       integrand *= integrand / var;
       integrand *= integrand;
       integrand *= f(lowBoundary);
    } while (integrand > epsilon);


    /// get such upper boundary 'x' that |integrand(x)| < epsilon
    double upperBoundary = mu;
    do {
       upperBoundary += var3;
       integrand = upperBoundary - mu;
       integrand *= integrand / var;
       integrand *= integrand;
       integrand *= f(upperBoundary);
    } while (integrand > epsilon);

    double integral = RandMath::integral([this, mu] (double x)
    {
        double aux = x - mu;
        aux *= aux;
        return aux * aux * f(x);
    },
    lowBoundary, upperBoundary);

    return integral / (var * var) - 3;
}

double ContinuousRand::likelihood(const QVector<double> &sample) const
{
    double res = 1.0;
    for (int i = 0; i != sample.size(); ++i)
        res *= f(sample[i]);
    return res;
}

double ContinuousRand::loglikelihood(const QVector<double> &sample) const
{
    double res = 0.0;
    for (int i = 0; i != sample.size(); ++i)
        res += std::log(f(sample[i]));
    return res;
}
