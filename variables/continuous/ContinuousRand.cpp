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
    /// use for distributions w/o explicit formula
    /// works good for unimodal and wide distributions

    double mu = E();
    if (std::isnan(mu) || std::isinf(mu))
        return NAN;

    double var = Var();
    if (std::isnan(var) || std::isinf(var))
        return NAN;

    double sigma = std::sqrt(var);

    /// get such low boundary 'x' that |integrand(x)| < epsilon
    static constexpr double epsilon = 1e-10;
    static constexpr int maxIter = 1000;
    int iter = 0;
    double lowBoundary = mu;
    double integrand = 0;
    do {
       lowBoundary -= var;
       double aux = (lowBoundary - mu) / sigma;
       integrand = aux * aux * aux;
       integrand *= f(lowBoundary);
    } while (-integrand > epsilon && ++iter < maxIter);

    if (iter == maxIter) /// can't take integral, integrand decreases too slow
        return NAN;

    /// get such upper boundary 'x' that |integrand(x)| < epsilon
    double upperBoundary = mu;
    iter = 0;
    do {
       upperBoundary += var;
       double aux = (upperBoundary - mu) / sigma;
       integrand = aux * aux * aux;
       integrand *= f(upperBoundary);
    } while (integrand > epsilon && ++iter < maxIter);

    if (iter == maxIter) /// can't take integral, integrand decreases too slow
        return NAN;

    double integral = RandMath::integral([this, mu] (double x)
    {
        double aux = x - mu;
        return aux * aux * aux * f(x);
    },
    lowBoundary, upperBoundary, epsilon);

    return integral / (sigma * var);
}

double ContinuousRand::ExcessKurtosis() const
{
    /// WARNING: attempt to calculate kurtosis by numerical method
    /// use for distributions w/o explicit formula
    /// works good for unimodal and wide distributions
    double mu = E();
    if (std::isnan(mu) || std::isinf(mu))
        return NAN;

    double var = Var();
    if (std::isnan(var) || std::isinf(var))
        return NAN;

    /// get such low boundary 'x' that |integrand(x)| < epsilon
    static constexpr double epsilon = 1e-10;
    double lowBoundary = mu;
    static constexpr int maxIter = 1000;
    int iter = 0;
    double integrand = 0;
    do {
       lowBoundary -= var;
       integrand = lowBoundary - mu;
       integrand *= integrand / var;
       integrand *= integrand;
       integrand *= f(lowBoundary);
    } while (integrand > epsilon && ++iter < maxIter);

    if (iter == maxIter) /// can't take integral, integrand decreases too slow
        return NAN;

    /// get such upper boundary 'x' that |integrand(x)| < epsilon
    double upperBoundary = mu;
    iter = 0;
    do {
       upperBoundary += var;
       integrand = upperBoundary - mu;
       integrand *= integrand / var;
       integrand *= integrand;
       integrand *= f(upperBoundary);
    } while (integrand > epsilon && ++iter < maxIter);

    if (iter == maxIter) /// can't take integral, integrand decreases too slow
        return NAN;

    double integral = RandMath::integral([this, mu] (double x)
    {
        double aux = x - mu;
        aux *= aux;
        return aux * aux * f(x);
    },
    lowBoundary, upperBoundary, epsilon);

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
