#include "ContinuousRand.h"

void ContinuousRand::pdf(const QVector<double> &x, QVector<double> &y) const
{
    size_t size = std::min(x.size(), y.size());
    for (size_t i = 0; i != size; ++i)
        y[i] = f(x[i]);
}

double ContinuousRand::ExpectedValue(const std::function<double (double)> &funPtr, double startPoint) const
{
    /// attempt to calculate expected value by numerical method
    /// use for distributions w/o explicit formula
    /// works good for unimodal and wide distributions

    /// get such low boundary 'x' that |integrand(x)| < epsilon
    static constexpr double epsilon = 1e-10;
    static constexpr int maxIter = 1000;
    int iter = 0;
    double lowBoundary = startPoint;
    double integrand = 0;
    /// WARNING: we use variance - so there can be deadlock if we don't define this function explicitly
    /// therefore function Variance() should stay pure and noone should calculate it by this function
    double var = Var();
    do {
       lowBoundary -= var;
       integrand = funPtr(lowBoundary) * f(lowBoundary);
    } while (std::fabs(integrand) > epsilon && ++iter < maxIter);

    if (iter == maxIter) /// can't take integral, integrand decreases too slow
        return NAN;

    /// get such upper boundary 'x' that |integrand(x)| < epsilon
    double upperBoundary = startPoint;
    iter = 0;
    do {
       upperBoundary += var;
       integrand = funPtr(upperBoundary) * f(upperBoundary);
    } while (std::fabs(integrand) > epsilon && ++iter < maxIter);

    if (iter == maxIter) /// can't take integral, integrand decreases too slow
        return NAN;

    return RandMath::integral([this, funPtr] (double x)
    {
        return funPtr(x) * f(x);
    },
    lowBoundary, upperBoundary, epsilon);
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
