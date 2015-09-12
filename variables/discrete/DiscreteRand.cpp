#include "DiscreteRand.h"

void DiscreteRand::pmf(const QVector<int> &x, QVector<double> &y) const
{
    if (x.size() != y.size())
        return;
    for (int i = 0; i != x.size(); ++i)
        y[i] = P(x[i]);
}

double DiscreteRand::Skewness() const
{
    /// Calculate skewness using Monte-Carlo method
    /// use only for distributions w/o explicit formula
    static constexpr int size = 1e7;
    long double skewness = 0.0;
    double mu = E();
    double var = Var();
    double sigma = std::sqrt(var);
    for (int i = 0; i != size; ++i)
    {
        double sample = variate() - mu;
        skewness += sample * sample * sample;
    }
    skewness /= (sigma * var);
    return skewness / size;
}

double DiscreteRand::ExcessKurtosis() const
{
    /// Calculate kurtosis using Monte-Carlo method
    /// use only for distributions w/o explicit formula
    static constexpr int size = 1e7;
    long double skewness = 0.0;
    double mu = E();
    double var = Var();
    for (int i = 0; i != size; ++i)
    {
        double sample = variate() - mu;
        sample *= sample;
        skewness += sample * sample;
    }
    skewness /= (var * var);
    return skewness / size - 3;
}

double DiscreteRand::likelihood(const QVector<int> &sample) const
{
    double res = 1.0;
    for (int i = 0; i != sample.size(); ++i)
        res *= std::log(P(sample[i]));
    return res;
}

double DiscreteRand::loglikelihood(const QVector<int> &sample) const
{
    double res = 0.0;
    for (int i = 0; i != sample.size(); ++i)
        res += std::log(P(sample[i]));
    return res;
}
