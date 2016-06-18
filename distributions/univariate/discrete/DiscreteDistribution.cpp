#include "DiscreteDistribution.h"

void DiscreteDistribution::ProbabilityMassFunction(const std::vector<int> &x, std::vector<double> &y) const
{
    size_t size = x.size();
    if (size > y.size())
        return;
    for (size_t i = 0; i != size; ++i)
        y[i] = P(x[i]);
}

double DiscreteDistribution::Quantile(double probability) const
{
    if (probability < 0 || probability > 1)
        return NAN;
    if (probability == 0.0)
        return -INFINITY;
    double mean = Mean();
    int down = static_cast<int>(std::floor(mean)), up = down + 1;
    double fu = F(up), fd = F(down);
    int iter = 0;
    static constexpr int maxIter = 1e4;
    /// go up
    while (fu < probability)
    {
        fd = fu;
        fu = F(++up);
        if (++iter > maxIter)
            return NAN;
    }
    down = up - 1;

    iter = 0;
    /// go down
    while (fd > probability)
    {
        fd = F(--down);
        if (++iter > maxIter)
            return NAN;
    }
    up = down + 1;

    /// if lower quantile is not equal probability, we return upper quantile
    return (fd < probability) ? up : down;
}

double DiscreteDistribution::Hazard(double x) const
{
    return P(x) / (1.0 - F(x));
}

double DiscreteDistribution::ExpectedValue(const std::function<double (double)> &funPtr, double startPoint) const
{
    static constexpr double epsilon = 1e-12;
    static constexpr int maxIter = 1e4;
    int iter = 0;
    long double sum = 0.0L;
    double addon = 0;
    double x = std::floor(startPoint);
    if (RandMath::areClose(x, startPoint, epsilon))
        --x;

    do {
        addon = funPtr(x);
        if (addon != 0.0)
            addon *= P(x);
        sum += addon;
        --x;
        if (++iter > maxIter) /// can't take sum, addon decreases too slow
            return INFINITY;
    } while (std::fabs(addon) > epsilon);

    iter = 0;
    x = std::ceil(startPoint);
    if (RandMath::areClose(x, startPoint, epsilon))
        ++x;
    do {
        addon = funPtr(x);
        if (addon != 0.0)
            addon *= P(x);
        addon *= P(x);
        sum += addon;
        ++x;
        if (++iter > maxIter) /// can't take sum, addon decreases too slow
            return INFINITY;
    } while (std::fabs(addon) > epsilon);

    return sum;
}

double DiscreteDistribution::Likelihood(const std::vector<int> &sample) const
{
    double res = 1.0;
    for (const int & var : sample )
        res *= P(var);
    return res;
}

double DiscreteDistribution::LogLikelihood(const std::vector<int> &sample) const
{
    double res = 0.0;
    for (const int & var : sample )
        res += std::log(P(var));
    return res;
}
