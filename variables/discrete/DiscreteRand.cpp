#include "DiscreteRand.h"

void DiscreteRand::pmf(const QVector<int> &x, QVector<double> &y) const
{
    if (x.size() != y.size())
        return;
    for (int i = 0; i != x.size(); ++i)
        y[i] = P(x[i]);
}

double DiscreteRand::Quantile(double probability) const
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

double DiscreteRand::Hazard(double x) const
{
    return P(x) / (1.0 - F(x));
}

double DiscreteRand::ExpectedValue(const std::function<double (double)> &funPtr, double startPoint) const
{
    static constexpr double epsilon = 1e-12;
    static constexpr int maxIter = 1e4; // why so small?
    int iter = 0;
    long double sum = 0.0L;
    double addon = 0;
    double x = std::floor(startPoint);
    if (RandMath::areEqual(x, startPoint, epsilon))
        --x;

    do {
        addon = funPtr(x);
        addon *= P(x);
        sum += addon;
        --x;
        if (++iter > maxIter) /// can't take sum, addon decreases too slow
            return INFINITY;
    } while (std::fabs(addon) > epsilon);

    //TODO: why we need next line?
    if (iter == maxIter) /// can't take sum, addon decreases too slow
        return INFINITY;

    iter = 0;
    x = std::ceil(startPoint);
    if (RandMath::areEqual(x, startPoint, epsilon))
        ++x;
    do {
        addon = funPtr(x);
        addon *= P(x);
        sum += addon;
        ++x;
        if (++iter > maxIter) /// can't take sum, addon decreases too slow
            return INFINITY;
    } while (std::fabs(addon) > epsilon);

    return sum;
}

double DiscreteRand::likelihood(const QVector<int> &sample) const
{
    double res = 1.0;
    for (const int & var : sample )
        res *= P(var);
    return res;
}

double DiscreteRand::loglikelihood(const QVector<int> &sample) const
{
    double res = 0.0;
    for (const int & var : sample )
        res += std::log(P(var));
    return res;
}
