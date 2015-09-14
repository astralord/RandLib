#include "DiscreteRand.h"

void DiscreteRand::pmf(const QVector<int> &x, QVector<double> &y) const
{
    if (x.size() != y.size())
        return;
    for (int i = 0; i != x.size(); ++i)
        y[i] = P(x[i]);
}

double DiscreteRand::ExpectedValue(const std::function<double (double)> &funPtr, double startPoint) const
{
    static constexpr double epsilon = 1e-12;
    static constexpr int maxIter = 1e4;
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
    } while (std::fabs(addon) > epsilon && ++iter < maxIter);

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
    } while (std::fabs(addon) > epsilon && ++iter < maxIter);

    if (iter == maxIter) /// can't take sum, addon decreases too slow
        return INFINITY;

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
