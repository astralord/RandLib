#include "DiscreteDistribution.h"

void DiscreteDistribution::ProbabilityMassFunction(const std::vector<int> &x, std::vector<double> &y) const
{
    size_t size = x.size();
    if (size > y.size())
        return;
    for (size_t i = 0; i != size; ++i)
        y[i] = P(x[i]);
}

int DiscreteDistribution::Mode() const
{
    /// Works only for unimodal and monotone from start point to mode distributions
    int x = Median();
    double prob = P(x), newProb = P(x + 1);
    if (prob < newProb) {
        do {
            ++x;
            prob = newProb;
            newProb = P(x + 1);
        } while (prob < newProb);
    }
    else {
        newProb = P(x - 1);
        while (prob < newProb) {
            --x;
            prob = newProb;
            newProb = P(x - 1);
        }
    }
    return x;
}

double DiscreteDistribution::QuantileImpl(double probability) const
{
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

double DiscreteDistribution::ExpectedValue(const std::function<double (double)> &funPtr, int minPoint, int maxPoint) const
{
    SUPPORT_TYPE suppType = supportType();
    int k = minPoint, upperBoundary = maxPoint;
    if (suppType == FINITE_T || suppType == RIGHTSEMIFINITE_T) {
        k = std::max(k, static_cast<int>(std::floor(MinValue())));
    }
    if (suppType == FINITE_T || suppType == LEFTSEMIFINITE_T) {
        upperBoundary = std::min(upperBoundary, static_cast<int>(std::ceil(MaxValue())));
    }

    double sum = 0;
    do {
        double addon = funPtr(k);
        if (addon != 0.0) {
            addon *= P(k);
            sum += addon;
        }
        ++k;
    } while (k <= upperBoundary);
    return sum;
}

double DiscreteDistribution::ExpectedValue(const std::function<double (double)> &funPtr, double startPoint) const
{
    SUPPORT_TYPE suppType = supportType();
    double sum = 0, addon = 0, prob = 0;
    if (suppType == FINITE_T) {
        int k = std::floor(MinValue());
        int upperBoundary = std::ceil(MaxValue());
        do {
            addon = funPtr(k);
            if (addon != 0.0) {
                prob = P(k);
                if (prob > 0.0) {
                    addon *= prob;
                    sum += addon;
                }
            }
            ++k;
        } while (k <= upperBoundary);
        return sum;
    }

    static constexpr double epsilon = 1e-12;
    static constexpr int maxIter = 1e4;
    int iter = 0;

    if (suppType == RIGHTSEMIFINITE_T) {
        int k = std::floor(MinValue());
        do {
            addon = funPtr(k);
            if (addon != 0.0) {
                prob = P(k);
                if (prob > 0.0) {
                    addon *= prob;
                    sum += addon;
                }
            }
            ++k;
            if (++iter > maxIter) /// can't take sum, addon decreases too slow
                return INFINITY;
        } while (prob > epsilon || std::fabs(addon) > epsilon || F(k) < 0.999);
        return sum;
    }

    if (suppType == LEFTSEMIFINITE_T) {
        int k = std::floor(MaxValue());
        do {
            addon = funPtr(k);
            if (addon != 0.0) {
                prob = P(k);
                if (prob > 0.0) {
                    addon *= prob;
                    sum += addon;
                }
            }
            --k;
            if (++iter > maxIter) /// can't take sum, addon decreases too slow
                return INFINITY;
        } while (prob > epsilon || std::fabs(addon) > epsilon || F(k) > 0.001);
        return sum;
    }

    // TODO: transform all this spagetti code into good one

    double x = std::floor(startPoint);
    if (RandMath::areClose(x, startPoint, epsilon))
        --x;
    int k = x;

    do {
        addon = funPtr(k);
        if (addon != 0.0) {
            prob = P(k);
            if (prob > 0.0) {
                addon *= prob;
                sum += addon;
            }
        }
        --k;
        if (++iter > maxIter) /// can't take sum, addon decreases too slow
            return INFINITY;
    } while (prob > epsilon || std::fabs(addon) > epsilon || F(k) < 0.999);

    iter = 0;
    prob = 0;
    x = std::ceil(startPoint);
    if (RandMath::areClose(x, startPoint, epsilon))
        ++x;
    k = x;
    do {
        addon = funPtr(k);
        if (addon != 0.0) {
            prob = P(k);
            if (prob > 0.0) {
                addon *= prob;
                sum += addon;
            }
        }
        ++k;
        if (++iter > maxIter) /// can't take sum, addon decreases too slow
            return INFINITY;
    } while (prob > epsilon || std::fabs(addon) > epsilon || F(k) > 0.001);

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
