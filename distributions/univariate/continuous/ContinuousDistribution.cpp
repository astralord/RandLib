#include "ContinuousDistribution.h"

void ContinuousDistribution::ProbabilityDensityFunction(const std::vector<double> &x, std::vector<double> &y) const
{
    size_t size = x.size();
    if (size > y.size())
        return;
    for (size_t i = 0; i != size; ++i)
        y[i] = f(x[i]);
}

double ContinuousDistribution::quantileImpl(double p) const
{
    double guess = 0.0;
    SUPPORT_TYPE supp = SupportType();
    if (supp == FINITE_T) {
        if (RandMath::findRoot([this, p] (double x)
        {
            return F(x) - p;
        }, MinValue(), MaxValue(), guess))
            return guess;
        return NAN;
    }

    /// We use quantile from sample as an initial guess
    static constexpr int SAMPLE_SIZE = 100;
    std::vector<double> sample(SAMPLE_SIZE);
    this->Sample(sample);
    std::sort(sample.begin(), sample.end());
    guess = sample[p * SAMPLE_SIZE];

    if (RandMath::findRoot([this, p] (double x)
    {
        double first = F(x) - p;
        double second = f(x);
        return DoublePair(first, second);
    }, guess))
        return guess;
    return NAN;
}

double ContinuousDistribution::quantileImpl1m(double p) const
{
    double guess = 0.0;
    SUPPORT_TYPE supp = SupportType();
    if (supp == FINITE_T) {
        if (RandMath::findRoot([this, p] (double x)
        {
            return S(x) - p;
        }, MinValue(), MaxValue(), guess))
            return guess;
        return NAN;
    }

    /// We use quantile from sample as an initial guess
    static constexpr int SAMPLE_SIZE = 100;
    std::vector<double> sample(SAMPLE_SIZE);
    this->Sample(sample);
    /// Sort in desceding order
    std::sort(sample.begin(), sample.end(), std::greater<>());
    guess = sample[p * SAMPLE_SIZE];

    if (RandMath::findRoot([this, p] (double x)
    {
        double first = p - S(x);
        double second = f(x);
        return DoublePair(first, second);
    }, guess))
        return guess;
    return NAN;
}

double ContinuousDistribution::Hazard(double x) const
{
    if (x < MinValue())
        return 0.0; /// 0/1
    if (x > MaxValue())
        return NAN; /// 0/0
    return f(x) / S(x);
}

double ContinuousDistribution::Mode() const
{
    double guess = Mean(); /// good starting point
    if (!std::isfinite(guess))
        guess = Median(); /// this shouldn't be nan or inf
    double root = 0;
    RandMath::findMin([this] (double x)
    {
        return -f(x);
    }, guess, root);
    return root;
}

double ContinuousDistribution::ExpectedValue(const std::function<double (double)> &funPtr, double minPoint, double maxPoint) const
{
    /// attempt to calculate expected value by numerical method
    /// use for distributions w/o explicit formula
    double lowerBoundary = minPoint, upperBoundary = maxPoint;
    if (isRightBounded())
        lowerBoundary = std::max(minPoint, lowerBoundary);
    if (isLeftBounded())
        upperBoundary = std::min(maxPoint, upperBoundary);

    if (lowerBoundary >= upperBoundary)
        return 0.0;

    double mode = this->Mode();
    bool isLeftBoundFinite = std::isfinite(lowerBoundary), isRightBoundFinite = std::isfinite(upperBoundary);

    double res = 0.0;
    /// Integrate from left boundary to mode
    if (isLeftBoundFinite) {
        res += RandMath::integral([this, funPtr] (double x)
        {
            return funPtr(x) * f(x);
        },
        lowerBoundary, mode);
    }
    else {
        res += RandMath::integral([this, funPtr, mode] (double x)
        {
            if (x <= 0.0)
                return 0.0;
            double t = mode - (1.0 - x) / x;
            double y = funPtr(t) * f(t);
            return y / (x * x);
        },
        0.0, 1.0);
    }

    /// Integrate from mode to right boundary
    if (isRightBoundFinite) {
        res += RandMath::integral([this, funPtr] (double x)
        {
            return funPtr(x) * f(x);
        },
        mode, upperBoundary);
    }
    else {
        res += RandMath::integral([this, funPtr, mode] (double x)
        {
            if (x >= 1.0)
                return 0.0;
            double denom = 1.0 - x;
            double t = mode + x / denom;
            double y = funPtr(t) * f(t);
            denom *= denom;
            return y / denom;
        },
        0.0, 1.0);
    }
    return res;
}

double ContinuousDistribution::Likelihood(const std::vector<double> &sample) const
{
    double res = 1.0;
    for (const double & var : sample )
        res *= f(var);
    return res;
}

double ContinuousDistribution::LogLikelihood(const std::vector<double> &sample) const
{
    double res = 0.0;
    for (const double & var : sample )
        res += std::log(f(var));
    return res;
}
