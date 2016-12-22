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
    double root = 0.0;
    SUPPORT_TYPE supp = SupportType();
    if (supp == FINITE_T) {
        if (RandMath::findRoot([this, p] (double x)
        {
            return F(x) - p;
        }, MinValue(), MaxValue(), root))
            return root;
        return NAN;
    }
    double mean = Mean();
    if (!std::isfinite(mean))
    {
        if (isLeftBounded()) {
            root = MinValue() + 1;
            while(f(root) <= MIN_POSITIVE)
                ++root;
        }
        else if (isRightBounded()) {
            root = MaxValue() - 1;
            while(f(root) <= MIN_POSITIVE)
                --root;
        }
        else {
            root = 0.0;
        }
    }
    else {
        root = mean;
    }

    if (RandMath::findRoot([this, p] (double x)
    {
        double first = F(x) - p;
        double second = f(x);
        return DoublePair(first, second);
    }, root))
        return root;
    return NAN;
}

double ContinuousDistribution::quantileImpl1m(double p) const
{
    double root = 0.0;
    SUPPORT_TYPE supp = SupportType();
    if (supp == FINITE_T) {
        if (RandMath::findRoot([this, p] (double x)
        {
            return F(x) - p;
        }, MinValue(), MaxValue(), root))
            return root;
        return NAN;
    }
    double mean = Mean();
    if (!std::isfinite(mean))
    {
        if (isLeftBounded()) {
            root = MinValue() + 1;
            while(f(root) <= MIN_POSITIVE)
                ++root;
        }
        else if (isRightBounded()) {
            root = MaxValue() - 1;
            while(f(root) <= MIN_POSITIVE)
                --root;
        }
        else {
            root = 0.0;
        }
    }
    else {
        root = mean;
    }

    if (RandMath::findRoot([this, p] (double x)
    {
        double first = F(x) - 1;
        first += p;
        double second = f(x);
        return DoublePair(first, second);
    }, root))
        return root;
    return NAN;
}

double ContinuousDistribution::Hazard(double x) const
{
    return f(x) / (1.0 - F(x));
}

double ContinuousDistribution::Mode() const
{
    double mu = Mean(); /// good starting point
    if (!std::isfinite(mu))
        mu = Median(); /// this shouldn't be nan or inf
    double root = 0;
    RandMath::findMin([this] (double x)
    {
        return -f(x);
    }, mu, root);
    return root;
}

double ContinuousDistribution::GetLeftLimit(double value, double epsilon) const
{
    if (!std::isfinite(f(value))) {
        if (std::fabs(value) < 1)
            value -= epsilon;
        else
            value *= 0.9999;
    }
    return value;
}

double ContinuousDistribution::GetRightLimit(double value, double epsilon) const
{
    if (!std::isfinite(f(value))) {
        if (std::fabs(value) < 1)
            value += epsilon;
        else
            value *= 1.0001;
    }
    return value;
}

double ContinuousDistribution::ExpectedValue(const std::function<double (double)> &funPtr, double minPoint, double maxPoint) const
{
    static constexpr double epsilon = 1e-10;
    double lowerBoundary = minPoint, upperBoundary = maxPoint;
    if (isRightBounded()) {
        lowerBoundary = std::max(minPoint, lowerBoundary);
        lowerBoundary = GetLeftLimit(lowerBoundary, epsilon);
    }
    if (isLeftBounded()) {
        upperBoundary = std::min(maxPoint, upperBoundary);
        upperBoundary = GetRightLimit(upperBoundary, epsilon);
    }

    if (lowerBoundary >= upperBoundary)
        return 0.0;

    return RandMath::integral([this, funPtr] (double x)
    {
        return funPtr(x) * f(x);
    },
    lowerBoundary, upperBoundary);
}


double ContinuousDistribution::ExpectedValue(const std::function<double (double)> &funPtr, double startPoint) const
{
    /// attempt to calculate expected value by numerical method
    /// use for distributions w/o explicit formula
    /// works good for unimodal distributions
    static constexpr double epsilon = 1e-10;

    // TODO: Get rid of startpoint
    double lowerBoundary = startPoint, upperBoundary = startPoint;
    SUPPORT_TYPE suppType = SupportType();
    if (suppType == FINITE_T)
    {
        lowerBoundary = GetLeftLimit(MinValue(), epsilon);
        upperBoundary = GetRightLimit(MaxValue(), epsilon);

        return RandMath::integral([this, funPtr] (double x)
        {
            return funPtr(x) * f(x);
        },
        lowerBoundary, upperBoundary);
    }
    else if (suppType == RIGHTSEMIFINITE_T) {
        lowerBoundary = GetLeftLimit(MinValue(), epsilon);
        return RandMath::integral([this, funPtr, lowerBoundary] (double x)
        {
            if (x >= 1.0)
                return 0.0;
            double denom = 1.0 - x;
            double t = lowerBoundary + x / denom;
            double y = funPtr(t) * f(t);
            denom *= denom;
            return y / denom;
        },
        0.0, 1.0);
    } else if (suppType == LEFTSEMIFINITE_T) {
        upperBoundary = GetRightLimit(MaxValue(), epsilon);
        return RandMath::integral([this, funPtr, upperBoundary] (double x)
        {
            if (x <= 0.0)
                return 0.0;
            double t = upperBoundary - (1.0 - x) / x;
            double y = funPtr(t) * f(t);
            return y / (x * x);
        },
        0.0, 1.0);
    }
    /// Infinite case
    return RandMath::integral([this, funPtr] (double x)
    {
        if (std::fabs(x) >= 1.0)
            return 0.0;
        double x2 = x * x;
        double denom = 1.0 - x2;
        double t = x / denom;
        double y = funPtr(t) * f(t);
        denom *= denom;
        return y * (1.0 + x2) / denom;
    },
    -1.0, 1.0);

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
