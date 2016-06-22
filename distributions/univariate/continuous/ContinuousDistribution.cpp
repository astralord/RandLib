#include "ContinuousDistribution.h"

void ContinuousDistribution::ProbabilityDensityFunction(const std::vector<double> &x, std::vector<double> &y) const
{
    size_t size = x.size();
    if (size > y.size())
        return;
    for (size_t i = 0; i != size; ++i)
        y[i] = f(x[i]);
}

double ContinuousDistribution::Quantile(double p) const
{
    if (p < 0 || p > 1)
        return NAN;

    double root = Mean(); /// good starting point
    if (!std::isfinite(root))
        root = 0.0;
    if (RandMath::findRoot([this, p] (double x)
    {
        return F(x) - p;
    },
    [this] (double x)
    {
        return f(x);
    }, root))
        return root;
    /// if we can't find quantile, then probably p == 1
    return INFINITY;
}

double ContinuousDistribution::Hazard(double x) const
{
    return f(x) / (1.0 - F(x));
}

double ContinuousDistribution::Mode() const
{
    // use only for unimodal distributions!

    double mu = Mean(); /// good starting point
    if (!std::isfinite(mu))
        mu = Median(); /// this shouldn't be nan or inf
    double step = 10 * Variance();
    if (!std::isfinite(step))
        step = 100; // dirty hack

    /// localization
    double a = mu - step;
    double b = mu + step;
    double fa = f(a), fb = f(b), fmu = f(mu);
    while (fa > fmu)
    {
       b = mu; fb = fmu;
       mu = a; fmu = fa;
       a -= step; fa = f(a);
    }
    while (fb > fmu)
    {
       a = mu;
       mu = b; fmu = fb;
       b += step; fb = f(b);
    }

    double root = 0;
    RandMath::findMin([this] (double x)
    {
        return -f(x);
    }, a, b, root);

    return root;
}

double ContinuousDistribution::getLeftLimit(double value, double epsilon) const
{
    if (!std::isfinite(f(value))) {
        if (std::fabs(value) < 1)
            value -= epsilon;
        else
            value *= 0.9999;
    }
    return value;
}

double ContinuousDistribution::getRightLimit(double value, double epsilon) const
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
        lowerBoundary = getLeftLimit(lowerBoundary, epsilon);
    }
    if (isLeftBounded()) {
        upperBoundary = std::min(maxPoint, upperBoundary);
        upperBoundary = getRightLimit(upperBoundary, epsilon);
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
    /// works good for unimodal and distributions
    static constexpr double epsilon1 = 1e-5;
    static constexpr double epsilon2 = 1e-10;
    static constexpr int maxIter = 1000;

    double lowerBoundary = startPoint, upperBoundary = startPoint;
    SUPPORT_TYPE suppType = supportType();
    if (suppType == FINITE_T)
    {
        lowerBoundary = getLeftLimit(MinValue(), epsilon2);
        upperBoundary = getRightLimit(MaxValue(), epsilon2);
    }
    else
    {
        int iter = 0;
        /// WARNING: we use variance - so there can be deadlock if we don't define this function explicitly
        /// therefore function Variance() should stay pure and noone should calculate it by this function
        double var = Variance();
        bool varIsInfinite = !std::isfinite(var);

        /// search lower boundary
        if (suppType == RIGHTSEMIFINITE_T) {
            lowerBoundary = getLeftLimit(MinValue(), epsilon2);
        }
        else if (varIsInfinite) {
            lowerBoundary = Quantile(0.001);
        }
        else
        {
            // TODO: check that lowerBoundary < upperBoundary
            /// get such lower boundary 'x' that f(x) < eps1 && |g(x)f(x)| < eps2 && F(x) < 0.001
            double fx = 1;
            do {
               lowerBoundary -= var;
               fx = f(lowerBoundary);
            } while ((fx > epsilon1 || std::fabs(funPtr(lowerBoundary)) * fx > epsilon2 || F(lowerBoundary) > 0.001) && ++iter < maxIter);

            if (iter == maxIter) /// can't take integral, integrand decreases too slow
                return NAN;
        }

        /// search upper boundary
        if (suppType == LEFTSEMIFINITE_T) {
            upperBoundary = getRightLimit(MaxValue(), epsilon2);
        }
        else if (varIsInfinite) {
            upperBoundary = Quantile(0.999);
        }
        else
        {
            /// get such upper boundary 'x' that f(x) < eps1 && |g(x)f(x)| < eps2 && F(x) > 0.999
            iter = 0;
            double fx = 1;
            do {
               upperBoundary += var;
               fx = f(upperBoundary);
            } while ((fx > epsilon1 || std::fabs(funPtr(upperBoundary)) * fx > epsilon2 || F(upperBoundary) < 0.999) && ++iter < maxIter);

            if (iter == maxIter) /// can't take integral, integrand decreases too slow
                return NAN;
        }
    }

    return RandMath::integral([this, funPtr] (double x)
    {
        return funPtr(x) * f(x);
    },
    lowerBoundary, upperBoundary);
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
