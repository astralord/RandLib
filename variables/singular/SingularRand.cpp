#include "SingularRand.h"

SingularRand::SingularRand()
{

}

double SingularRand::Quantile(double p) const
{
    if (p < 0 || p > 1)
        return NAN;

    double root = Mean(); /// good starting point
    if (std::isnan(root) || std::isinf(root))
        root = 0.0;
    if (RandMath::findRoot([this, p] (double x)
    {
        return F(x) - p;
    }, root))
        return root;
    /// if we can't find quantile, then probably p == 1
    return INFINITY;
}

double SingularRand::Hazard(double) const
{
    return NAN;
}

double SingularRand::ExpectedValue(const std::function<double (double)> &, double) const
{
    return NAN;
}

double SingularRand::Mode() const
{
    return NAN;
}

