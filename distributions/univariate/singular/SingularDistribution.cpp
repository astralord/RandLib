#include "SingularDistribution.h"

SingularDistribution::SingularDistribution()
{

}

double SingularDistribution::Hazard(double) const
{
    return NAN;
}

double SingularDistribution::ExpectedValue(const std::function<double (double)> &, double, double ) const
{
    return NAN;
}

double SingularDistribution::Mode() const
{
    return NAN;
}

double SingularDistribution::Likelihood(const std::vector<double> &) const
{
    return NAN;
}

double SingularDistribution::LogLikelihood(const std::vector<double> &) const
{
    return NAN;
}
