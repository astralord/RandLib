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

double SingularDistribution::LikelihoodFunction(const std::vector<double> &) const
{
    return NAN;
}

double SingularDistribution::LogLikelihoodFunction(const std::vector<double> &) const
{
    return NAN;
}
