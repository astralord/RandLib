#include "RademacherRand.h"
#include "BernoulliRand.h"

RademacherRand::RademacherRand()
{
}

std::string RademacherRand::name()
{
    return "Rademacher";
}

double RademacherRand::P(int k) const
{
    if (k == -1 || k == 1)
        return 0.5;
    return 0;
}

double RademacherRand::F(double x) const
{
    if (x < -1)
        return 0;
    return (x < 1) ? 0.5 : 1.0;
}

double RademacherRand::variate() const
{
    return BernoulliRand::standardVariate() ? 1.0 : -1.0;
}

std::complex<double> RademacherRand::CF(double t) const
{
    return std::cos(t);
}

double RademacherRand::Median() const
{
    return 0.0;
}

double RademacherRand::Mode() const
{
    /// any from {-1, 1}
    return variate();
}

double RademacherRand::Skewness() const
{
    return 0.0;
}

double RademacherRand::ExcessKurtosis() const
{
    return -2.0;
}
