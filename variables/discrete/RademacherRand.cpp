#include "RademacherRand.h"

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
    if (x < 1)
        return 0.5;
    return 1;
}

double RademacherRand::variate() const
{
    if ((signed)RandGenerator::variate() < 0)
        return -1;
    return 1;
}

std::complex<double> RademacherRand::CF(double t) const
{
    return std::cos(t);
}

double RademacherRand::quantile(double p) const
{
    if (p == 0)
        return -INFINITY;
    if (p == 0.5)
        return -1;
    if (p == 1)
        return 1;
    return NAN;
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
