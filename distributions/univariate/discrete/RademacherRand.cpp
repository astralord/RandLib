#include "RademacherRand.h"
#include "BernoulliRand.h"

RademacherRand::RademacherRand()
{
}

String RademacherRand::Name() const
{
    return "Rademacher";
}

double RademacherRand::P(const int & k) const
{
    return (k == 1 || k == -1) ? 0.5 : 0.0;
}

double RademacherRand::logP(const int & k) const
{
    return (k == 1 || k == -1) ? -M_LN2 : -INFINITY;
}

double RademacherRand::F(const int & k) const
{
    if (k < -1)
        return 0;
    return (k < 1) ? 0.5 : 1.0;
}

int RademacherRand::Variate() const
{
    return BernoulliRand::StandardVariate() ? 1 : -1;
}

double RademacherRand::Mean() const
{
    return 0;
}

double RademacherRand::Variance() const
{
    return 1;
}

std::complex<double> RademacherRand::CFImpl(double t) const
{
    return std::cos(t);
}

int RademacherRand::Median() const
{
    return -1;
}

int RademacherRand::Mode() const
{
    /// any from {-1, 1}
    return Variate();
}

double RademacherRand::Skewness() const
{
    return 0.0;
}

double RademacherRand::ExcessKurtosis() const
{
    return -2.0;
}

int RademacherRand::quantileImpl(double p) const
{
    return (p <= 0.5) ? -1 : 1;
}

int RademacherRand::quantileImpl1m(double p) const
{
    return (p >= 0.5) ? -1 : 1;
}

double RademacherRand::Entropy()
{
    return M_LN2;
}
